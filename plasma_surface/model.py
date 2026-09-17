"""Two-state surface model with identical mean-field and Gillespie kinetics.

Each column exposes one bare or modified site. A removal exposes a bare site.
No lateral interactions, transport, charging, shadowing, or chemical speciation.
"""
from dataclasses import asdict, dataclass
import math
import numpy as np


@dataclass(frozen=True)
class Parameters:
    site_density_m2: float = 7.0e18
    layer_nm: float = 0.136
    sticking: float = 0.25
    desorption_prefactor_s: float = 1.0e6
    desorption_barrier_ev: float = 0.65
    chemical_threshold_ev: float = 15.0
    physical_threshold_ev: float = 40.0
    chemical_yield_scale: float = 0.8
    physical_yield_scale: float = 0.15
    deposition_sticking: float = 0.2

    def __post_init__(self):
        for key, value in asdict(self).items():
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"{key} must be finite and nonnegative")
        if min(self.site_density_m2, self.layer_nm, self.chemical_threshold_ev,
               self.physical_threshold_ev) <= 0:
            raise ValueError("Density, layer thickness, and thresholds must be positive")
        if max(self.sticking, self.deposition_sticking) > 1:
            raise ValueError("Sticking probabilities must be <= 1")


@dataclass(frozen=True)
class Phase:
    name: str
    duration_s: float
    radical_flux_m2_s: float = 0.0
    ion_flux_m2_s: float = 0.0
    ion_energy_ev: float = 0.0
    precursor_flux_m2_s: float = 0.0
    temperature_k: float = 300.0

    def __post_init__(self):
        for key, value in asdict(self).items():
            if key != "name" and (not math.isfinite(value) or value < 0):
                raise ValueError(f"{key} must be finite and nonnegative")
        if self.temperature_k <= 0:
            raise ValueError("temperature_k must be positive")


def rates(phase: Phase, p: Parameters):
    """Per-site hazards [s^-1]: adsorption, desorption, chemistry, sputter, growth."""
    threshold = lambda e0, scale: scale * max(math.sqrt(phase.ion_energy_ev / e0) - 1, 0)
    return np.array([
        p.sticking * phase.radical_flux_m2_s / p.site_density_m2,
        p.desorption_prefactor_s * math.exp(-p.desorption_barrier_ev / (8.617333262e-5 * phase.temperature_k)),
        phase.ion_flux_m2_s / p.site_density_m2 * threshold(p.chemical_threshold_ev, p.chemical_yield_scale),
        phase.ion_flux_m2_s / p.site_density_m2 * threshold(p.physical_threshold_ev, p.physical_yield_scale),
        p.deposition_sticking * phase.precursor_flux_m2_s / p.site_density_m2,
    ])


def ale_recipe(energy_ev=28.0, dose_s=2.0, ion_s=2.0, temperature_k=300.0):
    return [
        Phase("modify", dose_s, radical_flux_m2_s=2e19, temperature_k=temperature_k),
        Phase("purge", 0.5, temperature_k=temperature_k),
        Phase("remove", ion_s, ion_flux_m2_s=2e19, ion_energy_ev=energy_ev, temperature_k=temperature_k),
        Phase("purge", 0.5, temperature_k=temperature_k),
    ]


def _check_run(phases, cycles):
    if not phases or not isinstance(cycles, int) or cycles < 1:
        raise ValueError("Provide at least one phase and a positive integer cycle count")


def mean_field(phases, p=None, cycles=5):
    """Exact integration; positive net removal means etching, negative means growth."""
    p = p or Parameters()
    _check_run(phases, cycles)
    theta = removed = deposited = elapsed = 0.0
    rows = []
    for cycle in range(1, cycles + 1):
        for phase in phases:
            a, d, c, s, g = rates(phase, p)
            k = a + d + c + s + g
            t = phase.duration_s
            if k:
                equilibrium = a / k
                relaxation = -math.expm1(-k * t) / k
                integral = equilibrium * t + (theta - equilibrium) * relaxation
                theta = equilibrium + (theta - equilibrium) * math.exp(-k * t)
            else:
                integral = theta * t
            removed += (c * integral + s * t) * p.layer_nm
            deposited += g * t * p.layer_nm
            elapsed += t
            rows.append(dict(cycle=cycle, phase=phase.name, time_s=elapsed,
                             coverage=float(theta), removed_nm=float(removed),
                             deposited_nm=float(deposited), net_removed_nm=float(removed-deposited)))
    return rows


def kmc(phases, p=None, cycles=5, sites=256, seed=42, max_events=2_000_000, observer=None):
    """Continuous-time Gillespie kMC on independent surface columns.

    Height variance is an independent-column statistic, not a feature profile.
    Optional observer receives copied states at phase starts, events, and ends;
    it draws no random numbers and does not change the simulation trajectory.
    """
    p = p or Parameters()
    _check_run(phases, cycles)
    if not isinstance(sites, int) or sites < 1 or max_events < 1:
        raise ValueError("sites and max_events must be positive integers")
    rng = np.random.default_rng(seed)
    modified = np.zeros(sites, dtype=bool)
    heights = np.zeros(sites, dtype=int)
    removed = deposited = events = 0
    elapsed = 0.0
    rows = []
    def emit(cycle, phase, time):
        if observer is not None:
            observer(dict(cycle=cycle, phase=phase.name, time_s=float(time),
                          modified=modified.copy(), heights=heights.copy(),
                          removed=removed, deposited=deposited, events=events))
    for cycle in range(1, cycles + 1):
        for phase in phases:
            r = rates(phase, p)
            t = 0.0
            emit(cycle, phase, elapsed)
            while t < phase.duration_s:
                n = int(modified.sum())
                hazards = r * np.array([sites-n, n, n, sites, sites])
                total = hazards.sum()
                if total == 0:
                    break
                t += rng.exponential(1 / total)
                if t > phase.duration_s:
                    break
                events += 1
                if events > max_events:
                    raise RuntimeError("Event budget exceeded; reduce flux, sites, or duration")
                event = int(np.searchsorted(np.cumsum(hazards), rng.random() * total, side="right"))
                if event < 3:
                    candidates = np.flatnonzero(~modified if event == 0 else modified)
                    site = int(rng.choice(candidates))
                else:
                    site = int(rng.integers(sites))
                modified[site] = event == 0
                if event in (2, 3):
                    heights[site] -= 1
                    removed += 1
                elif event == 4:
                    heights[site] += 1
                    deposited += 1
                emit(cycle, phase, elapsed+t)
            elapsed += phase.duration_s
            emit(cycle, phase, elapsed)
            rows.append(dict(cycle=cycle, phase=phase.name, time_s=elapsed,
                             coverage=float(modified.mean()), removed_nm=removed*p.layer_nm/sites,
                             deposited_nm=deposited*p.layer_nm/sites,
                             net_removed_nm=(removed-deposited)*p.layer_nm/sites,
                             roughness_nm=float(heights.std()*p.layer_nm), events=events))
    return rows
