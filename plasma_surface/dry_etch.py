"""Published F2/Si first-reaction kinetics and surface barrier sensitivity.

This module predicts conversion of pre-existing susceptible sites, not EPC.
No ion-assisted threshold or stoichiometric layer-removal yield is inferred.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import numpy as np

KB_EV = 8.617333262e-5
KB_J = 1.380649e-23
ROOT = Path(__file__).resolve().parents[1]


def f2_rate(law, temperature_K, pressure_Pa):
    if not math.isfinite(temperature_K) or not 298.15 <= temperature_K <= 1000:
        raise ValueError('Published fit applies only at 298.15 <= T <= 1000 K')
    if not math.isfinite(pressure_Pa) or pressure_Pa < 0:
        raise ValueError('Pressure must be finite and nonnegative')
    ratio = temperature_K / law['reference_temperature_K']
    activation = math.exp(-law['barrier_eV'] / (KB_EV * temperature_K))
    k = law['A0_m3_s'] * ratio ** law['n1'] * activation
    gamma = law['gamma0'] * ratio ** law['n2'] * activation
    density = pressure_Pa / (KB_J * temperature_K)
    return dict(k_m3_s=k, impact_probability=gamma, density_m3=density,
                hazard_s=k * density)


def first_event_kmc(hazard_s, sites, observation_times, seed):
    """Gillespie loss of susceptible sites; each site reacts at most once."""
    times = np.asarray(observation_times, dtype=float)
    if not np.isfinite(hazard_s) or hazard_s < 0 or not isinstance(sites, int) or sites < 1:
        raise ValueError('Require finite nonnegative hazard and positive integer site count')
    if not np.isfinite(times).all() or np.any(times < 0) or np.any(np.diff(times) < 0):
        raise ValueError('Observation times must be finite, nonnegative and ordered')
    rng = np.random.default_rng(seed)
    remaining = sites
    t = 0.0
    events = []
    if hazard_s:
        while remaining:
            t += rng.exponential(1 / (remaining * hazard_s))
            if len(times) == 0 or t > times[-1]:
                break
            events.append(t)
            remaining -= 1
    return np.searchsorted(events, times, side='right') / sites


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out = ROOT / 'docs/dry_etch_results'
    out.mkdir(parents=True, exist_ok=True)
    source = ROOT / 'data/literature/si_f2_kinetics.json'
    laws = json.loads(source.read_text())['laws']
    temps = np.linspace(298.15, 1000, 121)
    pressure = 133.32236842105263  # 1 Torr, matching the source illustration
    records = []
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
    checks = []
    for law in laws:
        rates = [f2_rate(law, t, pressure) for t in temps]
        axs[0].semilogy(temps, [r['hazard_s'] for r in rates], label=law['facet'])
        for t, r in zip(temps, rates):
            records.append(dict(facet=law['facet'], temperature_K=float(t), pressure_Pa=pressure, **r))
        rate = f2_rate(law, 500., pressure)['hazard_s']
        dimensionless_times = np.linspace(0, 4, 41)
        times = dimensionless_times / rate
        ensemble = np.array([first_event_kmc(rate, 512, times, seed=1000+i) for i in range(128)])
        mean = ensemble.mean(axis=0)
        expected = -np.expm1(-dimensionless_times)
        se = np.sqrt(expected*(1-expected)/(512*128))
        # Ignore deterministic t=0 in standardized residuals.
        z = float(np.max(np.abs(mean[1:]-expected[1:])/se[1:]))
        checks.append(dict(facet=law['facet'],temperature_K=500.,hazard_s=rate,
                           sites=512,replicates=128,max_standard_errors=z,passed=z < 5,
                           seed_policy='1000..1127; common random numbers across facets'))
        axs[1].plot(times[1:], mean[1:], label=law['facet'])
        axs[1].plot(times[1:], expected[1:], '--', lw=.8, color='black', alpha=.4)
    axs[0].set(xlabel='Surface/gas temperature (K)', ylabel='First-reaction hazard per site (1/s)', title='Published F2/Si rate law at 1 Torr F2')
    axs[1].set(xlabel='Exposure time (s)', ylabel='Fraction of initial sites reacted', title='500 K: Gillespie ensemble and exact expectation', xscale='log', ylim=(0,1.02))
    for ax in axs: ax.grid(alpha=.2);ax.legend(title='Si facet')
    fig.text(.5,.01,'First F2 reaction on susceptible F-terminated sites; no layer removal, ion bombardment or ALE cycling implied.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.04,1,1]);fig.savefig(out/'f2_surface_kinetics.png',dpi=160);plt.close(fig)
    with (out/'f2_rates.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    # HF values used only as conditional sensitivity with explicitly assumed prefactors.
    with (ROOT/'data/literature/sin_hf_pathways.csv').open() as f: hf=list(csv.DictReader(f))
    sensitivity=[]
    for row in hf:
        ea=float(row['activation_energy_eV'])
        for t in (300.,400.,500.,600.):
            for nu in (1e11,1e12,1e13):
                sensitivity.append(dict(pathway=row['pathway'],temperature_K=t,
                    assumed_prefactor_s=nu,barrier_eV=ea,conditional_hazard_s=nu*math.exp(-ea/(KB_EV*t)),
                    interpretation='Sensitivity of occupied reaction complex only; reference conventions/prefactors unresolved'))
    with (out/'sin_hf_sensitivity.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(sensitivity[0]));w.writeheader();w.writerows(sensitivity)
    manifest=dict(status='passed' if all(c['passed'] for c in checks) else 'failed',
        source=json.loads(source.read_text())['source'],source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),checks=checks,
        scope='Numerical reproduction of a published dry-etch first-event rate law; not new DFT, independent experimental validation or calibrated ALE EPC',
        f2_rate_rows=len(records),hf_sensitivity_rows=len(sensitivity))
    (out/'verification.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))
    if manifest['status']!='passed':raise RuntimeError('Kinetic verification failed')

if __name__=='__main__': main()
