"""Unit-explicit thermal-rate calculations; results do not bypass rate evidence gates."""
import math
import numpy as np
KB_EV=8.617333262145e-5
KB_J=1.380649e-23
H_EV_S=4.135667696e-15
AMU_KG=1.66053906660e-27


def positive(value,name,zero=False):
    if not math.isfinite(value) or (value<0 if zero else value<=0):
        raise ValueError(name+' must be finite and '+('nonnegative' if zero else 'positive'))
    return value


def arrhenius(barrier_eV,temperature_K,prefactor_s):
    """Conditional per-event hazard, not a gas-surface bimolecular coefficient."""
    positive(barrier_eV,'barrier',True);positive(temperature_K,'temperature');positive(prefactor_s,'prefactor',True)
    return prefactor_s*math.exp(-barrier_eV/(KB_EV*temperature_K))


def harmonic_free_energy(mode_energies_eV,temperature_K):
    """Quantum harmonic F including ZPE for explicitly selected positive real modes.

    Caller must treat translational/rotational/soft modes consistently. Imaginary
    and zero modes are rejected, not silently removed or replaced by a cutoff.
    """
    positive(temperature_K,'temperature')
    modes=np.asarray(mode_energies_eV)
    if np.iscomplexobj(modes) or modes.ndim!=1 or not np.isfinite(modes).all() or np.any(modes<=0):
        raise ValueError('Supply a one-dimensional list of positive real vibrational energies')
    modes=modes.astype(float);kT=KB_EV*temperature_K
    return float(np.sum(.5*modes+kT*np.log(-np.expm1(-modes/kT))))


def eyring(delta_free_energy_eV,temperature_K,transmission=1.):
    """Unimolecular same-standard-state TST rate; no gas pressure factor implied."""
    positive(temperature_K,'temperature')
    if not math.isfinite(delta_free_energy_eV):raise ValueError('Free energy must be finite')
    if not math.isfinite(transmission) or not 0<=transmission<=1:raise ValueError('Classical transmission must lie in [0,1]')
    return transmission*KB_EV*temperature_K/H_EV_S*math.exp(-delta_free_energy_eV/(KB_EV*temperature_K))


def harmonic_tst(electronic_barrier_eV,temperature_K,is_modes_eV,ts_real_modes_eV,ts_imaginary_modes,transmission=1.):
    """Fixed-volume surface HTST: Delta F from matched IS/TS stable modes.

    TS must have one unstable mode excluded, with otherwise matching degrees of
    freedom. Electronic minima/TS validation and applicability remain external.
    """
    positive(electronic_barrier_eV,'electronic barrier',True)
    if ts_imaginary_modes!=1 or len(is_modes_eV)!=len(ts_real_modes_eV)+1:
        raise ValueError('Require one TS imaginary mode and matched IS/TS mode counts')
    fi=harmonic_free_energy(is_modes_eV,temperature_K);ft=harmonic_free_energy(ts_real_modes_eV,temperature_K)
    barrier=electronic_barrier_eV+ft-fi
    return dict(electronic_barrier_eV=electronic_barrier_eV,delta_vibrational_free_energy_eV=ft-fi,activation_Helmholtz_free_energy_eV=barrier,rate_s=eyring(barrier,temperature_K,transmission),qualified_for_kmc=False)


def adsorption_hazard(pressure_Pa,gas_temperature_K,mass_amu,site_area_m2,sticking):
    """Ideal-gas impingement on a specified site area, with independent sticking."""
    positive(pressure_Pa,'pressure',True);positive(gas_temperature_K,'gas temperature');positive(mass_amu,'mass');positive(site_area_m2,'site area')
    if not math.isfinite(sticking) or not 0<=sticking<=1:raise ValueError('Sticking must lie in [0,1]')
    flux=pressure_Pa/math.sqrt(2*math.pi*mass_amu*AMU_KG*KB_J*gas_temperature_K)
    return dict(flux_m2_s=flux,hazard_s=flux*site_area_m2*sticking)


def isotropic_hop_diffusivity(per_neighbor_rate_s,hop_length_m,neighbors,dimension=2,correlation=1.):
    """D=f*z*k*l^2/(2*d); equivalent isotropic hops, rate is per neighbor."""
    positive(per_neighbor_rate_s,'hop rate',True);positive(hop_length_m,'hop length')
    if not isinstance(neighbors,int) or neighbors<1 or dimension not in (1,2,3):raise ValueError('Invalid topology')
    if not math.isfinite(correlation) or not 0<=correlation<=1:raise ValueError('Correlation factor must lie in [0,1]')
    return correlation*neighbors*per_neighbor_rate_s*hop_length_m**2/(2*dimension)


def md_count_rate(event_count,reactive_site_time_s,confidence=.95):
    """Count/exposure MLE and one-sided Poisson upper bound, not AIMD generation.

    Requires independent committed transitions and stationary eligible-site time;
    repeated barrier recrossings must not be counted as independent reactions.
    """
    from scipy.stats import chi2
    if not isinstance(event_count,int) or event_count<0:raise ValueError('Event count must be a nonnegative integer')
    positive(reactive_site_time_s,'eligible site-time')
    if not 0<confidence<1:raise ValueError('Confidence must lie in (0,1)')
    return dict(rate_s=event_count/reactive_site_time_s,upper_rate_s=float(chi2.ppf(confidence,2*(event_count+1))/(2*reactive_site_time_s)),one_sided_confidence=confidence,qualified_for_kmc=False)
