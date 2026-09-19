"""NIST Shomate ideal-gas thermochemistry with explicit reference and fit limits."""
from pathlib import Path
import json,math
R_J_MOL_K=8.31446261815324
KJ_MOL_PER_EV=96.48533212331002


def gas_thermo(species,temperature_K,partial_pressure_Pa,library=None):
    if library is None:
        library=json.loads((Path(__file__).resolve().parents[1]/'data/literature/nist_gas_shomate.json').read_text())
    r=library['species'][species];lo,hi=r['temperature_range_K']
    if not math.isfinite(temperature_K) or not lo<=temperature_K<=hi:raise ValueError('Temperature outside source fit range')
    if not math.isfinite(partial_pressure_Pa) or partial_pressure_Pa<=0:raise ValueError('Positive finite partial pressure required')
    A,B,C,D,E,F,G,H=r['coefficients'];t=temperature_K/1000.
    cp=A+B*t+C*t*t+D*t**3+E/t**2
    dh=A*t+B*t*t/2+C*t**3/3+D*t**4/4-E/t+F-H
    entropy=A*math.log(t)+B*t+C*t*t/2+D*t**3/3-E/(2*t*t)+G
    pressure=R_J_MOL_K*temperature_K*math.log(partial_pressure_Pa/library['standard_pressure_Pa'])/1000.
    correction=dh-temperature_K*entropy/1000.+pressure
    return dict(species=species,temperature_K=temperature_K,partial_pressure_Pa=partial_pressure_Pa,Cp_J_mol_K=cp,H_minus_H298_kJ_mol=dh,S_standard_J_mol_K=entropy,pressure_correction_eV=pressure/KJ_MOL_PER_EV,mu_minus_H298_eV=correction/KJ_MOL_PER_EV,mu_element_enthalpy_reference_eV=(H+correction)/KJ_MOL_PER_EV,source=r['url'],qualified_for_surface_rate=False)
