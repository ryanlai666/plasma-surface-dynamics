import json
from pathlib import Path
import math
import numpy as np
import pytest
from plasma_surface.dry_etch import f2_rate, first_event_kmc, KB_EV, KB_J

LAWS=json.loads((Path(__file__).resolve().parents[1]/'data/literature/si_f2_kinetics.json').read_text())['laws']

def test_source_reference_temperature_and_pressure_units():
    law=LAWS[0];t=298.15;p=133.32236842105263
    r=f2_rate(law,t,p)
    expected_k=8.82e-21*math.exp(-.13/(KB_EV*t))
    assert r['k_m3_s']==pytest.approx(expected_k)
    assert r['hazard_s']==pytest.approx(expected_k*p/(KB_J*t))
    assert f2_rate(law,t,2*p)['hazard_s']==pytest.approx(2*r['hazard_s'])
    assert f2_rate(law,t,0)['hazard_s']==0

def test_rejects_extrapolation_and_invalid_pressure():
    for t,p in [(200,1),(1001,1),(500,-1),(float('nan'),1)]:
        with pytest.raises(ValueError):f2_rate(LAWS[0],t,p)

def test_gillespie_first_event_conservation_and_limits():
    assert np.array_equal(first_event_kmc(0,50,[0,1,10],1),[0,0,0])
    x=first_event_kmc(2,50,[0,.01,.1,1,100],1)
    assert x[0]==0 and x[-1]==1 and np.all(np.diff(x)>=0)
    assert np.allclose(x*50,np.round(x*50))

def test_published_path_barriers_obey_same_saddle_reference():
    d=json.loads((Path(__file__).resolve().parents[1]/'data/literature/sicl4_surface_paths.json').read_text())
    for r in d['paths']:
        assert r['barrier_forward_eV']-r['barrier_reverse_eV']==pytest.approx(r['reaction_energy_eV'])
        assert r['barrier_forward_eV']>=0 and r['barrier_reverse_eV']>=0
