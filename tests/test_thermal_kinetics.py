import math
import numpy as np
import pytest
from plasma_surface.thermal_kinetics import (
    KB_EV,
    H_EV_S,
    arrhenius,
    harmonic_free_energy,
    harmonic_tst,
    eyring,
    adsorption_hazard,
    isotropic_hop_diffusivity,
    md_count_rate,
)


def test_harmonic_tst_classical_limit_and_detailed_balance():
    # kBT >> h nu: harmonic TST tends to the Vineyard frequency ratio.
    modes = [1e-7, 2e-7, 3e-7]
    ts = [1.5e-7, 2.5e-7]
    T = 400.0
    result = harmonic_tst(0.5, T, modes, ts, 1)
    expected = np.prod(modes) / np.prod(ts) / H_EV_S * math.exp(-0.5 / (KB_EV * T))
    assert result['rate_s'] == pytest.approx(expected, rel=1e-5)
    forward = eyring(0.7, T)
    reverse = eyring(0.5, T)
    assert forward / reverse == pytest.approx(math.exp(-0.2 / (KB_EV * T)))
    assert harmonic_free_energy([0.1, 0.2], 0.01) == pytest.approx(0.15)


def test_thermal_rates_reject_missing_or_unstable_mode_data():
    for modes in ([0, 0.1], [-0.1, 0.2], [complex(0, 0.1)], [float('nan')]):
        with pytest.raises(ValueError):
            harmonic_free_energy(modes, 400.0)
    with pytest.raises(ValueError):
        harmonic_tst(0.5, 400.0, [0.1, 0.2], [0.1], 0)
    with pytest.raises(ValueError):
        harmonic_tst(0.5, 400.0, [0.1, 0.2], [0.1, 0.2], 1)
    assert arrhenius(0.5, 400, 0) == 0


def test_flux_hop_conventions_and_zero_md_event_bound():
    a = adsorption_hazard(1.0, 400.0, 20.0, 1e-19, 0.2)
    assert a['hazard_s'] == pytest.approx(a['flux_m2_s'] * 2e-20)
    # Four-neighbor square lattice: total escape rate is 4*k, D=k*l^2.
    assert isotropic_hop_diffusivity(3.0, 2e-10, 4) == pytest.approx(1.2e-19)
    bound = md_count_rate(0, 1e-10)
    assert bound['rate_s'] == 0 and bound['upper_rate_s'] == pytest.approx(-math.log(0.05) / 1e-10)
    assert not bound['qualified_for_kmc']


def test_parameter_audit_matches_live_inputs_and_enabled_event_coverage():
    from pathlib import Path
    import csv, json, hashlib
    from plasma_surface.provenance import matches_recorded
    from plasma_surface.species_kmc import rate_constants

    root = Path(__file__).resolve().parents[1]
    folder = root / 'data/kinetic_audit'
    m = json.loads((folder / 'manifest.json').read_text())
    for path, sha in m['input_sha256'].items():
        assert matches_recorded(root / path, sha)
    assert matches_recorded(root / m['runner'], m['runner_sha256'])
    network = json.loads((root / 'configs/species_kmc_network.json').read_text())
    rows = list(csv.DictReader((folder / 'species_event_rates.csv').open()))
    for T in m['temperatures_K']:
        rr = [r for r in rows if float(r['temperature_K']) == T]
        assert [r['event_id'] for r in rr] == [e['id'] for e in network['events']]
        assert np.allclose([float(r['rate_s']) for r in rr], rate_constants(network, T))
        assert all(not r['activation_free_energy_eV'] for r in rr)
    params = list(csv.DictReader((folder / 'parameters.csv').open()))
    assert len(params) == m['parameter_rows']
    assert all(r['value'] == '' for r in params if r['model'] == 'missing')
