from pathlib import Path
import json, re, csv, hashlib
from plasma_surface.provenance import matches_recorded
from collections import Counter
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def xyz(path):
    lines = path.read_text().splitlines()
    n = int(lines[0])
    return (
        lines[1],
        np.array([[float(x) for x in l.split()[1:4]] for l in lines[2 : n + 2]]),
        Counter(l.split()[0] for l in lines[2 : n + 2]),
    )


def test_coadsorbate_provenance_stoichiometry_and_reference_energies():
    d = json.loads((ROOT / 'data/intermediate_campaign/refined_summary.json').read_text())
    for key in ('runner', 'helper', 'refinement_runner'):
        assert matches_recorded(ROOT / d[key], d[key + '_sha256'])
    assert len(d['results']) == 8
    for r in d['results']:
        folder = ROOT / r['folder']
        head, pos, composition = xyz(folder / 'final.extxyz')
        source_head, source_pos, parent_composition = xyz(ROOT / r['source'])
        expected = parent_composition + Counter(
            {'F': 1, 'H': 1} if r['coadsorbate'] == 'HF' else {'O': 1, 'H': 2}
        )
        assert composition == expected
        assert matches_recorded(folder / 'final.extxyz', r['final_sha256'])
        assert matches_recorded(folder / 'trajectory.extxyz', r['trajectory_sha256'])
        assert np.allclose(pos[r['fixed_indices']], source_pos[r['fixed_indices']], atol=1e-7)
        energy = float(re.search(r'\benergy=([^ ]+)', head).group(1))
        assert energy == pytest.approx(r['energy_eV'])
        assert r['incremental_association_energy_eV'] == pytest.approx(
            energy - r['reference']['energy_eV'] - r['gas_reference']['energy_eV']
        )
        assert r['converged'] and r['max_mobile_force_eV_A'] <= d['force_tolerance_eV_A']
        assert not r['rate_enabled']
        assert 'not a reaction path' in head


def test_candidate_network_retains_negative_curvature_and_has_no_fabricated_rates():
    d = json.loads((ROOT / 'data/intermediate_campaign/candidate_network.json').read_text())
    assert len(d['nodes']) == len(d['edges']) == 8
    assert any(n['curvature_screen']['minimum_eigenvalue_eV_A2'] < -0.02 for n in d['nodes'])
    for n in d['nodes']:
        assert len(n['symbols']) == len(n['positions_A'])
        assert all(max(i, j) < len(n['symbols']) for i, j in n['distance_graph'])
        assert matches_recorded(ROOT / n['structure'], n['structure_sha256'])
    assert all(
        e['rate_s'] is None and e['barrier_eV'] is None and not e['rate_enabled']
        for e in d['edges']
    )
