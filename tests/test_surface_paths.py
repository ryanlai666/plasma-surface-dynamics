"""Audit committed reaction-path artifacts without requiring ML model downloads."""

from pathlib import Path
import csv, hashlib, json, shlex
from plasma_surface.provenance import matches_recorded
import numpy as np
import pytest
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
SYSTEMS = [('Si100', 'F'), ('Si100', 'Cl'), ('Si111', 'F')]


def frames(path):
    lines = path.read_text().splitlines()
    result = []
    i = 0
    while i < len(lines):
        n = int(lines[i])
        meta = dict(t.split('=', 1) for t in shlex.split(lines[i + 1]) if '=' in t)
        atoms = [line.split() for line in lines[i + 2 : i + 2 + n]]
        result.append(
            (meta, [r[0] for r in atoms], np.array([[float(v) for v in r[1:4]] for r in atoms]))
        )
        i += n + 2
    return result


@pytest.mark.parametrize('surface,species', SYSTEMS)
def test_evaluated_path_energy_and_constraints(surface, species):
    p = ROOT / 'data/surface_paths' / surface / species / 'mace_neb'
    s = json.loads((p / 'summary.json').read_text())
    v = json.loads((p / 'curvature.json').read_text())
    fs = frames(p / 'images.extxyz')
    with (p / 'energies.csv').open() as f:
        rows = list(csv.DictReader(f))
    assert len(fs) == len(rows) == s['images']
    es = np.array([float(r['energy_eV']) for r in rows])
    for (meta, symbols, xyz), r in zip(fs, rows):
        assert symbols == ['Si'] * 36 + [species]
        assert meta['pbc'] == 'T T F'
        assert abs(float(meta['energy']) - float(r['energy_eV'])) < 1e-5
        assert np.array_equal(xyz[:-1], fs[0][2][:-1])
        assert np.isfinite(xyz).all()
    assert np.linalg.norm(fs[-1][2][-1] - fs[0][2][-1]) > 1
    assert abs(s['peak_above_IS_eV'] - (es.max() - es[0])) < 1e-8
    assert abs(s['FS_minus_IS_eV'] - (es[-1] - es[0])) < 1e-8
    assert 0 < s['peak_image'] < len(fs) - 1
    assert s['status'] == 'converged' and s['max_neb_force_eV_A'] <= 0.025
    assert v['re_evaluation_max_energy_error_eV'] < 1e-5
    assert v['endpoint_minima_in_mobile_subspace']
    assert matches_recorded(ROOT / s['command'][0], s['script_sha256'])


@pytest.mark.parametrize('surface,species', SYSTEMS)
def test_animation_has_only_calculated_frames(surface, species):
    p = ROOT / 'data/surface_paths' / surface / species / 'mace_neb'
    m = json.loads((p / 'render_manifest.json').read_text())
    with Image.open(p / 'path.gif') as im:
        assert im.n_frames == m['frame_count']
    assert (
        m['frame_image_indices'] == list(range(m['frame_count'])) and m['interpolated_frames'] == 0
    )
    for file, key in [
        ('energies.csv', 'energies_sha256'),
        ('images.extxyz', 'coordinates_sha256'),
        ('path.gif', 'gif_sha256'),
    ]:
        assert matches_recorded(p / file, m[key])


def test_source_slideshows_contain_only_three_stationary_points():
    p = ROOT / 'data/surface_paths/Si100_c4x2/SiCl4/published'
    manifest = json.loads((p / 'manifest.json').read_text())
    for r in manifest['paths']:
        with Image.open(p / (r['path'].lower() + '.gif')) as im:
            assert im.n_frames == 3
        assert r['frames'] == 3 and r['interpolation'].startswith('None:')


def test_remote_hf_energy_unit_conversion_and_missing_final_states():
    rows = json.loads((ROOT / 'data/literature/remote_hf_2020.json').read_text())
    assert len(rows) == 4
    for r in rows:
        assert r['activation_eV'] == pytest.approx(r['activation_E_over_kB_K'] * 8.617333262145e-5)
        assert r['FS_relative_eV'] is None
