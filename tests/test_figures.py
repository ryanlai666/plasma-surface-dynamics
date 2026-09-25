"""Presentation figures must trace back to the saved results they draw."""

import json
from pathlib import Path

import numpy as np
from PIL import Image

from plasma_surface.provenance import matches_recorded

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / 'docs/figures'
MANIFEST = json.loads((FIGURES / 'manifest.json').read_text(encoding='utf-8'))


def test_figures_were_drawn_by_the_current_code():
    assert matches_recorded('scripts/make_figures.py', MANIFEST['runner_sha256'])
    assert matches_recorded('plasma_surface/figstyle.py', MANIFEST['style_sha256'])


def test_figure_inputs_are_the_saved_results():
    assert MANIFEST['inputs_sha256']
    for path, sha in MANIFEST['inputs_sha256'].items():
        assert matches_recorded(path, sha), path


def test_every_figure_is_recorded_and_unchanged():
    on_disk = {f.name for f in FIGURES.iterdir() if f.suffix in ('.png', '.gif')}
    assert on_disk == set(MANIFEST['figures_sha256'])
    for name, sha in MANIFEST['figures_sha256'].items():
        assert matches_recorded(FIGURES / name, sha), name


def test_animations_show_every_saved_snapshot_and_nothing_else():
    lattice = np.load(ROOT / 'docs/species_kmc_results/lattice_trajectory.npz')
    with Image.open(FIGURES / 'species_kmc.gif') as gif:
        assert gif.n_frames == len(lattice['times'])
    trajectory = json.loads((ROOT / 'data/multilayer/demo_trajectory.json').read_text())
    with Image.open(FIGURES / 'multilayer_kmc.gif') as gif:
        assert gif.n_frames == len(trajectory['snapshots'])
