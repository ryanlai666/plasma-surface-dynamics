import importlib.util
from pathlib import Path
import numpy as np
import pytest

pytest.importorskip('ase')
from ase import Atoms

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    'fixed_frame', ROOT / 'scripts/fixed_frame_constraints.py'
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_position_projection_preserves_frame_and_distances():
    a = Atoms('SiNHF', positions=[[0, 0, 0], [0, 0, 2], [0.8, 0, 2.6], [2.5, 0, 0.5]])
    c = module.FixedFrameDistances([0], [(0, 1), (1, 2), (0, 3)], [2, 1, 2.5])
    new = a.positions.copy() + np.random.default_rng(1).normal(0, 0.1, (4, 3))
    c.adjust_positions(a, new)
    assert np.array_equal(new[0], a.positions[0])
    assert np.allclose(c.jacobian(new)[0], [2, 1, 2.5], atol=1e-10)


def test_force_projection_is_tangent_and_preserves_virtual_work():
    a = Atoms('SiNHF', positions=[[0, 0, 0], [0, 0, 2], [0.8, 0, 2.6], [2.5, 0, 0.5]])
    c = module.FixedFrameDistances([0], [(0, 1), (1, 2), (0, 3)], [2, 1, np.sqrt(6.5)])
    rng = np.random.default_rng(2)
    force = rng.normal(size=(4, 3))
    raw = force.copy()
    c.adjust_forces(a, force)
    _, j = c.jacobian(a.positions)
    assert np.allclose(force[0], 0) and np.allclose(j @ force.ravel(), 0, atol=1e-12)
    tangent = rng.normal(size=(4, 3))
    c.adjust_forces(a, tangent)
    assert np.sum(tangent * raw) == pytest.approx(np.sum(tangent * force), abs=1e-12)
