"""Exact fixed-frame, distance-constrained molecular minimization support."""

import numpy as np
from ase.constraints import FixConstraint


class FixedFrameDistances(FixConstraint):
    def __init__(self, fixed, pairs, targets):
        self.fixed = list(fixed)
        self.pairs = list(pairs)
        self.targets = np.asarray(targets, float)

    def get_removed_dof(self, atoms):
        return 3 * len(self.fixed) + len(self.pairs)

    def jacobian(self, positions):
        j = np.zeros((len(self.pairs), positions.size))
        values = []
        for k, (a, b) in enumerate(self.pairs):
            v = positions[a] - positions[b]
            r = np.linalg.norm(v)
            if r < 1e-8:
                raise ValueError('Coincident constrained atoms')
            values.append(r)
            if a not in self.fixed:
                j[k, 3 * a : 3 * a + 3] = v / r
            if b not in self.fixed:
                j[k, 3 * b : 3 * b + 3] = -v / r
        return np.array(values), j

    def adjust_positions(self, atoms, new):
        new[self.fixed] = atoms.positions[self.fixed]
        for _ in range(100):
            values, j = self.jacobian(new)
            residual = values - self.targets
            if abs(residual).max() < 1e-10:
                return
            new[:] -= (j.T @ np.linalg.solve(j @ j.T, residual)).reshape(-1, 3)
        raise RuntimeError('Distance projection did not converge')

    def adjust_forces(self, atoms, forces):
        forces[self.fixed] = 0.0
        _, j = self.jacobian(atoms.positions)
        forces[:] -= (j.T @ np.linalg.solve(j @ j.T, j @ forces.ravel())).reshape(-1, 3)
