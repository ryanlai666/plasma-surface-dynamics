"""Check actual path energies, endpoint/peak curvature and a second-model comparison."""

from pathlib import Path
import json, hashlib, csv
import numpy as np
import torch
from ase.io import read, write
from ase.calculators.singlepoint import SinglePointCalculator
from threadpoolctl import threadpool_limits
from mace.calculators import MACECalculator

ROOT = Path(__file__).resolve().parents[1]


def hessian_mobile(a, calc, delta=0.005):
    a = a.copy()
    a.calc = calc
    i = len(a) - 1
    h = np.zeros((3, 3))
    for j in range(3):
        x = a.positions[i, j]
        a.positions[i, j] = x + delta
        fp = a.get_forces()[i].copy()
        a.positions[i, j] = x - delta
        fm = a.get_forces()[i].copy()
        a.positions[i, j] = x
        h[:, j] = -(fp - fm) / (2 * delta)
    return np.linalg.eigvalsh((h + h.T) / 2).tolist()


def main():
    torch.set_num_threads(2)
    threadpool_limits(limits=1)
    paths = list(sorted((ROOT / 'data/surface_paths').glob('*/*/mace_neb/summary.json')))
    calc = MACECalculator(
        model_paths=str(ROOT / '.cache/mace/mp_0b2_small.model'),
        device='cpu',
        default_dtype='float64',
    )
    for path in paths:
        d = json.loads(path.read_text())
        folder = path.parent
        images = read(folder / 'images.extxyz', ':')
        curvature = {}
        with (folder / 'energies.csv').open() as f:
            expected = list(csv.DictReader(f))
        energy_errors = []
        # Shared ASE calculators must be frozen separately before multi-image export.
        # Re-evaluate each saved geometry, compare against the independently recorded CSV,
        # then attach an immutable calculator with that image's energy and raw forces.
        for a, row in zip(images, expected):
            a.calc = calc
            e = float(a.get_potential_energy())
            forces = a.get_forces(apply_constraint=False).copy()
            energy_errors.append(abs(e - float(row['energy_eV'])))
            a.calc = SinglePointCalculator(a, energy=e, forces=forces)
        if max(energy_errors) > 1e-5:
            raise ValueError('Re-evaluated path energy mismatch')
        write(folder / 'images.extxyz', images)
        write(folder / 'IS.extxyz', images[0])
        write(folder / 'FS.extxyz', images[-1])
        for name, i in [('IS', 0), ('peak', d['peak_image']), ('FS', len(images) - 1)]:
            curvature[name] = hessian_mobile(images[i], calc)
        peak = curvature['peak']
        endpoints_ok = all(min(curvature[s]) > 0 for s in ['IS', 'FS'])
        index_one = sum(v < -0.01 for v in peak) == 1 and sum(v > 0.01 for v in peak) == 2
        result = dict(
            re_evaluation_max_energy_error_eV=max(energy_errors),
            validator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            curvature_eigenvalues_eV_A2=curvature,
            finite_difference_A=0.005,
            endpoint_minima_in_mobile_subspace=endpoints_ok,
            peak_index_one_in_mobile_subspace=index_one,
            interpretation='Three-dimensional adsorbate-only Hessian with every Si fixed; does not validate full-surface vibrational stability or DFT accuracy.',
            eligible_constrained_saddle=d['status'] == 'converged' and endpoints_ok and index_one,
        )
        (folder / 'curvature.json').write_text(json.dumps(result, indent=2) + '\n')
        print(str(folder), result, flush=True)
    del calc
    model = ROOT / '.cache/mace/mpa_0_medium.model'
    calc = MACECalculator(model_paths=str(model), device='cpu', default_dtype='float64')
    for path in paths:
        folder = path.parent
        images = read(folder / 'images.extxyz', ':')
        es = []
        for a in images:
            a.calc = calc
            es.append(float(a.get_potential_energy()))
        rel = (np.array(es) - es[0]).tolist()
        (folder / 'second_model.json').write_text(
            json.dumps(
                dict(
                    model='MACE-MPA-0 medium',
                    model_sha256=hashlib.sha256(model.read_bytes()).hexdigest(),
                    energies_eV=es,
                    relative_energies_eV=rel,
                    sampled_peak_above_IS_eV=max(rel),
                    interpretation='Single-point evaluation on MP-0b2 images, not an independently optimized MPA-0 NEB or accuracy reference.',
                ),
                indent=2,
            )
            + '\n'
        )


if __name__ == '__main__':
    main()
