"""Screen all mobile host+adsorbate modes and escape unstable coadsorbate minima."""

from pathlib import Path
import json, hashlib, time
import numpy as np
import torch
from threadpoolctl import threadpool_limits
from ase.io import read, write
from ase.optimize import BFGS
from mace.calculators import MACECalculator
from molecular_surface_campaign import ROOT, independent, snapshot


def hessian(a, mobile, delta=0.01):
    xyz = a.positions.copy()
    h = np.zeros((3 * len(mobile), 3 * len(mobile)))
    for k, (i, j) in enumerate((i, j) for i in mobile for j in range(3)):
        p = xyz.copy()
        p[i, j] += delta
        a.set_positions(p)
        fp = a.get_forces(apply_constraint=False)[mobile].ravel()
        p = xyz.copy()
        p[i, j] -= delta
        a.set_positions(p)
        fm = a.get_forces(apply_constraint=False)[mobile].ravel()
        h[:, k] = -(fp - fm) / (2 * delta)
    a.set_positions(xyz)
    h = (h + h.T) / 2
    w, v = np.linalg.eigh(h)
    return w, v, h


def examine(path, calc):
    out = path.parent / 'full_stability'
    out.mkdir(exist_ok=False)
    a = read(path)
    a.calc = independent(calc)
    fixed = set(a.constraints[0].get_indices())
    mobile = [i for i in range(len(a)) if i not in fixed]
    start = time.perf_counter()
    attempts = []
    w, v, h = hessian(a, mobile)
    initial_w = w.tolist()
    for iteration in range(2):
        if w.min() >= -0.02:
            break
        choices = []
        for sign in (-1, 1):
            sub = out / f'escape_{iteration}_{sign}'
            sub.mkdir()
            b = a.copy()
            b.positions[mobile] += sign * 0.18 * v[:, 0].reshape(-1, 3)
            b.calc = independent(calc)
            frames = []
            opt = BFGS(b, logfile=str(sub / 'optimization.log'), maxstep=0.06)
            opt.attach(lambda: frames.append(snapshot(b)), interval=1)
            ok = bool(opt.run(fmax=0.01, steps=250))
            write(sub / 'trajectory.extxyz', frames)
            write(sub / 'final.extxyz', frames[-1])
            choices.append(
                (b.get_potential_energy(), b, ok, str(sub.relative_to(ROOT)).replace('\\', '/'))
            )
        valid = [x for x in choices if x[2]]
        attempts.append([dict(energy_eV=float(x[0]), converged=x[2], folder=x[3]) for x in choices])
        if not valid:
            break
        energy, a, ok, chosen = min(valid, key=lambda x: x[0])
        w, v, h = hessian(a, mobile)
    write(out / 'final.extxyz', snapshot(a))
    np.save(out / 'hessian_eV_A2.npy', h)
    result = dict(
        source=str(path.relative_to(ROOT)).replace('\\', '/'),
        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        mobile_indices=mobile,
        displacement_A=0.01,
        initial_eigenvalues_eV_A2=initial_w,
        final_eigenvalues_eV_A2=w.tolist(),
        minimum_eigenvalue_eV_A2=float(w.min()),
        appreciable_negative_modes=int(sum(w < -0.02)),
        max_mobile_force_eV_A=float(np.linalg.norm(a.get_forces(), axis=1).max()),
        energy_eV=float(a.get_potential_energy()),
        attempts=attempts,
        final_structure=str((out / 'final.extxyz').relative_to(ROOT)).replace('\\', '/'),
        final_sha256=hashlib.sha256((out / 'final.extxyz').read_bytes()).hexdigest(),
        hessian_sha256=hashlib.sha256((out / 'hessian_eV_A2.npy').read_bytes()).hexdigest(),
        elapsed_s=time.perf_counter() - start,
        scope='All mobile host and adsorbate Cartesian modes; fixed lower-host constraint retained; MACE only; finite-difference tolerance not a DFT minimum proof',
    )
    (out / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


def main():
    torch.set_num_threads(2)
    threadpool_limits(limits=1)
    data = ROOT / 'data/intermediate_campaign'
    d = json.loads((data / 'refined_summary.json').read_text())
    model = ROOT / '.cache/mace/mp_0b2_small.model'
    calc = MACECalculator(model_paths=str(model), device='cpu', default_dtype='float64')
    parents = {}
    results = []
    for r in d['results']:
        key = r['surface'] + f"_{r['start']}"
        if key not in parents:
            parents[key] = examine(ROOT / r['reference']['folder'] / 'final.extxyz', calc)
            print('parent', key, parents[key]['minimum_eigenvalue_eV_A2'], flush=True)
        checked = examine(ROOT / r['folder'] / 'final.extxyz', calc)
        checked.update(
            surface=r['surface'],
            coadsorbate=r['coadsorbate'],
            start=r['start'],
            parent_key=key,
            incremental_association_energy_eV=checked['energy_eV']
            - parents[key]['energy_eV']
            - r['gas_reference']['energy_eV'],
        )
        results.append(checked)
        manifest = dict(
            runner='scripts/validate_coadsorbate_minima.py',
            runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            model_sha256=d['model_sha256'],
            parents=parents,
            results=results,
            status='complete' if len(results) == 8 else 'running',
        )
        (data / 'full_stability.json').write_text(json.dumps(manifest, indent=2) + '\n')
        print(
            key,
            r['coadsorbate'],
            checked['minimum_eigenvalue_eV_A2'],
            checked['incremental_association_energy_eV'],
            flush=True,
        )


if __name__ == '__main__':
    main()
