"""Compute rigid-surface halogen diffusion with MACE and climbing-image NEB.
No interpolated frames are published: only optimized, energy-evaluated images.
"""

from pathlib import Path
import argparse, hashlib, json, time, importlib.metadata, copy, sys
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
import torch
from ase import Atom
from ase.build import diamond100, diamond111
from ase.constraints import FixAtoms
from ase.io import read, write
from ase.mep import NEB
from ase.optimize import BFGS
from threadpoolctl import threadpool_limits
from mace.calculators import MACECalculator

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--surface', choices=['Si100', 'Si111'], required=True)
    ap.add_argument('--species', choices=['F', 'Cl'], required=True)
    ap.add_argument('--steps', type=int, default=250)
    ap.add_argument('--initial-images', type=Path)
    args = ap.parse_args()
    out = ROOT / 'data/surface_paths' / args.surface / args.species / 'mace_neb'
    out.mkdir(parents=True, exist_ok=True)
    model = ROOT / '.cache/mace/mp_0b2_small.model'
    torch.set_num_threads(2)
    threadpool_limits(limits=1)
    calc = MACECalculator(model_paths=str(model), device='cpu', default_dtype='float64')
    builder = diamond100 if args.surface == 'Si100' else diamond111
    slab = builder('Si', size=(3, 2, 6), a=5.43, vacuum=10, periodic=False)
    slab.pbc = [True, True, False]
    top = np.where(slab.positions[:, 2] > slab.positions[:, 2].max() - 0.01)[0]
    # Two translationally equivalent atop guesses, separated by one surface lattice vector.
    order = sorted(top, key=lambda i: tuple(slab.positions[i, :2]))
    i = order[0]
    xy = slab.positions[i, :2]
    candidates = [
        j
        for j in top
        if slab.positions[j, 0] > xy[0] + 1 and abs(slab.positions[j, 1] - xy[1]) < 0.01
    ]
    if not candidates:
        raise ValueError('No neighboring equivalent surface site')
    j = min(candidates, key=lambda j: np.linalg.norm(slab.positions[j, :2] - xy))
    ends = []
    endpoint_checks = []
    for label, k in [('IS', i), ('FS', j)]:
        a = slab.copy()
        a.append(
            Atom(
                args.species,
                slab.positions[k] + [0.35, 0.25, 1.65 if args.species == 'F' else 2.15],
            )
        )
        a.set_constraint(FixAtoms(indices=range(len(slab))))
        a.calc = calc
        opt = BFGS(a, logfile=str(out / (label + '.log')), maxstep=0.15)
        ok = bool(opt.run(fmax=0.008, steps=250))
        from validate_surface_neb import hessian_mobile

        eig = hessian_mobile(a, calc)
        if not ok or min(eig) <= 0:
            raise RuntimeError(f'Endpoint {label} is not a stable minimum: {eig}')
        write(out / (label + '.extxyz'), a)
        endpoint_checks.append(
            dict(
                label=label,
                converged=ok,
                max_force=float(np.linalg.norm(a.get_forces(), axis=1).max()),
                curvature_eigenvalues_eV_A2=eig,
            )
        )
        ends.append(a)
    images = [ends[0]] + [ends[0].copy() for _ in range(7)] + [ends[1]]
    if args.initial_images:
        images = read(args.initial_images, ':')
        if len(images) != 9:
            raise ValueError('Nine restart images required')
        if any(
            np.max(abs(images[n].positions - ends[k].positions)) > 1e-4
            for n, k in [(0, 0), (-1, 1)]
        ):
            raise ValueError('Restart endpoints differ from checked minima')
    for a in images:
        # Independent ASE result caches; immutable model weights are shared, not duplicated.
        c = copy.copy(calc)
        c.atoms = None
        c.results = {}
        a.calc = c
    neb = NEB(images, climb=False, k=0.2, method='improvedtangent')
    if not args.initial_images:
        neb.interpolate(apply_constraint=True)
    start = time.perf_counter()
    opt = BFGS(neb, logfile=str(out / 'neb.log'), maxstep=0.15)
    opt.attach(lambda: write(out / 'checkpoint.extxyz', images), interval=5)
    opt.run(fmax=0.12, steps=60)
    neb.climb = True
    converged = bool(opt.run(fmax=0.025, steps=args.steps))
    forces = neb.get_forces()
    residual = float(np.linalg.norm(forces.reshape(-1, 3), axis=1).max())
    energies = [float(a.get_potential_energy()) for a in images]
    positions = np.array([a.positions for a in images])
    lengths = np.linalg.norm(np.diff(positions, axis=0).reshape(len(images) - 1, -1), axis=1)
    arc = np.r_[0, np.cumsum(lengths)]
    peak = int(np.argmax(energies))
    rows = []
    for n, a in enumerate(images):
        rows.append(
            dict(
                image=n,
                arc_length_A=float(arc[n]),
                energy_eV=energies[n],
                relative_energy_eV=energies[n] - energies[0],
                max_unconstrained_force_eV_A=float(np.linalg.norm(a.get_forces(), axis=1).max()),
            )
        )
    write(out / 'images.extxyz', images)
    import csv

    with (out / 'energies.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)
    summary = dict(
        command=sys.argv,
        initial_images_sha256=hashlib.sha256(args.initial_images.read_bytes()).hexdigest()
        if args.initial_images
        else None,
        optimizer='BFGS',
        surface=args.surface,
        species=args.species,
        reaction=f'{args.species}* (site A) -> {args.species}* (site B)',
        method='MACE-MP-0b2 small; rigid substrate CI-NEB',
        status='converged'
        if converged and all(x['converged'] for x in endpoint_checks)
        else 'unconverged',
        images=len(images),
        peak_image=peak,
        peak_above_IS_eV=max(energies) - energies[0],
        FS_minus_IS_eV=energies[-1] - energies[0],
        max_neb_force_eV_A=residual,
        endpoint_checks=endpoint_checks,
        seconds=time.perf_counter() - start,
        model_sha256=hashlib.sha256(model.read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        versions={
            p: importlib.metadata.version(p) for p in ['ase', 'mace-torch', 'torch', 'numpy']
        },
        assumptions=[
            'Ideal unreconstructed diamond surface, a=5.43 angstrom',
            '36 Si substrate atoms held fixed; one halogen atom mobile',
            '3 x 2 surface repeat, six atomic layers; 10 angstrom vacuum each side',
            'Periodic x/y only; no hydrogen passivation',
            'Single coverage, cell and image count; no convergence study',
            'Bulk-trained ML potential; no independent DFT accuracy validation',
            'NEB peak is a candidate saddle, not a frequency-verified full-surface TS',
        ],
        rendering='Nine evaluated final NEB images only; no interpolated animation frames',
    )
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
