"""Relaxed two-coordinate diagnostic for local cleavage, not a TS calculation."""

import os

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '2')
from pathlib import Path
import argparse, json, hashlib
import numpy as np
import torch
from ase.io import read, write
from ase.constraints import FixAtoms, FixBondLengths
from ase.optimize import BFGS
from cluster_reaction_paths_refined import clone, freeze, ROOT


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--checkpoint', type=Path, required=True)
    args = p.parse_args()
    from deepmd.calculator import DP

    torch.set_num_threads(2)
    calc = DP(model=str(args.checkpoint.resolve()), head='OMol25')
    src = ROOT / 'data/final_cleavage/SiF3_NH2_HF/omol25/IS.extxyz'
    out = src.parent.parent / 'local_coordinate_scan'
    out.mkdir(exist_ok=False)
    a = read(src)
    a.calc = clone(calc)
    start = a.get_distance(0, 4)
    nh = a.get_distance(4, 8)
    rows = []
    frames = []
    targets = list(zip(np.linspace(start, 3.3, 11), np.linspace(nh, 1.025, 11)))
    for direction, seq in [
        ('forward', list(enumerate(targets))),
        ('reverse', list(enumerate(targets))[::-1]),
    ]:
        for i, (sin, nhdist) in seq:
            # Translate N and its original H atoms together before imposing the new distances.
            a.set_constraint()
            v = a.positions[4] - a.positions[0]
            shift = v * (sin / np.linalg.norm(v) - 1)
            a.positions[[4, 5, 6]] += shift
            v = a.positions[8] - a.positions[4]
            a.positions[8] = a.positions[4] + v * nhdist / np.linalg.norm(v)
            a.set_constraint([FixAtoms(indices=[0, 1, 2, 3]), FixBondLengths([(0, 4), (4, 8)])])
            opt = BFGS(a, logfile=str(out / f'{direction}_{i:02d}.log'), maxstep=0.035)
            ok = bool(opt.run(fmax=0.04, steps=160))
            physical = a.get_forces(apply_constraint=False)
            projected = a.get_forces()
            b = freeze(a)
            b.info.update(
                calculation_kind='Evaluated relaxed two-coordinate scan; not a minimum-energy path or TS',
                direction=direction,
                target_SiN_A=float(sin),
                target_NH_A=float(nhdist),
            )
            frames.append(b)
            rows.append(
                dict(
                    direction=direction,
                    index=i,
                    energy_eV=float(a.get_potential_energy()),
                    target_SiN_A=float(sin),
                    target_NH_A=float(nhdist),
                    SiN_A=float(a.get_distance(0, 4)),
                    NH_A=float(a.get_distance(4, 8)),
                    SiF_A=float(a.get_distance(0, 7)),
                    HF_A=float(a.get_distance(7, 8)),
                    projected_force_eV_A=float(np.linalg.norm(projected, axis=1).max()),
                    mobile_physical_force_eV_A=float(np.linalg.norm(physical[4:], axis=1).max()),
                    converged=ok,
                )
            )
            write(out / 'evaluated_geometries.extxyz', frames)
            (out / 'summary.json').write_text(
                json.dumps(
                    dict(
                        runner='scripts/scan_final_cleavage_coordinates.py',
                        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                        source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),
                        model_sha256=hashlib.sha256(args.checkpoint.read_bytes()).hexdigest(),
                        status='complete' if len(rows) == 22 else 'running',
                        scope='Constrained local Si-N elongation and proton transfer; forward/reverse hysteresis diagnostic; no TS or kinetic barrier',
                        rate_enabled=False,
                        results=rows,
                    ),
                    indent=2,
                )
                + '\n'
            )
            print(direction, i, rows[-1], flush=True)


if __name__ == '__main__':
    main()
