"""BFGS continuation of the saved, unconverged final-cleavage band."""

import os

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '2')
from pathlib import Path
import argparse, json, hashlib, csv
import numpy as np
import torch
from ase.io import read, write
from ase.optimize import BFGS
from ase.mep import NEB
from cluster_reaction_paths_refined import clone, freeze, curvature, ROOT


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--checkpoint', type=Path, required=True)
    args = p.parse_args()
    from deepmd.calculator import DP

    torch.set_num_threads(2)
    calc = DP(model=str(args.checkpoint.resolve()), head='OMol25')
    src = ROOT / 'data/final_cleavage/SiF3_NH2_HF/omol25/images.extxyz'
    out = src.parent.parent / 'omol25_bfgs'
    out.mkdir(exist_ok=False)
    images = read(src, ':')
    for a in images:
        a.calc = clone(calc)
    neb = NEB(images, k=0.15, climb=False, method='improvedtangent')
    opt = BFGS(neb, logfile=str(out / 'neb.log'), maxstep=0.04)
    opt.attach(lambda: write(out / 'checkpoint.extxyz', [freeze(a) for a in images]), interval=20)
    opt.run(fmax=0.08, steps=250)
    neb.climb = True
    opt = BFGS(neb, logfile=str(out / 'climb.log'), maxstep=0.025)
    opt.attach(lambda: write(out / 'checkpoint.extxyz', [freeze(a) for a in images]), interval=20)
    ok = bool(opt.run(fmax=0.025, steps=350))
    saved = [freeze(a) for a in images]
    write(out / 'images.extxyz', saved)
    e = np.array([x.get_potential_energy() for x in saved])
    peak = int(e.argmax())
    mobile = list(range(4, 9))
    w, v = curvature(saved[peak], calc, mobile)
    r = dict(
        runner='scripts/refine_final_cleavage_neb.py',
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source=str(src.relative_to(ROOT)).replace('\\', '/'),
        source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),
        model_sha256=hashlib.sha256(args.checkpoint.read_bytes()).hexdigest(),
        status='converged_band_pending_saddle_checks' if ok else 'neb_unconverged',
        max_neb_force_eV_A=float(np.linalg.norm(neb.get_forces().reshape(-1, 3), axis=1).max()),
        peak_image=peak,
        peak_above_IS_eV=float(e.max() - e[0]),
        FS_minus_IS_eV=float(e[-1] - e[0]),
        peak_eigenvalues_eV_A2=w.tolist(),
        index_one_mobile_subspace=bool(sum(w < -0.02) == 1),
        rate_enabled=False,
        scope='Fixed SiF3 molecular frame; bare-N environments not represented; no DFT stationary-point validation',
    )
    with (out / 'energies.csv').open('w', newline='') as f:
        wr = csv.writer(f)
        wr.writerow(['image', 'energy_eV', 'relative_to_IS_eV'])
        wr.writerows((i, x, x - e[0]) for i, x in enumerate(e))
    write(out / 'DFT_check_geometries.extxyz', [saved[0], saved[peak], saved[-1]])
    (out / 'summary.json').write_text(json.dumps(r, indent=2) + '\n')
    print(r, flush=True)


if __name__ == '__main__':
    main()
