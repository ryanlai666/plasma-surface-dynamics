"""Resume saved molecular NEBs with FIRE and independent threaded calculators."""

from pathlib import Path
import argparse, csv, hashlib, json, time
import numpy as np
import torch
from threadpoolctl import threadpool_limits
from ase.io import read, write
from ase.mep import NEB
from ase.optimize import FIRE, BFGS
from ase.geometry import find_mic
from mace.calculators import MACECalculator
from molecular_surface_campaign import ROOT, independent, snapshot
from run_molecular_reaction_neb import curvature


def main():
    p = argparse.ArgumentParser()
    p.add_argument('species')
    args = p.parse_args()
    torch.set_num_threads(1)
    threadpool_limits(limits=1)
    source = ROOT / 'data/surface_paths/beta_Si3N4_001' / args.species / 'molecular_neb'
    out = source.parent / 'molecular_refined'
    out.mkdir(exist_ok=False)
    d = json.loads((source / 'summary.json').read_text())
    d['previous_runner'] = d['runner']
    d['previous_script_sha256'] = d['script_sha256']
    d.update(
        runner='scripts/refine_molecular_surface_paths.py',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        status='path_optimization',
        restart_source_sha256=hashlib.sha256(
            (source / 'checkpoint.extxyz').read_bytes()
        ).hexdigest(),
        optimizer='FIRE; independent ASE calculators, shared immutable MACE weights; threaded NEB',
    )

    def save():
        (out / 'summary.json').write_text(json.dumps(d, indent=2) + '\n')

    save()
    checkpoint = ROOT / '.cache/mace/mp_0b2_small.model'
    calc = MACECalculator(model_paths=str(checkpoint), device='cpu', default_dtype='float64')
    (out / 'restart.extxyz').write_bytes((source / 'checkpoint.extxyz').read_bytes())
    images = read(out / 'restart.extxyz', ':')
    mobile = d['mobile_indices']
    for a in images:
        a.calc = independent(calc)
    neb = NEB(images, k=0.25, climb=False, method='improvedtangent', parallel=True)
    opt = FIRE(neb, logfile=str(out / 'neb.log'), dt=0.03, maxstep=0.06)
    opt.attach(lambda: write(out / 'checkpoint.extxyz', images), interval=10)
    opt.run(fmax=0.12, steps=200)
    neb.climb = True
    opt = FIRE(neb, logfile=str(out / 'climb.log'), dt=0.025, maxstep=0.05)
    opt.attach(lambda: write(out / 'checkpoint.extxyz', images), interval=10)
    ok = bool(opt.run(fmax=0.035, steps=500))
    residual = float(np.linalg.norm(neb.get_forces().reshape(-1, 3), axis=1).max())
    saved = [snapshot(a) for a in images]
    write(out / 'images.extxyz', saved)
    es = np.array([a.get_potential_energy() for a in saved])
    peak = int(es.argmax())
    xyz = np.array([a.positions for a in images])
    arc = np.r_[
        0, np.cumsum(np.linalg.norm(np.diff(xyz, axis=0).reshape(len(images) - 1, -1), axis=1))
    ]
    rows = [
        dict(
            image=i,
            arc_length_A=float(arc[i]),
            energy_eV=float(e),
            relative_energy_eV=float(e - es[0]),
        )
        for i, e in enumerate(es)
    ]
    with (out / 'energies.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)
    d.update(
        status='neb_converged' if ok else 'neb_unconverged',
        images=len(images),
        peak_image=peak,
        peak_above_IS_eV=float(es.max() - es[0]),
        FS_minus_IS_eV=float(es[-1] - es[0]),
        max_neb_force_eV_A=residual,
    )
    save()
    if 0 < peak < len(images) - 1:
        eig, vec = curvature(saved[peak], calc, mobile)
        d['peak_hessian_eigenvalues_eV_A2'] = eig.tolist()
        d['index_one_in_mobile_subspace'] = bool(sum(eig < -0.02) == 1)
        save()
        if ok and d['index_one_in_mobile_subspace']:
            connections = []
            for sign in [-1, 1]:
                a = saved[peak].copy()
                a.positions[mobile] += sign * 0.12 * vec[:, 0].reshape(-1, 3)
                a.calc = independent(calc)
                opt = BFGS(a, logfile=str(out / f'downhill_{sign}.log'), maxstep=0.08)
                done = bool(opt.run(fmax=0.02, steps=220))
                write(out / f'downhill_{sign}.extxyz', snapshot(a))
                errors = [
                    float(
                        np.sqrt(
                            np.mean(
                                find_mic(a.positions[mobile] - b.positions[mobile], a.cell, a.pbc)[
                                    0
                                ]
                                ** 2
                            )
                        )
                    )
                    for b in [saved[0], saved[-1]]
                ]
                connections.append(dict(sign=sign, converged=done, endpoint_RMS_errors_A=errors))
            d['downhill_connections'] = connections
            matched = [
                int(np.argmin(c['endpoint_RMS_errors_A']))
                if min(c['endpoint_RMS_errors_A']) < 0.25
                else -1
                for c in connections
            ]
            d['endpoint_connectivity_confirmed'] = sorted(matched) == [0, 1] and all(
                c['converged'] for c in connections
            )
    save()
    print(json.dumps(d), flush=True)


if __name__ == '__main__':
    main()
