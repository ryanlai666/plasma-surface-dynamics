"""Search final fluorination in a constrained SiF3-NH2 molecular proxy."""

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "2")
from pathlib import Path
import argparse, json, hashlib, csv
import numpy as np
import torch

from ase.io import write
from ase.optimize import BFGS, FIRE
from ase.mep import NEB
from cluster_reaction_paths_refined import make_pair, clone, freeze, curvature, ROOT


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--checkpoint', type=Path, required=True)
    args = p.parse_args()
    from deepmd.calculator import DP

    torch.set_num_threads(2)
    calc = DP(model=str(args.checkpoint.resolve()), head='OMol25')
    out = ROOT / 'data/final_cleavage/SiF3_NH2_HF/omol25'
    out.mkdir(parents=True, exist_ok=False)
    initial, final, halide, proton = make_pair('SiN', 'HF')
    for a in (initial, final):
        symbols = a.get_chemical_symbols()
        symbols[1:4] = ['F'] * 3
        a.set_chemical_symbols(symbols)
        a.positions[1:4] *= 1.60 / 1.48
        a.info.update(
            calculation_kind='Constrained SiF3-NH2 final-cleavage proxy; not a periodic surface',
            charge_spin=np.array([0, 1]),
        )
    mobile = list(range(4, len(initial)))
    ends = []
    r = dict(
        reaction='SiF3NH2 + HF -> SiF4 + NH3',
        gap='A3',
        scope='Molecular proxy for terminal NH2 final Si-N cleavage; fixed SiF3 frame; not the graph environment or bare-N branch',
        model='DPA-3.3-1M OMol25',
        checkpoint_sha256=hashlib.sha256(args.checkpoint.read_bytes()).hexdigest(),
        runner='scripts/final_cleavage_proxy.py',
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        helper='scripts/cluster_reaction_paths_refined.py',
        helper_sha256=hashlib.sha256(
            (ROOT / 'scripts/cluster_reaction_paths_refined.py').read_bytes()
        ).hexdigest(),
        fixed_indices=[0, 1, 2, 3],
        mobile_indices=mobile,
        charge=0,
        multiplicity=1,
        rate_enabled=False,
        status='endpoint_search',
    )

    def save():
        (out / 'summary.json').write_text(json.dumps(r, indent=2) + '\n')

    save()
    for name, a in [('IS', initial), ('FS', final)]:
        a.calc = clone(calc)
        opt = BFGS(a, logfile=str(out / (name + '.log')), maxstep=0.06)
        ok = bool(opt.run(fmax=0.015, steps=300))
        eig, vec = curvature(a, calc, mobile)
        for attempt in range(2):
            if eig.min() >= -0.02:
                break
            trials = []
            for sign in [-1, 1]:
                b = a.copy()
                b.positions[mobile] += sign * 0.15 * vec[:, 0].reshape(-1, 3)
                b.calc = clone(calc)
                o = BFGS(
                    b, logfile=str(out / (name + f'_escape_{attempt}_{sign}.log')), maxstep=0.05
                )
                success = bool(o.run(fmax=0.015, steps=200))
                trials.append((b.get_potential_energy(), b, success))
            _, a, ok = min(trials, key=lambda x: x[0])
            eig, vec = curvature(a, calc, mobile)
        ends.append(freeze(a))
        write(out / (name + '.extxyz'), ends[-1])
        r[name] = dict(
            converged=ok,
            energy_eV=float(a.get_potential_energy()),
            Si_N_distance_A=float(a.get_distance(0, 4)),
            Si_incoming_F_distance_A=float(a.get_distance(0, halide)),
            HF_distance_A=float(a.get_distance(halide, proton)),
            N_incoming_H_distance_A=float(a.get_distance(4, proton)),
            max_mobile_force_eV_A=float(np.linalg.norm(a.get_forces(), axis=1).max()),
            hessian_eigenvalues_eV_A2=eig.tolist(),
        )
        save()
        print(name, r[name], flush=True)
    if not all(
        r[k]['converged'] and min(r[k]['hessian_eigenvalues_eV_A2']) >= -0.02 for k in ('IS', 'FS')
    ):
        r['status'] = 'endpoint_validation_failed'
        save()
        return
    if (
        r['IS']['Si_N_distance_A'] > 2.2
        or r['FS']['Si_N_distance_A'] < 2.8
        or r['FS']['N_incoming_H_distance_A'] > 1.3
    ):
        r['status'] = 'endpoint_identity_failed'
        save()
        return
    images = [ends[0].copy() for _ in range(8)] + [ends[1].copy()]
    for a in images:
        a.calc = clone(calc)
    neb = NEB(images, k=0.2, climb=False, method='improvedtangent')
    neb.interpolate(method='idpp', apply_constraint=True)
    opt = FIRE(neb, logfile=str(out / 'neb.log'), dt=0.03, maxstep=0.06)
    opt.attach(lambda: write(out / 'checkpoint.extxyz', [freeze(a) for a in images]), interval=25)
    opt.run(fmax=0.08, steps=350)
    neb.climb = True
    opt = FIRE(neb, logfile=str(out / 'climb.log'), dt=0.02, maxstep=0.04)
    opt.attach(lambda: write(out / 'checkpoint.extxyz', [freeze(a) for a in images]), interval=25)
    ok = bool(opt.run(fmax=0.03, steps=500))
    saved = [freeze(a) for a in images]
    write(out / 'images.extxyz', saved)
    es = np.array([a.get_potential_energy() for a in saved])
    peak = int(es.argmax())
    forces = float(np.linalg.norm(neb.get_forces().reshape(-1, 3), axis=1).max())
    r.update(
        status='neb_converged_pending_DFT_and_connectivity' if ok else 'neb_unconverged',
        peak_image=peak,
        peak_above_IS_eV=float(es.max() - es[0]),
        FS_minus_IS_eV=float(es[-1] - es[0]),
        max_neb_force_eV_A=forces,
    )
    with (out / 'energies.csv').open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(
            [
                'image',
                'energy_eV',
                'relative_to_IS_eV',
                'Si_N_A',
                'Si_incoming_F_A',
                'HF_A',
                'N_incoming_H_A',
            ]
        )
        w.writerows(
            (
                i,
                e,
                e - es[0],
                a.get_distance(0, 4),
                a.get_distance(0, halide),
                a.get_distance(halide, proton),
                a.get_distance(4, proton),
            )
            for i, (a, e) in enumerate(zip(saved, es))
        )
    if 0 < peak < len(images) - 1:
        eig, vec = curvature(saved[peak], calc, mobile)
        r['peak_curvature_eV_A2'] = eig.tolist()
        r['index_one_mobile_subspace'] = bool(sum(eig < -0.02) == 1)
    write(out / 'DFT_check_geometries.extxyz', [saved[0], saved[peak], saved[-1]])
    r['images_sha256'] = hashlib.sha256((out / 'images.extxyz').read_bytes()).hexdigest()
    save()
    print(r['status'], r.get('peak_above_IS_eV'), flush=True)


if __name__ == '__main__':
    main()
