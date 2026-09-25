"""Trace both sides of a constrained molecular saddle to the minima they reach.

These are evaluated minimization trajectories, not interpolated movies or an IRC.
"""

from pathlib import Path
import csv, hashlib, json
import numpy as np
import torch
from threadpoolctl import threadpool_limits
from ase.io import read, write
from ase.optimize import BFGS
from mace.calculators import MACECalculator
from molecular_surface_campaign import ROOT, independent, snapshot
from run_molecular_reaction_neb import curvature


def main():
    torch.set_num_threads(1)
    threadpool_limits(limits=1)
    source = ROOT / 'data/surface_paths/beta_Si3N4_001/H2O/molecular_refined'
    out = source.parent / 'mace_local_saddle'
    out.mkdir(exist_ok=False)
    prior = json.loads((source / 'summary.json').read_text())
    assert prior['status'] == 'neb_converged' and prior['index_one_in_mobile_subspace']
    calc = MACECalculator(
        model_paths=str(ROOT / '.cache/mace/mp_0b2_small.model'),
        device='cpu',
        default_dtype='float64',
    )
    peak = read(source / 'images.extxyz', prior['peak_image'])
    mobile = prior['mobile_indices']
    eig, vec = curvature(peak, calc, mobile)
    branches = []
    endpoints = []
    for sign in [1, -1]:
        a = peak.copy()
        a.positions[mobile] += sign * 0.12 * vec[:, 0].reshape(-1, 3)
        a.calc = independent(calc)
        saved = []

        def record():
            saved.append(snapshot(a))

        opt = BFGS(a, logfile=str(out / f'branch_{sign}.log'), maxstep=0.06)
        opt.attach(record, interval=4)
        ok = bool(opt.run(fmax=0.008, steps=250))
        record()
        ev, _ = curvature(a, calc, mobile)
        branches.append(saved)
        endpoints.append(
            dict(
                converged=ok,
                energy_eV=float(a.get_potential_energy()),
                hessian_eigenvalues_eV_A2=ev.tolist(),
                no_negative_curvature=bool(ev.min() > -0.02),
            )
        )
        write(out / f'branch_{sign}.extxyz', saved)
    peak.calc = independent(calc)
    images = list(reversed(branches[0])) + [snapshot(peak)] + branches[1]
    write(out / 'images.extxyz', images)
    es = np.array([a.get_potential_energy() for a in images])
    ts = len(branches[0])
    assert int(es.argmax()) == ts
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
            state='TS'
            if i == ts
            else 'IS'
            if i == 0
            else 'FS'
            if i == len(images) - 1
            else 'evaluated relaxation',
        )
        for i, e in enumerate(es)
    ]
    with (out / 'energies.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)
    confirmed = (
        all(e['converged'] and e['no_negative_curvature'] for e in endpoints)
        and np.linalg.norm(images[-1].positions - images[0].positions) > 0.3
        and sum(eig < -0.02) == 1
    )
    d = dict(
        surface='beta_Si3N4_001',
        species='H2O',
        status='constrained_saddle_connected' if confirmed else 'local_minimum_validation_failed',
        method='MACE-MP-0b2; evaluated downhill minimization branches from a CI-NEB saddle; not an IRC',
        source=str(source.relative_to(ROOT)),
        source_sha256=hashlib.sha256((source / 'images.extxyz').read_bytes()).hexdigest(),
        runner='scripts/trace_surface_saddle.py',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        images=len(images),
        peak_image=ts,
        peak_above_IS_eV=float(es.max() - es[0]),
        FS_minus_IS_eV=float(es[-1] - es[0]),
        IS=endpoints[0],
        FS=endpoints[1],
        mobile_indices=mobile,
        peak_hessian_eigenvalues_eV_A2=eig.tolist(),
        index_one_in_mobile_subspace=bool(sum(eig < -0.02) == 1),
        endpoint_connectivity_confirmed=bool(confirmed),
        original_precursor_connectivity_confirmed=False,
        assumptions=[
            'IS is the alternative intact-water minimum actually reached from this saddle, not the original intended precursor',
            'Only two substrate atoms and water are mobile',
            'Reversed minimization frames are a connectivity visualization, not a dynamical trajectory or a minimum-energy-path optimization',
            'Bulk-trained MACE is unvalidated for this surface chemistry',
        ],
    )
    (out / 'summary.json').write_text(json.dumps(d, indent=2) + '\n')
    print(json.dumps(d), flush=True)


if __name__ == '__main__':
    main()
