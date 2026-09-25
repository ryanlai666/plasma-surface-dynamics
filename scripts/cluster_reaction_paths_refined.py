"""Local Si-N/Si-O cleavage proxies: OMol25 NEB and optional direct DFT evaluation.

These capped molecular motifs supplement periodic surfaces; they are not bulk slabs.
"""

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import argparse, copy, csv, hashlib, json, time
import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from ase.io import read, write
from ase.optimize import BFGS, FIRE
from ase.mep import NEB
from ase.calculators.singlepoint import SinglePointCalculator

ROOT = Path(__file__).resolve().parents[1]


def clone(calc):
    c = copy.copy(calc)
    c.atoms = None
    c.results = {}
    return c


def freeze(a):
    e = float(a.get_potential_energy())
    f = a.get_forces(apply_constraint=False).copy()
    b = a.copy()
    b.calc = SinglePointCalculator(b, energy=e, forces=f)
    return b


def make_pair(motif, reactant):
    angles = np.arange(3) * 2 * np.pi / 3
    cap = [[0, 0, 0]] + [
        [1.48 * np.sqrt(8 / 9) * np.cos(x), 1.48 * np.sqrt(8 / 9) * np.sin(x), -1.48 / 3]
        for x in angles
    ]
    z = 1.74 if motif == 'SiN' else 1.64
    symbols = ['Si', 'H', 'H', 'H', 'N' if motif == 'SiN' else 'O']
    xyz = cap + [[0, 0, z]]
    if motif == 'SiN':
        symbols += ['H', 'H']
        xyz += [[0.82, 0, z + 0.58], [-0.82, 0, z + 0.58]]
    else:
        symbols += ['H']
        xyz += [[0.8, 0, z + 0.6]]
    halogen = 'F' if reactant == 'HF' else 'Cl'
    bond = 0.93 if halogen == 'F' else 1.29
    halide = len(symbols)
    proton = halide + 1
    symbols += [halogen, 'H']
    xyz += [[2.8, 0, z + 0.35], [2.8 - bond, 0, z + 0.35]]
    initial = Atoms(symbols, positions=xyz, cell=[24, 24, 24], pbc=False)
    final = initial.copy()
    shift = np.array([0, 0, 4.0 - z])
    final.positions[4:halide] += shift
    final.positions[halide] = [0, 0, 1.63 if halogen == 'F' else 2.08]
    final.positions[proton] = final.positions[4] + [0, 0.82, -0.58]
    for a in [initial, final]:
        a.set_constraint(FixAtoms(indices=[0, 1, 2, 3]))
        a.info['charge_spin'] = np.array([0, 1])
    return initial, final, halide, proton


def curvature(a, calc, indices, delta=0.008):
    b = a.copy()
    b.calc = clone(calc)
    h = np.zeros((len(indices) * 3, len(indices) * 3))
    for col, (i, j) in enumerate((i, j) for i in indices for j in range(3)):
        x = b.positions[i, j]
        b.positions[i, j] = x + delta
        fp = b.get_forces()[indices].reshape(-1).copy()
        b.positions[i, j] = x - delta
        fm = b.get_forces()[indices].reshape(-1).copy()
        b.positions[i, j] = x
        h[:, col] = -(fp - fm) / (2 * delta)
    return np.linalg.eigh((h + h.T) / 2)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--motif', choices=['SiN', 'SiO'], required=True)
    p.add_argument('--reactant', choices=['HF', 'HCl'], required=True)
    p.add_argument('--checkpoint', type=Path)
    p.add_argument('--evaluate-path', type=Path)
    args = p.parse_args()
    try:
        from threadpoolctl import threadpool_limits

        threadpool_limits(limits=1)
    except ImportError:
        pass
    if args.evaluate_path:
        from pyscf import gto, dft, lib
        from ase.calculators.calculator import Calculator, all_changes
        from ase.units import Hartree, Bohr

        lib.num_threads(2)

        class DFT(Calculator):
            implemented_properties = ['energy', 'forces']

            def calculate(self, atoms=None, properties=None, system_changes=all_changes):
                super().calculate(atoms, properties, system_changes)
                mol = gto.M(
                    atom=list(zip(atoms.get_chemical_symbols(), atoms.positions)),
                    basis='def2-svp',
                    charge=0,
                    spin=0,
                    unit='Angstrom',
                    verbose=0,
                )
                mf = dft.RKS(mol).density_fit()
                mf.xc = 'PBE'
                mf.grids.level = 3
                mf.conv_tol = 1e-9
                mf.max_cycle = 150
                e = mf.kernel()
                if not mf.converged:
                    raise RuntimeError('Unconverged DFT SCF')
                self.results = {
                    'energy': float(e * Hartree),
                    'forces': -mf.nuc_grad_method().kernel() * Hartree / Bohr,
                }

        calc = DFT()
        images = read(args.evaluate_path, ':')
        out = args.evaluate_path.parent / 'pbe_svp_on_omol_path'
        out.mkdir(exist_ok=False)
        es = []
        saved = []
        for i, a in enumerate(images):
            a.calc = clone(calc)
            b = freeze(a)
            saved.append(b)
            es.append(b.get_potential_energy())
            write(out / 'images.extxyz', saved)
            (out / 'energies.json').write_text(
                json.dumps(
                    dict(
                        energies_eV=es,
                        relative_energies_eV=(np.array(es) - es[0]).tolist(),
                        method='Direct density-fitted PBE/def2-SVP; grid3; SCF1e-9 Eh',
                        geometry_method='OMol25 path images, not DFT-optimized NEB',
                        all_completed_SCF_converged=True,
                    ),
                    indent=2,
                )
                + '\n'
            )
            print('DFT image', i, es[-1], flush=True)
        return
    import torch
    from deepmd.calculator import DP

    torch.set_num_threads(2)
    calc = DP(model=str(args.checkpoint.resolve()), head='OMol25')
    out = (
        ROOT
        / 'data/surface_paths'
        / (args.motif + '_capped_motif')
        / args.reactant
        / 'omol25_refined'
    )
    out.mkdir(parents=True, exist_ok=False)
    initial, final, halide, proton = make_pair(args.motif, args.reactant)
    mobile = list(range(4, len(initial)))
    original = out.parent / 'omol25_neb'
    initial = read(original / 'IS.extxyz')
    final = read(original / 'FS.extxyz')
    record = dict(
        motif=args.motif,
        reactant=args.reactant,
        status='endpoint_search',
        surface='Capped local bonding proxy; NOT a periodic nitride/oxide slab',
        reaction=(
            'H3SiNH2 + '
            + args.reactant
            + ' -> H3Si'
            + ('F' if args.reactant == 'HF' else 'Cl')
            + ' + NH3'
        )
        if args.motif == 'SiN'
        else (
            'H3SiOH + '
            + args.reactant
            + ' -> H3Si'
            + ('F' if args.reactant == 'HF' else 'Cl')
            + ' + H2O'
        ),
        mobile_indices=mobile,
        fixed_indices=[0, 1, 2, 3],
        model='DPA-3.3-1M OMol25',
        charge=0,
        multiplicity=1,
        checkpoint_sha256=hashlib.sha256(args.checkpoint.read_bytes()).hexdigest(),
        runner='scripts/cluster_reaction_paths_refined.py',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        assumptions=[
            'SiH3 cap fixed; surface-like local Si-N or Si-O cleavage proxy',
            'No bulk screening, periodicity, coverage or composition representation',
            'Not an independently calibrated SiN/SiO2 etch barrier',
        ],
    )

    def save():
        (out / 'summary.json').write_text(json.dumps(record, indent=2) + '\n')

    save()
    ends = []
    for name, a in [('IS', initial), ('FS', final)]:
        a.calc = clone(calc)
        opt = BFGS(a, logfile=str(out / (name + '.log')), maxstep=0.08)
        ok = bool(opt.run(fmax=0.015, steps=250))
        eig, vec = curvature(a, calc, mobile)
        for attempt in range(3):
            if eig.min() > -0.02:
                break
            trials = []
            for sign in [-1, 1]:
                b = a.copy()
                b.positions[mobile] += sign * 0.18 * vec[:, 0].reshape(-1, 3)
                b.calc = clone(calc)
                trial = BFGS(
                    b, logfile=str(out / (name + f'_escape_{attempt}_{sign}.log')), maxstep=0.06
                )
                converged = bool(trial.run(fmax=0.008, steps=200))
                trials.append((b.get_potential_energy(), b, converged))
            _, a, ok = min(trials, key=lambda x: x[0])
            eig, vec = curvature(a, calc, mobile)
        ends.append(freeze(a))
        write(out / (name + '.extxyz'), ends[-1])
        record[name] = dict(
            converged=ok,
            energy_eV=float(a.get_potential_energy()),
            Si_X_distance_A=float(a.get_distance(0, 4)),
            hessian_eigenvalues_eV_A2=eig.tolist(),
        )
        save()
    if not all(
        record[k]['converged'] and min(record[k]['hessian_eigenvalues_eV_A2']) > -0.02
        for k in ['IS', 'FS']
    ):
        record['status'] = 'endpoint_validation_failed'
        save()
        return
    if abs(record['IS']['Si_X_distance_A'] - record['FS']['Si_X_distance_A']) < 0.5:
        record['status'] = 'endpoints_collapsed'
        save()
        return
    images = [ends[0].copy()] + [ends[0].copy() for _ in range(9)] + [ends[1].copy()]
    for a in images:
        a.calc = clone(calc)
    neb = NEB(images, k=0.2, climb=False, method='improvedtangent')
    neb.interpolate(method='idpp', apply_constraint=True)
    opt = FIRE(neb, logfile=str(out / 'neb.log'), dt=0.025, maxstep=0.06)
    opt.attach(lambda: write(out / 'checkpoint.extxyz', images), interval=5)
    opt.run(fmax=0.10, steps=450)
    neb.climb = True
    opt = FIRE(neb, logfile=str(out / 'climb.log'), dt=0.02, maxstep=0.04)
    opt.attach(lambda: write(out / 'checkpoint.extxyz', images), interval=10)
    ok = bool(opt.run(fmax=0.025, steps=600))
    fmax = float(np.linalg.norm(neb.get_forces().reshape(-1, 3), axis=1).max())
    saved = [freeze(a) for a in images]
    write(out / 'images.extxyz', saved)
    es = np.array([a.get_potential_energy() for a in saved])
    peak = int(es.argmax())
    record.update(
        status='neb_converged' if ok else 'neb_unconverged',
        images=len(saved),
        peak_image=peak,
        peak_above_IS_eV=float(es.max() - es[0]),
        FS_minus_IS_eV=float(es[-1] - es[0]),
        max_neb_force_eV_A=fmax,
    )
    rows = [
        dict(
            image=i,
            energy_eV=float(e),
            relative_energy_eV=float(e - es[0]),
            Si_X_distance_A=float(a.get_distance(0, 4)),
            Si_halogen_distance_A=float(a.get_distance(0, halide)),
            halogen_H_distance_A=float(a.get_distance(halide, proton)),
        )
        for i, (e, a) in enumerate(zip(es, saved))
    ]
    with (out / 'energies.csv').open('w', newline='') as handle:
        w = csv.DictWriter(handle, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)
    if 0 < peak < len(images) - 1:
        eig, vec = curvature(saved[peak], calc, mobile)
        record['peak_hessian_eigenvalues_eV_A2'] = eig.tolist()
        record['index_one_in_mobile_subspace'] = bool(sum(eig < -0.02) == 1)
    save()
    print(json.dumps(record), flush=True)


if __name__ == '__main__':
    main()
