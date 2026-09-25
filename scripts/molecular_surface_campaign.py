"""Evaluate a molecular/surface matrix and retain every actual geometry and energy.

These are explicitly constrained approach curves, not transition-state searches.
Use run_molecular_reaction_neb.py for separately relaxed reactive paths.
"""

from pathlib import Path
import argparse
import copy
import csv
import hashlib
import importlib.metadata
import json
import time
import numpy as np
import torch
from threadpoolctl import threadpool_limits
from ase import Atoms
from ase.build import diamond100, diamond111, surface
from ase.constraints import FixAtoms
from ase.io import read, write
from ase.optimize import BFGS
from ase.calculators.singlepoint import SinglePointCalculator
from mace.calculators import MACECalculator

ROOT = Path(__file__).resolve().parents[1]
MOLECULES = ['HF', 'HCl', 'F2', 'Cl2', 'H2', 'H2O', 'CH3F', 'SiF4', 'SiCl4']
SURFACES = ['Si100', 'Si111', 'beta_Si3N4_001', 'alpha_quartz_001']
HEIGHTS = [6.5, 5.0, 4.0, 3.25, 2.75, 2.25, 1.9, 1.6, 1.3]


def molecule(name):
    diatomics = {
        'HF': ('FH', 0.92),
        'HCl': ('ClH', 1.28),
        'F2': ('F2', 1.42),
        'Cl2': ('Cl2', 1.99),
        'H2': ('H2', 0.75),
    }
    if name in diatomics:
        symbols, bond = diatomics[name]
        return Atoms(symbols, positions=[[0, 0, 0], [0, 0, bond]])
    if name == 'H2O':
        return Atoms('OH2', positions=[[0, 0, 0], [0.757, 0, 0.586], [-0.757, 0, 0.586]])
    angles = np.arange(3) * 2 * np.pi / 3
    if name == 'CH3F':
        xyz = [[0, 0, 0], [0, 0, 1.38]]
        xyz += [
            [1.09 * np.sqrt(8 / 9) * np.cos(a), 1.09 * np.sqrt(8 / 9) * np.sin(a), 1.38 + 1.09 / 3]
            for a in angles
        ]
        return Atoms('FCH3', positions=xyz)
    bond = 1.56 if name == 'SiF4' else 2.02
    xyz = [[0, 0, 0], [0, 0, -bond]]
    xyz += [
        [bond * np.sqrt(8 / 9) * np.cos(a), bond * np.sqrt(8 / 9) * np.sin(a), bond / 3]
        for a in angles
    ]
    return Atoms(name, positions=xyz)


def slab(name):
    if name in ('Si100', 'Si111'):
        builder = diamond100 if name == 'Si100' else diamond111
        atoms = builder('Si', size=(3, 2, 6), a=5.43, vacuum=12, periodic=False)
        source = {'construction': 'ASE ideal diamond cut, a=5.43 A; 3x2 repeat, six atomic layers'}
    else:
        cif = ROOT / (
            'docs/mace_results/2102550.cif'
            if name.startswith('beta')
            else 'data/reference/surface_bulk/9005018.cif'
        )
        bulk = read(cif)
        atoms = surface(bulk, (0, 0, 1), 2, vacuum=12)
        if name.startswith('alpha'):
            atoms = atoms.repeat((2, 1, 1))
        source = {
            'bulk_cif': str(cif.relative_to(ROOT)),
            'bulk_sha256': hashlib.sha256(cif.read_bytes()).hexdigest(),
            'source_url': 'https://www.crystallography.net/cod/'
            + ('2102550' if name.startswith('beta') else '9005018')
            + '.html',
            'construction': 'Unrelaxed (001) cleavage, two unit-cell layers; quartz doubled along a',
            'bulk_temperature_note': 'Quartz source structure measured at 398 K; this is not a finite-temperature simulation'
            if name.startswith('alpha')
            else 'Public beta-Si3N4 structure',
        }
    atoms.pbc = [True, True, False]
    atoms.set_constraint(FixAtoms(indices=range(len(atoms))))
    return atoms, source


def independent(calc):
    result = copy.copy(calc)
    result.atoms = None
    result.results = {}
    return result


def snapshot(atoms):
    e = float(atoms.get_potential_energy())
    f = atoms.get_forces(apply_constraint=False).copy()
    out = atoms.copy()
    out.calc = SinglePointCalculator(out, energy=e, forces=f)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--surfaces', nargs='+', choices=SURFACES, default=SURFACES)
    parser.add_argument('--species', nargs='+', choices=MOLECULES, default=MOLECULES)
    args = parser.parse_args()
    torch.set_num_threads(2)
    threadpool_limits(limits=1)
    checkpoint = ROOT / '.cache/mace/mp_0b2_small.model'
    calc = MACECalculator(model_paths=str(checkpoint), device='cpu', default_dtype='float64')
    gases = {}
    for name in args.species:
        a = molecule(name)
        a.calc = independent(calc)
        opt = BFGS(a, logfile=None, maxstep=0.08)
        ok = bool(opt.run(fmax=0.01, steps=150))
        if not ok:
            raise RuntimeError('Gas geometry did not converge: ' + name)
        a.positions[:, 2] -= a.positions[:, 2].min()
        gases[name] = a.copy()
    for surface_name in args.surfaces:
        base, source = slab(surface_name)
        base.calc = independent(calc)
        e_slab = float(base.get_potential_energy())
        top = base.positions[:, 2].max()
        si = [
            i for i, atom in enumerate(base) if atom.symbol == 'Si' and atom.position[2] > top - 2.0
        ]
        center = 0.5 * (base.cell[0] + base.cell[1])
        site = min(si, key=lambda i: np.linalg.norm(base.positions[i, :2] - center[:2]))
        for name in args.species:
            out = ROOT / 'data/surface_paths' / surface_name / name / 'mace_approach'
            if out.exists():
                raise FileExistsError(out)
            out.mkdir(parents=True)
            start = time.perf_counter()
            gas = gases[name].copy()
            gas.set_cell(base.cell)
            gas.pbc = base.pbc
            gas.calc = independent(calc)
            e_gas = float(gas.get_potential_energy())
            rows, images = [], []
            for n, height in enumerate(HEIGHTS):
                ads = gas.copy()
                ads.positions += [base.positions[site, 0], base.positions[site, 1], top + height]
                a = base + ads
                a.set_constraint(FixAtoms(indices=range(len(a))))
                a.calc = independent(calc)
                a.info.update(
                    surface=surface_name,
                    species=name,
                    height_A=height,
                    substrate_atoms=len(base),
                    calculation_kind='Rigid approach energy evaluation; not NEB or TS',
                )
                image = snapshot(a)
                e = image.get_potential_energy()
                f = image.get_forces(apply_constraint=False)
                rows.append(
                    dict(
                        image=n,
                        nearest_molecular_atom_height_A=height,
                        energy_eV=e,
                        interaction_energy_eV=e - e_slab - e_gas,
                        max_raw_force_eV_A=float(np.linalg.norm(f, axis=1).max()),
                    )
                )
                images.append(image)
            write(out / 'images.extxyz', images)
            write(out / 'slab.extxyz', snapshot(base))
            write(out / 'molecule.extxyz', snapshot(gas))
            with (out / 'energies.csv').open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=rows[0])
                writer.writeheader()
                writer.writerows(rows)
            minimum = min(rows, key=lambda row: row['interaction_energy_eV'])
            meta = dict(
                surface=surface_name,
                species=name,
                images=len(images),
                substrate_atoms=len(base),
                formula=images[0].get_chemical_formula(),
                surface_site_index=site,
                source=source,
                energy_reference='E(slab+molecule) - E(slab) - E(periodic molecular layer), same cell and fixed geometries',
                slab_energy_eV=e_slab,
                molecular_layer_energy_eV=e_gas,
                sampled_minimum_interaction_eV=minimum['interaction_energy_eV'],
                sampled_minimum_height_A=minimum['nearest_molecular_atom_height_A'],
                farthest_point_residual_eV=rows[0]['interaction_energy_eV'],
                model='MACE-MP-0b2 small',
                model_sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
                runner='scripts/molecular_surface_campaign.py',
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                versions={
                    p: importlib.metadata.version(p)
                    for p in ['ase', 'mace-torch', 'numpy', 'torch']
                },
                seconds=time.perf_counter() - start,
                TS_status='Not searched: scan maxima are NOT transition states',
                assumptions=[
                    'Ideal unpassivated rigid surface; no reconstruction, defects or thermal motion',
                    'One fixed molecule orientation and lateral site; finite periodic coverage',
                    'Gas molecule internally optimized first, then translated rigidly',
                    'No intermediate animation frame is synthesized: each displayed geometry has a model energy',
                    'Bulk-trained potential is unvalidated for these reactions; no ALE rates fitted',
                ],
            )
            (out / 'summary.json').write_text(json.dumps(meta, indent=2) + '\n')
            print(
                surface_name, name, 'sampled minimum', minimum['interaction_energy_eV'], flush=True
            )


if __name__ == '__main__':
    main()
