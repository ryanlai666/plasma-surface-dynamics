"""Direct PBE/def2-SVP energies and forces on immutable saved molecular geometries."""

from pathlib import Path
import argparse, hashlib, json, csv
import numpy as np
from ase.io import read, write
from ase.calculators.singlepoint import SinglePointCalculator
from ase.units import Hartree, Bohr
from pyscf import gto, dft, lib

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input', type=Path)
    p.add_argument('output', type=Path)
    args = p.parse_args()
    lib.num_threads(2)
    args.output.mkdir(parents=True, exist_ok=False)
    source = args.input.read_bytes()
    (args.output / 'input.extxyz').write_bytes(source)
    images = read(args.output / 'input.extxyz', ':')
    meta = dict(
        method='Density-fitted PBE/def2-SVP',
        grid_level=3,
        SCF_tolerance_Hartree=1e-9,
        charge=0,
        multiplicity=1,
        geometry_method='Fixed saved ML geometries; no DFT optimization',
        input=str(args.input),
        input_sha256=hashlib.sha256(source).hexdigest(),
        runner='scripts/dft_molecular_singlepoints.py',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        status='running',
        frames=[],
    )
    saved = []
    rows = []
    for i, a in enumerate(images):
        if a.pbc.any():
            raise ValueError('Molecular DFT only; periodic input rejected')
        mol = gto.M(
            atom=list(zip(a.get_chemical_symbols(), a.positions)),
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
        energy = float(mf.kernel() * Hartree)
        if not mf.converged:
            raise RuntimeError(f'SCF failed for image {i}')
        forces = -mf.nuc_grad_method().kernel() * Hartree / Bohr
        a.calc = SinglePointCalculator(a, energy=energy, forces=forces)
        saved.append(a)
        write(args.output / 'images.extxyz', saved)
        rows.append(
            dict(
                image=i,
                energy_eV=energy,
                relative_energy_eV=energy - rows[0]['energy_eV'] if rows else 0.0,
                max_raw_force_eV_A=float(np.linalg.norm(forces, axis=1).max()),
                SCF_converged=True,
            )
        )
        with (args.output / 'energies.csv').open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=rows[0])
            w.writeheader()
            w.writerows(rows)
        meta['frames'] = rows
        meta['status'] = 'complete' if len(saved) == len(images) else 'running'
        (args.output / 'summary.json').write_text(json.dumps(meta, indent=2) + '\n')
        print(i, energy, flush=True)


if __name__ == '__main__':
    main()
