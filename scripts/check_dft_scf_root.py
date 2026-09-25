"""Check alternative SCF starts and internal stability for a problematic image."""

from pathlib import Path
import json
import numpy as np
from pyscf import gto, dft, lib
from ase.io import read, write
from ase.units import Hartree, Bohr
from ase.calculators.singlepoint import SinglePointCalculator

ROOT = Path(__file__).resolve().parents[1]
out = ROOT / 'data/surface_paths/SiO_capped_motif/HF/dft_path'
a = read(out / 'input.extxyz', 5)
lib.num_threads(2)
mol = gto.M(
    atom=list(zip(a.get_chemical_symbols(), a.positions)),
    basis='def2-svp',
    charge=0,
    spin=0,
    unit='Angstrom',
    verbose=0,
)
results = []
for guess in ['atom', 'minao', 'huckel']:
    mf = dft.RKS(mol).density_fit()
    mf.xc = 'PBE'
    mf.grids.level = 3
    mf.conv_tol = 1e-9
    mf.max_cycle = 400
    mf.level_shift = 0.25
    mf.damp = 0.2
    mf.init_guess = guess
    mf.kernel()
    mf.level_shift = 0.0
    mf.damp = 0.0
    mf.kernel(dm0=mf.make_rdm1())
    stable = False
    if mf.converged:
        for k in range(3):
            orbitals, _, stable, _ = mf.stability(internal=True, external=False, return_status=True)
            if stable:
                break
            mf.kernel(dm0=mf.make_rdm1(orbitals, mf.mo_occ))
    record = dict(
        guess=guess,
        converged=bool(mf.converged),
        internally_stable=bool(stable),
        energy_eV=float(mf.e_tot * Hartree),
        electron_count=float(np.trace(mf.make_rdm1() @ mf.get_ovlp())),
    )
    if record['converged'] and record['internally_stable']:
        b = a.copy()
        b.calc = SinglePointCalculator(
            b, energy=record['energy_eV'], forces=-mf.nuc_grad_method().kernel() * Hartree / Bohr
        )
        write(out / ('image5_' + guess + '.extxyz'), b)
    results.append(record)
    (out / 'SCF_image5_audit.json').write_text(json.dumps(results, indent=2) + '\n')
    print(record, flush=True)
