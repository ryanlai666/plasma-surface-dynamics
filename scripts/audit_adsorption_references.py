"""Decompose relaxed adsorption energies to detect metastable clean-slab references."""

from pathlib import Path
import json, hashlib
import torch
from threadpoolctl import threadpool_limits
from ase.io import read, write
from mace.calculators import MACECalculator
from molecular_surface_campaign import ROOT, independent, snapshot


def main():
    torch.set_num_threads(2)
    threadpool_limits(limits=1)
    model = ROOT / '.cache/mace/mp_0b2_small.model'
    calc = MACECalculator(model_paths=str(model), device='cpu', default_dtype='float64')
    rows = []
    for p in (ROOT / 'data/surface_paths').glob('*/HF/relaxed_adsorption/summary.json'):
        d = json.loads(p.read_text())
        for r in d['starts']:
            folder = ROOT / Path(r['folder'].replace('\\', '/'))
            a = read(folder / 'trajectory.extxyz', -1)
            slab = a[:-2]
            mol = a[-2:]
            slab.calc = independent(calc)
            mol.calc = independent(calc)
            es = slab.get_potential_energy()
            em = mol.get_potential_energy()
            ec = r['energy_eV']
            deform = es - d['bare_slab']['energy_eV']
            molshift = em - d['gas_energy_eV']
            interaction = ec - es - em
            record = dict(
                surface=r['surface'],
                start=r['start'],
                combined_converged=r['converged'],
                bare_reference_converged=d['bare_slab']['converged'],
                frozen_fragment_interaction_eV=float(interaction),
                substrate_reference_shift_eV=float(deform),
                molecular_layer_reference_shift_eV=float(molshift),
                apparent_adsorption_difference_eV=float(
                    ec - d['bare_slab']['energy_eV'] - d['gas_energy_eV']
                ),
                closure_error_eV=float(
                    interaction
                    + deform
                    + molshift
                    - (ec - d['bare_slab']['energy_eV'] - d['gas_energy_eV'])
                ),
                lower_clean_slab_geometry_found=bool(deform < -0.1),
                screening_threshold_eV=-0.1,
                interpretation='Negative substrate shift below -0.1 eV flags a lower-energy clean-slab geometry than the relaxed reference; that threshold is diagnostic, not a convergence standard.',
                trajectory_sha256=hashlib.sha256(
                    (folder / 'trajectory.extxyz').read_bytes()
                ).hexdigest(),
                runner='scripts/audit_adsorption_references.py',
                runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            )
            write(folder / 'frozen_slab.extxyz', snapshot(slab))
            write(folder / 'frozen_molecular_layer.extxyz', snapshot(mol))
            (folder / 'reference_audit.json').write_text(json.dumps(record, indent=2) + '\n')
            rows.append(record)
    (ROOT / 'docs/dry_etch_results/adsorption_reference_audit.json').write_text(
        json.dumps(rows, indent=2) + '\n'
    )
    print(json.dumps(rows), flush=True)


if __name__ == '__main__':
    main()
