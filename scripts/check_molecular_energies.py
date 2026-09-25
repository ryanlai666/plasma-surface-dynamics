"""Numerically compare energy-only evaluation with ASE's energy/force evaluator."""

from pathlib import Path
import hashlib, json
import numpy as np
import torch
from threadpoolctl import threadpool_limits
from ase.io import read
from mace.calculators import MACECalculator
from screen_molecular_sites_energy import energy_snapshot
from molecular_surface_campaign import ROOT, independent


def main():
    torch.set_num_threads(1)
    threadpool_limits(limits=1)
    calc = MACECalculator(
        model_paths=str(ROOT / '.cache/mace/mp_0b2_small.model'),
        device='cpu',
        default_dtype='float64',
    )
    rows = []
    for surface in ['Si100', 'Si111', 'beta_Si3N4_001', 'alpha_quartz_001']:
        for species in ['HF', 'CH3F', 'SiCl4']:
            source = ROOT / 'data/surface_paths' / surface / species / 'mace_approach/images.extxyz'
            for index in [3, 8]:
                a = read(source, index)
                stored = a.get_potential_energy()
                a.calc = independent(calc)
                fast = energy_snapshot(a).get_potential_energy()
                reference = a.get_potential_energy()
                forces = a.get_forces(apply_constraint=False)
                assert (
                    abs(fast - reference) < 1e-7
                    and abs(stored - reference) < 1e-6
                    and np.isfinite(forces).all()
                )
                rows.append(
                    dict(
                        surface=surface,
                        species=species,
                        image=index,
                        energy_only_eV=fast,
                        regular_eV=reference,
                        stored_eV=stored,
                    )
                )
    result = dict(
        status='passed',
        geometries=len(rows),
        max_energy_only_error_eV=max(abs(r['energy_only_eV'] - r['regular_eV']) for r in rows),
        max_reproduction_error_eV=max(abs(r['stored_eV'] - r['regular_eV']) for r in rows),
        runner='scripts/check_molecular_energies.py',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        rows=rows,
    )
    (ROOT / 'docs/dry_etch_results/molecular_energy_checks.json').write_text(
        json.dumps(result, indent=2) + '\n'
    )
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}))


if __name__ == '__main__':
    main()
