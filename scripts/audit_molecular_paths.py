"""Reject geometrically pathological paths before interpreting numerical peaks."""

from pathlib import Path
import hashlib, json
import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[1]
for parent in (ROOT / 'data/surface_paths').glob('*_capped_motif/HF'):
    folder = parent / 'omol25_refined'
    meta = json.loads((folder / 'summary.json').read_text())
    images = read(folder / 'images.extxyz', ':')
    rows = []
    for i, a in enumerate(images):
        distances = a.get_all_distances()
        distances[np.diag_indices(len(a))] = np.inf
        ij = np.unravel_index(distances.argmin(), distances.shape)
        rows.append(
            dict(
                image=i,
                minimum_pair_distance_A=float(distances[ij]),
                pair_indices=[int(x) for x in ij],
                pair_species=[a[int(x)].symbol for x in ij],
                severe_overlap=bool(distances[ij] < 0.5),
            )
        )
    overlap = any(r['severe_overlap'] for r in rows)
    result = dict(
        status='rejected_as_reaction_path',
        reason='Severe atomic overlap and unconverged NEB'
        if overlap
        else 'NEB did not converge; no verified transition state',
        severe_overlap=overlap,
        rows=rows,
        ml_status=meta['status'],
        path_sha256=hashlib.sha256((folder / 'images.extxyz').read_bytes()).hexdigest(),
        runner='scripts/audit_molecular_paths.py',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        DFT_interpretation='Fixed-geometry diagnostic energies only; exclude this path peak from kinetic parameterization',
    )
    if (parent / 'dft_path/SCF_image5_audit.json').exists():
        result['rejected_DFT_images'] = [5]
        result['electronic_note'] = (
            'Anomalous Newton SCF solution; three independent alternative starts failed. Do not interpret its energy as a reaction-path point.'
        )
    (parent / 'path_quality.json').write_text(json.dumps(result, indent=2) + '\n')
    print(parent.parent.name, result['reason'], min(r['minimum_pair_distance_A'] for r in rows))
