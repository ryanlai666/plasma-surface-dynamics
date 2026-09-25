"""Fetch and extract published SiCl4/Si(100) stationary-point coordinates.

Requires ASE and PyMuPDF. Outputs are attributed source geometries, not new DFT.
"""

from pathlib import Path
import argparse
import hashlib
import json
import re
import urllib.request
import zipfile
import numpy as np
from ase import Atoms
from ase.io import write
import pymupdf

ROOT = Path(__file__).resolve().parents[1]
URL = 'https://mdpi-res.com/d_attachment/symmetry/symmetry-15-00213/article_deploy/symmetry-15-00213-s001.zip'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--archive', type=Path, default=ROOT / 'data/raw/sicl4_surface/supplement.zip'
    )
    args = parser.parse_args()
    if not args.archive.exists():
        args.archive.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(URL, timeout=90) as r:
            raw = r.read()
        args.archive.write_bytes(raw)
    raw = args.archive.read_bytes()
    target = ROOT / 'data/surface_paths/Si100_c4x2/SiCl4/published'
    expected = json.loads((target / 'provenance.json').read_text())
    if hashlib.sha256(raw).hexdigest() != expected['source_archive_sha256']:
        raise ValueError('Source archive checksum changed; review version before extraction')
    with zipfile.ZipFile(args.archive) as z:
        pdf = z.read('Supporting Information.pdf')
    doc = pymupdf.open(stream=pdf, filetype='pdf')
    text = '\n'.join(p.get_text() for p in doc)
    parts = re.split(r'^\d+\. ([^\n]+)\n', text, flags=re.M)
    inventory = []
    for i in range(1, len(parts), 2):
        name = parts[i].strip()
        key = name.replace('Si(100) surface.', 'clean_surface')
        rows = re.findall(r'\b(Si|Cl|H)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)', parts[i + 1])
        symbols = [r[0] for r in rows]
        xyz = np.array([[float(v) for v in r[1:]] for r in rows])
        count = 112 if key == 'clean_surface' else 5 if key == 'STC' else 117
        if len(rows) != count or not np.isfinite(xyz).all():
            raise ValueError('Incomplete coordinates: ' + name)
        a = Atoms(symbols, positions=xyz)
        a.info.update(
            source_doi='10.3390/sym15010213',
            source_label=name,
            coordinates_only=True,
            cell_not_supplied=True,
        )
        path = target / (key + '.extxyz')
        write(path, a)
        inventory.append(
            dict(
                label=key,
                natoms=len(a),
                formula=a.get_chemical_formula(),
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            )
        )
    if len(inventory) != 17:
        raise ValueError('Expected 17 source structures')
    (target / 'inventory.json').write_text(json.dumps(inventory, indent=2) + '\n')
    print('Verified source checksum and extracted 17 stationary structures')


if __name__ == '__main__':
    main()
