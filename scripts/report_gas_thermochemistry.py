from pathlib import Path
import json, csv, hashlib, sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from plasma_surface.gas_thermochemistry import gas_thermo


def main():
    library = ROOT / 'data/literature/nist_gas_shomate.json'
    data = json.loads(library.read_text())
    rows = []
    for species, r in data['species'].items():
        for T in (300.0, 400.0, 450.0, 500.0, 600.0):
            if not r['temperature_range_K'][0] <= T <= r['temperature_range_K'][1]:
                continue
            for pressure in (1.0, 133.32236842105263, 100000.0):
                rows.append(gas_thermo(species, T, pressure, data))
    out = ROOT / 'data/kinetic_audit'
    with (out / 'gas_chemical_potentials.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)
    (out / 'gas_thermochemistry_manifest.json').write_text(
        json.dumps(
            dict(
                rows=len(rows),
                source_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
                code_sha256=hashlib.sha256(
                    (ROOT / 'plasma_surface/gas_thermochemistry.py').read_bytes()
                ).hexdigest(),
                scope='Ideal-gas chemical-potential terms on source reference; no surface free energies or barriers inferred',
                water_below_500K='Excluded: not in retrieved Shomate fit range',
            ),
            indent=2,
        )
        + '\n'
    )
    print('Evaluated', len(rows), 'source-bounded gas states')


if __name__ == '__main__':
    main()
