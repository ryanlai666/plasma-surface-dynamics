"""Assign the connected saddle by actual bond distances, not its starting guess."""

from pathlib import Path
import hashlib, json
from ase.io import read

ROOT = Path(__file__).resolve().parents[1]
p = ROOT / 'data/surface_paths/beta_Si3N4_001/H2O/mace_local_saddle'
d = json.loads((p / 'summary.json').read_text())
images = read(p / 'images.extxyz', ':')
rows = []
for label, i in [('IS', 0), ('TS', d['peak_image']), ('FS', len(images) - 1)]:
    a = images[i]
    rows.append(
        dict(
            state=label,
            image=i,
            OH_bond_A=float(a.get_distance(28, 29, mic=True)),
            separated_O_H_A=float(a.get_distance(28, 30, mic=True)),
            H_N25_A=float(a.get_distance(30, 25, mic=True)),
            H_N26_A=float(a.get_distance(30, 26, mic=True)),
            O_Si15_A=float(a.get_distance(28, 15, mic=True)),
        )
    )
assert rows[0]['H_N25_A'] < 1.2 and rows[-1]['H_N26_A'] < 1.2
assert all(r['OH_bond_A'] < 1.2 and r['separated_O_H_A'] > 2 for r in rows)
audit = dict(
    event='H transfer N25-H + N26 -> N25 + N26-H; adsorbed Si15-OH spectator',
    origin='Water-derived OH + H surface configuration; NOT intact H2O dissociation',
    rows=rows,
    runner='scripts/audit_surface_saddle_identity.py',
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    path_sha256=hashlib.sha256((p / 'images.extxyz').read_bytes()).hexdigest(),
)
(p / 'identity_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
d['event_identity'] = audit['event']
d['original_chemical_assignment_rejected'] = (
    'Alternative intact-water IS: false; it is already OH + surface H.'
)
d['identity_audit'] = 'identity_audit.json'
d['assumptions'][0] = (
    'IS and FS are distinct H adsorption sites beside Si-OH; both are already dissociated water-derived configurations.'
)
(p / 'summary.json').write_text(json.dumps(d, indent=2) + '\n')
print(json.dumps(audit), flush=True)
