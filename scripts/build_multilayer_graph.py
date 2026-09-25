"""Construct a terminated Si3N4 bond graph from the public crystalline CIF."""

from pathlib import Path
import json, hashlib
import numpy as np
from ase.build import surface
from ase.io import read
from ase.neighborlist import neighbor_list

ROOT = Path(__file__).resolve().parents[1]


def main():
    source = ROOT / 'docs/mace_results/2102550.cif'
    bulk = read(source)
    atoms = surface(bulk, (0, 0, 1), 6, vacuum=6).repeat((2, 2, 1))
    atoms.pbc = [True, True, False]
    ii, jj = neighbor_list('ij', atoms, 2.05)
    bonds = sorted(
        {
            tuple(sorted([int(i), int(j)]))
            for i, j in zip(ii, jj)
            if atoms[i].symbol != atoms[j].symbol
        }
    )
    degree = np.zeros(len(atoms), int)
    for i, j in bonds:
        degree[i] += 1
        degree[j] += 1
    frac = atoms.get_scaled_positions(wrap=True)
    rng = np.random.default_rng(71)
    zmin = atoms.positions[:, 2].min()
    nodes = []
    for i, a in enumerate(atoms):
        target = 4 if a.symbol == 'Si' else 3
        missing = target - int(degree[i])
        assert missing >= 0
        termination = (
            ['H'] * missing
            if a.symbol == 'N'
            else rng.choice(['H', 'F', 'Cl'], size=missing, p=[0.2, 0.5, 0.3]).tolist()
        )
        nodes.append(
            dict(
                id=i,
                element=a.symbol,
                position_A=a.position.tolist(),
                column=int(frac[i, 0] * 6) + 6 * int(frac[i, 1] * 6),
                layer=int((a.position[2] - zmin) / bulk.cell[2, 2] + 1e-6),
                active=True,
                fixed=bool(a.position[2] < zmin + bulk.cell[2, 2] - 0.01),
                termination=termination,
            )
        )
    out = ROOT / 'data/multilayer'
    out.mkdir(exist_ok=True)
    d = dict(
        nodes=nodes,
        bonds=[list(x) for x in bonds],
        cell_A=atoms.cell.tolist(),
        pbc=[True, True, False],
        source=str(source.relative_to(ROOT)).replace('\\', '/'),
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        source_url='https://www.crystallography.net/cod/2102550.html',
        construction='beta-Si3N4 (001), six unit-cell layers, 2x2 lateral repeat; Si-N cutoff 2.05 A; fixed bottom unit-cell layer',
        accessibility='6x6 fractional-xy columns; current top atom exposed, finite penetration below moving local envelope',
        termination='Dangling N valences capped by H; dangling Si valences seeded H/F/Cl with probabilities 0.2/0.5/0.3; seed 71; illustrative unrelaxed chemistry',
        runner='scripts/build_multilayer_graph.py',
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    (out / 'sin_graph.json').write_text(json.dumps(d, indent=2) + '\n')
    print(len(nodes), 'substrate atoms', len(bonds), 'bonds')


if __name__ == '__main__':
    main()
