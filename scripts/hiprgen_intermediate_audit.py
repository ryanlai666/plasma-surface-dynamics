"""Run a limited HiPRGen bucketing/decision-tree pilot on explicit capped motifs.

The full MPI/OpenBabel/thermochemical pipeline is not run. No rates are assigned.
"""

from pathlib import Path
import argparse, ast, csv, hashlib, itertools, json, sqlite3
from enum import Enum
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
UPSTREAM = 'a0dddfedc21be0121745e5f33f27ad8aafe796ea'


class Terminal(Enum):
    KEEP = 1
    DISCARD = -1


def upstream_function(file, name, namespace):
    tree = ast.parse(file.read_text(encoding='utf8'))
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
    # Execute the complete unmodified upstream function, avoiding unrelated optional imports.
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(file), 'exec'), namespace)
    return namespace[name]


def library(ligand, halogen):
    entries = []
    for cap_h in range(5):
        for bound in range(5 - cap_h):
            count = 4 - cap_h - bound
            composition = {'Si': 1}
            if cap_h + (2 if ligand == 'NH2' else 1) * bound:
                composition['H'] = cap_h + (2 if ligand == 'NH2' else 1) * bound
            if bound:
                composition['N' if ligand == 'NH2' else 'O'] = bound
            if count:
                composition[halogen] = count
            name = (
                'Si'
                + (f'H{cap_h}' if cap_h else '')
                + (f'{halogen}{count}' if count else '')
                + (f'({ligand}){bound}' if bound else '')
            )
            entries.append(
                dict(
                    id=name,
                    composition=composition,
                    charge=0,
                    multiplicity=1,
                    phase='capped_molecular_proxy',
                    cap_H=cap_h,
                    bound_ligands=bound,
                    halogens=count,
                    evidence='candidate graph; energy and surface stability unverified',
                    fixed_cap_hydrogens=cap_h,
                )
            )
    for name, comp in [
        ('H' + halogen, {'H': 1, halogen: 1}),
        (
            'NH3' if ligand == 'NH2' else 'H2O',
            {'N': 1, 'H': 3} if ligand == 'NH2' else {'O': 1, 'H': 2},
        ),
    ]:
        entries.append(
            dict(
                id=name,
                composition=comp,
                charge=0,
                multiplicity=1,
                phase='gas',
                cap_H=None,
                bound_ligands=None,
                halogens=None,
                evidence='known molecular identity; no species free energy assigned',
            )
        )
    return entries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'data/reaction_network/hiprgen')
    args = parser.parse_args()
    snapshot = ROOT / 'third_party/hiprgen_snapshot'
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    bucket = upstream_function(
        snapshot / 'bucketing.py',
        'bucket',
        dict(
            sqlite3=sqlite3, combinations_with_replacement=itertools.combinations_with_replacement
        ),
    )
    decide = upstream_function(
        snapshot / 'reaction_questions.py', 'run_decision_tree', dict(Terminal=Terminal)
    )
    all_reactions = []
    all_species = []
    audit = []
    for ligand, halogen in itertools.product(['NH2', 'OH'], ['F', 'Cl']):
        family = ligand + '_' + halogen
        entries = library(ligand, halogen)
        mols = [
            SimpleNamespace(
                ind=i,
                species=[e for e, n in d['composition'].items() for _ in range(n)],
                metadata=d,
            )
            for i, d in enumerate(entries)
        ]
        db = out / (family + '_buckets.sqlite')
        bucket(mols, str(db))
        con = sqlite3.connect(db)
        groups = {}
        for a, b, c, g in con.execute('SELECT * FROM complexes'):
            groups.setdefault(c, []).append(tuple(i for i in (a, b) if i >= 0))
        con.close()
        stats = {}
        generated = 0
        kept = 0

        def shared(r, m, p):
            return bool(set(r['left']) & set(r['right']))

        def not_single_si(r, m, p):
            return any(
                sum(entries[i]['composition'].get('Si', 0) for i in r[k]) != 1
                for k in ['left', 'right']
            )

        def changing_cap(r, m, p):
            return next(
                entries[i]['cap_H'] for i in r['left'] if 'Si' in entries[i]['composition']
            ) != next(entries[i]['cap_H'] for i in r['right'] if 'Si' in entries[i]['composition'])

        def not_one_ligand(r, m, p):
            return (
                abs(
                    next(
                        entries[i]['bound_ligands']
                        for i in r['left']
                        if 'Si' in entries[i]['composition']
                    )
                    - next(
                        entries[i]['bound_ligands']
                        for i in r['right']
                        if 'Si' in entries[i]['composition']
                    )
                )
                != 1
            )

        def always(r, m, p):
            return True

        tree = [
            (shared, Terminal.DISCARD),
            (not_single_si, Terminal.DISCARD),
            (changing_cap, Terminal.DISCARD),
            (not_one_ligand, Terminal.DISCARD),
            (always, Terminal.KEEP),
        ]
        for sides in groups.values():
            for left, right in itertools.permutations(sides, 2):
                r = dict(
                    left=left,
                    right=right,
                    reactants=left + (-1,) * (2 - len(left)),
                    products=right + (-1,) * (2 - len(right)),
                    number_of_reactants=len(left),
                    number_of_products=len(right),
                )
                trace = []
                generated += 1
                accepted = decide(r, mols, {}, tree, trace)
                reason = trace[0].__name__
                stats[reason] = stats.get(reason, 0) + 1
                if not accepted:
                    continue
                kept += 1
                reactants = [entries[i]['id'] for i in left]
                products = [entries[i]['id'] for i in right]
                center = next(entries[i] for i in left if 'Si' in entries[i]['composition'])
                after = next(entries[i] for i in right if 'Si' in entries[i]['composition'])
                direction = (
                    'forward_substitution'
                    if after['halogens'] > center['halogens']
                    else 'reverse_substitution'
                )
                all_reactions.append(
                    dict(
                        id=family + '_' + str(kept),
                        family=family,
                        reactants=reactants,
                        products=products,
                        direction=direction,
                        cap_H=center['cap_H'],
                        initial_bound_ligands=center['bound_ligands'],
                        final_bound_ligands=after['bound_ligands'],
                        barrier_eV=None,
                        rate_s=None,
                        status='composition_and_fixed_cap_filter_passed_not_TS_verified',
                        needed=[
                            '3D structure and atom mapping',
                            'surface embedding and termination',
                            'minimum and TS searches',
                            'free energies/prefactors and gas chemical potential',
                        ],
                    )
                )
        all_species.extend(dict(family=family, **s) for s in entries)
        audit.append(
            dict(
                family=family,
                species=len(entries),
                composition_buckets=len(groups),
                directed_pairs=generated,
                retained=kept,
                filter_counts=stats,
            )
        )
    result = dict(
        upstream_repository='https://github.com/BlauGroup/HiPRGen',
        upstream_commit=UPSTREAM,
        execution_scope='Unmodified bucket() and run_decision_tree() functions; custom fixed-cap local-motif filters. Full MPI, molecule filtering, graph isomorphism and thermochemical pipeline NOT executed.',
        source_sha256={
            p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in snapshot.glob('*.py')
        },
        runner='scripts/hiprgen_intermediate_audit.py',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        families=audit,
        species=all_species,
        reactions=all_reactions,
        completeness='Finite supplied species library only; cannot prove mechanism completeness or kinetic accessibility',
    )
    (out / 'network.json').write_text(json.dumps(result, indent=2) + '\n')
    print(
        json.dumps(dict(families=audit, retained_directed_reactions=len(all_reactions))), flush=True
    )


if __name__ == '__main__':
    main()
