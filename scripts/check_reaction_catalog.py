"""Check elemental and charge accounting for proposed reactions; not kinetics."""
import json
from pathlib import Path


def validate(catalog):
    species=catalog['species']
    ids=set()
    for reaction in catalog['reactions']:
        if reaction['id'] in ids:raise ValueError('Duplicate reaction ID')
        ids.add(reaction['id'])
        totals=[]
        for side in ['reactants','products']:
            tally={'charge':0}
            for name,n in reaction[side].items():
                if not isinstance(n,int) or n<=0:raise ValueError('Stoichiometry must be a positive integer')
                s=species[name];tally['charge']+=n*s['charge']
                for element,count in s['composition'].items():tally[element]=tally.get(element,0)+n*count
            totals.append(tally)
        if totals[0]!=totals[1]:raise ValueError('Unbalanced reaction: '+reaction['id'])
    return len(ids)


if __name__=='__main__':
    root=Path(__file__).resolve().parents[1]
    c=json.loads((root/'configs/reaction_candidates.json').read_text())
    print(f'{validate(c)} candidate reactions pass elemental and charge balance; no kinetics validated.')
