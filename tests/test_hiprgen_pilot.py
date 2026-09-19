"""Check bounded HiPRGen pilot output without optional upstream dependencies."""
from collections import Counter
from pathlib import Path
import hashlib,json,sqlite3
ROOT=Path(__file__).resolve().parents[1]
def test_native_bucket_output_and_candidate_conservation():
    p=ROOT/'data/reaction_network/hiprgen';d=json.loads((p/'network.json').read_text())
    assert len(d['reactions'])==80 and len(d['families'])==4
    assert 'Full MPI' in d['execution_scope'] and 'NOT executed' in d['execution_scope']
    for f in d['families']:
        entries={s['id']:s for s in d['species'] if s['family']==f['family']}
        with sqlite3.connect(p/(f['family']+'_buckets.sqlite')) as con:assert con.execute('SELECT count(*) FROM complexes').fetchone()[0]==17+17*18//2
        for r in (r for r in d['reactions'] if r['family']==f['family']):
            sides=[]
            for side in ['reactants','products']:
                total=Counter()
                for name in r[side]:total.update(entries[name]['composition'])
                sides.append(total)
            assert sides[0]==sides[1]
            assert r['rate_s'] is None and r['barrier_eV'] is None
            assert abs(r['initial_bound_ligands']-r['final_bound_ligands'])==1
    assert sum(r['direction']=='forward_substitution' for r in d['reactions'])==40
    assert d['script_sha256']==hashlib.sha256((ROOT/d['runner']).read_bytes()).hexdigest()
    for name,sha in d['source_sha256'].items():assert hashlib.sha256((ROOT/'third_party/hiprgen_snapshot'/name).read_bytes()).hexdigest()==sha
