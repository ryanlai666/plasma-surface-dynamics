"""Run multilayer prototype, accessibility controls and evidence-gated control."""
from pathlib import Path
from collections import Counter
import json,hashlib,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from plasma_surface.multilayer import simulate

def main():
    source=ROOT/'data/multilayer/sin_graph.json';initial=json.loads(source.read_text());out=ROOT/'docs/multilayer_results';out.mkdir(exist_ok=True)
    r=simulate(initial,np.linspace(0,9,37),seed=71,policy='demonstration');(ROOT/'data/multilayer/demo_trajectory.json').write_text(json.dumps(r,separators=(',',':'))+'\n');runs=[]
    for depth in [0.,1.,2.,4.]:
        for seed in [71,72,73]:
            a=simulate(initial,[0,3,6,9],seed=seed,penetration_A=depth,policy='demonstration')
            runs.append(dict(penetration_A=depth,seed=seed,removed=dict(Counter(n['element'] for n in a['final']['nodes'] if not n['active'])),removed_by_initial_layer=dict(Counter(n['layer'] for n in a['final']['nodes'] if not n['active'])),events=len(a['events']),gas=a['snapshots'][-1]['gas'],conservation_exact=a['conservation_exact']))
    strict=simulate(initial,[0,9],policy='validated_only');assert not strict['events']
    d=dict(status='multilayer_connectivity_demonstration_not_validated_etch_prediction',substrate_atoms=len(initial['nodes']),initial_bonds=len(initial['bonds']),initial_layers=sorted(set(n['layer'] for n in initial['nodes'])),fixed_atoms=sum(n['fixed'] for n in initial['nodes']),baseline_events=len(r['events']),event_counts=dict(Counter(e['kind'] for e in r['events'])),removed=dict(Counter(n['element'] for n in r['final']['nodes'] if not n['active'])),removed_by_initial_layer=dict(Counter(n['layer'] for n in r['final']['nodes'] if not n['active'])),initial_layer_sizes=dict(Counter(n['layer'] for n in initial['nodes'])),newly_exposed_atoms=len({i for e in r['events'] for i in e['newly_exposed']}),initial_inventory=r['initial_inventory'],final_gas=r['snapshots'][-1]['gas'],conservation_exact=True,control_runs=runs,validated_only_events=len(strict['events']),parameters=r['parameters'],geometry_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),trajectory_sha256=hashlib.sha256((ROOT/'data/multilayer/demo_trajectory.json').read_bytes()).hexdigest(),code_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['plasma_surface/multilayer.py','plasma_surface/rate_evidence.py','scripts/run_multilayer.py']})
    (out/'summary.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d),flush=True)
if __name__=='__main__':main()
