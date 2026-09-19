"""Inventory encountered AND unexecuted candidate environments for targeted rate data."""
from pathlib import Path
from collections import Counter
from copy import deepcopy
import hashlib,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from plasma_surface.multilayer import candidates,apply
from plasma_surface.rate_evidence import environment

def main():
    trajectory=ROOT/'data/multilayer/demo_trajectory.json';d=json.loads(trajectory.read_text());records={}
    def add(graph,prec,event,executed):
        env=environment(graph,prec,event,d['parameters']['temperature_K']);key=env['key']
        if key not in records:
            ids=env['atom_ids'];records[key]=dict(**env,event_kind=event['kind'],reference_pathway=event['reference_pathway'],demonstration_barrier_eV=event['barrier_eV'],rate_status='missing_qualified_rate',observed_events=0,candidate_snapshots=0,requirements=['Build and relax matching terminated IS and FS, with HF precursor where appropriate','Verify charge/spin, boundary caps, reference energy, coverage and neighbor connectivity','Obtain connected TS, one imaginary mode, prefactor/free-energy and uncertainty','For adsorption: measured/calculated sticking, gas flux and site-area mapping instead of an arbitrary barrier'],example=dict(nodes=[deepcopy(graph['nodes'][i]) for i in ids],bonds=[b for b in graph['bonds'] if all(i in ids for i in b)],precursor_sites=[i for i in ids if prec[i]],event=event),structure_status='Connectivity request only; ligand/HF positions and relaxed electronic structure are not supplied')
        if 'FS_connectivity' not in records[key]['example']:
            after=deepcopy(graph);pafter=np.array(prec,copy=True);products={};product_formulas={}
            apply(after,pafter,event,products,product_formulas)
            records[key]['example'].update(cell_A=graph['cell_A'],pbc=graph['pbc'],FS_connectivity=dict(nodes=[after['nodes'][i] for i in env['atom_ids']],bonds=[b for b in after['bonds'] if all(i in env['atom_ids'] for i in b)],precursor_sites=[i for i in env['atom_ids'] if pafter[i]]),gas_delta=products,gas_formulas=product_formulas)
        records[key]['observed_events' if executed else 'candidate_snapshots']+=1
    graph=deepcopy(d['initial']);prec=np.zeros(len(graph['nodes']),bool);gas={};formulas={}
    for event in d['events']:
        add(graph,prec,event,True);apply(graph,prec,event,gas,formulas)
    for frame in d['snapshots']:
        graph=deepcopy(d['initial']);graph['bonds']=frame['bonds']
        for i,n in enumerate(graph['nodes']):n['active']=frame['active'][i];n['termination']=[e for e,count in frame['terminations'][i].items() for _ in range(count)]
        prec=np.array(frame['precursors'])
        for e in candidates(graph,prec,**dict(temperature=d['parameters']['temperature_K'],arrival=5.,desorption=100.,penetration_A=d['parameters']['penetration_A'],attenuation_A=1.,policy='demonstration')):add(graph,prec,e,False)
    ordered=sorted(records.values(),key=lambda r:(r['example']['event']['rate_status']=='disabled_missing_final_cleavage_barrier',r['observed_events']+r['candidate_snapshots']),reverse=True)
    for rank,r in enumerate(ordered,1):r['priority_rank']=rank;r['priority_basis']='Missing final-cleavage barriers first, then observed count + sampled candidate presence; priority heuristic, not degree of rate control'
    out=ROOT/'data/multilayer/rate_requests';out.mkdir(exist_ok=True)
    # Retain full missing-environment inventory; expand the top 12 reproducible requests separately.
    for r in ordered[:12]:
        p=out/r['key'][:16];p.mkdir(exist_ok=True);(p/'request.json').write_text(json.dumps(r,indent=2)+'\n')
    (out/'index.json').write_text(json.dumps(dict(requests=ordered,trajectory_sha256=hashlib.sha256(trajectory.read_bytes()).hexdigest(),runner='scripts/queue_multilayer_rates.py',runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')
    print(len(ordered),'unique missing rate environments;',sum(r['observed_events'] for r in ordered),'events covered')
if __name__=='__main__':main()
