"""Diagnose stalled multilayer recession without changing any event rates."""
from pathlib import Path
from collections import Counter
import gzip,hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plasma_surface.multilayer import simulate,candidates,accessibility,neighbors
from plasma_surface.rate_evidence import environment
ROOT=Path(__file__).resolve().parents[1]

def main():
    src=ROOT/'data/multilayer/demo_trajectory.json';d=json.loads(src.read_text());p=d['parameters']
    r=simulate(d['initial'],np.arange(0,36.001,.25),seed=d['seed'],temperature=p['temperature_K'],arrival=p['arrival_s'],desorption=p['desorption_s'],penetration_A=p['penetration_A'],attenuation_A=p['attenuation_A'],policy='demonstration',dose_s=p['dose_s'],purge_s=p['purge_s'])
    traj=ROOT/'data/multilayer/extended_12cycle_trajectory.json.gz';traj.write_bytes(gzip.compress(json.dumps(r,separators=(',',':')).encode(),mtime=0))
    cycles=[]
    for c in range(12):
        ee=[e for e in r['events'] if 3*c<=e['time_s']<3*(c+1)];removed=[e for e in ee if e['kind'].endswith('release')]
        cycles.append(dict(cycle=c+1,events=len(ee),removed_by_band=dict(Counter(r['initial']['nodes'][e['site']]['layer'] for e in removed)),removed_by_element=dict(Counter(r['initial']['nodes'][e['site']]['element'] for e in removed)),newly_exposed_ids=sorted({i for e in ee for i in e['newly_exposed']})))
    g=r['final'];depth,ex,re=accessibility(g,p['penetration_A']);nb=neighbors(g);blocked=[]
    # Probe each potential adsorbed-HF state separately; do not modify the run.
    for i,n in enumerate(g['nodes']):
        if not re[i] or n['element']!='Si':continue
        occ=np.zeros(len(g['nodes']),bool);occ[i]=True
        events=[e for e in candidates(g,occ,policy='demonstration') if e['site']==i and e['kind']=='SiN_cleavage']
        for e in events:
            j=e['partner'];blocked.append(dict(site=i,partner=j,initial_band=n['layer'],remaining_substrate_bonds=len(nb[i]),F_terminations=n['termination'].count('F'),partner_H=g['nodes'][j]['termination'].count('H'),depth_A=float(depth[i]),rate_s=e['rate_s'],rate_status=e['rate_status'],environment=environment(g,occ,e,p['temperature_K'])))
    sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
    report=dict(original_trajectory_sha256=sha(src),extended_trajectory=str(traj.relative_to(ROOT)).replace('\\','/'),extended_trajectory_sha256=sha(traj),parameters=r['parameters'],policy=r['policy'],cycles=cycles,blocked_final_cleavages=blocked,probe='One hypothetical adsorbed HF at a time on the unchanged final graph; these are candidate probes, not executed events',conclusion='No further removal after cycle 2 under unchanged assumptions; missing final-cleavage rates prevent depth progression. Not evidence of physical self-limiting ALE.',runner='scripts/analyze_multilayer_stall.py',runner_sha256=sha(Path(__file__)),engine_sha256=sha(ROOT/'plasma_surface/multilayer.py'))
    out=ROOT/'docs/multilayer_results';(out/'cycle_diagnosis.json').write_text(json.dumps(report,indent=2)+'\n')
    counts=np.array([[x['removed_by_band'].get(l,0) for x in cycles] for l in range(5,-1,-1)])
    fig,ax=plt.subplots(1,2,figsize=(12,4.4),gridspec_kw={'width_ratios':[2.2,1]},layout='constrained')
    ax[0].imshow(counts,cmap='YlGnBu',vmin=0,vmax=21,aspect='auto')
    for y in range(6):
        for x in range(12):ax[0].text(x,y,str(counts[y,x]),ha='center',va='center',color='white' if counts[y,x]>12 else '#263B50')
    ax[0].set(xticks=range(12),xticklabels=range(1,13),yticks=range(6),yticklabels=['B5 top','B4','B3','B2','B1','B0 fixed'],xlabel='Dose/purge cycle (same rates and access depth)',title='Actual substrate removals per cycle and initial band')
    totals=Counter('bare N' if x['partner_H']==0 else 'NH2' if x['partner_H']==2 else 'NH' for x in blocked)
    ax[1].bar(totals.keys(),totals.values(),color=['#CA7151','#7869A6'][:len(totals)])
    for i,v in enumerate(totals.values()):ax[1].text(i,v+.2,str(v),ha='center')
    ax[1].set(ylim=(0,16),ylabel='Blocked reachable SiF3 sites',xlabel='Remaining Si-N backbond partner',title='Why deeper removal stops')
    fig.suptitle('12-cycle diagnostic: a missing-rate plateau, not demonstrated ALE saturation',fontsize=14)
    fig.savefig(out/'cycle_depth_diagnosis.png',dpi=160);plt.close(fig)
    print(report['conclusion']);print('Blocked environments:',dict(totals))
if __name__=='__main__':main()
