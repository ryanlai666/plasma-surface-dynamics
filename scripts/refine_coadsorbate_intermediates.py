"""Refine coadsorbates and matched parent references with BFGS; retain initial attempts."""
from pathlib import Path
import csv,json,time,hashlib
import numpy as np
import torch
from threadpoolctl import threadpool_limits
from ase.io import read,write
from ase.optimize import BFGS
from mace.calculators import MACECalculator
from relax_coadsorbate_intermediates import ROOT, OUT, sha, adsorbate_curvature
from molecular_surface_campaign import independent,snapshot

def main():
    torch.set_num_threads(2);threadpool_limits(limits=1)
    model=ROOT/'.cache/mace/mp_0b2_small.model';calc=MACECalculator(model_paths=str(model),device='cpu',default_dtype='float64')
    original=OUT/'summary.json';d=json.loads(original.read_text());results=[];parents={}
    for r in d['results']:
        key=(r['surface'],r['start'])
        if key not in parents:
            parent_source=OUT/r['surface']/f"parent_{r['start']}"/'final.extxyz'
            parent_folder=parent_source.parent/'refinement';parent_folder.mkdir(exist_ok=True)
            parent=read(parent_source);parent.info=dict(surface=r['surface'],calculation_kind='Single-HF reference relaxation; not a reaction path')
            parent.calc=independent(calc);parent_frames=[];parent_opt=BFGS(parent,logfile=str(parent_folder/'optimization.log'),maxstep=.06)
            parent_opt.attach(lambda:parent_frames.append(snapshot(parent)),interval=1)
            parent_ok=bool(parent_opt.run(fmax=.02,steps=300));write(parent_folder/'trajectory.extxyz',parent_frames);write(parent_folder/'final.extxyz',parent_frames[-1])
            parents[key]=dict(converged=parent_ok,steps=parent_opt.nsteps,energy_eV=float(parent.get_potential_energy()),max_mobile_force_eV_A=float(np.linalg.norm(parent.get_forces(),axis=1).max()),folder=str(parent_folder.relative_to(ROOT)).replace('\\','/'),source_sha256=sha(parent_source),final_sha256=sha(parent_folder/'final.extxyz'))
            (parent_folder/'summary.json').write_text(json.dumps(parents[key],indent=2)+'\n')
            print('parent',key,parent_ok,parents[key]['max_mobile_force_eV_A'],flush=True)
        r=dict(r);r['reference']=parents[key]
        source=ROOT/r['folder'];target=source/'refinement';target.mkdir(exist_ok=True)
        a=read(source/'final.extxyz');a.info=dict(surface=r['surface'],coadsorbate=r['coadsorbate'],start=r['start'],calculation_kind='Local coadsorbate relaxation; not a reaction path',host_atoms=r['host_atoms'])
        a.calc=independent(calc);frames=[];opt=BFGS(a,logfile=str(target/'optimization.log'),maxstep=.06)
        start=time.perf_counter()
        # Refine every candidate to the same tighter tolerance, including failures.
        opt.attach(lambda:frames.append(snapshot(a)),interval=1)
        ok=bool(opt.run(fmax=.02,steps=300));write(target/'trajectory.extxyz',frames);write(target/'final.extxyz',frames[-1])
        energies=[float(x.get_potential_energy()) for x in frames]
        with (target/'energies.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['optimization_step','energy_eV','energy_relative_to_start_eV']);w.writerows((i,e,e-energies[0]) for i,e in enumerate(energies))
        n=r['host_atoms'];dist=a.get_all_distances(mic=True);result=dict(r)
        result.update(initial_attempt=r['folder'],initial_attempt_sha256=sha(source/'summary.json'),converged=ok,steps=opt.nsteps,frames=len(frames),energy_eV=energies[-1],max_mobile_force_eV_A=float(np.linalg.norm(a.get_forces(),axis=1).max()),elapsed_s=time.perf_counter()-start,folder=str(target.relative_to(ROOT)).replace('\\','/'),trajectory_sha256=sha(target/'trajectory.extxyz'),final_sha256=sha(target/'final.extxyz'),geometry=dict(parent_HF_distance_A=float(dist[n,n+1]),added_bond_distances_A=[float(dist[n+2,i]) for i in range(n+3,len(a))],closest_added_heavy_to_parent_H_A=float(dist[n+2,n+1]),minimum_added_to_host_A=float(dist[n+2:,:n].min())),adsorbate_curvature=adsorbate_curvature(a,n) if ok else None,incremental_association_energy_eV=float(energies[-1]-r['reference']['energy_eV']-r['gas_reference']['energy_eV']) if ok and r['reference']['converged'] and r['gas_reference']['converged'] else None,status='relaxed_candidate_not_TS_or_rate_evidence' if ok else 'unconverged_candidate')
        (target/'summary.json').write_text(json.dumps(result,indent=2)+'\n');results.append(result);print(r['surface'],r['coadsorbate'],r['start'],ok,result['max_mobile_force_eV_A'],result['incremental_association_energy_eV'],flush=True)
    report=dict(d);report.update(results=results,initial_campaign='data/intermediate_campaign/summary.json',initial_campaign_sha256=sha(original),refinement_runner='scripts/refine_coadsorbate_intermediates.py',refinement_runner_sha256=sha(Path(__file__)),method=d['method']+'; then BFGS refinement of every candidate',force_tolerance_eV_A=.02,maximum_steps=300,parent_reference_tolerance_eV_A=.02)
    (OUT/'refined_summary.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
