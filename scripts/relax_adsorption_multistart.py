"""Relax HF and the upper substrate from two distinct screened starts per surface."""
from pathlib import Path
import csv,hashlib,json,time
import numpy as np
import torch
from threadpoolctl import threadpool_limits
from ase.io import read,write
from ase.constraints import FixAtoms
from ase.optimize import FIRE
from mace.calculators import MACECalculator
from molecular_surface_campaign import ROOT,SURFACES,independent,snapshot

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def relax(a,calc,out):
    out.mkdir(parents=True,exist_ok=False)
    a.calc=independent(calc);frames=[]
    opt=FIRE(a,logfile=str(out/'optimization.log'),dt=.03,maxstep=.08)
    def record():
        frames.append(snapshot(a))
        write(out/'checkpoint.extxyz',frames[-1])
    opt.attach(record,interval=1)
    ok=bool(opt.run(fmax=.04,steps=300));write(out/'trajectory.extxyz',frames)
    energies=[float(f.get_potential_energy()) for f in frames]
    result=dict(converged=ok,steps=opt.nsteps,energy_eV=energies[-1],initial_energy_eV=energies[0],relaxation_energy_eV=energies[-1]-energies[0],max_mobile_force_eV_A=float(np.linalg.norm(a.get_forces(),axis=1).max()),frames=len(frames),trajectory_sha256=digest(out/'trajectory.extxyz'))
    with (out/'energies.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['step','energy_eV','relative_to_start_eV']);w.writerows((i,e,e-energies[0]) for i,e in enumerate(energies))
    return result,a

def main():
    torch.set_num_threads(2);threadpool_limits(limits=1)
    model=ROOT/'.cache/mace/mp_0b2_small.model';calc=MACECalculator(model_paths=str(model),device='cpu',default_dtype='float64')
    all_results=[]
    for surface in SURFACES:
        parent=ROOT/'data/surface_paths'/surface/'HF';out=parent/'relaxed_adsorption'
        if (out/'summary.json').exists():
            all_results.extend(json.loads((out/'summary.json').read_text())['starts']);continue
        out.mkdir(exist_ok=True)
        base=read(parent/'mace_approach/slab.extxyz');cut=float(np.median(base.positions[:,2]));fixed=np.flatnonzero(base.positions[:,2]<=cut).tolist();base.set_constraint(FixAtoms(indices=fixed))
        reference,relaxed_base=relax(base,calc,out/'bare_slab')
        gas=read(parent/'mace_approach/molecule.extxyz');gas.calc=independent(calc);gas_energy=float(gas.get_potential_energy())
        rows=list(csv.DictReader((parent/'mace_site_orientation/all_energies.csv').open()))
        rows=[r for r in rows if float(r['height_A'])>=1.9]
        rows.sort(key=lambda r:float(r['interaction_energy_eV']))
        chosen=[]
        for r in rows:
            if not chosen or (r['site']!=chosen[0]['site'] and r['orientation']!=chosen[0]['orientation']):chosen.append(r)
            if len(chosen)==2:break
        results=[]
        for i,r in enumerate(chosen):
            src=parent/'mace_site_orientation'/r['site']/r['orientation']/'images.extxyz';a=read(src,int(r['image']));a.set_constraint(FixAtoms(indices=fixed));initial=a.positions.copy()
            result,final=relax(a,calc,out/f'start_{i+1}')
            result.update(surface=surface,species='HF',start=i+1,site=r['site'],orientation=r['orientation'],start_height_A=float(r['height_A']),source=str(src.relative_to(ROOT)),source_sha256=digest(src),source_frame=int(r['image']),fixed_indices=fixed,mobile_substrate_atoms=len(base)-len(fixed),molecule_atoms=2,substrate_max_displacement_A=float(np.linalg.norm(final.positions[:len(base)]-initial[:len(base)],axis=1).max()),HF_distance_A=float(final.get_distance(len(base),len(base)+1,mic=True)),adsorption_energy_eV=float(result['energy_eV']-reference['energy_eV']-gas_energy) if result['converged'] and reference['converged'] else None,folder=str((out/f'start_{i+1}').relative_to(ROOT)))
            result['geometry_classification']='HF bond separated (>1.4 A)' if result['HF_distance_A']>1.4 else 'HF bond retained (<=1.4 A)'
            (out/f'start_{i+1}'/'summary.json').write_text(json.dumps(result,indent=2)+'\n');results.append(result);print(surface,i+1,result['converged'],result['max_mobile_force_eV_A'],flush=True)
        d=dict(method='MACE-MP-0b2 small; FIRE; upper substrate and HF mobile; lower half fixed',force_tolerance_eV_A=.04,max_steps=300,fixed_indices=fixed,bare_slab=reference,gas_energy_eV=gas_energy,gas_reference='Isolated neutral HF, same MACE model; one molecule per periodic surface supercell',model_sha256=digest(model),runner='scripts/relax_adsorption_multistart.py',script_sha256=digest(Path(__file__)),starts=results,limitations='No Hessian, TS or DFT validation. Ideal unpassivated slab; finite cell and fixed-bottom constraint. Local relaxations do not prove global minima or a reaction path.')
        (out/'summary.json').write_text(json.dumps(d,indent=2)+'\n');all_results.extend(results)
    (ROOT/'docs/dry_etch_results/relaxed_adsorption.json').write_text(json.dumps(all_results,indent=2)+'\n')
if __name__=='__main__':main()
