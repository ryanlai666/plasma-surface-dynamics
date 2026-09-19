"""Relaxed three-coordinate diagnostic for local cleavage, not a TS calculation."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','2')
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import torch
from ase.io import read,write
from ase.constraints import FixAtoms
from fixed_frame_constraints import FixedFrameDistances
from ase.optimize import BFGS
from cluster_reaction_paths_refined import clone,freeze,ROOT

def main():
    p=argparse.ArgumentParser();p.add_argument('--checkpoint',type=Path,required=True);args=p.parse_args()
    from deepmd.calculator import DP
    torch.set_num_threads(2);calc=DP(model=str(args.checkpoint.resolve()),head='OMol25')
    src=ROOT/'data/final_cleavage/SiF3_NH2_HF/omol25/IS.extxyz';out=src.parent.parent/'fixed_frame_three_coordinate_scan';out.mkdir(exist_ok=False)
    a=read(src);a.calc=clone(calc);start=a.get_distance(0,4);nh=a.get_distance(4,8);rows=[];frames=[]
    targets=list(zip(np.linspace(start,3.3,11),np.linspace(nh,1.025,11),np.linspace(a.get_distance(0,7),1.60,11)))
    fixed_xyz=a.positions[:4].copy()
    for direction,seq in [('forward',list(enumerate(targets))),('reverse',list(enumerate(targets))[::-1])]:
        for i,(sin,nhdist,sif) in seq:
            # Translate N and its original H atoms together before imposing the new distances.
            a.set_constraint();v=a.positions[4]-a.positions[0];shift=v*(sin/np.linalg.norm(v)-1);a.positions[[4,5,6]]+=shift
            v=a.positions[8]-a.positions[4];a.positions[8]=a.positions[4]+v*nhdist/np.linalg.norm(v)
            a.set_constraint(FixedFrameDistances([0,1,2,3],[(0,4),(4,8),(0,7)],[sin,nhdist,sif]));a.set_positions(a.positions.copy())
            opt=BFGS(a,logfile=str(out/f'{direction}_{i:02d}.log'),maxstep=.035);ok=bool(opt.run(fmax=.04,steps=160))
            physical=a.get_forces(apply_constraint=False);projected=a.get_forces();assert np.max(abs(a.positions[:4]-fixed_xyz))<1e-10;b=freeze(a);b.set_constraint(FixAtoms(indices=[0,1,2,3]));b.info.update(calculation_kind='Evaluated relaxed three-coordinate scan; not a minimum-energy path or TS',direction=direction,target_SiN_A=float(sin),target_NH_A=float(nhdist),target_SiF_A=float(sif));frames.append(b)
            rows.append(dict(direction=direction,index=i,energy_eV=float(a.get_potential_energy()),target_SiN_A=float(sin),target_NH_A=float(nhdist),target_SiF_A=float(sif),SiN_A=float(a.get_distance(0,4)),NH_A=float(a.get_distance(4,8)),SiF_A=float(a.get_distance(0,7)),HF_A=float(a.get_distance(7,8)),projected_force_eV_A=float(np.linalg.norm(projected,axis=1).max()),mobile_physical_force_eV_A=float(np.linalg.norm(physical[4:],axis=1).max()),converged=ok))
            write(out/'evaluated_geometries.extxyz',frames);(out/'summary.json').write_text(json.dumps(dict(runner='scripts/scan_final_cleavage_fixed_frame.py',runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),constraint_sha256=hashlib.sha256((ROOT/'scripts/fixed_frame_constraints.py').read_bytes()).hexdigest(),model_sha256=hashlib.sha256(args.checkpoint.read_bytes()).hexdigest(),status='complete' if len(rows)==22 else 'running',scope='Fixed-frame constrained Si-N elongation, proton transfer and Si-F approach; forward/reverse hysteresis diagnostic; no TS or kinetic barrier',rate_enabled=False,results=rows),indent=2)+'\n');print(direction,i,rows[-1],flush=True)
if __name__=='__main__':main()
