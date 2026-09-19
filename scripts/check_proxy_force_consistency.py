"""Finite-difference energy/force consistency at saved OMol25 proxy structures."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','2')
import argparse,json,hashlib
from pathlib import Path
import torch,numpy as np
from ase.io import read
from cluster_reaction_paths_refined import ROOT,clone

def main():
    p=argparse.ArgumentParser();p.add_argument('--checkpoint',type=Path,required=True);args=p.parse_args()
    from deepmd.calculator import DP
    torch.set_num_threads(2);calc=DP(model=str(args.checkpoint.resolve()),head='OMol25');src=ROOT/'data/final_cleavage/SiF3_NH2_HF/omol25/DFT_check_geometries.extxyz';results=[]
    for frame,a in enumerate(read(src,':')):
        a.set_constraint();a.calc=clone(calc);f=a.get_forces().copy();xyz=a.positions.copy();checks=[]
        for delta in (.001,.0003):
            for i in range(4,9):
                for j in range(3):
                    a.positions[:]=xyz;a.positions[i,j]+=delta;ep=a.get_potential_energy()
                    a.positions[:]=xyz;a.positions[i,j]-=delta;em=a.get_potential_energy()
                    fd=-(ep-em)/(2*delta);checks.append(dict(atom=i,axis=j,delta_A=delta,force_eV_A=float(f[i,j]),energy_difference_force_eV_A=float(fd),error_eV_A=float(fd-f[i,j])))
        results.append(dict(frame=frame,max_absolute_error_eV_A=max(abs(x['error_eV_A']) for x in checks),checks=checks));print(frame,results[-1]['max_absolute_error_eV_A'],flush=True)
    out=src.parent.parent/'force_consistency.json';out.write_text(json.dumps(dict(source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),runner='scripts/check_proxy_force_consistency.py',runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),model_sha256=hashlib.sha256(args.checkpoint.read_bytes()).hexdigest(),results=results,scope='Calculator consistency only; not accuracy against DFT'),indent=2)+'\n')
if __name__=='__main__':main()
