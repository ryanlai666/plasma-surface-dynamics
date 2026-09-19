"""Densify a saved NEB path and re-optimize every inserted image before publishing."""
from pathlib import Path
import argparse,copy,csv,hashlib,json,time,sys
import numpy as np
import torch
from threadpoolctl import threadpool_limits
from ase.io import read,write
from ase.mep import NEB
from ase.optimize import BFGS
from mace.calculators import MACECalculator
ROOT=Path(__file__).resolve().parents[1]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('folder',type=Path);args=ap.parse_args();p=args.folder
    torch.set_num_threads(2);threadpool_limits(limits=1)
    old=read(p/'images.extxyz',':');s=json.loads((p/'summary.json').read_text());attempt=p/'attempts/coarse_9';attempt.mkdir(parents=True,exist_ok=True)
    if (attempt/'summary.json').exists():raise ValueError('Refinement already recorded; choose a new explicit attempt')
    for name in ['images.extxyz','energies.csv','summary.json','curvature.json','second_model.json','neb.log']:
        src=p/name
        if src.exists():(attempt/name).write_bytes(src.read_bytes())
    images=[]
    for a,b in zip(old[:-1],old[1:]):
        images.append(a.copy());mid=a.copy();mid.positions=(a.positions+b.positions)/2;images.append(mid)
    images.append(old[-1].copy())
    calc=MACECalculator(model_paths=str(ROOT/'.cache/mace/mp_0b2_small.model'),device='cpu',default_dtype='float64')
    for a in images:c=copy.copy(calc);c.atoms=None;c.results={};a.calc=c
    neb=NEB(images,climb=False,k=.3,method='improvedtangent');opt=BFGS(neb,logfile=str(p/'refinement.log'),maxstep=.1);start=time.perf_counter()
    opt.attach(lambda:write(p/'checkpoint.extxyz',images),interval=5)
    opt.run(fmax=.07,steps=60);neb.climb=True;ok=bool(opt.run(fmax=.015,steps=160));res=float(np.linalg.norm(neb.get_forces().reshape(-1,3),axis=1).max())
    es=np.array([a.get_potential_energy() for a in images]);xyz=np.array([a.positions for a in images]);arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(xyz,axis=0).reshape(len(images)-1,-1),axis=1))]
    rows=[dict(image=i,arc_length_A=arc[i],energy_eV=e,relative_energy_eV=e-es[0],max_unconstrained_force_eV_A=np.linalg.norm(a.get_forces(),axis=1).max()) for i,(e,a) in enumerate(zip(es,images))]
    with (p/'energies.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    write(p/'images.extxyz',images)
    s.update(command=sys.argv,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),status='converged' if ok else 'unconverged',images=len(images),peak_image=int(es.argmax()),peak_above_IS_eV=float(es.max()-es[0]),FS_minus_IS_eV=float(es[-1]-es[0]),max_neb_force_eV_A=res,seconds=time.perf_counter()-start,refinement='17-image NEB from archived 9-image path; inserted starting guesses were force-optimized, not published as interpolated frames',rendering='Final evaluated images only',initial_images_sha256=hashlib.sha256((attempt/'images.extxyz').read_bytes()).hexdigest())
    (p/'summary.json').write_text(json.dumps(s,indent=2)+'\n');print(json.dumps(s),flush=True)
if __name__=='__main__':main()
