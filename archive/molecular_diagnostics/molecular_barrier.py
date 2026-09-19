"""Small NH3 inversion benchmark: direct DFT or an OMol25-trained ML potential.

This is a molecular saddle-point workflow check, not a surface etching barrier.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import time


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--backend',choices=['deepmd','pyscf'],required=True)
    parser.add_argument('--checkpoint',type=Path)
    parser.add_argument('--basis',default='def2-svp')
    parser.add_argument('--xc',default='PBE')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    import numpy as np
    from scipy.optimize import minimize
    from ase import Atoms
    from ase.io import write
    from ase.units import Hartree,Bohr
    start=time.perf_counter()
    meta=dict(started_utc=datetime.now(timezone.utc).isoformat(),backend=args.backend,
              formula='NH3',charge=0,multiplicity=1,periodic=False,threads=4,
              python=platform.python_version(),status='running',
              purpose='Gas-phase umbrella-inversion saddle diagnostic; not a SiN surface reaction',
              source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              method='Symmetry-reduced minimum and planar saddle optimization; relaxed umbrella scan; full Cartesian finite-difference Hessian')
    records=[]
    if args.backend=='deepmd':
        if args.checkpoint is None:parser.error('--checkpoint required for deepmd')
        import torch
        from deepmd.calculator import DP
        torch.set_num_threads(4)
        calc=DP(model=str(args.checkpoint.resolve()),head='OMol25')
        meta.update(model='DPA-3.3-1M',head='OMol25',label_kind='ML_prediction',
                    checkpoint_sha256=hashlib.sha256(args.checkpoint.read_bytes()).hexdigest(),
                    source='https://huggingface.co/deepmodelingcommunity/DPA-3.3-1M')
        packages=['deepmd-kit','torch','ase','numpy','scipy']
    else:
        from pyscf import gto,dft,lib
        lib.num_threads(4)
        meta.update(label_kind='new_DFT',functional=args.xc,basis=args.basis,
                    scf_tolerance_hartree=1e-10,grid_level=4)
        packages=['pyscf','ase','numpy','scipy']
    meta['packages']={p:importlib.metadata.version(p) for p in packages}
    def save(name,data):
        (args.output/name).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    save('manifest.json',meta)
    def evaluate(positions):
        positions=np.array(positions)
        if args.backend=='deepmd':
            atoms=Atoms('NH3',positions=positions,pbc=False,cell=[20,20,20])
            atoms.info['charge_spin']=np.array([0,1]);atoms.calc=calc
            energy=float(atoms.get_potential_energy());forces=atoms.get_forces()
        else:
            mol=gto.M(atom=list(zip(['N','H','H','H'],positions)),unit='Angstrom',
                      basis=args.basis,charge=0,spin=0,verbose=0)
            mf=dft.RKS(mol);mf.xc=args.xc;mf.grids.level=4;mf.conv_tol=1e-10;mf.max_cycle=100
            energy=float(mf.kernel()*Hartree)
            if not mf.converged:raise RuntimeError('DFT SCF did not converge')
            forces=-mf.nuc_grad_method().kernel()*Hartree/Bohr
        if not np.isfinite(energy) or not np.isfinite(forces).all():raise ValueError('Nonfinite calculation')
        records.append(dict(positions_A=positions.tolist(),energy_eV=energy,forces_eV_A=forces.tolist(),label_kind=meta['label_kind']))
        return energy,forces
    angles=np.arange(3)*2*np.pi/3
    radial=np.column_stack((np.cos(angles),np.sin(angles),np.zeros(3)))
    def geometry(r,q):return np.vstack(([0,0,q],radial*r))
    def objective(v):
        e,f=evaluate(geometry(*v))
        return e,np.array([-np.sum(f[1:]*radial),-f[0,2]])
    try:
        minimum=minimize(objective,[.94,.38],jac=True,method='L-BFGS-B',bounds=[(.7,1.2),(.05,.8)],options={'gtol':1e-5,'ftol':1e-13,'maxiter':80})
        r,q=minimum.x
        saddle=minimize(lambda v:(lambda e,g:(e,g[:1]))(*objective([v[0],0])),[r],jac=True,method='L-BFGS-B',bounds=[(.7,1.2)],options={'gtol':1e-5,'ftol':1e-13,'maxiter':60})
        ts=geometry(saddle.x[0],0);rs=geometry(r,q)
        emin,fmin=evaluate(rs);ets,fts=evaluate(ts)
        profile=[]
        for qi in np.linspace(-q,q,17):
            opt=minimize(lambda v:(lambda e,g:(e,g[:1]))(*objective([v[0],qi])),[r],jac=True,method='L-BFGS-B',bounds=[(.7,1.2)],options={'gtol':1e-5,'ftol':1e-13,'maxiter':40})
            e,f=evaluate(geometry(opt.x[0],qi))
            profile.append(dict(q_A=float(qi),radius_A=float(opt.x[0]),energy_eV=e,relative_energy_eV=e-emin,optimizer_success=bool(opt.success)))
        def hessian(positions,delta):
            h=np.empty((12,12))
            for i in range(12):
                plus=positions.copy().reshape(-1);minus=positions.copy().reshape(-1)
                plus[i]+=delta;minus[i]-=delta
                _,fp=evaluate(plus.reshape(4,3));_,fm=evaluate(minus.reshape(4,3))
                h[:,i]=-(fp-fm).reshape(-1)/(2*delta)
            return (h+h.T)/2
        # Two steps check that the negative curvature is not a finite-step artifact.
        eigs={}
        for delta in [.005,.01]:eigs[str(delta)]=np.linalg.eigvalsh(hessian(ts,delta)).tolist()
        mineigs=np.linalg.eigvalsh(hessian(rs,.005)).tolist()
        # Relax both signs of the unstable umbrella coordinate back to minima.
        connections=[]
        for sign in [-1,1]:
            bounds=[(.7,1.2),(-.8,-.01) if sign<0 else (.01,.8)]
            opt=minimize(objective,[saddle.x[0],sign*.08],jac=True,method='L-BFGS-B',bounds=bounds,options={'gtol':1e-5,'ftol':1e-13,'maxiter':80})
            connections.append(dict(q_A=float(opt.x[1]),radius_A=float(opt.x[0]),energy_eV=float(opt.fun),optimizer_success=bool(opt.success)))
        checks=dict(minimum_optimizer=bool(minimum.success),saddle_optimizer=bool(saddle.success),
                    minimum_stationary=bool(np.linalg.norm(fmin,axis=1).max()<.005),
                    saddle_stationary=bool(np.linalg.norm(fts,axis=1).max()<.005),
                    one_negative_mode_both_steps=all(sum(v<-.01 for v in es)==1 for es in eigs.values()),
                    minimum_no_negative_modes=sum(v<-.01 for v in mineigs)==0,
                    endpoint_connections=all(c['optimizer_success'] and abs(c['energy_eV']-emin)<1e-4 for c in connections) and connections[0]['q_A']*connections[1]['q_A']<0,
                    scan_optimizers=all(p['optimizer_success'] for p in profile))
        summary=dict(barrier_eV=ets-emin,minimum_energy_eV=emin,saddle_energy_eV=ets,
                     minimum_radius_A=float(r),minimum_q_A=float(q),saddle_radius_A=float(saddle.x[0]),
                     minimum_max_force_eV_A=float(np.linalg.norm(fmin,axis=1).max()),
                     saddle_max_force_eV_A=float(np.linalg.norm(fts,axis=1).max()),
                     saddle_cartesian_hessian_eigenvalues_eV_A2=eigs,minimum_hessian_eigenvalues_eV_A2=mineigs,
                     negative_mode_cutoff_eV_A2=-.01,connections=connections,checks=checks,
                     limitations=['No ZPE or finite-temperature free-energy correction','Symmetry-specific saddle search, not a general TS algorithm','Two-sided relaxation is a connectivity check, not a full IRC','NH3 inversion is not a surface diffusion or etching event','PBE reference differs from OMol25 training functional and basis'])
        save('summary.json',summary);save('profile.json',profile);save('evaluations.json',records)
        for name,pos in [('minimum',rs),('saddle',ts),('product',geometry(r,-q))]:write(args.output/(name+'.xyz'),Atoms('NH3',positions=pos))
        meta.update(status='verified' if all(checks.values()) else 'checks_failed',evaluations=len(records),elapsed_s=time.perf_counter()-start)
        print(json.dumps(summary,indent=2),flush=True)
    except Exception as exc:
        meta.update(status='failed',error=repr(exc));save('evaluations.json',records);raise
    finally:
        meta['finished_utc']=datetime.now(timezone.utc).isoformat();save('manifest.json',meta)


if __name__=='__main__':main()
