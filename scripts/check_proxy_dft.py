"""Independent DFT energy/force and SCF-stability checks on ML proxy geometries."""
from pathlib import Path
import argparse,json,hashlib,time
import numpy as np
from ase.io import read,write
from ase.calculators.singlepoint import SinglePointCalculator
from ase.units import Hartree,Bohr
from pyscf import gto,dft,lib
ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',type=Path);p.add_argument('output',type=Path);args=p.parse_args();lib.num_threads(2)
    args.output.mkdir(parents=True,exist_ok=False);raw=args.input.read_bytes();(args.output/'input.extxyz').write_bytes(raw);images=read(args.output/'input.extxyz',':');all_results=[]
    meta=dict(source=str(args.input),source_sha256=hashlib.sha256(raw).hexdigest(),runner='scripts/check_proxy_dft.py',runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='DFT single points and forces on saved ML structures; not DFT-optimized endpoints or TS',charge=0,multiplicity=1,results=all_results)
    for basis in ('def2-svp','def2-tzvp'):
        result=[];saved=[]
        for i,a in enumerate(images):
            if a.pbc.any():raise ValueError('Molecular proxy only')
            mol=gto.M(atom=list(zip(a.get_chemical_symbols(),a.positions)),basis=basis,charge=0,spin=0,verbose=0)
            mf=dft.RKS(mol).density_fit();mf.xc='PBE';mf.grids.level=3;mf.conv_tol=1e-9;mf.max_cycle=150
            t=time.perf_counter();e=mf.kernel();attempts=[];stable=False
            for trial in range(3):
                if not mf.converged:
                    mf=mf.newton();mf.max_cycle=100;e=mf.kernel()
                if not mf.converged:break
                mo,_,stable,_=mf.stability(internal=True,external=False,return_status=True)
                attempts.append(dict(energy_Hartree=float(e),internal_stable=bool(stable)))
                if stable:break
                dm=mf.make_rdm1(mo,mf.mo_occ);mf=mf.newton();mf.max_cycle=100;e=mf.kernel(dm0=dm)
            if not mf.converged or not stable:raise RuntimeError('SCF convergence or internal-stability check failed')
            forces=-mf.nuc_grad_method().kernel()*Hartree/Bohr;mlforce=a.get_forces(apply_constraint=False);mle=a.get_potential_energy();b=a.copy();b.info=dict(calculation_kind='PBE fixed-geometry diagnostic; not a DFT stationary point',basis=basis,charge=0,multiplicity=1);b.calc=SinglePointCalculator(b,energy=float(e*Hartree),forces=forces);saved.append(b)
            fixed=set(a.constraints[0].get_indices()) if a.constraints else set();mobile=[j for j in range(len(a)) if j not in fixed]
            result.append(dict(image=i,DFT_energy_eV=float(e*Hartree),ML_energy_eV=float(mle),max_mobile_DFT_force_eV_A=float(np.linalg.norm(forces[mobile],axis=1).max()),mobile_force_RMSE_eV_A=float(np.sqrt(np.mean((forces[mobile]-mlforce[mobile])**2))),SCF_converged=bool(mf.converged),internal_stable=bool(stable),external_spin_stability_checked=False,SCF_attempts=attempts,elapsed_s=time.perf_counter()-t))
            write(args.output/(basis+'_images.extxyz'),saved);print(basis,i,result[-1]['DFT_energy_eV'],result[-1]['max_mobile_DFT_force_eV_A'],flush=True)
            meta['in_progress_basis']=basis;meta['partial_results']=result;(args.output/'summary.json').write_text(json.dumps(meta,indent=2)+'\n')
        de=np.array([x['DFT_energy_eV'] for x in result]);me=np.array([x['ML_energy_eV'] for x in result])
        all_results.append(dict(basis=basis,functional='PBE',density_fitting=True,grid_level=3,SCF_tolerance_Hartree=1e-9,frames=result,relative_DFT_energies_eV=(de-de[0]).tolist(),relative_ML_energies_eV=(me-me[0]).tolist(),relative_energy_disagreement_eV=((me-me[0])-(de-de[0])).tolist(),energies_are_barriers=False))
        meta.pop('partial_results',None);meta.pop('in_progress_basis',None);(args.output/'summary.json').write_text(json.dumps(meta,indent=2)+'\n')
    meta['status']='complete';meta['limitations']=['Same closed-shell electronic sector only; no external spin-stability proof','Two basis sets at PBE, no functional-convergence claim','No DFT geometry optimization, saddle verification or thermochemistry','Fixed molecular SiF3 frame is not an embedded surface environment'];(args.output/'summary.json').write_text(json.dumps(meta,indent=2)+'\n')
if __name__=='__main__':main()
