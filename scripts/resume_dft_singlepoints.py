"""Resume failed fixed-geometry DFT with a documented Newton SCF fallback."""
from pathlib import Path
import csv,hashlib,json,sys
import numpy as np
from ase.io import read,write
from ase.calculators.singlepoint import SinglePointCalculator
from ase.units import Hartree,Bohr
from pyscf import gto,dft,lib
ROOT=Path(__file__).resolve().parents[1];out=Path(sys.argv[1]);lib.num_threads(2)
meta=json.loads((out/'summary.json').read_text());images=read(out/'input.extxyz',':');saved=read(out/'images.extxyz',':');rows=meta['frames']
meta['continuations']=[dict(runner='scripts/resume_dft_singlepoints.py',script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),start_image=len(rows),reason='Direct SCF failed at this image; use Newton fallback without changing geometry or functional/basis')]
for i in range(len(rows),len(images)):
    a=images[i];mol=gto.M(atom=list(zip(a.get_chemical_symbols(),a.positions)),basis='def2-svp',charge=0,spin=0,unit='Angstrom',verbose=0)
    mf=dft.RKS(mol).density_fit();mf.xc='PBE';mf.grids.level=3;mf.conv_tol=1e-9;mf.max_cycle=200
    mf.kernel();solver='DIIS'
    if not mf.converged:
        dm=mf.make_rdm1();mf=mf.newton();mf.max_cycle=180;mf.conv_tol=1e-9;mf.kernel(dm0=dm);solver='Newton after failed DIIS'
    if not mf.converged:
        meta['status']='SCF_failed';meta['failed_image']=i;(out/'summary.json').write_text(json.dumps(meta,indent=2)+'\n');raise RuntimeError('Newton SCF failed at '+str(i))
    energy=float(mf.e_tot*Hartree);forces=-mf.nuc_grad_method().kernel()*Hartree/Bohr;a.calc=SinglePointCalculator(a,energy=energy,forces=forces);saved.append(a);write(out/'images.extxyz',saved)
    rows.append(dict(image=i,energy_eV=energy,relative_energy_eV=energy-rows[0]['energy_eV'],max_raw_force_eV_A=float(np.linalg.norm(forces,axis=1).max()),SCF_converged=True,solver=solver))
    with (out/'energies.csv').open('w',newline='') as f:
        keys=list(rows[-1]);w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
    meta['frames']=rows;meta['status']='complete' if len(saved)==len(images) else 'running';meta['electronic_state_note']='Converged restricted-singlet SCF; no proof of lowest electronic or spin state';(out/'summary.json').write_text(json.dumps(meta,indent=2)+'\n');print(i,energy,solver,flush=True)
