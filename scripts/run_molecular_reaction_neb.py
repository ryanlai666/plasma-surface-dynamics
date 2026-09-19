"""Search molecular dissociation paths with a locally mobile surface pair.

A failed endpoint search is retained as a failed search, never turned into a TS.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import time
import numpy as np
import torch
from threadpoolctl import threadpool_limits
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.geometry import find_mic
from ase.mep import NEB
from ase.optimize import BFGS
from mace.calculators import MACECalculator
from molecular_surface_campaign import ROOT, independent, snapshot


def curvature(atoms, calc, mobile, delta=.005):
    a = atoms.copy()
    a.calc = independent(calc)
    h = np.zeros((3*len(mobile), 3*len(mobile)))
    for col, (i, axis) in enumerate((i,j) for i in mobile for j in range(3)):
        x = a.positions[i,axis]
        a.positions[i,axis] = x+delta
        plus = a.get_forces()[mobile].reshape(-1).copy()
        a.positions[i,axis] = x-delta
        minus = a.get_forces()[mobile].reshape(-1).copy()
        a.positions[i,axis] = x
        h[:,col] = -(plus-minus)/(2*delta)
    return np.linalg.eigh((h+h.T)/2)


def run(surface_name, species, calc, checkpoint):
    source = ROOT/'data/surface_paths'/surface_name/species/'mace_approach'
    out = source.parent/'molecular_neb'
    out.mkdir(exist_ok=False)
    base = read(source/'slab.extxyz')
    gas = read(source/'molecule.extxyz')
    parent = json.loads((source/'summary.json').read_text())
    site = parent['surface_site_index']
    preferred = 'N' if 'Si3N4' in surface_name else 'O' if 'quartz' in surface_name else 'Si'
    if species in ['F2','Cl2']:
        preferred = 'Si'
    near = [i for i,a in enumerate(base) if a.symbol==preferred and i!=site and a.position[2]>base.positions[:,2].max()-2.2]
    vectors = {i:find_mic(base.positions[i]-base.positions[site],base.cell,base.pbc)[0] for i in near}
    neighbor = min(near,key=lambda i:np.linalg.norm(vectors[i]))
    target = base.positions[site]+vectors[neighbor]
    away = vectors[neighbor].copy(); away[2]=0; away /= np.linalg.norm(away)
    gas.positions -= gas.positions[0]
    breaking = 2 if species=='H2O' else 1
    bond0 = float(np.linalg.norm(gas.positions[breaking]-gas.positions[0]))
    initial = base.copy()+gas.copy()
    final = base.copy()+gas.copy()
    n = len(base)
    anchor_height = 3.0 if species=='CH3F' else 2.5
    anchor = base.positions[site]+[0,0,anchor_height]
    rotated = gas.copy()
    rotated.rotate(rotated.positions[breaking],target+[0,0,1.1]-anchor,center=(0,0,0))
    initial.positions[n:] = rotated.positions+anchor
    first_height = 2.12 if gas[0].symbol=='Cl' else 1.75 if gas[0].symbol=='F' else 1.65
    final.positions[n:] = rotated.positions+base.positions[site]+[0,0,first_height]
    second_height = (1.52 if preferred=='Si' else 1.05) if gas[breaking].symbol=='H' else (2.12 if gas[breaking].symbol=='Cl' else 1.75 if gas[breaking].symbol=='F' else 1.48)
    destination = target+[0,0,second_height]+away*.35
    fragment = list(range(1,len(gas))) if species=='CH3F' else [breaking]
    final.positions[n+np.array(fragment)] += destination-final.positions[n+breaking]
    mobile = [site,neighbor]+list(range(n,len(initial)))
    fixed = [i for i in range(len(initial)) if i not in mobile]
    record = dict(surface=surface_name,species=species,status='endpoint_search',mobile_indices=mobile,
                  substrate_atoms=n,active_surface_pair=[site,neighbor],source=parent['source'],
                  proposed_reaction=f'{species} + reactive surface pair -> dissociated adsorbates',
                  breaking_bond_indices=[n,n+breaking],gas_bond_length_A=bond0,
                  method='MACE-MP-0b2; CI-NEB with two surface atoms and all molecular atoms mobile',
                  model_sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
                  runner='scripts/run_molecular_reaction_neb.py',script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  limitations=['Ideal unpassivated cleavage, one coverage and site; no slab convergence study',
                               'Only two substrate atoms move; bulk-trained ML model is unvalidated for surface chemistry',
                               'Electronic potential energies only; no excitation, ZPE or entropy',
                               'A highest image is not called a verified TS without curvature and connectivity checks'])
    def save():
        (out/'summary.json').write_text(json.dumps(record,indent=2)+'\n')
    save()
    start=time.perf_counter()
    endpoints=[]
    for label,a in [('IS',initial),('FS',final)]:
        a.set_constraint(FixAtoms(indices=fixed));a.calc=independent(calc)
        write(out/(label+'_guess.extxyz'),a,write_results=False)
        opt=BFGS(a,logfile=str(out/(label+'.log')),maxstep=.1)
        ok=bool(opt.run(fmax=.02,steps=220))
        saved=snapshot(a);write(out/(label+'.extxyz'),saved)
        eig,_=curvature(a,calc,mobile)
        endpoints.append(saved)
        record[label]=dict(converged=ok,energy_eV=float(a.get_potential_energy()),
                           max_mobile_force_eV_A=float(np.linalg.norm(a.get_forces()[mobile],axis=1).max()),
                           breaking_bond_A=float(a.get_distance(n,n+breaking,mic=True)),
                           surface_pair_distance_A=float(a.get_distance(site,neighbor,mic=True)),
                           hessian_eigenvalues_eV_A2=eig.tolist(),no_negative_curvature=bool(eig.min()>-.02))
        save()
    displacement=find_mic(endpoints[1].positions-endpoints[0].positions,base.cell,base.pbc)[0]
    record['endpoint_mobile_RMS_displacement_A']=float(np.sqrt(np.mean(displacement[mobile]**2)))
    if not all(record[k]['converged'] and record[k]['no_negative_curvature'] for k in ['IS','FS']):
        record['status']='endpoint_validation_failed';record['seconds']=time.perf_counter()-start;save();print(surface_name,species,record['status'],flush=True);return
    if abs(record['FS']['breaking_bond_A']-record['IS']['breaking_bond_A'])<.35:
        record['status']='no_distinct_dissociation_endpoints';record['seconds']=time.perf_counter()-start;save();print(surface_name,species,record['status'],flush=True);return
    images=[endpoints[0].copy()]+[endpoints[0].copy() for _ in range(9)]+[endpoints[1].copy()]
    for a in images:a.calc=independent(calc)
    neb=NEB(images,k=.25,climb=False,method='improvedtangent')
    neb.interpolate(method='idpp',apply_constraint=True)
    opt=BFGS(neb,logfile=str(out/'neb.log'),maxstep=.1)
    opt.attach(lambda:write(out/'checkpoint.extxyz',images),interval=5)
    opt.run(fmax=.15,steps=80);neb.climb=True
    ok=bool(opt.run(fmax=.035,steps=180))
    residual=float(np.linalg.norm(neb.get_forces().reshape(-1,3),axis=1).max())
    saved=[snapshot(a) for a in images];write(out/'images.extxyz',saved)
    es=np.array([a.get_potential_energy() for a in saved]);peak=int(es.argmax())
    xyz=np.array([a.positions for a in saved]);arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(xyz,axis=0).reshape(len(images)-1,-1),axis=1))]
    rows=[dict(image=i,arc_length_A=float(arc[i]),energy_eV=float(es[i]),relative_energy_eV=float(es[i]-es[0]),breaking_bond_A=float(a.get_distance(n,n+breaking,mic=True)),surface_pair_distance_A=float(a.get_distance(site,neighbor,mic=True))) for i,a in enumerate(saved)]
    with (out/'energies.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    record.update(status='neb_converged' if ok else 'neb_unconverged',images=len(images),peak_image=peak,
                  peak_above_IS_eV=float(es.max()-es[0]),FS_minus_IS_eV=float(es[-1]-es[0]),max_neb_force_eV_A=residual)
    if 0<peak<len(images)-1:
        eig,vec=curvature(saved[peak],calc,mobile)
        record['peak_hessian_eigenvalues_eV_A2']=eig.tolist()
        record['index_one_in_mobile_subspace']=bool(sum(eig<-.02)==1)
        if ok and record['index_one_in_mobile_subspace']:
            connections=[]
            for sign in [-1,1]:
                a=saved[peak].copy();a.positions[mobile]+=sign*.12*vec[:,0].reshape(-1,3);a.calc=independent(calc)
                opt=BFGS(a,logfile=str(out/f'downhill_{sign}.log'),maxstep=.1);done=bool(opt.run(fmax=.02,steps=180))
                e=float(a.get_potential_energy());write(out/f'downhill_{sign}.extxyz',snapshot(a))
                errors=[float(np.sqrt(np.mean(find_mic(a.positions[mobile]-b.positions[mobile],base.cell,base.pbc)[0]**2))) for b in endpoints]
                connections.append(dict(sign=sign,converged=done,energy_eV=e,endpoint_RMS_errors_A=errors,breaking_bond_A=float(a.get_distance(n,n+breaking,mic=True))))
            record['downhill_connections']=connections
            matched=[int(np.argmin(c['endpoint_RMS_errors_A'])) if min(c['endpoint_RMS_errors_A'])<.25 else -1 for c in connections]
            record['endpoint_connectivity_confirmed']=sorted(matched)==[0,1] and all(c['converged'] for c in connections)
    record['seconds']=time.perf_counter()-start;save();print(surface_name,species,record['status'],record.get('peak_above_IS_eV'),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases',nargs='+',required=True,help='surface/species pairs')
    args=parser.parse_args();torch.set_num_threads(2);threadpool_limits(limits=1)
    checkpoint=ROOT/'.cache/mace/mp_0b2_small.model'
    calc=MACECalculator(model_paths=str(checkpoint),device='cpu',default_dtype='float64')
    for case in args.cases:
        run(*case.split('/'),calc,checkpoint)


if __name__=='__main__':
    main()
