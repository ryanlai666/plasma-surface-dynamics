"""Evaluate three sites and chemically distinct orientations for each molecule/surface."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
import torch
from threadpoolctl import threadpool_limits
from ase.io import read,write
from ase.geometry import find_mic
from mace.calculators import MACECalculator
from molecular_surface_campaign import ROOT, SURFACES, MOLECULES, HEIGHTS, independent, snapshot


def sites(base):
    top_si=max(a.position[2] for a in base if a.symbol=='Si')
    silicon=[i for i,a in enumerate(base) if a.symbol=='Si' and a.position[2]>top_si-.02]
    center=.5*(base.cell[0]+base.cell[1])
    i=min(silicon,key=lambda j:np.linalg.norm(base.positions[j,:2]-center[:2]))
    others=[j for j,a in enumerate(base) if a.symbol!='Si' and a.position[2]>base.positions[:,2].max()-2]
    if others:
        vectors={j:find_mic(base.positions[j]-base.positions[i],base.cell,base.pbc)[0] for j in others}
        j=min(others,key=lambda j:np.linalg.norm(vectors[j]));v=vectors[j]
        return [('atop_Si',base.positions[i].copy(),[i]),('atop_'+base[j].symbol,base.positions[i]+v,[j]),('bridge_Si_'+base[j].symbol,base.positions[i]+.5*v,[i,j])]
    vectors={j:find_mic(base.positions[j]-base.positions[i],base.cell,base.pbc)[0] for j in silicon if j!=i}
    order=sorted(vectors,key=lambda j:np.linalg.norm(vectors[j]));j=order[0];v=vectors[j]
    k=next(k for k in order[1:] if abs(np.cross(v,vectors[k])[2])>1)
    angle=np.dot(v,vectors[k])/np.linalg.norm(v)/np.linalg.norm(vectors[k])
    # For square surface nets use a four-fold hollow, for triangular nets the triangle center.
    hollow=base.positions[i]+(v+vectors[k])/(2 if abs(angle)<.1 else 3)
    return [('atop_Si',base.positions[i].copy(),[i]),('bridge_Si_Si',base.positions[i]+.5*v,[i,j]),('hollow',hollow,[i,j,k])]


def main():
    torch.set_num_threads(2);threadpool_limits(limits=1)
    model=ROOT/'.cache/mace/mp_0b2_small.model'
    calc=MACECalculator(model_paths=str(model),device='cpu',default_dtype='float64')
    count=0
    for surface in SURFACES:
        for species in MOLECULES:
            parent=ROOT/'data/surface_paths'/surface/species
            base=read(parent/'mace_approach/slab.extxyz');gas=read(parent/'mace_approach/molecule.extxyz')
            root=parent/'mace_site_orientation';root.mkdir(exist_ok=False)
            base.calc=independent(calc);e_slab=float(base.get_potential_energy())
            rows=[]
            for site_name,site_position,site_atoms in sites(base):
                orientations=['upright','parallel'] if species in ['F2','Cl2','H2'] else ['upright','parallel','flipped']
                for orientation in orientations:
                    mol=gas.copy();mol.calc=None
                    if orientation=='parallel':mol.rotate(90,'y',center=(0,0,0))
                    if orientation=='flipped':mol.rotate(180,'y',center=(0,0,0))
                    mol.positions[:,:2]-=mol.positions[:,:2].mean(axis=0)
                    mol.positions[:,2]-=mol.positions[:,2].min()
                    mol.calc=independent(calc);e_mol=float(mol.get_potential_energy())
                    path=root/site_name/orientation;path.mkdir(parents=True)
                    point_rows=[];images=[]
                    for n,height in enumerate(HEIGHTS):
                        moved=mol.copy();moved.positions += [site_position[0],site_position[1],base.positions[:,2].max()+height]
                        a=base+moved;a.calc=independent(calc)
                        a.info.update(surface=surface,species=species,site=site_name,orientation=orientation,height_A=height,substrate_atoms=len(base),calculation_kind='Rigid evaluated approach scan; not NEB')
                        image=snapshot(a);images.append(image);e=float(image.get_potential_energy())
                        row=dict(image=n,height_A=height,energy_eV=e,interaction_energy_eV=e-e_slab-e_mol)
                        point_rows.append(row);rows.append(dict(site=site_name,orientation=orientation,**row))
                    write(path/'images.extxyz',images)
                    with (path/'energies.csv').open('w',newline='') as handle:
                        writer=csv.DictWriter(handle,fieldnames=point_rows[0]);writer.writeheader();writer.writerows(point_rows)
                    meta=dict(surface=surface,species=species,site=site_name,site_position_A=site_position.tolist(),site_atom_indices=site_atoms,
                              orientation=orientation,orientation_convention='Rotate the optimized gas geometry around y; center molecular centroid laterally over site; lowest atom defines height',
                              substrate_atoms=len(base),images=len(images),slab_energy_eV=e_slab,molecular_layer_energy_eV=e_mol,
                              minimum_interaction_energy_eV=min(r['interaction_energy_eV'] for r in point_rows),
                              peak_is_TS=False,method='Rigid evaluated scan, not a reaction path',
                              model_sha256=hashlib.sha256(model.read_bytes()).hexdigest(),
                              runner='scripts/screen_molecular_sites.py',script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
                    (path/'summary.json').write_text(json.dumps(meta,indent=2)+'\n');count+=1
            with (root/'all_energies.csv').open('w',newline='') as handle:
                writer=csv.DictWriter(handle,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
            print(surface,species,len(rows),'evaluations; cumulative scans',count,flush=True)


if __name__=='__main__':main()
