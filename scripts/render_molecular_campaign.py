"""Render evaluated molecular geometries; never interpolate displayed frames."""
from pathlib import Path
import csv,io,json,hashlib,argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from ase.io import read
from ase.data import covalent_radii,atomic_numbers
ROOT=Path(__file__).resolve().parents[1]
COLORS={'Si':'#819db7','N':'#466ddd','O':'#df594e','F':'#31b882','Cl':'#a2bc36','H':'#eeeeee','C':'#525867'}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def render(folder, title, kind, rows, images, curves=None):
    energies=np.array([float(r['interaction_energy_eV'] if kind=='scan' else r['relative_energy_eV']) for r in rows])
    axis=np.array([float(r['height_A']) if kind=='scan' else i for i,r in enumerate(rows)])
    xyz=np.array([a.positions for a in images]);lo=xyz.min(axis=(0,1));hi=xyz.max(axis=(0,1));center=(lo+hi)/2;span=max(hi-lo)
    path_meta=json.loads((folder/'summary.json').read_text())
    if folder.name=='mace_local_saddle':title='Si3N4 / OH + H: N-to-N H transfer'
    frames=[];frame_indices=[]
    for n,a in enumerate(images):
        if n and np.array_equal(a.positions,images[n-1].positions):continue
        frame_indices.append(n)
        assert abs(a.get_potential_energy()-float(rows[n]['energy_eV']))<1e-7
        fig=plt.figure(figsize=(10,4.8),dpi=95);ax=fig.add_subplot(121,projection='3d');ep=fig.add_subplot(122)
        pos=a.positions
        for i in range(len(a)):
            for j in range(i):
                cutoff=1.18*(covalent_radii[a.numbers[i]]+covalent_radii[a.numbers[j]])
                if .45<np.linalg.norm(pos[i]-pos[j])<cutoff:ax.plot(*pos[[i,j]].T,c='#8895a3',lw=1,alpha=.6)
        ax.scatter(*pos.T,s=[48 if x=='H' else 80 for x in a.get_chemical_symbols()],c=[COLORS[x] for x in a.get_chemical_symbols()],edgecolors='#38424b',linewidths=.45,depthshade=True)
        for setter,c in zip([ax.set_xlim,ax.set_ylim,ax.set_zlim],center):setter(c-span*.56,c+span*.56)
        ax.set_box_aspect([1,1,1]);ax.view_init(elev=24,azim=-58);ax.set_axis_off()
        state=''
        if kind!='scan':
            state='IS' if n==0 else 'FS' if n==len(images)-1 else ('TS (constrained)' if path_meta.get('endpoint_connectivity_confirmed') else 'candidate peak') if n==path_meta.get('peak_image') else f'image {n}'
        ax.set_title(title+(' / '+state if state else ''),fontsize=10)
        ep.plot(axis,energies,'o-',color='#177f93',ms=4,lw=1,label='ML evaluated energies')
        dft=folder.parent/'dft_path'
        if kind!='scan' and (dft/'summary.json').exists():
            ds=json.loads((dft/'summary.json').read_text())
            if ds.get('status')=='complete' and ds.get('input_sha256')==sha(folder/'images.extxyz'):
                de=np.array([r['energy_eV'] for r in ds['frames']]);de-=de[0];de=np.array([value if r.get('accepted_for_comparison',True) else np.nan for value,r in zip(de,ds['frames'])]);ep.plot(axis,de,'s--',color='#9864b1',ms=4,label='PBE/def2-SVP (gaps = rejected SCF)');ep.legend(fontsize=7)
        ep.scatter(axis[n],energies[n],s=85,c='#dc623c',zorder=5)
        ep.set(xlabel='Lowest molecular atom above slab (angstrom)' if kind=='scan' else 'Evaluated geometry index',ylabel='Interaction energy (eV)' if kind=='scan' else 'Energy relative to IS (eV)')
        if kind=='scan':ep.invert_xaxis()
        ep.grid(alpha=.2);ep.set_title(f'Evaluated image {n}: {energies[n]:+.3f} eV',fontsize=11)
        fig.suptitle('MACE rigid approach: NOT a transition-state path' if kind=='scan' else 'MACE: evaluated minimization branches, not IRC' if folder.name=='mace_local_saddle' else ('DPA-3.3 OMol25' if 'omol25' in folder.name else 'MACE-MP-0b2')+' reaction path: see status',fontsize=12)
        fig.text(.5,.025,'Actual evaluated geometries only | no interpolated animation frames',ha='center',fontsize=9)
        fig.subplots_adjust(left=.01,right=.98,bottom=.17,top=.83,wspace=.02)
        b=io.BytesIO();fig.savefig(b,format='png');b.seek(0);frames.append(Image.open(b).convert('RGB'));plt.close(fig)
    frames[0].save(folder/'path.gif',save_all=True,append_images=frames[1:],duration=800,loop=0)
    frames[frame_indices.index(int(energies.argmin() if kind=='scan' else energies.argmax()))].save(folder/'preview.png')
    (folder/'render_manifest.json').write_text(json.dumps(dict(frame_count=len(frames),frame_image_indices=frame_indices,interpolated_frames=0,coordinates_sha256=sha(folder/'images.extxyz'),energies_sha256=sha(folder/'energies.csv'),gif_sha256=sha(folder/'path.gif')),indent=2)+'\n')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--force',action='store_true');args=parser.parse_args()
    surfaces=['Si100','Si111','beta_Si3N4_001','alpha_quartz_001'];species=['HF','HCl','F2','Cl2','H2','H2O','CH3F','SiF4','SiCl4']
    overview=[]
    for surface in surfaces:
        for mol in species:
            parent=ROOT/'data/surface_paths'/surface/mol;root=parent/'mace_site_orientation'
            if not (root/'all_energies.csv').exists():continue
            rows=list(csv.DictReader((root/'all_energies.csv').open()));best=min(rows,key=lambda r:float(r['interaction_energy_eV']))
            folder=root/best['site']/best['orientation'];meta=json.loads((folder/'summary.json').read_text());curves={}
            for r in rows:curves.setdefault(r['site']+'/'+r['orientation'],[]).append(r)
            selected=list(csv.DictReader((folder/'energies.csv').open()))
            if args.force or not (folder/'path.gif').exists():render(folder,surface+' / '+mol+'\n'+best['site']+' / '+best['orientation'],'scan',selected,read(folder/'images.extxyz',':'),curves)
            fig,ax=plt.subplots(figsize=(7,4.3),layout='constrained')
            for label,rr in curves.items():ax.plot([float(r['height_A']) for r in rr],[float(r['interaction_energy_eV']) for r in rr],'o-',ms=3,label=label)
            ax.invert_xaxis();ax.grid(alpha=.2);ax.set(xlabel='Approach height (angstrom)',ylabel='Interaction energy (eV)',title=surface+' / '+mol+': evaluated sites and orientations');ax.legend(fontsize=7,ncol=2)
            fig.savefig(root/'site_comparison.png',dpi=130);plt.close(fig)
            rel=folder.relative_to(parent).as_posix()
            (parent/'README.md').write_text(f'# {surface} / {mol}\n\n![3D evaluated approach]({rel}/path.gif)\n\n![All sites and orientations](mace_site_orientation/site_comparison.png)\n\nThese are rigid approach scans, not NEB paths or transition states. Each GIF frame has its own evaluated energy. The displayed curve has the lowest sampled interaction energy among the tested starts; it is not a globally optimized adsorption energy.\n\n- Selected site/orientation: **{best["site"]} / {best["orientation"]}**.\n- Lowest sampled interaction energy: **{float(best["interaction_energy_eV"]):.5f} eV** at {best["height_A"]} angstrom.\n- Tested {len(curves)} site/orientation combinations and {len(rows)} geometries.\n- [All energies](mace_site_orientation/all_energies.csv), [selected coordinates]({rel}/images.extxyz), [selected provenance]({rel}/summary.json).\n- [Shared assumptions, equations and sources](../../../../docs/dry_etch_results/MOLECULAR_SURFACES.md).\n',encoding='utf8')
            if (parent/'relaxed_adsorption/summary.json').exists():
                with (parent/'README.md').open('a',encoding='utf8') as f:f.write('\n## Force-relaxed follow-up\n\n[Two starts with the upper substrate mobile and reference-state audit](../../../../docs/dry_etch_results/RELAXED_ADSORPTION.md).\n')
            overview.append(dict(surface=surface,species=mol,minimum_interaction_energy_eV=float(best['interaction_energy_eV']),site=best['site'],orientation=best['orientation'],height_A=float(best['height_A']),evaluations=len(rows),combinations=len(curves),folder=folder.relative_to(ROOT).as_posix()))
            print('Rendered',surface,mol,flush=True)
    dest=ROOT/'docs/dry_etch_results';(dest/'molecular_screening.json').write_text(json.dumps(overview,indent=2)+'\n')
    if len(overview)==36:
        data=np.array([[next(r['minimum_interaction_energy_eV'] for r in overview if r['surface']==s and r['species']==m) for s in surfaces] for m in species])
        fig,ax=plt.subplots(figsize=(8,6),layout='constrained');im=ax.imshow(data,cmap='viridis');fig.colorbar(im,ax=ax,label='Lowest sampled interaction energy (eV)')
        ax.set_xticks(range(4),['Si(100)','Si(111)','beta-Si3N4(001)','quartz(001)'],rotation=15);ax.set_yticks(range(9),species)
        for i in range(9):
            for j in range(4):ax.text(j,i,f'{data[i,j]:.2f}',ha='center',va='center',color='white' if data[i,j]<data.mean() else 'black')
        ax.set_title('Site/orientation screening: rigid scans, not barriers');fig.savefig(dest/'molecular_screening.png',dpi=170);plt.close(fig)
    for method in ['molecular_neb','molecular_refined','mace_local_saddle','omol25_neb','omol25_refined']:
        for folder in (ROOT/'data/surface_paths').glob('*/*/'+method):
            if not (folder/'energies.csv').exists():continue
            meta=json.loads((folder/'summary.json').read_text());rows=list(csv.DictReader((folder/'energies.csv').open()))
            render(folder,folder.parent.parent.name+' / '+folder.parent.name+'\n'+meta['status'],'neb',rows,read(folder/'images.extxyz',':'))
if __name__=='__main__':main()
