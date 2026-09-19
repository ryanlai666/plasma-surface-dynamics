"""Render IS -> TS -> FS from published surface stationary-point coordinates.

Only the three published stationary points are shown; no invented intermediate frames.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from ase.io import read
from matplotlib import font_manager
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/surface_paths/Si100_c4x2/SiCl4/published'
OUT=ROOT/'data/surface_paths/Si100_c4x2/SiCl4/published'
FONT=font_manager.findfont('DejaVu Sans')
def font(n):return ImageFont.truetype(FONT,n)

def render(xyz,symbols,progress,labels,energies,title):
    im=Image.new('RGB',(1000,620),'#f4f7fb');d=ImageDraw.Draw(im)
    d.text((25,16),title,fill='#15283b',font=font(23))
    d.text((25,49),'Published Si(100)-c(4x2) structures | reverse source pathway',fill='#44556b',font=font(15))
    # Fixed camera: perspective projection with equal spatial units.
    center=np.array([7.5,7.5,6.0]);v=xyz-center
    az=np.deg2rad(-36);el=np.deg2rad(24)
    right=np.array([np.cos(az),np.sin(az),0]);up=np.array([-np.sin(az)*np.sin(el),np.cos(az)*np.sin(el),np.cos(el)])
    depth=np.cross(right,up);z=v@depth;factor=32/(32-z)
    px=340+22*(v@right)*factor;py=295-22*(v@up)*factor
    radii={'Si':1.11,'Cl':1.02,'H':.31}
    colors={'Si':'#8ba5ba','Cl':'#44b989','H':'#f2f3f5'}
    # Pairs are visual covalent-distance guides, not computed bond orders.
    bonds=[]
    for i in range(len(xyz)):
        for j in range(i):
            dist=np.linalg.norm(xyz[i]-xyz[j])
            if .5<dist<1.18*(radii[symbols[i]]+radii[symbols[j]]):bonds.append(((z[i]+z[j])/2,i,j))
    items=[(float(z[i]),'atom',i) for i in range(len(xyz))]+[(float(a),'bond',(i,j)) for a,i,j in bonds]
    for _,kind,k in sorted(items,key=lambda x:x[0]):
        if kind=='bond':
            i,j=k;d.line((px[i],py[i],px[j],py[j]),fill='#657b8c',width=3)
        else:
            i=k;r=(8 if symbols[i]=='Si' else 9 if symbols[i]=='Cl' else 4)*factor[i]
            col='#e7ad49' if i==112 else colors[symbols[i]]
            d.ellipse((px[i]-r,py[i]-r,px[i]+r,py[i]+r),fill=col,outline='#344c60',width=1)
            d.ellipse((px[i]-.5*r,py[i]-.6*r,px[i]+.1*r,py[i]),fill='#ffffff')
    d.rounded_rectangle((682,100,978,452),radius=12,fill='white',outline='#c8d3df')
    d.text((700,115),'Published electronic energies',fill='#24384e',font=font(16))
    es=np.array(energies)-energies[0];low=min(es)-.15;high=max(es)+.3
    xp=np.array([720,827,940]);yp=375-(es-low)/(high-low)*195
    d.line(list(zip(xp,yp)),fill='#9aaabc',width=3)
    for i in range(3):
        d.ellipse((xp[i]-5,yp[i]-5,xp[i]+5,yp[i]+5),fill='#2d719b')
        d.text((xp[i]-24,yp[i]-25),f'{es[i]:.2f}',fill='#213d52',font=font(15))
        d.text((xp[i]-22,392),['IS','TS','FS'][i],fill='#213d52',font=font(16))
        d.text((xp[i]-24,417),labels[i],fill='#596f80',font=font(13))
    d.text((702,148),'eV relative to IS; lines are guides',fill='#52677b',font=font(13))
    if progress in (0,1,2):
        j=int(progress);d.ellipse((xp[j]-8,yp[j]-8,xp[j]+8,yp[j]+8),fill='#d55848')
        state=['IS','TS','FS'][j]+' / '+labels[j]+' (source geometry)'
    else:raise ValueError('Only original stationary points may be rendered')
    d.text((30,514),state,fill='#203c53',font=font(20))
    d.text((30,547),'Si: blue | reacting Si: gold | Cl: green | H: white',fill='#415c72',font=font(15))
    d.text((30,573),'Three published stationary points only; no connecting trajectory has been computed.',fill='#824b36',font=font(15))
    return im

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    spec=json.loads((ROOT/'data/literature/sicl4_surface_paths.json').read_text());manifest=[]
    for path in spec['paths']:
        labels=path['states'][::-1]
        atoms=[read(DATA/(x+'.extxyz')) for x in labels]
        assert all(a.get_chemical_symbols()==atoms[0].get_chemical_symbols() for a in atoms)
        symbols=atoms[0].get_chemical_symbols();xyz=[a.positions.copy() for a in atoms]
        # Remove only numerical rigid translation measured on the 32 bottom H atoms.
        bottom=np.array([s=='H' for s in symbols]);translations=[]
        for k in range(3):
            shift=(xyz[0][bottom]-xyz[k][bottom]).mean(axis=0);xyz[k]+=shift;translations.append(shift.tolist())
        view_shift=np.array([7.5,7.5,.643])-xyz[0][bottom].mean(axis=0)
        xyz=[x+view_shift for x in xyz]
        energies=np.array(path['relative_energies_kcal_mol'][::-1])/spec['conversion_kcal_mol_per_eV']
        title='SiCl4 recombination / '+path['id'] if path['id']!='IR_flip' else 'Si surface dimer flip / IR'
        frames=[];durations=[];minimum_separation=99.
        for q in (0,1,2):
            pos=xyz[q]
            delta=pos[:,None,:]-pos[None,:,:];dist=np.linalg.norm(delta,axis=-1);np.fill_diagonal(dist,np.inf);minimum_separation=min(minimum_separation,float(dist.min()))
            frames.append(render(pos,symbols,float(q),labels,energies,title));durations.append(1000 if q in (0,1,2) else 80)
        assert minimum_separation>.6,'Source atom overlap'
        name=path['id'].lower();frames[0].save(OUT/(name+'.gif'),save_all=True,append_images=frames[1:],duration=durations,loop=0,optimize=False)
        if path['id']=='OD':frames[1].save(OUT/'od_ts_preview.png')
        manifest.append(dict(path=path['id'],direction='reverse of source adsorption/flip',states=labels,natoms=len(symbols),frames=len(frames),barrier_eV=float(energies[1]-energies[0]),rigid_translation_A=translations,common_display_translation_A=view_shift.tolist(),min_pair_separation_A=minimum_separation,coordinate_mapping='Source atom order unchanged; identical symbols and counts verified',interpolation='None: three original stationary points only',file_sha256=hashlib.sha256((OUT/(name+'.gif')).read_bytes()).hexdigest()))
    (OUT/'manifest.json').write_text(json.dumps(dict(source='https://doi.org/10.3390/sym15010213',authors='Zhang, Zhu and Li',license='CC BY 4.0; adapted visualization',paths=manifest),indent=2)+'\n')
    print('Rendered',len(manifest),'source-based IS-TS-FS surface GIFs')
if __name__=='__main__':main()
