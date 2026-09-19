"""Render saved graph states with separate geometry, exposure and chemistry channels."""
from pathlib import Path
from collections import Counter
import hashlib, io, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from PIL import Image
ROOT = Path(__file__).resolve().parents[1]
COLORS = dict(Si='#546D91', N='#84B6D4', H='#9099A5', F='#179D83', Cl='#A45FB5', HF='#E6A51A', exposed='#D04A34', new='#12A6A3')


def frame_data(d, k):
    """Auditable plot masks; new means first exposed since preceding saved frame."""
    s = d['snapshots'][k]
    active = np.asarray(s['active'], bool)
    ever = set()
    fresh = set()
    previous = d['snapshots'][k-1]['time_s'] if k else -1.
    for e in d['events']:
        if e['time_s'] > s['time_s']:
            break
        for i in e['newly_exposed']:
            if i not in ever and e['time_s'] > previous:
                fresh.add(i)
            ever.add(i)
    return dict(active=active, exposed=active & np.asarray(s['exposed'], bool),
                fresh=np.array([i in fresh for i in range(len(active))]) & active,
                hf=active & np.asarray(s['precursors'], bool),
                first_exposed_ids=sorted(fresh), ever_exposed_ids=sorted(ever),
                term_counts={e:sum(t.get(e,0) for t,a in zip(s['terminations'],active) if a) for e in ('H','F','Cl')})


def main():
    src=ROOT/'data/multilayer/demo_trajectory.json'
    d=json.loads(src.read_text());out=ROOT/'docs/multilayer_results'
    nodes=d['initial']['nodes'];pos=np.array([n['position_A'] for n in nodes])
    columns=np.array([n['column'] for n in nodes]);fixed=np.array([n['fixed'] for n in nodes])
    layers=np.array([n['layer'] for n in nodes]);sizes=Counter(layers)
    # Fixed fractional-y row 2, chosen before inspecting trajectory removal.
    cut=columns//6==2
    colors=np.array([COLORS[n['element']] for n in nodes])
    zlo,zhi=pos[:,2].min()-.8,pos[:,2].max()+.8
    centers=[pos[layers==l,2].mean() for l in range(6)]
    bounds=[zlo]+[(a+b)/2 for a,b in zip(centers,centers[1:])]+[zhi]
    frames=[];records=[];removed_history=[]
    summary=json.loads((out/'summary.json').read_text())
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    for k,snap in enumerate(d['snapshots']):
        f=frame_data(d,k);active=f['active'];t=snap['time_s']
        counts=Counter(layers[~active]);removed=Counter(n['element'] for n,a in zip(nodes,active) if not a)
        fig=plt.figure(figsize=(15,9),dpi=100,facecolor='#F8FAFC')
        grid=fig.add_gridspec(2,3,left=.025,right=.975,bottom=.19,top=.79,width_ratios=[1.25,1,.8],height_ratios=[1,1],wspace=.28,hspace=.75)
        view=fig.add_subplot(grid[:,0],projection='3d');side=fig.add_subplot(grid[:,1]);chem=fig.add_subplot(grid[0,2]);bar=fig.add_subplot(grid[1,2])
        for ax in (view,side,chem,bar):ax.set_facecolor('#F8FAFC')
        local=[(a,b) for a,b in snap['bonds'] if np.linalg.norm(pos[a]-pos[b])<=2.1]
        view.add_collection3d(Line3DCollection([pos[[a,b]] for a,b in local],colors='#ABB8C7',linewidths=.7,alpha=.5))
        def dots(mask,**kw):
            if mask.any():view.scatter(*pos[mask].T,depthshade=False,**kw)
        dots(active & ~fixed,c=colors[active & ~fixed],s=24,alpha=.7)
        dots(active & fixed,c='#CED4DE',s=18,alpha=.6)
        dots(~active,facecolors='none',edgecolors='#B1B7C0',s=32,linewidths=.7)
        dots(f['exposed'],facecolors='none',edgecolors=COLORS['exposed'],s=62,linewidths=1.5)
        dots(f['fresh'],c=COLORS['new'],marker='*',s=120,edgecolors='white',linewidths=.4)
        dots(f['hf'],c=COLORS['HF'],marker='D',s=26)
        view.view_init(17,-62);view.set_box_aspect(np.ptp(pos,axis=0), zoom=1.22);view.set(zlim=(zlo,zhi));view.axis('off')
        view.set_title('01  Full substrate connectivity',loc='left',weight='bold',pad=12)
        view.text2D(.02,.02,'All 336 initial host sites\nGray base: fixed bottom band',transform=view.transAxes,fontsize=10,color='#56657A')
        for l in range(6):
            side.axhspan(bounds[l],bounds[l+1],color='#C3CDDB' if l==0 else '#DFEAF3',alpha=.6 if l%2==0 else .23,zorder=0)
        segments=[pos[[a,b]][:,[0,2]] for a,b in local if cut[a] and cut[b]]
        side.add_collection(LineCollection(segments,colors='#9BAABE',linewidths=1.1,zorder=1))
        def side_dots(mask,**kw):
            m=mask & cut
            side.scatter(pos[m,0],pos[m,2],**kw)
        side_dots(active,c=colors[active & cut],s=46,edgecolors='white',linewidths=.5,zorder=3)
        side_dots(~active,facecolors='none',edgecolors='#9FA9B7',s=55,linewidths=1,zorder=2)
        side_dots(f['exposed'],facecolors='none',edgecolors=COLORS['exposed'],s=100,linewidths=1.7,zorder=4)
        side_dots(f['fresh'],c=COLORS['new'],marker='*',s=160,zorder=6)
        side_dots(f['hf'],c=COLORS['HF'],marker='D',s=33,zorder=5)
        side.set(ylim=(zlo,zhi),xlim=(pos[cut,0].min()-1,pos[cut,0].max()+1),xlabel='x (angstrom)',ylabel='z (angstrom)')
        side.set_title('02  Fixed thin cross-section',loc='left',weight='bold',pad=12)
        side.text(0,1.015,'Fractional y in [2/6, 3/6); same atoms each frame',transform=side.transAxes,fontsize=8,color='#56657A')
        for l,z in enumerate(centers):side.text(1.015,z,f'B{l}',transform=side.get_yaxis_transform(),fontsize=9,color='#56657A',va='center')
        labels=['H caps','F caps','Cl caps','HF complexes'];values=[f['term_counts'][e] for e in ('H','F','Cl')]+[int(f['hf'].sum())]
        chem.barh(labels,values,color=[COLORS[e] for e in ('H','F','Cl','HF')],height=.58)
        chem.set_xlim(0,160);chem.invert_yaxis();chem.set_xlabel('Attached species count')
        chem.set_title('03  Surface chemistry',loc='left',weight='bold',pad=12)
        for j,v in enumerate(values):chem.text(v+2,j,str(v),va='center',fontsize=10)
        chem.text(0,-.36,'H/F/Cl counted separately, including mixed caps.\nHF occupancy is additional to these terminations.',transform=chem.transAxes,fontsize=8,color='#56657A')
        ls=list(range(5,-1,-1));bar.barh(range(6),[counts[l] for l in ls],color='#349D9C',height=.65)
        bar.set(yticks=range(6),yticklabels=[f'B{l}'+(' (fixed)' if l==0 else '') for l in ls],xlim=(0,56),xlabel='Removed / 56 initial hosts per band');bar.invert_yaxis()
        bar.set_title('04  Removal by initial depth',loc='left',weight='bold',pad=12)
        for j,l in enumerate(ls):bar.text(counts[l]+1,j,str(counts[l]),va='center',fontsize=9)
        phase='Complete' if t>=9 else ('HF dose' if t%3<2 else 'Purge')
        fig.text(.035,.952,'Multilayer plasma-surface bond graph',fontsize=22,weight='bold',color='#24354C')
        fig.text(.035,.914,f'{t:4.2f} s  |  {phase}  |  450 K  |  Unvalidated rate demonstration',fontsize=12,color='#56657A')
        fig.text(.965,.946,f"Removed: {removed['Si']} Si + {removed['N']} N",ha='right',fontsize=15,weight='bold',color='#24354C')
        fig.text(.965,.914,f"Newly exposed so far: {len(f['ever_exposed_ids'])}  |  Active bonds: {len(snap['bonds'])}",ha='right',fontsize=11,color='#56657A')
        timeline=fig.add_axes([.035,.865,.93,.018])
        for start in (0,3,6):
            timeline.axvspan(start,start+2,color='#DDEDE8');timeline.axvspan(start+2,start+3,color='#E5E8ED')
        timeline.axvline(t,color='#24354C',lw=2);timeline.set(xlim=(0,9),yticks=[],xticks=range(10));timeline.tick_params(labelsize=8,length=2)
        handles=[Line2D([],[],ls='',marker=m,markersize=8,markerfacecolor=fc,markeredgecolor=ec,label=label) for m,fc,ec,label in [('o',COLORS['Si'],COLORS['Si'],'Si host'),('o',COLORS['N'],COLORS['N'],'N host'),('o','none',COLORS['exposed'],'Currently exposed'),('*',COLORS['new'],COLORS['new'],'First exposed this interval'),('D',COLORS['HF'],COLORS['HF'],'Adsorbed HF'),('o','none','#9FA9B7','Removed host site')]]
        fig.legend(handles=handles,ncol=3,loc='lower center',bbox_to_anchor=(.5,.085),frameon=False,columnspacing=3,fontsize=10)
        fig.text(.5,.049,'Saved kMC states only; host coordinates fixed. B0-B5 are unit-cell depth bands, not atomic monolayers.\nPeriodic wrap bonds are omitted visually; the model retains them. Cap counts are not relaxed ligand coordinates.',ha='center',fontsize=9,color='#56657A')
        b=io.BytesIO();fig.savefig(b,format='png');b.seek(0);frames.append(Image.open(b).convert('RGB'));plt.close(fig)
        records.append(dict(snapshot_index=k,time_s=t,removed=dict(removed),first_exposed_ids=f['first_exposed_ids'],ever_exposed_ids=f['ever_exposed_ids'],termination_counts=f['term_counts'],HF_occupancy=int(f['hf'].sum()),drawn_bonds=len(local),periodic_bonds_omitted=len(snap['bonds'])-len(local)))
        removed_history.append([t,int((~active).sum())])
    # Reserve semantic colors explicitly: automatic GIF quantization can merge
    # small HF diamonds into the much more common red exposure outlines.
    from matplotlib.colors import to_rgb
    palette=frames[0].quantize(colors=240).getpalette()[:720]
    for color in COLORS.values():
        palette.extend(round(255*v) for v in to_rgb(color))
    palette=(palette+[255]*768)[:768]
    palette_image=Image.new('P',(1,1));palette_image.putpalette(palette)
    gif_frames=[im.quantize(palette=palette_image,dither=Image.Dither.NONE) for im in frames]
    gif_frames[0].save(out/'multilayer_kmc.gif',save_all=True,append_images=gif_frames[1:],duration=[350]*36+[1800],loop=0,optimize=False)
    frames[-1].save(out/'multilayer_preview.png')
    sheet=Image.new('RGB',(1500,900),'white')
    for j,k in enumerate((0,4,16,36)):
        thumb=frames[k].resize((750,450));sheet.paste(thumb,((j%2)*750,(j//2)*450))
    sheet.save(out/'multilayer_storyboard.png')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    (out/'animation_manifest.json').write_text(json.dumps(dict(frame_count=len(frames),snapshot_indices=list(range(len(frames))),interpolated_frames=0,trajectory_sha256=sha(src),gif_sha256=sha(out/'multilayer_kmc.gif'),renderer_sha256=sha(Path(__file__)),cross_section=dict(fractional_y_interval=[2/6,3/6],node_ids=np.flatnonzero(cut).tolist()),new_exposure_definition='First appearance in event newly_exposed before this snapshot, after preceding snapshot; stars only for hosts still active',frames=records),indent=2)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained');h=np.array(removed_history);axes[0].step(h[:,0],h[:,1],where='post');axes[0].set(xlabel='Time (s)',ylabel='Removed Si + N atoms',title='Three dose/purge cycles; fixed finite substrate')
    for t in [2,5,8]:axes[0].axvline(t,c='gray',ls=':',alpha=.6)
    for dep in [0,1,2,4]:
        rr=[sum(r['removed'].values()) for r in summary['control_runs'] if r['penetration_A']==dep];axes[1].scatter([dep]*len(rr),rr,color='#2A9C9D')
    axes[1].set(xlabel='Assumed access depth below moving envelope (A)',ylabel='Removed Si + N atoms at 9 s',title='Accessibility sensitivity: three seeds per depth');fig.savefig(out/'multilayer_controls.png',dpi=160);plt.close(fig)
    print('Rendered', len(frames), 'saved states')
if __name__=='__main__': main()
