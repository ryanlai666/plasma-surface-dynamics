"""Render actual multilayer graph snapshots, not interpolated atomic dynamics."""
from pathlib import Path
from collections import Counter
import hashlib,io,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]
COLORS=dict(Si='#95A5B6',N='#547ECC',H='#C7C9CB',F='#39AE91',Cl='#C681CE',HF='#F5C64C',fixed='#D9E0E6')

def main():
    src=ROOT/'data/multilayer/demo_trajectory.json';d=json.loads(src.read_text());out=ROOT/'docs/multilayer_results';nodes=d['initial']['nodes'];pos=np.array([n['position_A'] for n in nodes]);frames=[];layer_sizes=Counter(n['layer'] for n in nodes);summary=json.loads((out/'summary.json').read_text());removed_history=[]
    for snap in d['snapshots']:
        active=np.array(snap['active']);colors=[]
        for i,n in enumerate(nodes):
            term=snap['terminations'][i];key='fixed' if n['fixed'] else 'HF' if snap['precursors'][i] else 'Cl' if term.get('Cl',0) else 'F' if term.get('F',0) else 'H' if term.get('H',0) and n['element']=='Si' else n['element'];colors.append(COLORS[key])
        fig=plt.figure(figsize=(12.4,5.8),dpi=100);view=fig.add_subplot(131,projection='3d');side=fig.add_subplot(132);bar=fig.add_subplot(133)
        for a,b in snap['bonds']:
            if np.linalg.norm(pos[a]-pos[b])>2.1:continue
            view.plot(*pos[[a,b]].T,c='#8897A5',lw=.55,alpha=.45);side.plot(pos[[a,b],0],pos[[a,b],2],c='#CBD2D9',lw=.55,zorder=1)
        view.scatter(*pos[active].T,c=np.array(colors)[active],s=21,depthshade=False);view.view_init(20,-55);view.set(xlim=(pos[:,0].min()-1,pos[:,0].max()+1),ylim=(pos[:,1].min()-1,pos[:,1].max()+1),zlim=(pos[:,2].min()-1,pos[:,2].max()+1));view.set_box_aspect(np.ptp(pos,axis=0));view.axis('off');view.set_title('Remaining Si/N bond graph',fontsize=11)
        side.scatter(pos[active,0],pos[active,2],c=np.array(colors)[active],s=25,edgecolors='#607080',linewidths=.2,zorder=3)
        side.set(xlabel='x (angstrom)',ylabel='z (angstrom)',ylim=(pos[:,2].min()-1,pos[:,2].max()+1),title='Projection through all y positions');side.grid(alpha=.15)
        counts=Counter(n['layer'] for n,a in zip(nodes,active) if not a);layers=sorted(layer_sizes);bar.barh(layers,[counts[l]/layer_sizes[l] for l in layers],color='#289C9B');bar.set(xlim=(0,1),xlabel='Removed substrate-atom fraction',ylabel='Initial unit-cell depth band',yticks=layers,title='Depth-resolved removal');bar.text(.02,.98,'0 = fixed bottom; 5 = initial top',transform=bar.transAxes,va='top',fontsize=8)
        phase='3 cycles complete' if snap['time_s']==9 else 'HF dose' if snap['time_s']%3<2 else 'purge';fig.suptitle(f"Multilayer graph demonstration | 450 K | {snap['time_s']:.2f} s | {phase}",fontsize=15,weight='bold')
        handles=[Patch(fc=COLORS[k],label=label) for k,label in [('Si','Si'),('N','N'),('H','Si-H'),('F','F-bearing'),('Cl','Cl-bearing'),('HF','HF complex'),('fixed','fixed base')]]
        fig.legend(handles=handles,ncol=7,loc='lower center',bbox_to_anchor=(.5,.07),frameon=False,fontsize=9)
        fig.text(.5,.027,'Actual event snapshots; crystal host positions fixed. Color priority: HF > Cl > F > Si-H. Ligand counts are explicit in saved states.\nLocal rates and penetration are assumptions; this is not calibrated ALE, a relaxed atomistic movie, or a measured etch depth.',ha='center',fontsize=8.5)
        fig.subplots_adjust(left=.025,right=.98,bottom=.23,top=.84,wspace=.34);b=io.BytesIO();fig.savefig(b,format='png');b.seek(0);frames.append(Image.open(b).convert('RGB'));plt.close(fig);removed_history.append([snap['time_s'],sum(counts.values())])
    frames[0].save(out/'multilayer_kmc.gif',save_all=True,append_images=frames[1:],duration=260,loop=0);frames[-1].save(out/'multilayer_preview.png')
    (out/'animation_manifest.json').write_text(json.dumps(dict(frame_count=len(frames),snapshot_indices=list(range(len(frames))),interpolated_frames=0,trajectory_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),gif_sha256=hashlib.sha256((out/'multilayer_kmc.gif').read_bytes()).hexdigest(),renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained');h=np.array(removed_history);axes[0].step(h[:,0],h[:,1],where='post');axes[0].set(xlabel='Time (s)',ylabel='Removed Si + N atoms',title='Three dose/purge cycles; fixed finite substrate')
    for t in [2,5,8]:axes[0].axvline(t,c='gray',ls=':',alpha=.6)
    for dep in [0,1,2,4]:
        rr=[sum(r['removed'].values()) for r in summary['control_runs'] if r['penetration_A']==dep];axes[1].scatter([dep]*len(rr),rr,color='#2A9C9D')
    axes[1].set(xlabel='Assumed access depth below moving envelope (A)',ylabel='Removed Si + N atoms at 9 s',title='Accessibility sensitivity: three seeds per depth');fig.savefig(out/'multilayer_controls.png',dpi=160);plt.close(fig)
    print('Rendered',len(frames),'actual multilayer snapshots')
if __name__=='__main__':main()
