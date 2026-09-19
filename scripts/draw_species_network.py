"""Readable family panels generated from every enabled kMC event and state."""
from pathlib import Path
import json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch,Patch
ROOT=Path(__file__).resolve().parents[1]
COLORS=['#B5BEC9','#80CBA4','#2BA6A4','#A489CF']
HF='#F4C95D';REMOVED='#293E57';MISSING='#BE5078'
NAMES={'bridge_NH':'NH bridges','terminal_NH2':'Terminal NH2 ligands','bridge_N':'Bare N bridges','Si_H_rich':'Si-H-rich motif','pre_SiH2F':'Pre-existing SiH2F','pre_SiHF2':'Pre-existing SiHF2','Si_Si_backbond':'Si-Si backbond'}

def edge(ax,a,b,color,rad=0,style='-',shrink=19):
    ax.add_patch(FancyArrowPatch(a,b,connectionstyle=f'arc3,rad={rad}',arrowstyle='-|>',mutation_scale=12,lw=1.6,color=color,linestyle=style,shrinkA=shrink,shrinkB=shrink,zorder=1))
def save(fig,path):
    fig.savefig(path.with_suffix('.png'),dpi=170,facecolor='white')
    fig.savefig(path.with_suffix('.svg'),facecolor='white')
    p=path.with_suffix('.svg');p.write_text('\n'.join(x.rstrip() for x in p.read_text().splitlines())+'\n')
    plt.close(fig)
def main():
    source=ROOT/'configs/species_kmc_network.json';n=json.loads(source.read_text());out=ROOT/'docs/reaction_network';out.mkdir(exist_ok=True)
    states=n['states'];families=list(dict.fromkeys(s['family'] for s in states));reachable={e['target'] for e in n['events']};pos={}
    fig,axes=plt.subplots(len(families),1,figsize=(16,17));fig.subplots_adjust(left=.055,right=.975,top=.90,bottom=.055,hspace=.65)
    for row,(ax,family) in enumerate(zip(axes,families)):
        ids=[i for i,s in enumerate(states) if s['family']==family];events=[e for e in n['events'] if e['source'] in ids]
        stages=sorted(set(states[i]['fluorination_stage'] for i in ids if not states[i]['HF_complex']))
        for i in ids:
            s=states[i];stage=stages.index(s['fluorination_stage']);x=stage*2.1+(.78 if s['HF_complex'] else 0)
            pos[i]=(x,-.22 if s['HF_complex'] else .50)
        for e in events:
            a,b=pos[e['source']],pos[e['target']]
            if e['kind']=='reaction':
                edge(ax,a,b,'#B45928');gas=', '.join(k for k,v in e['gas_delta'].items() if v>0)
                ax.text((a[0]+b[0])/2,-.79,f"{e['pathway']}  |  {e['barrier_eV']:.2f} eV"+('\n+ '+gas if gas else ''),ha='center',va='top',fontsize=9,color='#843F1E')
            else:edge(ax,a,b,'#7B8794',rad=.23 if e['kind']=='desorption' else .08)
        for i in ids:
            s=states[i];x,y=pos[i];disabled=s['Si_removed'] and i not in reachable
            color='white' if disabled else HF if s['HF_complex'] else REMOVED if s['Si_removed'] else COLORS[min(s['fluorination_stage'],3)]
            label='unresolved\nrelease' if disabled else 'HF complex' if s['HF_complex'] else 'Si released' if s['Si_removed'] else 'F'+str(s['fluorination_stage'])
            ax.text(x,y,label,ha='center',va='center',fontsize=10,color='white' if s['Si_removed'] and not disabled else '#172B40',bbox=dict(boxstyle='round,pad=.48',fc=color,ec=MISSING if disabled else '#65768A',lw=1.3,linestyle='--' if disabled else '-'),zorder=4)
            if disabled:
                previous=next(j for j in ids if states[j]['HF_complex'] and states[j]['fluorination_stage']==3)
                edge(ax,pos[previous],pos[i],MISSING,style='--');ax.text((pos[previous][0]+x)/2,-.8,'No matched barrier\nDISABLED',ha='center',va='top',fontsize=9,color=MISSING)
        if family=='Si_Si_backbond':ax.text(3.4,.45,'Further chemistry unresolved',fontsize=10,color=MISSING,va='center')
        ax.set_xlim(-.6,9.1);ax.set_ylim(-1.4,1.15);ax.axis('off')
        title=NAMES.get(family,family.replace('_',' '));ax.text(-.5,1.04,f'{row+1}  {title}',fontsize=12,weight='bold',color='#172B40')
        ax.text(9,1.04,f'{len(ids)} states / {len(events)} enabled events',ha='right',fontsize=9,color='#65768A')
    fig.suptitle('HF / SiN:H - complete enabled kMC network',fontsize=23,weight='bold',y=.98,color='#172B40')
    fig.text(.5,.952,'Seven independent motif families, not a fully connected reaction graph | 45 states | 55 enabled events | 16 source pathways',ha='center',fontsize=11)
    handles=[Patch(fc=c,label='F'+str(i)) for i,c in enumerate(COLORS)]+[Patch(fc=HF,label='Adsorbed HF complex'),Patch(fc=REMOVED,label='Si released'),Patch(fc='white',ec=MISSING,ls='--',label='Missing / disabled')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.938),ncol=7,frameon=False,fontsize=10)
    fig.text(.5,.028,'Gray arrow pairs: HF adsorption / desorption. Brown arrows: chemical conversion, source Ea and gas products.\nF0-F3 count incorporated F at the central motif; adsorbed HF is additional. Rates remain conditional; no calibrated EPC.',ha='center',fontsize=10)
    save(fig,out/'species_kmc_network_compact')
    # Compact conceptual view is explicitly a repeated motif template, not an extra enabled graph.
    fig,ax=plt.subplots(figsize=(12,3.6));ax.axis('off');ax.set_xlim(-.7,6.7);ax.set_ylim(-1.7,1.4)
    nodes=[((0,0),'Available\nmotif',COLORS[0]),((2,0),'HF adsorption\ncomplex',HF),((4,0),'Converted\nmotif',COLORS[2]),((6,0),'Si released\nif parameterized',REMOVED)]
    for (x,y),label,c in nodes:ax.text(x,y,label,ha='center',va='center',fontsize=12,color='white' if c==REMOVED else '#172B40',bbox=dict(boxstyle='round,pad=.7',fc=c,ec='#65768A'),zorder=3)
    edge(ax,(0,0),(2,0),'#7B8794');edge(ax,(2,-.2),(0,-.2),'#7B8794',rad=-.6)
    edge(ax,(2,0),(4,0),'#B45928');edge(ax,(4,0),(6,0),'#B45928',style='--')
    for x,t in [(1,'+ HF'),(3,'bond conversion'),(5,'further steps')]:ax.text(x,.6,t,ha='center',fontsize=10)
    ax.text(1,-.9,'HF desorption',ha='center',fontsize=10,color='#65768A');ax.text(4.7,-.9,'F1 / F2 / F3 and residual states\nare expanded in the full network',ha='center',fontsize=10)
    fig.suptitle('Read one motif pathway from left to right',fontsize=17,weight='bold');fig.text(.5,.035,'Overview template only. The complete graph below specifies which steps are enabled or missing.',ha='center',fontsize=10)
    save(fig,out/'species_kmc_overview')
    (out/'render_manifest.json').write_text(json.dumps(dict(network_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),state_ids=[s['id'] for s in states],event_ids=[e['id'] for e in n['events']],families=families,state_colors=dict(F0=COLORS[0],F1=COLORS[1],F2=COLORS[2],F3=COLORS[3],HF_complex=HF,Si_released=REMOVED),interpretation='Full graph contains all configured states and events; overview is conceptual'),indent=2)+'\n')
    from draw_surface_states import main as draw_surfaces
    draw_surfaces()
    print('Rendered complete network and overview')
if __name__=='__main__':main()
