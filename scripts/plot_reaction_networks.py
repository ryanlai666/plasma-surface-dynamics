"""Plot candidate HiPRGen ladders and the actual enabled/disabled kMC network."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
ROOT=Path(__file__).resolve().parents[1]

def arrow(ax,a,b,color='#52687b',style='-',rad=0):
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=9,linewidth=1.1,color=color,linestyle=style,connectionstyle=f'arc3,rad={rad}',shrinkA=16,shrinkB=16))

def main():
    d=json.loads((ROOT/'data/reaction_network/hiprgen/network.json').read_text());out=ROOT/'docs/reaction_network';out.mkdir(exist_ok=True)
    fig,axes=plt.subplots(2,2,figsize=(15,10),layout='constrained')
    for ax,family in zip(axes.flat,[f['family'] for f in d['families']]):
        states=[s for s in d['species'] if s['family']==family and s['phase']=='capped_molecular_proxy'];positions={s['id']:(s['halogens'],-s['cap_H']) for s in states}
        for r in d['reactions']:
            if r['family']!=family or r['direction']!='forward_substitution':continue
            a=next(x for x in r['reactants'] if x in positions);b=next(x for x in r['products'] if x in positions);arrow(ax,positions[a],positions[b])
        for s in states:
            x,y=positions[s['id']];ax.scatter(x,y,s=1100,marker='s',c='#d5e9ef' if s['bound_ligands'] else '#dfe7c5',edgecolor='#455e70',zorder=3)
            ax.text(x,y,s['id'].replace('(','\n('),ha='center',va='center',fontsize=7,zorder=4)
        ax.set(xlim=(-.6,4.6),ylim=(-4.6,.7),xticks=range(5),yticks=[0,-1,-2,-3,-4],yticklabels=[0,1,2,3,4],xlabel='Halogen ligands on central Si',ylabel='Fixed H caps',title=family+': + H'+family.split('_')[1]+' / coproduct '+('NH3' if family.startswith('NH') else 'H2O'))
        ax.grid(alpha=.1)
    fig.suptitle('HiPRGen function-level pilot: 40 forward candidate substitutions + 40 reverse\nSupplied capped-motif library; no energy screening, surface stability or TS verification',fontsize=13)
    fig.savefig(out/'hiprgen_candidates.png',dpi=170);fig.savefig(out/'hiprgen_candidates.svg');plt.close(fig)
    n=json.loads((ROOT/'configs/species_kmc_network.json').read_text());families=list(dict.fromkeys(s['family'] for s in n['states']));positions={}
    fig,ax=plt.subplots(figsize=(17,10))
    for i,s in enumerate(n['states']):
        row=families.index(s['family']);x=s['fluorination_stage']+(.43 if s['HF_complex'] else 0)
        if s['family']=='Si_Si_backbond':x-=1
        if s['family'].startswith('pre_'):x=0.43 if s['HF_complex'] else 1 if s['Si_removed'] else 0
        positions[i]=(x,-row)
    for e in n['events']:
        a=positions[e['source']];b=positions[e['target']]
        if e['kind']=='reaction':
            arrow(ax,a,b,color='#bd633b');label=e['pathway']+' / '+str(e['barrier_eV'])+' eV';gas=', '.join(k for k,v in e['gas_delta'].items() if v>0)
            ax.text((a[0]+b[0])/2,a[1]+.26,label+('\n'+gas if gas else ''),ha='center',fontsize=7,color='#884321')
        else:arrow(ax,a,b,color='#91a0ac',rad=-.4 if e['kind']=='desorption' else 0)
    for i,s in enumerate(n['states']):
        x,y=positions[i];color='#f7df8e' if s['HF_complex'] else '#ccd9c0' if s['Si_removed'] else '#c8dce9'
        ax.scatter(x,y,s=820,marker='s',c=color,edgecolor='#536c7d',zorder=3)
        label='HF\ncomplex' if s['HF_complex'] else 'residual' if s['Si_removed'] else 'F'+str(s['fluorination_stage'])
        ax.text(x,y,label,ha='center',va='center',fontsize=7,zorder=4)
    for family in ['terminal_NH2','bridge_N']:
        row=families.index(family);arrow(ax,(3.43,-row),(4,-row),color='#9c5b85',style='--');ax.text(3.7,-row+.25,'missing barrier\ndisabled',ha='center',fontsize=7,color='#9c5b85')
    ax.text(1.8,-families.index('Si_Si_backbond'),'Further connectivity/steps unresolved',va='center',fontsize=9,color='#9c5b85')
    ax.set(xlim=(-1.4,4.65),ylim=(-6.65,.65));ax.set_axis_off()
    for row,family in enumerate(families):ax.text(-.75,-row,family.replace('_',' '),ha='center',va='center',fontsize=10)
    fig.suptitle('HF/SiN:H species-resolved kMC: 45 named states, 55 enabled events\nAll 16 source pathways represented; connecting motifs is a conditional topology assumption',fontsize=14)
    fig.text(.5,.045,'Gray arrows: assumed HF arrival (5/s during dose) and complex desorption (100/s) | Brown: source Ea with assumed prefactor 1e12/s\nDashed: unparameterized and disabled | Fixed neighboring anchors are implicit | No calibrated EPC or ion-assisted kinetics',ha='center',fontsize=10)
    fig.subplots_adjust(left=.02,right=.99,bottom=.13,top=.86);fig.savefig(out/'species_kmc_network.png',dpi=170);fig.savefig(out/'species_kmc_network.svg');plt.close(fig)
    print('Saved reaction network PNG/SVG plots')
if __name__=='__main__':main()
