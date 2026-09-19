"""Run and verify a conditional species-resolved HF/SiN:H sensitivity campaign."""
from pathlib import Path
import csv,hashlib,io,json,sys,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from plasma_surface.species_kmc import build_network,rate_constants,initial_counts,simulate,expectation
from plasma_surface.species_native import simulate as native

def write_csv(path,rows):
    with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)

def ledger(network,result,sites):
    elements=['Si','N','H','F'];S=np.array([[s['composition'].get(e,0) for e in elements] for s in network['states']]);G=np.array([[v.get(e,0) for e in elements] for v in network['gas_species'].values()]);expected=initial_counts(network,sites)@S
    assert np.array_equal(result['counts'].sum(axis=1),np.full(len(result['counts']),sites))
    assert np.all(result['counts']>=0)
    assert np.all(result['counts']@S+result['gas_counts']@G==expected)

def render_lattice(network,result,out):
    palette=['#aab3bd','#8dcfaf','#3ca49e','#8760b0','#283c51'];states=network['states'];n=result['lattice'].shape[1];side=int(np.sqrt(n));x,y=np.meshgrid(np.arange(side),np.arange(side));frames=[]
    for k in range(len(result['times'])):
        grid=result['lattice'][k];stage=np.array([states[i]['fluorination_stage'] for i in grid]);removed=np.array([states[i]['Si_removed'] for i in grid]);occupied=np.array([states[i]['HF_complex'] for i in grid]);colors=[palette[4] if gone else palette[min(j,3)] for j,gone in zip(stage,removed)]
        fig=plt.figure(figsize=(10,4.8),dpi=100);ax=fig.add_subplot(121);view=fig.add_subplot(122,projection='3d')
        ax.scatter(x.ravel(),y.ravel(),c=colors,s=70,marker='s');ax.scatter(x.ravel()[occupied],y.ravel()[occupied],s=12,c='#f5cb57');ax.set_aspect('equal');ax.set_axis_off();ax.set_title('Named surface motifs; yellow = HF complex')
        z=np.where(removed,0.,1.);view.scatter(x.ravel(),y.ravel(),z,c=colors,s=28,depthshade=True);view.view_init(28,-58);view.set_zlim(0,2);view.set_axis_off();view.set_title('Perspective of motif occupancy')
        fig.suptitle(f'HF/SiN:H conditional kMC | 400 K | t = {result["times"][k]:.2f} s',fontsize=13)
        fig.text(.5,.04,'Gray: F0 | green: F1 | teal: F2 | purple: F3 | dark: Si released\nSchematic reactive-site grid, not an atomistic crystal or physical film thickness',ha='center',fontsize=9)
        fig.subplots_adjust(left=.02,right=.98,bottom=.18,top=.82);b=io.BytesIO();fig.savefig(b,format='png');b.seek(0);frames.append(Image.open(b).convert('RGB'));plt.close(fig)
    frames[0].save(out/'species_kmc.gif',save_all=True,append_images=frames[1:],duration=250,loop=0)
    np.savez_compressed(out/'lattice_trajectory.npz',times=result['times'],states=result['lattice'],counts=result['counts'],gas_counts=result['gas_counts'])
    (out/'animation_manifest.json').write_text(json.dumps(dict(frame_count=len(frames),interpolated_frames=0,representation='Actual kMC site states on a schematic grid; no atomistic geometry or length scale',seed=741,trajectory_sha256=hashlib.sha256((out/'lattice_trajectory.npz').read_bytes()).hexdigest(),gif_sha256=hashlib.sha256((out/'species_kmc.gif').read_bytes()).hexdigest()),indent=2)+'\n')

def main():
    out=ROOT/'docs/species_kmc_results';out.mkdir(exist_ok=True);network=build_network();(ROOT/'configs/species_kmc_network.json').write_text(json.dumps(network,indent=2)+'\n')
    sites=1000;replicates=128;times=np.linspace(0,3,61);temperatures=[325.,350.,400.,450.];removed=np.array([s['Si_removed'] for s in network['states']]);gases=list(network['gas_species']);summaries=[];records=[];checks=[];rate_rows=[];saved={}
    fig,axs=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
    for temp in temperatures:
        rates=rate_constants(network,temperature=temp);samples=[native(network,rates,sites,times,seed=10000+i) for i in range(replicates)]
        for r in samples:ledger(network,r,sites)
        counts=np.array([r['counts'] for r in samples]);gas=np.array([r['gas_counts'] for r in samples]);mean=counts.mean(axis=0);gmean=gas.mean(axis=0);exact,egas=expectation(network,rates,sites,times)
        assert np.max(abs(exact.sum(axis=1)-sites))<1e-4
        se=counts.std(axis=0,ddof=1)/np.sqrt(replicates);mask=(se>.05)&(exact>2);z=float(np.max(abs(mean-exact)[mask]/se[mask]))
        assert z<6
        checks.append(dict(temperature_K=temp,max_standard_errors=z,conservation_exact=True,master_equation_count_error=float(np.max(abs(mean-exact)))))
        converted=mean[:,removed].sum(axis=1)/sites;axs[0].plot(times,converted,label=f'{temp:g} K');axs[0].plot(times,exact[:,removed].sum(axis=1)/sites,'k--',lw=.65,alpha=.4)
        for t,c,g in zip(times,mean,gmean):
            for i,s in enumerate(network['states']):records.append(dict(temperature_K=temp,time_s=float(t),species=s['id'],phase='surface',mean_count=float(c[i]),fraction_per_initial_motif=float(c[i]/sites)))
            for i,name in enumerate(gases):records.append(dict(temperature_K=temp,time_s=float(t),species=name,phase='reservoir_net' if name=='HF' else 'gas_product',mean_count=float(g[i]),fraction_per_initial_motif=float(g[i]/sites)))
        for e,k in zip(network['events'],rates):rate_rows.append(dict(temperature_K=temp,event=e['id'],kind=e['kind'],rate_s=float(k),source_pathway=e['pathway'] or '',barrier_eV=e['barrier_eV'],rate_status=e['rate_status']))
        summaries.append(dict(temperature_K=temp,Si_released_fraction=float(converted[-1]),products_per_initial_motif={name:float(gmean[-1,i]/sites) for i,name in enumerate(gases)},mean_events=float(np.mean([r['events'] for r in samples]))))
        saved[temp]=(mean,gmean,exact,egas)
    for gas in ['NH3','H2','SiF4','SiH2F2','SiHF3']:axs[1].plot(times,saved[400.][1][:,gases.index(gas)]/sites,label=gas)
    axs[0].set(xlabel='Time (s)',ylabel='Fraction of initial motifs releasing Si',title='Conditional HF exposure / purge; not EPC',ylim=(0,.55));axs[1].set(xlabel='Time (s)',ylabel='Products per initial reactive motif',title='Species-resolved products at 400 K')
    for ax in axs:ax.axvline(2,color='gray',ls=':',label='HF off / purge');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Literature Ea sensitivity with assumed arrival/desorption/prefactors',fontsize=12);fig.savefig(out/'species_kinetics.png',dpi=170);plt.close(fig)
    mean=saved[400.][0];families=list(dict.fromkeys(s['family'] for s in network['states']));fig,axes=plt.subplots(2,4,figsize=(14,6),layout='constrained')
    for ax,family in zip(axes.flat,families):
        for j in range(5):
            indices=[i for i,s in enumerate(network['states']) if s['family']==family and s['fluorination_stage']==j and not s['Si_removed']]
            if indices:ax.plot(times,mean[:,indices].sum(axis=1)/sites,label=f'F{j} (incl. HF complex)')
        indices=[i for i,s in enumerate(network['states']) if s['family']==family and s['Si_removed']]
        if indices:ax.plot(times,mean[:,indices].sum(axis=1)/sites,'k--',label='Si released')
        ax.set(title=family,xlabel='Time (s)',ylabel='Fraction of all initial motifs');ax.grid(alpha=.2);ax.legend(fontsize=6)
    axes.flat[-1].set_axis_off();fig.suptitle('Intermediate populations at 400 K: missing pathways remain blocked');fig.savefig(out/'intermediate_populations.png',dpi=160);plt.close(fig)
    # Independent Python SSA and native implementation agree statistically, not bitwise.
    rates=rate_constants(network);py=[simulate(network,rates,500,[0,2,3],seed=30000+i) for i in range(32)];cp=[native(network,rates,500,[0,2,3],seed=50000+i) for i in range(128)]
    for r in py:ledger(network,r,500)
    a=np.array([r['counts'][-1] for r in py]);b=np.array([r['counts'][-1] for r in cp]);se=np.sqrt(a.var(axis=0,ddof=1)/len(a)+b.var(axis=0,ddof=1)/len(b));mask=se>.1;backend_z=float(np.max(abs(a.mean(axis=0)-b.mean(axis=0))[mask]/se[mask]));assert backend_z<6
    benchmark={}
    for name,fn in [('python',simulate),('cpp',native)]:
        trials=[]
        for seed in range(661,666):
            start=time.perf_counter();r=fn(network,rates,1000,[0,2,3],seed=seed);elapsed=time.perf_counter()-start;trials.append(dict(seconds=elapsed,events=r['events']))
        benchmark[name]=dict(repetitions=len(trials),median_seconds=float(np.median([x['seconds'] for x in trials])),events_per_second=sum(x['events'] for x in trials)/sum(x['seconds'] for x in trials),trials=trials)
    sensitivity=[]
    for nu in [1e11,1e12,1e13]:
        for des in [10.,100.,1000.]:
            rr=rate_constants(network,prefactor=nu,desorption=des);cs,gs=expectation(network,rr,sites,times)
            sensitivity.append(dict(assumed_prefactor_s=nu,assumed_desorption_s=des,assumed_arrival_s=5.,temperature_K=400.,final_Si_released_fraction=float(cs[-1,removed].sum()/sites)))
    write_csv(out/'species_populations.csv',records);write_csv(out/'event_rates.csv',rate_rows);write_csv(out/'assumption_sensitivity.csv',sensitivity)
    write_csv(out/'event_counts_400K.csv',[dict(event=e['id'],count=int(c)) for e,c in zip(network['events'],native(network,rates,sites,times,seed=741)['event_counts'])])
    result=dict(status='conditional_sensitivity_not_validated_etch_prediction',states=len(network['states']),enabled_events=len(network['events']),distinct_source_pathways=16,sites=sites,replicates=replicates,exposure_s=2.,purge_s=1.,assumed_HF_arrival_s=5.,assumed_complex_desorption_s=100.,assumed_prefactor_s=1e12,initial_mixture=network['initial_motif_fractions'],temperatures=summaries,verification=checks,python_cpp_max_standard_errors=backend_z,benchmark=benchmark,seed_policy='C++ ensembles 10000..10127; Python check 30000..30031; independent C++ check 50000..50127; lattice 741',source_sha256=hashlib.sha256((ROOT/'data/literature/sin_hf_pathways.csv').read_bytes()).hexdigest(),code_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['plasma_surface/species_kmc.py','plasma_surface/species_native.py','cpp/species.cpp','scripts/run_species_kmc.py']},assumptions=network['assumptions'])
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');(out/'cpp_build.json').write_bytes((ROOT/'build/species_build_manifest.json').read_bytes())
    lattice=simulate(network,rates,400,np.linspace(0,3,41),seed=741,track_lattice=True);ledger(network,lattice,400);render_lattice(network,lattice,out)
    print(json.dumps(dict(verification=checks,backend_z=backend_z,benchmark=benchmark,results=summaries)),flush=True)
if __name__=='__main__':main()
