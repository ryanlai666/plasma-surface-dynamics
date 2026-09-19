"""Species-resolved conditional HF/SiN:H kinetics with explicit provenance.

This is a sensitivity model, not calibrated ALE or a film-removal prediction.
"""
from pathlib import Path
import csv,math
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
KB=8.617333262145e-5
GASES={'HF':{'H':1,'F':1},'H2':{'H':2},'NH3':{'N':1,'H':3},'SiF4':{'Si':1,'F':4},'SiH2F2':{'Si':1,'H':2,'F':2},'SiHF3':{'Si':1,'H':1,'F':3}}

def build_network():
    literature={r['pathway']:r for r in csv.DictReader((ROOT/'data/literature/sin_hf_pathways.csv').open())}
    states=[];events=[];lookup={};initial={}
    def state(name,comp,family,stage,occupied=False,removed=False):
        comp={e:int(n) for e,n in comp.items() if n};lookup[name]=len(states);states.append(dict(id=name,composition=comp,family=family,fluorination_stage=stage,HF_complex=occupied,Si_removed=removed));return name
    def event(name,src,dst,kind,gas=None,path=None):
        r=literature[path] if path else None
        events.append(dict(id=name,source=lookup[src],target=lookup[dst],kind=kind,gas_delta=gas or {},pathway=path,barrier_eV=float(r['activation_energy_eV']) if r else None,source_reaction=r['surface_motif_reaction'] if r else None,rate_status='conditional_Ea_and_assumed_prefactor' if r else 'assumed_adsorption_or_desorption_hazard',enabled=True))
    def precursor(name):
        s=states[lookup[name]];c=dict(s['composition']);c['H']=c.get('H',0)+1;c['F']=c.get('F',0)+1
        ads=state(name+'_HF',c,s['family'],s['fluorination_stage'],True)
        event('ads_'+name,name,ads,'adsorption',{'HF':-1});event('des_'+name,ads,name,'desorption',{'HF':1});return ads
    families=[('bridge_NH',.3,['P1c','P2c','P3c','P4']),('terminal_NH2',.15,['P1a','P2a','P3a',None]),('Si_H_rich',.1,['P1b','P2b','P3b','P4']),('bridge_N',.15,['P1d','P2d','P3d',None])]
    for family,fraction,paths in families:
        names=[]
        for j in range(4):
            comp={'Si':1,'N':4,'H':4+j,'F':j} if family=='bridge_NH' else {'Si':1,'N':4-j,'H':8-2*j,'F':j} if family=='terminal_NH2' else {'Si':1,'N':4,'H':j,'F':j} if family=='bridge_N' else {'Si':1,'N':1,'H':4-j,'F':j}
            names.append(state(family+'_F'+str(j),comp,family,j))
        initial[names[0]]=fraction
        residual=state(family+'_residual',{'N':4,'H':8} if family=='bridge_NH' else {'N':1,'H':2} if family=='Si_H_rich' else {'N':4,'H':4} if family=='bridge_N' else {},family,4,removed=True)
        for j,name in enumerate(names):
            ads=precursor(name)
            if paths[j] is None:continue
            coproduct={} if family in ['bridge_NH','bridge_N'] else {'NH3':1} if family=='terminal_NH2' else {'H2':1}
            if j==3:coproduct={'SiF4':1}
            event(paths[j]+'_'+family,ads,names[j+1] if j<3 else residual,'reaction',coproduct,paths[j])
    for family,fraction,path,h,f,product in [('pre_SiH2F',.1,'P2e',2,1,'SiH2F2'),('pre_SiHF2',.1,'P3e',1,2,'SiHF3')]:
        name=state(family,{'Si':1,'N':1,'H':h,'F':f},family,f);initial[name]=fraction;ads=precursor(name);res=state(family+'_residual',{'N':1,'H':1},family,f+1,removed=True);event(path,ads,res,'reaction',{product:1},path)
    name=state('Si_Si_F1',{'Si':2,'F':1},'Si_Si_backbond',1);initial[name]=.1;ads=precursor(name)
    after=state('Si_Si_cleaved_F2_H',{'Si':2,'H':1,'F':2},'Si_Si_backbond',2)
    event('P2f',ads,after,'reaction',{},'P2f')
    # Unparameterized candidate is recorded but never assigned a guessed rate.
    missing=[dict(id='bare_N_bridge_F3_release',status='disabled_missing_matching_final_step'),dict(id='Si_Si_cleavage_followup',status='disabled_missing_connected_followup_states'),dict(id='terminal_F3_cleavage_and_release',source='terminal_NH2_F3_HF',status='disabled_missing_matching_barrier_and_product_state',reason='P4 is a bridging NH event, not an NH2 terminal-ligand event'),dict(id='water_blocking_and_assisted_HF',status='disabled_missing_matched_rates'),dict(id='AFS_NH4F_retention_decomposition',status='disabled_missing_consistent_states_and_rates'),dict(id='ion_sputtering_and_ion_assisted_removal',status='disabled_missing_energy_angle_resolved_yields')]
    network=dict(status='conditional_literature_sensitivity_not_calibrated_ALE',states=states,events=events,gas_species=GASES,initial_motif_fractions=initial,unparameterized_candidates=missing,source='https://doi.org/10.1016/j.apsusc.2024.159414',assumptions=['Connecting stage-specific source motifs into chains is a reduced topology hypothesis; source amorphous environments are not identical','HF arrival and complex-desorption hazards and thermal prefactors are sensitivity assumptions','Ea references have not been resolved against full free-energy profiles; rates are conditional sensitivity only','Untracked neighboring substrate anchors are fixed and unchanged; each state represents one reactive Si-centered motif, not one atom','Initial motif mixture is assumed; no conversion to N/Si composition, thickness or EPC','No new-layer refill, periodic ALE cycle, lateral interactions, diffusion, ions or product readsorption','Terminal NH2 chain stalls at F3 because its final release barrier is missing; no P4 substitution across incompatible motifs'])
    validate(network);return network

def validate(network):
    states=network['states'];ids=set()
    for e in network['events']:
        if e['id'] in ids:raise ValueError('Duplicate event')
        ids.add(e['id']);balance=dict(states[e['target']]['composition'])
        for el,n in states[e['source']]['composition'].items():balance[el]=balance.get(el,0)-n
        for gas,nu in e['gas_delta'].items():
            for el,n in network['gas_species'][gas].items():balance[el]=balance.get(el,0)+nu*n
        if any(balance.values()):raise ValueError('Unbalanced '+e['id']+': '+str(balance))
    if abs(sum(network['initial_motif_fractions'].values())-1)>1e-12:raise ValueError('Initial fractions must sum to one')

def rate_constants(network,temperature=400.,prefactor=1e12,adsorption=5.,desorption=100.):
    if not all(math.isfinite(x) and x>0 for x in [temperature,prefactor]) or not all(math.isfinite(x) and x>=0 for x in [adsorption,desorption]):raise ValueError('Invalid rates or temperature')
    rates=[]
    for e in network['events']:
        if not e['enabled']:rates.append(0.);continue
        if e['kind']=='reaction':
            if e['barrier_eV'] is None:raise ValueError('Missing reaction barrier')
            rates.append(prefactor*math.exp(-e['barrier_eV']/(KB*temperature)))
        else:rates.append(adsorption if e['kind']=='adsorption' else desorption)
    return np.array(rates)

def initial_counts(network,sites):
    if not isinstance(sites,int) or sites<1:raise ValueError('Positive integer sites required')
    names=[s['id'] for s in network['states']];counts=np.zeros(len(names),dtype=np.int64)
    desired=[(names.index(k),sites*v) for k,v in network['initial_motif_fractions'].items()]
    for i,n in desired:counts[i]=int(n)
    for i,n in sorted(desired,key=lambda x:x[1]%1,reverse=True)[:sites-int(counts.sum())]:counts[i]+=1
    return counts

def simulate(network,rates,sites,times,seed=0,exposure_end=2.,track_lattice=False):
    if not math.isfinite(exposure_end) or exposure_end<0:raise ValueError('Invalid exposure duration')
    times=np.asarray(times,float)
    if not len(times) or np.any(~np.isfinite(times)) or np.any(np.diff(times)<0) or times[0]<0:raise ValueError('Invalid sample times')
    rates=np.asarray(rates,float)
    if rates.shape!=(len(network['events']),) or np.any(~np.isfinite(rates)) or np.any(rates<0):raise ValueError('Invalid rate vector')
    src=np.array([e['source'] for e in network['events']]);dst=np.array([e['target'] for e in network['events']]);ads=np.array([e['kind']=='adsorption' for e in network['events']]);gases=list(network['gas_species']);delta=np.array([[e['gas_delta'].get(g,0) for g in gases] for e in network['events']],dtype=np.int64)
    counts=initial_counts(network,sites);gas=np.zeros(len(gases),dtype=np.int64);hist=[];gh=[];frames=[];t=0.;sample=0;events=0;rng=np.random.default_rng(seed);event_counts=np.zeros(len(rates),dtype=np.int64)
    if track_lattice:
        lattice=np.repeat(np.arange(len(counts)),counts);rng.shuffle(lattice);pools=[list(np.flatnonzero(lattice==i)) for i in range(len(counts))]
    while sample<len(times):
        active=rates.copy()
        if t>=exposure_end:active[ads]=0
        hazards=active*counts[src];total=float(hazards.sum());candidate=t+rng.exponential(1/total) if total else math.inf
        boundary=exposure_end if t<exposure_end else math.inf;next_t=min(candidate,boundary)
        while sample<len(times) and times[sample]<=next_t:
            hist.append(counts.copy());gh.append(gas.copy())
            if track_lattice:frames.append(lattice.copy())
            sample+=1
        if sample==len(times):break
        if boundary<=candidate:t=boundary;continue
        event=int(np.searchsorted(np.cumsum(hazards),rng.random()*total,side='right'));event=min(event,len(rates)-1);i=src[event];j=dst[event]
        if counts[i]<=0:raise RuntimeError('Negative source inventory')
        counts[i]-=1;counts[j]+=1;gas+=delta[event];event_counts[event]+=1;events+=1;t=candidate
        if track_lattice:
            index=int(rng.integers(len(pools[i])));site=pools[i][index];pools[i][index]=pools[i][-1];pools[i].pop();pools[j].append(site);lattice[site]=j
    return dict(times=times,counts=np.array(hist),gas_counts=np.array(gh),lattice=np.array(frames) if track_lattice else None,event_counts=event_counts,events=events)

def expectation(network,rates,sites,times,exposure_end=2.):
    from scipy.integrate import solve_ivp
    n=len(network['states']);gases=list(network['gas_species']);times=np.asarray(times);y=np.r_[initial_counts(network,sites).astype(float),np.zeros(len(gases))];out=np.zeros((len(times),len(y)));out[times==0]=y
    for start,end,exposure in [(0.,min(exposure_end,float(times[-1])),True),(exposure_end,float(times[-1]),False)]:
        if end<=start:continue
        A=np.zeros((len(y),len(y)))
        for e,k in zip(network['events'],rates):
            if not exposure and e['kind']=='adsorption':continue
            i=e['source'];j=e['target'];A[i,i]-=k;A[j,i]+=k
            for g,nu in e['gas_delta'].items():A[n+gases.index(g),i]+=k*nu
        mask=(times>=start)&(times<=end);ts=times[mask]
        sol=solve_ivp(lambda t,v:A@v,(start,end),y,method='Radau',jac=A,rtol=2e-9,atol=1e-10,dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        out[mask]=sol.sol(ts).T;y=sol.y[:,-1]
    return out[:,:n],out[:,n:]
