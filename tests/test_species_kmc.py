"""Physical bookkeeping, phase switching, numerical reference and native parity."""
from copy import deepcopy
from pathlib import Path
import hashlib,json
import numpy as np
import pytest
from plasma_surface.species_kmc import build_network,rate_constants,initial_counts,simulate,expectation,validate
ROOT=Path(__file__).resolve().parents[1]

def atom_ledger(n,r,sites):
    elements=['Si','N','H','F'];S=np.array([[s['composition'].get(e,0) for e in elements] for s in n['states']]);G=np.array([[s.get(e,0) for e in elements] for s in n['gas_species'].values()])
    assert np.all(r['counts']@S+r['gas_counts']@G==initial_counts(n,sites)@S)
    assert np.all(r['counts'].sum(axis=1)==sites) and np.all(r['counts']>=0)

def test_all_source_events_have_balanced_named_states():
    n=build_network();validate(n);assert len(n['states'])==45 and len(n['events'])==55
    assert len({e['pathway'] for e in n['events'] if e['pathway']})==16
    n['events'][-1]['gas_delta']={'HF':1}
    with pytest.raises(ValueError,match='Unbalanced'):validate(n)

def test_missing_reaction_barrier_never_becomes_a_rate():
    n=build_network();next(e for e in n['events'] if e['kind']=='reaction')['barrier_eV']=None
    with pytest.raises(ValueError,match='Missing'):rate_constants(n)

def test_exposure_switch_conservation_and_actual_lattice():
    n=build_network();r=simulate(n,rate_constants(n),121,np.linspace(0,3,31),seed=171,track_lattice=True);atom_ledger(n,r,121)
    for grid,counts in zip(r['lattice'],r['counts']):assert np.array_equal(np.bincount(grid,minlength=len(n['states'])),counts)
    hf=list(n['gas_species']).index('HF');after=r['times']>=2
    assert np.all(np.diff(r['gas_counts'][after,hf])>=0)  # no new HF uptake after switch

def test_missing_terminal_release_remains_blocked():
    n=build_network();n['initial_motif_fractions']={'terminal_NH2_F0':1.};r=simulate(n,rate_constants(n,temperature=600.),100,[0,2,3],seed=55)
    assert r['gas_counts'][-1,list(n['gas_species']).index('SiF4')]==0
    ids=[i for i,s in enumerate(n['states']) if s['id'].startswith('terminal_NH2_F3')]
    assert r['counts'][-1,ids].sum()>95
    atom_ledger(n,r,100)

def test_zero_dose_and_zero_time_are_well_defined():
    n=build_network();initial=initial_counts(n,31);rates=rate_constants(n)
    for times,exposure in [([0],2.),([0,1,2],0.)]:
        r=simulate(n,rates,31,times,seed=0,exposure_end=exposure);assert np.all(r['counts']==initial) and r['events']==0
        mean,gas=expectation(n,rates,31,np.array(times),exposure_end=exposure);assert np.allclose(mean,initial) and np.allclose(gas,0)

def test_master_equation_predicts_independent_python_ensemble():
    n=build_network();rates=rate_constants(n,temperature=350.);times=np.array([0.,1.,2.,3.]);mean,gas=expectation(n,rates,300,times)
    samples=np.array([simulate(n,rates,300,times,seed=800+i)['counts'] for i in range(64)])
    se=samples.std(axis=0,ddof=1)/8;mask=(se>.1)&(mean>2)
    assert np.max(abs(samples.mean(axis=0)-mean)[mask]/se[mask])<6
    assert np.max(abs(mean.sum(axis=1)-300))<1e-5

def test_native_species_backend_and_budget():
    from plasma_surface.species_native import simulate as native,library
    try:library()
    except OSError:pytest.skip('Build optional species backend to exercise native checks')
    n=build_network();rates=rate_constants(n);r=native(n,rates,100,[0,1,2,3],seed=10);atom_ledger(n,r,100)
    with pytest.raises(RuntimeError,match='budget'):native(n,rates,100,[0,3],max_events=1)

def test_published_campaign_passed_and_source_hashes_match():
    d=json.loads((ROOT/'docs/species_kmc_results/summary.json').read_text());assert d['states']==45 and d['enabled_events']==55
    assert d['status'].startswith('conditional_') and d['distinct_source_pathways']==16
    assert d['python_cpp_max_standard_errors']<6 and all(c['conservation_exact'] and c['max_standard_errors']<6 for c in d['verification'])
    for p,sha in d['code_sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==sha
    assert d['source_sha256']==hashlib.sha256((ROOT/'data/literature/sin_hf_pathways.csv').read_bytes()).hexdigest()

def test_species_animation_is_backed_by_actual_saved_states():
    from PIL import Image
    p=ROOT/'docs/species_kmc_results';d=json.loads((p/'animation_manifest.json').read_text());traj=np.load(p/'lattice_trajectory.npz')
    assert d['interpolated_frames']==0 and len(traj['times'])==d['frame_count']==41
    with Image.open(p/'species_kmc.gif') as im:assert im.n_frames==41
    for grid,counts in zip(traj['states'],traj['counts']):assert np.array_equal(np.bincount(grid,minlength=45),counts)
    assert hashlib.sha256((p/'lattice_trajectory.npz').read_bytes()).hexdigest()==d['trajectory_sha256']
