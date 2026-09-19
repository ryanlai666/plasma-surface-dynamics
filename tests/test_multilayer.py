from pathlib import Path
from collections import Counter
from copy import deepcopy
import hashlib,json
import numpy as np
import pytest
from plasma_surface.multilayer import simulate,validate,accessibility,inventory,apply,candidates
from plasma_surface.rate_evidence import environment,qualified_rate
ROOT=Path(__file__).resolve().parents[1]

def initial():return json.loads((ROOT/'data/multilayer/sin_graph.json').read_text())

def test_crystal_graph_has_six_bands_and_valence_complete_terminations():
    g=initial();validate(g);assert len(g['nodes'])==336 and len(g['bonds'])==552
    assert {n['layer'] for n in g['nodes']}==set(range(6))
    assert {x for n in g['nodes'] for x in n['termination']}=={'H','F','Cl'}
    assert sum(n['fixed'] for n in g['nodes'])==56
    assert g['source_sha256']==hashlib.sha256((ROOT/g['source']).read_bytes()).hexdigest()

def test_every_saved_multilayer_frame_balances_all_elements_and_preserves_base():
    d=json.loads((ROOT/'data/multilayer/demo_trajectory.json').read_text());g=initial()
    for frame in d['snapshots']:
        g['bonds']=frame['bonds']
        for i,n in enumerate(g['nodes']):
            n['active']=frame['active'][i];n['termination']=[e for e,count in frame['terminations'][i].items() for _ in range(count)]
            if n['fixed']:assert n==d['initial']['nodes'][i]
        validate(g);total=inventory(g,frame['precursors'])
        for name,amount in frame['gas'].items():
            for element,count in d['gas_formulas'][name].items():total[element]+=amount*count
        assert all(total[e]==d['initial_inventory'].get(e,0) for e in set(total)|set(d['initial_inventory']))
    assert any(not n['active'] and n['layer']<5 for n in g['nodes'])

def test_fresh_subsurface_sites_are_exposed_by_real_removal():
    d=json.loads((ROOT/'data/multilayer/demo_trajectory.json').read_text());assert any(e['newly_exposed'] for e in d['events'] if e['kind'].endswith('release'))
    g=initial();before=accessibility(g,0)[1];i=int(np.flatnonzero(before)[0]);col=g['nodes'][i]['column'];g['nodes'][i]['active']=False;after=accessibility(g,0)[1]
    assert not after[i]
    assert any(a and not b and g['nodes'][j]['column']==col for j,(a,b) in enumerate(zip(after,before)))

def test_no_qualified_rates_means_no_fabricated_strict_chemistry():
    g=initial();d=simulate(g,[0,1,3],policy='validated_only');assert d['events']==[] and d['final']==g
    options=candidates(g,np.zeros(336,bool),policy='demonstration');event=options[0];env=environment(g,np.zeros(336,bool),event,450.)
    assert qualified_rate(None,env,event,{})==(None,'no_matched_record')
    assert qualified_rate({'rate_s':1},env,event,{})==(None,'incomplete_evidence')

def test_cannot_release_a_bound_substrate_atom():
    g=initial();i=next(n['id'] for n in g['nodes'] if n['element']=='Si' and not n['fixed']);event=dict(kind='Si_molecule_release',site=i,partner=None)
    with pytest.raises(ValueError,match='bonded'):apply(g,np.zeros(336,bool),event,{}, {})

def test_access_depth_is_relative_to_current_front_and_terminations_change_keys():
    g=initial();d,ex,re=accessibility(g,2.);assert np.all(re[d>2.]==False) and np.all(re[ex])
    prec=np.zeros(336,bool);event=candidates(g,prec,policy='demonstration')[0];a=environment(g,prec,event,450.)
    g['nodes'][event['site']]['termination'].append('Cl');b=environment(g,prec,event,450.);assert a['key']!=b['key']

def test_multilayer_artifact_provenance_and_rate_queue_cover_events():
    out=ROOT/'docs/multilayer_results';s=json.loads((out/'summary.json').read_text());assert s['validated_only_events']==0
    for path,sha in s['code_sha256'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha
    a=json.loads((out/'animation_manifest.json').read_text());assert a['interpolated_frames']==0 and a['frame_count']==37
    assert a['trajectory_sha256']==s['trajectory_sha256']==hashlib.sha256((ROOT/'data/multilayer/demo_trajectory.json').read_bytes()).hexdigest()
    q=json.loads((ROOT/'data/multilayer/rate_requests/index.json').read_text());assert q['trajectory_sha256']==s['trajectory_sha256']
    assert sum(r['observed_events'] for r in q['requests'])==s['baseline_events']
    assert any(r['observed_events']==0 and r['candidate_snapshots']>0 for r in q['requests'])

def test_surface_drawings_preserve_atom_inventory_without_invented_residuals():
    d=json.loads((ROOT/'docs/reaction_network/render_manifest.json').read_text());n=json.loads((ROOT/'configs/species_kmc_network.json').read_text());assert len(d['depictions'])==45
    for s,p in zip(n['states'],d['depictions']):
        assert s['id']==p['state_id']
        if p['unknown']:assert not p['nodes'] and s['Si_removed'];continue
        assert dict(Counter(a['element'] for a in p['nodes']))==s['composition']
        assert all(max(b['a'],b['b'])<len(p['nodes']) for b in p['bonds'])
    assert d['depiction_renderer_sha256']==hashlib.sha256((ROOT/d['depiction_renderer']).read_bytes()).hexdigest()


def test_unparameterized_fourth_fluorination_is_not_given_a_third_step_barrier():
    d=json.loads((ROOT/'data/multilayer/demo_trajectory.json').read_text());g=deepcopy(d['final']);prec=np.array(d['snapshots'][-1]['precursors'])
    # Explore all available precursor states without changing valence bookkeeping.
    for i,n in enumerate(g['nodes']):
        if n['active'] and n['element']=='Si' and n['termination'].count('F')==3:prec[i]=True
    ee=candidates(g,prec,penetration_A=100.,policy='demonstration')
    final=[e for e in ee if e['kind']=='SiN_cleavage' and g['nodes'][e['site']]['termination'].count('F')>=3]
    assert final
    for e in final:
        h=g['nodes'][e['partner']]['termination'].count('H')
        if h==1:assert e['reference_pathway']=='P4'
        else:assert e['rate_s']==0 and e['rate_status'].startswith('disabled_missing')


def test_rate_evidence_rejects_bad_context_artifacts_uncertainty_and_saddle(monkeypatch,tmp_path):
    import plasma_surface.rate_evidence as module
    monkeypatch.setattr(module,'ROOT',tmp_path)
    p=tmp_path/'fixture.txt';p.write_text('unit-test fixture only')
    event={'kind':'SiN_cleavage'};env={'key':'test'};context={'temperature_K':450.}
    r=dict(id='test',environment_key='test',reaction='SiN_cleavage',rate_s=1.,source='unit-test',method='DFT_TST',context=context,geometry_match_verified=True,applicability_reviewed=True,reference_state='same_adsorbed_or_surface_state',uncertainty={'basis':'unit-test interval','rate_s_interval':[.5,2.]},artifact_sha256={'fixture.txt':hashlib.sha256(p.read_bytes()).hexdigest()},checks=dict(imaginary_modes=1,IS_minimum=True,FS_minimum=True,TS_force_converged=True,connectivity=True,prefactor_or_free_energy=True,basis_cell_convergence=True,charge_spin_reviewed=True))
    assert qualified_rate(r,env,event,context)==(1.,'qualified_record')
    assert qualified_rate(r,env,event,{'temperature_K':400.})[0] is None
    rr=deepcopy(r);rr['checks']['connectivity']=False;assert qualified_rate(rr,env,event,context)[0] is None
    rr=deepcopy(r);rr['uncertainty']['rate_s_interval']=[2.,3.];assert qualified_rate(rr,env,event,context)[0] is None
    p.write_text('changed');assert qualified_rate(r,env,event,context)[0] is None


def test_multilayer_visualization_matches_saved_states_and_event_exposure():
    from PIL import Image
    d=json.loads((ROOT/'data/multilayer/demo_trajectory.json').read_text())
    out=ROOT/'docs/multilayer_results'
    m=json.loads((out/'animation_manifest.json').read_text())
    assert m['renderer_sha256']==hashlib.sha256((ROOT/'scripts/render_multilayer.py').read_bytes()).hexdigest()
    assert m['gif_sha256']==hashlib.sha256((out/'multilayer_kmc.gif').read_bytes()).hexdigest()
    assert Image.open(out/'multilayer_kmc.gif').n_frames==len(d['snapshots'])==len(m['frames'])
    gif=Image.open(out/'multilayer_kmc.gif');gif.seek(4)
    assert np.any(np.all(np.asarray(gif.convert('RGB'))==[230,165,26],axis=-1))  # HF gold survives GIF palette

    assert m['cross_section']['node_ids']==[n['id'] for n in d['initial']['nodes'] if n['column']//6==2]
    seen=set()
    for k,(s,r) in enumerate(zip(d['snapshots'],m['frames'])):
        expected={i for e in d['events'] if e['time_s']<=s['time_s'] for i in e['newly_exposed']}
        assert set(r['ever_exposed_ids'])==expected
        assert set(r['first_exposed_ids'])==expected-seen
        seen=expected
        assert r['removed']==dict(Counter(n['element'] for n,a in zip(d['initial']['nodes'],s['active']) if not a))
        assert r['HF_occupancy']==sum(p and a for p,a in zip(s['precursors'],s['active']))
        assert r['termination_counts']=={e:sum(t.get(e,0) for t,a in zip(s['terminations'],s['active']) if a) for e in ('H','F','Cl')}
        assert r['drawn_bonds']+r['periodic_bonds_omitted']==len(s['bonds'])
