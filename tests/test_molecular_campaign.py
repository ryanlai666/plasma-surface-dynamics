"""Audit published scientific artifacts without downloading ML checkpoints."""
from pathlib import Path
import csv,hashlib,json
import numpy as np
import pytest
from PIL import Image
from test_surface_paths import frames
ROOT=Path(__file__).resolve().parents[1]
SURFACES=['Si100','Si111','beta_Si3N4_001','alpha_quartz_001']
MOLECULES=['HF','HCl','F2','Cl2','H2','H2O','CH3F','SiF4','SiCl4']

@pytest.mark.parametrize('surface',SURFACES)
def test_site_orientation_matrix_preserves_geometries_and_energy_references(surface):
    total=0
    for species in MOLECULES:
        parent=ROOT/'data/surface_paths'/surface/species
        base=frames(parent/'mace_approach/slab.extxyz')[0][2]
        gas=frames(parent/'mace_approach/molecule.extxyz')[0][2]
        gas_dist=np.linalg.norm(gas[:,None,:]-gas[None,:,:],axis=2)
        paths=list((parent/'mace_site_orientation').glob('*/*/summary.json'))
        assert len(paths)==(6 if species in ['F2','Cl2','H2'] else 9)
        for path in paths:
            p=path.parent;s=json.loads(path.read_text());fs=frames(p/'images.extxyz');rows=list(csv.DictReader((p/'energies.csv').open()))
            assert len(fs)==len(rows)==9 and s['peak_is_TS'] is False
            assert s['script_sha256']==hashlib.sha256((ROOT/s['runner']).read_bytes()).hexdigest()
            for (meta,symbols,xyz),r in zip(fs,rows):
                assert np.isfinite(xyz).all() and meta['pbc']=='T T F'
                assert np.allclose(xyz[:len(base)],base,atol=1e-8)
                molecular=xyz[len(base):]
                distances=np.linalg.norm(molecular[:,None,:]-molecular[None,:,:],axis=2)
                assert np.allclose(distances,gas_dist,atol=3e-8)
                assert molecular[:,2].min()-base[:,2].max()==pytest.approx(float(r['height_A']),abs=2e-8)
                assert np.allclose(molecular[:,:2].mean(axis=0),s['site_position_A'][:2],atol=2e-8)
                assert float(meta['energy'])==pytest.approx(float(r['energy_eV']),abs=1e-7)
                reference=float(r['energy_eV'])-s['slab_energy_eV']-s['molecular_layer_energy_eV']
                assert reference==pytest.approx(float(r['interaction_energy_eV']),abs=1e-8)
                total+=1
        all_rows=list(csv.DictReader((parent/'mace_site_orientation/all_energies.csv').open()))
        assert len(all_rows)==len(paths)*9
    assert total==648


def test_molecular_gifs_use_only_evaluated_frames():
    overview=json.loads((ROOT/'docs/dry_etch_results/molecular_screening.json').read_text())
    assert len(overview)==36
    for r in overview:
        p=ROOT/r['folder'];m=json.loads((p/'render_manifest.json').read_text())
        assert m['frame_count']==9 and m['interpolated_frames']==0
        with Image.open(p/'path.gif') as im:assert im.n_frames==9
        for file,key in [('images.extxyz','coordinates_sha256'),('energies.csv','energies_sha256'),('path.gif','gif_sha256')]:
            assert hashlib.sha256((p/file).read_bytes()).hexdigest()==m[key]


def test_reaction_labels_do_not_promote_failed_searches_to_verified_ts():
    summaries=list((ROOT/'data/surface_paths').glob('*/*/molecular_neb/summary.json'))
    assert len(summaries)==7
    for p in summaries:
        s=json.loads(p.read_text());assert s['status']!='endpoint_search'
        assert s['script_sha256']==hashlib.sha256((ROOT/s['runner']).read_bytes()).hexdigest()
        if s.get('endpoint_connectivity_confirmed'):
            assert s['status']=='neb_converged' and s['index_one_in_mobile_subspace']
            assert all(s[k]['converged'] and s[k]['no_negative_curvature'] for k in ['IS','FS'])
        if 'peak_above_IS_eV' in s:
            rows=list(csv.DictReader((p.parent/'energies.csv').open()));es=np.array([float(r['energy_eV']) for r in rows])
            assert s['peak_above_IS_eV']==pytest.approx(es.max()-es[0],abs=1e-8)


def test_direct_dft_results_preserve_geometry_and_reject_unresolved_root():
    paths=list((ROOT/'data/surface_paths').glob('*_capped_motif/HF/dft_*/summary.json'))
    assert len(paths)==4
    for p in paths:
        s=json.loads(p.read_text());assert s['status']=='complete'
        assert s['script_sha256']==hashlib.sha256((ROOT/s['runner']).read_bytes()).hexdigest()
        assert s['input_sha256']==hashlib.sha256((p.parent/'input.extxyz').read_bytes()).hexdigest()
        if p.parents[2].name=='SiO_capped_motif' and p.parent.name=='dft_path':
            assert s['accepted_frame_count']==10
            assert s['frames'][5]['accepted_for_comparison'] is False
            assert sum(r['accepted_for_comparison'] for r in s['frames'])==10
            assert 'rejected' in s['electronic_quality']
            assert s['SCF_audit_script_sha256']==hashlib.sha256((ROOT/s['SCF_audit_runner']).read_bytes()).hexdigest()
        before=frames(p.parent/'input.extxyz');after=frames(p.parent/'images.extxyz')
        assert len(before)==len(after)==len(s['frames'])
        for a,b,r in zip(before,after,s['frames']):
            assert np.array_equal(a[2],b[2]) and a[1]==b[1]
            assert r['SCF_converged']
            assert float(b[0]['energy'])==pytest.approx(r['energy_eV'],abs=1e-8)


def test_energy_only_evaluator_matches_independent_regular_evaluations():
    d=json.loads((ROOT/'docs/dry_etch_results/molecular_energy_checks.json').read_text())
    assert d['status']=='passed' and d['geometries']==24
    assert d['max_energy_only_error_eV']<1e-7 and d['max_reproduction_error_eV']<1e-6
    assert d['script_sha256']==hashlib.sha256((ROOT/d['runner']).read_bytes()).hexdigest()


def test_refined_reaction_artifacts_and_ts_classification():
    paths=list((ROOT/'data/surface_paths').glob('*/*/molecular_refined/summary.json'))
    assert len(paths)==2
    for p in paths:
        d=json.loads(p.read_text());assert d['status'] in ['neb_converged','neb_unconverged']
        assert d['script_sha256']==hashlib.sha256((ROOT/d['runner']).read_bytes()).hexdigest()
        fs=frames(p.parent/'images.extxyz');rows=list(csv.DictReader((p.parent/'energies.csv').open()))
        assert len(fs)==len(rows)==d['images']
        fixed=[i for i in range(len(fs[0][1])) if i not in d['mobile_indices']]
        for f,r in zip(fs,rows):
            assert np.array_equal(f[2][fixed],fs[0][2][fixed])
            assert float(f[0]['energy'])==pytest.approx(float(r['energy_eV']),abs=1e-7)
        if d.get('endpoint_connectivity_confirmed'):
            assert d['status']=='neb_converged' and d['index_one_in_mobile_subspace']
            assert d['max_neb_force_eV_A']<=.035
            assert sum(x<-.02 for x in d['peak_hessian_eigenvalues_eV_A2'])==1
            assert all(d[k]['converged'] and d[k]['no_negative_curvature'] for k in ['IS','FS'])


def test_local_saddle_uses_the_minima_actually_connected():
    p=ROOT/'data/surface_paths/beta_Si3N4_001/H2O/mace_local_saddle'
    d=json.loads((p/'summary.json').read_text());assert d['status']=='constrained_saddle_connected'
    assert d['endpoint_connectivity_confirmed'] and not d['original_precursor_connectivity_confirmed']
    assert all(d[k]['converged'] and min(d[k]['hessian_eigenvalues_eV_A2'])>0 for k in ['IS','FS'])
    assert sum(x<-.02 for x in d['peak_hessian_eigenvalues_eV_A2'])==1
    fs=frames(p/'images.extxyz');es=np.array([float(f[0]['energy']) for f in fs])
    assert d['peak_image']==int(es.argmax()) and d['peak_above_IS_eV']==pytest.approx(es.max()-es[0],abs=1e-8)
    audit=json.loads((p/'identity_audit.json').read_text());assert 'N25' in audit['event'] and 'N26' in audit['event']
    assert all(r['OH_bond_A']<1.2 and r['separated_O_H_A']>2 for r in audit['rows'])
    assert audit['rows'][0]['H_N25_A']<1.2 and audit['rows'][-1]['H_N26_A']<1.2
    assert audit['path_sha256']==hashlib.sha256((p/'images.extxyz').read_bytes()).hexdigest()
    assert d['script_sha256']==hashlib.sha256((ROOT/d['runner']).read_bytes()).hexdigest()
