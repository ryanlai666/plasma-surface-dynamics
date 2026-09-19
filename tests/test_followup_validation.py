from pathlib import Path
import json,hashlib
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]


def test_full_host_stability_excludes_unstable_references():
    d=json.loads((ROOT/'data/intermediate_campaign/full_stability.json').read_text())
    q=json.loads((ROOT/'data/intermediate_campaign/stability_qualification.json').read_text())
    assert d['status']=='complete' and len(d['results'])==8 and len(d['parents'])==4
    assert hashlib.sha256((ROOT/d['runner']).read_bytes()).hexdigest()==d['runner_sha256']
    for r in list(d['parents'].values())+d['results']:
        assert hashlib.sha256((ROOT/r['final_structure']).read_bytes()).hexdigest()==r['final_sha256']
        h=np.load(ROOT/Path(r['final_structure']).parent/'hessian_eV_A2.npy')
        assert np.linalg.eigvalsh(h)[0]==pytest.approx(r['minimum_eigenvalue_eV_A2'])
        assert h.shape==(3*len(r['mobile_indices']),)*2
    assert sum(r['screened_association_energy_eV'] is None for r in q['results'])>=2
    for r in q['results']:
        if r['parent_minimum_curvature_eV_A2']<-.02 or r['candidate_minimum_curvature_eV_A2']<-.02:
            assert r['screened_association_energy_eV'] is None
        assert not r['rate_enabled']


def test_direct_dft_diagnostics_remain_excluded_from_rates():
    base=ROOT/'data/final_cleavage/SiF3_NH2_HF'
    d=json.loads((base/'dft_initial_path/summary.json').read_text())
    assert d['status']=='complete' and len(d['results'])==2
    assert hashlib.sha256((ROOT/d['runner']).read_bytes()).hexdigest()==d['runner_sha256']
    for r in d['results']:
        e=np.array([x['DFT_energy_eV'] for x in r['frames']])
        assert np.allclose(e-e[0],r['relative_DFT_energies_eV'])
        assert not r['energies_are_barriers']
        assert all(x['SCF_converged'] and x['internal_stable'] for x in r['frames'])
        assert r['frames'][0]['max_mobile_DFT_force_eV_A']>.02
    for name in ('omol25','omol25_bfgs'):
        r=json.loads((base/name/'summary.json').read_text());assert r['status']=='neb_unconverged' and not r['rate_enabled']
    scan=json.loads((base/'fixed_frame_three_coordinate_scan/summary.json').read_text())
    assert scan['status']=='complete' and len(scan['results'])==22 and not scan['rate_enabled']
    for r in scan['results']:
        assert r['SiN_A']==pytest.approx(r['target_SiN_A'],abs=1e-6)
        assert r['NH_A']==pytest.approx(r['target_NH_A'],abs=1e-6)

        assert r['SiF_A']==pytest.approx(r['target_SiF_A'],abs=1e-6)


def test_corrected_scan_frame_and_invalid_run_are_distinguished():
    base=ROOT/'data/final_cleavage/SiF3_NH2_HF'
    def frames(path):
        lines=path.read_text().splitlines();result=[];i=0
        while i<len(lines):
            n=int(lines[i]);result.append(np.array([[float(x) for x in line.split()[1:4]] for line in lines[i+2:i+2+n]]));i+=n+2
        return result
    initial=frames(base/'omol25/IS.extxyz')[0]
    for xyz in frames(base/'fixed_frame_three_coordinate_scan/evaluated_geometries.extxyz'):
        assert np.allclose(xyz[:4],initial[:4],atol=1e-8)
    bad=json.loads((base/'local_coordinate_scan/INVALIDATION.json').read_text())
    assert bad['maximum_fixed_atom_displacement_A']>1 and not bad['rates_enabled']
