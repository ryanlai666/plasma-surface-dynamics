"""Physical bookkeeping, phase switching, numerical reference and native parity."""

from copy import deepcopy
from pathlib import Path
import hashlib, json
from plasma_surface.provenance import matches_recorded
import numpy as np
import pytest
from plasma_surface.species_kmc import (
    build_network,
    rate_constants,
    initial_counts,
    simulate,
    expectation,
    validate,
)

ROOT = Path(__file__).resolve().parents[1]


def atom_ledger(n, r, sites):
    elements = ['Si', 'N', 'H', 'F']
    S = np.array([[s['composition'].get(e, 0) for e in elements] for s in n['states']])
    G = np.array([[s.get(e, 0) for e in elements] for s in n['gas_species'].values()])
    assert np.all(r['counts'] @ S + r['gas_counts'] @ G == initial_counts(n, sites) @ S)
    assert np.all(r['counts'].sum(axis=1) == sites) and np.all(r['counts'] >= 0)


def test_all_source_events_have_balanced_named_states():
    n = build_network()
    validate(n)
    assert len(n['states']) == 45 and len(n['events']) == 55
    assert len({e['pathway'] for e in n['events'] if e['pathway']}) == 16
    n['events'][-1]['gas_delta'] = {'HF': 1}
    with pytest.raises(ValueError, match='Unbalanced'):
        validate(n)


def test_missing_reaction_barrier_never_becomes_a_rate():
    n = build_network()
    next(e for e in n['events'] if e['kind'] == 'reaction')['barrier_eV'] = None
    with pytest.raises(ValueError, match='Missing'):
        rate_constants(n)


def test_exposure_switch_conservation_and_actual_lattice():
    n = build_network()
    r = simulate(n, rate_constants(n), 121, np.linspace(0, 3, 31), seed=171, track_lattice=True)
    atom_ledger(n, r, 121)
    for grid, counts in zip(r['lattice'], r['counts']):
        assert np.array_equal(np.bincount(grid, minlength=len(n['states'])), counts)
    hf = list(n['gas_species']).index('HF')
    after = r['times'] >= 2
    assert np.all(np.diff(r['gas_counts'][after, hf]) >= 0)  # no new HF uptake after switch


def test_missing_terminal_release_remains_blocked():
    n = build_network()
    n['initial_motif_fractions'] = {'terminal_NH2_F0': 1.0}
    r = simulate(n, rate_constants(n, temperature=600.0), 100, [0, 2, 3], seed=55)
    assert r['gas_counts'][-1, list(n['gas_species']).index('SiF4')] == 0
    ids = [i for i, s in enumerate(n['states']) if s['id'].startswith('terminal_NH2_F3')]
    assert r['counts'][-1, ids].sum() > 95
    atom_ledger(n, r, 100)


def test_zero_dose_and_zero_time_are_well_defined():
    n = build_network()
    initial = initial_counts(n, 31)
    rates = rate_constants(n)
    for times, exposure in [([0], 2.0), ([0, 1, 2], 0.0)]:
        r = simulate(n, rates, 31, times, seed=0, exposure_end=exposure)
        assert np.all(r['counts'] == initial) and r['events'] == 0
        mean, gas = expectation(n, rates, 31, np.array(times), exposure_end=exposure)
        assert np.allclose(mean, initial) and np.allclose(gas, 0)


def test_master_equation_predicts_independent_python_ensemble():
    n = build_network()
    rates = rate_constants(n, temperature=350.0)
    times = np.array([0.0, 1.0, 2.0, 3.0])
    mean, gas = expectation(n, rates, 300, times)
    samples = np.array([simulate(n, rates, 300, times, seed=800 + i)['counts'] for i in range(64)])
    se = samples.std(axis=0, ddof=1) / 8
    mask = (se > 0.1) & (mean > 2)
    assert np.max(abs(samples.mean(axis=0) - mean)[mask] / se[mask]) < 6
    assert np.max(abs(mean.sum(axis=1) - 300)) < 1e-5


def test_native_species_backend_and_budget():
    from plasma_surface.species_native import simulate as native, library

    try:
        library()
    except OSError:
        pytest.skip('Build optional species backend to exercise native checks')
    n = build_network()
    rates = rate_constants(n)
    r = native(n, rates, 100, [0, 1, 2, 3], seed=10)
    atom_ledger(n, r, 100)
    with pytest.raises(RuntimeError, match='budget'):
        native(n, rates, 100, [0, 3], max_events=1)


def test_published_campaign_passed_and_source_hashes_match():
    d = json.loads((ROOT / 'docs/species_kmc_results/summary.json').read_text())
    assert d['states'] == 45 and d['enabled_events'] == 55
    assert d['status'].startswith('conditional_') and d['distinct_source_pathways'] == 16
    assert d['python_cpp_max_standard_errors'] < 6 and all(
        c['conservation_exact'] and c['max_standard_errors'] < 6 for c in d['verification']
    )
    for p, sha in d['code_sha256'].items():
        assert matches_recorded(ROOT / p, sha)
    assert matches_recorded(ROOT / 'data/literature/sin_hf_pathways.csv', d['source_sha256'])


def test_species_animation_is_backed_by_actual_saved_states():
    from PIL import Image

    p = ROOT / 'docs/species_kmc_results'
    d = json.loads((p / 'animation_manifest.json').read_text())
    traj = np.load(p / 'lattice_trajectory.npz')
    assert d['interpolated_frames'] == 0 and len(traj['times']) == d['frame_count'] == 41
    with Image.open(p / 'species_kmc.gif') as im:
        assert im.n_frames == 41
    for grid, counts in zip(traj['states'], traj['counts']):
        assert np.array_equal(np.bincount(grid, minlength=45), counts)
    assert matches_recorded(p / 'lattice_trajectory.npz', d['trajectory_sha256'])


def test_network_plot_includes_all_events_and_matches_animation_palette():
    p = ROOT / 'docs/reaction_network'
    d = json.loads((p / 'render_manifest.json').read_text())
    n = build_network()
    assert d['state_ids'] == [s['id'] for s in n['states']]
    assert d['event_ids'] == [e['id'] for e in n['events']]
    assert matches_recorded(ROOT / 'configs/species_kmc_network.json', d['network_sha256'])
    a = json.loads((ROOT / 'docs/species_kmc_results/animation_manifest.json').read_text())
    assert d['state_colors'] == a['state_colors']
    assert matches_recorded(ROOT / a['renderer'], a['renderer_sha256'])


def test_released_animation_replays_exactly_with_recorded_seed():
    p = ROOT / 'docs/species_kmc_results'
    traj = np.load(p / 'lattice_trajectory.npz')
    n = build_network()
    r = simulate(n, rate_constants(n), 400, traj['times'], seed=741, track_lattice=True)
    for saved, key in [('states', 'lattice'), ('counts', 'counts'), ('gas_counts', 'gas_counts')]:
        assert np.array_equal(traj[saved], r[key])
    atom_ledger(n, r, 400)
    # This is a single inventory: every site stays in its initial motif family.
    family = np.array([s['family'] for s in n['states']])
    assert np.all(family[traj['states']] == family[traj['states'][0]])
    removed = np.array([s['Si_removed'] for s in n['states']])
    assert np.all(np.diff(removed[traj['states']].astype(int), axis=0) >= 0)


def test_sensitivity_mixture_identity_and_topological_ceiling():
    p = ROOT / 'docs/species_kmc_results/kinetic_priorities.json'
    d = json.loads(p.read_text())
    n = build_network()
    assert d['mixture_linearity_error'] < 1e-7
    assert matches_recorded(ROOT / d['runner'], d['runner_sha256'])
    total = 0.0
    names = [s['id'] for s in n['states']]
    for name, weight in n['initial_motif_fractions'].items():
        reached = {names.index(name)}
        while True:
            new = reached | {
                e['target'] for e in n['events'] if e['source'] in reached and e['enabled']
            }
            if new == reached:
                break
            reached = new
        if any(n['states'][i]['Si_removed'] for i in reached):
            total += weight
    assert total == pytest.approx(0.6)
    for r in d['sensitivities']:
        expected = np.log(r['higher_rate_yield'] / r['lower_rate_yield']) / (
            2 * np.log(r['relative_step'])
        )
        assert expected == pytest.approx(r['log_yield_sensitivity'], abs=1e-12)


def test_separate_literature_channels_preserve_source_units():
    p = ROOT / 'docs/species_kmc_results/literature_channel_rates.json'
    d = json.loads(p.read_text())
    assert matches_recorded(
        (ROOT / 'data/literature/jung2020_coadsorption_rates.json'), d['source_sha256']
    )
    for r in d['rows']:
        assert r['rate_s'] == pytest.approx(
            r['A_s'] * np.exp(-r['Ea_over_R_K'] / r['temperature_K'])
        )
        assert r['barrier_eV'] == pytest.approx(r['Ea_over_R_K'] * 8.617333262145e-5)
