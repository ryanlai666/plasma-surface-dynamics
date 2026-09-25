import numpy as np
import pytest
from plasma_surface.model import Parameters, Phase, ale_recipe, kmc
from plasma_surface.animate import sample_trajectory


def test_observation_preserves_seeded_trajectory_and_copies_states():
    phases = ale_recipe(35, 2, 4)
    snapshots = []
    expected = kmc(phases, cycles=2, sites=64, seed=4)
    actual = kmc(phases, cycles=2, sites=64, seed=4, observer=snapshots.append)
    assert actual == expected
    assert snapshots[0]['events'] == 0
    assert not snapshots[0]['modified'].any()
    assert not snapshots[0]['heights'].any()
    for s in snapshots:
        assert s['heights'].sum() == s['deposited'] - s['removed']
    assert all(a['time_s'] <= b['time_s'] for a, b in zip(snapshots, snapshots[1:]))
    assert snapshots[-1]['modified'].mean() == actual[-1]['coverage']
    assert -snapshots[-1]['heights'].mean() * Parameters().layer_nm == pytest.approx(
        actual[-1]['net_removed_nm']
    )


def test_observer_mutation_does_not_change_solver():
    def mutate(s):
        s['modified'][:] = False
        s['heights'][:] = 999

    assert kmc(ale_recipe(), sites=16, observer=mutate) == kmc(ale_recipe(), sites=16)


def test_uniform_sampling_keeps_idle_and_endpoints():
    history, snapshots = sample_trajectory(
        [Phase('idle', 1)], Parameters(), cycles=1, side=4, frames=11
    )
    assert [s['time_s'] for s in snapshots] == pytest.approx(np.linspace(0, 1, 11))
    assert all(s['events'] == 0 for s in snapshots)
    assert snapshots[-1]['time_s'] == history[-1]['time_s']
    assert not snapshots[-1]['heights'].any()
