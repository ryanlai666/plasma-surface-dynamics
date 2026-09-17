import numpy as np
import pytest
from plasma_surface import model, native

pytestmark = pytest.mark.skipif(not native.library_path().is_file(), reason="Optional C++ backend not built")


@pytest.mark.parametrize("phases", [model.ale_recipe(), model.ale_recipe(75, .5, 3, 450),
    [model.Phase("growth", 2, precursor_flux_m2_s=1e19)],
    [model.Phase("mixed", 2, radical_flux_m2_s=2e19, ion_flux_m2_s=2e19, ion_energy_ev=60, precursor_flux_m2_s=1e19)],
    [model.Phase("idle", 2)]])
def test_native_mean_field_matches_python(phases):
    reference = model.mean_field(phases)
    result = native.mean_field(phases)
    for expected, actual in zip(reference, result):
        for field in ["coverage", "net_removed_nm", "removed_nm", "deposited_nm", "time_s"]:
            assert actual[field] == pytest.approx(expected[field], abs=1e-12)


@pytest.mark.parametrize("phases", [model.ale_recipe(),
    [model.Phase("mixed", 2, radical_flux_m2_s=2e19, ion_flux_m2_s=2e19, ion_energy_ev=60, precursor_flux_m2_s=1e19)]])
def test_native_kmc_ensemble_matches_exact_solution(phases):
    target = model.mean_field(phases, cycles=2)[-1]
    result = [native.kmc(phases, cycles=2, sites=256, seed=i)[-1] for i in range(96)]
    for field in ["coverage", "net_removed_nm", "removed_nm", "deposited_nm"]:
        values = np.array([r[field] for r in result])
        assert abs(values.mean()-target[field]) < 5*values.std(ddof=1)/np.sqrt(len(values))+1e-6
    for row in result:
        assert 0 <= row["coverage"] <= 1
        assert row["net_removed_nm"] == pytest.approx(row["removed_nm"]-row["deposited_nm"])


def test_native_seed_and_event_limit():
    assert native.kmc(model.ale_recipe(), seed=12) == native.kmc(model.ale_recipe(), seed=12)
    with pytest.raises(RuntimeError, match="Event budget"):
        native.kmc(model.ale_recipe(), max_events=1)
    with pytest.raises(ValueError):
        native.kmc(model.ale_recipe(), seed=-1)
