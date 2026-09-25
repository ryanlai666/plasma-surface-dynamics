from dataclasses import replace
import numpy as np
import pytest
from plasma_surface.model import Parameters, Phase, ale_recipe, mean_field, kmc
from plasma_surface.workflows import campaign, calibrate, train_surrogate
from plasma_surface.io import write_csv


def test_adsorption_analytic_limit():
    p = replace(Parameters(), desorption_prefactor_s=0)
    phase = Phase("adsorb", 2, radical_flux_m2_s=p.site_density_m2)
    result = mean_field([phase], p, 1)[-1]
    assert result["coverage"] == pytest.approx(1 - np.exp(-p.sticking * 2))
    assert result["net_removed_nm"] == 0


def test_no_flux_and_subthreshold_removal():
    assert mean_field([Phase("idle", 10)], cycles=1)[-1]["net_removed_nm"] == 0
    assert mean_field(ale_recipe(energy_ev=10))[-1]["net_removed_nm"] == 0


def test_growth_sign_and_balance():
    p = Parameters()
    phase = Phase("growth", 3, precursor_flux_m2_s=p.site_density_m2)
    result = mean_field([phase], p, 1)[-1]
    assert result["net_removed_nm"] == pytest.approx(-3 * p.deposition_sticking * p.layer_nm)
    stochastic = kmc([phase], p, 1, seed=12)[-1]
    assert stochastic["net_removed_nm"] == stochastic["removed_nm"] - stochastic["deposited_nm"]


def test_kmc_ensemble_matches_mean_field():
    phases = ale_recipe()
    reference = mean_field(phases, cycles=2)[-1]
    samples = [kmc(phases, cycles=2, sites=128, seed=i)[-1] for i in range(48)]
    for name in ["net_removed_nm", "coverage"]:
        values = np.array([s[name] for s in samples])
        assert (
            abs(values.mean() - reference[name])
            < 5 * values.std(ddof=1) / np.sqrt(len(values)) + 1e-5
        )


def test_reproducibility_and_event_budget():
    assert kmc(ale_recipe(), sites=32, seed=7) == kmc(ale_recipe(), sites=32, seed=7)
    with pytest.raises(RuntimeError):
        kmc(ale_recipe(), max_events=1)


@pytest.mark.parametrize(
    "kwargs", [{"temperature_k": 0}, {"ion_energy_ev": -1}, {"duration_s": float("nan")}]
)
def test_invalid_phase(kwargs):
    with pytest.raises(ValueError):
        Phase(**(dict(name="invalid", duration_s=1) | kwargs))


def test_campaign_shards_match_serial():
    serial = campaign(32)
    parallel = sorted(
        campaign(32, shard=0, shards=2) + campaign(32, shard=1, shards=2),
        key=lambda r: r["recipe_id"],
    )
    assert serial == parallel


def test_surrogate_split_and_outputs(tmp_path):
    report = train_surrogate(campaign(64), tmp_path)
    assert not set(report["train_recipe_ids"]) & set(report["test_recipe_ids"])
    assert np.isfinite(report["mae_nm"])
    assert (tmp_path / "surrogate_predictions.csv").is_file()


def test_calibration_recovers_synthetic_parameters(tmp_path):
    p = replace(Parameters(), sticking=0.38, chemical_yield_scale=1.1)
    rows = []
    for dose in [0.3, 1, 3]:
        for energy in [22, 35, 60]:
            trace = mean_field(ale_recipe(energy, dose), p, cycles=2)
            rows.append(
                dict(
                    energy_ev=energy,
                    dose_s=dose,
                    ion_s=2,
                    temperature_k=300,
                    cycle=2,
                    epc_nm=trace[-1]["net_removed_nm"] - trace[-5]["net_removed_nm"],
                    sigma_nm=0.001,
                    source_id="synthetic_parameter_recovery_test",
                )
            )
    path = tmp_path / "measurements.csv"
    write_csv(path, rows)
    fit = calibrate(path, tmp_path / "fit.json")
    assert fit["success"]
    assert fit["jacobian_rank"] == 2
    assert fit["sticking"] == pytest.approx(0.38, abs=1e-5)
    assert fit["chemical_yield_scale"] == pytest.approx(1.1, abs=1e-5)
