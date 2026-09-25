"""Reproducible campaigns, surrogate evaluation, and experimental calibration."""

from dataclasses import asdict, replace
import csv
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import qmc
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import GroupShuffleSplit
from .model import Parameters, ale_recipe, mean_field, kmc
from .io import write_csv, write_json, manifest

FEATURES = ["energy_ev", "dose_s", "ion_s", "temperature_k"]


def campaign(n=128, seed=42, shard=0, shards=1, backend="python"):
    if n < 16 or not 0 <= shard < shards:
        raise ValueError("n >= 16 and 0 <= shard < shards required")
    solver = mean_field
    if backend == "cpp":
        from .native import mean_field as solver
    elif backend != "python":
        raise ValueError("backend must be python or cpp")
    points = qmc.scale(
        qmc.LatinHypercube(4, seed=seed).random(n), [5, 0.1, 0.1, 280], [90, 6, 6, 500]
    )
    rows = []
    for i, x in enumerate(points):
        if i % shards != shard:
            continue
        history = solver(ale_recipe(*x), cycles=5)
        epc = history[-1]["net_removed_nm"] - history[-5]["net_removed_nm"]
        rows.append(
            dict(
                recipe_id=i,
                **dict(zip(FEATURES, map(float, x))),
                epc_nm=float(epc),
                data_origin="synthetic_uncalibrated_model",
            )
        )
    return rows


def train_surrogate(rows, output, seed=42):
    """Held-out recipe evaluation, plus a separate high-energy extrapolation audit."""
    x = np.array([[r[f] for f in FEATURES] for r in rows], dtype=float)
    y = np.array([r["epc_nm"] for r in rows], dtype=float)
    groups = np.array([r["recipe_id"] for r in rows])
    train, test = next(
        GroupShuffleSplit(n_splits=1, test_size=0.25, random_state=seed).split(x, y, groups)
    )
    model = ExtraTreesRegressor(n_estimators=160, min_samples_leaf=2, random_state=seed, n_jobs=1)
    model.fit(x[train], y[train])
    prediction = model.predict(x[test])
    spread = np.std([tree.predict(x[test]) for tree in model.estimators_], axis=0)
    low = x[:, 0] < 65
    extrapolation = ExtraTreesRegressor(
        n_estimators=160, min_samples_leaf=2, random_state=seed, n_jobs=1
    )
    extrapolation.fit(x[low], y[low])
    metrics = dict(
        mae_nm=float(mean_absolute_error(y[test], prediction)),
        rmse_nm=float(np.sqrt(mean_squared_error(y[test], prediction))),
        constant_baseline_mae_nm=float(
            mean_absolute_error(y[test], np.full(len(test), y[train].mean()))
        ),
        high_energy_holdout_mae_nm=float(
            mean_absolute_error(y[~low], extrapolation.predict(x[~low]))
        ),
        train_recipe_ids=groups[train].tolist(),
        test_recipe_ids=groups[test].tolist(),
        uncertainty_note="Tree spread is heuristic disagreement, not calibrated confidence or physical uncertainty.",
        target="Cycle-5 net EPC from the illustrative mean-field model; no experimental validation.",
    )
    write_json(Path(output) / "surrogate_metrics.json", metrics)
    write_csv(
        Path(output) / "surrogate_predictions.csv",
        [
            dict(
                recipe_id=int(groups[i]),
                true_epc_nm=float(y[i]),
                predicted_epc_nm=float(pred),
                tree_spread_nm=float(std),
            )
            for i, pred, std in zip(test, prediction, spread)
        ],
    )
    return metrics


def demo(output, n=128, seed=42, backend="python"):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output = Path(output)
    p = Parameters()
    phases = ale_recipe()
    rom_solver, kmc_solver = mean_field, kmc
    if backend == "cpp":
        from .native import mean_field as rom_solver, kmc as kmc_solver
    elif backend != "python":
        raise ValueError("backend must be python or cpp")
    rom = rom_solver(phases, p)
    stochastic = kmc_solver(phases, p, sites=512, seed=seed)
    write_csv(output / "mean_field.csv", rom)
    write_csv(output / "kmc.csv", stochastic)
    rows = campaign(n, seed, backend=backend)
    write_csv(output / "campaign.csv", rows)
    metrics = train_surrogate(rows, output, seed)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    for data, label in [(rom, "Mean field"), (stochastic, "kMC (512 sites)")]:
        axes[0].plot(
            [r["time_s"] for r in data], [r["net_removed_nm"] for r in data], ".-", label=label
        )
    axes[0].set(xlabel="Time (s)", ylabel="Net removal (nm)")
    axes[0].legend()
    energies = np.linspace(0, 100, 101)
    for dose in [0.5, 2, 6]:
        values = []
        for energy in energies:
            trace = mean_field(ale_recipe(energy, dose))
            values.append(trace[-1]["net_removed_nm"] - trace[-5]["net_removed_nm"])
        axes[1].plot(energies, values, label=f"Modify {dose:g} s")
    axes[1].set(xlabel="Ion energy (eV)", ylabel="Cycle-5 net EPC (nm/cycle)")
    axes[1].legend()
    fig.suptitle("Illustrative silicon ALE model - uncalibrated synthetic results")
    fig.savefig(output / "overview.png", dpi=170)
    plt.close(fig)
    write_json(
        output / "manifest.json",
        manifest(
            dict(
                parameters=asdict(p),
                phases=[asdict(x) for x in phases],
                n=n,
                seed=seed,
                cycles=5,
                kmc_sites=512,
                backend=backend,
            )
        ),
    )
    return metrics


def calibrate(csv_path, output):
    """Fit ONLY sticking and chemical yield to compatible cycle-resolved ALE data."""
    with Path(csv_path).open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    required = FEATURES + ["cycle", "epc_nm", "sigma_nm", "source_id"]
    if len(rows) < 4 or any(any(not r.get(k) for k in required) for r in rows):
        raise ValueError(f"At least four complete rows required: {required}")
    values = np.array(
        [[float(r[k]) for k in FEATURES + ["cycle", "epc_nm", "sigma_nm"]] for r in rows]
    )
    if not np.isfinite(values).all() or (values[:, -1] <= 0).any():
        raise ValueError("Finite values and positive measurement uncertainties required")
    if any(float(r["cycle"]) < 1 or not float(r["cycle"]).is_integer() for r in rows):
        raise ValueError("cycle must be a positive integer")

    def residual(x):
        p = replace(Parameters(), sticking=x[0], chemical_yield_scale=x[1])
        predictions = []
        for row in rows:
            trace = mean_field(
                ale_recipe(*(float(row[f]) for f in FEATURES)), p, int(float(row["cycle"]))
            )
            previous = trace[-5]["net_removed_nm"] if len(trace) > 4 else 0.0
            predictions.append(trace[-1]["net_removed_nm"] - previous)
        return (np.array(predictions) - values[:, -2]) / values[:, -1]

    result = least_squares(residual, [0.25, 0.8], bounds=([0, 0], [1, 5]))
    report = dict(
        sticking=float(result.x[0]),
        chemical_yield_scale=float(result.x[1]),
        success=bool(result.success),
        weighted_residual_sum_squares=float(2 * result.cost),
        jacobian_rank=int(np.linalg.matrix_rank(result.jac)),
        jacobian_singular_values=np.linalg.svd(result.jac, compute_uv=False).tolist(),
        sources=sorted({r["source_id"] for r in rows}),
        limitation="In-sample fit only. Check identifiability, residuals, and independent experiments before predicting.",
    )
    import hashlib

    report["input_sha256"] = hashlib.sha256(Path(csv_path).read_bytes()).hexdigest()
    write_json(output, report)
    return report
