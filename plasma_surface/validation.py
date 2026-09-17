"""Reproducible local simulation/verification campaign (no experimental claims).

Run from the project root: python -m plasma_surface.validation
"""
from dataclasses import asdict, replace
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import uuid

import numpy as np

from . import model, native
from .benchmark import benchmark
from .datasets import EXPECTED_MD5, inventory
from .io import manifest, write_csv, write_json
from .materials import compare
from .workflows import campaign, demo, train_surrogate


def run(output="docs/results"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    root = Path(__file__).resolve().parents[1]
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    print("Running regression tests", flush=True)
    tests = subprocess.run([sys.executable, "-m", "pytest", "-q",
                            "-p", "no:cacheprovider",
                            f"--basetemp=outputs/pytest-campaign-{uuid.uuid4().hex}"], cwd=root,
                           capture_output=True, text=True)
    (output/"tests.txt").write_text(tests.stdout+tests.stderr, encoding="utf-8")
    if tests.returncode:
        raise RuntimeError(f"Tests failed; see {output/'tests.txt'}")
    if not native.library_path().exists():
        raise RuntimeError("Build C++ first: python scripts/build_cpp.py")

    print("Running Python and C++ demos and SiNx scenarios", flush=True)
    demo(output/"python", backend="python")
    demo(output/"cpp", backend="cpp")
    compare(root/"configs/materials.json", output/"materials", "cpp")

    print("Evaluating a 4096-recipe campaign and surrogate", flush=True)
    rows = campaign(4096, seed=42, backend="cpp")
    reference = campaign(4096, seed=42, backend="python")
    write_csv(output/"campaign.csv", rows)
    max_backend_error = max(abs(a["epc_nm"]-b["epc_nm"]) for a, b in zip(rows, reference))
    metrics = train_surrogate(rows, output/"ml", seed=42)

    print("Checking finite-site stochastic ensembles", flush=True)
    phases = model.ale_recipe()
    target = model.mean_field(phases, cycles=5)[-1]["net_removed_nm"]
    convergence, samples = [], []
    for backend, solver, sizes, repeats in [
        ("cpp", native.kmc, [128, 512, 2048, 8192], 128),
        ("python", model.kmc, [128, 512], 32),
    ]:
        for sites in sizes:
            values = []
            for seed in range(repeats):
                value = solver(phases, cycles=5, sites=sites, seed=seed)[-1]["net_removed_nm"]
                values.append(value)
                samples.append(dict(backend=backend, sites=sites, seed=seed, net_removed_nm=value))
            std = float(np.std(values, ddof=1))
            mean = float(np.mean(values))
            se = std/np.sqrt(repeats)
            convergence.append(dict(backend=backend, sites=sites, repeats=repeats,
                                    mean_nm=mean, std_nm=std, standard_error_nm=float(se),
                                    reference_nm=target, absolute_bias_nm=abs(mean-target),
                                    within_five_se=bool(abs(mean-target) <= 5*se+1e-12)))
            print(f"  {backend}, {sites} sites: bias={abs(mean-target):.3g} nm, SE={se:.3g} nm", flush=True)
    write_csv(output/"convergence.csv", convergence)
    write_csv(output/"ensemble_samples.csv", samples)

    # Mechanism controls: no chemical flux vs no ions, isolated from cycle history.
    p = model.Parameters()
    def epc(phases):
        history = native.mean_field(phases, p, cycles=10)
        return history[-1]["net_removed_nm"]-history[-1-len(phases)]["net_removed_nm"]
    controls = []
    for energy in [10., 28., 75.]:
        recipe = model.ale_recipe(energy_ev=energy)
        for label, sequence in [
            ("combined", recipe),
            ("ions_only", [replace(s, radical_flux_m2_s=0) for s in recipe]),
            ("chemistry_only", [replace(s, ion_flux_m2_s=0) for s in recipe]),
        ]:
            controls.append(dict(energy_ev=energy, control=label, cycle10_epc_nm=epc(sequence)))
    write_csv(output/"half_cycle_controls.csv", controls)
    saturation = [dict(dose_s=float(dose), cycle10_epc_nm=epc(model.ale_recipe(28, float(dose), 20)))
                  for dose in [0, .1, .25, .5, 1, 2, 4, 8, 16, 32]]
    write_csv(output/"dose_saturation.csv", saturation)

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
    cpp = [r for r in convergence if r["backend"] == "cpp"]
    sizes = np.array([r["sites"] for r in cpp])
    deviations = np.array([r["std_nm"] for r in cpp])
    axes[0].loglog(sizes, deviations, "o-", label="C++ kMC ensemble SD")
    axes[0].loglog(sizes, deviations[0]*np.sqrt(sizes[0]/sizes), "--", label="N^(-1/2) guide")
    axes[0].set(xlabel="Number of sites", ylabel="SD of 5-cycle removal (nm)")
    axes[0].legend()
    axes[1].plot([r["dose_s"] for r in saturation], [r["cycle10_epc_nm"] for r in saturation], "o-")
    axes[1].set(xlabel="Modification duration (s)", ylabel="Cycle-10 EPC (nm/cycle)", title="28 eV; removal duration 20 s")
    figure.suptitle("Numerical verification of an uncalibrated surface model")
    figure.savefig(output/"verification.png", dpi=170)
    plt.close(figure)

    print("Benchmarking both backends", flush=True)
    timing = benchmark(output/"benchmark.json", repeats=5)
    public = dict(status="not_present", experimentally_validated=False)
    archive = root/"data/raw/si_hcl/data.zip"
    if archive.exists():
        if hashlib.md5(archive.read_bytes()).hexdigest() != EXPECTED_MD5:
            raise RuntimeError("Public data checksum mismatch")
        public = inventory(archive, output/"public_structure_inventory.csv")
        public.update(status="checksum_verified", doi="10.5281/zenodo.10211009",
                      experimentally_validated=False, note="Geometry inventory is not etch-rate validation")

    checks = dict(backend_agreement=max_backend_error < 1e-12,
                  stochastic_ensembles=all(r["within_five_se"] for r in convergence),
                  chemistry_only_zero=all(r["cycle10_epc_nm"] == 0 for r in controls if r["control"] == "chemistry_only"),
                  below_threshold_zero=all(r["cycle10_epc_nm"] == 0 for r in controls if r["energy_ev"] == 10),
                  dose_saturates=abs(saturation[-1]["cycle10_epc_nm"]-saturation[-2]["cycle10_epc_nm"]) < 1e-5)
    report = dict(timestamp_utc=datetime.now(timezone.utc).isoformat(), checks=checks,
                  max_python_cpp_epc_difference_nm=max_backend_error,
                  reference_five_cycle_removal_nm=target, convergence=convergence,
                  surrogate=metrics, speedup=timing["wall_time_speedup"], public_data=public,
                  experimental_validation="Not performed: no compatible measured EPC dataset is available locally.")
    write_json(output/"summary.json", report)
    run_manifest = manifest(dict(backend="cpp", seed=42, samples=4096, parameters=asdict(p),
                                 phases=[asdict(s) for s in phases], reference_cycles=5,
                                 cpp_ensemble_repeats=128, python_ensemble_repeats=32))
    run_manifest["validation_driver_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    run_manifest["artifacts_sha256"] = {str(f.relative_to(output)):hashlib.sha256(f.read_bytes()).hexdigest()
                                        for f in sorted(output.rglob("*"))
                                        if f.is_file() and f != output/"manifest.json" and f != output/"REPORT.md"}
    write_json(output/"manifest.json", run_manifest)
    passed = all(checks.values())
    (output/"REPORT.md").write_text(f"""# Simulation and numerical verification report

Generated {report['timestamp_utc']} using `python -m plasma_surface.validation`.

**Numerical checks: {'PASS' if passed else 'FAIL'}.** This verifies implementation and internal consistency, not agreement with measured plasma etching.

- Regression test output: [tests.txt](tests.txt).
- 4,096 recipes evaluated by both backends; maximum EPC difference: {max_backend_error:.3g} nm/cycle.
- Six stochastic ensembles: all means within five standard errors of the exact mean-field expectation: {checks['stochastic_ensembles']}.
- Subthreshold/chemical-only controls and dose saturation: {all(checks[k] for k in ['chemistry_only_zero', 'below_threshold_zero', 'dose_saturates'])}.
- Synthetic surrogate held-out MAE: {metrics['mae_nm']:.5f} nm/cycle; constant baseline: {metrics['constant_baseline_mae_nm']:.5f}; high-energy holdout: {metrics['high_energy_holdout_mae_nm']:.5f}.
- Median kMC times: Python {timing['python']['median_seconds']:.6f} s, C++ {timing['cpp']['median_seconds']:.6f} s; workload-specific speedup {timing['wall_time_speedup']:.1f}×, including the improved site-selection algorithm.
- Public DFT archive: {public['status']}. Structures are not experimental etch-rate labels.

![Convergence and saturation](verification.png)

## Results to inspect

- [Python ALE figure](python/overview.png) and [C++ ALE figure](cpp/overview.png).
- [Si/SiNx hypothetical scenarios](materials/comparison.png).
- [Stochastic convergence](convergence.csv), [raw ensemble samples](ensemble_samples.csv), and [half-cycle controls](half_cycle_controls.csv).
- [Surrogate metrics](ml/surrogate_metrics.json) and [predictions](ml/surrogate_predictions.csv).
- [Benchmark](benchmark.json), [summary](summary.json), and [provenance](manifest.json).

SiNx composition cards remain hypothetical. Default rates are uncalibrated. Statistical tolerances were fixed at five standard errors before evaluating ensembles; no seeds were discarded. Finite-site scatter should scale approximately as N^(-1/2), but this is not proof of experimental accuracy. Public geometry checks do not validate reaction barriers. No experimental validation or HPC scaling was performed.
""", encoding="utf-8")
    if not passed:
        raise RuntimeError("Verification failed; inspect summary.json")
    print(f"PASS: {output/'REPORT.md'}", flush=True)


if __name__ == "__main__":
    run()
