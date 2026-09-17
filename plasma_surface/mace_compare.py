"""Optional CPU comparison of two public MACE checkpoints. No DFT accuracy claim.

Run with the separate MACE environment from the repository root:
    .venv-mace/Scripts/python -m plasma_surface.mace_compare
"""
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import time
import urllib.request

import numpy as np
from .io import write_csv, write_json

MODELS = {
    "mp_0b2_small": "https://github.com/ACEsuit/mace-foundations/releases/download/mace_mp_0b2/mace-small-density-agnesi-stress.model",
    "mpa_0_medium": "https://github.com/ACEsuit/mace-foundations/releases/download/mace_mpa_0/mace-mpa-0-medium.model",
}


def run(output="docs/mace_results"):
    import torch
    from ase.build import bulk
    from ase.io import read, write
    from mace.calculators import mace_mp
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    torch.set_num_threads(2)
    torch.set_num_interop_threads(1)
    root = Path(__file__).resolve().parents[1]
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    cache = root/".cache/mace"
    cache.mkdir(parents=True, exist_ok=True)
    cif = root/"data/raw/cod/2102550.cif"
    if not cif.exists():
        cif.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen("https://www.crystallography.net/cod/2102550.cif", timeout=60) as r:
            cif.write_bytes(r.read())
    nitride = read(cif)
    assert nitride.get_chemical_formula() == "N8Si6", nitride.get_chemical_formula()
    structures = {"Si_diamond": bulk("Si", "diamond", a=5.43, cubic=True),
                  "beta_Si3N4": nitride}
    nitrogen = [i for i, a in enumerate(nitride) if a.symbol == "N"]
    for vacancies in [2, 3]:
        atoms = nitride.copy()
        del atoms[nitrogen[-vacancies:]]
        structures[f"Si6N{8-vacancies}_unrelaxed_vacancies"] = atoms
    cases = []
    for name, atoms in structures.items():
        for strain in [-.03, 0., .03]:
            case = atoms.copy()
            case.set_cell(atoms.cell*(1+strain), scale_atoms=True)
            case.info.update(case_name=name, linear_strain=strain,
                             origin="COD_2102550_derived" if name != "Si_diamond" else "ASE_diamond_a_5.43_A")
            cases.append(case)
    write(output/"input_structures.extxyz", cases)
    (output/"2102550.cif").write_bytes(cif.read_bytes())
    predictions, checks, provenance, forces = [], [], {}, {}
    for name, url in MODELS.items():
        path = cache/(name+".model")
        if not path.exists():
            print(f"Downloading {name}", flush=True)
            with urllib.request.urlopen(url, timeout=120) as response:
                payload = response.read(200_000_001)
            if len(payload) > 200_000_000:
                raise ValueError("Checkpoint exceeds the local comparison download budget")
            path.write_bytes(payload)
        provenance[name] = dict(url=url, sha256=hashlib.sha256(path.read_bytes()).hexdigest(), bytes=path.stat().st_size,
                                license="MIT; see upstream model repository")
        print(f"Loading {name} on CPU", flush=True)
        calc = mace_mp(model=str(path), device="cpu", default_dtype="float64", dispersion=False)
        for i, case in enumerate(cases):
            atoms = case.copy()
            atoms.calc = calc
            start = time.perf_counter()
            energy = float(atoms.get_potential_energy())
            force = atoms.get_forces()
            elapsed = time.perf_counter()-start
            if not np.isfinite(force).all() or not np.isfinite(energy):
                raise ValueError("Nonfinite prediction")
            forces[name, i] = force.copy()
            symbols = atoms.get_chemical_symbols()
            predictions.append(dict(model=name, case_id=i, structure=case.info["case_name"],
                strain=case.info["linear_strain"], atoms=len(atoms), n_si_ratio=symbols.count("N")/symbols.count("Si"),
                energy_ev=energy, energy_ev_per_atom=energy/len(atoms),
                force_rms_ev_A=float(np.sqrt(np.mean(force**2))),
                max_force_ev_A=float(np.max(np.linalg.norm(force,axis=1))), seconds=elapsed))
            print(f"  {case.info['case_name']}, strain {case.info['linear_strain']:+.2f}: {energy/len(atoms):.6f} eV/atom", flush=True)
        # Central finite difference at a displaced beta-Si3N4 configuration.
        probe = nitride.copy()
        probe.positions[0, 0] += .04
        probe.calc = calc
        analytical = float(probe.get_forces()[0,0])
        step = 1e-4
        probe.positions[0,0] += step
        plus = probe.get_potential_energy()
        probe.positions[0,0] -= 2*step
        minus = probe.get_potential_energy()
        numerical = float(-(plus-minus)/(2*step))
        error = abs(analytical-numerical)
        checks.append(dict(model=name, analytic_force_ev_A=analytical, finite_difference_force_ev_A=numerical,
                           abs_difference_ev_A=error, tolerance_ev_A=1e-4, passed=error < 1e-4))
    disagreement = []
    for i, case in enumerate(cases):
        disagreement.append(dict(case_id=i, structure=case.info["case_name"], strain=case.info["linear_strain"],
            force_component_rmse_between_models_ev_A=float(np.sqrt(np.mean((forces["mp_0b2_small",i]-forces["mpa_0_medium",i])**2)))))
    write_csv(output/"predictions.csv", predictions)
    write_csv(output/"model_disagreement.csv", disagreement)
    write_json(output/"force_gradient_checks.json", checks)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    for ax, label in zip(axes.flat, structures):
        for name in MODELS:
            subset = [r for r in predictions if r["model"] == name and r["structure"] == label]
            zero = next(r["energy_ev_per_atom"] for r in subset if r["strain"] == 0)
            ax.plot([r["strain"]*100 for r in subset], [(r["energy_ev_per_atom"]-zero)*1000 for r in subset], "o-", label=name)
        ax.set(title=label, xlabel="Isotropic linear strain (%)", ylabel="E - E(unstrained) (meV/atom)")
        ax.legend(fontsize=8)
    fig.suptitle("MACE checkpoint comparison: structural probes, not etch validation")
    fig.savefig(output/"comparison.png",dpi=160)
    plt.close(fig)
    report = dict(timestamp_utc=datetime.now(timezone.utc).isoformat(), models=provenance,
                  packages={k:importlib.metadata.version(k) for k in ["mace-torch","torch","ase","e3nn","numpy"]},
                  device="cpu", threads=2, dtype="float64", cases=len(cases), predictions=len(predictions),
                  structure_source="https://www.crystallography.net/cod/2102550.html", structure_license="CC0",
                  structure_sha256=hashlib.sha256(cif.read_bytes()).hexdigest(),
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  all_force_gradient_checks_passed=all(c["passed"] for c in checks),
                  limitations=["No independent DFT energy/force labels: disagreement is not prediction error.",
                               "Nitrogen-vacancy structures are unrelaxed stress probes, not realistic amorphous SiNx films.",
                               "No collision trajectories, plasma chemistry, surface charging, or etch yields evaluated.",
                               "Energies referenced separately within each composition; no cross-composition stability ranking."])
    write_json(output/"manifest.json",report)
    (output/"REPORT.md").write_text(f"""# MACE force-field comparison

Evaluated {len(cases)} structural probes using two official checkpoints on CPU (24 energy/force evaluations).

- MACE-MP-0b2 small and MACE-MPA-0 medium.
- Diamond Si and public beta-Si3N4 (COD 2102550).
- Derived Si6N6 and Si6N5 cells with unrelaxed nitrogen vacancies.
- Each structure at -3%, 0%, and +3% isotropic linear strain.
- Energy-gradient finite-difference checks passed: {report['all_force_gradient_checks_passed']} (absolute tolerance 1e-4 eV/angstrom).

![Comparison](comparison.png)

[Predictions](predictions.csv), [force disagreement](model_disagreement.csv), [gradient checks](force_gradient_checks.json), and [checkpoint/environment provenance](manifest.json).

Energy offsets are removed separately for each structure and model. Do not use these plots to rank compositions thermodynamically. No reference DFT forces or energies are available for these probes, so model disagreement is not an accuracy score. The finite-difference check verifies energy/force consistency, not agreement with experiments. Vacancy cells are deliberately unrelaxed structural tests, not synthesized or equilibrated SiNx films. This comparison does not validate plasma-impact chemistry or ALE rates.

The public crystal data are CC0: du Boulay et al., Acta Crystallographica B 60 (2004), 388-405, [COD 2102550](https://www.crystallography.net/cod/2102550.html). Models and citations: [MACE foundation-model repository](https://github.com/ACEsuit/mace-foundations). Model checkpoints remain cached locally and are not committed.
""", encoding="utf-8")
    if not report["all_force_gradient_checks_passed"]:
        raise RuntimeError("Energy/force consistency check failed; inspect report")
    print(f"Completed: {output/'REPORT.md'}", flush=True)


if __name__ == "__main__":
    run()
