"""Explicit composition-tagged scenarios, never an inferred composition-rate law."""

from dataclasses import asdict, replace
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from .model import Parameters, ale_recipe, mean_field
from .io import write_csv, write_json, manifest


def compare(cards_path, output, backend="python"):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    solver = mean_field
    if backend == "cpp":
        from .native import mean_field as solver
    elif backend != "python":
        raise ValueError("backend must be python or cpp")
    config = json.loads(Path(cards_path).read_text(encoding="utf-8"))
    if not config.get("cards"):
        raise ValueError("At least one material card required")
    output = Path(output)
    rows, effective = [], []
    fig, ax = plt.subplots(figsize=(8, 5), layout="constrained")
    for card in config["cards"]:
        x, hydrogen = card["n_si_ratio"], card["hydrogen_atomic_fraction"]
        if not math.isfinite(x) or x < 0 or not math.isfinite(hydrogen) or not 0 <= hydrogen < 1:
            raise ValueError("Invalid composition metadata")
        p = replace(Parameters(), **card["parameters"])
        effective.append(dict(**card, effective_parameters=asdict(p)))
        values = []
        for energy in np.linspace(0, 100, 101):
            history = solver(ale_recipe(float(energy)), p)
            epc = history[-1]["net_removed_nm"] - history[-5]["net_removed_nm"]
            values.append(epc)
            rows.append(
                dict(
                    material=card["name"],
                    n_si_ratio=x,
                    hydrogen_atomic_fraction=hydrogen,
                    energy_ev=float(energy),
                    epc_nm=epc,
                    data_origin=config["status"],
                )
            )
        ax.plot(np.linspace(0, 100, 101), values, label=card["name"])
    ax.set(
        xlabel="Ion energy (eV)",
        ylabel="Cycle-5 net EPC (nm/cycle)",
        title="Si / SiNx hypothetical parameter sensitivity\nNot measured composition dependence",
    )
    ax.legend()
    write_csv(output / "comparison.csv", rows)
    fig.savefig(output / "comparison.png", dpi=170)
    plt.close(fig)
    write_json(
        output / "manifest.json",
        manifest(
            dict(
                cards=effective,
                status=config["status"],
                backend=backend,
                cards_sha256=hashlib.sha256(Path(cards_path).read_bytes()).hexdigest(),
            )
        ),
    )
    return dict(rows=len(rows), output=str(output), status=config["status"])
