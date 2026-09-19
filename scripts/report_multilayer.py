"""Publish computed multilayer status and the rate-evidence work queue."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
def main():
    out=ROOT/'docs/multilayer_results';d=json.loads((out/'summary.json').read_text());q=json.loads((ROOT/'data/multilayer/rate_requests/index.json').read_text())['requests']
    text="""# Bond-resolved multilayer etching prototype

**Implemented and run; not validated predictive ALE.** This Python engine is separate from the original 45-state Python/C++ motif model. The new bond-graph engine has not been ported to C++. Rate and access assumptions must be validated before interpreting physical etch rates.

![Actual multilayer graph snapshots](multilayer_kmc.gif)

## What changed physically in the model

| Requirement | Implemented behavior | Remaining limitation |
|---|---|---|
| Subsurface connectivity | 336 Si/N substrate atoms, 552 initial Si-N bonds; public beta-Si3N4 crystal, six unit-cell depth bands, periodic lateral connectivity | Ideal crystalline positions stay fixed; no amorphous reconstruction or force relaxation |
| H/F/Cl termination | Every dangling valence has an explicit attached-species entry; every active Si has four bonds/ligands and every N three | Initial Si cap proportions are assumed; no fitted coverage, N-F chemistry or electronic charge assignment |
| Bond changes and removal | HF cleavage deletes a specific Si-N/Si-Si edge, attaches F to Si and H to its partner; gas release requires all substrate backbonds broken | Event set excludes ions, fluorocarbon films and salt retention |
| Exposure below the surface | Removal updates a moving local height envelope; deeper atoms become exposed or enter a finite access zone | 6x6 fractional-xy columns approximate accessibility; no pore diffusion, shadowing or transport solver |
| Depth-dependent modification | HF arrival attenuates with depth below the moving local envelope; buried atoms beyond the access cutoff do not adsorb HF | Penetration and attenuation lengths are assumptions, not measured values |
| Local-environment rates | Rules distinguish N hydrogenation, incorporated F, Si-H and Si-Cl exchange; exact environment evidence can override rates in strict mode | Literature seeds are unvalidated transfers; no qualified graph-specific records yet |
| Conservation | Si, N, H, F and Cl conserved with the signed gas reservoir after every event | Conservation does not establish correct kinetics |

Initial dangling N bonds are H capped. Dangling Si caps are sampled H/F/Cl with probabilities 0.2/0.5/0.3 and seed 71. The bottom unit-cell band is fixed. A band contains 56 substrate atoms and is **not one atomic monolayer**. Hatched supports in the separate motif drawings are schematic; the depth bands here come from explicit crystal coordinates.

For column c, the access depth is `d_i = z_top(c) - z_i`. Exposed nodes have d=0; an assumed penetrating reactant can reach d <= d_max. During a dose, `a_ads,i = 5 exp(-d_i / 1 angstrom) s^-1` in the demonstration; arrival is zero during purge. Local precursor loss is assumed 100 s^-1. The default strict policy admits only qualified rate records, so it does not silently use these hazards.

Si-N cleavage uses stage/motif values from the [existing HF/SiN:H source extract](../../data/literature/sin_hf_pathways.csv) only as **demonstration seeds**. These source environments differ from the crystal graph. Si-H exchange and Si-Si cleavage are also unvalidated transfers. Fourth Si-N cleavage is seeded by P4 only for an NH bridge; unmatched final bare-N or NH2 cleavage is assigned zero rate even in the demonstration, rather than borrowing a third-fluorination barrier. Si-Cl exchange (0.75 eV), detached Si-molecule release (0.53 eV), and NH3 release (0.45 eV) are explicitly assumed barriers. All use an assumed 1e12 s^-1 prefactor. In particular, molecular desorption is not mislabeled as the source P4 bond-cleavage event. No graph-specific barrier or local vibrational prefactor has been calculated here.

The network can emit stoichiometrically complete SiHhFfClc molecules with h+f+c=4 once no substrate bonds remain, plus NH3, H2 and HCl where applicable. Mixed-halide release is a demonstration channel, not established product yield. NHx retained on neighboring N is distinct from released NH3.

## Actual run and controls

Three 2 s HF doses plus 1 s purges were run at an assumed 450 K. This is pulsed reactive etching, not a demonstrated self-limiting ALE sequence. The finite fixed-bottom domain imposes another possible saturation mechanism.

"""
    text+=f"The baseline contains **{d['baseline_events']} events**, removes **{d['removed'].get('Si',0)} Si and {d['removed'].get('N',0)} N**, and exposes **{d['newly_exposed_atoms']} previously covered substrate atoms**. Element conservation is exact. These are simulated demonstration counts, not measurements.\n\n| Initial depth band | Removed / initial substrate atoms |\n|---|---:|\n"
    for band in range(5,-1,-1):text+=f"| {band} ({'initial top' if band==5 else 'fixed bottom' if band==0 else 'subsurface'}) | {d['removed_by_initial_layer'].get(str(band),0)} / {d['initial_layer_sizes'][str(band)]} |\n"
    text+="""
![Access-depth control runs](multilayer_controls.png)

Twelve controls vary the allowed access depth (0, 1, 2, 4 angstrom) over three random seeds. They test how strongly this modeling assumption affects removal, not a physical uncertainty interval. With the empty qualified-rate library, the `validated_only` control produces **zero events**.

[Run metadata and control results](summary.json), [full event history and snapshots](../../data/multilayer/demo_trajectory.json), [animation provenance](animation_manifest.json), [initial graph](../../data/multilayer/sin_graph.json).

## Improving rates instead of transferring one barrier everywhere

1. **Identify the actual environment.** Each event gets a two-shell descriptor containing reacting-atom roles, Si/N connectivity, H/F/Cl counts, HF occupancy, substrate bond lengths and temperature. This is a candidate-grouping key, not a unique relaxed adsorbate geometry or electronic state.
2. **Collect IS/FS connectivity requests.** Record both local graphs and gas stoichiometry, including candidate events that never occurred. Construct and relax ligands, precursor placement, boundary caps and electronic states before running atomistic calculations. The requests are not ready-made DFT inputs.
3. **Compute matched evidence.** For chemistry, require stable IS/FS, a force-converged index-one saddle, connected endpoints, basis/cell and charge/spin checks, and prefactors/free energies with consistent references. Adsorption instead needs flux, sticking and site-area evidence. Experimental rates need matching film/recipe and held-out validation.
4. **Store applicability and uncertainty.** A qualified record must match the environment and operating context, identify its source, reference local artifacts with verified SHA256 hashes, document an applicability review and supply a rate interval with its basis. The software checks these records; it cannot substitute for the underlying scientific review.
5. **Select expensive work efficiently.** The current queue prioritizes missing final-cleavage barriers, then ranks observed frequency plus candidate presence. This is a coverage heuristic, not rate control. Combine it with the independent motif-model sensitivity results and model disagreement before allocating DFT/NEB calculations. Split future ML fits by whole environment/reaction, not adjacent path frames.
6. **Calibrate transfer before interpolation.** Benchmark small representative motifs with one consistent DFT method, compute barrier/reference corrections against ML paths, and validate on withheld neighboring environments. Extrapolation or an unreviewed descriptor match stays disabled; no universal barrier offsets are installed.

"""
    text+=f"The saved campaign produced **{len(q)} distinct missing-rate environments**, including unexecuted candidates. Twelve high-priority requests have individual folders; the full index retains all cases. A large count is evidence of unresolved parameterization, not validated mechanism completeness.\n\n| Rank | Event | Observed executions | Sampled candidate occurrences | Request |\n|---|---|---:|---:|---|\n"
    for r in q[:12]:text+=f"| {r['priority_rank']} | {r['event_kind']} | {r['observed_events']} | {r['candidate_snapshots']} | [IS/FS connectivity](../../data/multilayer/rate_requests/{r['key'][:16]}/request.json) |\n"
    text+="""
[All missing environments](../../data/multilayer/rate_requests/index.json) | [Empty qualified library](../../configs/multilayer_rate_library.json) | [Evidence gate](../../plasma_surface/rate_evidence.py) | [Broader literature and scientific critique](../SCIENTIFIC_REVIEW.md).

The present source barriers describe hydrogenated amorphous nitride, whereas this graph starts from a crystal. This material mismatch alone prevents declaring the transferred rates validated. Coadsorbed HF/H2O, salt formation, ion-created defects and charge/spin-dependent paths require additional, separately qualified branches.

## Reproduction

```sh
# ASE environment needed only to reconstruct the public-CIF graph:
python scripts/build_multilayer_graph.py
# Explicit demonstration run plus strict-policy and access-depth controls:
python scripts/run_multilayer.py
python scripts/queue_multilayer_rates.py
python scripts/render_multilayer.py
python scripts/report_multilayer.py
python -m pytest -q
```

Library calls default to `policy='validated_only'`; use `policy='demonstration'` explicitly to reproduce the unvalidated mechanism exercise. Do not convert these atom counts into EPC or claim an ALE window. The remaining work is environment-specific energetics and matched validation, not additional animation layers.
"""
    (out/'REPORT.md').write_text(text,encoding='utf-8');print('Reported multilayer run and',len(q),'rate gaps')
if __name__=='__main__':main()
