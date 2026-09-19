# Collecting DFT, transition-state, and diffusion data

This collection separates **published DFT labels**, **published geometry-only structures**, **literature barriers**, **reaction-specific literature rates**, and **ML predictions**. None is silently substituted for a calibrated plasma etch rate. Source register checked 2026-09-18: [machine-readable catalog](../data/literature/atomistic_sources.json).

## What is now available

| Collection | Local artifact | Content and interpretation |
|---|---|---|
| Public Si/O/C/F DFT subset | [70 frames with energy and force labels](../data/reference/atomistic/dft_subset.extxyz), [index](../data/reference/atomistic/frame_index.json) | 45 evenly sampled frames from nine bulk/etching members plus all 25 configurations from six quasi-static drag groups |
| Public Si-H-Cl geometries | [22 geometries](../data/reference/atomistic/hcl_geometries.extxyz), [index](../data/reference/atomistic/hcl_geometry_index.json) | Seven adsorbates and 15 structures named TS by their authors; no source energy/force/Hessian labels in these CONTCAR files |
| Surface diffusion reference | [Two diffusion barriers](../data/literature/diffusion_barriers.csv) | Cl hopping and SiCl-complex diffusion on the specified Si(111)-(5x5) surface; no prefactors or path coordinates acquired |
| Existing SiN:H/HF reference | [Original two DFT pathway entries; expanded to [all 16](../data/literature/sin_hf_pathways.csv)](../data/literature/sin_hf_barriers.csv) | Chemistry-specific author-poster values; kept separate from chlorine chemistry |
| Dry-etch reaction parameters and paths | [Surface TS and rates](DRY_ETCH_PARAMETERS.md) | Four F2/Si rate laws, 16 SiN:H/HF fluorination pathways, six matched SiCl4 surface triplets, and relevant impact/salt data |

### Published labels and provenance

The Si/O/C/F source is Oh, You, Kim, Lee, Han and Kang, [Model and Data from: A Lightweight Universal Machine-Learning Interatomic Potential via Knowledge Distillation for Scalable Atomistic Simulations](https://zenodo.org/records/19491140), DOI 10.5281/zenodo.19491140, **CC BY 4.0**. The 89.7 MB `dft.tar` already existed in `MLIP_benchmark`; its provider MD5 was verified before reuse. The full archive contains 6,163 frames across ten extended-XYZ members; the repository retains a 70-frame starter subset and an inventory for expanding it. The source README identifies etching configurations as ML-generated trajectories subsequently labeled by DFT, not AIMD. Its `qsd` data are quasi-static drag configurations with `group` and `dist` metadata. Some sample distances are strongly repulsive: a scan's highest energy is **not** automatically a reaction barrier. We preserve source frame text, including `energy` and `free_energy`; the index uses `energy`. All 70 selected energy/force arrays are finite; the 45 non-QSD frames carry `converged=True`, while the 25 QSD frames do not supply that flag. Missing convergence metadata is not treated as proof of SCF convergence. Audit the full XC/pseudopotential settings before pooling these labels with other calculations.

The Si-H-Cl source is Martinez i Diaz, Li, Prats and Sklenard, [Data for: Atomistic description of Si etching with HCl](https://zenodo.org/records/10211009), DOI 10.5281/zenodo.10211009, **CC BY 4.0**. These are the previously downloaded structures, now converted into a convenient indexed ASE collection. They are not 22 newly computed DFT results. Use the associated [reaction-mechanism paper](https://doi.org/10.1016/j.apsusc.2024.159836) to map t1-t20 labels to reactions and consistent reactant reservoirs before assigning barriers. Do not calculate an activation energy by subtracting unmatched slabs with different atom counts.

[Provenance](../data/reference/atomistic/provenance.json) records original archive hashes, source attribution, selection rules, and output hashes. [Archive inventory](../data/reference/atomistic/archive_inventory.json) records the count and hash of every DFT trajectory member. Selected source frames are unchanged; HCl geometries were converted from VASP to extended XYZ with ASE. Upstream licenses and author credit apply to the data independently of project authorship.

## Reproduce the collection

Use an environment with ASE and NumPy, such as this project's `.venv-mace` or the existing benchmark environment. All outputs go into this repository; the sibling project is read only.

```powershell
.venv-mace\Scripts\python.exe scripts/collect_atomistic.py `
  --etch-archive ..\MLIP_benchmark\data\downloads\etch_dft.tar
```

On a new computer, obtain `dft.tar` directly from the [versioned Zenodo record](https://zenodo.org/records/19491140), and recreate the HCl archive with `python -m plasma_surface.cli fetch-data`. Pass the local DFT tar path with `--etch-archive`. The curator rejects archives whose published MD5 differs and reads regular members without extracting arbitrary paths.

For training, split whole trajectories and parent configurations before subsampling. Several QSD groups come from the same parent impact configuration; grouping only by individual drag path would still leak closely related structures. This subset is for smoke tests and method development, not a statistically representative training set or unbiased model ranking.

## Where to obtain more data

| Priority | Source | Concrete next acquisition | Scope and limits |
|---|---|---|---|
| 1 | [Si-HCl DFT/TS archive and paper](https://zenodo.org/records/10211009) | Recover source transition-state/endpoint correspondence and electronic energies, or recompute matched geometries with one protocol | Direct Si surface relevance; present archive has geometries only |
| 1 | [SiN:H + HF DFT study](https://doi.org/10.1016/j.apsusc.2024.159414), [public author poster](https://avssymposium.org/ALD2024/Sessions/SupplementalDocumentDownload/80689?sessionId=79078) | Expand reaction-specific barrier entries with exact energy references and acquire corresponding structures where available | Amorphous hydrogenated nitride; not transferable to a zero-H chlorine model |
| 2 | [SAIT SiN benchmark](https://github.com/SAITPublic/MLFF-Framework) | Begin with sampled `SiN.tar`; inventory N/Si, H content, split identifiers, energies and forces | Valuable composition-resolved structural DFT; do not assume halogen or transition-state coverage. Provider Drive link identified; archive not downloaded here |
| 2 | [Transition1x](https://www.nature.com/articles/s41597-022-01870-w), [author code](https://gitlab.com/matschreiner/Transition1x) | Select small complete reactant/path/TS/product cases for testing a generic TS workflow | Organic reaction-path benchmark; full HDF5 is about 6.6 GB. Audit overlap with OMol25/model training before calling it an independent test |
| 2 | [OMol25](https://fair-chem.github.io/omol25/) | Filter a small molecule/reaction subset with explicit element, charge, multiplicity and provenance requirements | Molecular DFT energy/force data at wB97M-V/def2-TZVPD; not a periodic semiconductor surface dataset. Follow Hugging Face access terms |
| 3 | [OMol25 electronic structures](https://fair-chem.github.io/omol25-elec/) | Acquire selected raw ORCA outputs for inspected molecule IDs through the provider | Useful for auditing calculations and wavefunctions; the full collection is far too large for this local project |

Public availability does not guarantee redistributable full-text PDFs or supplementary files. The collection includes attributed factual values and licensed geometry subsets; papers retain their own rights. Source summaries come from available publisher/author records, not a claim to have audited every full paper.

## Diffusion: what the published numbers mean

Asari, Nara and Ohno, [First-principles study of chlorine atom diffusion on Si(111)-(5x5) reconstructed surfaces](https://doi.org/10.1016/j.apsusc.2012.07.015), report **1.73 eV** for Cl hopping and **1.34 eV** for SiCl-complex diffusion in their [publisher abstract](https://www.sciencedirect.com/author/6603600410/jun-nara). These are different microscopic mechanisms. They do not justify using one universal diffusion barrier for Si(001), SiNx, different coverages, or charged surfaces. The numerical extract retains the source DOI, context and missing-prefactor status; coordinates and path verification remain unavailable here.

For a future lattice model, define directional hops with

$$
k_{i\to j}=\nu_{ij}\exp[-(E^\ddagger_{ij}-E_i)/(k_BT)].
$$

Use the same saddle reference for reverse hops and check detailed balance under equilibrium conditions. For equivalent uncorrelated hops of length a in d dimensions, with **total** escape rate Gamma summed over neighbors,

$$
D=\frac{\Gamma a^2}{2d}.
$$

Do not substitute a per-neighbor rate for Gamma without its coordination factor. Anisotropic pathways require a diffusion tensor; interacting lattices need correlation effects. The present independent-column ALE solver has no hops, so these values are reference data rather than active kMC parameters.

## Computing missing barriers locally

The local inventory found DPA-3.3-1M with an **OMol25 head**, not UMA. [The checkpoint's official model card](https://huggingface.co/deepmodelingcommunity/DPA-3.3-1M) specifies molecular task selection and `[charge, multiplicity]` inputs. Existing MACE-MP and NequIP materials checkpoints should not be relabeled as molecular OMol25 potentials. No existing benchmark environment or checkpoint was replaced.

Start with a defined dry-etch surface event and matched IS/TS/FS, not a gas-product conformational change. The current [parameter audit and source-specific calculations](DRY_ETCH_PARAMETERS.md) identify what can actually enter a kinetic model. Obtain a consistent cell, constraints, coverage and charge/spin state, relax matching endpoints, run CI-NEB and saddle refinement, then verify the mobile-subspace Hessian and connectivity. Check cluster-size/termination effects if a cluster is used. Compute vibrational corrections/prefactors before promoting electronic barriers to production rates.

A short AIMD trajectory does not replace a rare-event search. For SiNx, begin with several amorphous configurations at each N/Si and H fraction; a single crystal or isolated cluster cannot establish the material's barrier distribution. A full 177-181 atom HCl slab DFT/NEB calculation is substantially more expensive than a molecular calculation and has not been run here.
