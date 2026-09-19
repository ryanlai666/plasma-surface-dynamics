# Public data register

Checked 2026-09-18. Source availability and reuse rights are separate questions. Preserve authors, version, license, retrieved date, original units, and checksum for every acquisition. The downloader stores the full source metadata.

| Source | Contents and planned use | Limits / acquisition |
|---|---|---|
| [Martinez i Diaz et al., Si–HCl DFT data](https://zenodo.org/records/10211009), DOI 10.5281/zenodo.10211009 | Minima and transition-state CONTCAR geometries; atomistic reaction starting points | CC BY 4.0; 74,637-byte ZIP downloaded and verified; 22 structures inventoried. No barrier/etch-rate labels inferred from geometry. HCl differs from the illustrative Cl₂/Ar-inspired cycle. |
| [SAIT MLFF Framework](https://github.com/SAITPublic/MLFF-Framework), [NeurIPS paper](https://papers.neurips.cc/paper_files/paper/2023/hash/a1859debfb3b59d094f3504d5ebb6c25-Abstract-Datasets_and_Benchmarks.html) | Public SiN DFT datasets in extended XYZ, including varied structures and compositions; useful for atomistic benchmarks | Provider links sampled and raw archives via Google Drive. Not downloaded in this starter. Review archive-specific license and size before redistribution. No halogen reaction coverage should be assumed. |
| [Hong et al., HF etching of amorphous Si₃N₄](https://pubs.acs.org/doi/10.1021/acsami.4c07949) | Published atomistic etching study; supporting information includes yield and coarse-grained reaction-probability material | Public SI inspected: six impact yields and seven state coefficients extracted with provenance; CC BY-NC 4.0 source. No checkpoint, training frames or matched TS coordinates provided. Not fitted to the generic cards. |
| [Vella and Graves, Si–Cl₂–Ar⁺ ALE study](https://figshare.com/articles/media/29218103) | MD/reduced-order context for energy and fluence dependence | Record is media, ~47 MB, CC BY-NC 4.0; do not assume it is a tabular training dataset. Not downloaded. |
| [Daly et al., fluorocarbon ICP diagnostics](https://zenodo.org/records/7704879) | Optical emission, images, and tool logs for plasma-diagnostics ML | ~48.8 GB split archive; defer on a local-machine first pass. Not direct EPC labels and not downloaded. Inspect license metadata before use. |
| [LXCat redistribution policy](https://master.lxcat.net/instructions/redistribution.php) | Future electron/ion collision inputs for plasma kinetics | Obtain directly from the provider with database-specific references; third-party redistribution is restricted. No LXCat files bundled. |

## Downloaded DFT archive

Credit: Biel Martinez i Diaz, Jing Li, Hector Prats, Benoit Sklénard (2023), *Data for: Atomistic description of Si etching with HCl*, Zenodo, DOI 10.5281/zenodo.10211009. License: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

Original archive MD5: `34644b3c55b44441dc64492759cb6d31`. Local raw data are excluded from source control by default. `record.json` contains provider metadata; `provenance.json` contains retrieval time, hashes, and creators; `structure_inventory.csv` adds derived member metadata without modifying structures. Recreate the inventory offline with:

```powershell
python -m plasma_surface.cli inventory data/raw/si_hcl/data.zip
```

## A local-computer SiN data path

Begin with the provider's sampled SiN archive rather than raw trajectories. Inspect available composition counts and reference energy/force conventions. Use a small deterministic subset for analysis; keep all snapshots from a trajectory together. Preserve the provider's out-of-distribution split and add composition-held-out evaluation. If labels lack reactive species, restrict the claim to Si/N structural energetics. A bulk cohesive energy cannot be substituted for an ion-assisted etching barrier.

For digitized literature curves, retain figure/table identifiers, digitizer settings, axis calibration, inferred uncertainty, process conditions, and license notes. Do not silently combine thermal ALE, continuous etching, ion beams, and plasma ALE into one target.

## Collected atomistic subset and surface pathways

The [atomistic acquisition guide](ATOMISTIC_DATA.md) now links a 70-frame public DFT energy/force subset, indexed Si-H-Cl geometries including 15 source-named transition states, and two cited Cl diffusion barriers. See the [dry-etch reaction report](dry_etch_results/REPORT.md) for published-rate reproduction and surface pathways and [candidate fragments and experiment sets](REACTION_CANDIDATES.md) for follow-up research.

The [dry-etch acquisition and parameter audit](DRY_ETCH_PARAMETERS.md) adds four F2/Si rate laws, all 16 SiN:H/HF fluorination pathways, salt energetics and 17 [SiCl4 surface structures](../data/reference/sicl4_surface/README.md) forming six matched IS/TS/FS paths.
