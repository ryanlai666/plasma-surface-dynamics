# kMC surface-lattice animations

## Visualization choice

The implemented model has independent surface columns with bare/modified states and integer height offsets. A square **site map plus a fixed-row height cross-section** represents those variables directly. Tile colors encode state; a pink outline marks columns removed since the preceding sampled frame. No bonds or Si/N atom colors are added, because the kernel does not track them. The highlighted map row identifies exactly which columns appear in the cross-section.

A ball-and-stick crystal or randomly colored Si/N lattice would imply structure and species-resolved reactions that the current model does not contain. The square arrangement is a display layout: site adjacency has no effect on rates.

The design follows the general separation of simulation data and visual glyphs described in the [OVITO visual-element documentation](https://www.ovito.org/manual/reference/pipelines/visual_elements/index.html), using property-based state colors as in [OVITO particle visualization](https://ovito.org/manual/reference/pipelines/visual_elements/particles.html). [kmos](https://github.com/mhoffman/kmos) is a relevant future lattice-kMC framework if lateral interactions and explicit site types are added. No OVITO/kmos dependency is needed to view or reproduce these GIFs.

## Reproduce

```powershell
python -m plasma_surface.animate
```

Pillow (installed with the current plotting dependencies) writes GIFs without FFmpeg or a GPU. The command reads `configs/materials.json` and writes four top-view/cross-section animations and four 3D perspective animations, representative PNG frames, sampled metrics, and a manifest to `docs/animations/`.

| Setting | Value |
|---|---|
| Systems | Si, SiN0.8, SiN1.0, Si3N4 |
| Sites | 16 x 16 independent columns |
| Random seed | 42 for each scenario |
| ALE cycles | 3 |
| Modification / purge / removal / purge | 2 / 0.5 / 4 / 0.5 s |
| Ion energy during removal | 35 eV |
| Neutral and ion flux in their respective phases | 2e19 m^-2 s^-1 |
| Temperature | 300 K |
| Simulated duration | 21 s |
| Sampling | 121 equally spaced times, including endpoints |
| GIF playback | 100 ms per frame; 12.1 s per loop |

The same phase schedule and display scales apply to all materials. The 35 eV energy is above all four illustrative chemical thresholds and below all four physical thresholds, allowing modified-layer removal in every scenario without the physical sputtering channel. It is a demonstration choice, not an experimentally established operating window.

## State fidelity and tests

The optional observer in `model.kmc` receives copied states at phase starts, accepted events, and phase ends. It does not draw random numbers. Uniform sampling selects the last state at or before each sample time; it never interpolates positions or invents events. Multiple events may occur between displayed frames. Purges retain their true durations. At the loop boundary the trajectory restarts visibly from the initial state.

Tests verify identical phase summaries with/without the observer, immunity to callback mutation, nondecreasing timestamps, zero-flux sampling, and the column-height/material-count balance. The generator also checks seeded parity for each material, the final balance, GIF frame counts, and playback duration. Numerical verification does not establish experimental accuracy.

SiNx labels specify the scenario's bulk-composition metadata. They do not represent spatially resolved Si/N identities, preferential nitrogen loss, crystallographic orientation, or equilibrated amorphous material. All nitride kinetic parameters remain hypothetical; the inherited 0.136 nm increment is an illustrative conversion. The cross-section is an independent-column statistic, not a predicted trench profile or validated AFM roughness.

## 3D perspective view

Files ending in `_3d.gif` project the same sampled states as the original site-map movies. A fixed pinhole camera draws exposed column tops and visible side faces. Surface colors and pink removal outlines have the same meaning in both views. The height reference is projected with the columns; zero denotes the initial surface. There is no interpolation, smoothing, or camera motion. All scenarios use identical camera and height scales.

Lateral grid spacing is schematic and vertical coordinates are exaggerated 1.5x for readability. Solid sides visualize the column height offsets; they do not introduce bulk atom positions or a crystallographic lattice. Occlusion can hide sites in perspective, so use the paired top view to inspect every site. This is a visualization of independent-column kinetics, not a morphology prediction.

The manifest records the camera offset, focal distance in pixels, exaggeration, and hashes of both GIFs for each material. Both views contain 121 frames at 100 ms each and use the same CSV time series. Reproduce all views with the command above.
