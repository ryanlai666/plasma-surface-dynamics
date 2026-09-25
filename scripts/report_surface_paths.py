"""Build README and path-report tables from checked calculation outputs."""

from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]


def main():
    rows = []
    galleries = []
    summary_rows = []
    for p in sorted((ROOT / 'data/surface_paths').glob('*/*/mace_neb/summary.json')):
        s = json.loads(p.read_text())
        v = json.loads(p.with_name('curvature.json').read_text())
        m = json.loads(p.with_name('second_model.json').read_text())
        rel = p.parent.relative_to(ROOT).as_posix()
        state = (
            'constrained saddle checked'
            if v['eligible_constrained_saddle']
            else 'candidate only; see checks'
        )
        rows.append(
            f"| {s['surface']} / {s['species']} | {s['reaction']} | 0 | {s['peak_above_IS_eV']:.4f} | {s['FS_minus_IS_eV']:+.6f} | {s['max_neb_force_eV_A']:.4f} | {state} |"
        )
        summary_rows.append(
            f"| {s['surface']} / {s['species']} migration | {s['peak_above_IS_eV']:.4f} | MACE CI-NEB; {s['images']} evaluated images; {state} |"
        )
        galleries.append(
            f"### {s['surface']} / {s['species']}\n\n![{s['species']} surface migration with calculated image energies](../../{rel}/path.gif)\n\n[Coordinates, energies and checks](../../{rel}/README.md). MPA-0 single-point sampled peak on this path: **{m['sampled_peak_above_IS_eV']:.4f} eV**; this is a model-sensitivity result, not an independent TS calculation.\n"
        )
    if len(rows) != 3:
        raise ValueError('Three completed path folders required')
    table = (
        '| Surface / species | Reaction | IS (eV) | Peak above IS (eV) | FS minus IS (eV) | Max NEB force (eV/A) | Validation |\n|---|---|---:|---:|---:|---:|---|\n'
        + '\n'.join(rows)
    )
    p = ROOT / 'docs/dry_etch_results/SURFACE_PATHS.md'
    s = p.read_text(encoding='utf8')
    start = '<!-- CALCULATED_TABLE -->'
    end = '<!-- END CALCULATED_TABLE -->'
    content = start + '\n\n' + table + '\n\n' + '\n'.join(galleries) + '\n' + end
    if end in s:
        s = s[: s.index(start)] + content + s[s.index(end) + len(end) :]
    else:
        s = s.replace(start, content)
    p.write_text(s, encoding='utf8')
    p = ROOT / 'README.md'
    s = p.read_text(encoding='utf8')
    a = (
        s.index('### Surface paths: calculated images and energies')
        if '### Surface paths: calculated images and energies' in s
        else s.index('### Surface IS -> TS -> FS animation')
    )
    b = s.index('### Literature-parameterized dry-etch kinetics', a)
    replacement = (
        r"""### Surface paths: calculated images and energies

Each new GIF shows **force-optimized, energy-evaluated NEB images** (9 for each F path; 17 for Cl). There are no geometrically interpolated animation frames. Separate folders hold each surface/reactant combination, its coordinates, pointwise energies, checks and source information.

| Surface / reactant and event | Peak above reactant (eV) | Evidence |
|---|---:|---|
"""
        + '\n'.join(summary_rows)
        + r"""
| Fluorinated Si3N4 / HF + H2O removal | 0.6775 | Published cluster DFT, R6 |
| Fluorinated SiO2 / HF + H2O removal | 0.6376 | Published cluster DFT, R6 |
| Fluorinated Si3N4 / HF + HF(v=1) | 0.1513 effective | DFT barrier with source vibrational-energy adjustment |
| Fluorinated SiO2 / HF + HF(v=1) | 0.6946 | Source model assumes excitation quenching |
| Reconstructed Si(100) / SiCl4 recombination, OD | 2.4024 | Published DFT stationary points only |

The calculated migration paths use MACE-MP-0b2 on **ideal rigid silicon slabs**. Curvature checks apply only to the mobile halogen; these are not validated DFT barriers or complete ALE cycles. A second MACE model evaluates every image for model sensitivity. The HF-assisted values come from [Jung et al. (2020)](https://doi.org/10.1116/1.5125569); missing FS energies and raw paths remain unfilled.

![Three calculated surface energy profiles](docs/dry_etch_results/calculated_path_comparison.png)

![Calculated F migration on Si100](data/surface_paths/Si100/F/mace_neb/path.gif)

<details>
<summary>Cl on Si(100), F on Si(111), and published SiCl4 stationary points</summary>

![Calculated Cl migration on Si100](data/surface_paths/Si100/Cl/mace_neb/path.gif)

![Calculated F migration on Si111](data/surface_paths/Si111/F/mace_neb/path.gif)

![Three published SiCl4 stationary points](data/surface_paths/Si100_c4x2/SiCl4/published/od.gif)

The SiCl4 slideshow contains only the published IS, TS and FS. It illustrates SiCl3* + Cl* -> SiCl4(physisorbed), a recombination channel; no connecting trajectory or substrate-Si removal is claimed.

</details>

[Full curves, TS table, reaction/rate equations and hypotheses](docs/dry_etch_results/SURFACE_PATHS.md) | [Surface/reactant folder index](data/surface_paths/README.md) | [All six published stationary-point triplets](docs/dry_etch_results/TRANSITION_STATES.md).

"""
    )
    s = s[:a] + replacement + s[b:]
    p.write_text(s, encoding='utf8')


if __name__ == '__main__':
    main()
