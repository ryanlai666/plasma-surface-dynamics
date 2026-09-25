"""Render only evaluated NEB images, with the energy of each displayed geometry."""

from pathlib import Path
import csv, io, json, hashlib
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from ase.io import read

ROOT = Path(__file__).resolve().parents[1]


def main():
    summaries = []
    for folder in sorted((ROOT / 'data/surface_paths').glob('*/*/mace_neb')):
        if not (folder / 'summary.json').exists():
            continue
        meta = json.loads((folder / 'summary.json').read_text())
        atoms = read(folder / 'images.extxyz', ':')
        with (folder / 'energies.csv').open() as f:
            rows = list(csv.DictReader(f))
        energies = np.array([float(r['relative_energy_eV']) for r in rows])
        arc = np.array([float(r['arc_length_A']) for r in rows])
        assert len(atoms) == len(rows) == meta['images']
        frames = []
        xyz = np.array([a.positions for a in atoms])
        lo = xyz.min(axis=(0, 1))
        hi = xyz.max(axis=(0, 1))
        span = max(hi - lo)
        for n, a in enumerate(atoms):
            assert abs(a.get_potential_energy() - float(rows[n]['energy_eV'])) < 1e-8
            fig = plt.figure(figsize=(11, 5.5), dpi=115)
            ax = fig.add_subplot(121, projection='3d')
            ep = fig.add_subplot(122)
            x = a.positions
            colors = ['#8da6ba'] * (len(a) - 1) + [
                '#14a778' if meta['species'] == 'F' else '#d19023'
            ]
            ax.scatter(
                *x.T,
                s=[48] * (len(a) - 1) + [155],
                c=colors,
                depthshade=True,
                edgecolors='#32465a',
                linewidths=0.4,
            )
            for i in range(len(a)):
                for j in range(i):
                    d = np.linalg.norm(x[i] - x[j])
                    cutoff = 2.65 if i < len(a) - 1 else (2.1 if meta['species'] == 'F' else 2.65)
                    if 0.6 < d < cutoff:
                        ax.plot(*x[[i, j]].T, c='#758799', lw=1, alpha=0.65)
            center = (lo + hi) / 2
            for setter, c in zip([ax.set_xlim, ax.set_ylim, ax.set_zlim], center):
                setter(c - span * 0.55, c + span * 0.55)
            ax.set_box_aspect([1, 1, 1])
            ax.view_init(elev=25, azim=-55)
            ax.set_axis_off()
            label = (
                'IS'
                if n == 0
                else 'FS'
                if n == len(atoms) - 1
                else 'candidate TS'
                if n == meta['peak_image']
                else f'image {n}'
            )
            ax.set_title(f'{meta["species"]} on {meta["surface"]}: {label}', fontsize=13)
            ep.plot(arc, energies, 'o-', c='#397da2', lw=1, ms=5, label='Evaluated NEB images')
            ep.scatter([arc[n]], [energies[n]], s=105, c='#d44f49', zorder=5)
            ep.annotate(
                f'peak {max(energies):.3f} eV',
                (arc[meta['peak_image']], max(energies)),
                xytext=(0, 14),
                textcoords='offset points',
                ha='center',
                fontsize=10,
            )
            ep.set(
                xlabel='Configuration-space path length (angstrom)',
                ylabel='Energy relative to IS (eV)',
            )
            ep.set_title(f'Image {n}: {energies[n]:+.4f} eV', fontsize=13)
            ep.margins(y=0.22)
            ep.grid(alpha=0.18)
            fig.suptitle('MACE-MP-0b2 / rigid-surface climbing-image NEB', fontsize=15)
            fig.text(
                0.5,
                0.025,
                'Actual optimized images only | lines join evaluated energies | no motion interpolation',
                ha='center',
                fontsize=10,
            )
            fig.subplots_adjust(left=0.01, right=0.97, bottom=0.16, top=0.84, wspace=0.03)
            bio = io.BytesIO()
            fig.savefig(bio, format='png')
            bio.seek(0)
            frame = Image.open(bio).convert('RGB')
            frames.append(frame)
            plt.close(fig)
        frames[0].save(
            folder / 'path.gif',
            save_all=True,
            append_images=frames[1:],
            duration=[
                1100 if n in (0, meta['peak_image'], len(frames) - 1) else 650
                for n in range(len(frames))
            ],
            loop=0,
        )
        frames[meta['peak_image']].save(folder / 'peak.png')
        (folder / 'render_manifest.json').write_text(
            json.dumps(
                dict(
                    frame_count=len(frames),
                    frame_image_indices=list(range(len(frames))),
                    interpolated_frames=0,
                    energies_sha256=hashlib.sha256(
                        (folder / 'energies.csv').read_bytes()
                    ).hexdigest(),
                    coordinates_sha256=hashlib.sha256(
                        (folder / 'images.extxyz').read_bytes()
                    ).hexdigest(),
                    gif_sha256=hashlib.sha256((folder / 'path.gif').read_bytes()).hexdigest(),
                ),
                indent=2,
            )
            + '\n'
        )
        summaries.append(meta)
        (folder / 'README.md').write_text(
            f"# {meta['surface']} / {meta['species']} migration\n\n![Calculated path](path.gif)\n\n- Method: {meta['method']}\n- Numerical status: {meta['status']}; maximum NEB force {meta['max_neb_force_eV_A']:.5f} eV/angstrom.\n- Peak above IS: {meta['peak_above_IS_eV']:.6f} eV; FS minus IS: {meta['FS_minus_IS_eV']:.6f} eV.\n- [Every image and energy](energies.csv), [coordinates and forces](images.extxyz), [full provenance](summary.json).\n- Animation has {len(frames)} evaluated images, with no interpolated frames.\n\n## Model assumptions\n\n"
            + '\n'.join('- ' + x for x in meta['assumptions'])
            + '\n',
            encoding='utf8',
        )
    fig, axes = plt.subplots(1, len(summaries), figsize=(12, 3.8), layout='constrained')
    for ax, meta in zip(np.atleast_1d(axes), summaries):
        folder = ROOT / 'data/surface_paths' / meta['surface'] / meta['species'] / 'mace_neb'
        with (folder / 'energies.csv').open() as f:
            rows = list(csv.DictReader(f))
        arc = [float(r['arc_length_A']) for r in rows]
        e = [float(r['relative_energy_eV']) for r in rows]
        second = json.loads((folder / 'second_model.json').read_text())
        ax.plot(arc, e, 'o-', ms=4, label='MP-0b2 optimized path')
        ax.plot(arc, second['relative_energies_eV'], 's--', ms=3, label='MPA-0 on same images')
        ax.set(
            title=f"{meta['species']} / {meta['surface']}",
            xlabel='Path length (angstrom)',
            ylabel='Energy relative to own IS (eV)',
        )
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
    fig.suptitle('Rigid-surface migration: evaluated energies, not DFT validation', fontsize=13)
    fig.savefig(ROOT / 'docs/dry_etch_results/calculated_path_comparison.png', dpi=180)
    plt.close(fig)
    print('Rendered', len(summaries), 'calculated paths')


if __name__ == '__main__':
    main()
