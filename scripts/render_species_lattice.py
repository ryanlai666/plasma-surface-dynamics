"""State-colored kMC rendering; all frames come from recorded SSA snapshots."""

from pathlib import Path
import io, json, hashlib
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from PIL import Image
from draw_species_network import COLORS, HF, REMOVED


def render_lattice(network, result, out):
    states = network['states']
    n = result['lattice'].shape[1]
    side = int(np.sqrt(n))
    assert side * side == n
    x, y = np.meshgrid(np.arange(side), np.arange(side))
    frames = []
    for k, t in enumerate(result['times']):
        grid = result['lattice'][k]
        removed = np.array([states[i]['Si_removed'] for i in grid])
        occupied = np.array([states[i]['HF_complex'] for i in grid])
        stage_colors = [COLORS[min(states[i]['fluorination_stage'], 3)] for i in grid]
        colors = [
            REMOVED if r else HF if h else c for r, h, c in zip(removed, occupied, stage_colors)
        ]
        fig = plt.figure(figsize=(11, 5.5), dpi=110)
        ax = fig.add_subplot(121)
        view = fig.add_subplot(122, projection='3d')
        ax.scatter(
            x.ravel(),
            y.ravel(),
            c=colors,
            edgecolors=stage_colors,
            s=82,
            linewidths=1.1,
            marker='s',
        )
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title('Top view: actual site states', fontsize=12)
        z = np.where(removed, 0.0, 1.0)
        view.scatter(x.ravel(), y.ravel(), z, c=colors, s=35, depthshade=False)
        view.view_init(28, -58)
        view.set_zlim(0, 2)
        view.axis('off')
        view.set_title('Perspective: released sites lowered', fontsize=12)
        phase = 'HF exposure' if t < 2 else 'Purge: HF arrival off'
        fig.suptitle(f'HF / SiN:H | 400 K | {t:.2f} s | {phase}', fontsize=16, weight='bold')
        handles = [Patch(fc=c, label=f'F{i}') for i, c in enumerate(COLORS)] + [
            Patch(fc=HF, label='HF complex'),
            Patch(fc=REMOVED, label='Si released'),
        ]
        fig.legend(
            handles=handles,
            loc='lower center',
            bbox_to_anchor=(0.5, 0.10),
            ncol=6,
            frameon=False,
            fontsize=10,
        )
        fig.text(
            0.5,
            0.065,
            f'HF complexes: {occupied.sum()}/{n}   |   Si released: {removed.sum()}/{n}   |   Gold fill = adsorbed HF; border = underlying F stage',
            ha='center',
            fontsize=9,
        )
        fig.text(
            0.5,
            0.025,
            'Schematic motif grid; heights are state indicators, not a crystal structure or physical film thickness.',
            ha='center',
            fontsize=9,
            color='#586779',
        )
        fig.subplots_adjust(left=0.02, right=0.98, bottom=0.21, top=0.84)
        b = io.BytesIO()
        fig.savefig(b, format='png')
        b.seek(0)
        frames.append(Image.open(b).convert('RGB'))
        plt.close(fig)
    frames[0].save(
        out / 'species_kmc.gif', save_all=True, append_images=frames[1:], duration=250, loop=0
    )
    frames[len(frames) // 2].save(out / 'species_kmc_preview.png')
    np.savez_compressed(
        out / 'lattice_trajectory.npz',
        times=result['times'],
        states=result['lattice'],
        counts=result['counts'],
        gas_counts=result['gas_counts'],
    )
    (out / 'animation_manifest.json').write_text(
        json.dumps(
            dict(
                frame_count=len(frames),
                interpolated_frames=0,
                representation='Actual kMC site states on a schematic grid; no atomistic geometry or length scale',
                seed=741,
                state_colors=dict(
                    F0=COLORS[0],
                    F1=COLORS[1],
                    F2=COLORS[2],
                    F3=COLORS[3],
                    HF_complex=HF,
                    Si_released=REMOVED,
                ),
                renderer='scripts/render_species_lattice.py',
                renderer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                trajectory_sha256=hashlib.sha256(
                    (out / 'lattice_trajectory.npz').read_bytes()
                ).hexdigest(),
                gif_sha256=hashlib.sha256((out / 'species_kmc.gif').read_bytes()).hexdigest(),
            ),
            indent=2,
        )
        + '\n'
    )


if __name__ == '__main__':
    root = Path(__file__).resolve().parents[1]
    out = root / 'docs/species_kmc_results'
    d = np.load(out / 'lattice_trajectory.npz')
    render_lattice(
        json.loads((root / 'configs/species_kmc_network.json').read_text()),
        dict(times=d['times'], lattice=d['states'], counts=d['counts'], gas_counts=d['gas_counts']),
        out,
    )
