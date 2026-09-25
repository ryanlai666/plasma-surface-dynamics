"""Shared Matplotlib style for presentation figures.

The colours are a validated categorical palette (checked for colour-vision
deficiency and normal-vision separation) plus one ordinal blue ramp for ordered
quantities such as temperature or fluorination stage.  Text always uses ink
colours; series colours only mark data.
"""

from __future__ import annotations

import matplotlib as mpl
import numpy as np

# Categorical slots, used in this order and never cycled.
BLUE = '#2a78d6'
ORANGE = '#eb6834'
AQUA = '#1baf7a'
YELLOW = '#eda100'
RED = '#e34948'
# Ordinal ramp: light -> dark, one hue (validated with --ordinal).
RAMP4 = ['#86b6ef', '#3987e5', '#1c5cab', '#0d366b']

SURFACE = '#fcfcfb'
INK = '#0b0b0b'
INK2 = '#52514e'
MUTED = '#898781'
GRID = '#e1e0d9'
AXIS = '#c3c2b7'
PHASE_FILL = '#f0efec'
EMPTY = '#ffffff'


def apply() -> None:
    mpl.rcParams.update(
        {
            'font.family': 'sans-serif',
            'font.sans-serif': ['Segoe UI', 'Helvetica Neue', 'Arial', 'DejaVu Sans'],
            'font.size': 10.5,
            'text.color': INK,
            'axes.labelcolor': INK2,
            'axes.edgecolor': AXIS,
            'axes.linewidth': 0.8,
            'axes.facecolor': SURFACE,
            'axes.grid': True,
            'axes.grid.axis': 'y',
            'axes.axisbelow': True,
            'axes.spines.top': False,
            'axes.spines.right': False,
            'axes.titlesize': 11.5,
            'axes.titleweight': 'bold',
            'axes.titlelocation': 'left',
            'axes.titlepad': 10,
            'grid.color': GRID,
            'grid.linewidth': 0.7,
            'xtick.color': MUTED,
            'ytick.color': MUTED,
            'xtick.labelcolor': INK2,
            'ytick.labelcolor': INK2,
            'lines.linewidth': 2.0,
            'lines.solid_capstyle': 'round',
            'legend.frameon': False,
            'figure.facecolor': SURFACE,
            'savefig.facecolor': SURFACE,
            'savefig.dpi': 150,
        }
    )


def header(fig, title: str, subtitle: str | None = None, top: float = 0.97) -> None:
    """Finding-as-title at top left, with an optional muted subtitle."""
    fig.text(0.012, top, title, fontsize=15, weight='bold', color=INK, va='top')
    if subtitle:
        fig.text(0.012, top - 0.055, subtitle, fontsize=10, color=INK2, va='top', wrap=True)


def footnote(fig, text: str) -> None:
    fig.text(0.012, 0.012, text, fontsize=8.5, color=MUTED, va='bottom')


def shade_phases(ax, spans, label_y: float = 1.0) -> None:
    """Shade (start, end, label) spans and name them above the plot area."""
    for i, (start, end, label) in enumerate(spans):
        if i % 2 == 0:
            ax.axvspan(start, end, color=PHASE_FILL, lw=0, zorder=0)
        ax.text(
            (start + end) / 2,
            label_y,
            label,
            transform=ax.get_xaxis_transform(),
            ha='center',
            va='bottom',
            fontsize=8.5,
            color=MUTED,
        )


def spread(values, min_gap: float) -> list[float]:
    """Nudge label positions apart so no two are closer than ``min_gap``."""
    order = np.argsort(values)
    placed = np.array(values, dtype=float)[order]
    for i in range(1, len(placed)):
        placed[i] = max(placed[i], placed[i - 1] + min_gap)
    shift = max(0.0, placed[-1] - max(values)) / 2 if len(placed) else 0.0
    placed -= shift
    out = np.empty_like(placed)
    out[order] = placed
    return out.tolist()


def end_labels(ax, x: float, ys, texts, colors, min_gap: float, dx: float = 0.0) -> None:
    """Direct labels at line ends: coloured tick mark, ink text."""
    for y0, y, text, color in zip(ys, spread(ys, min_gap), texts, colors):
        ax.annotate(
            text,
            xy=(x, y0),
            xytext=(x + dx, y),
            textcoords='data',
            va='center',
            ha='left',
            fontsize=9.5,
            color=INK,
            annotation_clip=False,
            arrowprops=dict(arrowstyle='-', color=color, lw=1.2, shrinkA=0, shrinkB=2),
        )


def bar_value_labels(ax, bars, fmt: str = '{:.0f}', pad: float = 3) -> None:
    for bar in bars:
        width = bar.get_width()
        ax.annotate(
            fmt.format(width),
            xy=(max(width, 0), bar.get_y() + bar.get_height() / 2),
            xytext=(pad, 0),
            textcoords='offset points',
            va='center',
            fontsize=9,
            color=INK2,
        )
