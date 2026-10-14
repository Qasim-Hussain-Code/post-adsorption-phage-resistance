"""Figure style for the manuscript figures.

The palette is the validated default from the project's visualisation
guidelines. Both the three-slot categorical set and the five-step ordinal blue
ramp used for the plaque score were checked with the palette validator against
the light chart surface before being adopted:

  categorical #2a78d6, #eb6834, #1baf7a
    lightness band, chroma floor, all-pairs colour vision deficiency
    separation (worst 9.2) and normal-vision separation (worst 24.0) all pass.
    The aqua slot sits at 2.74:1 against the surface, below the 3:1 bar, so
    every figure that uses it carries visible numeric labels and has a
    machine-readable table in results/tables as its accessible twin.

  ordinal #86b6ef, #5598e7, #2a78d6, #1c5cab, #104281
    monotone in lightness, adjacent lightness gaps above the floor, light end
    at 2.06:1, single hue. Used only for the plaque score, which is an ordered
    scale, never for nominal categories.

These are print figures, so there is no dark variant and no hover layer. Series
identity is carried by marker shape as well as colour, so the figures survive
greyscale reproduction.
"""

from __future__ import annotations

import matplotlib
import matplotlib.pyplot as plt

SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
AXIS_RULE = "#c3c2b7"

# Categorical slots, assigned in fixed order and never cycled.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
          "#e87ba4", "#008300", "#6250d6", "#e34948"]
BLUE, ORANGE, AQUA = SERIES[0], SERIES[1], SERIES[2]

# Ordered ramp for the 0 to 4 plaque score.
ORDINAL_BLUE = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]

# Status steps, reserved for meaning and never used as a series colour.
STATUS_GOOD = "#0ca30c"
STATUS_CRITICAL = "#d03b3b"


def apply_style() -> None:
    matplotlib.rcParams.update({
        "figure.facecolor": SURFACE,
        "figure.edgecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "savefig.edgecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "axes.edgecolor": AXIS_RULE,
        "axes.linewidth": 0.8,
        "axes.labelcolor": INK_SECONDARY,
        "axes.titlecolor": INK_PRIMARY,
        "axes.titlesize": 11,
        "axes.titleweight": "semibold",
        "axes.labelsize": 9.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": GRIDLINE,
        "grid.linewidth": 0.6,
        "grid.linestyle": "-",
        "xtick.color": INK_MUTED,
        "ytick.color": INK_MUTED,
        "xtick.labelcolor": INK_SECONDARY,
        "ytick.labelcolor": INK_SECONDARY,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.frameon": False,
        "legend.fontsize": 8.5,
        "legend.labelcolor": INK_SECONDARY,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Segoe UI", "Arial"],
        "font.size": 9.5,
        "text.color": INK_PRIMARY,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.dpi": 110,
    })


def save(fig, outdir, stem: str, dpi: int = 300) -> list:
    """Write a raster copy for review and a vector copy for typesetting.

    The PDF writer stamps a creation date by default, which would make an
    otherwise identical rerun produce a different file. Suppressing it keeps
    reruns byte-identical, so a change in a figure always means a change in the
    data behind it.
    """
    written = []
    for ext in ("png", "pdf"):
        path = outdir / f"{stem}.{ext}"
        metadata = {"CreationDate": None} if ext == "pdf" else None
        fig.savefig(path, dpi=dpi, bbox_inches="tight",
                    facecolor=SURFACE, edgecolor=SURFACE, metadata=metadata)
        written.append(path)
    plt.close(fig)
    return written


def tidy(ax, xgrid: bool = True, ygrid: bool = False) -> None:
    """Hairline grid on one axis only, so the marks stay dominant."""
    ax.grid(xgrid, axis="x")
    ax.grid(ygrid, axis="y")
    ax.tick_params(length=3, width=0.8)
