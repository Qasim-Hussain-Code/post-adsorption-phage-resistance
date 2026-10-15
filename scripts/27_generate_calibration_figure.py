"""Figure for the association study: why the null had to be rebuilt twice.

The calibration argument behind the accessory genome scan is the most involved
piece of reasoning in this project, and in the tables it is four numbers per
null. Shown as distributions it can be checked by eye: a null that is correct
puts the p-values on the diagonal, and the two that are not miss it in
opposite directions.

Palette and style are those validated for script 06; no new colour slots are
introduced, so no re-validation is required.

Outputs:
    figures/fig9_association_calibration.{png,pdf}
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from papr import config, paths, plotting as P

NULLS = [
    ("p_surface_adjusted_simulated_null", "matched simulation", P.ORANGE),
    ("p_surface_adjusted", "binned empirical", P.AQUA),
    ("p_surface_adjusted_adaptive", "adaptive, normal tail", P.SERIES[3]),
    ("p_surface_adjusted_fitted_tail", "adaptive, fitted tail", P.BLUE),
]


def main() -> int:
    paths.bootstrap()
    P.apply_style()

    loci = pd.read_csv(paths.TABLES / "table48_gwas_all_loci.csv")
    cal = pd.read_csv(paths.TABLES / "table46_gwas_calibration.csv")

    fig = plt.figure(figsize=(14.2, 11.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.3, wspace=0.26)

    # a: the whole argument in one panel. Uniform p-values track the diagonal.
    ax = fig.add_subplot(gs[0, 0])
    for col, label, colour in NULLS:
        v = np.sort(loci[col].dropna().values)
        if not v.size:
            continue
        expected = (np.arange(1, v.size + 1) - 0.5) / v.size
        ax.plot(expected, v, color=colour, linewidth=2.0, label=label,
                zorder=3, solid_capstyle="round")
    ax.plot([0, 1], [0, 1], color=P.INK_MUTED, linewidth=1.0, linestyle="--",
            zorder=2)
    ax.annotate("a null that is correct\nlies on the diagonal", xy=(0.70, 0.70),
                xytext=(0.52, 0.34), fontsize=8.6, color=P.INK_SECONDARY,
                ha="center",
                arrowprops=dict(arrowstyle="-", color=P.INK_MUTED,
                                linewidth=0.8))
    ax.annotate("too liberal: p-values\ncrowd towards zero", xy=(0.42, 0.14),
                xytext=(0.78, 0.09), fontsize=8.6, color=P.INK_SECONDARY,
                ha="center",
                arrowprops=dict(arrowstyle="-", color=P.INK_MUTED,
                                linewidth=0.8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("expected quantile under a correct null")
    ax.set_ylabel("observed p-value")
    ax.set_title("a  Four nulls, one diagonal", loc="left")
    ax.legend(frameon=False, fontsize=8.6, loc="upper left")
    P.tidy(ax, xgrid=True, ygrid=True)

    # b: the tail is what FDR control depends on, so show it on its own scale.
    ax = fig.add_subplot(gs[0, 1])
    rows = cal[cal["outcome"] == "surface_adjusted"].set_index("null")
    order = ["simulated", "empirical", "adaptive", "adaptive_fitted_tail"]
    names = ["matched\nsimulation", "binned\nempirical", "adaptive,\nnormal tail",
             "adaptive,\nfitted tail"]
    colours = [P.ORANGE, P.AQUA, P.SERIES[3], P.BLUE]
    x = np.arange(len(order))
    vals = [100 * float(rows.loc[k, "frac_below_01"]) for k in order]
    ax.bar(x, vals, width=0.58, color=colours, edgecolor=P.SURFACE,
           linewidth=1.1, zorder=3)
    ax.axhline(1.0, color=P.INK_PRIMARY, linewidth=1.1, linestyle="--", zorder=4)
    ax.annotate("1 percent is what a correct null gives",
                xy=(len(order) - 0.5, 1.0), xytext=(len(order) - 0.55, 3.4),
                fontsize=8.6, color=P.INK_SECONDARY, ha="right",
                arrowprops=dict(arrowstyle="-", color=P.INK_MUTED,
                                linewidth=0.8))
    for xi, v in zip(x, vals):
        ax.text(xi, v + 0.16, f"{v:.2f}", ha="center", va="bottom",
                fontsize=8.8, color=P.INK_PRIMARY)
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=8.6)
    ax.set_ylabel("percent of loci with p <= 0.01")
    ax.set_title("b  Behaviour in the tail that matters", loc="left")
    P.tidy(ax, xgrid=False, ygrid=True)

    # c: why the fixed bins failed. Heavy tails sit where they dropped loci.
    ax = fig.add_subplot(gs[1, 0])
    binned = loci["p_surface_adjusted"].notna()
    edges = np.linspace(2, 50, 25)
    # Share of each group, not counts: the groups differ more than twofold in
    # size, and the claim is about the shape of the two distributions.
    for mask, label, colour in (
            (binned, f"kept by the fixed bins (n = {int(binned.sum())})", P.AQUA),
            (~binned, f"dropped by them (n = {int((~binned).sum())})", P.ORANGE)):
        v = loci.loc[mask, "local_tail_df"].dropna().values
        if not v.size:
            continue
        ax.hist(v, bins=edges, weights=np.full(v.size, 100.0 / v.size),
                histtype="step", color=colour, linewidth=2.0, label=label,
                zorder=3)
    ax.set_xlim(2, 52)
    ax.set_ylim(0, 68)
    ax.set_xlabel("fitted degrees of freedom of the local tail "
                  "(50 is the cap, and is effectively normal)")
    ax.set_ylabel("percent of the group")
    ax.set_title("c  Loci the fixed bins dropped have heavier tails", loc="left")
    med_kept = float(loci.loc[binned, "local_tail_df"].median())
    med_drop = float(loci.loc[~binned, "local_tail_df"].median())
    ax.annotate("low df means a heavy tail,\nwhere a normal reference\nover-rejects",
                xy=(3.1, 29.0), xytext=(11.0, 41.0), fontsize=8.6,
                color=P.INK_SECONDARY, ha="left",
                arrowprops=dict(arrowstyle="-", color=P.INK_MUTED,
                                linewidth=0.8))
    ax.text(0.97, 0.47, f"median df {med_kept:.0f} against {med_drop:.1f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8.8,
            color=P.INK_PRIMARY)
    ax.legend(frameon=False, fontsize=8.6, loc="upper right")
    P.tidy(ax, xgrid=True, ygrid=True)

    # d: the three loci the full scan flags, against the diagnostic that
    # disqualifies them. A collapsed, clumped neighbourhood is not a null.
    ax = fig.add_subplot(gs[1, 1])
    flagged = ((loci["q_surface_adjusted_fitted_tail"].fillna(1) <= config.FDR_ALPHA)
               | (loci["q_lysis_fitted_tail"].fillna(1) <= config.FDR_ALPHA))
    ax.scatter(loci.loc[~flagged, "neighbour_clump_fraction"],
               loci.loc[~flagged, "local_scale_ratio"], s=7,
               facecolor=P.INK_MUTED, alpha=0.26, edgecolor="none",
               label="locus not flagged", zorder=3)
    ax.scatter(loci.loc[flagged, "neighbour_clump_fraction"],
               loci.loc[flagged, "local_scale_ratio"], s=86,
               facecolor=P.STATUS_CRITICAL, edgecolor=P.SURFACE, linewidth=1.1,
               label="flagged by the fitted tail", zorder=5)
    ax.set_xlabel("fraction of the neighbourhood piled at its centre")
    ax.set_ylabel("local spread, as a fraction of the spread over all loci")
    ax.set_title("d  Every locus the full scan flags sits in a collapsed null",
                 loc="left")
    ax.annotate("a null estimated from one set of\nnear-identical clade markers "
                "says\nnothing about a different clade",
                xy=(float(loci.loc[flagged, "neighbour_clump_fraction"].mean()),
                    float(loci.loc[flagged, "local_scale_ratio"].mean())),
                xytext=(0.26, 0.60), textcoords="axes fraction", fontsize=8.6,
                color=P.INK_SECONDARY, ha="left",
                arrowprops=dict(arrowstyle="-", color=P.INK_MUTED,
                                linewidth=0.8))
    ax.legend(frameon=False, fontsize=8.6, loc="upper right",
              markerscale=1.0, handletextpad=0.6)
    P.tidy(ax, xgrid=True, ygrid=True)

    written = P.save(fig, paths.FIGURES, "fig9_association_calibration")
    for path in written:
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
