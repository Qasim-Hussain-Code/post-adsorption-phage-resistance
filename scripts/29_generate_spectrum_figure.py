"""Figure for the phage-resolved analysis: what averaging over a panel costs.

Palette and style are those validated for script 06; no new colour slots are
introduced, so no re-validation is required.

Outputs:
    figures/fig10_phage_resolved_spectrum.{png,pdf}
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from papr import config, paths, plotting as P

Z_THRESHOLD = 2.5758293035489004


def main() -> int:
    paths.bootstrap()
    P.apply_style()

    spec = pd.read_csv(paths.TABLES / "table50_spectrum_tests.csv")
    per_phage = pd.read_csv(paths.TABLES / "table49_phage_resolved_effects.csv")
    structure = pd.read_csv(paths.TABLES / "table54_trait_structure.csv")
    survivors = set(spec.loc[spec["q_protective"] <= config.FDR_ALPHA, "system"])

    fig = plt.figure(figsize=(14.2, 17.0))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1], hspace=0.34,
                          wspace=0.26)

    # a: the panel average against the effect on the phages a system acts on.
    ax = fig.add_subplot(gs[0, 0])
    d = spec.dropna(subset=["spectrum_difference"])
    other = d[~d["system"].isin(survivors)]
    won = d[d["system"].isin(survivors)]
    ax.axhline(0, color=P.AXIS_RULE, linewidth=0.9, zorder=2)
    ax.axvline(0, color=P.AXIS_RULE, linewidth=0.9, zorder=2)
    ax.plot([-0.35, 0.05], [-0.35, 0.05], color=P.INK_MUTED, linewidth=1.0,
            linestyle="--", zorder=2)
    ax.scatter(other["panel_difference"], other["spectrum_difference"], s=38,
               facecolor=P.INK_MUTED, alpha=0.55, edgecolor=P.SURFACE,
               linewidth=0.8, zorder=3, label="does not clear both nulls")
    ax.scatter(won["panel_difference"], won["spectrum_difference"], s=120,
               facecolor=P.STATUS_CRITICAL, edgecolor=P.SURFACE, linewidth=1.2,
               zorder=5, label="clears both nulls")
    for r in won.itertuples(index=False):
        ax.annotate(r.system, xy=(r.panel_difference, r.spectrum_difference),
                    xytext=(14, 6), textcoords="offset points", fontsize=9,
                    color=P.INK_PRIMARY, fontweight="semibold")
    ax.set_xlabel("difference in lysis rate, averaged over the whole panel")
    ax.set_ylabel("difference on the phages the system acts on")
    ax.set_title("a  What averaging over 94 phages costs", loc="left")
    ax.annotate("on the diagonal a system\nacts on the panel evenly",
                xy=(-0.17, -0.17), xytext=(0.33, 0.93),
                textcoords="axes fraction", fontsize=8.6,
                color=P.INK_SECONDARY, ha="left",
                arrowprops=dict(arrowstyle="-", color=P.INK_MUTED,
                                linewidth=0.8))
    ax.legend(frameon=False, fontsize=8.6, loc="lower left")
    P.tidy(ax, xgrid=True, ygrid=True)

    # b: the surviving system across every phage in the panel.
    ax = fig.add_subplot(gs[0, 1])
    name = sorted(survivors)[0] if survivors else spec.iloc[0]["system"]
    sub = per_phage[per_phage["system"] == name].sort_values("score_z")
    vals = sub["score_z"].to_numpy()
    colours = np.where(vals < -Z_THRESHOLD, P.STATUS_CRITICAL, P.INK_MUTED)
    ax.bar(np.arange(vals.size), vals, width=0.9, color=colours,
           edgecolor="none", zorder=3)
    ax.axhline(-Z_THRESHOLD, color=P.INK_PRIMARY, linewidth=1.1,
               linestyle="--", zorder=4)
    ax.axhline(0, color=P.AXIS_RULE, linewidth=0.9, zorder=2)
    n_hit = int((vals < -Z_THRESHOLD).sum())
    ax.annotate(f"{n_hit} phages beyond the threshold,\nspanning "
                f"{int(spec.loc[spec.system == name, 'n_groups_protective'].iloc[0])}"
                " distinct phage groups",
                xy=(n_hit / 2, vals.min() * 0.75), xytext=(vals.size * 0.33,
                                                           vals.min() * 0.55),
                fontsize=8.6, color=P.INK_SECONDARY, ha="left",
                arrowprops=dict(arrowstyle="-", color=P.INK_MUTED,
                                linewidth=0.8))
    ax.set_xlabel(f"the {vals.size} phages, ordered by effect")
    ax.set_ylabel("score against lysis (negative is protective)")
    ax.set_title(f"b  {name} across the panel", loc="left")
    # The stability check belongs beside the result, not in a later section.
    drop = paths.TABLES / "table56_matched_drop.csv"
    if drop.exists():
        md = pd.read_csv(drop)
        md = md[md["system"] == name]
        if not md.empty:
            r0 = md.iloc[0]
            n_carr = int(spec.loc[spec.system == name, "n_carriers"].iloc[0])
            ax.text(0.04, 0.97,
                    f"but one clonal cluster holds "
                    f"{int(r0.carriers_dropped)} of the {n_carr} carriers,\n"
                    f"and removing it takes this to "
                    f"p = {r0.p_dropping_cluster:.3f},\nwhere random drops of "
                    f"the same size leave a median\nof "
                    f"p = {r0.median_p_random_drop:.4f}. Those "
                    f"{int(r0.carriers_dropped)} strains carry it.",
                    transform=ax.transAxes, ha="left", va="top",
                    fontsize=8.4, color=P.STATUS_CRITICAL)
    P.tidy(ax, xgrid=False, ygrid=True)

    # c: why one null was not enough.
    ax = fig.add_subplot(gs[1, 0])
    lo = min(structure["coordinate_null_concordance"].min(),
             structure["observed_concordance"].min()) - 0.02
    hi = max(structure["observed_concordance"].max(), 1.0) + 0.01
    ax.plot([lo, hi], [lo, hi], color=P.INK_MUTED, linewidth=1.0,
            linestyle="--", zorder=2)
    ax.scatter(structure["coordinate_null_concordance"],
               structure["observed_concordance"], s=44, facecolor=P.ORANGE,
               edgecolor=P.SURFACE, linewidth=0.9, zorder=4,
               label="against the coordinate null")
    ax.scatter(structure["clonal_null_concordance"],
               structure["observed_concordance"], s=44, facecolor=P.BLUE,
               edgecolor=P.SURFACE, linewidth=0.9, zorder=5,
               label="against the clonal null")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("clonal concordance of the null traits")
    ax.set_ylabel("clonal concordance of the real system")
    ax.set_title("c  The coordinate null is not clonal enough", loc="left")
    ax.annotate("every system sits above the line:\nreal systems track lineages "
                "more\ntightly than that null can",
                xy=(structure["coordinate_null_concordance"].median(),
                    structure["observed_concordance"].median()),
                xytext=(0.04, 0.72), textcoords="axes fraction", fontsize=8.6,
                color=P.INK_SECONDARY, ha="left",
                arrowprops=dict(arrowstyle="-", color=P.INK_MUTED,
                                linewidth=0.8))
    ax.legend(frameon=False, fontsize=8.6, loc="lower right")
    P.tidy(ax, xgrid=True, ygrid=True)

    # d: what each null admits, and what survives both.
    ax = fig.add_subplot(gs[1, 1])
    floor = spec[["p_protective_coordinate", "p_protective_clonal"]].min().min()
    jitter = floor * 0.6
    x = spec["p_protective_coordinate"].clip(lower=floor) + jitter
    y = spec["p_protective_clonal"].clip(lower=floor) + jitter
    keep = spec["system"].isin(survivors)
    ax.axvspan(floor * 0.5, config.FDR_ALPHA, color=P.ORANGE, alpha=0.08,
               zorder=1)
    ax.axhspan(floor * 0.5, config.FDR_ALPHA, color=P.BLUE, alpha=0.08,
               zorder=1)
    ax.scatter(x[~keep], y[~keep], s=40, facecolor=P.INK_MUTED, alpha=0.6,
               edgecolor=P.SURFACE, linewidth=0.8, zorder=3)
    ax.scatter(x[keep], y[keep], s=120, facecolor=P.STATUS_CRITICAL,
               edgecolor=P.SURFACE, linewidth=1.2, zorder=5)
    for r, xi, yi in zip(spec.itertuples(index=False), x, y):
        if r.system in survivors or r.p_protective_coordinate <= config.FDR_ALPHA:
            ax.annotate(r.system, xy=(xi, yi), xytext=(7, 4),
                        textcoords="offset points", fontsize=8.2,
                        color=P.INK_PRIMARY)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("p against the coordinate null")
    ax.set_ylabel("p against the clonal null")
    n_coord = int((spec["p_protective_coordinate"] <= config.FDR_ALPHA).sum())
    n_both = int(((spec["p_protective_coordinate"] <= config.FDR_ALPHA)
                  & (spec["p_protective_clonal"] <= config.FDR_ALPHA)).sum())
    ax.set_title(f"d  {n_coord} clear one null, {n_both} clear both, "
                 f"{len(survivors)} survives the correction", loc="left")
    ax.annotate(f"shading marks p <= {config.FDR_ALPHA} against each null.\n"
                f"Only the filled point also clears a\n"
                f"Benjamini-Hochberg correction for the\n"
                f"{len(spec)} systems tested.",
                xy=(0.03, 0.03), xycoords="axes fraction", fontsize=8.4,
                color=P.INK_SECONDARY, ha="left", va="bottom")
    ax.legend(handles=[
        plt.Line2D([], [], marker="o", linestyle="none", markersize=9,
                   markerfacecolor=P.STATUS_CRITICAL, markeredgecolor=P.SURFACE,
                   label="clears both nulls and the correction"),
        plt.Line2D([], [], marker="o", linestyle="none", markersize=7,
                   markerfacecolor=P.INK_MUTED, markeredgecolor=P.SURFACE,
                   label="does not")],
        frameon=False, fontsize=8.6, loc="upper right")
    P.tidy(ax, xgrid=True, ygrid=True)

    # e: the direction no defense mechanism explains is better supported.
    ax = fig.add_subplot(gs[2, 0])
    floor_q = min(spec[["q_protective", "q_increased"]].min().min(), 0.02)
    qx = spec["q_protective"].clip(lower=floor_q * 0.8)
    qy = spec["q_increased"].clip(lower=floor_q * 0.8)
    prot = spec["q_protective"] <= config.FDR_ALPHA
    incr = spec["q_increased"] <= config.FDR_ALPHA
    ax.axvline(config.FDR_ALPHA, color=P.INK_PRIMARY, linewidth=1.0,
               linestyle="--", zorder=4)
    ax.axhline(config.FDR_ALPHA, color=P.INK_PRIMARY, linewidth=1.0,
               linestyle="--", zorder=4)
    ax.scatter(qx[~prot & ~incr], qy[~prot & ~incr], s=38,
               facecolor=P.INK_MUTED, alpha=0.55, edgecolor=P.SURFACE,
               linewidth=0.8, zorder=3)
    ax.scatter(qx[prot], qy[prot], s=120, facecolor=P.BLUE,
               edgecolor=P.SURFACE, linewidth=1.2, zorder=5,
               label="protective, clears both nulls")
    ax.scatter(qx[incr], qy[incr], s=120, facecolor=P.ORANGE,
               edgecolor=P.SURFACE, linewidth=1.2, zorder=5,
               label="increased lysis, clears both nulls")
    # The flagged systems cluster tightly, so labels are fanned out rather
    # than stacked on one another.
    flagged = [(r, a, b) for r, a, b in zip(spec.itertuples(index=False), qx, qy)
               if r.q_protective <= config.FDR_ALPHA
               or r.q_increased <= config.FDR_ALPHA]
    flagged.sort(key=lambda t: t[1])
    rank = 0
    for r, a, b in flagged:
        if r.q_protective <= config.FDR_ALPHA:
            ax.annotate(r.system, xy=(a, b), xytext=(9, 5),
                        textcoords="offset points", fontsize=8.2, ha="left",
                        color=P.INK_PRIMARY)
            continue
        # The increased-direction systems sit on top of one another at the
        # bottom of the panel, so their labels fan up and to the left.
        ax.annotate(r.system, xy=(a, b), xytext=(-26, 16 + 17 * rank),
                    textcoords="offset points", fontsize=8.2, ha="right",
                    color=P.INK_PRIMARY,
                    arrowprops=dict(arrowstyle="-", color=P.INK_MUTED,
                                    linewidth=0.7))
        rank += 1
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(top=3.0)
    ax.set_xlabel("q for protecting against lysis")
    ax.set_ylabel("q for increasing lysis")
    ax.set_title(f"e  {int(incr.sum())} systems in the impossible direction, "
                 f"{int(prot.sum())} in the sensible one", loc="left")
    ax.legend(frameon=False, fontsize=8.6, loc="lower left")
    P.tidy(ax, xgrid=True, ygrid=True)

    # f: the check that separates them.
    ax = fig.add_subplot(gs[2, 1])
    stab_path = paths.TABLES / "table55_leave_cluster_out.csv"
    if stab_path.exists():
        st = pd.read_csv(stab_path)
        order = (st.groupby(["system", "direction"])["p_value"].max()
                 .sort_values().reset_index())
        order = order.rename(columns={"p_value": "p_value"})
        for i, row in enumerate(order.itertuples(index=False)):
            part = st[(st.system == row.system)
                      & (st.direction == row.direction)]
            colour = P.BLUE if row.direction == "protective" else P.ORANGE
            jitter = (np.arange(len(part)) % 7 - 3) * 0.022
            ax.scatter(part["p_value"], np.full(len(part), i) + jitter, s=26,
                       facecolor=colour, alpha=0.5, edgecolor="none", zorder=3)
            ax.scatter([part["p_value"].max()], [i], s=130, marker="D",
                       facecolor=colour, edgecolor=P.SURFACE, linewidth=1.1,
                       zorder=5)
        ax.axvline(config.FDR_ALPHA, color=P.INK_PRIMARY, linewidth=1.1,
                   linestyle="--", zorder=4)
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels([f"{r.system}\n({r.direction})"
                            for r in order.itertuples(index=False)],
                           fontsize=8.4)
        ax.set_ylim(-0.6, len(order) - 0.4)
        ax.set_xscale("log")
        ax.set_xlabel("p after dropping one clonal cluster of carriers")
        broke = int((order["p_value"] > config.FDR_ALPHA).sum())
        ax.set_title(f"f  {broke} of {len(order)} rest on one group of "
                     "relatives", loc="left")
        ax.annotate("diamonds mark the worst case.\nEverything right of the "
                    "line\nfails when one group of\nrelatives is removed.",
                    xy=(0.97, 0.06), xycoords="axes fraction", fontsize=8.4,
                    color=P.INK_SECONDARY, ha="right", va="bottom")
        P.tidy(ax, xgrid=True, ygrid=False)

    written = P.save(fig, paths.FIGURES, "fig10_phage_resolved_spectrum")
    for path in written:
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
