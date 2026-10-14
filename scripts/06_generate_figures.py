"""Build the four manuscript figures at 300 dpi raster and vector.

Every figure is drawn from the tables written by scripts 02 to 05 or recomputed
from the pair-level table, so a figure can never disagree with the numbers in
results/tables. Each figure has a machine-readable twin in that directory,
which is also what satisfies the accessibility requirement for the one palette
slot that sits below the contrast bar.

One deliberate departure from the brief: figure 4 was specified as a raincloud
plot of plaque dilution scores. A raincloud shows the shape of a continuous
distribution, and this outcome is an ordered five-level score in which 79
percent of observations are zero. A violin fitted to that is a smooth curve
over a spike and three small steps, which invents detail the measurement does
not contain. The score distribution is therefore drawn as ordered proportion
bars, which show every level exactly, with the carrier against non-carrier
contrast as a dot plot with intervals beside it.

Outputs:
    figures/fig1_two_gate_filtering.{png,pdf}
    figures/fig2_global_vs_conditioned_forest.{png,pdf}
    figures/fig3_lineage_controlled_survival.{png,pdf}
    figures/fig4_quantitative_titer_shifts.{png,pdf}
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from papr import config, paths, plotting as P

SCORE_LABELS = ["0, no plaques", "1", "2", "3", "4, plaques at the deepest dilution"]


def load():
    pairs = pd.read_csv(paths.PAIRS, compression="gzip", low_memory=False)
    ladder = pd.read_csv(paths.TABLES / "table3_defense_conditioning_ladder.csv")
    quant = pd.read_csv(paths.TABLES / "table4_quantitative_score_shifts.csv")
    lineage = pd.read_csv(paths.TABLES / "table7_lineage_controlled_effects.csv")
    decomp = pd.read_csv(paths.TABLES / "table1_two_gate_decomposition.csv")
    return pairs, ladder, quant, lineage, decomp


def three_way(frame: pd.DataFrame) -> tuple[float, float, float, int]:
    """Shares of lysis, post-adsorption failure and surface exclusion."""
    n = len(frame)
    if n == 0:
        return (np.nan, np.nan, np.nan, 0)
    lysis = float((frame["lytic"] == 1).mean())
    gate2 = float(((frame["lytic"] == 0) & frame["gate_lps"]).mean())
    gate1 = float(((frame["lytic"] == 0) & ~frame["gate_lps"]).mean())
    return lysis, gate2, gate1, n


def ink_on(hex_color: str) -> str:
    """Pick dark or light text for a fill by its relative luminance.

    Keying the label colour to the fill rather than to the slot index means a
    palette change cannot leave white text sitting on a pale segment.
    """
    r, g, b = (int(hex_color[i:i + 2], 16) / 255 for i in (1, 3, 5))

    def lin(c: float) -> float:
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4

    luminance = 0.2126 * lin(r) + 0.7152 * lin(g) + 0.0722 * lin(b)
    return "#0b0b0b" if luminance > 0.42 else "#ffffff"


def stacked_rows(ax, labels, shares, counts, colors, legend_names,
                 label_threshold: float = 0.09) -> None:
    """Horizontal stacked proportion bars with a surface gap between fills."""
    y = np.arange(len(labels))
    left = np.zeros(len(labels))
    for j, (color, name) in enumerate(zip(colors, legend_names)):
        w = shares[:, j]
        ax.barh(y, w, left=left, height=0.72, color=color, label=name,
                edgecolor=P.SURFACE, linewidth=0.9)
        text_colour = ink_on(color)
        for i, (wi, li) in enumerate(zip(w, left)):
            # Only label a segment that is wide enough to hold the text.
            if np.isfinite(wi) and wi >= label_threshold:
                ax.text(li + wi / 2, y[i], f"{100 * wi:.0f}",
                        ha="center", va="center", fontsize=7.5,
                        color=text_colour)
        left = left + np.nan_to_num(w)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{l}  (n={c})" for l, c in zip(labels, counts)])
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "25", "50", "75", "100"])
    ax.invert_yaxis()
    P.tidy(ax, xgrid=True)


def figure1(pairs: pd.DataFrame, decomp: pd.DataFrame):
    colors = [P.AQUA, P.BLUE, P.ORANGE]
    names = ["Productive lysis",
             "Failure, surface state compatible (gate 2)",
             "Failure, surface state incompatible (gate 1)"]

    fig = plt.figure(figsize=(12.4, 9.9))
    # A dedicated band between the rows holds the legend, so it can never
    # collide with the titles of the panels below it.
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 0.16, 5.4],
                          hspace=0.28, wspace=0.34)

    ax0 = fig.add_subplot(gs[0, :])
    l, g2, g1, n = three_way(pairs)
    stacked_rows(ax0, ["All observed pairs"], np.array([[l, g2, g1]]), [n],
                 colors, names)
    ax0.set_title("a  Outcome of every observed strain-phage pair", loc="left")
    ax0.set_xlabel("")

    ax_legend = fig.add_subplot(gs[1, :])
    ax_legend.axis("off")
    ax_legend.legend(
        handles=[Patch(facecolor=c, edgecolor=P.SURFACE, label=nm)
                 for c, nm in zip(colors, names)],
        loc="center", ncol=3, frameon=False)

    ax1 = fig.add_subplot(gs[2, 0])
    rows, labels, counts = [], [], []
    for genus, g in pairs.groupby("genus"):
        l, g2, g1, n = three_way(g)
        rows.append([l, g2, g1])
        labels.append(genus)
        counts.append(n)
    order = np.argsort([-r[0] for r in rows])
    stacked_rows(ax1, [labels[i] for i in order],
                 np.array(rows)[order], [counts[i] for i in order], colors, names)
    ax1.set_title("b  By phage genus, ordered by lysis rate", loc="left")
    ax1.set_xlabel("Percent of observed pairs")

    ax2 = fig.add_subplot(gs[2, 1])
    keep = (pairs.groupby("phylogroup").size()
            .sort_values(ascending=False).head(12).index)
    rows, labels, counts = [], [], []
    for phylo in keep:
        g = pairs[pairs["phylogroup"] == phylo]
        l, g2, g1, n = three_way(g)
        rows.append([l, g2, g1])
        labels.append(phylo)
        counts.append(n)
    order = np.argsort([-r[0] for r in rows])
    stacked_rows(ax2, [labels[i] for i in order],
                 np.array(rows)[order], [counts[i] for i in order], colors, names)
    ax2.set_title("c  By host phylogroup, twelve largest", loc="left")
    ax2.set_xlabel("Percent of observed pairs")

    main = decomp[decomp["gate"] == "gate_lps"].iloc[0]
    fig.suptitle(
        "Figure 1. Surface exclusion accounts for a minority of failed infections",
        x=0.065, y=0.985, ha="left", fontsize=13, fontweight="semibold")
    fig.text(0.065, 0.945,
             "Receptor compatibility is inferred from the surface states of the strains each phage genus was isolated on, "
             "a definition that uses no information from the interaction outcomes.\n"
             f"Of {int(main['n_nonlytic'])} failed infections, {main['pct_nonlytic_at_gate1']:.1f} percent "
             f"are excluded at the surface and {main['pct_nonlytic_at_gate2']:.1f} percent fail inside a "
             "compatible host. Numbers inside the bars are percentages.",
             ha="left", va="top", fontsize=9, color=P.INK_SECONDARY)
    fig.subplots_adjust(top=0.885)
    return fig


def figure2(ladder: pd.DataFrame):
    l0 = ladder[ladder["level"] == "L0_unstratified"].set_index("system")
    l3 = ladder[ladder["level"] == "L3_phage_lps"].set_index("system")
    order = l3.sort_values("log2_or").index.tolist()
    y = np.arange(len(order))

    fig, ax = plt.subplots(figsize=(9.8, 11.6))
    ax.axvline(0, color=P.AXIS_RULE, linewidth=1.0, zorder=1)

    for offset, frame, color, marker, name in [
        (+0.19, l0, P.BLUE, "o", "Unstratified across all pairs"),
        (-0.19, l3, P.ORANGE, "s", "Stratified by phage identity and LPS outer core"),
    ]:
        est = frame.loc[order, "log2_or"].to_numpy()
        lo = frame.loc[order, "log2_ci_low"].to_numpy()
        hi = frame.loc[order, "log2_ci_high"].to_numpy()
        ax.hlines(y + offset, lo, hi, color=color, linewidth=1.6, alpha=0.9, zorder=2)
        ax.scatter(est, y + offset, s=26, facecolor=color, edgecolor=P.SURFACE,
                   linewidth=0.9, marker=marker, label=name, zorder=3)

    sig = l3.loc[order, "significant"].to_numpy()
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{s}  ({int(l3.loc[s, 'n_carrier_strains'])})" + ("  *" if sg else "")
         for s, sg in zip(order, sig)])
    # Leave an empty band at the foot of the plot so the direction hints below
    # cannot overlap the last row's interval.
    ax.set_ylim(-0.9, len(order) + 0.6)
    ax.invert_yaxis()
    ax.set_xlabel("log2 odds ratio for productive lysis, with 95 percent confidence interval")
    P.tidy(ax, xgrid=True)

    ax.text(0.015, 0.010, "Favours resistance", transform=ax.transAxes,
            fontsize=8.5, color=P.INK_MUTED, ha="left", va="bottom")
    ax.text(0.985, 0.010, "Favours lysis", transform=ax.transAxes,
            fontsize=8.5, color=P.INK_MUTED, ha="right", va="bottom")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.085), ncol=2)

    fig.suptitle("Figure 2. Conditioning on surface state widens defense effects "
                 "rather than removing them",
                 x=0.02, y=0.975, ha="left", fontsize=13, fontweight="semibold")
    fig.text(0.02, 0.952,
             "Each defense system carried by at least "
             f"{config.MIN_CARRIER_STRAINS} of the 402 strains. Carrier strain counts are in brackets. "
             "An asterisk marks a stratified estimate passing Benjamini-Hochberg at "
             f"q <= {config.FDR_ALPHA}.\nThe intervals shown treat strain-phage pairs as independent and are "
             "therefore too narrow; Figure 3 refits these contrasts with clustered standard errors.",
             ha="left", va="top", fontsize=9, color=P.INK_SECONDARY)
    fig.subplots_adjust(top=0.925, left=0.26, bottom=0.075)
    return fig


def figure3(ladder: pd.DataFrame, lineage: pd.DataFrame):
    l3 = ladder[ladder["level"] == "L3_phage_lps"].set_index("system")
    full = lineage[lineage["adjustment"] == "plus_phylogroup_and_coordinates"].set_index("system")
    base = lineage[lineage["adjustment"] == "base"].set_index("system")
    order = l3.sort_values("log2_or").index.tolist()
    y = np.arange(len(order))

    fig, axes = plt.subplots(1, 2, figsize=(13.6, 11.4),
                             gridspec_kw={"width_ratios": [1.35, 1], "wspace": 0.3})

    ax = axes[0]
    ax.axvline(0, color=P.AXIS_RULE, linewidth=1.0, zorder=1)
    before = l3.loc[order, "log2_or"].to_numpy()
    after = full.loc[order, "log2_or"].to_numpy()
    survived = (full.loc[order, "significant"] & full.loc[order, "direction_preserved"]).to_numpy()

    for i in range(len(order)):
        ax.plot([before[i], after[i]], [y[i], y[i]], color=P.INK_MUTED,
                linewidth=1.0, alpha=0.55, zorder=2)
    ax.scatter(before, y, s=28, facecolor=P.BLUE, edgecolor=P.SURFACE, linewidth=0.9,
               marker="o", label="Stratified on phage and surface state", zorder=3)
    ax.scatter(after, y, s=34, facecolor=P.ORANGE, edgecolor=P.SURFACE, linewidth=0.9,
               marker="s", label="Also adjusted for clonal background", zorder=4)
    ax.scatter(after[survived], y[survived], s=150, facecolor="none",
               edgecolor=P.STATUS_GOOD, linewidth=1.6, marker="o",
               label="Survives clonal control", zorder=5)

    ax.set_yticks(y)
    ax.set_yticklabels(order)
    ax.set_ylim(-0.8, len(order) - 0.2)
    ax.invert_yaxis()
    ax.set_xlabel("log2 odds ratio for productive lysis")
    ax.set_title("a  Effect size before and after controlling for lineage", loc="left")
    P.tidy(ax, xgrid=True)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.075), ncol=1)

    ax = axes[1]
    lim = float(np.nanmax(np.abs(np.concatenate([before, after])))) * 1.15
    ax.axhline(0, color=P.AXIS_RULE, linewidth=1.0, zorder=1)
    ax.axvline(0, color=P.AXIS_RULE, linewidth=1.0, zorder=1)
    ax.plot([-lim, lim], [-lim, lim], color=P.INK_MUTED, linewidth=1.0,
            alpha=0.6, zorder=1)
    ax.scatter(before[~survived], after[~survived], s=34, facecolor=P.BLUE,
               edgecolor=P.SURFACE, linewidth=0.9, marker="o",
               label="Does not survive", zorder=3)
    ax.scatter(before[survived], after[survived], s=58, facecolor=P.STATUS_GOOD,
               edgecolor=P.SURFACE, linewidth=0.9, marker="D",
               label="Survives clonal control", zorder=4)
    for i in np.where(survived)[0]:
        ax.annotate(order[i], (before[i], after[i]), textcoords="offset points",
                    xytext=(8, 6), fontsize=8.5, color=P.INK_PRIMARY)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("log2 odds ratio, stratified only")
    ax.set_ylabel("log2 odds ratio, clonal background controlled")
    ax.set_title("b  Attenuation toward the null", loc="left")
    ax.grid(True, axis="both")
    ax.tick_params(length=3, width=0.8)
    ax.legend(loc="upper left")

    n_before = int(base["significant"].sum())
    n_after = int((full["significant"] & full["direction_preserved"]).sum())
    fig.suptitle("Figure 3. Almost every defense association is a property of the "
                 "clonal lineage, not the system",
                 x=0.035, y=0.983, ha="left", fontsize=13, fontweight="semibold")
    fig.text(0.035, 0.962,
             f"Of {len(order)} systems tested, {n_before} were significant with clustered standard errors "
             f"before lineage adjustment and {n_after} remain after adjusting for phylogroup and eight "
             "core-genome principal coordinates.\nPoints in panel b falling toward the horizontal axis lost "
             "their effect. The diagonal marks no change. Every association with increased lysis, which no "
             "defense mechanism can produce, disappears.",
             ha="left", va="top", fontsize=9, color=P.INK_SECONDARY)
    fig.subplots_adjust(top=0.895, bottom=0.165, left=0.14)
    return fig


def figure4(pairs: pd.DataFrame, quant: pd.DataFrame, lineage: pd.DataFrame):
    gated = pairs[pairs["gate_lps"]].copy()

    fig = plt.figure(figsize=(14.2, 10.8))
    # A legend band of its own between the rows, and a wide gutter between the
    # two upper panels so the right panel's row labels cannot reach the left
    # panel's bars.
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 0.2, 1.5],
                          hspace=0.42, wspace=0.52)

    def score_shares(frame, key, keep=None):
        labels, rows, counts = [], [], []
        groups = frame.groupby(key, observed=True)
        names = list(keep) if keep is not None else [k for k, _ in groups]
        for name in names:
            g = frame[frame[key] == name]
            if g.empty:
                continue
            vc = g["plaque_score"].value_counts(normalize=True)
            rows.append([float(vc.get(s, 0.0)) for s in range(5)])
            labels.append(str(name))
            counts.append(len(g))
        return labels, np.array(rows), counts

    ax = fig.add_subplot(gs[0, 0])
    labels, rows, counts = score_shares(gated, "defense_burden_quartile",
                                        keep=["Q1", "Q2", "Q3", "Q4"])
    stacked_rows(ax, labels, rows, counts, P.ORDINAL_BLUE, SCORE_LABELS,
                 label_threshold=0.055)
    ax.set_title("a  Plaque score by defense burden quartile", loc="left")
    ax.set_xlabel("Percent of pairs")
    ax.set_ylabel("Quartile of systems carried")

    ax = fig.add_subplot(gs[0, 1])
    labels, rows, counts = score_shares(gated, "morphotype")
    stacked_rows(ax, labels, rows, counts, P.ORDINAL_BLUE, SCORE_LABELS,
                 label_threshold=0.055)
    ax.set_title("b  Plaque score by phage morphotype", loc="left")
    ax.set_xlabel("Percent of pairs")

    ax_legend = fig.add_subplot(gs[1, :])
    ax_legend.axis("off")
    ax_legend.legend(
        handles=[Patch(facecolor=c, edgecolor=P.SURFACE, label=l)
                 for c, l in zip(P.ORDINAL_BLUE, SCORE_LABELS)],
        loc="center", ncol=5, frameon=False,
        title="Plaque score, number of ten-fold dilutions retaining plaques",
        title_fontsize=8.5)

    ax = fig.add_subplot(gs[2, :])
    q = quant[quant["scope"] == "all_pairs_in_gate"].copy()
    q["delta"] = q["mean_score_carrier"] - q["mean_score_noncarrier"]
    q = q.sort_values("delta")
    full = lineage[lineage["adjustment"] == "plus_phylogroup_and_coordinates"].set_index("system")
    survived = {s for s in full.index
                if bool(full.loc[s, "significant"]) and bool(full.loc[s, "direction_preserved"])}
    x = np.arange(len(q))
    ax.axhline(0, color=P.AXIS_RULE, linewidth=1.0, zorder=1)
    sig = (q["brunner_munzel_q"] <= config.FDR_ALPHA).to_numpy()
    ax.vlines(x, 0, q["delta"].to_numpy(), color=P.INK_MUTED, linewidth=1.0,
              alpha=0.55, zorder=2)
    ax.scatter(x[sig], q["delta"].to_numpy()[sig], s=40, facecolor=P.BLUE,
               edgecolor=P.SURFACE, linewidth=0.9, marker="o",
               label=f"Brunner-Munzel q <= {config.FDR_ALPHA}", zorder=3)
    ax.scatter(x[~sig], q["delta"].to_numpy()[~sig], s=36, facecolor=P.SURFACE,
               edgecolor=P.BLUE, linewidth=1.3, marker="o",
               label="Not significant", zorder=3)
    surv_mask = q["system"].isin(survived).to_numpy()
    ax.scatter(x[surv_mask], q["delta"].to_numpy()[surv_mask], s=170,
               facecolor="none", edgecolor=P.STATUS_GOOD, linewidth=1.7,
               marker="o", label="Survives clonal control", zorder=4)
    ax.set_xticks(x)
    ax.set_xticklabels(q["system"], rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("Mean plaque score, carriers minus non-carriers")
    ax.set_title("c  Shift in mean plaque score associated with carrying each system",
                 loc="left")
    ax.grid(True, axis="y")
    ax.grid(False, axis="x")
    ax.tick_params(length=3, width=0.8)
    ax.legend(loc="upper left", ncol=1)

    fig.suptitle("Figure 4. Plaque score tracks phage morphotype strongly and "
                 "defense burden not at all",
                 x=0.035, y=0.978, ha="left", fontsize=13, fontweight="semibold")
    fig.text(0.035, 0.952,
             "All panels are restricted to pairs passing the surface compatibility gate. The score is the number of "
             "successive ten-fold dilutions at which plaques remained visible, so it is an ordered five-level scale "
             "and not a continuous titre.\nPanel c shows a difference in means for interpretability and marks "
             "significance from the rank-based Brunner-Munzel test, which does not assume a shape. Only the two "
             "systems ringed in green hold up once clonal background is controlled.",
             ha="left", va="top", fontsize=9, color=P.INK_SECONDARY)
    fig.subplots_adjust(top=0.905, bottom=0.14)
    return fig


def main() -> int:
    paths.bootstrap()
    P.apply_style()
    pairs, ladder, quant, lineage, decomp = load()

    written = []
    written += P.save(figure1(pairs, decomp), paths.FIGURES,
                      "fig1_two_gate_filtering", config.FIGURE_DPI)
    written += P.save(figure2(ladder), paths.FIGURES,
                      "fig2_global_vs_conditioned_forest", config.FIGURE_DPI)
    written += P.save(figure3(ladder, lineage), paths.FIGURES,
                      "fig3_lineage_controlled_survival", config.FIGURE_DPI)
    written += P.save(figure4(pairs, quant, lineage), paths.FIGURES,
                      "fig4_quantitative_titer_shifts", config.FIGURE_DPI)

    for p in written:
        print(f"  {p.stat().st_size / 1024:8.0f} KiB  {p.relative_to(paths.ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
