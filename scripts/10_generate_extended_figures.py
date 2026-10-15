"""Figures for the direct phage annotation and the robustness analyses.

Kept separate from script 06 so that the core pipeline, steps 01 to 06, stays
runnable without the heavier annotation step that steps 07 to 09 depend on.

Palette and style are those validated for script 06; no new colour slots are
introduced, so no re-validation is required.

Outputs:
    figures/fig5_direct_antidefense_annotation.{png,pdf}
    figures/fig6_robustness_of_surviving_effects.{png,pdf}
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

from papr import config, paths, plotting as P


def figure5(adf: pd.DataFrame, validation: pd.DataFrame,
            comparison: pd.DataFrame):
    families = [c for c in adf.columns
                if c.startswith("Anti_") or c in ("NADP", "Other")]
    counts = adf[families].sum().sort_values()
    counts = counts[counts > 0]

    fig = plt.figure(figsize=(14.0, 9.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.42, wspace=0.3)

    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(len(counts))
    ax.barh(y, counts.to_numpy(), height=0.7, color=P.BLUE,
            edgecolor=P.SURFACE, linewidth=0.9)
    for i, v in enumerate(counts.to_numpy()):
        ax.text(v + 0.6, i, f"{v}  ({100 * v / len(adf):.0f}%)",
                va="center", fontsize=8, color=P.INK_SECONDARY)
    ax.set_yticks(y)
    ax.set_yticklabels(counts.index)
    ax.set_xlim(0, counts.max() * 1.32)
    ax.set_xlabel(f"Phages carrying the family, of {len(adf)}")
    ax.set_title("a  Anti-defense families found in the 96 phage genomes", loc="left")
    P.tidy(ax, xgrid=True)

    ax = fig.add_subplot(gs[0, 1])
    v = validation.dropna(subset=["prevalence_refseq_imputed"]).copy()
    ax.plot([0, 1], [0, 1], color=P.INK_MUTED, linewidth=1.0, alpha=0.7, zorder=1)
    sizes = 18 + 9 * v["n_phages_here"].to_numpy()
    ax.scatter(v["prevalence_refseq_imputed"], v["prevalence_direct"],
               s=sizes, facecolor=P.BLUE, edgecolor=P.SURFACE, linewidth=0.9,
               alpha=0.9, zorder=3)
    for r in v.itertuples(index=False):
        if abs(r.prevalence_direct - r.prevalence_refseq_imputed) > 0.12:
            ax.annotate(r.genus, (r.prevalence_refseq_imputed, r.prevalence_direct),
                        textcoords="offset points", xytext=(8, -3),
                        fontsize=8, color=P.INK_PRIMARY)
    rho = v[["prevalence_direct", "prevalence_refseq_imputed"]].corr(
        method="spearman").iloc[0, 1]
    ax.set_xlim(-0.05, 1.08)
    ax.set_ylim(-0.05, 1.08)
    ax.set_xlabel("Prevalence imputed from RefSeq relatives of the genus")
    ax.set_ylabel("Prevalence measured in these phages")
    ax.set_title("b  Direct annotation against the published survey", loc="left")
    # Bottom right is the empty quadrant here: genera cluster on the diagonal
    # and near the top right, so a caption there cannot collide with a label.
    ax.text(0.97, 0.03, f"Spearman {rho:+.3f} over {len(v)} genera\n"
                        "point size is phages per genus",
            transform=ax.transAxes, fontsize=8.5, ha="right", va="bottom",
            color=P.INK_SECONDARY)
    ax.grid(True, axis="both")
    ax.tick_params(length=3, width=0.8)

    ax = fig.add_subplot(gs[1, :])
    comp = comparison.dropna(
        subset=["direct_interaction_or", "imputed_interaction_or"]).copy()
    comp = comp.sort_values("direct_interaction_or")
    y = np.arange(len(comp))
    ax.axvline(0, color=P.AXIS_RULE, linewidth=1.0, zorder=1)
    d = np.log2(comp["direct_interaction_or"].to_numpy())
    lo = np.log2(comp["direct_ci_low"].to_numpy())
    hi = np.log2(comp["direct_ci_high"].to_numpy())
    imp = np.log2(comp["imputed_interaction_or"].to_numpy())
    ax.hlines(y - 0.16, lo, hi, color=P.BLUE, linewidth=1.6, alpha=0.9, zorder=2)
    ax.scatter(d, y - 0.16, s=30, facecolor=P.BLUE, edgecolor=P.SURFACE,
               linewidth=0.9, marker="o", label="Direct annotation of each genome",
               zorder=3)
    ax.scatter(imp, y + 0.16, s=30, facecolor=P.ORANGE, edgecolor=P.SURFACE,
               linewidth=0.9, marker="s", label="Imputed from the phage genus",
               zorder=3)
    sig = comp["direct_significant"].to_numpy()
    ax.set_yticks(y)
    ax.set_yticklabels([f"{s}  ({f.replace('Anti_', '')})" + ("  *" if g else "")
                        for s, f, g in zip(comp["system"],
                                           comp["antidefense_family"], sig)])
    ax.set_ylim(-0.8, len(comp) + 0.4)
    ax.invert_yaxis()
    ax.set_xlabel("log2 interaction odds ratio: host system against phages carrying "
                  "the matching anti-defense family")
    ax.set_title("c  Rescue interactions, measured annotation against the genus prior",
                 loc="left")
    P.tidy(ax, xgrid=True)
    ax.text(0.012, 0.012, "System still protects", transform=ax.transAxes,
            fontsize=8.5, color=P.INK_MUTED, ha="left", va="bottom")
    ax.text(0.988, 0.012, "System protects less, as rescue predicts",
            transform=ax.transAxes, fontsize=8.5, color=P.INK_MUTED,
            ha="right", va="bottom")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2)

    fig.suptitle("Figure 5. Annotating the phages directly sharpens the test and "
                 "leaves its answer unchanged",
                 x=0.035, y=0.982, ha="left", fontsize=13, fontweight="semibold")
    fig.text(0.035, 0.958,
             "Genes were called with Prodigal and searched against the 178 curated anti-defense profiles at their own "
             "gathering cutoffs, the same criterion AntiDefenseFinder applies.\n"
             "Panel c compares the interaction estimated from these measured calls with the one estimated from "
             "genus-level prevalence; an asterisk marks significance at q <= "
             f"{config.FDR_ALPHA} under the direct annotation.",
             ha="left", va="top", fontsize=9, color=P.INK_SECONDARY)
    fig.subplots_adjust(top=0.895, bottom=0.135, left=0.17)
    return fig


def figure6(mde: pd.DataFrame, sweep: pd.DataFrame, within: pd.DataFrame,
            het: pd.DataFrame, survivors: list[str]):
    fig = plt.figure(figsize=(14.0, 13.2))
    # Panel a carries one row per system, so it needs roughly twice the height
    # of the lower row for the labels to clear one another.
    gs = fig.add_gridspec(2, 2, height_ratios=[2.0, 1], hspace=0.3, wspace=0.3)

    ax = fig.add_subplot(gs[0, :])
    m = mde.copy().sort_values("odds_ratio")
    y = np.arange(len(m))
    ax.axvline(0, color=P.AXIS_RULE, linewidth=1.0, zorder=1)
    est = np.log2(m["odds_ratio"].to_numpy())
    lo = np.log2(m["ci_low"].to_numpy())
    hi = np.log2(m["ci_high"].to_numpy())
    sig = m["significant"].to_numpy()
    ax.hlines(y, lo, hi, color=P.BLUE, linewidth=1.5, alpha=0.85, zorder=2)
    ax.scatter(est[~sig], y[~sig], s=28, facecolor=P.SURFACE, edgecolor=P.BLUE,
               linewidth=1.3, marker="o", label="Not significant", zorder=3)
    ax.scatter(est[sig], y[sig], s=46, facecolor=P.STATUS_GOOD,
               edgecolor=P.SURFACE, linewidth=0.9, marker="D",
               label="Survives clonal control", zorder=4)
    det = np.log2(m["min_detectable_or_protective"].to_numpy())
    ax.scatter(det, y, s=34, facecolor=P.ORANGE, edgecolor=P.SURFACE,
               linewidth=0.9, marker="|", zorder=3,
               label="Protective effect detectable at 80 percent power")
    ax.set_yticks(y)
    ax.set_yticklabels(m["system"], fontsize=8)
    ax.set_ylim(-0.8, len(m) + 0.4)
    ax.invert_yaxis()
    ax.set_xlabel("log2 odds ratio for productive lysis")
    ax.set_title("a  What each test found, and what it was able to find", loc="left")
    P.tidy(ax, xgrid=True)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.055), ncol=3)

    ax = fig.add_subplot(gs[1, 0])
    ax.plot(sweep["n_clusters"], sweep["n_protective"], marker="o",
            color=P.BLUE, linewidth=1.8, markersize=6,
            markeredgecolor=P.SURFACE, markeredgewidth=0.9,
            label="Associated with reduced lysis")
    ax.plot(sweep["n_clusters"], sweep["n_increased_lysis"], marker="s",
            color=P.ORANGE, linewidth=1.8, markersize=6,
            markeredgecolor=P.SURFACE, markeredgewidth=0.9,
            label="Associated with increased lysis")
    top = max(3, int(sweep[["n_protective", "n_increased_lysis"]].to_numpy().max()))
    ax.set_ylim(-0.35, top + 1.4)
    used = sweep[np.isclose(sweep["threshold"], config.CLONAL_DISTANCE_THRESHOLD)]
    if len(used):
        ax.axvline(float(used["n_clusters"].iloc[0]), color=P.INK_MUTED,
                   linewidth=1.0, alpha=0.8)
        ax.annotate("threshold used", (float(used["n_clusters"].iloc[0]), -0.2),
                    textcoords="offset points", xytext=(6, 0),
                    fontsize=8, color=P.INK_MUTED, va="bottom")
    ax.set_xlabel("Clonal clusters at the threshold")
    ax.set_ylabel("Systems significant at q <= 0.05")
    ax.set_title("b  Sensitivity to the clonal distance threshold", loc="left")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.grid(True, axis="both")
    ax.tick_params(length=3, width=0.8)
    ax.legend(loc="upper right")

    ax = fig.add_subplot(gs[1, 1])
    if not het.empty:
        offsets = {s: (i - (len(survivors) - 1) / 2) * 0.22
                   for i, s in enumerate(survivors)}
        colours = {s: c for s, c in zip(survivors, [P.BLUE, P.ORANGE, P.AQUA])}
        markers = {s: m for s, m in zip(survivors, ["o", "s", "D"])}
        genera = (het.groupby("genus")["odds_ratio"].median()
                  .sort_values().index.tolist())
        ypos = {g: i for i, g in enumerate(genera)}
        ax.axvline(0, color=P.AXIS_RULE, linewidth=1.0, zorder=1)
        for s in survivors:
            sub = het[het["system"] == s]
            if sub.empty:
                continue
            yy = np.array([ypos[g] for g in sub["genus"]]) + offsets[s]
            ax.hlines(yy, np.log2(sub["ci_low"]), np.log2(sub["ci_high"]),
                      color=colours[s], linewidth=1.4, alpha=0.85, zorder=2)
            ax.scatter(np.log2(sub["odds_ratio"]), yy, s=26,
                       facecolor=colours[s], edgecolor=P.SURFACE, linewidth=0.9,
                       marker=markers[s], label=s, zorder=3)
        ax.set_yticks(range(len(genera)))
        ax.set_yticklabels(genera, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("log2 odds ratio within the genus")
        ax.legend(loc="lower left", fontsize=8)
    ax.set_title("c  Do the survivors protect across phage genera", loc="left")
    P.tidy(ax, xgrid=True)

    n_sig = int(mde["significant"].sum())
    fig.suptitle("Figure 6. Neither surviving association is robust to every "
                 "design, and most of the nulls are uninformative",
                 x=0.035, y=0.985, ha="left", fontsize=13, fontweight="semibold")
    fig.text(0.035, 0.966,
             f"Panel a places each system's estimate beside the protective odds ratio its own standard error would "
             "have detected at 80 percent power. Where that marker sits close to the estimate, a null excludes "
             "little, which is the case for most systems here.\n"
             f"{n_sig} of {len(mde)} systems survive clonal control and that count is stable across every clonal "
             "distance threshold tried, but neither survivor passes all of the stricter tests in "
             "results/robustness_report.txt. Panel c shows both are protective in direction across most phage "
             "genera.",
             ha="left", va="top", fontsize=9, color=P.INK_SECONDARY)
    fig.subplots_adjust(top=0.925, bottom=0.07, left=0.16)
    return fig


def main() -> int:
    paths.bootstrap()
    P.apply_style()

    adf = pd.read_csv(paths.PHAGE_ADF_SYSTEMS)
    validation = pd.read_csv(paths.TABLES / "table15_antidefense_direct_vs_imputed.csv")
    comparison = pd.read_csv(paths.TABLES / "table23_direct_vs_imputed_rescue.csv")

    written = []
    written += P.save(figure5(adf, validation, comparison), paths.FIGURES,
                      "fig5_direct_antidefense_annotation", config.FIGURE_DPI)

    robust = paths.TABLES / "table18_detectable_effects.csv"
    if robust.exists():
        mde = pd.read_csv(robust)
        sweep = pd.read_csv(paths.TABLES / "table19_threshold_sensitivity.csv")
        within = pd.read_csv(paths.TABLES / "table16_within_cluster_effects.csv")
        het = pd.read_csv(paths.TABLES / "table21_genus_heterogeneity.csv")
        survivors = sorted(mde[mde["significant"]]["system"])
        written += P.save(figure6(mde, sweep, within, het, survivors),
                          paths.FIGURES, "fig6_robustness_of_surviving_effects",
                          config.FIGURE_DPI)
    else:
        print("  robustness tables not present; skipping figure 6")

    for p in written:
        print(f"  {p.stat().st_size / 1024:8.0f} KiB  {p.relative_to(paths.ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
