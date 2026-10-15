"""The summary figure for the extended analysis.

Three things the earlier figures do not show: the corrected partition of
explained variation, what shrinkage does to 79 system estimates at once, and
the consolidated candidate list with the evidence behind each entry.

Palette and style are those validated for script 06; no new colour slots are
introduced, so no re-validation is required.

Outputs:
    figures/fig7_consolidated_evidence.{png,pdf}
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from papr import config, paths, plotting as P


def main() -> int:
    paths.bootstrap()
    P.apply_style()

    shap = pd.read_csv(paths.TABLES / "table25_shapley_variance_partition.csv")
    allsys = pd.read_csv(paths.TABLES / "table29_all_system_effects.csv")
    summary = pd.read_csv(paths.TABLES / "table30_shrinkage_summary.csv").iloc[0]
    cand = pd.read_csv(paths.TABLES / "table33_consolidated_candidates.csv")

    fig = plt.figure(figsize=(14.2, 12.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.3], hspace=0.3, wspace=0.28)

    # Panel a: Shapley partition of explained variation.
    ax = fig.add_subplot(gs[0, 0])
    s = shap.sort_values("shapley_pseudo_r2")
    y = np.arange(len(s))
    ax.barh(y, s["shapley_pseudo_r2"], height=0.65, color=P.BLUE,
            edgecolor=P.SURFACE, linewidth=0.9)
    for i, r in enumerate(s.itertuples(index=False)):
        ax.text(r.shapley_pseudo_r2 + 0.004, i,
                f"{r.shapley_pseudo_r2:.3f}  ({100 * r.share_of_explained:.0f}%, "
                f"{r.n_terms} terms)",
                va="center", fontsize=8, color=P.INK_SECONDARY)
    ax.set_yticks(y)
    ax.set_yticklabels(s["block"])
    ax.set_xlim(0, s["shapley_pseudo_r2"].max() * 1.75)
    ax.set_xlabel("Shapley share of McFadden pseudo R2")
    ax.set_title("a  What explains whether a phage lyses a strain", loc="left")
    P.tidy(ax, xgrid=True)

    # Panel b: shrinkage of 79 estimates toward the pooled mean.
    ax = fig.add_subplot(gs[0, 1])
    d = allsys.dropna(subset=["posterior_odds_ratio", "raw_odds_ratio"]).copy()
    raw = np.log2(d["raw_odds_ratio"].to_numpy())
    shr = np.log2(d["posterior_odds_ratio"].to_numpy())
    excl = d["posterior_excludes_null"].fillna(False).to_numpy().astype(bool)
    lim = float(np.nanmax(np.abs(raw))) * 1.1
    ax.plot([-lim, lim], [-lim, lim], color=P.INK_MUTED, linewidth=1.0,
            alpha=0.6, zorder=1)
    ax.axhline(np.log2(summary["pooled_odds_ratio"]), color=P.ORANGE,
               linewidth=1.2, zorder=2)
    ax.axvline(0, color=P.AXIS_RULE, linewidth=1.0, zorder=1)
    ax.axhline(0, color=P.AXIS_RULE, linewidth=1.0, zorder=1)
    sizes = 14 + 26 * (d["n_carrier_strains"] / d["n_carrier_strains"].max())
    ax.scatter(raw[~excl], shr[~excl], s=sizes[~excl], facecolor=P.SURFACE,
               edgecolor=P.BLUE, linewidth=1.1, marker="o",
               label="Shrunk interval includes 1", zorder=3)
    ax.scatter(raw[excl], shr[excl], s=sizes[excl] + 20, facecolor=P.BLUE,
               edgecolor=P.SURFACE, linewidth=0.9, marker="D",
               label="Shrunk interval excludes 1", zorder=4)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("log2 odds ratio, system fitted on its own")
    ax.set_ylabel("log2 odds ratio after shrinkage")
    ax.set_title(f"b  {len(d)} systems estimated jointly", loc="left")
    ax.text(0.03, 0.97,
            f"pooled OR {summary['pooled_odds_ratio']:.3f}\n"
            f"between-system SD {summary['tau']:.2f}\n"
            "point size is carrier strains",
            transform=ax.transAxes, fontsize=8, va="top", color=P.INK_SECONDARY)
    ax.grid(True, axis="both")
    ax.tick_params(length=3, width=0.8)
    ax.legend(loc="lower right", fontsize=8)

    # Panel c: consolidated candidates, three estimators side by side.
    ax = fig.add_subplot(gs[1, :])
    c = cand.sort_values(["checks_passed", "genus_consistency"],
                         ascending=[True, True]).reset_index(drop=True)
    y = np.arange(len(c))
    ax.axvline(0, color=P.AXIS_RULE, linewidth=1.0, zorder=1)
    series = [
        ("adjusted_odds_ratio", P.BLUE, "o", "Lineage-adjusted regression", +0.24),
        ("shrunk_odds_ratio", P.ORANGE, "s", "After shrinkage across systems", 0.0),
        ("within_cluster_odds_ratio", P.AQUA, "D", "Within clonal cluster", -0.24),
    ]
    for col, colour, marker, label, off in series:
        v = np.log2(c[col].to_numpy(dtype=float))
        ok = np.isfinite(v)
        ax.scatter(v[ok], y[ok] + off, s=42, facecolor=colour,
                   edgecolor=P.SURFACE, linewidth=0.9, marker=marker,
                   label=label, zorder=3)
    # Anchor the note to the left edge of the panel rather than to a position
    # on the odds ratio scale, where it would read as an estimate.
    from matplotlib.transforms import blended_transform_factory
    edge = blended_transform_factory(ax.transAxes, ax.transData)
    for i, r in enumerate(c.itertuples(index=False)):
        if not np.isfinite(r.within_cluster_odds_ratio):
            ax.text(0.012, i - 0.24, "within-cluster not estimable",
                    fontsize=7.5, va="center", ha="left",
                    color=P.INK_MUTED, transform=edge)
    ax.set_yticks(y)
    ax.set_yticklabels([
        f"{r.system}   ({r.n_carrier_strains} strains, {r.checks_passed}/5 checks, "
        f"{r.genera_below_one}/{r.genera_estimable} genera)"
        for r in c.itertuples(index=False)], fontsize=8.5)
    ax.set_ylim(-0.7, len(c) - 0.3)
    ax.set_xlabel("log2 odds ratio for productive lysis")
    ax.set_title("c  Every candidate the study produced, on one scale", loc="left")
    P.tidy(ax, xgrid=True)
    ax.text(0.012, 0.015, "Favours resistance", transform=ax.transAxes,
            fontsize=8.5, color=P.INK_MUTED, ha="left", va="bottom")
    ax.text(0.988, 0.015, "Favours lysis", transform=ax.transAxes,
            fontsize=8.5, color=P.INK_MUTED, ha="right", va="bottom")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=3)

    surf = shap.set_index("block").loc["host surface", "shapley_pseudo_r2"]
    defn = shap.set_index("block").loc["defense systems", "shapley_pseudo_r2"]
    top = int(c["checks_passed"].max())
    fig.suptitle("Figure 7. Host surface outweighs defense, and no candidate "
                 "system survives every check",
                 x=0.035, y=0.982, ha="left", fontsize=13, fontweight="semibold")
    fig.text(0.035, 0.961,
             f"Panel a allocates explained variation by exact Shapley value, so it does not depend on the order the "
             f"blocks are entered: host surface takes {surf:.3f} against {defn:.3f} for all 33 defense systems, "
             "reversing the impression given by the ordered decomposition in script 03.\n"
             "Panel b shows each system pulled toward the pooled mean in orange by an amount set by its own "
             f"precision. Panel c collects every candidate nominated anywhere in the study; the best pass {top} of 5 "
             "independent checks and none passes all five.",
             ha="left", va="top", fontsize=9, color=P.INK_SECONDARY)
    fig.subplots_adjust(top=0.9, bottom=0.085, left=0.2)

    written = P.save(fig, paths.FIGURES, "fig7_consolidated_evidence",
                     config.FIGURE_DPI)
    for p in written:
        print(f"  {p.stat().st_size / 1024:8.0f} KiB  {p.relative_to(paths.ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
