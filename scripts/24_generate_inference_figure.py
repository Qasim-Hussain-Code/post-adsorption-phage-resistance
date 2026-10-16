"""Figure for the extended inference: structure, nulls, posteriors, verdict.

Palette and style are those validated for script 06; no new colour slots are
introduced, so no re-validation is required.

Outputs:
    figures/fig8_inference_and_verdict.{png,pdf}
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


def main() -> int:
    paths.bootstrap()
    P.apply_style()

    phylo = pd.read_csv(paths.TABLES / "table40_phylogenetic_null_tests.csv")
    bayes = pd.read_csv(paths.TABLES / "table41_bayesian_pooled_effects.csv")
    cand = pd.read_csv(paths.TABLES / "table33_consolidated_candidates.csv")

    fig = plt.figure(figsize=(14.2, 12.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.26, wspace=0.3)

    # a: how lineage-confined a system is, against how much that inflates the null
    ax = fig.add_subplot(gs[0, 0])
    ax.scatter(phylo["phylogenetic_r2"], phylo["inflation_factor"], s=34,
               facecolor=P.BLUE, edgecolor=P.SURFACE, linewidth=0.9, zorder=3)
    # Alternate the label offsets so that neighbouring points do not collide.
    labelled = [r for r in phylo.itertuples(index=False)
                if r.phylogenetic_r2 > 0.15 or r.inflation_factor > 1.3]
    labelled.sort(key=lambda r: r.phylogenetic_r2)
    for i, r in enumerate(labelled):
        dy = 9 if i % 2 == 0 else -11
        ax.annotate(r.system, (r.phylogenetic_r2, r.inflation_factor),
                    textcoords="offset points", xytext=(8, dy), fontsize=7.5,
                    color=P.INK_PRIMARY)
    ax.axhline(1.0, color=P.AXIS_RULE, linewidth=1.0, zorder=1)
    ax.set_xlabel("Share of carriage explained by core-genome position")
    ax.set_ylabel("Inflation of the null standard deviation")
    ax.set_title("a  The size of the population structure problem", loc="left")
    ax.grid(True, axis="both")
    ax.tick_params(length=3, width=0.8)

    # b: phylogenetically matched null, before and after removing surface type
    ax = fig.add_subplot(gs[0, 1])
    lo = 1e-4
    x = np.maximum(phylo["q_structured"].to_numpy(), lo)
    y = np.maximum(phylo["q_structured_surface_adjusted"].to_numpy(), lo)
    prot = (phylo["observed_difference_surface_adjusted"] < 0).to_numpy()
    ax.plot([lo, 1], [lo, 1], color=P.INK_MUTED, linewidth=1.0, alpha=0.6, zorder=1)
    for mask, colour, marker, label in [
        (prot, P.BLUE, "o", "Associated with less lysis"),
        (~prot, P.ORANGE, "s", "Associated with more lysis"),
    ]:
        ax.scatter(x[mask], y[mask], s=34, facecolor=colour,
                   edgecolor=P.SURFACE, linewidth=0.9, marker=marker,
                   label=label, zorder=3)
    marked = [r for r in phylo.itertuples(index=False)
              if min(r.q_structured, r.q_structured_surface_adjusted) <= 0.05]
    marked.sort(key=lambda r: r.q_structured)
    for i, r in enumerate(marked):
        dx, dy = (9, 7) if i % 2 == 0 else (9, -12)
        ax.annotate(r.system,
                    (max(r.q_structured, lo),
                     max(r.q_structured_surface_adjusted, lo)),
                    textcoords="offset points", xytext=(dx, dy), fontsize=7.5,
                    color=P.INK_PRIMARY)
    for v in (0.05,):
        ax.axhline(v, color=P.AXIS_RULE, linewidth=0.9, zorder=1)
        ax.axvline(v, color=P.AXIS_RULE, linewidth=0.9, zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, 1.3)
    ax.set_ylim(lo, 1.3)
    ax.set_xlabel("q against a phylogenetically matched null")
    ax.set_ylabel("q with host surface type removed")
    ax.set_title("b  Phylogeny alone is not enough", loc="left")
    ax.grid(True, axis="both")
    ax.tick_params(length=3, width=0.8)
    ax.legend(loc="lower right", fontsize=8)

    # c: consolidated verdict, checks passed with the posterior beside it
    ax = fig.add_subplot(gs[1, :])
    c = cand.sort_values(["checks_passed", "prob_protective"],
                         ascending=[True, True]).reset_index(drop=True)
    y = np.arange(len(c))
    ax.barh(y, c["checks_passed"], height=0.62, color=P.BLUE,
            edgecolor=P.SURFACE, linewidth=0.9, zorder=2)
    for i, r in enumerate(c.itertuples(index=False)):
        pp = ("" if not np.isfinite(r.prob_protective)
              else f"   P(protective) {r.prob_protective:.3f}")
        wc = ("within-cluster not estimable"
              if not np.isfinite(r.within_cluster_odds_ratio)
              else f"within-cluster OR {r.within_cluster_odds_ratio:.2f} "
                   f"on {r.within_cluster_clusters} clusters")
        ax.text(r.checks_passed + 0.08, i, f"{wc}{pp}", va="center",
                fontsize=7.8, color=P.INK_SECONDARY)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.system}  ({r.n_carrier_strains} strains)"
                        for r in c.itertuples(index=False)], fontsize=8.5)
    ax.set_xlim(0, 7.9)
    ax.set_xlabel("Independent checks passed, of seven")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_title("c  Every candidate the study produced, and how much supports it",
                 loc="left")
    P.tidy(ax, xgrid=True)

    best = int(cand["checks_passed"].max())
    n_struct = int((phylo["q_structured"] <= 0.05).sum())
    n_surf = int((phylo["q_structured_surface_adjusted"] <= 0.05).sum())
    fig.suptitle("Figure 8. Harder nulls, proper posteriors, and a candidate "
                 "list that no single test supports",
                 x=0.035, y=0.983, ha="left", fontsize=13, fontweight="semibold")
    fig.text(0.035, 0.962,
             "Panel a shows that a system confined to one lineage widens the null for any test that ignores ancestry, "
             f"by up to a factor of {phylo['inflation_factor'].max():.1f} for "
             f"{phylo.loc[phylo['inflation_factor'].idxmax(), 'system']}.\n"
             f"Panel b compares a null matched on ancestry ({n_struct} systems significant) with the same null after "
             f"host surface type is removed from the outcome ({n_surf} significant); the associations with more lysis, "
             "which no defense can cause, are the ones that move. Panel c collects every candidate; the best passes "
             f"{best} of seven checks.",
             ha="left", va="top", fontsize=9, color=P.INK_SECONDARY)
    fig.subplots_adjust(top=0.9, bottom=0.07, left=0.22)

    written = P.save(fig, paths.FIGURES, "fig8_inference_and_verdict",
                     config.FIGURE_DPI)
    for p in written:
        print(f"  {p.stat().st_size / 1024:8.0f} KiB  {p.relative_to(paths.ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
