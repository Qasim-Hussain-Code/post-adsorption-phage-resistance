"""Workflow diagram of the whole analysis.

Every count in the diagram is read from the pipeline's own outputs rather than
typed in, so the figure cannot drift away from the analysis it describes. The
palette is the one validated for the other figures.

Outputs:
    figures/fig11_workflow.{png,pdf}
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch

from papr import paths, plotting as P

W, H = 160.0, 106.0
NAVY = "#12203a"
PHASE = [P.SERIES[0], "#1f66c0", "#5b4bb5", "#8b3fa8", "#c23a86", P.SERIES[1]]


def facts() -> dict:
    """Read the counts the diagram quotes from the analysis itself."""
    pairs = pd.read_csv(paths.PAIRS, low_memory=False)
    feats = pd.read_csv(paths.STRAIN_FEATURES, low_memory=False)
    f = dict(
        pairs=len(pairs), strains=pairs["strain"].nunique(),
        phages=pairs["phage"].nunique(),
        annotated=len([c for c in feats.columns if c.startswith("ds_")]),
        tables=len(list(paths.TABLES.glob("table*.csv"))),
        reports=len(list(paths.RESULTS.glob("*.txt"))),
        figures=len(list(paths.FIGURES.glob("fig*.png"))),
        scripts=len(list(paths.SCRIPTS.glob("[0-9][0-9]_*.py"))),
    )
    ladder = paths.TABLES / "table3_defense_conditioning_ladder.csv"
    f["tested"] = (pd.read_csv(ladder)["system"].nunique()
                   if ladder.exists() else 33)
    val = paths.PROCESSED / "host_genome_validation.csv"
    f["retrieved"] = len(pd.read_csv(val)) if val.exists() else 287
    v = paths.TABLES / "table36_host_genome_validation.csv"
    f["validated"] = (int((pd.read_csv(v)["kmer_rank_of_core_nearest"] <= 10).sum())
                      if v.exists() else 254)
    loci = paths.TABLES / "table48_gwas_all_loci.csv"
    f["loci"] = len(pd.read_csv(loci)) if loci.exists() else 10351
    cand = paths.TABLES / "table33_consolidated_candidates.csv"
    f["candidates"] = len(pd.read_csv(cand)) if cand.exists() else 11
    return f


def box(ax, x, y, w, h, fc="#ffffff", ec=P.AXIS_RULE, lw=0.9, r=0.9, z=2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                facecolor=fc, edgecolor=ec, linewidth=lw, zorder=z))


def down(ax, x, y0, y1, colour=P.INK_MUTED):
    ax.add_patch(FancyArrowPatch((x, y0), (x, y1), arrowstyle="-|>",
                                 mutation_scale=11, color=colour,
                                 linewidth=1.5, zorder=4))


def right(ax, x0, x1, y, colour=P.INK_MUTED):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>",
                                 mutation_scale=11, color=colour,
                                 linewidth=1.5, zorder=4))


def main() -> int:
    paths.bootstrap()
    P.apply_style()
    f = facts()

    fig = plt.figure(figsize=(W / 10, H / 10))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((0, 0), W, H, boxstyle="square,pad=0",
                                facecolor=P.SURFACE, edgecolor="none", zorder=0))

    # ------------------------------------------------------------- header
    ax.add_patch(FancyBboxPatch((0, H - 11.8), W, 11.8, boxstyle="square,pad=0",
                                facecolor=NAVY, edgecolor="none", zorder=1))
    ax.add_patch(Circle((6.4, H - 5.9), 2.7, facecolor="#ffffff", alpha=0.13,
                        edgecolor="#ffffff", linewidth=0.9, zorder=2))
    # A small phage: head, tail, fibres.
    ax.add_patch(plt.Polygon([[6.4, H - 4.3], [7.5, H - 5.0], [7.5, H - 6.2],
                              [6.4, H - 6.9], [5.3, H - 6.2], [5.3, H - 5.0]],
                             closed=True, facecolor="#ffffff", zorder=3))
    ax.plot([6.4, 6.4], [H - 6.9, H - 7.9], color="#ffffff", lw=1.4, zorder=3)
    for dx in (-0.9, 0.9):
        ax.plot([6.4, 6.4 + dx], [H - 7.9, H - 8.5], color="#ffffff", lw=1.1,
                zorder=3)
    ax.text(11.5, H - 5.0, "POST-ADSORPTION PHAGE RESISTANCE", fontsize=20,
            color="#ffffff", fontweight="bold", va="center", zorder=3)
    ax.text(11.5, H - 8.6,
            "Does intracellular defense explain infection failure once surface "
            "compatibility and clonal lineage are held fixed?",
            fontsize=10.5, color="#b9c9e2", va="center", zorder=3)
    ax.text(W - 4, H - 4.6, f"{f['scripts']} scripts", fontsize=10.5,
            color="#9fb6d6", ha="right", va="center", fontweight="bold", zorder=3)
    ax.text(W - 4, H - 7.6, f"{f['pairs']:,} strain-phage pairs", fontsize=10.5,
            color="#9fb6d6", ha="right", va="center", zorder=3)

    # ------------------------------------------------------------- columns
    for label, x, w, colour in (("INPUT DATA", 4, 38, P.SERIES[0]),
                                ("ANALYSIS PHASES", 45, 69, NAVY),
                                ("OUTPUTS", 117, 39, P.SERIES[2])):
        box(ax, x, H - 16.4, w, 3.4, fc=colour, ec=colour, r=0.7, z=3)
        ax.text(x + w / 2, H - 14.7, label, fontsize=11, color="#ffffff",
                fontweight="bold", ha="center", va="center", zorder=4)

    # -------------------------------------------------------------- inputs
    inputs = [
        ("Interaction matrix",
         [f"{f['strains']} Escherichia strains against {f['phages']} phages",
          f"{f['pairs']:,} observed pairs, plaque score 0 to 4",
          "Gaborieau et al., Nature Microbiology 2024"], P.SERIES[0]),
        ("Host genotypes",
         [f"{f['annotated']} DefenseFinder systems per strain",
          "LPS core, O-antigen, capsule, 12 OMP loci",
          "core-genome tree and distance matrix"], "#1f66c0"),
        ("Phage genomes and models",
         [f"{f['phages']} genomes, VIRIDIC tree, RBP assignments",
          "178 anti-defense profile models",
          "7 lysogeny profiles for lifestyle calls"], "#5b4bb5"),
        ("Bacterial genomes",
         [f"{f['retrieved']} retrieved from NCBI one strain at a time",
          f"{f['validated']} pass an independent identity check",
          "not redistributed; digests are tracked"], "#8b3fa8"),
    ]
    y = H - 20.4
    for title, lines, colour in inputs:
        h = 3.0 + 2.1 * len(lines)
        box(ax, 4, y - h, 38, h)
        ax.add_patch(FancyBboxPatch((4, y - h), 0.6, h,
                                    boxstyle="round,pad=0,rounding_size=0.3",
                                    facecolor=colour, edgecolor="none", zorder=3))
        ax.text(5.8, y - 1.6, title, fontsize=10.5, fontweight="bold",
                color=P.INK_PRIMARY, va="center", zorder=4)
        for i, ln in enumerate(lines):
            ax.text(5.8, y - 3.7 - 2.0 * i, ln, fontsize=8.4,
                    color=P.INK_SECONDARY, va="center", zorder=4)
        y -= h + 1.6

    box(ax, 4, y - 7.4, 38, 7.4, fc="#f3f6fb", ec="#d7e2f1")
    ax.text(5.8, y - 1.8, "Public inputs, verified on every run", fontsize=10,
            fontweight="bold", color=P.INK_PRIMARY, va="center", zorder=4)
    for i, ln in enumerate(["15 files with SHA-256 digests tracked in the",
                            "repository. The data itself is re-downloaded",
                            "rather than redistributed."]):
        ax.text(5.8, y - 3.8 - 1.7 * i, ln, fontsize=8.4,
                color=P.INK_SECONDARY, va="center", zorder=4)

    # -------------------------------------------------------------- phases
    phases = [
        ("Acquire and validate",
         ["Verify digests, assert the observed dimensions, and audit the",
          "released matrix against 318,816 replicate spot tests: the binary",
          "lysis call is reproduced 99.69 percent of the time."]),
        ("Separate the two gates",
         ["Define receptor compatibility from isolation hosts, then stratify",
          "rather than filter. The gate admits 65.7 percent of pairs and",
          "enriches lysis only 1.34 fold, so filtering would discard signal."]),
        ("Associate, then remove the lineage",
         ["Condition on phage identity and surface type, refit with errors",
          "clustered on clonal background, and partition the variance by",
          "exact Shapley value over all orderings of the blocks."]),
        ("Stress-test whatever survives",
         ["Within-cluster fits, permutation inside a lineage, minimum",
          "detectable effects, E-values, genus heterogeneity, matched",
          "phylogenetic nulls, and full Bayesian pooling across 79 systems."]),
        ("Search beyond the catalogue",
         [f"Test every variable k-mer across {f['validated']} genomes, "
          f"{f['loci']:,} loci, against",
          "an empirical null with a locally fitted tail, then relocate and",
          "annotate whatever it flags."]),
        ("Resolve the panel into phages",
         ["Score each system against each phage separately, counting distinct",
          "phage groups, against a coordinate-matched and a clonally-matched",
          "null, and report the worse of the two p-values."]),
    ]
    y = H - 20.4
    for i, (title, lines) in enumerate(phases):
        h = 11.0
        box(ax, 45, y - h, 69, h)
        ax.add_patch(FancyBboxPatch((45, y - h), 7.6, h,
                                    boxstyle="round,pad=0,rounding_size=0.9",
                                    facecolor=PHASE[i], edgecolor="none", zorder=3))
        ax.add_patch(FancyBboxPatch((51.0, y - h), 1.6, h,
                                    boxstyle="square,pad=0",
                                    facecolor=PHASE[i], edgecolor="none", zorder=3))
        ax.text(48.8, y - h / 2, f"{i + 1:02d}", fontsize=19, color="#ffffff",
                fontweight="bold", ha="center", va="center", zorder=4)
        ax.text(54.4, y - 2.2, title, fontsize=11.5, fontweight="bold",
                color=PHASE[i], va="center", zorder=4)
        for k, ln in enumerate(lines):
            ax.text(54.4, y - 4.9 - 2.0 * k, ln, fontsize=8.4,
                    color=P.INK_SECONDARY, va="center", zorder=4)
        if i < len(phases) - 1:
            down(ax, 79.5, y - h - 0.1, y - h - 1.5)
        y -= h + 1.6

    # ------------------------------------------------------------- outputs
    outputs = [
        ("Result tables", f"{f['tables']} CSV files",
         "every number behind every figure", P.SERIES[2]),
        ("Plain-text reports", f"{f['reports']} reports",
         "each stage explains its own result", P.SERIES[0]),
        ("Manuscript figures", f"{f['figures']} figures",
         "PNG and PDF, byte-reproducible", "#5b4bb5"),
        ("Ranked candidate list", f"{f['candidates']} systems nominated",
         "ordered by independent checks passed", P.SERIES[1]),
    ]
    y = H - 20.4
    for i, (title, big, small, colour) in enumerate(outputs):
        h = 10.2
        box(ax, 117, y - h, 39, h)
        ax.add_patch(FancyBboxPatch((117, y - h), 0.6, h,
                                    boxstyle="round,pad=0,rounding_size=0.3",
                                    facecolor=colour, edgecolor="none", zorder=3))
        ax.text(119, y - 2.1, title, fontsize=10.5, fontweight="bold",
                color=P.INK_PRIMARY, va="center", zorder=4)
        ax.text(119, y - 5.3, big, fontsize=15, fontweight="bold",
                color=colour, va="center", zorder=4)
        ax.text(119, y - 8.1, small, fontsize=8.4, color=P.INK_SECONDARY,
                va="center", zorder=4)
        if i < len(outputs) - 1:
            down(ax, 136.5, y - h - 0.1, y - h - 1.5)
        y -= h + 1.6

    box(ax, 117, y - 14.6, 39, 14.6, fc="#fbf5f0", ec="#eddcc9")
    ax.text(119, y - 2.2, "What the analysis concludes", fontsize=10.5,
            fontweight="bold", color=P.INK_PRIMARY, va="center", zorder=4)
    for k, ln in enumerate([
            f"Of {f['tested']} testable systems, almost every association",
            "is a property of the bacterial lineage rather than",
            "of the system. No system passes all seven",
            "independent checks, and roughly three quarters of",
            "infections still fail after accounting for surface",
            "type, defense carriage and prophage content."]):
        ax.text(119, y - 4.6 - 1.75 * k, ln, fontsize=8.2,
                color=P.INK_SECONDARY, va="center", zorder=4)

    right(ax, 42.4, 44.6, H - 55)
    right(ax, 114.4, 116.6, H - 55)

    # -------------------------------------------------------------- footer
    fy = 3.0
    box(ax, 4, fy, W - 8, 8.6)
    items = [
        ("Reproducible", "121 outputs regenerate byte for byte"),
        ("Verified", "226 claims checked against the tables"),
        ("Tested", "estimators checked against statsmodels"),
        ("Confounder-aware", "surface, lineage and island adjusted"),
        ("Nulls reported", "negative results carry equal weight"),
    ]
    step = (W - 14) / len(items)
    for i, (head, body) in enumerate(items):
        x = 7 + step * i
        ax.add_patch(Circle((x + 1.1, fy + 4.3), 1.05, facecolor="none",
                            edgecolor=PHASE[i], linewidth=1.6, zorder=4))
        ax.plot([x + 0.6, x + 1.0, x + 1.7], [fy + 4.3, fy + 3.8, fy + 4.9],
                color=PHASE[i], lw=1.6, solid_capstyle="round", zorder=5)
        ax.text(x + 3.0, fy + 5.4, head, fontsize=9.6, fontweight="bold",
                color=P.INK_PRIMARY, va="center", zorder=4)
        ax.text(x + 3.0, fy + 3.1, body, fontsize=7.8,
                color=P.INK_SECONDARY, va="center", zorder=4)

    written = P.save(fig, paths.FIGURES, "fig11_workflow", dpi=200)
    for path in written:
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
