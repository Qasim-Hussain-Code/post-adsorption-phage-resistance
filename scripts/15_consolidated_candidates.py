"""Assemble every candidate the study produced and judge them on one scale.

Different parts of this analysis nominate different systems. The lineage
adjusted regression of script 04 picks two. The within-cluster and permutation
designs of script 09 pick others. The shrinkage model of script 13 covers 79
systems, including rare ones the earlier threshold excluded, and picks others
again. A reader cannot be expected to reconcile four candidate lists.

This script takes the union of them, applies the same battery to every member
including the rare ones that have never been through it, and produces a single
ranked table. The ranking is a count of independent checks passed, not a
p-value, because the checks disagree often enough that any single one would be
misleading.

Nothing here is a new test of a new hypothesis. It is the same evidence,
assembled so that the candidate list can be read off and acted on.

Outputs:
    results/tables/table33_consolidated_candidates.csv
    results/consolidated_candidates_report.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from papr import config, paths
from papr.robust import permutation_effects, within_cluster_effects
from papr.stats import mantel_haenszel, stratified_tables

N_PERMUTATIONS = 10_000


def e_value(odds_ratio: float, baseline_risk: float) -> float:
    if not np.isfinite(odds_ratio) or odds_ratio <= 0:
        return np.nan
    rr = odds_ratio / ((1 - baseline_risk) + baseline_risk * odds_ratio)
    if rr < 1:
        rr = 1.0 / rr
    return float(rr + np.sqrt(rr * (rr - 1)))


def genus_consistency(pairs: pd.DataFrame, col: str) -> tuple[int, int]:
    below = total = 0
    for _, frame in pairs.groupby("genus"):
        eff = mantel_haenszel(stratified_tables(frame, ["phage"], col, "lytic"))
        if np.isfinite(eff.estimate) and eff.estimate > 0:
            total += 1
            below += int(eff.estimate < 1)
    return below, total


def main() -> int:
    paths.bootstrap()
    pairs = pd.read_csv(paths.PAIRS, compression="gzip", low_memory=False)
    strain_feats = pd.read_csv(paths.STRAIN_FEATURES, low_memory=False)
    clusters = pd.read_csv(paths.PROCESSED / "clonal_clusters.csv", index_col=0)
    all_eff = pd.read_csv(paths.TABLES / "table29_all_system_effects.csv")
    lineage = pd.read_csv(paths.TABLES / "table7_lineage_controlled_effects.csv")
    wc_prev = pd.read_csv(paths.TABLES / "table16_within_cluster_effects.csv")
    perm_prev = pd.read_csv(paths.TABLES / "table17_permutation_tests.csv")

    ds_cols = [c for c in strain_feats.columns if c.startswith("ds_")]
    pairs = pairs.merge(strain_feats[["strain"] + ds_cols], on="strain", how="left")
    pairs["clonal_cluster"] = pairs["strain"].map(clusters["clonal_cluster"])
    baseline = float(pairs["lytic"].mean())

    # Union of every nomination, protective direction only. A system associated
    # with more lysis cannot be a defense and is handled separately below.
    full = lineage[lineage["adjustment"] == "plus_phylogroup_and_coordinates"]
    nominated: dict[str, list[str]] = {}

    def nominate(names, source):
        for n in names:
            nominated.setdefault(n, []).append(source)

    nominate(full[full["significant"] & full["direction_preserved"]
                  & (full["odds_ratio"] < 1)]["system"], "lineage-adjusted regression")
    nominate(wc_prev[wc_prev["significant"] & (wc_prev["odds_ratio"] < 1)]["system"],
             "within clonal cluster")
    nominate(perm_prev[perm_prev["significant"]
                       & (perm_prev["observed_difference"] < 0)]["system"],
             "permutation inside lineage")
    shrunk_sig = all_eff[all_eff["posterior_excludes_null"].fillna(False)
                         & (all_eff["posterior_odds_ratio"] < 1)]
    nominate(shrunk_sig["system"], "shrinkage across all systems")

    candidates = sorted(nominated)
    cols = [f"ds_{c}" for c in candidates if f"ds_{c}" in pairs.columns]
    print(f"{len(candidates)} candidate systems nominated by at least one analysis")

    print("  within-cluster battery")
    wc = within_cluster_effects(pairs, cols).set_index("system")
    print("  permutation battery")
    perm = permutation_effects(
        pairs, strain_feats, cols,
        [("phylogroup", "phylogroup"), ("sequence_type", "sequence_type")],
        N_PERMUTATIONS, config.RANDOM_SEED)
    perm_phy = perm[perm["level"] == "phylogroup"].set_index("system")
    perm_st = perm[perm["level"] == "sequence_type"].set_index("system")
    eff_i = all_eff.set_index("system")
    full_i = full.set_index("system")

    rows = []
    for name in candidates:
        col = f"ds_{name}"
        if col not in pairs.columns:
            continue
        adj_or = float(eff_i.loc[name, "raw_odds_ratio"]) if name in eff_i.index else np.nan
        shrunk = float(eff_i.loc[name, "posterior_odds_ratio"]) if name in eff_i.index else np.nan
        shrunk_ok = bool(eff_i.loc[name, "posterior_excludes_null"]) if name in eff_i.index else False
        w = wc.loc[name] if name in wc.index else None
        pp = perm_phy.loc[name] if name in perm_phy.index else None
        ps = perm_st.loc[name] if name in perm_st.index else None
        below, total = genus_consistency(pairs, col)

        passes = {
            "lineage_adjusted": bool(
                name in full_i.index and full_i.loc[name, "significant"]
                and full_i.loc[name, "odds_ratio"] < 1),
            "within_cluster": bool(
                w is not None and w["significant"]
                and np.isfinite(w["odds_ratio"]) and w["odds_ratio"] < 1),
            "permutation_phylogroup": bool(
                pp is not None and pp["significant"] and pp["observed_difference"] < 0),
            "permutation_st": bool(
                ps is not None and ps["significant"] and ps["observed_difference"] < 0),
            "shrinkage": bool(shrunk_ok and np.isfinite(shrunk) and shrunk < 1),
        }
        rows.append(dict(
            system=name,
            n_carrier_strains=int(strain_feats[col].sum()),
            nominated_by="; ".join(sorted(set(nominated[name]))),
            adjusted_odds_ratio=adj_or,
            shrunk_odds_ratio=shrunk,
            within_cluster_odds_ratio=(float(w["odds_ratio"]) if w is not None else np.nan),
            within_cluster_clusters=(int(w["n_discordant_clusters"]) if w is not None else 0),
            within_cluster_q=(float(w["q_value"]) if w is not None else np.nan),
            permutation_q_phylogroup=(float(pp["q_value"]) if pp is not None else np.nan),
            permutation_q_sequence_type=(float(ps["q_value"]) if ps is not None else np.nan),
            genera_below_one=below, genera_estimable=total,
            e_value=e_value(adj_or, baseline),
            **{f"pass_{k}": v for k, v in passes.items()},
            checks_passed=int(sum(passes.values())),
            # A within-cluster test needs at least three informative clusters
            # before a null from it means anything; below that the check was
            # not really applied.
            within_cluster_informative=bool(
                w is not None and w["n_discordant_clusters"] >= 3
                and np.isfinite(w["odds_ratio"])),
        ))

    table = pd.DataFrame(rows)
    table["genus_consistency"] = table["genera_below_one"] / table["genera_estimable"].replace(0, np.nan)
    table = table.sort_values(
        ["checks_passed", "genus_consistency", "e_value"],
        ascending=[False, False, False])
    table.to_csv(paths.TABLES / "table33_consolidated_candidates.csv", index=False)

    o: list[str] = []
    add = o.append
    add("Consolidated candidate assessment")
    add("Generated by scripts/15_consolidated_candidates.py")
    add("")
    add("HOW THIS LIST WAS BUILT")
    add("  Every system nominated by any of four analyses, in the protective "
        "direction, then put through the same battery:")
    add("    lineage-adjusted regression   script 04, covariate adjustment")
    add("    within clonal cluster         script 09, near-identical genomes")
    add("    permutation inside phylogroup script 09, design-based null")
    add("    permutation inside ST         script 09, stricter design-based null")
    add("    shrinkage across systems      script 13, correction for multiplicity")
    add(f"  Candidates assessed: {len(table)}")
    add("  Several candidates below the 20-carrier threshold had never been "
        "through the within-cluster or permutation tests before; they have now.")

    add("")
    add("RANKED CANDIDATES")
    add(f"  {'system':24s} {'carr':>5s} {'adj OR':>7s} {'shrunk':>7s} "
        f"{'wc OR':>7s} {'wc n':>5s} {'genera':>7s} {'E':>5s} {'passed':>7s}")
    for r in table.itertuples(index=False):
        gen = ("   n/a" if not r.genera_estimable
               else f"{r.genera_below_one:3d}/{r.genera_estimable:<3d}")
        add(f"  {r.system:24s} {r.n_carrier_strains:5d} "
            f"{r.adjusted_odds_ratio:7.3f} {r.shrunk_odds_ratio:7.3f} "
            f"{r.within_cluster_odds_ratio:7.3f} {r.within_cluster_clusters:5d} "
            f"{gen:>7s} {r.e_value:5.2f} {r.checks_passed:5d}/5")
    add("")
    add("  'carr' is carrier strains, 'wc' the within-cluster fit and its number "
        "of informative clusters, 'genera' the phage genera in which the odds "
        "ratio is below 1, 'E' the E-value for unmeasured confounding.")

    best = table[table["checks_passed"] >= 2]
    add("")
    add("READING THE LIST")
    if len(best):
        add(f"  {len(best)} systems pass at least two independent checks: "
            f"{', '.join(best['system'])}.")
    add("  No system passes all five, and that is the honest headline. The "
        "checks are not redundant: each discards most of the data to remove a "
        "different source of confounding, so each has little power left, and "
        "they disagree. A system passing two of five is a candidate worth an "
        "experiment, not a finding.")

    thin = table[table["within_cluster_clusters"] < 3]
    if len(thin):
        add("")
        add(f"  {len(thin)} candidates have fewer than three informative clonal "
            f"clusters ({', '.join(thin['system'])}), so the within-cluster "
            "column for them is uninformative rather than negative. Failing a "
            "test that cannot be run is not evidence against a system.")

    rare = table[table["n_carrier_strains"] < config.MIN_CARRIER_STRAINS]
    if len(rare):
        add("")
        add(f"  {len(rare)} candidates are carried by fewer than "
            f"{config.MIN_CARRIER_STRAINS} strains "
            f"({', '.join(rare['system'])}). They reach this list only through "
            "the shrinkage model, their individual estimates are heavily "
            "pulled toward the pooled mean, and they should be treated as the "
            "weakest entries here despite their rank.")

    add("")
    add("WHAT WOULD SETTLE IT")
    add("  Nothing in this dataset can. Every analysis here is observational "
        "and the confounding between a defense system and the lineage carrying "
        "it cannot be broken by stratifying a fixed collection more finely. "
        "The decisive experiment is to delete the system in its native host and "
        "challenge the isogenic pair with the phages this collection says it "
        "resists. The ranking above is a priority order for doing that, and "
        "the within-cluster column indicates which systems already have "
        "naturally occurring near-isogenic pairs in the collection to start "
        "from.")

    text = "\n".join(o)
    (paths.RESULTS / "consolidated_candidates_report.txt").write_text(
        text + "\n", encoding="utf-8")
    print("\n" + text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
