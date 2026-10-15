"""Does any pair of defense systems protect when neither does alone.

Carrying more defense systems does not make a strain harder to lyse: script 03
found no monotone relationship between the number of systems and the lysis
rate, which is one of the cleaner null results in this analysis. That null
assumes systems act independently and additively. It would be missed entirely
if protection came from specific combinations rather than from counting, so
this script tests the pairs directly.

Each pair is fitted with both main effects, their interaction, phage identity,
host surface type and the clonal background covariates, with standard errors
clustered on clonal cluster. A pair is only tested when all four carriage
combinations occur in at least ten strains, since an interaction estimated from
an almost empty cell is not an estimate.

The script also asks a softer version of the same question: whether carrying
systems from several distinct mechanistic classes, rather than many systems of
one kind, is what matters.

Outputs:
    results/tables/table31_defense_pair_interactions.csv
    results/tables/table32_mechanistic_class_burden.csv
    results/defense_combinations_report.txt
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
from scipy import stats

from papr import config, paths
from papr.models import fit_logit_cluster, term_effect
from papr.stats import benjamini_hochberg

# Every one of the four carriage combinations needs this many strains.
MIN_CELL_STRAINS = 10

# Coarse mechanistic grouping of the annotated systems. Membership is decided
# by name prefix, which follows DefenseFinder's own naming.
CLASS_RULES: list[tuple[str, tuple[str, ...]]] = [
    ("restriction_modification", ("RM_Type", "Dnd", "Dpd", "BREX", "Dsr")),
    ("crispr", ("CAS_", "Cas")),
    ("abortive_infection", ("Abi", "Rex", "Rst", "Lit", "PrrC", "RloC", "Old_")),
    ("toxin_antitoxin", ("MazEF", "PsyrTA", "SanaTA", "RosmerTA", "ShosTA",
                         "MqsRAC", "DarTG", "Mok_Hok_Sok", "Retron")),
    ("nucleotide_signalling", ("CBASS", "Pycsar", "Thoeris", "Mokosh", "DRT",
                               "dGTPase", "dCTPdeaminase")),
    ("atpase_nlr_like", ("Avs_", "NLR_like", "Lamassu", "Septu", "Gabija",
                         "Hachiman", "Shedu", "Zorya", "Wadjet", "Kiwa",
                         "Druantia", "Menshen", "Shango", "Borvo", "Eleos",
                         "Olokun", "Dazbog", "Hna", "Gao_", "radar", "SoFic")),
]


def classify(system: str) -> str:
    for label, prefixes in CLASS_RULES:
        if any(system.startswith(p) for p in prefixes):
            return label
    return "other"


def main() -> int:
    paths.bootstrap()
    pairs = pd.read_csv(paths.PAIRS, compression="gzip", low_memory=False)
    strain_feats = pd.read_csv(paths.STRAIN_FEATURES, low_memory=False)
    pcoa = pd.read_csv(paths.PCOA, index_col=0)
    clusters = pd.read_csv(paths.PROCESSED / "clonal_clusters.csv", index_col=0)

    ds_cols = [c for c in strain_feats.columns if c.startswith("ds_")]
    tested = strain_feats[strain_feats["strain"].isin(set(pairs["strain"]))]
    systems = sorted([c for c in ds_cols
                      if tested[c].sum() >= config.MIN_CARRIER_STRAINS])
    pairs = pairs.merge(strain_feats[["strain"] + ds_cols], on="strain", how="left")
    pairs["clonal_cluster"] = pairs["strain"].map(clusters["clonal_cluster"])

    y = pairs["lytic"].to_numpy(dtype=float)
    groups = pairs["clonal_cluster"].to_numpy()
    coords = pcoa.reindex(pairs["strain"]).to_numpy(dtype=float)
    coords = (coords - coords.mean(axis=0)) / coords.std(axis=0)
    base = np.hstack([
        np.ones((len(pairs), 1)),
        pd.get_dummies(pairs["phage"].astype("string"), drop_first=True,
                       dtype=float).to_numpy(),
        pd.get_dummies(pairs["lps_core"].astype("string"), drop_first=True,
                       dtype=float).to_numpy(),
        pd.get_dummies(pairs["phylogroup"].astype("string"), drop_first=True,
                       dtype=float).to_numpy(),
        coords,
    ])

    candidates = []
    for a, b in itertools.combinations(systems, 2):
        cells = [int(((tested[a] == x) & (tested[b] == z)).sum())
                 for x, z in [(1, 1), (1, 0), (0, 1), (0, 0)]]
        if min(cells) >= MIN_CELL_STRAINS:
            candidates.append((a, b, cells))
    print(f"testing {len(candidates)} of "
          f"{len(systems) * (len(systems) - 1) // 2} pairs")

    rows = []
    for i, (a, b, cells) in enumerate(candidates, 1):
        xa = pairs[a].to_numpy(dtype=float)
        xb = pairs[b].to_numpy(dtype=float)
        x = np.hstack([(xa * xb).reshape(-1, 1), xa.reshape(-1, 1),
                       xb.reshape(-1, 1), base])
        try:
            res = fit_logit_cluster(y, x, groups)
            inter = term_effect(res, 0)
            ea = term_effect(res, 1)
            eb = term_effect(res, 2)
        except Exception:
            continue
        rows.append(dict(
            system_a=a[3:], system_b=b[3:],
            n_both=cells[0], n_a_only=cells[1], n_b_only=cells[2],
            n_neither=cells[3],
            or_a_alone=ea["odds_ratio"], or_b_alone=eb["odds_ratio"],
            interaction_odds_ratio=inter["odds_ratio"],
            interaction_ci_low=inter["ci_low"],
            interaction_ci_high=inter["ci_high"],
            interaction_p_value=inter["p_value"],
            combined_odds_ratio=float(
                ea["odds_ratio"] * eb["odds_ratio"] * inter["odds_ratio"])))
        if i % 25 == 0 or i == len(candidates):
            print(f"  {i}/{len(candidates)}")

    inter_tab = pd.DataFrame(rows)
    inter_tab["interaction_q_value"] = benjamini_hochberg(
        inter_tab["interaction_p_value"].to_numpy())
    inter_tab["significant"] = inter_tab["interaction_q_value"] <= config.FDR_ALPHA
    inter_tab = inter_tab.sort_values("interaction_odds_ratio")
    inter_tab.to_csv(paths.TABLES / "table31_defense_pair_interactions.csv",
                     index=False)

    # Mechanistic class burden against raw system count.
    class_of = {c: classify(c[3:]) for c in ds_cols}
    classes = sorted(set(class_of.values()))
    strain_class = pd.DataFrame(index=tested["strain"])
    for cl in classes:
        cols = [c for c in ds_cols if class_of[c] == cl]
        strain_class[cl] = (tested.set_index("strain")[cols].sum(axis=1) > 0).astype(int)
    strain_class["n_classes"] = strain_class[classes].sum(axis=1)
    strain_class["n_systems"] = tested.set_index("strain")[ds_cols].sum(axis=1)

    per_strain = pairs.groupby("strain")["lytic"].mean()
    joined = strain_class.join(per_strain.rename("lysis_rate"))
    rho_classes = stats.spearmanr(joined["n_classes"], joined["lysis_rate"])
    rho_systems = stats.spearmanr(joined["n_systems"], joined["lysis_rate"])

    class_rows = []
    for cl in classes:
        cols = [c for c in ds_cols if class_of[c] == cl]
        carriers = joined[joined[cl] == 1]
        others = joined[joined[cl] == 0]
        if len(carriers) < 5 or len(others) < 5:
            continue
        u = stats.mannwhitneyu(carriers["lysis_rate"], others["lysis_rate"],
                               alternative="two-sided")
        class_rows.append(dict(
            mechanistic_class=cl, n_systems_in_class=len(cols),
            n_strains_carrying=int(len(carriers)),
            mean_lysis_carriers=float(carriers["lysis_rate"].mean()),
            mean_lysis_noncarriers=float(others["lysis_rate"].mean()),
            mann_whitney_p=float(u.pvalue)))
    class_tab = pd.DataFrame(class_rows)
    if len(class_tab):
        class_tab["q_value"] = benjamini_hochberg(
            class_tab["mann_whitney_p"].to_numpy())
        class_tab["significant"] = class_tab["q_value"] <= config.FDR_ALPHA
    class_tab.to_csv(paths.TABLES / "table32_mechanistic_class_burden.csv",
                     index=False)

    o: list[str] = []
    add = o.append
    add("Defense system combinations report")
    add("Generated by scripts/14_defense_combinations.py")
    add("")
    add("PAIRWISE INTERACTIONS")
    add(f"  systems considered                 : {len(systems)}")
    add(f"  possible pairs                     : "
        f"{len(systems) * (len(systems) - 1) // 2}")
    add(f"  pairs with at least {MIN_CELL_STRAINS} strains in each of the four "
        f"carriage cells: {len(candidates)}")
    add(f"  pairs fitted                       : {len(inter_tab)}")
    nsig = int(inter_tab["significant"].sum())
    add(f"  interactions significant at q <= {config.FDR_ALPHA} : {nsig}")
    add("")
    if nsig:
        add(f"  {'pair':44s} {'both':>5s} {'interact':>8s} {'95 pct CI':>16s} "
            f"{'q':>9s}")
        for r in inter_tab[inter_tab["significant"]].itertuples(index=False):
            add(f"  {r.system_a + ' x ' + r.system_b:44s} {r.n_both:5d} "
                f"{r.interaction_odds_ratio:8.3f} "
                f"[{r.interaction_ci_low:6.3f}, {r.interaction_ci_high:6.3f}] "
                f"{r.interaction_q_value:9.2e}")
    else:
        add("  No pair of defense systems shows a detectable interaction. "
            "Carrying two systems is, as far as this analysis can tell, no "
            "better and no worse than what the two contribute separately.")
    add("")
    add("  Strongest interactions in each direction, significant or not:")
    add(f"  {'pair':44s} {'both':>5s} {'interact':>8s} {'q':>9s}")
    for r in pd.concat([inter_tab.head(4), inter_tab.tail(4)]).itertuples(index=False):
        add(f"  {r.system_a + ' x ' + r.system_b:44s} {r.n_both:5d} "
            f"{r.interaction_odds_ratio:8.3f} {r.interaction_q_value:9.2e}")

    add("")
    add("  An interaction below 1 would mean two systems together block lysis "
        "more than their separate effects predict, which is what synergy "
        "between a sensing system and an abortive infection system would look "
        "like. An interaction above 1 would mean redundancy, the second system "
        "adding nothing once the first is present.")

    add("")
    add("DOES BREADTH MATTER MORE THAN COUNT")
    add(f"  mechanistic classes defined        : {len(classes)}")
    add(f"  correlation of per-strain lysis rate with number of systems : "
        f"rho = {rho_systems.statistic:+.4f} (p = {rho_systems.pvalue:.3g})")
    add(f"  correlation of per-strain lysis rate with number of classes : "
        f"rho = {rho_classes.statistic:+.4f} (p = {rho_classes.pvalue:.3g})")
    add("")
    if len(class_tab):
        add(f"  {'class':26s} {'systems':>7s} {'carriers':>8s} "
            f"{'lysis carr':>10s} {'lysis non':>9s} {'q':>9s}")
        for r in class_tab.sort_values("mean_lysis_carriers").itertuples(index=False):
            mark = " *" if r.significant else "  "
            add(f"  {r.mechanistic_class:26s} {r.n_systems_in_class:7d} "
                f"{r.n_strains_carrying:8d} "
                f"{100 * r.mean_lysis_carriers:9.2f}% "
                f"{100 * r.mean_lysis_noncarriers:8.2f}% {r.q_value:9.2e}{mark}")
        nclass = int(class_tab["significant"].sum())
        add("")
        add(f"  Mechanistic classes associated with per-strain lysis rate at "
            f"q <= {config.FDR_ALPHA}: {nclass} of {len(class_tab)}. These are "
            "unadjusted strain-level comparisons and inherit every confound "
            "that script 04 removed, so they are descriptive only.")

    add("")
    add("CONCLUSION")
    if nsig == 0:
        add("  The null on defense burden is not hiding combination effects. "
            "Of the "
            f"{len(inter_tab)} pairs with enough strains in every carriage "
            "combination to support an interaction term, none shows one. "
            "Counting systems fails to predict lysis, and so does pairing "
            "them.")
    else:
        add(f"  {nsig} pairs show an interaction, so the burden null is not the "
            "whole story. Each would need the same clonal-background scrutiny "
            "applied to the single systems in script 09 before being treated "
            "as more than a candidate.")
    add("  Neither test reaches three-way or higher combinations. With 402 "
        "strains and this much correlation between defense repertoires, those "
        "cells are too thinly populated to estimate, and no amount of "
        "modelling recovers information the collection does not contain.")

    text = "\n".join(o)
    (paths.RESULTS / "defense_combinations_report.txt").write_text(
        text + "\n", encoding="utf-8")
    print("\n" + text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
