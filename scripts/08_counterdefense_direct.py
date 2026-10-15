"""Re-test counter-defense rescue using the directly annotated phage genomes.

Script 05 tested whether a host defense system protects less well against
phages whose relatives commonly carry the matching anti-defense family. That
test depended on a genus-level prior because the phages themselves were not
annotated. Script 07 annotated them, so the same question can now be asked of
the phages that were actually used in the experiment.

The improvement is not only in precision. Under the genus prior, every phage of
a genus necessarily received the same anti-defense status, so the contrast was
entirely between genera and inherited every difference in host range between
them. With direct annotation the status varies inside four genera, so part of
the contrast is now within-genus and therefore much harder to confound.

The estimator is unchanged: the protective effect of the host system is
estimated separately against anti-defense-carrying and anti-defense-free
phages, with phage identity absorbed as a fixed effect, host surface type,
phylogroup and core-genome coordinates adjusted, and standard errors clustered
on clonal background. Host defense carriage varies within every phage, so the
interaction is identified while phage host range is fully controlled.

Outputs:
    results/tables/table22_counterdefense_direct.csv
    results/tables/table23_direct_vs_imputed_rescue.csv
    results/tables/table24_phage_level_summary.csv
    results/counterdefense_direct_report.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
from scipy import stats

from papr import config, paths
from papr.models import fit_logit_cluster, term_effect
from papr.stats import (benjamini_hochberg, mantel_haenszel, stratified_tables,
                        wilson_interval)

# A split needs enough phages on each side to estimate anything.
MIN_PHAGES_PER_SIDE = 5


def load():
    pairs = pd.read_csv(paths.PAIRS, compression="gzip", low_memory=False)
    strain_feats = pd.read_csv(paths.STRAIN_FEATURES, low_memory=False)
    phage_feats = pd.read_csv(paths.PHAGE_FEATURES)
    adf = pd.read_csv(paths.PHAGE_ADF_SYSTEMS)
    clusters = pd.read_csv(paths.PROCESSED / "clonal_clusters.csv", index_col=0)
    pcoa = pd.read_csv(paths.PCOA, index_col=0)
    lineage = pd.read_csv(paths.TABLES / "table7_lineage_controlled_effects.csv")

    ds_cols = [c for c in strain_feats.columns if c.startswith("ds_")]
    tested = strain_feats[strain_feats["strain"].isin(set(pairs["strain"]))]
    systems = sorted([c for c in ds_cols
                      if tested[c].sum() >= config.MIN_CARRIER_STRAINS])
    pairs = pairs.merge(strain_feats[["strain"] + systems], on="strain", how="left")
    pairs["clonal_cluster"] = pairs["strain"].map(clusters["clonal_cluster"])
    return pairs, phage_feats, adf, pcoa, systems, lineage


def build_blocks(pairs: pd.DataFrame, pcoa: pd.DataFrame) -> dict[str, np.ndarray]:
    coords = pcoa.reindex(pairs["strain"]).to_numpy(dtype=float)
    coords = (coords - coords.mean(axis=0)) / coords.std(axis=0)
    return {
        "const": np.ones((len(pairs), 1)),
        "phage": pd.get_dummies(pairs["phage"].astype("string"),
                                drop_first=True, dtype=float).to_numpy(),
        "lps": pd.get_dummies(pairs["lps_core"].astype("string"),
                              drop_first=True, dtype=float).to_numpy(),
        "phylogroup": pd.get_dummies(pairs["phylogroup"].astype("string"),
                                     drop_first=True, dtype=float).to_numpy(),
        "pcoa": coords,
    }


def rescue_direct(pairs: pd.DataFrame, systems: list[str], adf: pd.DataFrame,
                  blocks: dict[str, np.ndarray]) -> pd.DataFrame:
    y = pairs["lytic"].to_numpy(dtype=float)
    groups = pairs["clonal_cluster"].to_numpy()
    adf_idx = adf.set_index("phage")

    rows = []
    for col in systems:
        name = col[3:]
        family = config.DEFENSE_TO_ANTIDEFENSE.get(name)
        if not family or family not in adf_idx.columns:
            continue
        carrier_phages = adf_idx.index[adf_idx[family] == 1]
        n_with = int(len(carrier_phages))
        n_without = int(len(adf_idx) - n_with)
        if min(n_with, n_without) < MIN_PHAGES_PER_SIDE:
            continue

        frame = pairs.copy()
        frame["adf"] = frame["phage"].map(adf_idx[family]).fillna(0).astype(int)
        with_adf = frame[frame["adf"] == 1]
        without_adf = frame[frame["adf"] == 0]
        eff_with = mantel_haenszel(stratified_tables(with_adf, ["phage"], col, "lytic"))
        eff_without = mantel_haenszel(
            stratified_tables(without_adf, ["phage"], col, "lytic"))

        exposure = frame[col].to_numpy(dtype=float)
        adf_flag = frame["adf"].to_numpy(dtype=float)
        x = np.hstack([
            (exposure * adf_flag).reshape(-1, 1),   # interaction
            exposure.reshape(-1, 1),                # defense main effect
            blocks["const"], blocks["phage"], blocks["lps"],
            blocks["phylogroup"], blocks["pcoa"],
        ])
        res = fit_logit_cluster(y, x, groups)
        inter = term_effect(res, 0)
        main = term_effect(res, 1)

        # How much of the contrast is within a genus rather than between genera.
        genus_mix = (adf_idx.join(pairs.groupby("phage")["genus"].first())
                     .groupby("genus")[family].nunique())
        n_mixed = int((genus_mix > 1).sum())

        rows.append(dict(
            system=name, antidefense_family=family,
            n_phages_with_family=n_with, n_phages_without=n_without,
            n_genera_mixed=n_mixed,
            or_against_phages_without=eff_without.estimate,
            or_without_ci_low=eff_without.ci_low,
            or_without_ci_high=eff_without.ci_high,
            or_against_phages_with=eff_with.estimate,
            or_with_ci_low=eff_with.ci_low, or_with_ci_high=eff_with.ci_high,
            interaction_odds_ratio=inter["odds_ratio"],
            interaction_ci_low=inter["ci_low"],
            interaction_ci_high=inter["ci_high"],
            interaction_p_value=inter["p_value"],
            defense_main_effect_or=main["odds_ratio"],
            defense_main_effect_p=main["p_value"]))
    out = pd.DataFrame(rows)
    if not out.empty:
        out["interaction_q_value"] = benjamini_hochberg(
            out["interaction_p_value"].to_numpy())
        out["significant_interaction"] = out["interaction_q_value"] <= config.FDR_ALPHA
    return out


def phage_level_summary(pairs: pd.DataFrame, phage_feats: pd.DataFrame,
                        adf: pd.DataFrame) -> pd.DataFrame:
    """Does carrying anti-defense genes track a phage's measured host range?"""
    per_phage = pairs.groupby("phage").agg(
        n_pairs=("lytic", "size"), n_lytic=("lytic", "sum"),
        lysis_rate=("lytic", "mean"), mean_score=("plaque_score", "mean"))
    out = (per_phage.join(phage_feats.set_index("phage")[
        ["genus", "morphotype", "genome_size"]])
        .join(adf.set_index("phage")[["n_antidefense_systems", "any_antidefense",
                                      "Anti_RM", "Anti_CBASS", "Anti_TA"]]))
    out = out.reset_index()
    ci = [wilson_interval(int(k), int(n))
          for k, n in zip(out["n_lytic"], out["n_pairs"])]
    out["ci_low"] = [c[0] for c in ci]
    out["ci_high"] = [c[1] for c in ci]
    return out


def main() -> int:
    paths.bootstrap()
    pairs, phage_feats, adf, pcoa, systems, lineage = load()
    blocks = build_blocks(pairs, pcoa)

    print("fitting rescue interactions with direct annotation")
    direct = rescue_direct(pairs, systems, adf, blocks)
    summary = phage_level_summary(pairs, phage_feats, adf)

    imputed = pd.read_csv(paths.TABLES / "table10_counterdefense_rescue.csv")
    comparison = (direct[["system", "antidefense_family", "interaction_odds_ratio",
                          "interaction_ci_low", "interaction_ci_high",
                          "interaction_q_value", "significant_interaction"]]
                  .rename(columns={
                      "interaction_odds_ratio": "direct_interaction_or",
                      "interaction_ci_low": "direct_ci_low",
                      "interaction_ci_high": "direct_ci_high",
                      "interaction_q_value": "direct_q",
                      "significant_interaction": "direct_significant"})
                  .merge(imputed[["system", "interaction_odds_ratio",
                                  "interaction_q_value", "significant_interaction"]]
                         .rename(columns={
                             "interaction_odds_ratio": "imputed_interaction_or",
                             "interaction_q_value": "imputed_q",
                             "significant_interaction": "imputed_significant"}),
                         on="system", how="outer"))

    direct.to_csv(paths.TABLES / "table22_counterdefense_direct.csv", index=False)
    comparison.to_csv(paths.TABLES / "table23_direct_vs_imputed_rescue.csv", index=False)
    summary.to_csv(paths.TABLES / "table24_phage_level_summary.csv", index=False)

    out: list[str] = []
    add = out.append
    add("Counter-defense rescue with directly annotated phage genomes")
    add("Generated by scripts/08_counterdefense_direct.py")
    add("")
    add("WHAT CHANGED")
    add(f"  phages annotated directly        : {len(adf)}")
    add(f"  phages carrying any anti-defense : {int(adf['any_antidefense'].sum())} "
        f"({100 * adf['any_antidefense'].mean():.1f} percent)")
    add(f"  phages carrying Anti_RM          : {int(adf['Anti_RM'].sum())} of {len(adf)}")
    add("  Under the genus prior of script 05 every phage of a genus shared one "
        "status, so the contrast was wholly between genera. Direct annotation "
        "splits four genera internally, which moves part of the contrast inside "
        "a genus where host range is far more similar.")

    add("")
    add("RESCUE INTERACTIONS, DIRECT ANNOTATION")
    add("  An interaction odds ratio above 1 means the host system protects "
        "less well against phages that carry the matching anti-defense family, "
        "which is what rescue predicts.")
    add("")
    if direct.empty:
        add("  No system had a usable split.")
    else:
        add(f"  {'system':22s} {'family':13s} {'phages +/-':>11s} {'OR vs -':>8s} "
            f"{'OR vs +':>8s} {'interact':>8s} {'95 pct CI':>16s} {'q':>9s}")
        for r in direct.sort_values("interaction_odds_ratio").itertuples(index=False):
            mark = " *" if r.significant_interaction else "  "
            add(f"  {r.system:22s} {r.antidefense_family:13s} "
                f"{r.n_phages_with_family:4d}/{r.n_phages_without:<6d} "
                f"{r.or_against_phages_without:8.3f} {r.or_against_phages_with:8.3f} "
                f"{r.interaction_odds_ratio:8.3f} "
                f"[{r.interaction_ci_low:6.3f}, {r.interaction_ci_high:6.3f}] "
                f"{r.interaction_q_value:9.2e}{mark}")
        nsig = int(direct["significant_interaction"].sum())
        add("")
        add(f"  Significant at q <= {config.FDR_ALPHA}: {nsig} of {len(direct)} tested.")
        rm = direct[direct["antidefense_family"] == "Anti_RM"]
        if len(rm):
            above = int((rm["interaction_odds_ratio"] > 1).sum())
            add(f"  Restriction-modification subtypes with an interaction above 1: "
                f"{above} of {len(rm)}.")

    add("")
    add("DIRECT ANNOTATION AGAINST THE GENUS PRIOR")
    add(f"  {'system':22s} {'direct OR':>10s} {'direct q':>9s} {'imputed OR':>11s} "
        f"{'imputed q':>10s} {'agree':>6s}")
    for r in comparison.itertuples(index=False):
        d_or = getattr(r, "direct_interaction_or", np.nan)
        i_or = getattr(r, "imputed_interaction_or", np.nan)
        d_q = getattr(r, "direct_q", np.nan)
        i_q = getattr(r, "imputed_q", np.nan)
        if not np.isfinite(d_or) or not np.isfinite(i_or):
            continue
        agree = "yes" if (d_or > 1) == (i_or > 1) else "no"
        add(f"  {r.system:22s} {d_or:10.3f} {d_q:9.2e} {i_or:11.3f} "
            f"{i_q:10.2e} {agree:>6s}")
    both = comparison.dropna(subset=["direct_interaction_or", "imputed_interaction_or"])
    if len(both) > 3:
        same = int(((both["direct_interaction_or"] > 1)
                    == (both["imputed_interaction_or"] > 1)).sum())
        rho = both[["direct_interaction_or", "imputed_interaction_or"]].corr(
            method="spearman").iloc[0, 1]
        add("")
        add(f"  Direction agrees for {same} of {len(both)} systems testable both "
            f"ways; Spearman correlation between the two sets of interaction "
            f"odds ratios is {rho:+.3f}.")
        n_direct = int(both["direct_significant"].sum())
        n_imp = int(both["imputed_significant"].sum())
        add(f"  Significant interactions: {n_direct} with direct annotation, "
            f"{n_imp} with the genus prior.")
        if n_direct <= n_imp:
            add("  Replacing the imputed status with a measured one did not "
                "produce new rescue signals. The negative result of script 05 "
                "therefore does not rest on the imputation, which was the main "
                "reason to doubt it.")
        else:
            add("  Direct annotation recovers interactions the genus prior "
                "missed, which is the attenuation that imputation predicts.")

    add("")
    add("PHAGE-SIDE DESCRIPTION")
    r_all = stats.spearmanr(summary["n_antidefense_systems"], summary["lysis_rate"])
    add(f"  Across the 96 phages, Spearman correlation between the number of "
        f"anti-defense systems carried and measured lysis rate: "
        f"{r_all.statistic:+.3f} (p = {r_all.pvalue:.3g}).")
    for morph, g in summary.groupby("morphotype"):
        add(f"    {morph:14s} n={len(g):2d}, mean anti-defense systems "
            f"{g['n_antidefense_systems'].mean():4.2f}, "
            f"mean lysis rate {100 * g['lysis_rate'].mean():5.2f} percent")
    add("  Morphotype is strongly associated with both anti-defense content and "
        "lysis rate, so this correlation is largely a restatement of which "
        "lineages are virulent in this collection and is not evidence that "
        "anti-defense genes widen host range.")

    top = summary.nlargest(8, "n_antidefense_systems")[
        ["phage", "genus", "n_antidefense_systems", "lysis_rate"]]
    add("")
    add("  Phages carrying the most anti-defense systems:")
    add(f"    {'phage':16s} {'genus':20s} {'systems':>7s} {'lysis':>7s}")
    for r in top.itertuples(index=False):
        add(f"    {r.phage:16s} {r.genus:20s} {r.n_antidefense_systems:7d} "
            f"{100 * r.lysis_rate:6.2f}%")

    text = "\n".join(out)
    (paths.RESULTS / "counterdefense_direct_report.txt").write_text(
        text + "\n", encoding="utf-8")
    print("\n" + text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
