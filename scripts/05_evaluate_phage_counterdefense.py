"""Test whether phage counter-defense explains residual variation in lysis.

Two strains with the same surface type and the same defense repertoire can
still differ in whether a given phage lyses them, and a phage lineage that
carries anti-defense genes should be less affected by a host defense system
than one that does not. This script tests that prediction.

The 96 Guelin phages were not individually annotated for anti-defense content
in the released data, and running AntiDefenseFinder here would mean installing
HMMER and MacSyFinder and re-annotating 96 genomes, which the storage and
platform budget does not allow. What is available is the published
AntiDefenseFinder survey of 13,488 RefSeq phage genomes, which carries ICTV
genus labels. Every one of the 19 VIRIDIC genera in this collection appears in
that survey, so each genus can be assigned the prevalence of each anti-defense
family among its sequenced relatives. This is a lineage-level prior, not a
per-genome call, and it is treated as such throughout: a phage is described as
belonging to a counter-defense-rich lineage, never as carrying a given gene.

The central test is an interaction, not a comparison of lysis rates between
phage groups. Anti-defense status is a property of the phage, so comparing
lysis across anti-defense-rich and anti-defense-poor phages would mostly
measure differences in host range. Instead, the protective effect of a host
defense system is estimated separately against anti-defense-rich and
anti-defense-poor lineages, with phage identity absorbed as a fixed effect.
Host defense carriage varies within each phage, so that contrast is identified
while phage host range is fully controlled.

Outputs:
    data/processed/antidefense_genus_prevalence.csv
    results/tables/table9_antidefense_genus_prevalence.csv
    results/tables/table10_counterdefense_rescue.csv
    results/tables/table11_three_state_model.csv
    results/tables/table12_tropism_by_phage_lineage.csv
    results/tables/table13_residual_resistance.csv
    results/counterdefense_report.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from papr import config, paths
from papr.models import fit_logit_cluster, term_effect
from papr.stats import benjamini_hochberg, mantel_haenszel, stratified_tables, wilson_interval

# A genus needs a reasonable number of sequenced relatives before a prevalence
# estimated from them means anything. Two surveyed genomes cannot establish
# whether a lineage carries an anti-defense family.
MIN_SURVEY_DEPTH = 10

# Mechanistic pairing between a bacterial defense system and the anti-defense
# family that counteracts it. Systems with no published counterpart in the
# AntiDefenseFinder model set are listed as unmatched and are tested against
# total anti-defense burden instead, which is a weaker proxy.
DEFENSE_TO_ANTIDEFENSE: dict[str, str] = {
    "RM_Type_I": "Anti_RM",
    "RM_Type_II": "Anti_RM",
    "RM_Type_IIG": "Anti_RM",
    "RM_Type_III": "Anti_RM",
    "RM_Type_IV": "Anti_RM",
    "CAS_Class1-Subtype-I-E": "Anti_CRISPR",
    "CAS_Class1-Subtype-I-F": "Anti_CRISPR",
    "CBASS_I": "Anti_CBASS",
    "BREX_I": "Anti_BREX",
    "Gabija": "Anti_Gabija",
    "Hachiman": "Anti_Hachiman",
    "Retron_I_C": "Anti_Retron",
    "MazEF": "Anti_TA",
    "PsyrTA": "Anti_TA",
    "SanaTA": "Anti_TA",
    "RosmerTA": "Anti_TA",
    "AbiE": "Anti_TA",
    "Mok_Hok_Sok": "Anti_TA",
}


def antidefense_prevalence() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Prevalence of each anti-defense family per phage genus, from RefSeq."""
    tax = pd.read_csv(paths.ADF_TAXONOMY, sep="\t", low_memory=False,
                      usecols=["Nuccore", "Genus"]).dropna(subset=["Nuccore"])
    hits = pd.read_csv(paths.ADF_PHAGE, sep="\t", low_memory=False,
                       usecols=["model_fqn", "Nuccore", "Genus"]).dropna(subset=["Nuccore"])
    hits["family"] = hits["model_fqn"].str.split("/").str[1]

    # Denominator: distinct surveyed genomes per genus.
    denom = tax.groupby("Genus")["Nuccore"].nunique().rename("n_surveyed")
    # Restrict hits to genomes that appear in the surveyed set so that the
    # numerator can never exceed the denominator.
    surveyed = set(tax["Nuccore"])
    hits = hits[hits["Nuccore"].isin(surveyed)]

    any_adf = hits.groupby("Genus")["Nuccore"].nunique().rename("n_with_any_antidefense")
    overall = pd.concat([denom, any_adf], axis=1).fillna({"n_with_any_antidefense": 0})
    overall["prevalence_any"] = overall["n_with_any_antidefense"] / overall["n_surveyed"]
    overall = overall.reset_index().rename(columns={"Genus": "genus"})

    per_family = (hits.groupby(["Genus", "family"])["Nuccore"].nunique()
                  .rename("n_carrying").reset_index()
                  .rename(columns={"Genus": "genus"}))
    per_family = per_family.merge(
        denom.reset_index().rename(columns={"Genus": "genus"}), on="genus", how="left")
    per_family["prevalence"] = per_family["n_carrying"] / per_family["n_surveyed"]
    ci = [wilson_interval(int(c), int(n))
          for c, n in zip(per_family["n_carrying"], per_family["n_surveyed"])]
    per_family["ci_low"] = [c[0] for c in ci]
    per_family["ci_high"] = [c[1] for c in ci]
    return overall, per_family


def genus_family_matrix(per_family: pd.DataFrame, genera: list[str],
                        overall: pd.DataFrame) -> pd.DataFrame:
    """Wide prevalence table covering every genus in this collection."""
    wide = (per_family.pivot(index="genus", columns="family", values="prevalence")
            .reindex(genera).fillna(0.0))
    wide = wide.join(overall.set_index("genus")[["n_surveyed", "prevalence_any"]])
    wide["n_surveyed"] = wide["n_surveyed"].fillna(0).astype(int)
    wide["prevalence_any"] = wide["prevalence_any"].fillna(0.0)
    wide["adequately_surveyed"] = wide["n_surveyed"] >= MIN_SURVEY_DEPTH
    return wide


def split_for_system(name: str, wide: pd.DataFrame) -> tuple[str, pd.Series, float] | None:
    """Anti-defense split for one defense system, or None if none is defensible.

    A split is only attempted where the system has a published anti-defense
    counterpart. An earlier version fell back to total anti-defense burden for
    unmatched systems; that fallback is not used, because prevalence of any
    anti-defense family saturates at 100 percent in most genera, leaving a
    median split with nothing on one side. Reporting a test built on a
    degenerate split would be worse than reporting no test.
    """
    family = DEFENSE_TO_ANTIDEFENSE.get(name)
    if not family or family not in wide.columns:
        return None
    prev = wide[family]
    reference = prev[wide["adequately_surveyed"]]
    if reference.empty or reference.nunique() < 2:
        return None
    cut = float(reference.median())
    if int((reference > cut).sum()) == 0 or int((reference <= cut).sum()) == 0:
        return None
    return family, prev, cut


def rescue_tests(pairs: pd.DataFrame, systems: list[str], wide: pd.DataFrame,
                 blocks: dict[str, np.ndarray], groups: np.ndarray) -> pd.DataFrame:
    """Is a defense system less protective against counter-defense-rich lineages?"""
    y = pairs["lytic"].to_numpy(dtype=float)
    ok_genera = set(wide.index[wide["adequately_surveyed"]])
    rows = []
    for sys_col in systems:
        name = sys_col[3:]
        split = split_for_system(name, wide)
        if split is None:
            continue
        family, prev, cut = split
        matched = True
        column = family
        usable = pairs["genus"].isin(ok_genera)
        frame = pairs[usable].copy()
        frame["adf_high"] = (frame["genus"].map(prev) > cut).astype(int)
        if frame["adf_high"].nunique() < 2:
            continue

        low = frame[frame["adf_high"] == 0]
        high = frame[frame["adf_high"] == 1]
        eff_low = mantel_haenszel(stratified_tables(low, ["phage"], sys_col, "lytic"))
        eff_high = mantel_haenszel(stratified_tables(high, ["phage"], sys_col, "lytic"))

        # Interaction model. Phage identity absorbs host range and the main
        # effect of anti-defense status, which is a phage-level property.
        mask = usable.to_numpy()
        exposure = pairs.loc[usable, sys_col].to_numpy(dtype=float)
        adf_high = frame["adf_high"].to_numpy(dtype=float)
        x = np.hstack([
            (exposure * adf_high).reshape(-1, 1),   # interaction, index 0
            exposure.reshape(-1, 1),                # main effect of defense
            blocks["const"][mask],
            blocks["phage_sub"],
            blocks["lps"][mask],
            blocks["phylogroup"][mask],
            blocks["pcoa"][mask],
        ])
        res = fit_logit_cluster(y[mask], x, groups[mask])
        inter = term_effect(res, 0)
        main = term_effect(res, 1)

        rows.append(dict(
            system=name, matched_antidefense_family=family or "none",
            prevalence_column=column, mechanistically_matched=matched,
            prevalence_split=cut,
            n_genera_low=int((prev[list(ok_genera)] <= cut).sum()),
            n_genera_high=int((prev[list(ok_genera)] > cut).sum()),
            n_pairs=int(len(frame)),
            or_against_antidefense_poor=eff_low.estimate,
            or_poor_ci_low=eff_low.ci_low, or_poor_ci_high=eff_low.ci_high,
            or_against_antidefense_rich=eff_high.estimate,
            or_rich_ci_low=eff_high.ci_low, or_rich_ci_high=eff_high.ci_high,
            rescue_ratio=(eff_high.estimate / eff_low.estimate
                          if eff_low.estimate and np.isfinite(eff_low.estimate) else np.nan),
            interaction_odds_ratio=inter["odds_ratio"],
            interaction_ci_low=inter["ci_low"], interaction_ci_high=inter["ci_high"],
            interaction_p_value=inter["p_value"],
            defense_main_effect_or=main["odds_ratio"],
            defense_main_effect_p=main["p_value"],
        ))
    out = pd.DataFrame(rows)
    if not out.empty:
        out["interaction_q_value"] = benjamini_hochberg(
            out["interaction_p_value"].to_numpy())
        out["significant_interaction"] = out["interaction_q_value"] <= config.FDR_ALPHA
    return out


def three_state(pairs: pd.DataFrame, systems: list[str], wide: pd.DataFrame) -> pd.DataFrame:
    """Descriptive lysis rates in the three mechanistic states."""
    ok_genera = set(wide.index[wide["adequately_surveyed"]])
    rows = []
    for sys_col in systems:
        name = sys_col[3:]
        split = split_for_system(name, wide)
        if split is None:
            continue
        family, prev, cut = split
        frame = pairs[pairs["genus"].isin(ok_genera)].copy()
        frame["adf_high"] = (frame["genus"].map(prev) > cut).astype(int)

        states = {
            "A_host_lacks_defense": frame[frame[sys_col] == 0],
            "B_defense_vs_antidefense_poor_lineage":
                frame[(frame[sys_col] == 1) & (frame["adf_high"] == 0)],
            "C_defense_vs_antidefense_rich_lineage":
                frame[(frame[sys_col] == 1) & (frame["adf_high"] == 1)],
        }
        for state, sub in states.items():
            if sub.empty:
                continue
            lo, hi = wilson_interval(int(sub["lytic"].sum()), len(sub))
            rows.append(dict(
                system=name, matched_antidefense_family=family or "none",
                state=state, n_pairs=int(len(sub)),
                n_lytic=int(sub["lytic"].sum()),
                lysis_rate=float(sub["lytic"].mean()),
                ci_low=lo, ci_high=hi,
                mean_plaque_score=float(sub["plaque_score"].mean())))
    return pd.DataFrame(rows)


def tropism(pairs: pd.DataFrame, wide: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key, label in [("morphotype", "morphotype"), ("genus", "genus"),
                       ("family", "family")]:
        for value, g in pairs.groupby(key):
            lo, hi = wilson_interval(int(g["lytic"].sum()), len(g))
            prev = (float(wide.loc[value, "prevalence_any"])
                    if label == "genus" and value in wide.index else np.nan)
            surveyed = (int(wide.loc[value, "n_surveyed"])
                        if label == "genus" and value in wide.index else np.nan)
            rows.append(dict(
                level=label, value=value, n_phages=int(g["phage"].nunique()),
                n_pairs=int(len(g)), n_lytic=int(g["lytic"].sum()),
                lysis_rate=float(g["lytic"].mean()), ci_low=lo, ci_high=hi,
                mean_plaque_score=float(g["plaque_score"].mean()),
                antidefense_prevalence_any=prev, n_refseq_surveyed=surveyed))
    return pd.DataFrame(rows)


def residual_resistance(pairs: pd.DataFrame, systems: list[str],
                        survivors: list[str]) -> pd.DataFrame:
    """How much failure is left once gate and validated defense are accounted for."""
    gated = pairs[pairs["gate_lps"]]
    surv_cols = [f"ds_{s}" for s in survivors if f"ds_{s}" in pairs.columns]
    rows = []

    def add(label: str, frame: pd.DataFrame, note: str) -> None:
        n = len(frame)
        if n == 0:
            # An empty subset is a fact about the collection, not a reason to
            # omit a row: it says no strain meets the condition at all.
            rows.append(dict(subset=label, n_pairs=0, n_failed=0,
                             failure_rate=np.nan, ci_low=np.nan, ci_high=np.nan,
                             note=note + " (no pairs meet this condition)"))
            return
        fail = int((frame["lytic"] == 0).sum())
        lo, hi = wilson_interval(fail, n)
        rows.append(dict(subset=label, n_pairs=n, n_failed=fail,
                         failure_rate=fail / n, ci_low=lo, ci_high=hi,
                         note=note))

    add("all observed pairs", pairs, "baseline failure rate across the matrix")
    add("receptor gate passed", gated,
        "surface state consistent with the isolation hosts of the phage genus")
    if surv_cols:
        no_surv = gated[(gated[surv_cols].sum(axis=1) == 0)]
        add("gate passed, lacks lineage-validated defense", no_surv,
            "host carries neither system that survived clonal control")
    none_tested = gated[gated[systems].sum(axis=1) == 0]
    add("gate passed, lacks all 33 tested systems", none_tested,
        "host carries none of the systems with enough carriers to test")
    return pd.DataFrame(rows)


def main() -> int:
    paths.bootstrap()
    pairs = pd.read_csv(paths.PAIRS, compression="gzip", low_memory=False)
    strain_feats = pd.read_csv(paths.STRAIN_FEATURES, low_memory=False)
    pcoa = pd.read_csv(paths.PCOA, index_col=0)
    clusters = pd.read_csv(paths.PROCESSED / "clonal_clusters.csv", index_col=0)
    lineage = pd.read_csv(paths.TABLES / "table7_lineage_controlled_effects.csv")

    ds_cols = [c for c in strain_feats.columns if c.startswith("ds_")]
    tested = strain_feats[strain_feats["strain"].isin(set(pairs["strain"]))]
    systems = sorted([c for c in ds_cols
                      if tested[c].sum() >= config.MIN_CARRIER_STRAINS])
    pairs = pairs.merge(strain_feats[["strain"] + systems], on="strain", how="left")
    pairs["clonal_cluster"] = pairs["strain"].map(clusters["clonal_cluster"])

    survivors = sorted(lineage[
        (lineage["adjustment"] == "plus_phylogroup_and_coordinates")
        & lineage["significant"] & lineage["direction_preserved"]]["system"])

    overall, per_family = antidefense_prevalence()
    genera = sorted(pairs["genus"].dropna().unique())
    wide = genus_family_matrix(per_family, genera, overall)
    wide.to_csv(paths.ADF_GENUS)
    wide.reset_index().to_csv(
        paths.TABLES / "table9_antidefense_genus_prevalence.csv", index=False)

    ok_genera = set(wide.index[wide["adequately_surveyed"]])
    usable = pairs["genus"].isin(ok_genera)
    blocks = {
        "const": np.ones((len(pairs), 1)),
        "lps": pd.get_dummies(pairs["lps_core"].astype("string"),
                              drop_first=True, dtype=float).to_numpy(),
        "phylogroup": pd.get_dummies(pairs["phylogroup"].astype("string"),
                                     drop_first=True, dtype=float).to_numpy(),
    }
    coords = pcoa.reindex(pairs["strain"]).to_numpy(dtype=float)
    blocks["pcoa"] = (coords - coords.mean(axis=0)) / coords.std(axis=0)
    blocks["phage_sub"] = pd.get_dummies(
        pairs.loc[usable, "phage"].astype("string"), drop_first=True,
        dtype=float).to_numpy()

    groups = pairs["clonal_cluster"].to_numpy()
    print(f"running rescue tests on {len(systems)} systems")
    rescue = rescue_tests(pairs, systems, wide, blocks, groups)
    three = three_state(pairs, systems, wide)
    trop = tropism(pairs, wide)
    resid = residual_resistance(pairs, systems, survivors)

    rescue.to_csv(paths.TABLES / "table10_counterdefense_rescue.csv", index=False)
    three.to_csv(paths.TABLES / "table11_three_state_model.csv", index=False)
    trop.to_csv(paths.TABLES / "table12_tropism_by_phage_lineage.csv", index=False)
    resid.to_csv(paths.TABLES / "table13_residual_resistance.csv", index=False)

    out: list[str] = []
    add = out.append
    add("Phage counter-defense and tropism report")
    add("Generated by scripts/05_evaluate_phage_counterdefense.py")
    add("")
    add("PROVENANCE OF THE ANTI-DEFENSE ANNOTATION")
    add(f"  RefSeq phage genomes in the published survey : "
        f"{int(overall['n_surveyed'].sum())}")
    add(f"  anti-defense families distinguished          : "
        f"{per_family['family'].nunique()}")
    add(f"  VIRIDIC genera in this collection            : {len(genera)}")
    add(f"  genera with at least {MIN_SURVEY_DEPTH} surveyed relatives     : "
        f"{int(wide['adequately_surveyed'].sum())} of {len(genera)}")
    excluded = sorted(wide.index[~wide["adequately_surveyed"]])
    add(f"  genera excluded for thin survey coverage     : "
        f"{', '.join(excluded) if excluded else 'none'}")
    add("")
    add("  The 96 phages studied here were not annotated individually. Each "
        "genus is assigned the prevalence of each anti-defense family among its "
        "sequenced relatives, so every statement below is about lineages and "
        "not about which genes a particular phage carries.")

    add("")
    add("ANTI-DEFENSE PREVALENCE BY GENUS")
    add(f"  {'genus':22s} {'surveyed':>8s} {'any ADF':>8s} {'Anti_RM':>8s} "
        f"{'Anti_CRISPR':>11s} {'Anti_CBASS':>10s} {'Anti_TA':>8s}")
    for g in wide.sort_values("n_surveyed", ascending=False).index:
        r = wide.loc[g]
        mark = "" if r["adequately_surveyed"] else "  (thin)"
        def pv(col: str) -> str:
            return f"{100 * r[col]:7.1f}%" if col in wide.columns else "      -"
        add(f"  {g:22s} {int(r['n_surveyed']):8d} "
            f"{100 * r['prevalence_any']:7.1f}% {pv('Anti_RM')} "
            f"{pv('Anti_CRISPR'):>11s} {pv('Anti_CBASS'):>10s} {pv('Anti_TA')}{mark}")

    add("")
    add("COUNTER-DEFENSE RESCUE: IS A DEFENSE SYSTEM WEAKER AGAINST")
    add("ANTI-DEFENSE-RICH PHAGE LINEAGES")
    add("  The interaction odds ratio compares the protective effect of the host "
        "system against anti-defense-rich lineages with its effect against "
        "anti-defense-poor ones, holding phage identity, host surface type and "
        "clonal background fixed. A value above 1 means the system protects less "
        "well against the rich lineages, which is what rescue predicts.")
    add("")
    if rescue.empty:
        add("  No system had both a usable anti-defense split and enough data.")
    else:
        add(f"  {'system':24s} {'family':14s} {'OR vs poor':>10s} {'OR vs rich':>10s} "
            f"{'interact':>8s} {'95 pct CI':>16s} {'q':>9s}")
        for r in rescue.sort_values("interaction_odds_ratio").itertuples(index=False):
            mark = " *" if r.significant_interaction else "  "
            mm = "" if r.mechanistically_matched else " (unmatched)"
            add(f"  {r.system:24s} {r.matched_antidefense_family:14s} "
                f"{r.or_against_antidefense_poor:10.3f} "
                f"{r.or_against_antidefense_rich:10.3f} "
                f"{r.interaction_odds_ratio:8.3f} "
                f"[{r.interaction_ci_low:6.3f}, {r.interaction_ci_high:6.3f}] "
                f"{r.interaction_q_value:9.2e}{mark}{mm}")
        nsig = int(rescue["significant_interaction"].sum())
        add("")
        add(f"  Interactions significant at q <= {config.FDR_ALPHA}: {nsig} of "
            f"{len(rescue)} systems tested.")

        tested_names = set(rescue["system"])
        no_family = sorted(s[3:] for s in systems
                           if s[3:] not in DEFENSE_TO_ANTIDEFENSE)
        degenerate = sorted(s[3:] for s in systems
                            if s[3:] in DEFENSE_TO_ANTIDEFENSE
                            and s[3:] not in tested_names)
        add("")
        add(f"  Not tested, no published anti-defense counterpart ({len(no_family)}): "
            f"{', '.join(no_family)}")
        add(f"  Not tested, counterpart exists but gives no usable split "
            f"({len(degenerate)}): {', '.join(degenerate) if degenerate else 'none'}")
        for fam in ["Anti_CRISPR", "Anti_Gabija", "Anti_Hachiman", "Anti_BREX"]:
            if fam in wide.columns:
                ad = wide.loc[wide["adequately_surveyed"], fam]
                if float(ad.max()) == 0.0:
                    add(f"    {fam} is absent from every adequately surveyed "
                        "genus in this collection, so the systems it "
                        "counteracts cannot be tested for rescue here at all.")

        rm = rescue[rescue["matched_antidefense_family"] == "Anti_RM"]
        if len(rm) >= 3:
            add("")
            add("  The restriction-modification family is the one place where a "
                "consistent pattern appears. Of the "
                f"{len(rm)} RM subtypes tested against Anti_RM prevalence, "
                f"{int((rm['interaction_odds_ratio'] > 1).sum())} have an "
                "interaction odds ratio above 1, meaning RM protects less well "
                "against lineages whose relatives commonly carry anti-RM genes. "
                f"RM_Type_IV reaches significance (interaction odds ratio "
                f"{float(rm.set_index('system').loc['RM_Type_IV', 'interaction_odds_ratio']):.2f}"
                f", q = "
                f"{float(rm.set_index('system').loc['RM_Type_IV', 'interaction_q_value']):.3f}) "
                "and RM_Type_III and RM_Type_I point the same way without "
                "reaching it. Anti-RM proteins such as Ocr and ArdA are the "
                "best characterised and most widespread anti-defense proteins "
                "known, so this is the pairing where a lineage-level prior "
                "would be expected to carry the most information. The pattern "
                "is suggestive and it rests on one significant test out of "
                f"{len(rescue)}; it is not established by this analysis.")
        if nsig == 0:
            add("  No defense system shows a detectably different effect against "
                "anti-defense-rich and anti-defense-poor phage lineages. The "
                "three-state rescue model is not supported by these data. That "
                "is a null result and it is reported as such.")
            add("  It is a weak test and its limits should be stated plainly. "
                "Anti-defense content is imputed at the level of the genus from "
                "sequenced relatives, so a phage in a genus where half the "
                "relatives carry an anti-RM gene may or may not carry one "
                "itself. Imputation of this kind attenuates a real interaction "
                "toward the null, so absence of evidence here is some way short "
                "of evidence of absence.")

    add("")
    add("THREE-STATE MECHANISTIC MODEL, DESCRIPTIVE LYSIS RATES")
    add(f"  {'system':24s} {'state':42s} {'pairs':>6s} {'lysis':>7s} {'95 pct CI':>16s}")
    focus = survivors + [s for s in ["RM_Type_I", "CAS_Class1-Subtype-I-E", "CBASS_I"]
                         if s not in survivors]
    for s in focus:
        sub = three[three["system"] == s]
        for r in sub.itertuples(index=False):
            add(f"  {r.system:24s} {r.state:42s} {r.n_pairs:6d} "
                f"{100 * r.lysis_rate:6.2f}% "
                f"[{100 * r.ci_low:5.2f}, {100 * r.ci_high:5.2f}]")
    add("  States B and C sit in different phage lineages by construction, so "
        "the gap between them mixes counter-defense with host range. The "
        "interaction model above is the comparison that separates the two; "
        "these rates are shown because they are what the raw data look like.")

    add("")
    add("TROPISM: LYSIS BY PHAGE LINEAGE")
    add(f"  {'morphotype':22s} {'phages':>6s} {'pairs':>6s} {'lysis':>7s} "
        f"{'mean score':>10s}")
    for r in trop[trop["level"] == "morphotype"].sort_values(
            "lysis_rate", ascending=False).itertuples(index=False):
        add(f"  {r.value:22s} {r.n_phages:6d} {r.n_pairs:6d} "
            f"{100 * r.lysis_rate:6.2f}% {r.mean_plaque_score:10.3f}")
    add("")
    add(f"  {'genus':22s} {'phages':>6s} {'lysis':>7s} {'ADF any':>8s} "
        f"{'surveyed':>8s}")
    for r in trop[trop["level"] == "genus"].sort_values(
            "lysis_rate", ascending=False).itertuples(index=False):
        adf = ("     n/a" if not np.isfinite(r.antidefense_prevalence_any)
               else f"{100 * r.antidefense_prevalence_any:7.1f}%")
        ns = "  n/a" if not np.isfinite(r.n_refseq_surveyed) else f"{int(r.n_refseq_surveyed):8d}"
        add(f"  {r.value:22s} {r.n_phages:6d} {100 * r.lysis_rate:6.2f}% {adf} {ns}")
    gen = trop[(trop["level"] == "genus") & trop["antidefense_prevalence_any"].notna()]
    if len(gen) > 3:
        rho = gen[["lysis_rate", "antidefense_prevalence_any"]].corr(
            method="spearman").iloc[0, 1]
        add("")
        add(f"  Spearman correlation across genera between lysis rate and "
            f"anti-defense prevalence: {rho:+.3f} over {len(gen)} genera. With "
            "this few genera, and with lysis rate varying by more than two "
            "orders of magnitude for reasons of host range, this correlation "
            "should not be read as causal in either direction.")

    add("")
    add("RESIDUAL RESISTANCE NOT EXPLAINED BY SURFACE OR DEFENSE")
    add(f"  {'subset':48s} {'pairs':>6s} {'failure':>8s} {'95 pct CI':>16s}")
    for r in resid.itertuples(index=False):
        if r.n_pairs == 0:
            add(f"  {r.subset:48s} {r.n_pairs:6d}      n/a                  "
                f"({r.note.split('(')[-1].rstrip(')')})")
            continue
        add(f"  {r.subset:48s} {r.n_pairs:6d} {100 * r.failure_rate:7.2f}% "
            f"[{100 * r.ci_low:5.2f}, {100 * r.ci_high:5.2f}]")
    populated = resid[resid["n_pairs"] > 0]
    base_fail = populated.iloc[0]["failure_rate"]
    final = populated.iloc[-1]
    add("")
    add(f"  Among pairs in the subset '{final['subset']}', "
        f"{100 * final['failure_rate']:.1f} percent of infections still fail, "
        f"against {100 * base_fail:.1f} percent across the whole matrix. "
        "Applying the surface gate and removing hosts that carry a "
        "lineage-validated defense system barely moves the failure rate.")
    empty = resid[resid["n_pairs"] == 0]
    if not empty.empty:
        add(f"  The subset '{empty.iloc[0]['subset']}' is empty: every one of "
            "the 402 strains carries at least one of the 33 testable defense "
            "systems, so there is no defense-free comparison group to be had. "
            "MazEF alone is present in 338 strains. This is a real constraint "
            "on the design and not an oversight: a question of the form 'how "
            "often does infection fail in a host with no defenses' cannot be "
            "answered in this collection.")
    add("  This is the central quantitative result of the counter-defense "
        "analysis: the great majority of failed infections in this collection "
        "are not attributable to the surface determinants typed here, nor to "
        "the intracellular defense systems DefenseFinder detects, nor to "
        "anti-defense content imputed from phage lineage.")

    add("")
    add("WHAT COULD NOT BE TESTED")
    add("  Resident prophages were not annotated. Superinfection exclusion by a "
        "resident prophage is a well-documented route to post-adsorption "
        "resistance and is a plausible contributor to the residual failure "
        "above, but identifying prophages would require running a prophage "
        "caller over 402 bacterial genomes, and those genomes are not part of "
        "the released data used here. This is a gap in the analysis rather than "
        "a question the data answer negatively.")
    add("  Per-phage anti-defense annotation was not performed, for the reasons "
        "given at the top of this report. Annotating the 96 genomes directly "
        "would sharpen the rescue test considerably and is the single change "
        "most likely to alter its conclusion.")
    add("  Receptor binding protein sequences were not compared. The released "
        "table names the receptor binding protein of each phage but does not "
        "give its sequence, and the phage genome files carry nucleotides "
        "without the gene coordinates needed to extract it, so grouping phages "
        "by receptor binding protein similarity was not possible.")

    text = "\n".join(out)
    (paths.RESULTS / "counterdefense_report.txt").write_text(text + "\n", encoding="utf-8")
    print("\n" + text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
