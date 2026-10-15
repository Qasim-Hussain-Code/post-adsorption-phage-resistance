"""Stress the two surviving associations and quantify what the nulls can exclude.

Script 04 controlled for clonal background by putting phylogroup and
core-genome coordinates in a regression. That is a modelling assumption: it
works only if the covariates capture the lineage structure that matters. This
script attacks the same question from directions that do not rely on that
assumption, and reports what each one can and cannot establish.

  within-cluster       Mantel-Haenszel restricted to strata that contain both
                       carriers and non-carriers of the same near-identical
                       core genome. Between-lineage confounding is removed by
                       construction rather than by adjustment.
  permutation          carrier labels shuffled among strains inside a
                       phylogroup, and inside a Sequence Type, preserving the
                       lineage structure and the number of carriers. Gives an
                       empirical null with no distributional assumption.
  detectable effect    for every system, the odds ratio that the observed
                       standard error would have detected at 80 percent power.
                       A null is only informative if it excludes something.
  threshold sweep      the headline count re-derived at every clonal distance
                       threshold, so the one chosen is not load-bearing.
  E-value              how strong an unmeasured lineage-linked confounder would
                       have to be to explain a surviving association away.
  heterogeneity        whether the surviving systems protect broadly across
                       phage genera or only against a few.

Outputs:
    results/tables/table16_within_cluster_effects.csv
    results/tables/table17_permutation_tests.csv
    results/tables/table18_detectable_effects.csv
    results/tables/table19_threshold_sensitivity.csv
    results/tables/table20_evalues.csv
    results/tables/table21_genus_heterogeneity.csv
    results/robustness_report.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

from papr import config, paths
from papr.models import fit_logit_cluster
from papr.stats import benjamini_hochberg, mantel_haenszel, stratified_tables

N_PERMUTATIONS = 10_000


def load():
    pairs = pd.read_csv(paths.PAIRS, compression="gzip", low_memory=False)
    strain_feats = pd.read_csv(paths.STRAIN_FEATURES, low_memory=False)
    clusters = pd.read_csv(paths.PROCESSED / "clonal_clusters.csv", index_col=0)
    pcoa = pd.read_csv(paths.PCOA, index_col=0)
    lineage = pd.read_csv(paths.TABLES / "table7_lineage_controlled_effects.csv")
    ds_cols = [c for c in strain_feats.columns if c.startswith("ds_")]
    tested = strain_feats[strain_feats["strain"].isin(set(pairs["strain"]))]
    systems = sorted([c for c in ds_cols
                      if tested[c].sum() >= config.MIN_CARRIER_STRAINS])
    pairs = pairs.merge(strain_feats[["strain"] + systems], on="strain", how="left")
    pairs["clonal_cluster"] = pairs["strain"].map(clusters["clonal_cluster"])
    return pairs, strain_feats, systems, pcoa, lineage


def within_cluster(pairs: pd.DataFrame, systems: list[str]) -> pd.DataFrame:
    """Compare carriers with non-carriers inside the same clonal cluster.

    Strata are phage by clonal cluster. A stratum is informative only if the
    same near-identical genomic background appears both with and without the
    system, so this uses a small slice of the data and answers a sharper
    question: among bacteria that are otherwise nearly the same, does carrying
    the system change the outcome.
    """
    rows = []
    for col in systems:
        name = col[3:]
        # Count informative clusters before testing, since that is the binding
        # constraint on what this design can see.
        per_cluster = pairs.groupby("clonal_cluster")[col].agg(["mean", "size"])
        discordant = per_cluster[(per_cluster["mean"] > 0) & (per_cluster["mean"] < 1)]
        tabs = stratified_tables(pairs, ["phage", "clonal_cluster"], col, "lytic")
        eff = mantel_haenszel(tabs)
        rows.append(dict(
            system=name,
            n_discordant_clusters=int(len(discordant)),
            n_strains_in_discordant=int(discordant["size"].sum() / 96) if len(discordant) else 0,
            odds_ratio=eff.estimate, ci_low=eff.ci_low, ci_high=eff.ci_high,
            p_value=eff.p_value, n_strata_used=eff.n_strata_used))
    out = pd.DataFrame(rows)
    out["q_value"] = benjamini_hochberg(out["p_value"].to_numpy())
    out["significant"] = out["q_value"] <= config.FDR_ALPHA
    out["log2_or"] = np.log2(out["odds_ratio"])
    return out


def permutation_tests(pairs: pd.DataFrame, strain_feats: pd.DataFrame,
                      systems: list[str]) -> pd.DataFrame:
    """Shuffle carrier labels inside a lineage and rebuild the null.

    The statistic is the difference in mean per-strain lysis rate between
    carriers and non-carriers. Every strain was challenged with the same 96
    phages, so per-strain lysis rates are directly comparable and no adjustment
    for phage composition is needed. Permuting inside a phylogroup, and more
    strictly inside a Sequence Type, keeps the lineage structure intact, so a
    system that merely marks a lineage cannot produce a low p-value.
    """
    rng = np.random.default_rng(config.RANDOM_SEED)
    per_strain = pairs.groupby("strain")["lytic"].mean()
    feats = strain_feats.set_index("strain").loc[per_strain.index]
    y = per_strain.to_numpy()

    rows = []
    for level_name, key in [("phylogroup", "phylogroup"), ("sequence_type", "sequence_type")]:
        groups = feats[key].astype("string").fillna("unknown").to_numpy()
        uniq = np.unique(groups)
        index_by_group = {g: np.where(groups == g)[0] for g in uniq}
        for col in systems:
            name = col[3:]
            x = feats[col].to_numpy().astype(int)
            if x.sum() == 0 or x.sum() == len(x):
                continue
            observed = y[x == 1].mean() - y[x == 0].mean()

            # Only groups containing both states can contribute variation.
            usable = [g for g in uniq
                      if 0 < x[index_by_group[g]].sum() < len(index_by_group[g])]
            n_usable_strains = int(sum(len(index_by_group[g]) for g in usable))
            if not usable:
                rows.append(dict(level=level_name, system=name,
                                 observed_difference=observed,
                                 n_usable_groups=0, n_usable_strains=0,
                                 p_value=np.nan, null_mean=np.nan, null_sd=np.nan))
                continue

            null = np.empty(N_PERMUTATIONS)
            xperm = x.copy()
            for b in range(N_PERMUTATIONS):
                for g in usable:
                    idx = index_by_group[g]
                    xperm[idx] = rng.permutation(x[idx])
                null[b] = y[xperm == 1].mean() - y[xperm == 0].mean()
            # Two-sided empirical p-value with the observed value included.
            p = (np.sum(np.abs(null - null.mean()) >= abs(observed - null.mean())) + 1) \
                / (N_PERMUTATIONS + 1)
            rows.append(dict(level=level_name, system=name,
                             observed_difference=float(observed),
                             n_usable_groups=len(usable),
                             n_usable_strains=n_usable_strains,
                             p_value=float(p), null_mean=float(null.mean()),
                             null_sd=float(null.std(ddof=1))))
    out = pd.DataFrame(rows)
    out["q_value"] = np.nan
    for level in out["level"].unique():
        m = out["level"] == level
        out.loc[m, "q_value"] = benjamini_hochberg(out.loc[m, "p_value"].to_numpy())
    out["significant"] = out["q_value"] <= config.FDR_ALPHA
    return out


def detectable_effects(lineage: pd.DataFrame) -> pd.DataFrame:
    """The odds ratio each test would have detected at 80 percent power.

    Derived from the standard error the fit actually produced, so it reflects
    this design and this amount of clustering rather than an idealised
    calculation. A system reported as null with a detectable effect close to 1
    is genuinely uninformative about moderate effects.
    """
    full = lineage[lineage["adjustment"] == "plus_phylogroup_and_coordinates"].copy()
    z_alpha, z_power = 1.959964, 0.841621
    full["min_detectable_or"] = np.exp((z_alpha + z_power) * full["std_error"])
    full["min_detectable_or_protective"] = 1.0 / full["min_detectable_or"]
    full["excludes_or_below"] = full["min_detectable_or_protective"]
    return full[["system", "odds_ratio", "ci_low", "ci_high", "std_error",
                 "p_value", "q_value", "significant",
                 "min_detectable_or", "min_detectable_or_protective"]]


def threshold_sweep(pairs: pd.DataFrame, strains: list[str], pcoa: pd.DataFrame,
                    systems: list[str]) -> pd.DataFrame:
    """Re-derive the headline count at every clonal distance threshold."""
    dist = pd.read_csv(paths.CORE_DISTANCES, sep="\t", index_col=0)
    dist.index = dist.index.astype(str)
    dist.columns = dist.columns.astype(str)
    arr = dist.loc[strains, strains].to_numpy(dtype=float)
    arr = 0.5 * (arr + arr.T)
    np.fill_diagonal(arr, 0.0)
    link = hierarchy.linkage(squareform(arr, checks=False), method="single")

    const = np.ones((len(pairs), 1))
    phage = pd.get_dummies(pairs["phage"].astype("string"),
                           drop_first=True, dtype=float).to_numpy()
    lps = pd.get_dummies(pairs["lps_core"].astype("string"),
                         drop_first=True, dtype=float).to_numpy()
    phylo = pd.get_dummies(pairs["phylogroup"].astype("string"),
                           drop_first=True, dtype=float).to_numpy()
    coords = pcoa.reindex(pairs["strain"]).to_numpy(dtype=float)
    coords = (coords - coords.mean(axis=0)) / coords.std(axis=0)
    y = pairs["lytic"].to_numpy(dtype=float)
    base = np.hstack([const, phage, lps, phylo, coords])

    rows = []
    for t in [1e-5, 5e-5, 1e-4, 2.5e-4, 5e-4, 1e-3, 2.5e-3, 5e-3, 1e-2]:
        labels = hierarchy.fcluster(link, t=t, criterion="distance")
        mapping = dict(zip(strains, labels))
        groups = pairs["strain"].map(mapping).to_numpy()
        pvals, ors = [], []
        for col in systems:
            x = np.hstack([pairs[col].to_numpy(dtype=float).reshape(-1, 1), base])
            res = fit_logit_cluster(y, x, groups)
            ors.append(float(np.exp(res.params[0])))
            pvals.append(float(res.pvalues[0]))
        q = benjamini_hochberg(np.array(pvals))
        sig = q <= config.FDR_ALPHA
        names = [c[3:] for c in systems]
        prot = [n for n, s, o in zip(names, sig, ors) if s and o < 1]
        inv = [n for n, s, o in zip(names, sig, ors) if s and o > 1]
        rows.append(dict(threshold=t, n_clusters=int(len(set(labels))),
                         n_significant=int(sig.sum()),
                         n_protective=len(prot), n_increased_lysis=len(inv),
                         protective_systems=",".join(sorted(prot)),
                         increased_lysis_systems=",".join(sorted(inv))))
    return pd.DataFrame(rows)


def e_values(lineage: pd.DataFrame, baseline_risk: float) -> pd.DataFrame:
    """How strong an unmeasured confounder would need to be.

    The E-value of VanderWeele and Ding: the minimum association, on the risk
    ratio scale, that an unmeasured confounder would need with both the
    exposure and the outcome to fully account for an observed association. The
    odds ratio is first converted to an approximate risk ratio, because the
    outcome here is common at about 21 percent and an odds ratio would overstate
    the risk ratio.
    """
    def to_rr(odds_ratio: float) -> float:
        return odds_ratio / ((1 - baseline_risk) + baseline_risk * odds_ratio)

    def e(rr: float) -> float:
        if not np.isfinite(rr) or rr <= 0:
            return np.nan
        if rr < 1:
            rr = 1.0 / rr
        return rr + np.sqrt(rr * (rr - 1))

    full = lineage[lineage["adjustment"] == "plus_phylogroup_and_coordinates"].copy()
    rows = []
    for r in full.itertuples(index=False):
        rr = to_rr(r.odds_ratio)
        # The confidence limit nearer the null is the one that must be moved.
        bound = r.ci_high if r.odds_ratio < 1 else r.ci_low
        rr_bound = to_rr(bound)
        rows.append(dict(system=r.system, odds_ratio=r.odds_ratio,
                         approx_risk_ratio=rr, e_value_estimate=e(rr),
                         ci_bound_nearer_null=bound,
                         e_value_ci_bound=(np.nan if (r.odds_ratio < 1) == (rr_bound > 1)
                                           else e(rr_bound)),
                         significant=bool(r.significant)))
    return pd.DataFrame(rows).sort_values("e_value_estimate", ascending=False)


def genus_heterogeneity(pairs: pd.DataFrame, survivors: list[str]) -> pd.DataFrame:
    """Do the surviving systems protect across genera, or only against some?"""
    rows = []
    for name in survivors:
        col = f"ds_{name}"
        logs, weights, details = [], [], []
        for g, frame in pairs.groupby("genus"):
            eff = mantel_haenszel(stratified_tables(frame, ["phage"], col, "lytic"))
            if not np.isfinite(eff.estimate) or eff.estimate <= 0:
                continue
            se = (np.log(eff.ci_high) - np.log(eff.ci_low)) / (2 * 1.959964)
            if not np.isfinite(se) or se <= 0:
                continue
            logs.append(np.log(eff.estimate))
            weights.append(1.0 / se**2)
            details.append(dict(system=name, genus=g, odds_ratio=eff.estimate,
                                ci_low=eff.ci_low, ci_high=eff.ci_high,
                                n_strata_used=eff.n_strata_used,
                                n_pairs=int(len(frame))))
        logs_a, w = np.array(logs), np.array(weights)
        pooled = float(np.sum(w * logs_a) / np.sum(w))
        q_stat = float(np.sum(w * (logs_a - pooled) ** 2))
        df = len(logs_a) - 1
        q_p = float(stats.chi2.sf(q_stat, df)) if df > 0 else np.nan
        i2 = float(max(0.0, (q_stat - df) / q_stat)) if q_stat > 0 else 0.0
        for d in details:
            d.update(pooled_odds_ratio=float(np.exp(pooled)), cochran_q=q_stat,
                     q_df=df, q_p_value=q_p, i_squared=i2,
                     n_genera_with_estimate=len(logs_a))
            rows.append(d)
    return pd.DataFrame(rows)


def main() -> int:
    paths.bootstrap()
    pairs, strain_feats, systems, pcoa, lineage = load()
    strains = sorted(pairs["strain"].unique())
    full = lineage[lineage["adjustment"] == "plus_phylogroup_and_coordinates"]
    survivors = sorted(full[full["significant"] & full["direction_preserved"]]["system"])
    baseline = float(pairs["lytic"].mean())

    print("within-cluster analysis")
    wc = within_cluster(pairs, systems)
    print("permutation tests")
    perm = permutation_tests(pairs, strain_feats, systems)
    print("detectable effects")
    mde = detectable_effects(lineage)
    print("clonal threshold sweep")
    sweep = threshold_sweep(pairs, strains, pcoa, systems)
    print("E-values")
    ev = e_values(lineage, baseline)
    print("genus heterogeneity")
    het = genus_heterogeneity(pairs, survivors)

    wc.to_csv(paths.TABLES / "table16_within_cluster_effects.csv", index=False)
    perm.to_csv(paths.TABLES / "table17_permutation_tests.csv", index=False)
    mde.to_csv(paths.TABLES / "table18_detectable_effects.csv", index=False)
    sweep.to_csv(paths.TABLES / "table19_threshold_sensitivity.csv", index=False)
    ev.to_csv(paths.TABLES / "table20_evalues.csv", index=False)
    het.to_csv(paths.TABLES / "table21_genus_heterogeneity.csv", index=False)

    out: list[str] = []
    add = out.append
    add("Robustness and sensitivity report")
    add("Generated by scripts/09_robustness_and_sensitivity.py")
    add("")
    add(f"  systems carried forward: {len(systems)}")
    add(f"  systems surviving lineage adjustment in script 04: "
        f"{', '.join(survivors) if survivors else 'none'}")

    add("")
    add("1. WITHIN CLONAL CLUSTER, COMPARING NEAR-IDENTICAL GENOMES")
    add("  Strata are phage by clonal cluster, so only bacteria sharing a "
        "near-identical core genome are compared with one another.")
    add(f"  {'system':24s} {'disc clusters':>13s} {'OR':>7s} {'95 pct CI':>16s} "
        f"{'q':>9s} {'strata':>7s}")
    for r in wc.sort_values("odds_ratio").itertuples(index=False):
        mark = " *" if r.significant else "  "
        if not np.isfinite(r.odds_ratio):
            add(f"  {r.system:24s} {r.n_discordant_clusters:13d} "
                f"{'n/a':>7s} {'not estimable':>16s} {'':>9s} {r.n_strata_used:7d}")
            continue
        add(f"  {r.system:24s} {r.n_discordant_clusters:13d} {r.odds_ratio:7.3f} "
            f"[{r.ci_low:6.3f}, {r.ci_high:6.3f}] {r.q_value:9.2e} "
            f"{r.n_strata_used:7d}{mark}")
    wc_sig = wc[wc["significant"]]
    add("")
    add(f"  Significant at q <= {config.FDR_ALPHA}: {len(wc_sig)} of {len(wc)} "
        f"({', '.join(sorted(wc_sig['system'])) if len(wc_sig) else 'none'})")
    for name in survivors:
        if name in set(wc["system"]):
            r = wc.set_index("system").loc[name]
            add(f"  {name}: {int(r['n_discordant_clusters'])} discordant clusters, "
                f"OR {r['odds_ratio']:.3f} [{r['ci_low']:.3f}, {r['ci_high']:.3f}], "
                f"q = {r['q_value']:.2e}")
    add("  This design is the strictest available and also the least powerful: "
        "most clonal clusters are small and uniform in defense content, so few "
        "strata are informative. Treat an agreement here as corroboration and a "
        "null here as weak evidence against.")

    add("")
    add("2. PERMUTATION INSIDE A LINEAGE")
    add("  Carrier labels shuffled among strains of the same phylogroup, then "
        "the same Sequence Type, holding the number of carriers fixed. The "
        "statistic is the difference in mean per-strain lysis rate.")
    for level in ["phylogroup", "sequence_type"]:
        sub = perm[perm["level"] == level]
        sig = sub[sub["significant"]]
        add(f"  permuting within {level}: {len(sig)} of {len(sub)} systems at "
            f"q <= {config.FDR_ALPHA} "
            f"({', '.join(sorted(sig['system'])) if len(sig) else 'none'})")
    add("")
    add(f"  {'system':24s} {'level':15s} {'obs diff':>9s} {'null sd':>8s} "
        f"{'groups':>7s} {'p':>9s} {'q':>9s}")
    focus = sorted(set(survivors) | set(perm[perm["significant"]]["system"]))
    for name in focus:
        for level in ["phylogroup", "sequence_type"]:
            sub = perm[(perm["system"] == name) & (perm["level"] == level)]
            if sub.empty:
                continue
            r = sub.iloc[0]
            if not np.isfinite(r["p_value"]):
                add(f"  {name:24s} {level:15s} {'no usable groups':>35s}")
                continue
            mark = " *" if r["significant"] else "  "
            add(f"  {name:24s} {level:15s} {r['observed_difference']:+9.4f} "
                f"{r['null_sd']:8.4f} {int(r['n_usable_groups']):7d} "
                f"{r['p_value']:9.2e} {r['q_value']:9.2e}{mark}")
    add("  A system whose association is purely a marker of lineage should look "
        "ordinary against this null, because the shuffling happens inside the "
        "lineage and cannot break a lineage-level pattern.")

    add("")
    add("3. WHAT THE NULLS CAN EXCLUDE")
    add("  For each system, the odds ratio its own standard error would have "
        "detected at 80 percent power.")
    nulls = mde[~mde["significant"]].sort_values("min_detectable_or_protective")
    add(f"  {'system':24s} {'OR':>7s} {'detectable protective OR':>25s}")
    for r in nulls.itertuples(index=False):
        add(f"  {r.system:24s} {r.odds_ratio:7.3f} "
            f"{r.min_detectable_or_protective:25.3f}")
    worst = nulls["min_detectable_or_protective"].max()
    best = nulls["min_detectable_or_protective"].min()
    add("")
    add(f"  Across the {len(nulls)} systems reported as null, the detectable "
        f"protective odds ratio ranges from {best:.2f} to {worst:.2f}. A system "
        f"whose figure is near {worst:.2f} is only able to exclude very large "
        "protective effects, so calling it null says little about a moderate "
        "one. This is the honest reading of most of the nulls in this analysis.")

    add("")
    add("4. SENSITIVITY TO THE CLONAL DISTANCE THRESHOLD")
    add(f"  {'threshold':>10s} {'clusters':>9s} {'significant':>12s} "
        f"{'protective':>11s} {'increased':>10s}  systems")
    for r in sweep.itertuples(index=False):
        mark = " <-" if abs(r.threshold - config.CLONAL_DISTANCE_THRESHOLD) < 1e-12 else "   "
        add(f"  {r.threshold:10.5f} {r.n_clusters:9d} {r.n_significant:12d} "
            f"{r.n_protective:11d} {r.n_increased_lysis:10d}{mark} "
            f"{r.protective_systems}")
    stable = sweep["protective_systems"].nunique()
    add("")
    add(f"  Distinct sets of protective systems across the nine thresholds: "
        f"{stable}. The threshold changes only how the standard errors are "
        "clustered, not the point estimates, so this sweep tests whether the "
        "conclusion depends on that choice.")

    add("")
    add("5. HOW STRONG AN UNMEASURED CONFOUNDER WOULD HAVE TO BE")
    add(f"  Outcome risk in the collection: {100 * baseline:.2f} percent, which "
        "is common enough that odds ratios are converted to approximate risk "
        "ratios before the E-value is computed.")
    add(f"  {'system':24s} {'OR':>7s} {'approx RR':>10s} {'E-value':>8s} "
        f"{'E-value at CI':>14s}")
    for r in ev.head(8).itertuples(index=False):
        eci = "   n/a" if not np.isfinite(r.e_value_ci_bound) else f"{r.e_value_ci_bound:14.2f}"
        mark = " *" if r.significant else "  "
        add(f"  {r.system:24s} {r.odds_ratio:7.3f} {r.approx_risk_ratio:10.3f} "
            f"{r.e_value_estimate:8.2f} {eci}{mark}")
    for name in survivors:
        r = ev.set_index("system").loc[name]
        add(f"  {name}: an unmeasured confounder would need to be associated "
            f"with both carriage and lysis by a risk ratio of at least "
            f"{r['e_value_estimate']:.2f} to explain this away, and at least "
            f"{r['e_value_ci_bound']:.2f} to move the confidence interval to "
            "the null.")
    add("  Phylogroup and core-genome position are already adjusted for, so any "
        "such confounder would have to act beyond the lineage structure those "
        "capture.")

    add("")
    add("6. ARE THE SURVIVORS BROADLY PROTECTIVE")
    for name in survivors:
        sub = het[het["system"] == name]
        if sub.empty:
            continue
        first = sub.iloc[0]
        below = int((sub["odds_ratio"] < 1).sum())
        add("")
        add(f"  {name}: estimable in {int(first['n_genera_with_estimate'])} genera, "
            f"odds ratio below 1 in {below} of them. "
            f"Pooled {first['pooled_odds_ratio']:.3f}, Cochran Q = "
            f"{first['cochran_q']:.1f} on {int(first['q_df'])} df "
            f"(p = {first['q_p_value']:.3g}), I squared = "
            f"{100 * first['i_squared']:.0f} percent.")
        add(f"    {'genus':22s} {'OR':>7s} {'95 pct CI':>16s} {'pairs':>6s}")
        for r in sub.sort_values("odds_ratio").itertuples(index=False):
            add(f"    {r.genus:22s} {r.odds_ratio:7.3f} "
                f"[{r.ci_low:6.3f}, {r.ci_high:6.3f}] {r.n_pairs:6d}")
    add("")
    add("  High heterogeneity would mean a system protects against some phage "
        "lineages and not others, which is what a mechanistic defense should "
        "look like. Low heterogeneity with a consistent direction is also "
        "compatible with a general effect on infection efficiency.")

    text = "\n".join(out)
    (paths.RESULTS / "robustness_report.txt").write_text(text + "\n", encoding="utf-8")
    print("\n" + text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
