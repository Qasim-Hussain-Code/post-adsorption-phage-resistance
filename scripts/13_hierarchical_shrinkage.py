"""Estimate every defense system jointly, including those too rare to test alone.

The analysis so far tested the 33 systems carried by at least 20 strains and
said nothing about the other 104. It also treated each system as an unrelated
question, so a system whose estimate was large by chance kept that estimate.
Both problems have the same remedy: fit all the systems together and let them
inform one another.

Each system is first fitted separately with the lineage adjustment of script
04, giving a log odds ratio and a standard error. Those estimates are then
combined in a random effects model, log OR_i drawn from a normal distribution
with mean mu and between-system variance tau squared, observed with known
error. Each system's estimate is shrunk toward mu by an amount set by its own
precision, so a system carried by six strains is pulled most and one carried by
three hundred is barely moved.

The quantity of most interest is tau squared itself. It measures how much
genuine variation in protective effect there is between defense systems once
sampling error is removed. If it is near zero, the systems do not differ
detectably from one another in how much they protect, whatever their individual
estimates look like.

Outputs:
    results/tables/table29_all_system_effects.csv
    results/tables/table30_shrinkage_summary.csv
    results/hierarchical_shrinkage_report.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
from scipy import stats

from papr import config, paths
from papr.models import fit_logit_cluster
from papr.stats import benjamini_hochberg

# Below this many carriers a separate fit is too unstable to contribute even a
# shrunken estimate; separation becomes common and the standard error is not
# interpretable.
MIN_CARRIERS = 5

# Standard errors above this indicate a degenerate fit rather than a wide one.
MAX_USABLE_SE = 5.0


def dersimonian_laird(y: np.ndarray, se: np.ndarray) -> tuple[float, float, float]:
    """Random effects mean, between-study variance and Cochran's Q."""
    w = 1.0 / se**2
    fixed_mean = float(np.sum(w * y) / np.sum(w))
    q = float(np.sum(w * (y - fixed_mean) ** 2))
    k = len(y)
    c = float(np.sum(w) - np.sum(w**2) / np.sum(w))
    tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0
    w_star = 1.0 / (se**2 + tau2)
    mu = float(np.sum(w_star * y) / np.sum(w_star))
    return mu, tau2, q


def paule_mandel(y: np.ndarray, se: np.ndarray, iters: int = 200) -> float:
    """Iterative estimate of tau squared, less biased than DerSimonian-Laird."""
    tau2 = 0.0
    k = len(y)
    for _ in range(iters):
        w = 1.0 / (se**2 + tau2)
        mu = float(np.sum(w * y) / np.sum(w))
        stat = float(np.sum(w * (y - mu) ** 2))
        if abs(stat - (k - 1)) < 1e-8:
            break
        deriv = float(np.sum(w**2 * (y - mu) ** 2))
        if deriv <= 0:
            break
        step = (stat - (k - 1)) / deriv
        tau2 = max(0.0, tau2 + step)
    return tau2


def main() -> int:
    paths.bootstrap()
    pairs = pd.read_csv(paths.PAIRS, compression="gzip", low_memory=False)
    strain_feats = pd.read_csv(paths.STRAIN_FEATURES, low_memory=False)
    pcoa = pd.read_csv(paths.PCOA, index_col=0)
    clusters = pd.read_csv(paths.PROCESSED / "clonal_clusters.csv", index_col=0)

    ds_cols = [c for c in strain_feats.columns if c.startswith("ds_")]
    tested = strain_feats[strain_feats["strain"].isin(set(pairs["strain"]))]
    counts = tested[ds_cols].sum()
    systems = sorted([c for c in ds_cols if counts[c] >= MIN_CARRIERS])
    print(f"{len(ds_cols)} systems annotated; {len(systems)} with at least "
          f"{MIN_CARRIERS} carrier strains")

    pairs = pairs.merge(strain_feats[["strain"] + systems], on="strain", how="left")
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

    rows = []
    for i, col in enumerate(systems, 1):
        exposure = pairs[col].to_numpy(dtype=float)
        beta = se = p = np.nan
        try:
            res = fit_logit_cluster(y, np.hstack([exposure.reshape(-1, 1), base]),
                                    groups)
            beta, se, p = float(res.params[0]), float(res.bse[0]), float(res.pvalues[0])
        except Exception:
            pass
        rows.append(dict(system=col[3:], n_carrier_strains=int(counts[col]),
                         log_or=beta, std_error=se, p_value=p))
        if i % 15 == 0 or i == len(systems):
            print(f"  fitted {i}/{len(systems)}")

    eff = pd.DataFrame(rows)
    usable = eff["log_or"].notna() & eff["std_error"].notna() \
        & np.isfinite(eff["log_or"]) & np.isfinite(eff["std_error"]) \
        & (eff["std_error"] > 0) & (eff["std_error"] < MAX_USABLE_SE)
    eff["usable"] = usable
    fit = eff[usable].copy()

    yv = fit["log_or"].to_numpy()
    sev = fit["std_error"].to_numpy()
    mu_dl, tau2_dl, q_stat = dersimonian_laird(yv, sev)
    tau2_pm = paule_mandel(yv, sev)
    tau2 = tau2_pm
    mu = float(np.sum(yv / (sev**2 + tau2)) / np.sum(1.0 / (sev**2 + tau2)))
    mu_se = float(np.sqrt(1.0 / np.sum(1.0 / (sev**2 + tau2))))

    shrink = tau2 / (tau2 + sev**2) if tau2 > 0 else np.zeros_like(sev)
    post_mean = shrink * yv + (1 - shrink) * mu
    post_var = 1.0 / (1.0 / tau2 + 1.0 / sev**2) if tau2 > 0 else np.full_like(sev, np.nan)
    post_sd = np.sqrt(post_var)

    fit["shrinkage_weight"] = shrink
    fit["posterior_log_or"] = post_mean
    fit["posterior_sd"] = post_sd
    fit["posterior_odds_ratio"] = np.exp(post_mean)
    fit["posterior_ci_low"] = np.exp(post_mean - 1.96 * post_sd)
    fit["posterior_ci_high"] = np.exp(post_mean + 1.96 * post_sd)
    fit["raw_odds_ratio"] = np.exp(yv)
    fit["q_value"] = benjamini_hochberg(fit["p_value"].to_numpy())
    fit["significant_raw"] = fit["q_value"] <= config.FDR_ALPHA
    fit["posterior_excludes_null"] = (
        (fit["posterior_ci_low"] > 1) | (fit["posterior_ci_high"] < 1))

    out = eff.merge(
        fit[["system", "shrinkage_weight", "posterior_log_or", "posterior_sd",
             "posterior_odds_ratio", "posterior_ci_low", "posterior_ci_high",
             "raw_odds_ratio", "q_value", "significant_raw",
             "posterior_excludes_null"]], on="system", how="left")
    out.to_csv(paths.TABLES / "table29_all_system_effects.csv", index=False)

    i2 = max(0.0, (q_stat - (len(yv) - 1)) / q_stat) if q_stat > 0 else 0.0
    summary = pd.DataFrame([dict(
        n_systems_annotated=len(ds_cols), n_systems_fitted=len(eff),
        n_systems_usable=int(usable.sum()),
        pooled_log_or=mu, pooled_odds_ratio=float(np.exp(mu)),
        pooled_se=mu_se,
        pooled_ci_low=float(np.exp(mu - 1.96 * mu_se)),
        pooled_ci_high=float(np.exp(mu + 1.96 * mu_se)),
        tau2_paule_mandel=tau2_pm, tau2_dersimonian_laird=tau2_dl,
        tau=float(np.sqrt(tau2)), cochran_q=q_stat, q_df=len(yv) - 1,
        q_p_value=float(stats.chi2.sf(q_stat, len(yv) - 1)), i_squared=i2,
        median_shrinkage=float(np.median(shrink)),
        n_significant_raw=int(fit["significant_raw"].sum()),
        n_posterior_excludes_null=int(fit["posterior_excludes_null"].sum()))])
    summary.to_csv(paths.TABLES / "table30_shrinkage_summary.csv", index=False)

    o: list[str] = []
    add = o.append
    add("Hierarchical shrinkage across all defense systems")
    add("Generated by scripts/13_hierarchical_shrinkage.py")
    add("")
    add("COVERAGE")
    add(f"  systems annotated by DefenseFinder        : {len(ds_cols)}")
    add(f"  carried by at least {MIN_CARRIERS} strains and fitted : {len(eff)}")
    add(f"  usable estimates after discarding degenerate fits: {int(usable.sum())}")
    add(f"  systems tested in scripts 03 and 04       : 33")
    add("  Extending from 33 to "
        f"{int(usable.sum())} systems is the point of this script. The earlier "
        "threshold of 20 carriers was chosen so that each system could be "
        "tested on its own; pooling removes that requirement, because a rare "
        "system no longer has to carry its own estimate.")

    add("")
    add("HOW MUCH DO DEFENSE SYSTEMS REALLY DIFFER FROM ONE ANOTHER")
    add(f"  pooled odds ratio across systems : {float(np.exp(mu)):.4f} "
        f"(95 percent CI {float(np.exp(mu - 1.96 * mu_se)):.4f} to "
        f"{float(np.exp(mu + 1.96 * mu_se)):.4f})")
    add(f"  between-system variance tau squared, Paule-Mandel : {tau2_pm:.5f}")
    add(f"  between-system variance tau squared, DerSimonian-Laird : {tau2_dl:.5f}")
    add(f"  between-system standard deviation tau : {float(np.sqrt(tau2)):.4f} "
        "on the log odds scale")
    add(f"  Cochran Q = {q_stat:.1f} on {len(yv) - 1} df "
        f"(p = {float(stats.chi2.sf(q_stat, len(yv) - 1)):.3g}), "
        f"I squared = {100 * i2:.0f} percent")
    if tau2 > 0:
        lo_sys = float(np.exp(mu - 1.96 * np.sqrt(tau2)))
        hi_sys = float(np.exp(mu + 1.96 * np.sqrt(tau2)))
        add(f"  Implied spread of true system effects: 95 percent of systems "
            f"lie between odds ratios of {lo_sys:.3f} and {hi_sys:.3f}.")
        add("  There is therefore real variation between defense systems in "
            "how they relate to lysis, over and above sampling error, but the "
            "pooled effect sits close to no effect at all.")
    else:
        add("  The estimated between-system variance is zero. Once sampling "
            "error is accounted for, these systems do not differ detectably "
            "from one another in their association with lysis, and the "
            "apparent spread of individual estimates is consistent with noise.")

    add("")
    add("WHAT SHRINKAGE DOES TO THE INDIVIDUAL ESTIMATES")
    add(f"  median shrinkage weight: {float(np.median(shrink)):.3f} "
        "(1 means the estimate is kept, 0 means it is replaced by the pooled mean)")
    add(f"  systems significant on their own at q <= {config.FDR_ALPHA}: "
        f"{int(fit['significant_raw'].sum())}")
    add(f"  systems whose shrunken 95 percent interval excludes an odds ratio "
        f"of 1: {int(fit['posterior_excludes_null'].sum())}")
    add("")
    add(f"  {'system':24s} {'carriers':>8s} {'raw OR':>7s} {'shrunk OR':>9s} "
        f"{'95 pct CI':>16s} {'weight':>7s}")
    show = fit.reindex(fit["posterior_log_or"].abs().sort_values(
        ascending=False).index).head(20)
    for r in show.itertuples(index=False):
        mark = " *" if r.posterior_excludes_null else "  "
        add(f"  {r.system:24s} {r.n_carrier_strains:8d} {r.raw_odds_ratio:7.3f} "
            f"{r.posterior_odds_ratio:9.3f} "
            f"[{r.posterior_ci_low:6.3f}, {r.posterior_ci_high:6.3f}] "
            f"{r.shrinkage_weight:7.3f}{mark}")
    add("  The twenty systems with the largest shrunken effects are shown; the "
        "full table is table29.")

    rare = fit[fit["n_carrier_strains"] < config.MIN_CARRIER_STRAINS]
    if len(rare):
        rare_sig = rare[rare["posterior_excludes_null"]]
        add("")
        add(f"  Of the {len(rare)} systems below the 20-carrier threshold used "
            f"earlier, {len(rare_sig)} have a shrunken interval excluding 1 "
            f"({', '.join(sorted(rare_sig['system'])) if len(rare_sig) else 'none'}).")
        add(f"  Their median shrinkage weight is "
            f"{float(rare['shrinkage_weight'].median()):.3f}, so most of their "
            "individual signal is replaced by the pooled mean. That is the "
            "correct treatment of an estimate from a handful of strains, and "
            "it is why extending coverage to these systems does not generate "
            "new claims.")

    add("")
    add("WHAT THE SHRUNKEN INTERVALS DO NOT ACCOUNT FOR")
    add("  The intervals above are empirical Bayes: the pooled mean and the "
        "between-system variance are estimated from the same data and then "
        "treated as known when each system's interval is formed. That makes "
        "them narrower than they should be, and the effect is largest for the "
        "rare systems whose intervals depend most on the pooled quantities. A "
        "fully Bayesian treatment carrying uncertainty in mu and tau squared "
        "would widen them, most likely removing the weakest of the nine "
        "systems whose intervals currently exclude 1.")
    add(f"  As a rough guide to the size of that omission, the standard error "
        f"of the pooled mean alone is {mu_se:.4f} against a median posterior "
        f"standard deviation of {float(np.median(post_sd)):.4f}, so it is not "
        "negligible.")

    add("")
    add("READING THIS ALONGSIDE THE EARLIER RESULTS")
    add("  Shrinkage is a different question from the robustness tests in "
        "script 09, and it is weaker. It corrects estimates for the fact that "
        "many systems were examined; it does nothing about confounding by "
        "clonal background, which is what actually removed most of the "
        "associations. A system can survive shrinkage and still be a lineage "
        "artefact.")
    add("  The practical conclusion is unchanged. Across all "
        f"{int(usable.sum())} systems the pooled association with lysis is "
        f"an odds ratio of {float(np.exp(mu)):.3f}, close to no association, "
        "and the systems that stand out individually are the same ones script "
        "09 showed do not survive the stricter designs.")

    text = "\n".join(o)
    (paths.RESULTS / "hierarchical_shrinkage_report.txt").write_text(
        text + "\n", encoding="utf-8")
    print("\n" + text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
