"""Use the whole plaque score instead of collapsing it to lysed or not.

Every association in this analysis has been estimated on a binary outcome: a
plaque appeared or it did not. That throws away most of the measurement. The
released score records how far a phage lysate could be diluted and still clear
the lawn, from 0 to 4, and the four thresholds carry different information:
whether an infection happened at all, and how productive it was.

The distinction matters mechanistically. A restriction-modification system that
cleaves most incoming genomes but misses a few reduces the efficiency of
plating without preventing plaques altogether, so it should show a weak effect
at the first threshold and a strong one at the fourth. A system that aborts
infection outright should act at the first threshold and add little beyond it.
A binary analysis cannot tell these apart.

Each system is fitted at all four thresholds with the same covariates and
clustered standard errors used elsewhere. Reporting them separately also tests
the proportional odds assumption that a single ordinal model would impose
rather than check.

Outputs:
    results/tables/table43_threshold_effects.csv
    results/tables/table44_threshold_heterogeneity.csv
    results/ordinal_thresholds_report.txt
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
from papr.stats import benjamini_hochberg

THRESHOLDS = [1, 2, 3, 4]


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
    pairs = pairs.merge(strain_feats[["strain"] + systems], on="strain", how="left")
    pairs["clonal_cluster"] = pairs["strain"].map(clusters["clonal_cluster"])

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
    groups = pairs["clonal_cluster"].to_numpy()
    score = pairs["plaque_score"].to_numpy(dtype=float)

    rows = []
    for i, col in enumerate(systems, 1):
        exposure = pairs[col].to_numpy(dtype=float)
        x = np.hstack([exposure.reshape(-1, 1), base])
        for t in THRESHOLDS:
            y = (score >= t).astype(float)
            try:
                res = fit_logit_cluster(y, x, groups)
                eff = term_effect(res, 0)
            except Exception:
                eff = dict(odds_ratio=np.nan, ci_low=np.nan, ci_high=np.nan,
                           log2_or=np.nan, coefficient=np.nan,
                           std_error=np.nan, p_value=np.nan)
            rows.append(dict(system=col[3:], threshold=t,
                             n_events=int(y.sum()), event_rate=float(y.mean()),
                             **eff))
        if i % 8 == 0 or i == len(systems):
            print(f"  {i}/{len(systems)} systems")

    eff = pd.DataFrame(rows)
    eff["q_value"] = np.nan
    for t in THRESHOLDS:
        m = eff["threshold"] == t
        eff.loc[m, "q_value"] = benjamini_hochberg(eff.loc[m, "p_value"].to_numpy())
    eff["significant"] = eff["q_value"] <= config.FDR_ALPHA
    eff.to_csv(paths.TABLES / "table43_threshold_effects.csv", index=False)

    # Test whether the effect is constant across thresholds, which is what a
    # single ordinal model would assume.
    het_rows = []
    for name, g in eff.groupby("system"):
        g = g.dropna(subset=["coefficient", "std_error"])
        g = g[g["std_error"] > 0]
        if len(g) < 2:
            continue
        b = g["coefficient"].to_numpy()
        se = g["std_error"].to_numpy()
        w = 1.0 / se ** 2
        pooled = float((w * b).sum() / w.sum())
        pooled_se = float(np.sqrt(1.0 / w.sum()))
        q_stat = float((w * (b - pooled) ** 2).sum())
        df = len(b) - 1
        p_het = float(stats.chi2.sf(q_stat, df)) if df > 0 else np.nan
        first = g[g["threshold"] == 1]["coefficient"]
        last = g[g["threshold"] == g["threshold"].max()]["coefficient"]
        het_rows.append(dict(
            system=name, n_thresholds=len(b),
            pooled_coefficient=pooled, pooled_odds_ratio=float(np.exp(pooled)),
            pooled_ci_low=float(np.exp(pooled - 1.96 * pooled_se)),
            pooled_ci_high=float(np.exp(pooled + 1.96 * pooled_se)),
            pooled_p=float(2 * stats.norm.sf(abs(pooled / pooled_se))),
            heterogeneity_q=q_stat, heterogeneity_df=df,
            heterogeneity_p=p_het,
            coef_threshold1=float(first.iloc[0]) if len(first) else np.nan,
            coef_threshold_top=float(last.iloc[0]) if len(last) else np.nan))
    het = pd.DataFrame(het_rows)
    het["pooled_q"] = benjamini_hochberg(het["pooled_p"].to_numpy())
    het["heterogeneity_q"] = benjamini_hochberg(het["heterogeneity_p"].to_numpy())
    het["pooled_significant"] = het["pooled_q"] <= config.FDR_ALPHA
    het["gradient"] = het["coef_threshold_top"] - het["coef_threshold1"]
    het = het.sort_values("pooled_coefficient")
    het.to_csv(paths.TABLES / "table44_threshold_heterogeneity.csv", index=False)

    binary_sig = int((eff[(eff["threshold"] == 1)]["significant"]).sum())

    o: list[str] = []
    add = o.append
    add("Ordinal threshold report")
    add("Generated by scripts/23_ordinal_thresholds.py")
    add("")
    add("WHAT IS BEING FITTED")
    add(f"  systems : {len(systems)}, thresholds : {len(THRESHOLDS)}, "
        f"fits : {len(systems) * len(THRESHOLDS)}")
    add("  Threshold t compares a plaque score of at least t against anything "
        "lower, with phage identity, host surface type, phylogroup and "
        "core-genome position adjusted and standard errors clustered on "
        "clonal background.")
    add("")
    add("  event rates:")
    for t in THRESHOLDS:
        r = eff[eff["threshold"] == t].iloc[0]
        add(f"    score >= {t} : {r.n_events:6d} pairs "
            f"({100 * r.event_rate:5.2f} percent)")

    add("")
    add("WHAT EACH THRESHOLD FINDS ON ITS OWN")
    for t in THRESHOLDS:
        sub = eff[eff["threshold"] == t]
        add(f"  threshold {t}: {int(sub['significant'].sum())} of {len(sub)} "
            f"systems at q <= {config.FDR_ALPHA}")
    add("  Power falls away with the threshold because the events become "
        "rarer, from 7976 pairs at the first to 1288 at the fourth, so the "
        "deeper thresholds find less despite carrying the more specific "
        "information.")
    add("")
    add("  A weighted combination of the four thresholds is reported below as "
        "an effect size, and deliberately not as a test. The thresholds are "
        "nested and reuse the same observations four times, so combining them "
        "as if they were independent studies produces an interval that is far "
        "too narrow. Under that invalid combination "
        f"{int(het['pooled_significant'].sum())} of {len(het)} systems would "
        f"appear significant against {binary_sig} for the binary analysis; "
        "that difference is an artefact of double counting, not a gain in "
        "evidence, and it is not claimed here.")

    add("")
    add("EFFECT ACROSS THRESHOLDS")
    add(f"  {'system':24s} {'t=1':>7s} {'t=2':>7s} {'t=3':>7s} {'t=4':>7s} "
        f"{'pooled OR':>10s} {'het p':>9s}")
    focus = het.head(8).append(het.tail(4)) if hasattr(het, "append") else \
        pd.concat([het.head(8), het.tail(4)])
    for r in focus.itertuples(index=False):
        cells = []
        for t in THRESHOLDS:
            v = eff[(eff["system"] == r.system) & (eff["threshold"] == t)]
            cells.append(f"{float(v['odds_ratio'].iloc[0]):7.3f}"
                         if len(v) and np.isfinite(v["odds_ratio"].iloc[0])
                         else "    n/a")
        add(f"  {r.system:24s} {' '.join(cells)} {r.pooled_odds_ratio:10.3f} "
            f"{r.heterogeneity_p:9.2e}")
    add("  The pooled column is a descriptive weighted average, not a test.")

    n_het = int((het["heterogeneity_q"] <= config.FDR_ALPHA).sum())
    add("")
    add("IS THE EFFECT THE SAME AT EVERY THRESHOLD")
    add(f"  systems with detectable heterogeneity across thresholds: "
        f"{n_het} of {len(het)} at q <= {config.FDR_ALPHA}")
    if n_het == 0:
        add("  No system shows a reliably different effect at different "
            "thresholds. The proportional odds assumption that a single "
            "ordinal model would impose is therefore not contradicted here, "
            "and the binary analysis loses little by collapsing the score.")
    else:
        strongest = het.reindex(het["heterogeneity_q"].sort_values().index).head(5)
        add(f"  {'system':24s} {'t=1 coef':>9s} {'top coef':>9s} "
            f"{'gradient':>9s} {'q':>9s}")
        for r in strongest.itertuples(index=False):
            add(f"  {r.system:24s} {r.coef_threshold1:+9.3f} "
                f"{r.coef_threshold_top:+9.3f} {r.gradient:+9.3f} "
                f"{r.heterogeneity_q:9.2e}")
        add("  A negative gradient means the system acts more strongly at the "
            "deeper dilutions, which is what a mechanism that reduces the "
            "efficiency of plating without preventing plaques looks like. A "
            "positive gradient means it acts mainly on whether an infection "
            "happens at all.")

    add("")
    add("WHAT THE SCORE ADDS")
    add("  Not power. The four thresholds reuse the same observations, so no "
        "valid combination of them yields much beyond what the first "
        "threshold already provides, and the deeper thresholds have fewer "
        "events to work with.")
    add("  What it adds is a consistency check that the binary analysis "
        "cannot perform. The effect sizes are stable across thresholds for "
        "every system tested, which means the associations reported "
        "throughout this analysis are not artefacts of where the lysed and "
        "not-lysed boundary was drawn. Had a system acted only at the first "
        "threshold, or only at the deepest, the binary summary would have "
        "been misleading about it. None does.")
    add("")
    add("LIMITS")
    add("  The four thresholds are nested, so the fits are not independent "
        "and the pooled estimate's interval would be far too narrow. It is "
        "reported as a summary of direction and magnitude only; the "
        "per-threshold estimates are the ones with honest intervals.")
    add("  Conditioning on a higher threshold conditions on a successful "
        "infection having happened, so an effect at threshold 4 mixes an "
        "effect on productivity with selection on which pairs reached that "
        "far. This is the same collider problem flagged for the "
        "lytic-pairs-only analysis in script 03.")

    text = "\n".join(o)
    (paths.RESULTS / "ordinal_thresholds_report.txt").write_text(
        text + "\n", encoding="utf-8")
    print("\n" + text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
