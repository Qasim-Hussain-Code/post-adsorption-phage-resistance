"""Design-based checks that do not rely on covariate adjustment.

Shared by the robustness sweep over all tested systems (script 09) and the
consolidated candidate assessment (script 15), so both apply an identical
estimator and a difference between them can only come from the systems chosen.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .models import fit_logit_cluster
from .stats import benjamini_hochberg, mantel_haenszel, stratified_tables


def within_cluster_effects(pairs: pd.DataFrame, systems: list[str]) -> pd.DataFrame:
    """Carriers against non-carriers inside the same clonal cluster.

    Only clusters containing a system both present and absent carry information
    about it, so the fit is restricted to those and a cluster fixed effect
    absorbs the shared genomic background. Standard errors are clustered on the
    same variable. A Mantel-Haenszel estimate over the same strata is returned
    alongside; it treats pairs as independent and is miscalibrated, which the
    caller can demonstrate by counting impossible results.
    """
    y_all = pairs["lytic"].to_numpy(dtype=float)
    rows = []
    for col in systems:
        name = col[3:]
        per_cluster = pairs.groupby("clonal_cluster")[col].agg(["mean"])
        discordant = per_cluster[(per_cluster["mean"] > 0) & (per_cluster["mean"] < 1)]
        keep = set(discordant.index)

        eff = mantel_haenszel(
            stratified_tables(pairs, ["phage", "clonal_cluster"], col, "lytic"))

        mask = pairs["clonal_cluster"].isin(keep).to_numpy()
        n_strains = int(pairs.loc[mask, "strain"].nunique())
        or_reg = lo_reg = hi_reg = p_reg = np.nan
        if len(keep) >= 2 and mask.sum() > 0:
            sub = pairs.loc[mask]
            exposure = sub[col].to_numpy(dtype=float)
            if 0 < exposure.sum() < len(exposure):
                x = np.hstack([
                    exposure.reshape(-1, 1),
                    np.ones((int(mask.sum()), 1)),
                    pd.get_dummies(sub["phage"].astype("string"),
                                   drop_first=True, dtype=float).to_numpy(),
                    pd.get_dummies(sub["clonal_cluster"].astype("string"),
                                   drop_first=True, dtype=float).to_numpy(),
                ])
                try:
                    res = fit_logit_cluster(
                        y_all[mask], x, sub["clonal_cluster"].to_numpy())
                    beta, se = float(res.params[0]), float(res.bse[0])
                    if np.isfinite(beta) and np.isfinite(se) and se > 0:
                        or_reg = float(np.exp(beta))
                        lo_reg = float(np.exp(beta - 1.96 * se))
                        hi_reg = float(np.exp(beta + 1.96 * se))
                        p_reg = float(res.pvalues[0])
                except Exception:
                    pass

        rows.append(dict(
            system=name,
            n_discordant_clusters=int(len(discordant)),
            n_strains_in_discordant=n_strains,
            mh_odds_ratio=eff.estimate, mh_ci_low=eff.ci_low,
            mh_ci_high=eff.ci_high, mh_p_value=eff.p_value,
            n_strata_used=eff.n_strata_used,
            odds_ratio=or_reg, ci_low=lo_reg, ci_high=hi_reg, p_value=p_reg))

    out = pd.DataFrame(rows)
    out["q_value"] = benjamini_hochberg(out["p_value"].to_numpy())
    out["mh_q_value"] = benjamini_hochberg(out["mh_p_value"].to_numpy())
    out["significant"] = out["q_value"] <= 0.05
    out["mh_significant"] = out["mh_q_value"] <= 0.05
    out["log2_or"] = np.log2(out["odds_ratio"])
    return out


def permutation_effects(pairs: pd.DataFrame, strain_feats: pd.DataFrame,
                        systems: list[str], levels: list[tuple[str, str]],
                        n_permutations: int, seed: int) -> pd.DataFrame:
    """Shuffle carrier labels inside a lineage and rebuild the null.

    The statistic is the difference in mean per-strain lysis rate between
    carriers and non-carriers. Every strain was challenged with the same 96
    phages, so per-strain rates are comparable without adjustment. Shuffling
    inside the lineage means a system that merely marks a lineage cannot
    produce a low p-value.
    """
    rng = np.random.default_rng(seed)
    per_strain = pairs.groupby("strain")["lytic"].mean()
    feats = strain_feats.set_index("strain").loc[per_strain.index]
    y = per_strain.to_numpy()

    rows = []
    for level_name, key in levels:
        groups = feats[key].astype("string").fillna("unknown").to_numpy()
        uniq = np.unique(groups)
        index_by_group = {g: np.where(groups == g)[0] for g in uniq}
        for col in systems:
            name = col[3:]
            x = feats[col].to_numpy().astype(int)
            if x.sum() == 0 or x.sum() == len(x):
                continue
            observed = y[x == 1].mean() - y[x == 0].mean()

            usable = [g for g in uniq
                      if 0 < x[index_by_group[g]].sum() < len(index_by_group[g])]
            n_usable_strains = int(sum(len(index_by_group[g]) for g in usable))
            if not usable:
                rows.append(dict(level=level_name, system=name,
                                 observed_difference=float(observed),
                                 n_usable_groups=0, n_usable_strains=0,
                                 p_value=np.nan, null_mean=np.nan,
                                 null_sd=np.nan))
                continue

            null = np.empty(n_permutations)
            xperm = x.copy()
            for b in range(n_permutations):
                for g in usable:
                    idx = index_by_group[g]
                    xperm[idx] = rng.permutation(x[idx])
                null[b] = y[xperm == 1].mean() - y[xperm == 0].mean()
            p = (np.sum(np.abs(null - null.mean()) >= abs(observed - null.mean())) + 1) \
                / (n_permutations + 1)
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
    out["significant"] = out["q_value"] <= 0.05
    return out
