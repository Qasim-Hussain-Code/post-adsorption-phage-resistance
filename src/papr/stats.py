"""Closed-form and rank-based estimators used throughout the analysis.

Every function here is deterministic given its inputs, except the two that take
an explicit random seed for resampling. Nothing in this module fits a
predictive model.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd
from scipy import stats


@dataclass
class Effect:
    """An odds ratio with an interval, a p-value and the data behind it."""

    estimate: float
    ci_low: float
    ci_high: float
    p_value: float
    n_exposed: int
    n_unexposed: int
    n_events_exposed: int
    n_events_unexposed: int
    n_strata_used: int = 1
    method: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def _counts(exposure: np.ndarray, outcome: np.ndarray) -> tuple[int, int, int, int]:
    e = np.asarray(exposure).astype(bool)
    y = np.asarray(outcome).astype(bool)
    a = int(np.sum(e & y))
    b = int(np.sum(e & ~y))
    c = int(np.sum(~e & y))
    d = int(np.sum(~e & ~y))
    return a, b, c, d


def fisher_odds_ratio(exposure, outcome, method: str = "fisher_exact") -> Effect:
    """Unstratified odds ratio with a conditional exact p-value.

    The interval is the Woolf (log) interval with a Haldane-Anscombe
    continuity correction of 0.5 applied only when a cell is empty, so that an
    interval is still reported for a sparse table instead of an infinity.
    """
    a, b, c, d = _counts(exposure, outcome)
    table = [[a, b], [c, d]]
    _, p = stats.fisher_exact(table)
    if min(a, b, c, d) == 0:
        aa, bb, cc, dd = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    else:
        aa, bb, cc, dd = float(a), float(b), float(c), float(d)
    or_hat = (aa * dd) / (bb * cc)
    se = np.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return Effect(
        estimate=float(or_hat),
        ci_low=float(np.exp(np.log(or_hat) - 1.96 * se)),
        ci_high=float(np.exp(np.log(or_hat) + 1.96 * se)),
        p_value=float(p),
        n_exposed=a + b,
        n_unexposed=c + d,
        n_events_exposed=a,
        n_events_unexposed=c,
        n_strata_used=1,
        method=method,
    )


def mantel_haenszel(tables: list[tuple[int, int, int, int]],
                    method: str = "cochran_mantel_haenszel") -> Effect:
    """Cochran-Mantel-Haenszel pooled odds ratio.

    The standard error is the Robins, Breslow and Greenland estimator, which is
    consistent both when strata are few and large and when they are many and
    sparse. The p-value is the Cochran-Mantel-Haenszel chi-squared test with a
    continuity correction. Strata in which either the exposure or the outcome
    does not vary carry no information about the common odds ratio and are
    dropped; the number retained is reported so that the reader can see how
    much of the data each estimate rests on.
    """
    R = S = 0.0
    vr = vs = vrs = 0.0
    num = den = 0.0
    used = 0
    tot = [0, 0, 0, 0]
    for a, b, c, d in tables:
        n = a + b + c + d
        tot[0] += a
        tot[1] += b
        tot[2] += c
        tot[3] += d
        if n == 0 or (a + b) == 0 or (c + d) == 0 or (a + c) == 0 or (b + d) == 0:
            continue
        used += 1
        R += a * d / n
        S += b * c / n
        # Robins, Breslow and Greenland write the variance of the pooled log
        # odds ratio as
        #     sum(P R) / (2 R^2) + sum(P S + Q R) / (2 R S) + sum(Q S) / (2 S^2)
        # with P = (a+d)/n and Q = (b+c)/n. The three sums are accumulated
        # here and divided below. The middle sum carries no factor of one half
        # of its own: writing one here would understate the standard error by
        # about a tenth and give intervals that are correspondingly too narrow.
        vr += (a + d) * a * d / n**2
        vs += (b + c) * b * c / n**2
        vrs += ((a + d) * b * c + (b + c) * a * d) / n**2
        num += a - (a + b) * (a + c) / n
        den += (a + b) * (c + d) * (a + c) * (b + d) / (n**2 * (n - 1)) if n > 1 else 0.0

    a_t, b_t, c_t, d_t = tot
    if used == 0 or R == 0 or S == 0 or den == 0:
        return Effect(np.nan, np.nan, np.nan, np.nan, a_t + b_t, c_t + d_t,
                      a_t, c_t, used, method)
    or_mh = R / S
    se = np.sqrt(vr / (2 * R**2) + vrs / (2 * R * S) + vs / (2 * S**2))
    chi2 = (abs(num) - 0.5) ** 2 / den
    p = float(stats.chi2.sf(chi2, df=1))
    return Effect(
        estimate=float(or_mh),
        ci_low=float(np.exp(np.log(or_mh) - 1.96 * se)),
        ci_high=float(np.exp(np.log(or_mh) + 1.96 * se)),
        p_value=p,
        n_exposed=a_t + b_t,
        n_unexposed=c_t + d_t,
        n_events_exposed=a_t,
        n_events_unexposed=c_t,
        n_strata_used=used,
        method=method,
    )


def stratified_tables(frame: pd.DataFrame, strata: list[str],
                      exposure: str, outcome: str) -> list[tuple[int, int, int, int]]:
    """Two-by-two tables of exposure against outcome, one per stratum."""
    out: list[tuple[int, int, int, int]] = []
    for _, g in frame.groupby(strata, observed=True, sort=False):
        out.append(_counts(g[exposure].to_numpy(), g[outcome].to_numpy()))
    return out


def benjamini_hochberg(p_values) -> np.ndarray:
    """Benjamini-Hochberg step-up q-values, order preserved, NaNs passed through."""
    p = np.asarray(p_values, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


def mann_whitney(a, b) -> dict:
    """Two-sided Mann-Whitney U with a rank-biserial effect size."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if a.size == 0 or b.size == 0:
        return dict(statistic=np.nan, p_value=np.nan, rank_biserial=np.nan,
                    n_a=int(a.size), n_b=int(b.size),
                    median_a=np.nan, median_b=np.nan, method="mann_whitney_u")
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    # Rank-biserial correlation: the difference between the probability that a
    # random draw from a exceeds one from b and the reverse.
    rb = 2.0 * u / (a.size * b.size) - 1.0
    return dict(statistic=float(u), p_value=float(p), rank_biserial=float(rb),
                n_a=int(a.size), n_b=int(b.size),
                median_a=float(np.median(a)), median_b=float(np.median(b)),
                method="mann_whitney_u")


def brunner_munzel(a, b) -> dict:
    """Brunner-Munzel test, which does not assume equal variances or shapes.

    Reports the stochastic superiority p_hat, the probability that a random
    draw from b exceeds a random draw from a, with ties counted as half. Unlike
    Mann-Whitney this remains valid when the two distributions have different
    shapes, which matters here because defense carriage can compress the upper
    tail of the plaque score without shifting the median.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if a.size < 2 or b.size < 2:
        return dict(statistic=np.nan, p_value=np.nan, p_hat=np.nan,
                    n_a=int(a.size), n_b=int(b.size), method="brunner_munzel")
    try:
        res = stats.brunnermunzel(a, b, alternative="two-sided")
        stat, p = float(res.statistic), float(res.pvalue)
    except ValueError:
        stat, p = np.nan, np.nan
    # p_hat computed directly so that it is defined even when the test is not.
    gt = 0.0
    chunk = 4096
    for i in range(0, a.size, chunk):
        blk = a[i:i + chunk][:, None]
        gt += float(np.sum(b[None, :] > blk) + 0.5 * np.sum(b[None, :] == blk))
    p_hat = gt / (a.size * b.size)
    return dict(statistic=stat, p_value=p, p_hat=float(p_hat),
                n_a=int(a.size), n_b=int(b.size), method="brunner_munzel")


def wilson_interval(successes: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval for a proportion, stable when the count is 0 or n."""
    if n == 0:
        return (np.nan, np.nan)
    p = successes / n
    d = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / d
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return (float(max(0.0, centre - half)), float(min(1.0, centre + half)))


def cluster_bootstrap_ci(frame: pd.DataFrame, cluster: str, statistic,
                         n_boot: int = 2000, seed: int = 20260101,
                         alpha: float = 0.05) -> tuple[float, float]:
    """Percentile interval from resampling whole clusters with replacement.

    Resampling clusters rather than rows keeps the dependence between pairs
    that share a bacterial genome, which is the unit that actually repeats in
    this design.
    """
    rng = np.random.default_rng(seed)
    groups = [g for _, g in frame.groupby(cluster, observed=True, sort=True)]
    k = len(groups)
    if k < 3:
        return (np.nan, np.nan)
    vals: list[float] = []
    for _ in range(n_boot):
        pick = rng.integers(0, k, size=k)
        sample = pd.concat([groups[i] for i in pick], ignore_index=True)
        v = statistic(sample)
        if v is not None and np.isfinite(v):
            vals.append(float(v))
    if len(vals) < 50:
        return (np.nan, np.nan)
    lo = float(np.quantile(vals, alpha / 2))
    hi = float(np.quantile(vals, 1 - alpha / 2))
    return (lo, hi)
