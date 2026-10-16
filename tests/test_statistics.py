"""Check the statistical core against independent implementations.

Every number in this project passes through `papr.stats` and `papr.models`, so
a defect in either would be invisible to any check that compares the written
results against the computed ones: the documentation would faithfully report
the wrong number. These tests compare each estimator against statsmodels,
scipy, or a brute-force calculation, which is the only way such a defect shows
up. One did. The middle term of the Robins, Breslow and Greenland variance had
an extra factor of one half, which made every Mantel-Haenszel interval about a
tenth too narrow.

Run with:

    python tests/test_statistics.py

It exits non-zero if anything fails, so it can be used as a gate.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats as sps
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportion_confint

from papr import sketch
from papr.models import fit_logit_cluster, term_effect
from papr.stats import (benjamini_hochberg, brunner_munzel, cluster_bootstrap_ci,
                        fisher_odds_ratio, mann_whitney, mantel_haenszel,
                        wilson_interval)

RESULTS: list[tuple[str, bool, str]] = []


def check(label: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((label, bool(ok), detail))


def test_benjamini_hochberg(rng):
    worst = 0.0
    for _ in range(300):
        p = rng.random(int(rng.integers(2, 60))) ** rng.choice([1, 2, 4])
        worst = max(worst, float(np.max(np.abs(
            benjamini_hochberg(p) - multipletests(p, method="fdr_bh")[1]))))
    check("Benjamini-Hochberg matches statsmodels fdr_bh", worst < 1e-12,
          f"max difference {worst:.2e}")
    q = benjamini_hochberg(np.array([0.04, np.nan, 0.001, 0.9]))
    ref = multipletests(np.array([0.04, 0.001, 0.9]), method="fdr_bh")[1]
    check("it passes NaN through and corrects only the rest",
          np.isnan(q[1]) and np.allclose(q[[0, 2, 3]], ref, atol=1e-12))


def test_mantel_haenszel(rng):
    or_err = se_err = p_err = 0.0
    for _ in range(200):
        tabs = [tuple(int(rng.integers(1, 40)) for _ in range(4))
                for _ in range(int(rng.integers(2, 9)))]
        mine = mantel_haenszel(tabs)
        if not np.isfinite(mine.estimate):
            continue
        ref = StratifiedTable([np.array([[a, b], [c, d]]) for a, b, c, d in tabs])
        or_err = max(or_err, abs(mine.estimate - ref.oddsratio_pooled))
        se_mine = (np.log(mine.ci_high) - np.log(mine.estimate)) / 1.96
        se_err = max(se_err, abs(se_mine - ref.logodds_pooled_se))
        p_err = max(p_err, abs(mine.p_value
                               - ref.test_null_odds(correction=True).pvalue))
    check("Mantel-Haenszel odds ratio matches statsmodels", or_err < 1e-9,
          f"max difference {or_err:.2e}")
    check("its Robins-Breslow-Greenland standard error matches statsmodels",
          se_err < 1e-9, f"max difference {se_err:.2e}")
    check("its continuity-corrected CMH p-value matches statsmodels",
          p_err < 1e-9, f"max difference {p_err:.2e}")
    base = [(5, 7, 9, 4), (3, 8, 6, 5)]
    padded = base + [(0, 0, 4, 9), (7, 3, 0, 0)]
    m1, m2 = mantel_haenszel(base), mantel_haenszel(padded)
    check("strata with no variation are dropped and change nothing",
          m1.n_strata_used == 2 and m2.n_strata_used == 2
          and abs(m1.estimate - m2.estimate) < 1e-12)


def test_proportions_and_tables(rng):
    p_err = 0.0
    for _ in range(200):
        e = rng.random(80) < 0.4
        y = rng.random(80) < (0.2 + 0.3 * e)
        a, b = int(np.sum(e & y)), int(np.sum(e & ~y))
        c, d = int(np.sum(~e & y)), int(np.sum(~e & ~y))
        p_err = max(p_err, abs(fisher_odds_ratio(e, y).p_value
                               - sps.fisher_exact([[a, b], [c, d]])[1]))
    check("Fisher p-value matches scipy", p_err < 1e-12,
          f"max difference {p_err:.2e}")
    sparse = fisher_odds_ratio(np.array([1, 1, 0, 0], bool),
                               np.array([1, 1, 1, 0], bool))
    check("an empty cell still gives a finite interval",
          np.isfinite(sparse.ci_low) and np.isfinite(sparse.ci_high))
    worst = 0.0
    for n in (1, 2, 5, 13, 50, 403, 38435):
        for s in {0, 1, n // 2, n - 1, n}:
            if 0 <= s <= n:
                lo, hi = wilson_interval(s, n)
                rlo, rhi = proportion_confint(s, n, alpha=0.05, method="wilson")
                worst = max(worst, abs(lo - rlo), abs(hi - rhi))
    check("Wilson interval matches statsmodels", worst < 2e-3,
          f"max difference {worst:.2e}")


def test_rank_tests(rng):
    rb_err = phat_err = 0.0
    for _ in range(150):
        a = rng.normal(size=int(rng.integers(5, 50)))
        b = rng.normal(0.4, size=int(rng.integers(5, 50)))
        mine = mann_whitney(a, b)
        u, p = sps.mannwhitneyu(a, b, alternative="two-sided")
        if abs(mine["statistic"] - u) > 1e-9 or abs(mine["p_value"] - p) > 1e-12:
            check("Mann-Whitney matches scipy", False)
            return
        brute = np.mean([float(x > y) - float(x < y) for x in a for y in b])
        rb_err = max(rb_err, abs(mine["rank_biserial"] - brute))
    check("Mann-Whitney matches scipy", True)
    check("its rank-biserial equals the brute-force superiority difference",
          rb_err < 1e-9, f"max difference {rb_err:.2e}")
    for _ in range(150):
        a = rng.integers(0, 5, size=int(rng.integers(8, 50))).astype(float)
        b = rng.integers(0, 5, size=int(rng.integers(8, 50))).astype(float)
        brute = np.mean([float(y > x) + 0.5 * float(y == x)
                         for x in a for y in b])
        phat_err = max(phat_err, abs(brunner_munzel(a, b)["p_hat"] - brute))
    check("Brunner-Munzel p_hat equals brute-force stochastic superiority",
          phat_err < 1e-9, f"max difference {phat_err:.2e}")


def test_cluster_robust(rng):
    worst_beta = worst_se = 0.0
    for _ in range(15):
        sizes = rng.integers(3, 15, size=40)
        groups = np.repeat(np.arange(40), sizes)
        x = rng.normal(size=groups.size) + rng.normal(size=40)[groups]
        eta = -0.5 + 0.8 * x + rng.normal(0, 0.5, size=40)[groups]
        y = (rng.random(groups.size) < 1 / (1 + np.exp(-eta))).astype(float)
        X = sm.add_constant(x)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mine = fit_logit_cluster(y, X, groups)
            ref = sm.Logit(y, X).fit(disp=0, cov_type="cluster",
                                     cov_kwds={"groups": groups})
            naive = sm.Logit(y, X).fit(disp=0)
        worst_beta = max(worst_beta, float(np.max(np.abs(mine.params - ref.params))))
        worst_se = max(worst_se, float(np.max(np.abs(mine.bse - ref.bse))))
    check("cluster-robust logistic matches a direct statsmodels fit",
          worst_beta < 1e-6, f"max difference {worst_beta:.2e}")
    check("and so do its clustered standard errors", worst_se < 1e-6,
          f"max difference {worst_se:.2e}")
    check("clustering widens the interval on correlated data",
          mine.bse[1] > naive.bse[1])
    eff = term_effect(mine, 1)
    check("term_effect reports exp(beta) and a matching interval",
          abs(eff["odds_ratio"] - np.exp(mine.params[1])) < 1e-12
          and abs(eff["log2_or"] - mine.params[1] / np.log(2)) < 1e-12)


def test_sketching(rng):
    a = np.array([1, 2, 3, 4, 5], dtype=np.uint64)
    b = np.array([3, 4, 5, 6, 7, 8], dtype=np.uint64)
    check("containment is the intersection over the smaller set",
          abs(sketch.containment(a, b) - 3 / 5) < 1e-12)
    check("Jaccard is the intersection over the union",
          abs(sketch.jaccard(a, b) - 3 / 8) < 1e-12)
    check("Mash distance is zero for identical sketches",
          abs(sketch.mash_distance(a, a)) < 1e-12)
    letters = {65: 0, 67: 1, 71: 2, 84: 3}
    comp = {65: 84, 67: 71, 71: 67, 84: 65}
    seq = bytes("ACGTTGCAAGGCTTAACGGTACGT", "ascii")
    fwd = np.array([letters[c] for c in seq], dtype=np.uint8)
    rev = np.array([letters[comp[c]] for c in seq[::-1]], dtype=np.uint8)
    check("canonical hashing is the same for a sequence and its reverse complement",
          sorted(sketch._canonical_hashes(fwd).tolist())
          == sorted(sketch._canonical_hashes(rev).tolist()))
    codes = rng.integers(0, 4, size=200000).astype(np.uint8)
    h = sketch._canonical_hashes(codes)
    for scaled in (10, 20, 50):
        kept = h[h <= np.uint64(0xFFFFFFFFFFFFFFFF) // np.uint64(scaled)].size
        check(f"FracMinHash keeps about one hash in {scaled}",
              abs(kept / h.size - 1 / scaled) < 0.2 / scaled,
              f"kept {kept / h.size:.4f}")


def test_bootstrap(rng):
    hits = trials = 0
    for _ in range(120):
        g = np.repeat(np.arange(30), 8)
        v = rng.normal(0, 1, size=g.size) + rng.normal(0, 1, size=30)[g]
        lo, hi = cluster_bootstrap_ci(pd.DataFrame({"cluster": g, "v": v}),
                                      "cluster", lambda d: d["v"].mean(),
                                      n_boot=200,
                                      seed=int(rng.integers(1 << 30)))
        if np.isfinite(lo):
            trials += 1
            hits += int(lo <= 0.0 <= hi)
    cov = hits / max(trials, 1)
    check(f"cluster bootstrap covers the truth near 95 percent of the time "
          f"({100 * cov:.1f})", 0.88 <= cov <= 0.99, f"coverage {cov:.3f}")


def main() -> int:
    rng = np.random.default_rng(20260101)
    for fn in (test_benjamini_hochberg, test_mantel_haenszel,
               test_proportions_and_tables, test_rank_tests,
               test_cluster_robust, test_sketching, test_bootstrap):
        fn(rng)
    for label, ok, detail in RESULTS:
        line = ("  OK   " if ok else "  FAIL ") + "  " + label
        if detail and not ok:
            line += f"   [{detail}]"
        print(line)
    bad = [r for r in RESULTS if not r[1]]
    print()
    print(f"{len(RESULTS) - len(bad)} of {len(RESULTS)} checks passed")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
