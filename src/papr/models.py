"""Binomial regression with cluster-robust covariance.

Shared by the population structure control and the counter-defense analysis so
that both use an identical estimator.
"""

from __future__ import annotations

import warnings

import numpy as np
import statsmodels.api as sm


def fit_logit_cluster(y: np.ndarray, x: np.ndarray, groups: np.ndarray):
    """Logistic fit with standard errors clustered on `groups`.

    A Newton solve is used rather than the iteratively reweighted least squares
    path that statsmodels takes for GLM, because the latter calls a
    singular-value least squares routine on the full n by k weighted design and
    exhausts the LAPACK workspace at this sample size. Newton forms the k by k
    Hessian instead, which is small and well conditioned here. If Newton fails
    on a badly scaled design the fit falls back to a quasi-Newton search.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = sm.Logit(y, x)
        try:
            return model.fit(method="newton", maxiter=200, disp=0,
                             cov_type="cluster", cov_kwds={"groups": groups})
        except (np.linalg.LinAlgError, MemoryError, ValueError):
            return model.fit(method="bfgs", maxiter=600, disp=0,
                             cov_type="cluster", cov_kwds={"groups": groups})


def term_effect(res, index: int) -> dict:
    """Odds ratio, interval and p-value for one coefficient of a fitted model."""
    beta = float(res.params[index])
    se = float(res.bse[index])
    return dict(
        odds_ratio=float(np.exp(beta)),
        ci_low=float(np.exp(beta - 1.96 * se)),
        ci_high=float(np.exp(beta + 1.96 * se)),
        log2_or=float(beta / np.log(2.0)),
        coefficient=beta,
        std_error=se,
        p_value=float(res.pvalues[index]),
    )
