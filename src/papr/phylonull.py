"""Phylogenetically matched null traits.

A defense system is not scattered at random across a collection: it sits in
defense islands that are inherited vertically, so its carriers are relatives.
Any test that treats carriage as independent of ancestry will therefore find
associations that belong to the lineage rather than the system. The null used
throughout this project answers that by simulating traits with the same
prevalence and the same degree of phylogenetic confinement as the real one,
and asking whether the real system does anything a matched arbitrary trait
would not.

Shared by the phylogenetic null tests and the phage-resolved analysis so that
both generate the null the same way.
"""

from __future__ import annotations

import numpy as np


def phylogenetic_r2(carrier: np.ndarray, axes: np.ndarray) -> float:
    """Share of carriage variance explained by core-genome position."""
    x = np.hstack([np.ones((axes.shape[0], 1)), axes])
    beta, *_ = np.linalg.lstsq(x, carrier.astype(float), rcond=None)
    fitted = x @ beta
    ss_res = float(np.sum((carrier - fitted) ** 2))
    ss_tot = float(np.sum((carrier - carrier.mean()) ** 2))
    return 0.0 if ss_tot <= 0 else max(0.0, min(1.0, 1.0 - ss_res / ss_tot))


def simulate_traits(axes: np.ndarray, n_carriers: int, r2: float,
                    n_sim: int, rng: np.random.Generator) -> np.ndarray:
    """Null carriage vectors with matched prevalence and phylogenetic signal.

    A random direction through the principal coordinates supplies the
    structured part, white noise the rest, mixed so that the structured part
    accounts for the requested share of variance. Thresholding at the observed
    number of carriers fixes prevalence exactly.
    """
    n = axes.shape[0]
    b = rng.standard_normal((axes.shape[1], n_sim))
    signal = axes @ b
    signal -= signal.mean(axis=0, keepdims=True)
    sd = signal.std(axis=0, keepdims=True)
    sd[sd == 0] = 1.0
    signal /= sd
    noise = rng.standard_normal((n, n_sim))
    z = np.sqrt(r2) * signal + np.sqrt(max(0.0, 1.0 - r2)) * noise
    # Top n_carriers entries of each column become carriers.
    cut = np.partition(z, n - n_carriers, axis=0)[n - n_carriers]
    return (z >= cut).astype(np.int8)


def cluster_block_traits(carrier: np.ndarray, labels: np.ndarray,
                         n_sim: int, rng: np.random.Generator) -> np.ndarray:
    """Null traits that keep the trait's clonal distribution and move it.

    Matching prevalence and the variance explained by a handful of principal
    coordinates is not enough for a defense system. Systems sit in islands
    inherited by near-identical relatives, which is structure finer than those
    coordinates can carry, and a null without it is too easy to beat: carriage
    that tracks clonal lineages more tightly than the null does will score
    against a phenotype that is itself clonally structured.

    This null preserves the trait exactly as a clonal object. It takes the
    number of carriers in each clonal cluster and permutes those counts among
    clusters of the same size, so prevalence and the whole within-cluster
    carriage pattern are unchanged while which lineages carry the trait is
    randomised. The question it poses is the right one: does this system do
    more than a trait of identical prevalence and identical clonal shape that
    happens to sit on other lineages?
    """
    uniq, inv = np.unique(labels, return_inverse=True)
    sizes = np.bincount(inv)
    counts = np.bincount(inv, weights=carrier.astype(float)).round().astype(int)
    members = [np.flatnonzero(inv == c) for c in range(uniq.size)]
    # Clusters are only interchangeable with others of their own size.
    by_size: dict[int, np.ndarray] = {}
    for size in np.unique(sizes):
        by_size[int(size)] = np.flatnonzero(sizes == size)

    out = np.zeros((carrier.size, n_sim), dtype=np.int8)
    for size, idx in by_size.items():
        # One independent permutation of this size class per simulation.
        order = rng.random((n_sim, idx.size)).argsort(axis=1)
        drawn = counts[idx][order]
        if size == 1:
            out[np.concatenate([members[c] for c in idx]), :] = drawn.T
            continue
        for position, c in enumerate(idx):
            k = drawn[:, position]
            # Rank a fresh random draw within the cluster and take the lowest
            # k of them, which picks k of its members uniformly per simulation.
            ranks = rng.random((n_sim, size)).argsort(axis=1)
            out[members[c], :] = (ranks < k[:, None]).T
    return out
