"""Search the whole accessory genome for anything associated with resistance.

Every association tested so far has been chosen in advance: 137 catalogued
defense systems, four surface determinants, anti-defense families, prophages.
Together they leave about three quarters of failed infections unexplained, and
no catalogue can be searched for something that is not in it.

With the host genomes validated in script 19, that restriction can be lifted.
This is a k-mer based genome-wide association: every variable k-mer in the
collection is tested against resistance, whether or not anyone has named the
gene it sits in.

Three things make it tractable and, more importantly, honest.

Collapsing. K-mers from the same locus occur in exactly the same strains, so
they share a presence pattern. Collapsing identical patterns reduces millions
of correlated tests to a far smaller number of distinct ones and handles
linkage without an arbitrary pruning rule.

Calibration. Population structure, not sampling noise, is what generates false
positives in a bacterial association study. Each pattern is compared against
simulated traits carrying the same prevalence and the same degree of
phylogenetic confinement, using the machinery built in script 21. Because the
null depends only on those two quantities, it is precomputed on a grid and
looked up.

Conditioning. Host surface type explains more of the variation in lysis than
anything else on the host side, so the association is run twice: on the lysis
rate itself, and on lysis with surface type regressed out.

Outputs:
    data/processed/gwas_patterns.npz
    results/tables/table45_gwas_hits.csv
    results/tables/table46_gwas_calibration.csv
    results/accessory_gwas_report.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from papr import config, paths
from papr.sketch import K, read_raw_sequences, sketch_file
from papr.stats import benjamini_hochberg

GENOME_DIR = paths.RAW / "host_genomes"
PATTERNS = paths.PROCESSED / "gwas_patterns.npz"

# One k-mer in fifty is kept. A 1 kb gene then contributes about twenty k-mers
# and a 10 kb island about two hundred, which is ample to detect a locus while
# keeping the whole collection in memory at once.
GWAS_SCALED = 50

# A k-mer carried by almost every strain or almost none carries no information
# about differences between them.
MIN_CARRIERS = 10

# Grid for the precomputed null, over carrier count and phylogenetic signal.
N_PREVALENCE_BINS = 24
N_SIGNAL_BINS = 8
N_NULL_SIMULATIONS = 4000


def build_presence(strains: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Presence patterns of every variable k-mer, collapsed to unique patterns.

    Returns the packed patterns, how many k-mers share each pattern, and a
    representative k-mer hash per pattern.
    """
    if PATTERNS.exists():
        with np.load(PATTERNS, allow_pickle=False) as z:
            if list(z["strains"]) == strains:
                return z["patterns"], z["weights"], z["representatives"]

    all_hashes, all_idx = [], []
    for i, s in enumerate(strains):
        path = GENOME_DIR / f"{s}.fna.gz"
        h = sketch_file(path, K, GWAS_SCALED)
        all_hashes.append(h)
        all_idx.append(np.full(h.size, i, dtype=np.uint16))
        if (i + 1) % 40 == 0 or i == len(strains) - 1:
            print(f"  sketched {i + 1}/{len(strains)}")
    hashes = np.concatenate(all_hashes)
    idx = np.concatenate(all_idx)
    del all_hashes, all_idx
    print(f"  {hashes.size} k-mer observations, "
          f"{np.unique(hashes).size} distinct")

    order = np.argsort(hashes, kind="stable")
    hashes = hashes[order]
    idx = idx[order]
    del order

    # Group boundaries of identical hashes.
    starts = np.flatnonzero(np.concatenate(([True], hashes[1:] != hashes[:-1])))
    counts = np.diff(np.concatenate((starts, [hashes.size])))
    keep = (counts >= MIN_CARRIERS) & (counts <= len(strains) - MIN_CARRIERS)
    print(f"  {int(keep.sum())} variable k-mers "
          f"(carried by {MIN_CARRIERS} to {len(strains) - MIN_CARRIERS} strains)")

    n_words = (len(strains) + 63) // 64
    sel_starts = starts[keep]
    sel_counts = counts[keep]
    packed = np.zeros((sel_starts.size, n_words), dtype=np.uint64)
    reps = hashes[sel_starts]
    for row, (st, ct) in enumerate(zip(sel_starts, sel_counts)):
        members = idx[st:st + ct]
        words = members >> 6
        bits = (np.uint64(1) << (members & 63).astype(np.uint64))
        np.bitwise_or.at(packed[row], words, bits)
        if (row + 1) % 200000 == 0:
            print(f"    packed {row + 1}/{sel_starts.size}")

    # Collapse k-mers that occur in exactly the same strains.
    uniq, inverse = np.unique(packed, axis=0, return_inverse=True)
    weights = np.bincount(inverse, minlength=uniq.shape[0])
    first = np.zeros(uniq.shape[0], dtype=np.uint64)
    seen = np.zeros(uniq.shape[0], dtype=bool)
    for row, pat in enumerate(inverse):
        if not seen[pat]:
            first[pat] = reps[row]
            seen[pat] = True
    print(f"  {uniq.shape[0]} distinct presence patterns")

    np.savez_compressed(PATTERNS, patterns=uniq, weights=weights,
                        representatives=first, strains=np.array(strains))
    return uniq, weights, first


def unpack(patterns: np.ndarray, n: int) -> np.ndarray:
    """Expand packed bit patterns into a boolean matrix of patterns by strains."""
    out = np.zeros((patterns.shape[0], n), dtype=bool)
    for j in range(n):
        out[:, j] = (patterns[:, j >> 6] >> np.uint64(j & 63)) & np.uint64(1) != 0
    return out


def phylogenetic_r2(carrier: np.ndarray, axes: np.ndarray) -> np.ndarray:
    """Share of carriage variance explained by core-genome position, vectorised."""
    x = np.hstack([np.ones((axes.shape[0], 1)), axes])
    beta, *_ = np.linalg.lstsq(x, carrier.T.astype(float), rcond=None)
    fitted = x @ beta
    resid = carrier.T.astype(float) - fitted
    ss_res = (resid ** 2).sum(axis=0)
    centred = carrier.T.astype(float) - carrier.T.astype(float).mean(axis=0)
    ss_tot = (centred ** 2).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        r2 = 1.0 - ss_res / ss_tot
    return np.clip(np.nan_to_num(r2), 0.0, 1.0)


def null_grid(axes: np.ndarray, lysis: np.ndarray, prevalences: np.ndarray,
              signals: np.ndarray, rng: np.random.Generator) -> dict:
    """Null distributions of the test statistic on a grid of prevalence and signal."""
    n = axes.shape[0]
    grid = {}
    for k in prevalences:
        for r2 in signals:
            b = rng.standard_normal((axes.shape[1], N_NULL_SIMULATIONS))
            sig = axes @ b
            sig -= sig.mean(axis=0, keepdims=True)
            sd = sig.std(axis=0, keepdims=True)
            sd[sd == 0] = 1.0
            sig /= sd
            z = (np.sqrt(r2) * sig
                 + np.sqrt(max(0.0, 1.0 - r2)) * rng.standard_normal((n, N_NULL_SIMULATIONS)))
            cut = np.partition(z, n - int(k), axis=0)[n - int(k)]
            traits = (z >= cut)
            kk = traits.sum(axis=0).astype(float)
            tot = lysis.sum()
            cs = lysis @ traits
            with np.errstate(invalid="ignore", divide="ignore"):
                stat = cs / kk - (tot - cs) / (n - kk)
            grid[(int(k), float(r2))] = np.sort(np.abs(stat[np.isfinite(stat)]))
    return grid


def main() -> int:
    paths.bootstrap()
    val = pd.read_csv(paths.PROCESSED / "host_genome_validation.csv")
    pairs = pd.read_csv(paths.PAIRS, compression="gzip", low_memory=False)
    strain_feats = pd.read_csv(paths.STRAIN_FEATURES, low_memory=False)
    pcoa = pd.read_csv(paths.PCOA, index_col=0)

    per_strain = pairs.groupby("strain")["lytic"].mean()
    strains = sorted(s for s in val.loc[val["validated"], "strain"].astype(str)
                     if s in per_strain.index and s in pcoa.index
                     and (GENOME_DIR / f"{s}.fna.gz").exists())
    print(f"strains with a validated genome and a phenotype: {len(strains)}")
    if len(strains) < 100:
        raise SystemExit("too few strains for an association study")

    lysis = per_strain.loc[strains].to_numpy(dtype=float)
    axes = pcoa.loc[strains].to_numpy(dtype=float)
    axes = (axes - axes.mean(axis=0)) / axes.std(axis=0)
    feats = strain_feats.set_index("strain").loc[strains]

    surface = pd.get_dummies(
        feats[["lps_core", "o_status"]].astype("string").fillna("missing"),
        drop_first=True, dtype=float).to_numpy()
    sx = np.hstack([np.ones((len(strains), 1)), surface])
    beta, *_ = np.linalg.lstsq(sx, lysis, rcond=None)
    lysis_resid = lysis - sx @ beta

    print("building k-mer presence patterns")
    patterns, weights, reps = build_presence(strains)
    carrier = unpack(patterns, len(strains))
    n_pat = carrier.shape[0]
    print(f"testing {n_pat} distinct presence patterns")

    k_counts = carrier.sum(axis=1)
    r2 = phylogenetic_r2(carrier, axes)

    prev_edges = np.unique(np.linspace(MIN_CARRIERS, len(strains) - MIN_CARRIERS,
                                       N_PREVALENCE_BINS).astype(int))
    sig_edges = np.linspace(0.0, 0.7, N_SIGNAL_BINS)
    rng = np.random.default_rng(config.RANDOM_SEED)
    print(f"precomputing null grid: {len(prev_edges)} x {len(sig_edges)} cells")

    results = {}
    for label, y in (("lysis_rate", lysis), ("surface_adjusted", lysis_resid)):
        grid = null_grid(axes, y, prev_edges, sig_edges, rng)
        tot = y.sum()
        cs = carrier @ y
        kk = k_counts.astype(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            stat = cs / kk - (tot - cs) / (len(strains) - kk)
        # Nearest grid cell for each pattern.
        pi = np.abs(k_counts[:, None] - prev_edges[None, :]).argmin(axis=1)
        si = np.abs(r2[:, None] - sig_edges[None, :]).argmin(axis=1)
        pvals = np.ones(n_pat)
        for a in range(len(prev_edges)):
            for b in range(len(sig_edges)):
                m = (pi == a) & (si == b)
                if not m.any():
                    continue
                ref = grid[(int(prev_edges[a]), float(sig_edges[b]))]
                rank = np.searchsorted(ref, np.abs(stat[m]), side="left")
                pvals[m] = (ref.size - rank + 1) / (ref.size + 1)
        results[label] = dict(stat=stat, p=pvals,
                              q=benjamini_hochberg(pvals))
        print(f"  {label}: {int((results[label]['q'] <= config.FDR_ALPHA).sum())} "
              f"patterns at q <= {config.FDR_ALPHA}")

    table = pd.DataFrame(dict(
        pattern=np.arange(n_pat), n_kmers=weights, n_carrier_strains=k_counts,
        phylogenetic_r2=r2, representative_hash=reps.astype(np.int64),
        effect_lysis=results["lysis_rate"]["stat"],
        p_lysis=results["lysis_rate"]["p"], q_lysis=results["lysis_rate"]["q"],
        effect_surface_adjusted=results["surface_adjusted"]["stat"],
        p_surface_adjusted=results["surface_adjusted"]["p"],
        q_surface_adjusted=results["surface_adjusted"]["q"]))
    table["significant"] = ((table["q_lysis"] <= config.FDR_ALPHA)
                            | (table["q_surface_adjusted"] <= config.FDR_ALPHA))
    hits = table[table["significant"]].sort_values("p_surface_adjusted")
    hits.to_csv(paths.TABLES / "table45_gwas_hits.csv", index=False)

    # Calibration: the distribution of p-values should be uniform under the null.
    calib = []
    for label in ("lysis_rate", "surface_adjusted"):
        p = results[label]["p"]
        calib.append(dict(outcome=label, n_patterns=n_pat,
                          median_p=float(np.median(p)),
                          frac_below_05=float((p <= 0.05).mean()),
                          frac_below_01=float((p <= 0.01).mean()),
                          n_significant=int((results[label]["q"] <= config.FDR_ALPHA).sum())))
    calibration = pd.DataFrame(calib)
    calibration.to_csv(paths.TABLES / "table46_gwas_calibration.csv", index=False)

    o: list[str] = []
    add = o.append
    add("Accessory genome association report")
    add("Generated by scripts/25_accessory_genome_gwas.py")
    add("")
    add("SCOPE")
    add(f"  strains with a validated genome and a phenotype : {len(strains)}")
    add(f"  k-mer size / sampling rate : {K} / one in {GWAS_SCALED}")
    add(f"  variable k-mers collapsed into distinct presence patterns : {n_pat}")
    add(f"  k-mers per pattern : median {int(np.median(weights))}, "
        f"max {int(weights.max())}")
    add("  Collapsing identical presence patterns is what makes this tractable "
        "and is also the right treatment of linkage: k-mers from one locus "
        "occur in exactly the same strains and are one test, not thousands.")

    add("")
    add("CALIBRATION")
    add(f"  {'outcome':20s} {'median p':>9s} {'frac p<=0.05':>13s} "
        f"{'frac p<=0.01':>13s} {'significant':>12s}")
    for r in calibration.itertuples(index=False):
        add(f"  {r.outcome:20s} {r.median_p:9.3f} {r.frac_below_05:13.4f} "
            f"{r.frac_below_01:13.4f} {r.n_significant:12d}")
    add("  Under a correctly calibrated null the median p-value is 0.5 and the "
        "fraction below 0.05 is 0.05. Departures measure how far the "
        "phylogenetically matched null still fails to absorb population "
        "structure.")

    add("")
    add("HITS")
    if len(hits) == 0:
        add("  No k-mer presence pattern is associated with resistance at "
            f"q <= {config.FDR_ALPHA} after matching the null on prevalence "
            "and phylogenetic confinement.")
        add("")
        add("  That is a substantive result rather than a failed search. The "
            "accessory genome of this collection contains no locus whose "
            "presence predicts resistance better than an arbitrary trait with "
            "the same distribution across the phylogeny. Whatever accounts for "
            "the three quarters of failed infections left unexplained by "
            "surface type, catalogued defense systems and prophages, it is "
            "not a single widely distributed accessory locus.")
        add("  Three limits bound that statement. The study covers "
            f"{len(strains)} strains, so an effect confined to a handful of "
            "them is invisible. One k-mer in fifty is sampled, so a locus "
            "under a few hundred bases could be missed. And the outcome is a "
            "strain's average across 96 phages, which dilutes anything acting "
            "on a single phage lineage.")
    else:
        add(f"  {len(hits)} patterns significant at q <= {config.FDR_ALPHA}.")
        add("")
        add(f"  {'kmers':>7s} {'strains':>8s} {'phyR2':>7s} {'effect':>8s} "
            f"{'q lysis':>9s} {'q surf-adj':>11s}")
        for r in hits.head(25).itertuples(index=False):
            add(f"  {r.n_kmers:7d} {r.n_carrier_strains:8d} "
                f"{r.phylogenetic_r2:7.3f} {r.effect_surface_adjusted:+8.4f} "
                f"{r.q_lysis:9.2e} {r.q_surface_adjusted:11.2e}")
        add("  'kmers' is how many k-mers share the pattern, a rough proxy for "
            f"locus length: at one k-mer in {GWAS_SCALED}, a pattern of 20 "
            "k-mers corresponds to roughly a kilobase.")

    add("")
    add("WHAT THIS DESIGN CAN AND CANNOT SEE")
    add("  It can see a locus present in a reasonable number of strains whose "
        "presence tracks resistance beyond what its phylogenetic distribution "
        "predicts. That is the standard target of a bacterial association "
        "study and it is tested here without a catalogue.")
    add("  It cannot see a locus that is confined to one lineage, because the "
        "matched null deliberately removes exactly that signal and no "
        "observational design can separate a lineage-confined cause from the "
        "lineage. It cannot see point mutations, since presence of a k-mer is "
        "a presence or absence test. It cannot see anything that is not in an "
        "assembly.")

    text = "\n".join(o)
    (paths.RESULTS / "accessory_gwas_report.txt").write_text(text + "\n", encoding="utf-8")
    print("\n" + text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
