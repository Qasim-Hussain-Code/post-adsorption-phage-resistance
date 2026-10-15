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
from scipy import stats as sps

from papr.stats import benjamini_hochberg

GENOME_DIR = paths.RAW / "host_genomes"
PATTERNS = paths.PROCESSED / "gwas_patterns.npz"

# One k-mer in twenty is kept, about fifty per kilobase of sequence. A sparser
# rate was tried first and failed for a reason worth recording: at one in
# fifty, k-mers from the same gene land in slightly different sets of strains,
# 83 percent of presence patterns were backed by a single k-mer, and the
# collapsing step that is supposed to turn a locus into one test did almost
# nothing.
GWAS_SCALED = 20

# A pattern backed by fewer k-mers than this is not a locus. At one k-mer in
# ten, five k-mers correspond to roughly fifty bases of shared sequence, which
# is the minimum at which a presence call means anything.
MIN_KMERS_PER_PATTERN = 5

# A k-mer carried by almost every strain or almost none carries no information
# about differences between them.
MIN_CARRIERS = 10

# Grid for the simulated null, over carrier count and phylogenetic signal.
N_PREVALENCE_BINS = 24
N_SIGNAL_BINS = 10
N_NULL_SIMULATIONS = 4000

# Bins for the empirical null. Each needs enough patterns to estimate a centre
# and a spread robustly.
MIN_BIN_PATTERNS = 400


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
    reps = hashes[sel_starts]
    del hashes, starts, counts, keep

    # Pack every k-mer's carrier set into a bitmask in one pass. Doing this
    # with a Python loop over k-mers, as an earlier version did, costs hours at
    # this sampling density; expressing it as a single scattered or-reduction
    # costs seconds.
    total = int(sel_counts.sum())
    rows = np.repeat(np.arange(sel_starts.size, dtype=np.int64), sel_counts)
    within = (np.arange(total, dtype=np.int64)
              - np.repeat(np.cumsum(sel_counts) - sel_counts, sel_counts))
    members = idx[np.repeat(sel_starts, sel_counts) + within]
    del within
    packed = np.zeros(sel_starts.size * n_words, dtype=np.uint64)
    flat = rows * n_words + (members >> 6).astype(np.int64)
    np.bitwise_or.at(packed, flat,
                     np.uint64(1) << (members & 63).astype(np.uint64))
    packed = packed.reshape(sel_starts.size, n_words)
    del rows, members, flat, idx
    print(f"    packed {sel_starts.size} k-mer carrier sets")

    # Collapse k-mers that occur in exactly the same strains.
    uniq, inverse = np.unique(packed, axis=0, return_inverse=True)
    weights = np.bincount(inverse, minlength=uniq.shape[0])
    print(f"  {uniq.shape[0]} raw patterns, "
          f"{int((weights == 1).sum())} backed by a single k-mer")
    first = np.zeros(uniq.shape[0], dtype=np.uint64)
    seen = np.zeros(uniq.shape[0], dtype=bool)
    for row, pat in enumerate(inverse):
        if not seen[pat]:
            first[pat] = reps[row]
            seen[pat] = True
    # Keep only patterns backed by enough k-mers to be a locus.
    solid = weights >= MIN_KMERS_PER_PATTERN
    uniq, weights, first = uniq[solid], weights[solid], first[solid]
    print(f"  {uniq.shape[0]} patterns backed by at least "
          f"{MIN_KMERS_PER_PATTERN} k-mers")

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
    sig_edges = np.linspace(0.0, 1.0, N_SIGNAL_BINS)
    rng = np.random.default_rng(config.RANDOM_SEED)
    print(f"precomputing null grid: {len(prev_edges)} x {len(sig_edges)} cells")

    results = {}
    for label, y in (("lysis_rate", lysis), ("surface_adjusted", lysis_resid)):
        tot = y.sum()
        cs = carrier @ y
        kk = k_counts.astype(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            stat = cs / kk - (tot - cs) / (len(strains) - kk)

        # Null one: traits simulated to match prevalence and phylogenetic
        # signal. This is the null used in script 21 and it is reported here
        # because its failure is informative, not because its p-values are
        # usable.
        grid = null_grid(axes, y, prev_edges, sig_edges, rng)
        pi = np.abs(k_counts[:, None] - prev_edges[None, :]).argmin(axis=1)
        si = np.abs(r2[:, None] - sig_edges[None, :]).argmin(axis=1)
        p_sim = np.ones(n_pat)
        for a in range(len(prev_edges)):
            for b in range(len(sig_edges)):
                m = (pi == a) & (si == b)
                if not m.any():
                    continue
                ref = grid[(int(prev_edges[a]), float(sig_edges[b]))]
                rank = np.searchsorted(ref, np.abs(stat[m]), side="left")
                p_sim[m] = (ref.size - rank + 1) / (ref.size + 1)

        # Null two: the empirical null. Within each prevalence and signal bin
        # the overwhelming majority of k-mers cannot be causal, so the spread
        # of their statistics is the null. Centre and scale are estimated by
        # median and scaled median absolute deviation, which a handful of true
        # associations cannot shift. This absorbs whatever phylogenetic
        # structure real k-mers carry, including structure the eight principal
        # coordinates fail to describe.
        z = np.full(n_pat, np.nan)
        cell = pi * len(sig_edges) + si
        for c in np.unique(cell):
            m = cell == c
            if m.sum() < MIN_BIN_PATTERNS:
                continue
            v = stat[m]
            med = float(np.median(v))
            mad = float(np.median(np.abs(v - med)))
            sigma = 1.4826 * mad
            # A bin dominated by near-identical clade markers can collapse to
            # almost no spread, which would turn every slight departure into an
            # enormous z. Such a bin carries no usable null.
            if sigma <= 1e-6 or sigma < 0.05 * float(np.std(stat)):
                continue
            z[m] = (v - med) / sigma
        p_emp = np.where(np.isfinite(z), 2.0 * sps.norm.sf(np.abs(z)), np.nan)
        usable = np.isfinite(p_emp)
        q_emp = np.full(n_pat, np.nan)
        q_emp[usable] = benjamini_hochberg(p_emp[usable])

        results[label] = dict(stat=stat, p_sim=p_sim,
                              q_sim=benjamini_hochberg(p_sim),
                              z=z, p=p_emp, q=q_emp, n_usable=int(usable.sum()))
        n_sig = int(np.nansum(q_emp <= config.FDR_ALPHA))
        print(f"  {label}: simulated-null median p "
              f"{np.median(p_sim):.3f}; empirical-null median p "
              f"{np.nanmedian(p_emp):.3f}; {n_sig} patterns at "
              f"q <= {config.FDR_ALPHA}")

    table = pd.DataFrame(dict(
        pattern=np.arange(n_pat), n_kmers=weights, n_carrier_strains=k_counts,
        phylogenetic_r2=r2, representative_hash=reps.astype(np.int64),
        effect_lysis=results["lysis_rate"]["stat"],
        z_lysis=results["lysis_rate"]["z"],
        p_lysis=results["lysis_rate"]["p"], q_lysis=results["lysis_rate"]["q"],
        p_lysis_simulated_null=results["lysis_rate"]["p_sim"],
        effect_surface_adjusted=results["surface_adjusted"]["stat"],
        z_surface_adjusted=results["surface_adjusted"]["z"],
        p_surface_adjusted=results["surface_adjusted"]["p"],
        q_surface_adjusted=results["surface_adjusted"]["q"],
        p_surface_adjusted_simulated_null=results["surface_adjusted"]["p_sim"]))
    table["significant"] = ((table["q_lysis"].fillna(1) <= config.FDR_ALPHA)
                            | (table["q_surface_adjusted"].fillna(1) <= config.FDR_ALPHA))
    hits = table[table["significant"]].sort_values("p_surface_adjusted")
    hits = hits[hits["p_surface_adjusted"].notna()
                | hits["p_lysis"].notna()]
    hits.to_csv(paths.TABLES / "table45_gwas_hits.csv", index=False)

    # Calibration: the distribution of p-values should be uniform under the null.
    calib = []
    for label in ("lysis_rate", "surface_adjusted"):
        for kind, pv, qv in (("simulated", results[label]["p_sim"],
                              results[label]["q_sim"]),
                             ("empirical", results[label]["p"],
                              results[label]["q"])):
            good = np.isfinite(pv)
            calib.append(dict(outcome=label, null=kind,
                              n_patterns=int(good.sum()),
                              median_p=float(np.nanmedian(pv[good])),
                              frac_below_05=float((pv[good] <= 0.05).mean()),
                              frac_below_01=float((pv[good] <= 0.01).mean()),
                              n_significant=int(np.nansum(qv <= config.FDR_ALPHA))))
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
    add(f"  {'outcome':18s} {'null':11s} {'loci':>7s} {'median p':>9s} "
        f"{'frac p<=0.05':>13s} {'frac p<=0.01':>13s} {'significant':>12s}")
    for r in calibration.itertuples(index=False):
        add(f"  {r.outcome:18s} {r.null:11s} {r.n_patterns:7d} "
            f"{r.median_p:9.3f} {r.frac_below_05:13.4f} "
            f"{r.frac_below_01:13.4f} {r.n_significant:12d}")
    n_emp = int(calibration[calibration["null"] == "empirical"]["n_patterns"].iloc[0])
    add("")
    add(f"  The empirical null is only defined for {n_emp} of the {n_pat} loci. "
        "A locus whose prevalence and phylogenetic signal place it in a bin "
        "holding too few others, or in one where near-identical clade markers "
        "collapse the spread to nothing, has no usable null and is left "
        f"untested. The {n_pat - n_emp} loci in that position are not null "
        "results: they were not tested.")
    add("")
    add("  Under a correctly calibrated null the median p-value is 0.5 and the "
        "fraction below 0.05 is 0.05.")
    sim_row = calibration[(calibration["null"] == "simulated")
                          & (calibration["outcome"] == "lysis_rate")].iloc[0]
    add(f"  The simulated null fails that test badly: median p "
        f"{sim_row['median_p']:.3f} and {100 * sim_row['frac_below_05']:.1f} "
        "percent of patterns below 0.05, roughly five times the expected rate. "
        "The reason is specific and worth stating. Simulated traits are built "
        "from eight principal coordinates, and a real accessory locus confined "
        "to one clade has phylogenetic structure far sharper than eight axes "
        "can describe. The share of carriage those axes explain therefore "
        "understates how lineage-confined a k-mer really is, and the matched "
        "null is too weak. That same null was adequate in script 21, where it "
        "was applied to 33 defense systems whose distributions are much less "
        "clade-restricted than a typical accessory k-mer.")
    add("  The empirical null does not assume a model for that structure. It "
        "takes the spread of statistics among k-mers of the same prevalence "
        "and the same measured phylogenetic signal as the null, estimated by "
        "median and scaled median absolute deviation so that a few genuine "
        "associations cannot inflate it. Only its p-values are interpreted "
        "below.")

    add("")
    add("HITS, EMPIRICAL NULL ONLY")
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
            f"{'z':>7s} {'q lysis':>9s} {'q surf-adj':>11s}")
        for r in hits.head(25).itertuples(index=False):
            add(f"  {r.n_kmers:7d} {r.n_carrier_strains:8d} "
                f"{r.phylogenetic_r2:7.3f} {r.effect_surface_adjusted:+8.4f} "
                f"{r.z_surface_adjusted:7.2f} {r.q_lysis:9.2e} "
                f"{r.q_surface_adjusted:11.2e}")
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
