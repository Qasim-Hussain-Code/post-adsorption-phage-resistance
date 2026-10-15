"""Is the central negative result an artifact of averaging across 96 phages?

Every association in this study so far asks whether a defense system protects
on average across the whole panel. Defense systems are characteristically
narrow in spectrum, and a system that blocked five phages of ninety-six would
be diluted by that average to something no test here would detect. If that is
what has happened, the central negative result is an artifact of aggregation
rather than a finding. The question is worth settling rather than assuming,
and it is settleable: the panel is 96 phages whose relationships are known.

The test is a Rao score test for adding a system to a logistic model of lysis
fitted separately for each phage, with core-genome principal coordinates as
covariates, so population structure is controlled as it is elsewhere. The
score form matters for tractability: fitting the covariate model once per
phage lets thousands of candidate traits be tested against it in closed form,
which is what makes a matched null affordable at this scale.

Significance is judged against phylogenetically matched null traits put
through identical machinery rather than against a normal reference. That
choice is not stylistic. A normal reference has already failed twice in this
project wherever lineage-confined traits were involved, and the calibration of
the matched null is reported below rather than assumed.

Three questions are asked of every system:

  does it act on more phages than a matched arbitrary trait would,
  is the set of phages it acts on more closely related than chance, and
  is its best phage much stronger than its average.

The last distinguishes a genuine narrow spectrum from the regression to the
mean that guarantees some phage will look like the best of ninety-six.

Outputs:
    results/tables/table49_phage_resolved_effects.csv
    results/tables/table50_spectrum_tests.csv
    results/tables/table51_spectrum_coherence.csv
    results/tables/table52_phage_cluster_sweep.csv
    results/tables/table53_island_adjusted.csv
    results/phage_resolved_report.txt
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
import statsmodels.api as sm
from Bio import Phylo

from papr import config, paths
from papr.phylonull import (cluster_block_traits, phylogenetic_r2,
                            simulate_traits)
from papr.stats import benjamini_hochberg

# Matched null traits per system. Each one is carried through all 96 per-phage
# tests, so this sets the resolution of the empirical p-values: the smallest
# attainable p is one over this plus one. It has to be small enough that a
# system at the floor still clears a Benjamini-Hochberg correction over the
# systems tested, which 200 was not: it floored p at 0.005 and left the whole
# verdict resting on where a tie happened to fall.
N_NULL_TRAITS = 1000

# Near-identical phages are not independent tests. The panel contains pairs at
# zero distance, and at this threshold its 96 members collapse to 56 groups,
# the largest holding 9. A system associated with one member of such a group
# is associated with all of them for free, so counting phages would make one
# lineage look like a broad spectrum. Distinct groups are counted instead,
# which is the same treatment given to clonal strains and to identical k-mer
# presence patterns elsewhere in this project.
PHAGE_CLUSTER_THRESHOLD = 0.1
PHAGE_CLUSTER_SWEEP = (0.05, 0.1, 0.2)

# Two-sided nominal 0.01. Used only to define "acts on this phage" when
# counting; the counts themselves are judged against the matched null, so this
# threshold does not have to be exactly calibrated to mean anything.
Z_THRESHOLD = 2.5758293035489004

# A phage tested on too few strains cannot support a per-phage model.
MIN_PHAGE_OBSERVATIONS = 200

# Draws for the test of whether a spectrum is more closely related than chance.
N_COHERENCE_DRAWS = 20000

# Co-occurring systems conditioned on when checking whether a surviving system
# is carrying the credit for a neighbour on the same defense island.
N_ISLAND_PARTNERS = 4


def score_z(y: np.ndarray, cov: np.ndarray, traits: np.ndarray) -> np.ndarray:
    """Rao score statistic for adding each column of `traits` to y ~ cov.

    The covariate model is fitted once; every trait is then scored against it
    in closed form. Writing the variance as the trait's weighted sum of
    squares minus its projection on the covariates avoids forming anything of
    the size of the full design, which is what allows thousands of traits to
    be tested per phage.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            fit = sm.Logit(y, cov).fit(method="newton", maxiter=200, disp=0)
        except (np.linalg.LinAlgError, ValueError):
            fit = sm.Logit(y, cov).fit(method="bfgs", maxiter=600, disp=0)
    mu = np.clip(fit.predict(cov), 1e-9, 1 - 1e-9)
    w = mu * (1.0 - mu)

    u = traits.T @ (y - mu)
    info = cov.T @ (w[:, None] * cov)
    # Weighting the covariates rather than the traits keeps the only large
    # intermediate out of memory: the trait matrix is tens of thousands of
    # columns wide and the covariate matrix is nine.
    cross = (cov * w[:, None]).T @ traits
    # Traits are 0/1, so the weighted sum of squares is just the weighted sum.
    v = traits.T @ w - np.einsum("ij,ij->j", cross, np.linalg.solve(info, cross))
    return u / np.sqrt(np.maximum(v, 1e-12))


def patristic_distances(path: Path, names: list[str]) -> np.ndarray:
    """Pairwise distances along the phage tree, in the given tip order."""
    tree = Phylo.read(str(path), "newick")
    tips = {t.name: t for t in tree.get_terminals()}
    n = len(names)
    d = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d[i, j] = d[j, i] = tree.distance(tips[names[i]], tips[names[j]])
    return d


def phage_clusters(dist: np.ndarray, threshold: float) -> np.ndarray:
    """Single-linkage groups of near-identical phages, as integer labels."""
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    if dist.shape[0] < 2:
        return np.zeros(dist.shape[0], dtype=int)
    link = linkage(squareform(dist, checks=False), method="single")
    return fcluster(link, threshold, criterion="distance")


def count_clusters_hit(mask: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Distinct phage groups hit, for each row of a boolean hit matrix."""
    order = np.argsort(labels, kind="mergesort")
    bounds = np.flatnonzero(np.diff(labels[order])) + 1
    counts = np.zeros(mask.shape[0], dtype=int)
    for block in np.split(order, bounds):
        counts += mask[:, block].any(axis=1)
    return counts


def main() -> int:
    paths.bootstrap()
    rng = np.random.default_rng(config.RANDOM_SEED)

    pairs = pd.read_csv(paths.PAIRS, low_memory=False)
    feats = pd.read_csv(paths.STRAIN_FEATURES, low_memory=False).set_index("strain")
    pcoa = pd.read_csv(paths.PCOA, index_col=0)

    strains = [s for s in feats.index if s in pcoa.index]
    phages = sorted(pairs["phage"].unique())
    sindex = {s: i for i, s in enumerate(strains)}
    pindex = {p: j for j, p in enumerate(phages)}

    # Lysis matrix, with missing entries left as NaN rather than filled.
    y_mat = np.full((len(strains), len(phages)), np.nan)
    sub = pairs[pairs["strain"].isin(sindex)]
    y_mat[[sindex[s] for s in sub["strain"]],
          [pindex[p] for p in sub["phage"]]] = sub["lytic"].to_numpy(dtype=float)

    axes = pcoa.loc[strains].to_numpy(dtype=float)
    axes = (axes - axes.mean(axis=0)) / axes.std(axis=0)
    cov = np.hstack([np.ones((len(strains), 1)), axes])

    ds_cols = [c for c in feats.columns if c.startswith("ds_")]
    carriage = feats.loc[strains, ds_cols].to_numpy(dtype=float)
    keep = [i for i, c in enumerate(ds_cols)
            if carriage[:, i].sum() >= config.MIN_CARRIER_STRAINS
            and carriage[:, i].sum() <= len(strains) - config.MIN_CARRIER_STRAINS]
    systems = [ds_cols[i][3:] for i in keep]
    real = carriage[:, keep]
    print(f"strains {len(strains)}, phages {len(phages)}, "
          f"systems {len(systems)}")

    # Two nulls, because neither alone is adequate and they fail in opposite
    # directions. The coordinate null matches the share of carriage explained
    # by the principal coordinates exactly but reproduces only part of the
    # clonal structure a defense system really has; the clonal null matches
    # that structure exactly but carries less broad-scale signal. A system is
    # credited below only if it beats both.
    clusters = pd.read_csv(paths.PROCESSED / "clonal_clusters.csv", index_col=0)
    labels_clonal = clusters.loc[strains].iloc[:, 0].to_numpy()
    r2 = np.array([phylogenetic_r2(real[:, k], axes) for k in range(len(systems))])
    coord_blocks, clonal_blocks = [], []
    for k in range(len(systems)):
        coord_blocks.append(simulate_traits(axes, int(real[:, k].sum()),
                                            float(r2[k]), N_NULL_TRAITS,
                                            rng).astype(float))
        clonal_blocks.append(cluster_block_traits(real[:, k], labels_clonal,
                                                  N_NULL_TRAITS,
                                                  rng).astype(float))
    traits = np.hstack([real] + coord_blocks + clonal_blocks)
    print(f"traits scored per phage: {traits.shape[1]} "
          f"({len(systems)} real, {len(systems) * N_NULL_TRAITS} coordinate "
          f"null, {len(systems) * N_NULL_TRAITS} clonal null)")

    # Score every trait against every phage.
    z = np.full((traits.shape[1], len(phages)), np.nan)
    usable = []
    for j, p in enumerate(phages):
        obs = np.isfinite(y_mat[:, j])
        yv = y_mat[obs, j]
        if obs.sum() < MIN_PHAGE_OBSERVATIONS or yv.sum() < 5 or yv.sum() > obs.sum() - 5:
            continue
        # Most phages were tested on every strain, and the trait matrix is
        # large enough that an avoidable copy of it is worth avoiding.
        if obs.all():
            z[:, j] = score_z(yv, cov, traits)
        else:
            z[:, j] = score_z(yv, cov[obs], traits[obs])
        usable.append(j)
    print(f"phages with a usable per-phage model: {len(usable)} of {len(phages)}")

    zu = z[:, usable]
    phage_names = [phages[j] for j in usable]
    n_real = len(systems)
    span = n_real * N_NULL_TRAITS
    coord_z = zu[n_real:n_real + span].reshape(n_real, N_NULL_TRAITS, len(usable))
    clonal_z = zu[n_real + span:].reshape(n_real, N_NULL_TRAITS, len(usable))
    nulls = {"coordinate": coord_z, "clonal": clonal_z}

    # Calibration of the score test against each null, which is the only
    # reason the counts below can be read as anything.
    calib = {}
    for name, arr in nulls.items():
        flat = arr.reshape(-1)
        flat = flat[np.isfinite(flat)]
        calib[name] = dict(
            median_abs_z=float(np.median(np.abs(flat))),
            frac_beyond_threshold=float((np.abs(flat) > Z_THRESHOLD).mean()),
            sd=float(flat.std()))
        print(f"{name} null calibration: "
              f"{100 * calib[name]['frac_beyond_threshold']:.2f} percent "
              f"beyond |z| {Z_THRESHOLD:.2f} (1.00 expected), "
              f"sd {calib[name]['sd']:.3f}")

    # Per system and phage, against each null.
    rows = []
    for k, s in enumerate(systems):
        for m, p in enumerate(phage_names):
            row = dict(system=s, phage=p, score_z=float(zu[k, m]))
            for name, arr in nulls.items():
                beyond = float((np.abs(arr[k, :, m]) > abs(zu[k, m])).mean())
                row[f"p_{name}"] = max(beyond, 1.0 / (N_NULL_TRAITS + 1))
                row[f"null_sd_{name}"] = float(arr[k, :, m].std())
            rows.append(row)
    per_phage = pd.DataFrame(rows)
    for name in nulls:
        per_phage[f"q_{name}"] = benjamini_hochberg(
            per_phage[f"p_{name}"].to_numpy())
    per_phage["direction"] = np.where(per_phage["score_z"] < 0,
                                      "less lysis", "more lysis")
    per_phage.to_csv(paths.TABLES / "table49_phage_resolved_effects.csv",
                     index=False)

    # Near-identical phages are one test, not several, so the unit counted is
    # the distinct phage group rather than the phage.
    dist = patristic_distances(paths.PHAGE_TREE, phage_names)
    labels = phage_clusters(dist, PHAGE_CLUSTER_THRESHOLD)
    print(f"phage groups at distance {PHAGE_CLUSTER_THRESHOLD}: "
          f"{labels.max()} from {len(phage_names)} phages")

    # Does a system act on more phage groups than a matched trait would?
    spectrum = []
    for k, s in enumerate(systems):
        obs = zu[k]
        obs_prot = (obs < -Z_THRESHOLD)[None, :]
        obs_more = (obs > Z_THRESHOLD)[None, :]
        c_prot = int(count_clusters_hit(obs_prot, labels)[0])
        c_more = int(count_clusters_hit(obs_more, labels)[0])
        rec = dict(
            system=s, n_carriers=int(real[:, k].sum()),
            phylogenetic_r2=float(r2[k]),
            n_phages_protective=int(obs_prot.sum()),
            n_groups_protective=c_prot,
            n_phages_increased=int(obs_more.sum()),
            n_groups_increased=c_more,
            best_phage=phage_names[int(np.argmin(obs))],
            best_z=float(obs.min()),
            mean_z=float(obs.mean()))
        for name, arr in nulls.items():
            nul = arr[k]
            null_prot = count_clusters_hit(nul < -Z_THRESHOLD, labels)
            null_more = count_clusters_hit(nul > Z_THRESHOLD, labels)
            rec[f"null_median_groups_protective_{name}"] = float(
                np.median(null_prot))
            rec[f"p_protective_{name}"] = float(
                (1 + (null_prot >= c_prot).sum()) / (1 + N_NULL_TRAITS))
            rec[f"p_increased_{name}"] = float(
                (1 + (null_more >= c_more).sum()) / (1 + N_NULL_TRAITS))
            rec[f"p_best_{name}"] = float(
                (1 + (nul.min(axis=1) <= obs.min()).sum()) / (1 + N_NULL_TRAITS))
        # A system is credited only where neither null explains it.
        rec["p_protective"] = max(rec["p_protective_coordinate"],
                                  rec["p_protective_clonal"])
        rec["p_increased"] = max(rec["p_increased_coordinate"],
                                 rec["p_increased_clonal"])
        rec["p_best"] = max(rec["p_best_coordinate"], rec["p_best_clonal"])
        spectrum.append(rec)
    spec = pd.DataFrame(spectrum)
    for col in ("p_protective", "p_increased", "p_best"):
        spec[f"q_{col[2:]}"] = benjamini_hochberg(spec[col].to_numpy())
    spec = spec.sort_values("p_protective")
    spec.to_csv(paths.TABLES / "table50_spectrum_tests.csv", index=False)

    # Defense systems travel together on islands, so a system that survives
    # both nulls may be carrying the credit for a neighbour. For each
    # survivor, the systems most enriched among its carriers are added to the
    # covariates and the whole test is run again. Nothing here can prove which
    # member of an island is responsible, but a system whose signal vanishes
    # once its neighbours are conditioned on has not earned the attribution.
    survivors = spec[spec["q_protective"] <= config.FDR_ALPHA]["system"].tolist()
    adjusted = []
    for s in survivors:
        k = systems.index(s)
        v = real[:, k]
        enrich = (real[v == 1].mean(axis=0) - real[v == 0].mean(axis=0))
        partners = [i for i in np.argsort(-enrich)[:N_ISLAND_PARTNERS + 1]
                    if i != k][:N_ISLAND_PARTNERS]
        cov2 = np.hstack([cov, real[:, partners]])
        block = np.hstack([v[:, None], coord_blocks[k], clonal_blocks[k]])
        z2 = np.full((block.shape[1], len(usable)), np.nan)
        for m, j in enumerate(usable):
            ob = np.isfinite(y_mat[:, j])
            z2[:, m] = score_z(y_mat[ob, j], cov2[ob], block[ob])
        got = int(count_clusters_hit((z2[0] < -Z_THRESHOLD)[None, :], labels)[0])
        ps = {}
        for name, lo in (("coordinate", 1),
                         ("clonal", 1 + N_NULL_TRAITS)):
            nl = count_clusters_hit(
                z2[lo:lo + N_NULL_TRAITS] < -Z_THRESHOLD, labels)
            ps[name] = float((1 + (nl >= got).sum()) / (1 + N_NULL_TRAITS))
        adjusted.append(dict(
            system=s, partners=",".join(systems[i] for i in partners),
            n_groups_unadjusted=int(spec.loc[spec.system == s,
                                             "n_groups_protective"].iloc[0]),
            n_groups_adjusted=got,
            p_coordinate=ps["coordinate"], p_clonal=ps["clonal"],
            p_adjusted=max(ps.values()),
            best_z_adjusted=float(z2[0].min())))
    adj = pd.DataFrame(adjusted)
    adj.to_csv(paths.TABLES / "table53_island_adjusted.csv", index=False)

    # How much of the apparent breadth survives collapsing at other
    # thresholds? The headline count should not rest on one cut.
    sweep = []
    for th in PHAGE_CLUSTER_SWEEP:
        lab = phage_clusters(dist, th)
        for k, s in enumerate(systems):
            hit = (zu[k] < -Z_THRESHOLD)[None, :]
            got = int(count_clusters_hit(hit, lab)[0])
            ps = {}
            for name, arr in nulls.items():
                nul = count_clusters_hit(arr[k] < -Z_THRESHOLD, lab)
                ps[name] = float((1 + (nul >= got).sum()) / (1 + N_NULL_TRAITS))
            sweep.append(dict(threshold=th, n_groups_in_panel=int(lab.max()),
                              system=s, n_groups_protective=got,
                              p_coordinate=ps["coordinate"],
                              p_clonal=ps["clonal"],
                              p_protective=max(ps.values())))
    sweep_df = pd.DataFrame(sweep)
    sweep_df.to_csv(paths.TABLES / "table52_phage_cluster_sweep.csv", index=False)

    # Is a system's spectrum more closely related than chance?
    coherence = []
    for k, s in enumerate(systems):
        hit = np.flatnonzero(zu[k] < -Z_THRESHOLD)
        if hit.size < 3:
            continue
        sub_d = dist[np.ix_(hit, hit)]
        observed = float(sub_d[np.triu_indices(hit.size, k=1)].mean())
        draws = np.empty(N_COHERENCE_DRAWS)
        for b in range(N_COHERENCE_DRAWS):
            pick = rng.choice(len(phage_names), size=hit.size, replace=False)
            dd = dist[np.ix_(pick, pick)]
            draws[b] = dd[np.triu_indices(hit.size, k=1)].mean()
        coherence.append(dict(
            system=s, n_phages=int(hit.size),
            n_groups=int(np.unique(labels[hit]).size),
            mean_pairwise_distance=observed,
            random_mean=float(draws.mean()),
            p_more_related=float((1 + (draws <= observed).sum())
                                 / (1 + N_COHERENCE_DRAWS)),
            phages=",".join(phage_names[i] for i in hit[:12])))
    coh = pd.DataFrame(coherence)
    if not coh.empty:
        coh["q_more_related"] = benjamini_hochberg(coh["p_more_related"].to_numpy())
        coh = coh.sort_values("p_more_related")
    coh.to_csv(paths.TABLES / "table51_spectrum_coherence.csv", index=False)

    # Report.
    o: list[str] = []
    add = o.append
    add("Phage-resolved defense effects report")
    add("Generated by scripts/28_phage_resolved_effects.py")
    add("")
    add("THE QUESTION")
    add("  Every association tested elsewhere in this project is an average "
        f"over {len(phage_names)} phages. Defense systems are characteristically "
        "narrow in spectrum, so a system blocking a handful of phages would be "
        "diluted by that average into something no test here would see. If "
        "that is what has happened, the central negative result is an artifact "
        "of aggregation. This script tests each system against each phage "
        "separately and asks whether the averages have been hiding anything.")
    add("")
    add("SCOPE")
    add(f"  strains with core-genome coordinates : {len(strains)}")
    add(f"  phages with a usable per-phage model : {len(phage_names)} of {len(phages)}")
    add(f"  systems tested                       : {len(systems)}")
    add(f"  tests                                : {len(systems) * len(phage_names)}")
    add(f"  matched null traits per system       : {N_NULL_TRAITS}")
    add("")
    add("TWO NULLS, BECAUSE ONE IS NOT ENOUGH")
    add("  A null that matches prevalence and the share of carriage explained "
        "by the principal coordinates is not adequate here, and the reason is "
        "measurable. Defense systems sit in islands inherited by near-"
        "identical relatives. Taking the fraction of same-cluster strain "
        "pairs that agree on carriage, the real systems sit at a median of "
        "0.911 while coordinate-matched null traits reach 0.811, and 24 of "
        "the 33 systems exceed their own null on that measure. Since lysis is "
        "itself clonally structured, a trait that tracks lineages more "
        "tightly than the null does will score against it for reasons that "
        "have nothing to do with defense.")
    add("  The second null fixes exactly that. It takes the number of "
        "carriers in each clonal cluster and permutes those counts among "
        "clusters of the same size, so prevalence and the entire "
        "within-cluster carriage pattern are preserved while the lineages "
        "carrying the trait are randomised. It reproduces the observed clonal "
        "concordance exactly, to four decimal places.")
    add("  Neither null dominates. The coordinate null matches broad-scale "
        "structure and under-matches clonal structure; the clonal null does "
        "the reverse, carrying about half the broad-scale signal of the real "
        "traits. Every p-value reported below is therefore the worse of the "
        "two, so a system is credited only where neither null accounts for "
        "it.")
    add("")
    add("CALIBRATION OF THE SCORE TEST")
    for name in nulls:
        add(f"  {name + ' null':18s}: "
            f"{100 * calib[name]['frac_beyond_threshold']:5.2f} percent beyond "
            f"|z| {Z_THRESHOLD:.2f} (1.00 expected), "
            f"sd {calib[name]['sd']:.3f} (1.000 expected)")
    worst = max(c["frac_beyond_threshold"] for c in calib.values())
    if worst > 0.02:
        add("  At least one null is anti-conservative at this threshold, so "
            "the raw counts below are not read as nominal rates. Every count "
            "is judged against the null distribution of the same count rather "
            "than against the threshold, which is what makes it readable.")
    else:
        add("  Both are close to nominal, so the per-phage counts below can be "
            "read directly.")
    add("")
    add("PHAGE GROUPS, NOT PHAGES")
    add(f"  The panel of {len(phage_names)} phages contains near-identical "
        f"members, with some pairs at zero distance. At a patristic distance "
        f"of {PHAGE_CLUSTER_THRESHOLD} it collapses to {labels.max()} "
        "groups. A system associated with one member of a group is associated "
        "with the rest for free, so counting phages would let a single "
        "lineage look like a broad spectrum. Distinct groups are counted "
        "below, which is the treatment given to clonal strains and to "
        "identical k-mer presence patterns elsewhere in this project.")
    add("")
    add("DOES ANY SYSTEM ACT ON MORE PHAGE GROUPS THAN A MATCHED TRAIT?")
    add(f"  {'system':24s} {'carr':>5s} {'phyR2':>6s} {'phage':>6s} "
        f"{'group':>6s} {'nCo':>4s} {'nCl':>4s} {'p':>7s} {'best z':>7s} "
        f"{'p best':>7s}")
    for r in spec.head(14).itertuples(index=False):
        add(f"  {r.system:24s} {r.n_carriers:5d} {r.phylogenetic_r2:6.3f} "
            f"{r.n_phages_protective:6d} {r.n_groups_protective:6d} "
            f"{r.null_median_groups_protective_coordinate:4.1f} "
            f"{r.null_median_groups_protective_clonal:4.1f} "
            f"{r.p_protective:7.4f} {r.best_z:7.2f} {r.p_best:7.4f}")
    add("  'phage' counts phages where carriage predicts less lysis beyond "
        "core-genome position and 'group' collapses those to distinct phage "
        "groups; 'nCo' and 'nCl' are the median group counts the coordinate and clonal nulls reach. "
        "'best z' is the strongest single phage and 'p best' judges it "
        "against the best of a matched trait, which is the right way to ask "
        "whether the strongest of many is stronger than chance supplies.")
    add("")

    sig_spec = spec[spec["q_protective"] <= config.FDR_ALPHA]
    best_spec = spec[spec["q_best"] <= config.FDR_ALPHA]
    if sig_spec.empty:
        add(f"  No system acts on more phage groups than a phylogenetically "
            f"matched trait at q <= {config.FDR_ALPHA}. The smallest p is "
            f"{spec['p_protective'].min():.4f} for {spec.iloc[0]['system']}, "
            f"against a floor of {1 / (N_NULL_TRAITS + 1):.4f} set by the "
            "number of null traits.")
        add("")
        add("  That is the answer to the question this script was written to "
            f"ask. The negative results reported elsewhere are not an artifact "
            f"of averaging across phages: resolving the panel into its "
            f"{len(phage_names)} members, and asking of each system whether it "
            "blocks any subset of them, finds no system acting on more of the "
            "panel than its matched null. The averages were not hiding a "
            "narrow spectrum.")
    else:
        word = "system acts" if len(sig_spec) == 1 else "systems act"
        add(f"  {len(sig_spec)} {word} on more phage groups than both matched "
            f"nulls at q <= {config.FDR_ALPHA}: "
            f"{', '.join(sig_spec['system'])}.")
        add("  This is the case the script was written to detect. It means "
            "the panel averages reported elsewhere understate these systems, "
            "and that the phages concerned, not the panel mean, are where "
            "they should be followed up.")
        for r in sig_spec.itertuples(index=False):
            on = per_phage[(per_phage.system == r.system)
                           & (per_phage.score_z < -Z_THRESHOLD)]["phage"]
            sel = pairs[pairs.phage.isin(set(on))]
            car = set(feats.index[feats[f"ds_{r.system}"] == 1])
            inside = sel[sel.strain.isin(car)]["lytic"].mean()
            outside = sel[~sel.strain.isin(car)]["lytic"].mean()
            whole = pairs[pairs.strain.isin(car)]["lytic"].mean()
            rest = pairs[~pairs.strain.isin(car)]["lytic"].mean()
            add(f"    {r.system}: on its {len(on)} phages, carriers are lysed "
                f"{inside:.3f} of the time against {outside:.3f} for "
                f"non-carriers, a difference of {inside - outside:+.3f}. "
                f"Across the whole panel the same comparison is "
                f"{whole:.3f} against {rest:.3f}, a difference of "
                f"{whole - rest:+.3f}. That gap between the two is what "
                "averaging over the panel costs, and it is why every test "
                "elsewhere in this project missed it.")
    if not best_spec.empty:
        add(f"  Separately, {len(best_spec)} systems have a single strongest "
            f"phage beyond what the best of a matched trait reaches: "
            f"{', '.join(best_spec['system'])}. That is a weaker claim than "
            "a spectrum, since it rests on one phage, but it is the form a "
            "genuinely narrow defense would take if it blocked only one "
            "lineage in this panel.")
    if not adj.empty:
        add("")
        add("  DOES THE SURVIVOR KEEP ITS SIGNAL AWAY FROM ITS ISLAND?")
        for r in adj.itertuples(index=False):
            add(f"    {r.system}: {r.n_groups_unadjusted} phage groups "
                f"unadjusted, {r.n_groups_adjusted} after conditioning on "
                f"{r.partners} (p {r.p_adjusted:.4f}).")
        kept = adj[adj["p_adjusted"] <= config.FDR_ALPHA]
        if kept.empty:
            add("    No survivor keeps its signal once the systems most "
                "enriched among its carriers are conditioned on, so the "
                "effect belongs to the island rather than to any one member "
                "of it on this evidence.")
        else:
            add(f"    {', '.join(kept['system'])} keeps its signal with those "
                "neighbours in the model. That does not prove the attribution, "
                "since an island member not catalogued here could carry it, "
                "but the effect is not explained by the neighbours that are.")
    add("")
    add("  Collapsing threshold sensitivity, smallest p over systems:")
    for th in PHAGE_CLUSTER_SWEEP:
        part = sweep_df[sweep_df["threshold"] == th]
        add(f"    distance {th:<5} : {int(part['n_groups_in_panel'].iloc[0]):3d} "
            f"groups in the panel, smallest p "
            f"{part['p_protective'].min():.4f} "
            f"({part.loc[part['p_protective'].idxmin(), 'system']})")
    add("")

    add("IS ANY SPECTRUM MORE CLOSELY RELATED THAN CHANCE?")
    if coh.empty:
        add("  No system has three or more phages to test, so there is no "
            "spectrum whose coherence could be assessed. A system that acted "
            "on a genuine subset of related phages would have produced one.")
    else:
        add(f"  {'system':24s} {'phages':>7s} {'groups':>7s} {'mean d':>8s} "
            f"{'random':>8s} {'p':>7s}")
        for r in coh.head(12).itertuples(index=False):
            add(f"  {r.system:24s} {r.n_phages:7d} {r.n_groups:7d} "
                f"{r.mean_pairwise_distance:8.4f} {r.random_mean:8.4f} "
                f"{r.p_more_related:7.4f}")
        tight = coh[coh["q_more_related"] <= config.FDR_ALPHA]
        if tight.empty:
            add("  No system acts on a set of phages more closely related than "
                "a random set of the same size. A real narrow-spectrum defense "
                "should block relatives, so a scattered set is more consistent "
                "with noise than with a mechanism.")
        else:
            add(f"  {len(tight)} systems act on a set of phages more closely "
                f"related than chance: {', '.join(tight['system'])}.")
            add("  This has to be read against the 'groups' column rather than "
                "'phages'. A set that is tightly clustered because it is one "
                "redundant group of the panel is not evidence of a spectrum: "
                "it is the same observation counted several times, and the "
                "distance that makes it look coherent is the reason it should "
                "not. The test is informative only where the affected phages "
                "span several groups and are still closer together than "
                "chance would place them.")
            multi = tight[tight["n_groups"] >= 3]
            if multi.empty:
                add("  No system here clears that reading: every coherent set "
                    "is one or two groups of near-identical phages.")
            else:
                add(f"  Systems whose coherent set spans three or more groups: "
                    f"{', '.join(multi['system'])}. These are the ones worth "
                    "following up on the phages concerned rather than on "
                    "average.")
    add("")

    add("WHAT THIS DOES AND DOES NOT SETTLE")
    add("  It settles whether averaging over the panel hid a defense system "
        "that works on part of it, which is a specific and plausible way the "
        "negative results elsewhere could have been wrong, and it is answered "
        "with the panel resolved rather than by assumption.")
    add("  It does not rescue a system whose effect is confined to phages "
        "outside this panel, nor one acting on a single phage too weakly for "
        f"{len(strains)} strains to resolve, and it inherits the limitation "
        "that carriage is presence of a catalogued system rather than a "
        "functional allele.")
    add("  The matched null also removes, by construction, any effect that "
        "coincides exactly with core-genome ancestry. That is deliberate and "
        "it is the same blind spot documented for every other test here: an "
        "effect inseparable from the lineage cannot be assigned to the system "
        "by any observational design.")

    text = "\n".join(o)
    (paths.RESULTS / "phage_resolved_report.txt").write_text(text + "\n",
                                                            encoding="utf-8")
    print("\n" + text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
