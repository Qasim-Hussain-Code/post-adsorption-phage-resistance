"""Test which defense associations survive control for bacterial population structure.

Defense islands are inherited vertically and co-segregate with clonal lineages,
so a system that happens to be common in a lineage with an unusual surface
chemistry or an unusual growth habit will show an association with lysis that
has nothing to do with its biochemical activity. Script 03 found fourteen
systems associated with more lysis rather than less, which no defense mechanism
predicts, so this confounding is demonstrably present and not hypothetical.

Clonal background enters in three ways:
  phylogroup       a coarse 15-level Clermont assignment
  principal coords eight axes of classical multidimensional scaling on the
                   core-genome distance matrix, which resolve structure within
                   a phylogroup
  clonal clusters  single-linkage groups of near-identical core genomes, used
                   as the clustering unit for the robust standard errors

Standard errors are clustered on clonal cluster throughout. A strain
contributes up to 96 correlated observations, and strains in a clonal cluster
are near-replicates of one another, so treating pairs as independent would
produce intervals that are far too narrow. A generalized linear mixed model
with a random intercept per strain was considered; a cluster-robust binomial
GLM is reported instead because it makes no assumption about the shape of the
between-lineage distribution and is stable at this number of fixed effects.

Outputs:
    data/processed/core_genome_pcoa.csv
    data/processed/clonal_clusters.csv
    results/tables/table7_lineage_controlled_effects.csv
    results/tables/table8_joint_defense_model.csv
    results/population_structure_report.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

from papr import config, paths
from papr.models import fit_logit_cluster
from papr.stats import benjamini_hochberg

ADJUSTMENTS: dict[str, list[str]] = {
    "base": [],
    "plus_phylogroup": ["phylogroup"],
    "plus_principal_coordinates": ["pcoa"],
    "plus_phylogroup_and_coordinates": ["phylogroup", "pcoa"],
}


def load_distance_matrix(strains: list[str]) -> pd.DataFrame:
    d = pd.read_csv(paths.CORE_DISTANCES, sep="\t", index_col=0)
    d.index = d.index.astype(str)
    d.columns = d.columns.astype(str)
    missing = [s for s in strains if s not in d.index]
    if missing:
        raise SystemExit(f"distance matrix is missing {len(missing)} strains: {missing[:5]}")
    d = d.loc[strains, strains]
    arr = d.to_numpy(dtype=float)
    # Enforce exact symmetry and a zero diagonal before the spectral step.
    arr = 0.5 * (arr + arr.T)
    np.fill_diagonal(arr, 0.0)
    return pd.DataFrame(arr, index=strains, columns=strains)


def principal_coordinates(dist: pd.DataFrame, n_axes: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Classical multidimensional scaling (Gower's principal coordinates).

    This is an eigendecomposition of a doubly centred matrix of squared
    distances. It is a deterministic linear projection with no fitted
    parameters and no objective function to optimise, which is why it is used
    here in place of a nonlinear embedding.
    """
    d2 = dist.to_numpy(dtype=float) ** 2
    n = d2.shape[0]
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ d2 @ j
    b = 0.5 * (b + b.T)
    vals, vecs = np.linalg.eigh(b)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    pos = vals > 0
    total = float(vals[pos].sum())
    k = min(n_axes, int(pos.sum()))
    coords = vecs[:, :k] * np.sqrt(vals[:k])
    axes = pd.DataFrame(coords, index=dist.index,
                        columns=[f"PCo{i + 1}" for i in range(k)])
    axes.index.name = "strain"
    variance = pd.DataFrame({
        "axis": [f"PCo{i + 1}" for i in range(k)],
        "eigenvalue": vals[:k],
        "variance_explained": vals[:k] / total,
        "cumulative_variance_explained": np.cumsum(vals[:k]) / total,
    })
    return axes, variance


def clonal_clusters(dist: pd.DataFrame, threshold: float) -> tuple[pd.Series, pd.DataFrame]:
    """Group near-identical core genomes by single-linkage distance cut."""
    condensed = squareform(dist.to_numpy(dtype=float), checks=False)
    link = hierarchy.linkage(condensed, method="single")
    labels = hierarchy.fcluster(link, t=threshold, criterion="distance")
    series = pd.Series(labels, index=dist.index, name="clonal_cluster")

    sweep = []
    for t in [1e-5, 5e-5, 1e-4, 2.5e-4, 5e-4, 1e-3, 2.5e-3, 5e-3, 1e-2]:
        lab = hierarchy.fcluster(link, t=t, criterion="distance")
        sizes = pd.Series(lab).value_counts()
        sweep.append(dict(threshold=t, n_clusters=int(sizes.size),
                          largest_cluster=int(sizes.iloc[0]),
                          n_singletons=int((sizes == 1).sum())))
    return series, pd.DataFrame(sweep)


def design_blocks(pairs: pd.DataFrame, pcoa: pd.DataFrame) -> dict[str, np.ndarray]:
    """Fixed-effect blocks, built once and reused across every model fit."""
    const = np.ones((len(pairs), 1))
    phage = pd.get_dummies(pairs["phage"].astype("string"),
                           drop_first=True, dtype=float).to_numpy()
    lps = pd.get_dummies(pairs["lps_core"].astype("string"),
                         drop_first=True, dtype=float).to_numpy()
    phylo = pd.get_dummies(pairs["phylogroup"].astype("string"),
                           drop_first=True, dtype=float).to_numpy()
    coords = pcoa.reindex(pairs["strain"]).to_numpy(dtype=float)
    # Standardise the coordinates so that the eight axes enter on a common
    # scale and the optimiser is well conditioned.
    coords = (coords - coords.mean(axis=0)) / coords.std(axis=0)
    return {"const": const, "phage": phage, "lps": lps,
            "phylogroup": phylo, "pcoa": coords}


def fit_one(y: np.ndarray, exposure: np.ndarray, blocks: dict[str, np.ndarray],
            extras: list[str], groups: np.ndarray) -> dict:
    parts = [exposure.reshape(-1, 1), blocks["const"], blocks["phage"], blocks["lps"]]
    parts += [blocks[k] for k in extras]
    x = np.hstack(parts)
    res = fit_logit_cluster(y, x, groups)
    beta = float(res.params[0])
    se = float(res.bse[0])
    p = float(res.pvalues[0])
    return dict(odds_ratio=float(np.exp(beta)),
                ci_low=float(np.exp(beta - 1.96 * se)),
                ci_high=float(np.exp(beta + 1.96 * se)),
                log2_or=float(beta / np.log(2.0)),
                std_error=se, p_value=p, n_parameters=int(x.shape[1]),
                converged=bool(getattr(res, "mle_retvals", {}).get("converged", True)))


def joint_model(y: np.ndarray, pairs: pd.DataFrame, systems: list[str],
                blocks: dict[str, np.ndarray], groups: np.ndarray) -> pd.DataFrame:
    """All defense systems together, mutually adjusted and lineage adjusted."""
    exposures = pairs[systems].to_numpy(dtype=float)
    x = np.hstack([exposures, blocks["const"], blocks["phage"], blocks["lps"],
                   blocks["phylogroup"], blocks["pcoa"]])
    res = fit_logit_cluster(y, x, groups)
    k = len(systems)
    beta = np.asarray(res.params[:k], dtype=float)
    se = np.asarray(res.bse[:k], dtype=float)
    p = np.asarray(res.pvalues[:k], dtype=float)
    out = pd.DataFrame({
        "system": [s[3:] for s in systems],
        "odds_ratio": np.exp(beta),
        "ci_low": np.exp(beta - 1.96 * se),
        "ci_high": np.exp(beta + 1.96 * se),
        "log2_or": beta / np.log(2.0),
        "std_error": se, "p_value": p,
    })
    out["q_value"] = benjamini_hochberg(out["p_value"].to_numpy())
    out["significant"] = out["q_value"] <= config.FDR_ALPHA
    return out.sort_values("odds_ratio")


def main() -> int:
    paths.bootstrap()
    pairs = pd.read_csv(paths.PAIRS, compression="gzip", low_memory=False)
    strain_feats = pd.read_csv(paths.STRAIN_FEATURES, low_memory=False)
    ladder = pd.read_csv(paths.TABLES / "table3_defense_conditioning_ladder.csv")

    ds_cols = [c for c in strain_feats.columns if c.startswith("ds_")]
    tested = strain_feats[strain_feats["strain"].isin(set(pairs["strain"]))]
    systems = sorted([c for c in ds_cols
                      if tested[c].sum() >= config.MIN_CARRIER_STRAINS])
    pairs = pairs.merge(strain_feats[["strain"] + systems], on="strain", how="left")

    strains = sorted(pairs["strain"].unique())
    dist = load_distance_matrix(strains)
    pcoa, variance = principal_coordinates(dist, config.N_PCOA_AXES)
    clusters, sweep = clonal_clusters(dist, config.CLONAL_DISTANCE_THRESHOLD)

    pcoa.to_csv(paths.PCOA)
    clusters.to_frame().to_csv(paths.PROCESSED / "clonal_clusters.csv")
    variance.to_csv(paths.TABLES / "table7b_pcoa_variance.csv", index=False)

    pairs["clonal_cluster"] = pairs["strain"].map(clusters)
    groups = pairs["clonal_cluster"].to_numpy()
    y = pairs["lytic"].to_numpy(dtype=float)
    blocks = design_blocks(pairs, pcoa)

    print(f"fitting {len(systems)} systems x {len(ADJUSTMENTS)} adjustments")
    rows = []
    for i, sys_col in enumerate(systems, 1):
        exposure = pairs[sys_col].to_numpy(dtype=float)
        for label, extras in ADJUSTMENTS.items():
            res = fit_one(y, exposure, blocks, extras, groups)
            rows.append(dict(system=sys_col[3:], adjustment=label, **res))
        if i % 8 == 0 or i == len(systems):
            print(f"  {i}/{len(systems)}")
    eff = pd.DataFrame(rows)
    eff["q_value"] = np.nan
    for label in eff["adjustment"].unique():
        m = eff["adjustment"] == label
        eff.loc[m, "q_value"] = benjamini_hochberg(eff.loc[m, "p_value"].to_numpy())
    eff["significant"] = eff["q_value"] <= config.FDR_ALPHA

    # Attenuation relative to the unadjusted stratified estimate from script 03.
    l3 = (ladder[ladder["level"] == "L3_phage_lps"]
          .set_index("system")[["estimate", "log2_or"]]
          .rename(columns={"estimate": "cmh_odds_ratio", "log2_or": "cmh_log2_or"}))
    eff = eff.join(l3, on="system")
    # The retained fraction is only meaningful when there was an effect to
    # retain. Dividing by a stratified log2 odds ratio that is already within
    # noise of zero produces arbitrarily large ratios, so those are left
    # undefined rather than reported as though they meant something.
    MIN_INTERPRETABLE_LOG2_OR = 0.1
    interpretable = eff["cmh_log2_or"].abs() >= MIN_INTERPRETABLE_LOG2_OR
    eff["retained_fraction_of_log2_effect"] = np.where(
        interpretable, eff["log2_or"] / eff["cmh_log2_or"], np.nan)
    eff["direction_preserved"] = (
        np.sign(eff["log2_or"]) == np.sign(eff["cmh_log2_or"]))
    eff.to_csv(paths.TABLES / "table7_lineage_controlled_effects.csv", index=False)

    joint = joint_model(y, pairs, systems, blocks, groups)
    joint.to_csv(paths.TABLES / "table8_joint_defense_model.csv", index=False)

    full = eff[eff["adjustment"] == "plus_phylogroup_and_coordinates"].set_index("system")
    base = eff[eff["adjustment"] == "base"].set_index("system")

    out: list[str] = []
    add = out.append
    add("Population structure control report")
    add("Generated by scripts/04_control_population_structure.py")
    add("")
    add("CORE-GENOME PRINCIPAL COORDINATES")
    add(f"  strains embedded: {len(strains)}")
    for r in variance.itertuples(index=False):
        add(f"  {r.axis:6s} eigenvalue {r.eigenvalue:10.5f}  "
            f"variance {100 * r.variance_explained:5.2f} percent  "
            f"cumulative {100 * r.cumulative_variance_explained:5.2f} percent")
    add(f"  The {config.N_PCOA_AXES} retained axes carry "
        f"{100 * variance['cumulative_variance_explained'].iloc[-1]:.1f} percent "
        "of the variation in core-genome distance.")

    add("")
    add("CLONAL CLUSTERS USED FOR THE ROBUST STANDARD ERRORS")
    add(f"  {'threshold':>10s} {'clusters':>9s} {'largest':>8s} {'singletons':>11s}")
    for r in sweep.itertuples(index=False):
        mark = "  <- used" if abs(r.threshold - config.CLONAL_DISTANCE_THRESHOLD) < 1e-12 else ""
        add(f"  {r.threshold:10.5f} {r.n_clusters:9d} {r.largest_cluster:8d} "
            f"{r.n_singletons:11d}{mark}")
    n_clusters = int(clusters.nunique())
    add(f"  Clusters at the chosen threshold: {n_clusters}, against "
        f"{tested['sequence_type'].nunique()} Sequence Types and "
        f"{tested['phylogroup'].nunique()} phylogroups among the same strains.")

    add("")
    add("EFFECT OF LINEAGE ADJUSTMENT ON THE DEFENSE ASSOCIATIONS")
    for label in ADJUSTMENTS:
        sub = eff[eff["adjustment"] == label]
        add(f"  {label:34s} {int(sub['significant'].sum()):2d} of {len(sub)} "
            f"significant at q <= {config.FDR_ALPHA}, "
            f"median |log2 OR| = {sub['log2_or'].abs().median():.3f}, "
            f"mean parameters = {sub['n_parameters'].mean():.0f}")
    add("")
    add("  The first row is the same contrast as the receptor-conditioned "
        "estimate in script 03 but fitted as a regression with standard errors "
        "clustered on clonal background. Comparing it with the Mantel-Haenszel "
        "intervals shows how much of the apparent precision in script 03 came "
        "from treating correlated pairs as independent.")

    add("")
    add("PER-SYSTEM COMPARISON: STRATIFIED, THEN FULLY LINEAGE ADJUSTED")
    add(f"  {'system':26s} {'CMH OR':>7s} {'base OR':>8s} {'adj OR':>7s} "
        f"{'95 pct CI':>16s} {'adj q':>9s} {'retained':>9s}")
    for s in full.sort_values("odds_ratio").index:
        r = full.loc[s]
        b = base.loc[s]
        mark = " *" if r["significant"] else "  "
        ret = r["retained_fraction_of_log2_effect"]
        ret_s = "      n/a" if not np.isfinite(ret) else f"{ret:9.2f}"
        add(f"  {s:26s} {r['cmh_odds_ratio']:7.3f} {b['odds_ratio']:8.3f} "
            f"{r['odds_ratio']:7.3f} [{r['ci_low']:6.3f}, {r['ci_high']:6.3f}] "
            f"{r['q_value']:9.2e} {ret_s}{mark}")
    add("  'retained' is the ratio of the lineage-adjusted log2 odds ratio to "
        "the stratified one. A value near 1 means lineage adjustment left the "
        "effect intact; near 0 means the association was almost entirely a "
        "property of the lineage; a negative value means the effect reversed. "
        "It is shown as n/a where the stratified effect was already within "
        "noise of zero, since a ratio to approximately zero carries no "
        "information.")

    surv = full[full["significant"] & full["direction_preserved"]]
    prot = surv[surv["odds_ratio"] < 1]
    inv = surv[surv["odds_ratio"] > 1]
    lost = base[base["significant"]].index.difference(surv.index)
    add("")
    add("WHAT SURVIVES")
    add(f"  systems significant before lineage adjustment : "
        f"{int(base['significant'].sum())} of {len(base)}")
    add(f"  systems significant after full adjustment     : "
        f"{int(full['significant'].sum())} of {len(full)}")
    add(f"  surviving with direction preserved            : {len(surv)}")
    add(f"    associated with reduced lysis               : {len(prot)} "
        f"({', '.join(sorted(prot.index)) if len(prot) else 'none'})")
    add(f"    associated with increased lysis             : {len(inv)} "
        f"({', '.join(sorted(inv.index)) if len(inv) else 'none'})")
    add(f"  lost on lineage adjustment                    : {len(lost)} "
        f"({', '.join(sorted(lost)) if len(lost) else 'none'})")
    add("")
    med_ret = full.loc[full["direction_preserved"], "retained_fraction_of_log2_effect"].median()
    add(f"  Median retained fraction of the log2 effect among systems that keep "
        f"their direction: {med_ret:.2f}.")
    add("")
    n_inverse_before = int(((base["odds_ratio"] > 1) & base["significant"]).sum())
    if len(inv) == 0:
        add("  The associations pointing the wrong way provide an internal check "
            "on the adjustment, and it passes. A defense system cannot make its "
            "host easier to lyse, so every significant association with "
            "increased lysis must be confounding. Before lineage adjustment "
            f"there were {n_inverse_before} such associations in this "
            "regression and 14 under the stratified estimates of script 03. "
            "After adjusting for phylogroup and core-genome position there are "
            "none. The covariates are therefore removing the kind of signal "
            "they were introduced to remove, which makes the attenuation of "
            "the protective associations harder to dismiss as overcorrection.")
    else:
        add(f"  {len(inv)} systems remain associated with increased lysis after "
            "adjusting for phylogroup and core-genome position. A defense "
            "system cannot make a host easier to lyse, so these are either "
            "residual confounding that these covariates do not capture, or "
            "markers for something else in the genome, for example a mobile "
            "element carrying both the defense locus and a determinant of "
            "phage susceptibility. They are reported as an unresolved feature "
            "of the data rather than interpreted as biology.")

    add("")
    add("JOINT MODEL: ALL SYSTEMS MUTUALLY ADJUSTED")
    add(f"  {'system':26s} {'OR':>7s} {'95 pct CI':>16s} {'q':>9s}")
    for r in joint.itertuples(index=False):
        mark = " *" if r.significant else "  "
        add(f"  {r.system:26s} {r.odds_ratio:7.3f} [{r.ci_low:6.3f}, {r.ci_high:6.3f}] "
            f"{r.q_value:9.2e}{mark}")
    add(f"  Significant in the joint model: {int(joint['significant'].sum())} "
        f"of {len(joint)}. Defense systems co-occur, so a coefficient here is "
        "the association of one system holding the other 32 fixed, which is a "
        "stricter and more collinear contrast than the per-system fits above.")

    text = "\n".join(out)
    (paths.RESULTS / "population_structure_report.txt").write_text(
        text + "\n", encoding="utf-8")
    print("\n" + text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
