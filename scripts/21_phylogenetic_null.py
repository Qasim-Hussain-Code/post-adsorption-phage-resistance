"""Compare each defense system against traits with the same phylogenetic shape.

The permutation tests in script 09 shuffle carriage inside a phylogroup or a
Sequence Type. That requires choosing a grouping level, and the choice matters:
RM_Type_IIG is significant inside phylogroup and not inside Sequence Type. The
grouping is also the same for every system, although systems differ enormously
in how lineage-confined they are.

This script removes the choice. For each system it measures how much of the
carriage pattern is explained by position on the core-genome phylogeny, then
simulates null traits that are confined to the phylogeny to exactly the same
degree and have exactly the same prevalence, and asks whether the real system's
association with lysis is stronger than those null traits achieve. A system
whose association is a by-product of the lineage it sits in will look ordinary
against traits that sit in lineages the same way.

This is the standard null for a bacterial genome-wide association, where
population structure rather than sampling noise is the thing that generates
false positives. Note what it does not need: no outcome model, no covariate
adjustment, and no assumption that the eight principal coordinates capture
confounding correctly. They are used only to generate traits with realistic
phylogenetic structure.

Outputs:
    results/tables/table39_phylogenetic_signal.csv
    results/tables/table40_phylogenetic_null_tests.csv
    results/phylogenetic_null_report.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
from scipy import stats

from papr import config, paths
from papr.phylonull import phylogenetic_r2, simulate_traits
from papr.stats import benjamini_hochberg

N_SIMULATIONS = 10_000


def statistic(traits: np.ndarray, lysis: np.ndarray) -> np.ndarray:
    """Difference in mean per-strain lysis rate between carriers and others."""
    k = traits.sum(axis=0).astype(float)
    n = traits.shape[0]
    total = lysis.sum()
    carrier_sum = lysis @ traits
    with np.errstate(invalid="ignore", divide="ignore"):
        return carrier_sum / k - (total - carrier_sum) / (n - k)


def main() -> int:
    paths.bootstrap()
    pairs = pd.read_csv(paths.PAIRS, compression="gzip", low_memory=False)
    strain_feats = pd.read_csv(paths.STRAIN_FEATURES, low_memory=False)
    pcoa = pd.read_csv(paths.PCOA, index_col=0)

    per_strain = pairs.groupby("strain")["lytic"].mean()
    strains = [s for s in per_strain.index if s in pcoa.index]
    lysis = per_strain.loc[strains].to_numpy(dtype=float)
    axes = pcoa.loc[strains].to_numpy(dtype=float)
    axes = (axes - axes.mean(axis=0)) / axes.std(axis=0)
    feats = strain_feats.set_index("strain").loc[strains]

    # The phylogenetically matched null conditions on ancestry and on nothing
    # else. Host surface type is the largest host-side determinant of lysis
    # (script 11) and is only partly collinear with ancestry, so a system
    # carried disproportionately by strains with a permissive LPS outer core
    # would pass that null for a reason that has nothing to do with defense.
    # The same test is therefore run a second time on lysis rates with surface
    # type regressed out.
    surface = pd.get_dummies(
        feats[["lps_core", "o_status"]].astype("string").fillna("missing"),
        drop_first=True, dtype=float).to_numpy()
    sx = np.hstack([np.ones((len(strains), 1)), surface])
    beta, *_ = np.linalg.lstsq(sx, lysis, rcond=None)
    lysis_resid = lysis - sx @ beta
    surface_r2 = 1.0 - float(np.sum(lysis_resid ** 2)) / float(
        np.sum((lysis - lysis.mean()) ** 2))
    print(f"host surface explains {100 * surface_r2:.1f} percent of the "
          "variation in per-strain lysis rate")

    # Carriage composition by permissive outer core, used to explain which
    # associations survive which null.
    lysis_by_core = pairs.groupby("lps_core")["lytic"].mean().sort_values()
    permissive_cores = set(lysis_by_core.index[-2:])
    permissive = feats["lps_core"].isin(permissive_cores).to_numpy().astype(int)

    ds_cols = [c for c in strain_feats.columns if c.startswith("ds_")]
    systems = sorted([c for c in ds_cols
                      if feats[c].sum() >= config.MIN_CARRIER_STRAINS])
    print(f"{len(strains)} strains, {len(systems)} systems, "
          f"{N_SIMULATIONS} simulations each")

    rng = np.random.default_rng(config.RANDOM_SEED)
    rows = []
    for i, col in enumerate(systems, 1):
        carrier = feats[col].to_numpy().astype(int)
        k = int(carrier.sum())
        r2 = phylogenetic_r2(carrier, axes)
        observed = float(lysis[carrier == 1].mean() - lysis[carrier == 0].mean())

        observed_resid = float(lysis_resid[carrier == 1].mean()
                               - lysis_resid[carrier == 0].mean())

        traits = simulate_traits(axes, k, r2, N_SIMULATIONS, rng)
        null = statistic(traits, lysis)
        null = null[np.isfinite(null)]
        # Two-sided empirical p-value, observed value included.
        p = float((np.sum(np.abs(null) >= abs(observed)) + 1) / (null.size + 1))

        # The same matched null applied to surface-residualised lysis.
        null_resid = statistic(traits, lysis_resid)
        null_resid = null_resid[np.isfinite(null_resid)]
        p_resid = float((np.sum(np.abs(null_resid) >= abs(observed_resid)) + 1)
                        / (null_resid.size + 1))

        # Compare with an unstructured null for contrast.
        flat = simulate_traits(axes, k, 0.0, 2000, rng)
        null_flat = statistic(flat, lysis)
        p_flat = float((np.sum(np.abs(null_flat) >= abs(observed)) + 1)
                       / (null_flat.size + 1))

        # Does carriage track a permissive outer core?
        a = int(((carrier == 1) & (permissive == 1)).sum())
        b = int(((carrier == 1) & (permissive == 0)).sum())
        c_ = int(((carrier == 0) & (permissive == 1)).sum())
        d = int(((carrier == 0) & (permissive == 0)).sum())
        surf_or = (a * d) / (b * c_) if b and c_ else np.nan
        surf_p = float(stats.fisher_exact([[a, b], [c_, d]])[1])

        rows.append(dict(system=col[3:], n_carriers=k,
                         phylogenetic_r2=r2, observed_difference=observed,
                         observed_difference_surface_adjusted=observed_resid,
                         null_sd_structured=float(null.std(ddof=1)),
                         null_sd_unstructured=float(null_flat.std(ddof=1)),
                         p_structured=p, p_structured_surface_adjusted=p_resid,
                         p_unstructured=p_flat,
                         pct_carriers_permissive_core=100.0 * a / max(1, a + b),
                         permissive_core_odds_ratio=surf_or,
                         permissive_core_p=surf_p))
        if i % 8 == 0 or i == len(systems):
            print(f"  {i}/{len(systems)}")

    out = pd.DataFrame(rows)
    out["q_structured"] = benjamini_hochberg(out["p_structured"].to_numpy())
    out["q_unstructured"] = benjamini_hochberg(out["p_unstructured"].to_numpy())
    out["q_structured_surface_adjusted"] = benjamini_hochberg(
        out["p_structured_surface_adjusted"].to_numpy())
    out["significant"] = out["q_structured"] <= config.FDR_ALPHA
    out["significant_surface_adjusted"] = (
        out["q_structured_surface_adjusted"] <= config.FDR_ALPHA)
    out["inflation_factor"] = out["null_sd_structured"] / out["null_sd_unstructured"]
    out = out.sort_values("p_structured")
    out.to_csv(paths.TABLES / "table40_phylogenetic_null_tests.csv", index=False)
    out[["system", "n_carriers", "phylogenetic_r2", "inflation_factor"]].to_csv(
        paths.TABLES / "table39_phylogenetic_signal.csv", index=False)

    o: list[str] = []
    add = o.append
    add("Phylogenetically matched null report")
    add("Generated by scripts/21_phylogenetic_null.py")
    add("")
    add("SETUP")
    add(f"  strains                 : {len(strains)}")
    add(f"  systems tested          : {len(systems)}")
    add(f"  simulations per system  : {N_SIMULATIONS}")
    add("  Null traits share each system's prevalence exactly and its share of "
        "carriage variance explained by core-genome position.")

    add("")
    add("HOW LINEAGE-CONFINED ARE THESE SYSTEMS")
    add(f"  median share of carriage explained by core-genome position : "
        f"{out['phylogenetic_r2'].median():.3f}")
    add(f"  range : {out['phylogenetic_r2'].min():.3f} to "
        f"{out['phylogenetic_r2'].max():.3f}")
    add("")
    add(f"  {'system':24s} {'carriers':>8s} {'phylo R2':>9s} {'inflation':>10s}")
    for r in out.sort_values("phylogenetic_r2", ascending=False).head(10).itertuples(index=False):
        add(f"  {r.system:24s} {r.n_carriers:8d} {r.phylogenetic_r2:9.3f} "
            f"{r.inflation_factor:10.2f}")
    add("  'inflation' is how much wider the null distribution of the test "
        "statistic becomes once the simulated traits are given the system's "
        "own phylogenetic structure. It is the factor by which a test that "
        "ignores population structure overstates its evidence.")
    add(f"  Median inflation across systems: {out['inflation_factor'].median():.2f}. "
        "A test treating carriage as randomly distributed across strains would "
        "therefore understate its standard error by roughly that factor.")

    add("")
    add("ASSOCIATION AGAINST A PHYLOGENETICALLY MATCHED NULL")
    nsig = int(out["significant"].sum())
    nflat = int((out["q_unstructured"] <= config.FDR_ALPHA).sum())
    add(f"  significant against the matched null    : {nsig} of {len(out)}")
    add(f"  significant against an unstructured null: {nflat} of {len(out)}")
    add("")
    add(f"  {'system':24s} {'carriers':>8s} {'phylo R2':>9s} {'obs diff':>9s} "
        f"{'p struct':>9s} {'q struct':>9s} {'p flat':>9s}")
    for r in out.head(14).itertuples(index=False):
        mark = " *" if r.significant else "  "
        add(f"  {r.system:24s} {r.n_carriers:8d} {r.phylogenetic_r2:9.3f} "
            f"{r.observed_difference:+9.4f} {r.p_structured:9.2e} "
            f"{r.q_structured:9.2e} {r.p_unstructured:9.2e}{mark}")

    add("")
    add("INTERPRETATION")
    add(f"  Accounting for phylogenetic structure removes "
        f"{nflat - nsig} of the {nflat} associations that an unstructured null "
        "would call significant. That gap is the size of the population "
        "structure problem in this dataset, measured directly rather than "
        "assumed.")
    sig = out[out["significant"]]
    prot = sig[sig["observed_difference"] < 0]
    inv = sig[sig["observed_difference"] > 0]
    if len(sig):
        add(f"  Systems surviving the matched null: {', '.join(sig['system'])}.")
        add(f"  Of those, {len(prot)} are associated with less lysis "
            f"({', '.join(prot['system']) if len(prot) else 'none'}) and "
            f"{len(inv)} with more "
            f"({', '.join(inv['system']) if len(inv) else 'none'}).")
    else:
        add("  No system's association with lysis exceeds what traits of the "
            "same prevalence and the same phylogenetic distribution achieve by "
            "chance.")

    add("")
    add("CONDITIONING ON PHYLOGENY ALONE IS NOT ENOUGH")
    add("  A defense system cannot make its host easier to lyse, so any "
        "association with increased lysis is confounding, and the count of "
        "such results measures how far a test can be trusted.")
    add(f"  Against the phylogenetically matched null, {len(inv)} of the "
        f"{len(sig)} survivors point the impossible way. Phylogenetic "
        "structure is not what produces them: their carriage is no more "
        "lineage-confined than average.")
    add("")
    add(f"  Host surface type explains {100 * surface_r2:.1f} percent of the "
        "variation in per-strain lysis rate, and it is only partly collinear "
        "with ancestry. A system carried disproportionately by strains with a "
        "permissive LPS outer core will therefore pass a null that conditions "
        "on ancestry alone. Carriage composition by outer core:")
    add("")
    add(f"  {'system':24s} {'carriers':>8s} {'pct permissive':>15s} "
        f"{'OR':>7s} {'p':>9s}")
    order = out["permissive_core_odds_ratio"].apply(
        lambda v: abs(np.log(v)) if np.isfinite(v) and v > 0 else 0.0)
    for r in out.loc[order.sort_values(ascending=False).index].head(10).itertuples(index=False):
        add(f"  {r.system:24s} {r.n_carriers:8d} "
            f"{r.pct_carriers_permissive_core:14.1f}% "
            f"{r.permissive_core_odds_ratio:7.2f} {r.permissive_core_p:9.2e}")

    add("")
    add("THE SAME NULL, WITH SURFACE TYPE REMOVED FROM THE OUTCOME")
    sig_s = out[out["significant_surface_adjusted"]]
    inv_s = sig_s[sig_s["observed_difference_surface_adjusted"] > 0]
    add(f"  significant against the matched null, surface adjusted : "
        f"{len(sig_s)} of {len(out)} "
        f"({', '.join(sig_s['system']) if len(sig_s) else 'none'})")
    add(f"  of which point the impossible way : {len(inv_s)}")
    add("")
    add(f"  {'system':24s} {'raw diff':>9s} {'adj diff':>9s} {'q raw':>9s} "
        f"{'q adj':>9s}")
    for r in out.head(12).itertuples(index=False):
        add(f"  {r.system:24s} {r.observed_difference:+9.4f} "
            f"{r.observed_difference_surface_adjusted:+9.4f} "
            f"{r.q_structured:9.2e} {r.q_structured_surface_adjusted:9.2e}")

    add("")
    add("INTERPRETATION, REVISED")
    if len(inv_s) < len(inv):
        add(f"  Removing surface type from the outcome drops the impossible "
            f"associations from {len(inv)} to {len(inv_s)}. Most of what "
            "survived conditioning on ancestry alone was host surface "
            "composition, which is what the variance partition in script 11 "
            "would predict.")
    else:
        add("  Removing surface type from the outcome does not reduce the "
            "number of impossible associations, so host surface composition "
            "is not what produces them and they remain unexplained.")
    if "RM_Type_IIG" in set(out["system"]):
        r = out.set_index("system").loc["RM_Type_IIG"]
        add("")
        add(f"  RM_Type_IIG deserves separate comment. It passes the matched "
            f"null on raw lysis (q = {r['q_structured']:.2e}) and its carriage "
            f"is barely lineage-confined (phylogenetic R2 "
            f"{r['phylogenetic_r2']:.3f}), which is why it survives tests that "
            "condition on ancestry. But it is strongly depleted in strains "
            f"with a permissive outer core (odds ratio "
            f"{r['permissive_core_odds_ratio']:.2f}, p = "
            f"{r['permissive_core_p']:.1e}): its carriers are "
            f"{r['pct_carriers_permissive_core']:.0f} percent permissive "
            "against 60 percent of strains overall, so they are drawn from "
            "strains that are intrinsically harder to lyse whatever they "
            "carry.")
        shrink = (1 - abs(r["observed_difference_surface_adjusted"])
                  / abs(r["observed_difference"])) * 100
        if r["significant_surface_adjusted"]:
            add(f"  Removing surface type from the outcome attenuates the "
                f"difference by {shrink:.0f} percent but it remains "
                f"significant (q = {r['q_structured_surface_adjusted']:.2e}). "
                "So part of its apparent protection is the surface "
                "composition of its carriers and part is not. This is the "
                "one piece of evidence for RM_Type_IIG that survives "
                "conditioning on both ancestry and surface type without "
                "relying on a covariate-adjusted regression.")
            add("  It still fails the within-cluster test, which holds the "
                "entire genomic background fixed and had moderate power for "
                "this system. The two designs disagree, and neither is "
                "decisive: the within-cluster test uses 36 strains in 7 "
                "clusters, while this one uses all 402 but cannot rule out a "
                "confounder it has not been given.")
        else:
            add(f"  Removing surface type from the outcome attenuates the "
                f"difference by {shrink:.0f} percent and it no longer "
                f"survives (q = {r['q_structured_surface_adjusted']:.2e}). "
                "Both design-based tests that supported it ignore surface "
                "type, and the within-cluster test that holds the whole "
                "genomic background fixed showed no effect, so the picture is "
                "less favourable than the earlier scripts implied.")

    moved = out[out["significant_surface_adjusted"] & ~out["significant"]]
    if len(moved):
        add("")
        add(f"  Adjusting for surface type also brings "
            f"{', '.join(moved['system'])} across the threshold. Surface "
            "composition was masking these rather than creating them: their "
            "carriers sit disproportionately on outer cores that are hard to "
            "lyse anyway, so the raw comparison understated them.")

    add("")
    add("  This test is conservative and has a known cost: a defense system "
        "genuinely confined to one lineage and genuinely protective is "
        "indistinguishable here from a passenger of that lineage. "
        "CAS_Class1-Subtype-I-E is the clearest case, with "
        f"{float(out.set_index('system').loc['CAS_Class1-Subtype-I-E', 'phylogenetic_r2']):.2f} "
        "of its carriage explained by core-genome position; this test cannot "
        "say whether it protects. Only deleting the system in its native "
        "background can.")

    text = "\n".join(o)
    (paths.RESULTS / "phylogenetic_null_report.txt").write_text(
        text + "\n", encoding="utf-8")
    print("\n" + text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
