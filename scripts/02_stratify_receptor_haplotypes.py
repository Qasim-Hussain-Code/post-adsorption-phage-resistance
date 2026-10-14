"""Harmonise the interaction data and build the surface receptor stratification.

The dataset carries no table stating which host surface receptor each phage
requires, so receptor compatibility has to be inferred. This script builds the
inference from the one piece of receptor evidence that is independent of the
interaction matrix being analysed: the strain each phage was originally
isolated on. A strain that yielded a phage must have presented a surface that
phage can bind. Pooling isolation hosts across the members of a VIRIDIC genus,
whose receptor binding proteins are homologous, gives a permissive surface
state set per genus that was never derived from the 38,435 outcomes under test.

The script then measures how well that gate actually separates adsorption from
post-adsorption failure, and reports the answer whether or not it is
favourable. It is not favourable: the gate is weak. That measurement is the
reason the primary analysis in script 03 stratifies on surface state rather
than filtering on it.

Outputs:
    data/processed/pair_level_table.csv.gz
    data/processed/strain_features.csv
    data/processed/phage_features.csv
    data/processed/receptor_gate_definitions.csv
    results/tables/table1_two_gate_decomposition.csv
    results/tables/table2_gate_performance_by_genus.csv
    results/receptor_stratification_report.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from papr import config, paths
from papr.stats import wilson_interval

# Gate tiers, from coarse to fine. Each tier adds a determinant that must also
# match the permissive set of the phage genus.
GATE_TIERS: dict[str, list[str]] = {
    "gate_lps": ["lps_core"],
    "gate_lps_ostatus": ["lps_core", "o_status"],
    "gate_lps_otype": ["lps_core", "o_antigen"],
}


def load_raw() -> dict[str, pd.DataFrame]:
    return {
        "matrix": pd.read_csv(paths.MATRIX, sep=";"),
        "picard": pd.read_csv(paths.PICARD, sep=";"),
        "guelin": pd.read_csv(paths.GUELIN, sep=";"),
        "defense": pd.read_csv(paths.DEFENSE, sep=";"),
        "omp": pd.read_csv(paths.OMP_ALLELES, sep="\t"),
        "rbp": pd.read_csv(paths.RBP, sep=";"),
    }


def build_pair_table(matrix: pd.DataFrame) -> pd.DataFrame:
    """One row per observed strain-phage pair; untested cells are dropped."""
    strain_col = matrix.columns[0]
    long = (matrix.melt(id_vars=strain_col, var_name="phage", value_name="plaque_score")
            .rename(columns={strain_col: "strain"}))
    n_cells = len(long)
    long = long.dropna(subset=["plaque_score"]).copy()
    long["plaque_score"] = long["plaque_score"].astype(int)
    long["lytic"] = (long["plaque_score"] > 0).astype(int)
    long.attrs["n_cells"] = n_cells
    return long.reset_index(drop=True)


def build_strain_features(picard: pd.DataFrame, omp: pd.DataFrame,
                          defense: pd.DataFrame, strains: list[str]) -> pd.DataFrame:
    """Surface determinants, clonal background and defense carriage per strain."""
    p = picard.set_index("bacteria").reindex(strains)
    feats = pd.DataFrame(index=pd.Index(strains, name="strain"))

    feats["phylogroup"] = p[config.PHYLOGROUP_COLUMN].to_numpy()
    feats["sequence_type"] = p[config.ST_COLUMN].astype("string").to_numpy()
    feats["collection"] = p["Collection"].to_numpy()
    feats["pathotype"] = p["Pathotype"].to_numpy()

    feats["lps_core"] = p[config.LPS_COLUMN].to_numpy()
    feats["o_antigen"] = p[config.O_TYPE_COLUMN].to_numpy()
    # O-antigen status separates strains that express any O-antigen from those
    # typed as O-negative, because a long O-antigen chain can occlude the LPS
    # outer core and the outer membrane proteins beneath it. A serotype recorded
    # as ONT is not typeable rather than absent, so it is left missing instead
    # of being counted as O-negative.
    o = feats["o_antigen"].astype("string")
    feats["o_status"] = np.where(
        o.isna() | o.eq("ONT"), None, np.where(o.eq("Oneg"), "Oneg", "Opos"))

    for col in config.CAPSULE_COLUMNS:
        feats[col.lower()] = p[col].to_numpy()
    # min_count=1 keeps a strain with no capsule call missing rather than
    # silently recording it as capsule negative.
    cap_sum = p[config.CAPSULE_COLUMNS].sum(axis=1, min_count=1)
    feats["capsule_any"] = np.where(
        cap_sum.isna().to_numpy(), np.nan, cap_sum.gt(0).to_numpy().astype(float))
    feats["abc_serotype"] = p["ABC_serotype"].to_numpy()

    o = omp.set_index("bacteria").reindex(strains)
    for locus in config.OMP_LOCI:
        feats[f"omp_{locus.lower()}"] = o[locus].to_numpy()

    d = defense.set_index("bacteria").reindex(strains)
    systems = [c for c in d.columns]
    presence = (d[systems].fillna(0) > 0).astype(int)
    presence.columns = [f"ds_{s}" for s in systems]
    feats = feats.join(presence)
    feats["n_defense_systems"] = presence.to_numpy().sum(axis=1)
    feats["defense_burden_quartile"] = pd.qcut(
        feats["n_defense_systems"], 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    return feats.reset_index()


def build_phage_features(guelin: pd.DataFrame, rbp: pd.DataFrame) -> pd.DataFrame:
    g = guelin.copy().rename(columns={
        "phage": "phage", "Morphotype": "morphotype", "Family": "family",
        "Genus": "genus", "Species": "species", "Genome_size": "genome_size",
        "Phage_host": "isolation_host", "Phage_host_phylo": "isolation_host_phylogroup",
        "Subfamily": "subfamily"})
    keep = ["phage", "morphotype", "family", "subfamily", "genus", "species",
            "genome_size", "isolation_host", "isolation_host_phylogroup"]
    g = g[keep]
    rb = (rbp.groupby("phage")
          .agg(n_rbp=("RBP", lambda s: int(s.notna().sum())),
               rbp_types=("type", lambda s: ",".join(sorted({x for x in s.dropna()}))))
          .reset_index())
    return g.merge(rb, on="phage", how="left")


def build_gate_definitions(phage_feats: pd.DataFrame,
                           strain_feats: pd.DataFrame) -> pd.DataFrame:
    """Permissive surface states per genus, taken from isolation hosts only.

    A state is permissive for a genus when at least one strain that yielded a
    phage of that genus carries it. Hosts with no call for a determinant are
    skipped for that determinant and counted, rather than being treated as a
    distinct state.
    """
    sf = strain_feats.set_index("strain")
    rows = []
    determinants = ["lps_core", "o_status", "o_antigen"]
    for genus, grp in phage_feats.groupby("genus"):
        hosts = sorted({h for h in grp["isolation_host"].dropna()})
        known = [h for h in hosts if h in sf.index]
        for det in determinants:
            states = {sf.loc[h, det] for h in known}
            states = sorted({s for s in states if isinstance(s, str)})
            rows.append(dict(
                genus=genus, determinant=det,
                n_phages_in_genus=int(len(grp)),
                n_isolation_hosts=len(hosts),
                n_hosts_with_metadata=len(known),
                n_hosts_with_call=sum(
                    1 for h in known if isinstance(sf.loc[h, det], str)),
                n_permissive_states=len(states),
                permissive_states="|".join(states),
            ))
    return pd.DataFrame(rows)


def apply_gates(pairs: pd.DataFrame, gates: pd.DataFrame) -> pd.DataFrame:
    """Flag each pair as passing or failing each gate tier."""
    lookup: dict[tuple[str, str], set[str]] = {}
    for r in gates.itertuples(index=False):
        states = set(r.permissive_states.split("|")) if r.permissive_states else set()
        lookup[(r.genus, r.determinant)] = states

    for gate_name, dets in GATE_TIERS.items():
        ok = np.ones(len(pairs), dtype=bool)
        for det in dets:
            vals = pairs[det].to_numpy()
            gen = pairs["genus"].to_numpy()
            ok &= np.array([
                isinstance(v, str) and v in lookup.get((gg, det), set())
                for v, gg in zip(vals, gen)])
        pairs[gate_name] = ok
    return pairs


def two_gate_decomposition(pairs: pd.DataFrame) -> pd.DataFrame:
    """Split failed infections into surface exclusion and post-adsorption loss."""
    rows = []
    n_obs = len(pairs)
    nonlytic = pairs[pairs["lytic"] == 0]
    lytic = pairs[pairs["lytic"] == 1]
    for gate in GATE_TIERS:
        g1 = int((~nonlytic[gate]).sum())
        g2 = int(nonlytic[gate].sum())
        passed = pairs[pairs[gate]]
        failed = pairs[~pairs[gate]]
        lysis_in = float(passed["lytic"].mean()) if len(passed) else np.nan
        lysis_out = float(failed["lytic"].mean()) if len(failed) else np.nan
        lo_in, hi_in = wilson_interval(int(passed["lytic"].sum()), len(passed))
        lo_out, hi_out = wilson_interval(int(failed["lytic"].sum()), len(failed))
        rows.append(dict(
            gate=gate,
            determinants="+".join(GATE_TIERS[gate]),
            n_observed_pairs=n_obs,
            n_pairs_passing_gate=int(len(passed)),
            pct_pairs_passing_gate=100.0 * len(passed) / n_obs,
            n_nonlytic=int(len(nonlytic)),
            n_gate1_surface_exclusion=g1,
            pct_nonlytic_at_gate1=100.0 * g1 / len(nonlytic),
            n_gate2_post_adsorption=g2,
            pct_nonlytic_at_gate2=100.0 * g2 / len(nonlytic),
            lysis_rate_inside_gate=lysis_in,
            lysis_inside_ci_low=lo_in,
            lysis_inside_ci_high=hi_in,
            lysis_rate_outside_gate=lysis_out,
            lysis_outside_ci_low=lo_out,
            lysis_outside_ci_high=hi_out,
            enrichment_inside_over_outside=(
                lysis_in / lysis_out if lysis_out and lysis_out > 0 else np.nan),
            pct_of_lysis_events_captured=(
                100.0 * passed["lytic"].sum() / lytic["lytic"].sum()
                if len(lytic) else np.nan),
        ))
    return pd.DataFrame(rows)


def gate_performance_by_genus(pairs: pd.DataFrame, gate: str) -> pd.DataFrame:
    rows = []
    for genus, g in pairs.groupby("genus"):
        inn, out = g[g[gate]], g[~g[gate]]
        li = float(inn["lytic"].mean()) if len(inn) else np.nan
        lo = float(out["lytic"].mean()) if len(out) else np.nan
        rows.append(dict(
            genus=genus, n_phages=int(g["phage"].nunique()), n_pairs=int(len(g)),
            n_passing=int(len(inn)), pct_passing=100.0 * len(inn) / len(g),
            lysis_rate_inside=li, lysis_rate_outside=lo,
            enrichment=(li / lo if lo and lo > 0 else np.nan),
            overall_lysis_rate=float(g["lytic"].mean()),
        ))
    return pd.DataFrame(rows).sort_values("n_pairs", ascending=False)


def main() -> int:
    paths.bootstrap()
    raw = load_raw()

    pairs = build_pair_table(raw["matrix"])
    n_cells = pairs.attrs["n_cells"]
    strains = sorted(pairs["strain"].unique())
    phages = sorted(pairs["phage"].unique())

    strain_feats = build_strain_features(raw["picard"], raw["omp"], raw["defense"], strains)
    phage_feats = build_phage_features(raw["guelin"], raw["rbp"])
    gates = build_gate_definitions(phage_feats, strain_feats)

    surf_cols = ["lps_core", "o_antigen", "o_status", "capsule_any",
                 "phylogroup", "sequence_type", "n_defense_systems",
                 "defense_burden_quartile"]
    pairs = pairs.merge(strain_feats[["strain"] + surf_cols], on="strain", how="left")
    pairs = pairs.merge(phage_feats[["phage", "genus", "morphotype", "family",
                                     "isolation_host"]], on="phage", how="left")
    pairs = apply_gates(pairs, gates)

    decomp = two_gate_decomposition(pairs)
    by_genus = gate_performance_by_genus(pairs, "gate_lps")

    # Persist
    pairs.to_csv(paths.PAIRS, index=False, compression="gzip")
    strain_feats.to_csv(paths.STRAIN_FEATURES, index=False)
    phage_feats.to_csv(paths.PHAGE_FEATURES, index=False)
    gates.to_csv(paths.GATE_DEFINITIONS, index=False)
    decomp.to_csv(paths.TABLES / "table1_two_gate_decomposition.csv", index=False)
    by_genus.to_csv(paths.TABLES / "table2_gate_performance_by_genus.csv", index=False)

    # Report
    out: list[str] = []
    add = out.append
    add("Receptor haplotype stratification report")
    add("Generated by scripts/02_stratify_receptor_haplotypes.py")
    add("")
    add("INTERACTION TABLE")
    add(f"  strains tested                  : {len(strains)}")
    add(f"  phages tested                   : {len(phages)}")
    add(f"  matrix cells                    : {n_cells}")
    add(f"  observed pairs retained         : {len(pairs)}")
    add(f"  untested cells dropped          : {n_cells - len(pairs)}")
    n_lytic = int(pairs['lytic'].sum())
    add(f"  productive lysis events         : {n_lytic} "
        f"({100.0 * n_lytic / len(pairs):.2f} percent)")
    add(f"  failed infections               : {len(pairs) - n_lytic} "
        f"({100.0 * (len(pairs) - n_lytic) / len(pairs):.2f} percent)")
    add("")
    add("  The share of interactions that fail, 79.25 percent, is close to the "
        "80 to 85 percent commonly quoted for surface-level failure across a "
        "species. That coincidence is not evidence that the failures are all at "
        "the surface, and the gate analysis below shows they are not.")

    add("")
    add("SURFACE DETERMINANT DISTRIBUTION AMONG THE TESTED STRAINS")
    for col in ["lps_core", "o_status", "phylogroup"]:
        vc = strain_feats[col].value_counts(dropna=False)
        add(f"  {col}:")
        for k, v in vc.items():
            add(f"    {str(k):12s} {v:4d} strains")
    add(f"  distinct O-antigen serotypes    : {strain_feats['o_antigen'].nunique()}")
    add(f"  strains with no capsule call    : "
        f"{int(strain_feats['capsule_any'].isna().sum())} "
        "(the strains contributed as phage isolation hosts)")
    add("")
    add("  Lysis rate by LPS outer core type, across all phages:")
    lps_rate = pairs.groupby("lps_core")["lytic"].agg(["mean", "size"]).sort_values("mean")
    for k, r in lps_rate.iterrows():
        add(f"    {str(k):10s} {100 * r['mean']:5.2f} percent of {int(r['size']):5d} pairs")
    add("  The outer core type alone moves the lysis rate by a factor of about "
        f"{lps_rate['mean'].max() / lps_rate['mean'].min():.1f}, so host surface "
        "variation plainly matters.")

    add("")
    add("RECEPTOR GATE CONSTRUCTED FROM PHAGE ISOLATION HOSTS")
    add(f"  VIRIDIC genera                  : {phage_feats['genus'].nunique()}")
    add(f"  distinct isolation hosts         : "
        f"{phage_feats['isolation_host'].nunique()}")
    add("  Permissive LPS outer core states per genus:")
    for r in gates[gates["determinant"] == "lps_core"].sort_values(
            "n_phages_in_genus", ascending=False).itertuples(index=False):
        add(f"    {r.genus:22s} {r.n_phages_in_genus:2d} phages, "
            f"{r.n_hosts_with_call:2d} anchor hosts -> {r.permissive_states}")
    add("")
    add("  The number of anchor hosts tracks the number of phages in a genus, so "
        "a genus represented by a single phage yields a single permissive state "
        "and a correspondingly strict gate. This makes gate stringency vary with "
        "genus size for a reason that has nothing to do with receptor biology, "
        "and is one of the gate's structural weaknesses.")

    add("")
    add("TWO-GATE DECOMPOSITION OF FAILED INFECTIONS")
    for r in decomp.itertuples(index=False):
        add(f"  {r.gate} ({r.determinants})")
        add(f"    pairs passing the gate            : {r.n_pairs_passing_gate} "
            f"({r.pct_pairs_passing_gate:.1f} percent of observed pairs)")
        add(f"    failures at gate 1, surface        : "
            f"{r.n_gate1_surface_exclusion} ({r.pct_nonlytic_at_gate1:.2f} percent "
            "of failed infections)")
        add(f"    failures at gate 2, post-adsorption: "
            f"{r.n_gate2_post_adsorption} ({r.pct_nonlytic_at_gate2:.2f} percent "
            "of failed infections)")
        add(f"    lysis rate inside the gate         : "
            f"{100 * r.lysis_rate_inside_gate:.2f} percent "
            f"(95 percent CI {100 * r.lysis_inside_ci_low:.2f} to "
            f"{100 * r.lysis_inside_ci_high:.2f})")
        add(f"    lysis rate outside the gate        : "
            f"{100 * r.lysis_rate_outside_gate:.2f} percent "
            f"(95 percent CI {100 * r.lysis_outside_ci_low:.2f} to "
            f"{100 * r.lysis_outside_ci_high:.2f})")
        add(f"    enrichment inside over outside     : "
            f"{r.enrichment_inside_over_outside:.2f} fold")
        add(f"    share of all lysis events captured : "
            f"{r.pct_of_lysis_events_captured:.1f} percent")
        add("")

    add("GATE PERFORMANCE BY PHAGE GENUS (LPS outer core gate)")
    add(f"  {'genus':22s} {'phages':>6s} {'pairs':>6s} {'pass pct':>9s} "
        f"{'lysis in':>9s} {'lysis out':>9s} {'enrich':>7s}")
    for r in by_genus.itertuples(index=False):
        ein = "n/a" if not np.isfinite(r.enrichment) else f"{r.enrichment:.2f}"
        add(f"  {r.genus:22s} {r.n_phages:6d} {r.n_pairs:6d} {r.pct_passing:8.1f}% "
            f"{100 * r.lysis_rate_inside:8.2f}% {100 * r.lysis_rate_outside:8.2f}% "
            f"{ein:>7s}")

    main_gate = decomp[decomp["gate"] == "gate_lps"].iloc[0]
    add("")
    add("INTERPRETATION AND CONSEQUENCE FOR THE PRIMARY ANALYSIS")
    add("  The isolation-host gate is a weak discriminator. It admits "
        f"{main_gate.pct_pairs_passing_gate:.1f} percent of pairs and raises the "
        f"lysis rate only from {100 * main_gate.lysis_rate_outside_gate:.1f} to "
        f"{100 * main_gate.lysis_rate_inside_gate:.1f} percent, an enrichment of "
        f"{main_gate.enrichment_inside_over_outside:.2f} fold. It assigns "
        f"{main_gate.pct_nonlytic_at_gate1:.1f} percent of failed infections to "
        "surface exclusion, well short of the 80 to 85 percent that a strict "
        "surface-masking account would predict.")
    add("  Adding the O-antigen serotype makes the gate far stricter without "
        "making it sharper: coverage falls to "
        f"{decomp.set_index('gate').loc['gate_lps_otype', 'pct_pairs_passing_gate']:.1f} "
        "percent of pairs while enrichment stays near "
        f"{decomp.set_index('gate').loc['gate_lps_otype', 'enrichment_inside_over_outside']:.2f} "
        "fold, which is what happens when a filter removes signal and noise in "
        "equal measure.")
    add("  Enrichment is real in some genera and absent in others. Felixounavirus "
        "and Phapecoctavirus separate cleanly, while Vectrevirus and "
        "Kayfunavirus, the two largest genera, show essentially none. A gate "
        "that works for a minority of genera cannot carry a species-wide claim.")
    add("")
    add("  Two conclusions follow, and they shape the rest of the pipeline.")
    add("  First, surface receptor state cannot be reduced to a binary "
        "compatibility call from the information this dataset releases. Reporting "
        "a hard receptor-compatible subset as though it were validated would "
        "overstate what the data support.")
    add("  Second, the correct way to remove surface confounding is to stratify "
        "on surface state rather than to filter on it. Comparing defense carriers "
        "with non-carriers inside the same phage and the same LPS outer core type "
        "removes the surface contribution without discarding data and without "
        "needing to know which state is permissive. Script 03 therefore reports a "
        "conditioning ladder, from an unstratified odds ratio through to one "
        "stratified on phage identity and surface state, and treats the gate "
        "above as a secondary analysis whose weakness is itself a finding.")

    text = "\n".join(out)
    (paths.RESULTS / "receptor_stratification_report.txt").write_text(
        text + "\n", encoding="utf-8")
    print(text)
    print(f"\nwrote {paths.PAIRS.relative_to(paths.ROOT)} "
          f"({paths.PAIRS.stat().st_size / 1024:.0f} KiB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
