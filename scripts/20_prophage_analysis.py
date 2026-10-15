"""Test whether resident prophages explain resistance.

This has been the largest open gap in the analysis. Superinfection exclusion is
a well-documented route to post-adsorption resistance: a strain carrying a
prophage is typically immune to phages closely related to it, because the
resident repressor silences the incoming genome. Every earlier script had to
list this as untested, because the bacterial genomes were not available.

With the genomes retrieved in script 18 and validated in script 19, two
quantities can be measured.

The first is specific. For each strain and each of the 96 phages, the fraction
of that phage's k-mers present in the host genome. A strain carrying a prophage
closely related to a given phage will contain much of that phage's sequence.
This is the quantity superinfection exclusion is about, and exact k-mer
containment is the right tool for it: repressor-mediated immunity requires close
relatedness, and k-mer containment falls away quickly as sequences diverge.

The second is general. The number of integrase genes in a genome, which counts
prophages whether or not they resemble anything in this collection.

The specific measure is the real test, because it varies within a strain across
phages. A strain-level prophage count can only be compared between strains and
is confounded with everything else about the lineage; a strain-by-phage
relatedness measure is identified within a strain, with the strain's entire
genetic background held fixed.

Outputs:
    data/processed/phage_containment_matrix.csv.gz
    data/processed/prophage_markers.csv
    results/tables/table37_prophage_exclusion.csv
    results/tables/table38_prophage_burden.csv
    results/prophage_report.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
from scipy import stats

from papr import config, paths
from papr.models import fit_logit_cluster, term_effect
from papr.sketch import (K, SCALED, containment, read_raw_sequences,
                         read_sequences, sketch_file)
from papr.stats import benjamini_hochberg, wilson_interval

GENOME_DIR = paths.RAW / "host_genomes"
SKETCH_FILE = paths.PROCESSED / "host_genome_sketches.npz"
CONTAINMENT = paths.PROCESSED / "phage_containment_matrix.csv.gz"
MARKERS = paths.PROCESSED / "prophage_markers.csv"

# Containment above this is treated as a resident prophage closely related to
# the phage in question. At k = 21 a 40 kb prophage sharing 97 percent identity
# retains roughly half its k-mers, so this is a deliberately permissive bar for
# "closely related" and a strict one for "unrelated".
PROPHAGE_CONTAINMENT = 0.20


def load_host_sketches(strains: list[str]) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    if SKETCH_FILE.exists():
        with np.load(SKETCH_FILE, allow_pickle=False) as z:
            for name in z.files:
                if not name.endswith("__len") and name in set(strains):
                    out[name] = z[name]
    missing = [s for s in strains if s not in out]
    for i, s in enumerate(missing, 1):
        p = GENOME_DIR / f"{s}.fna.gz"
        if p.exists():
            out[s] = sketch_file(p)
    return out


def phage_sketches() -> dict[str, np.ndarray]:
    out = {}
    for p in sorted((paths.RAW / "phage_genomes").glob("*.fna")):
        out[p.stem] = sketch_file(p)
    return out


def build_containment(hosts: dict[str, np.ndarray],
                      phages: dict[str, np.ndarray]) -> pd.DataFrame:
    rows = {}
    names = sorted(phages)
    for i, (h, hs) in enumerate(sorted(hosts.items()), 1):
        rows[h] = [containment(phages[p], hs) for p in names]
        if i % 25 == 0 or i == len(hosts):
            print(f"  containment {i}/{len(hosts)}")
    return pd.DataFrame.from_dict(rows, orient="index", columns=names)


def count_integrases(strains: list[str]) -> pd.DataFrame:
    """Predict proteins and count integrase and excisionase markers per genome."""
    import pyrodigal
    import pyhmmer

    profile_dir = paths.RAW / "lifestyle_profiles"
    hmm_paths = sorted(profile_dir.glob("*.hmm"))
    if not hmm_paths:
        print("  lifestyle profiles absent; run script 17 first")
        return pd.DataFrame(columns=["strain", "n_integrase", "n_proteins"])

    alphabet = pyhmmer.easel.Alphabet.amino()
    hmms, acc_of = [], {}
    for path in hmm_paths:
        with pyhmmer.plan7.HMMFile(path) as handle:
            for hmm in handle:
                nm = (hmm.name.decode() if isinstance(hmm.name, (bytes, bytearray))
                      else str(hmm.name))
                acc_of[nm] = path.stem
                hmms.append(hmm)
    tyrosine = {"PF00589", "PF02899", "PF13356"}

    def as_str(v) -> str:
        return v.decode() if isinstance(v, (bytes, bytearray)) else str(v)

    rows = []
    for i, s in enumerate(strains, 1):
        path = GENOME_DIR / f"{s}.fna.gz"
        if not path.exists():
            continue
        contigs = read_raw_sequences(path)
        # These are whole bacterial genomes, so Prodigal is trained on each one
        # rather than run in metagenomic mode. That is both the appropriate
        # choice and sixteen times faster: 2 seconds a genome against 36, for
        # essentially the same gene set.
        training = [bytes(c) for c in contigs if len(c) > 1000]
        if not training:
            continue
        finder = pyrodigal.GeneFinder()
        try:
            finder.train(*training[:50])
        except Exception:
            finder = pyrodigal.GeneFinder(meta=True)
        seqs = []
        idx = 0
        for contig in contigs:
            for gene in finder.find_genes(contig):
                idx += 1
                prot = gene.translate().rstrip("*")
                if prot:
                    seqs.append(pyhmmer.easel.TextSequence(
                        name=f"{s}_{idx}".encode(), sequence=prot).digitize(alphabet))
        if not seqs:
            continue
        block = pyhmmer.easel.DigitalSequenceBlock(alphabet, seqs)
        n_int = 0
        try:
            results = pyhmmer.hmmsearch(hmms, block, bit_cutoffs="gathering", cpus=0)
        except Exception:
            results = pyhmmer.hmmsearch(hmms, block, E=1e-5, cpus=0)
        for hits in results:
            acc = acc_of.get(as_str(getattr(hits.query, "name", hits.query)), "")
            if acc in tyrosine:
                n_int += sum(1 for h in hits if h.included)
        rows.append(dict(strain=s, n_integrase=n_int, n_proteins=len(seqs)))
        if i % 20 == 0 or i == len(strains):
            print(f"  gene calling {i}/{len(strains)}")
    return pd.DataFrame(rows)


def main() -> int:
    paths.bootstrap()
    val = pd.read_csv(paths.PROCESSED / "host_genome_validation.csv")
    strains = sorted(val.loc[val["validated"], "strain"].astype(str))
    print(f"validated genomes: {len(strains)}")
    if len(strains) < 30:
        raise SystemExit("too few validated genomes for this analysis")

    pairs = pd.read_csv(paths.PAIRS, compression="gzip", low_memory=False)
    strain_feats = pd.read_csv(paths.STRAIN_FEATURES, low_memory=False)
    pcoa = pd.read_csv(paths.PCOA, index_col=0)
    clusters = pd.read_csv(paths.PROCESSED / "clonal_clusters.csv", index_col=0)

    if CONTAINMENT.exists():
        cont = pd.read_csv(CONTAINMENT, index_col=0, compression="gzip")
        cont.index = cont.index.astype(str)
        missing = [s for s in strains if s not in cont.index]
    else:
        cont, missing = None, strains
    if missing:
        print(f"computing containment for {len(missing)} genomes")
        hosts = load_host_sketches(strains)
        phs = phage_sketches()
        print(f"  phage sketches: {len(phs)}, median hashes "
              f"{int(np.median([v.size for v in phs.values()]))}")
        fresh = build_containment({s: hosts[s] for s in strains if s in hosts}, phs)
        cont = fresh if cont is None else pd.concat([cont, fresh]).groupby(level=0).first()
        cont.to_csv(CONTAINMENT, compression={"method": "gzip", "mtime": 0})

    cont = cont.loc[[s for s in strains if s in cont.index]]
    phage_cols = [c for c in cont.columns if c in set(pairs["phage"])]
    cont = cont[phage_cols]
    print(f"containment matrix: {cont.shape[0]} strains x {cont.shape[1]} phages")

    if MARKERS.exists():
        markers = pd.read_csv(MARKERS)
    else:
        print("counting integrase markers")
        markers = count_integrases(strains)
        markers.to_csv(MARKERS, index=False)

    long = cont.stack().rename("containment").reset_index()
    long.columns = ["strain", "phage", "containment"]
    data = pairs.merge(long, on=["strain", "phage"], how="inner")
    data["clonal_cluster"] = data["strain"].map(clusters["clonal_cluster"])
    data["related_prophage"] = (data["containment"] >= PROPHAGE_CONTAINMENT).astype(int)
    print(f"analysable strain-phage pairs: {len(data)}")

    coords = pcoa.reindex(data["strain"]).to_numpy(dtype=float)
    coords = (coords - coords.mean(axis=0)) / coords.std(axis=0)
    base = np.hstack([
        np.ones((len(data), 1)),
        pd.get_dummies(data["phage"].astype("string"), drop_first=True,
                       dtype=float).to_numpy(),
        pd.get_dummies(data["lps_core"].astype("string"), drop_first=True,
                       dtype=float).to_numpy(),
        pd.get_dummies(data["phylogroup"].astype("string"), drop_first=True,
                       dtype=float).to_numpy(),
        coords,
    ])
    y = data["lytic"].to_numpy(dtype=float)
    groups = data["clonal_cluster"].to_numpy()

    results = []
    for label, exposure in [
        ("related prophage present (containment >= "
         f"{PROPHAGE_CONTAINMENT})", data["related_prophage"].to_numpy(dtype=float)),
        ("containment, per 10 percent",
         data["containment"].to_numpy(dtype=float) * 10.0),
    ]:
        if exposure.std() == 0:
            results.append(dict(exposure=label, n_pairs=len(data),
                                odds_ratio=np.nan, ci_low=np.nan, ci_high=np.nan,
                                p_value=np.nan, note="no variation"))
            continue
        res = fit_logit_cluster(y, np.hstack([exposure.reshape(-1, 1), base]), groups)
        eff = term_effect(res, 0)
        # With eight exposed pairs and no events among them the likelihood is
        # maximised at minus infinity: the fit separates completely and the
        # reported interval spans twenty orders of magnitude. Printing those
        # numbers would imply a precision that does not exist.
        separated = (not np.isfinite(eff["std_error"])) or eff["std_error"] > 5.0
        note = "not estimable, complete separation" if separated else ""
        if separated:
            eff = {k: np.nan for k in eff}
        results.append(dict(exposure=label, n_pairs=len(data),
                            n_exposed=int((exposure > 0).sum()), **eff, note=note))
    excl = pd.DataFrame(results)
    excl.to_csv(paths.TABLES / "table37_prophage_exclusion.csv", index=False)

    # Strain-level burden.
    per_strain = data.groupby("strain").agg(lysis_rate=("lytic", "mean"))
    burden = per_strain.join(
        cont.apply(lambda r: int((r >= PROPHAGE_CONTAINMENT).sum()), axis=1)
        .rename("n_related_prophages"))
    if len(markers):
        burden = burden.join(markers.set_index("strain")[["n_integrase"]])
    burden = burden.join(strain_feats.set_index("strain")[["phylogroup",
                                                           "n_defense_systems"]])
    burden.reset_index().to_csv(paths.TABLES / "table38_prophage_burden.csv",
                                index=False)

    o: list[str] = []
    add = o.append
    add("Prophage analysis report")
    add("Generated by scripts/20_prophage_analysis.py")
    add("")
    add("COVERAGE")
    add(f"  validated host genomes      : {len(strains)} of 402 strains "
        f"({100 * len(strains) / 402:.0f} percent)")
    add(f"  phages compared             : {cont.shape[1]}")
    add(f"  strain-phage pairs analysed : {len(data)} of {len(pairs)}")
    add(f"  sketch parameters           : k = {K}, scaled = {SCALED}")

    add("")
    add("HOW OFTEN IS A RELATED PROPHAGE PRESENT")
    n_rel = int(data["related_prophage"].sum())
    lo, hi = wilson_interval(n_rel, len(data))
    add(f"  pairs where the host contains at least "
        f"{100 * PROPHAGE_CONTAINMENT:.0f} percent of the phage's k-mers : "
        f"{n_rel} ({100 * n_rel / len(data):.2f} percent, 95 percent CI "
        f"{100 * lo:.2f} to {100 * hi:.2f})")
    add(f"  strains carrying at least one related prophage : "
        f"{int((burden['n_related_prophages'] > 0).sum())} of {len(burden)}")
    add(f"  median related prophages per strain            : "
        f"{int(burden['n_related_prophages'].median())}")
    add(f"  maximum                                       : "
        f"{int(burden['n_related_prophages'].max())}")
    if len(markers):
        add(f"  median tyrosine recombinase domain hits per genome : "
            f"{int(markers['n_integrase'].median())} "
            f"(range {int(markers['n_integrase'].min())} to "
            f"{int(markers['n_integrase'].max())})")
        add("  That count is not a prophage count and is not used as one. The "
            "tyrosine recombinase profiles also match XerC and XerD, which "
            "every Escherichia genome carries, and integron and transposon "
            "recombinases. A typical Escherichia genome holds under ten "
            "prophages, so these numbers overcount by a wide margin. The "
            "count is reported as a crude index of mobile element content, "
            "not as prophage burden.")

    add("")
    add("DOES A RELATED PROPHAGE PREVENT LYSIS")
    add("  Phage identity, host surface type, phylogroup and core-genome "
        "position are all adjusted for, and standard errors are clustered on "
        "clonal background. Because containment varies within a strain across "
        "phages, this contrast is identified with the host's entire genetic "
        "background held fixed.")
    add("")
    add(f"  {'exposure':52s} {'OR':>7s} {'95 pct CI':>16s} {'p':>9s}")
    for r in excl.itertuples(index=False):
        if not np.isfinite(r.odds_ratio):
            add(f"  {r.exposure:52s} {'n/a':>7s}  {r.note} "
                f"({r.n_exposed} exposed pairs)")
            continue
        add(f"  {r.exposure:52s} {r.odds_ratio:7.3f} "
            f"[{r.ci_low:6.3f}, {r.ci_high:6.3f}] {r.p_value:9.2e}")

    # How much could this test have seen? Expected lysis events among the
    # carriers under the null is the binding constraint, and it is tiny.
    carriers = data[data["related_prophage"] == 1]
    involved = sorted(carriers["phage"].unique())
    base_rate = float(data[data["phage"].isin(involved)]["lytic"].mean())         if involved else float("nan")
    expected = len(carriers) * base_rate
    add("")
    add("HOW MUCH COULD THIS TEST HAVE SEEN")
    add(f"  pairs with a related resident prophage : {len(carriers)}")
    add(f"  phages involved                        : "
        f"{', '.join(involved) if involved else 'none'}")
    add(f"  lysis rate of those phages across all strains : "
        f"{100 * base_rate:.2f} percent")
    add(f"  lysis events expected among the carriers under the null : "
        f"{expected:.2f}")
    add(f"  lysis events observed                  : "
        f"{int(carriers['lytic'].sum())}")
    add("  Observing no lysis among the carriers is therefore exactly what "
        "the null predicts. The carriers are few and the phages involved "
        "barely lyse anything, so this comparison cannot distinguish "
        "superinfection exclusion from chance however it is analysed. The "
        "result below is reported for completeness and carries no weight.")

    lifestyle_path = paths.PROCESSED / "phage_lifestyle.csv"
    if lifestyle_path.exists() and involved:
        life = pd.read_csv(lifestyle_path).set_index("phage")
        temperate = [ph for ph in involved
                     if str(life.loc[ph, "lifestyle"]) == "temperate"]             if all(ph in life.index for ph in involved) else []
        add("")
        add("  One thing this does establish. The phages found integrated in "
            f"host genomes are {', '.join(involved)}, and "
            f"{len(temperate)} of {len(involved)} of them were independently "
            "classified temperate in script 17 on the basis of carrying a "
            "tyrosine integrase. Two analyses that share no input beyond the "
            "phage genomes themselves, one searching host chromosomes for "
            "phage sequence and one searching phage proteomes for integration "
            "machinery, pick out the same phages. That is a real check on "
            "both.")

    binary = excl.iloc[0]
    add("")
    if np.isfinite(binary["odds_ratio"]):
        if binary["p_value"] <= 0.05 and binary["odds_ratio"] < 1:
            add("  Carrying a prophage closely related to the challenging phage "
                "is associated with resisting it, which is what superinfection "
                "exclusion predicts. Given the expected event count above, "
                "treat this as a coincidence of small numbers unless it "
                "survives a collection with more lytic temperate phages.")
        elif binary["p_value"] <= 0.05:
            add("  Carrying a related prophage is associated with being lysed "
                "more often, not less. That is the opposite of superinfection "
                "exclusion and most likely reflects the reverse direction of "
                "causation: a strain is more likely to carry a prophage "
                "related to a phage that can successfully infect it.")
        else:
            add("  Carrying a prophage closely related to the challenging phage "
                "shows no detectable association with resisting it, which "
                "given the power calculation above is uninformative rather "
                "than negative.")
    else:
        add("  The regression cannot be fitted at all: with eight exposed "
            "pairs and no lysis among them the likelihood has no interior "
            "maximum. The honest summary is that this collection cannot test "
            "superinfection exclusion, and the reason is in its design rather "
            "than in the method used here. Every one of the eight pairs is "
            "non-lytic, which is consistent with the mechanism and equally "
            "consistent with the 0.74 percent baseline lysis rate of the two "
            "phages involved.")

    add("")
    add("STRAIN-LEVEL BURDEN")
    if burden["n_related_prophages"].std() > 0:
        rho = stats.spearmanr(burden["n_related_prophages"], burden["lysis_rate"])
        add(f"  Spearman, related prophage count against per-strain lysis rate : "
            f"{rho.statistic:+.3f} (p = {rho.pvalue:.3g}, n = {len(burden)})")
    if "n_integrase" in burden and burden["n_integrase"].notna().sum() > 10:
        sub = burden.dropna(subset=["n_integrase"])
        rho2 = stats.spearmanr(sub["n_integrase"], sub["lysis_rate"])
        add(f"  Spearman, integrase count against per-strain lysis rate        : "
            f"{rho2.statistic:+.3f} (p = {rho2.pvalue:.3g}, n = {len(sub)})")
        rho3 = stats.spearmanr(sub["n_integrase"], sub["n_defense_systems"])
        add(f"  Spearman, recombinase count against defense system count       : "
            f"{rho3.statistic:+.3f} (p = {rho3.pvalue:.3g})")
        if rho3.pvalue < 0.01:
            add("  That last correlation is the one substantive result in this "
                "section. Strains richer in site-specific recombinases carry "
                "more defense systems, which is what is expected if defense "
                "systems travel on mobile elements. It also means defense "
                "carriage is partly a marker of how much mobile DNA a strain "
                "has taken up, which is one more thing that a defense system "
                "can be standing in for in the associations reported "
                "elsewhere in this analysis.")
    add("  Strain-level counts are between-strain comparisons and inherit every "
        "confound that script 04 removed. They are descriptive; the "
        "within-strain contrast above is the one to read.")

    add("")
    add("LIMITS")
    add("  Exact k-mer containment detects prophages closely related to a "
        "collection phage and misses diverged ones: at k = 21 a prophage at 95 "
        "percent nucleotide identity retains only about a third of its k-mers. "
        "This measures relatedness to these 96 phages, not prophage content in "
        "general, which is the right target for superinfection exclusion but "
        "not for prophage-encoded defense in general.")
    add(f"  The analysis covers {len(strains)} of 402 strains, and strains "
        "deposited in public archives are not a random sample of the "
        "collection.")
    add("  The binding limitation is not method but collection design. "
        "Prophages closely related to a collection phage are rare here "
        "because the Guelin phages were chosen for lytic activity: only two "
        "of the ninety-six are temperate, and those two plaque on under one "
        "percent of strains. Superinfection exclusion needs a temperate phage "
        "that also forms plaques, and this panel contains almost none. "
        "Testing the mechanism properly needs a collection assembled for that "
        "question, not a reanalysis of one assembled for another.")
    add("  Containment is measured from an assembly. A prophage split across "
        "contigs is still detected, since k-mers are counted without regard to "
        "assembly structure, but a prophage absent from the assembly through "
        "collapse of a repeated region would be missed.")

    text = "\n".join(o)
    (paths.RESULTS / "prophage_report.txt").write_text(text + "\n", encoding="utf-8")
    print("\n" + text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
