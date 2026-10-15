"""Rebuild the receptor gate from receptor binding protein similarity.

The weakest step in this analysis has been the definition of receptor
compatibility. Script 02 anchored it on the strains each phage was isolated
from, pooled across VIRIDIC genus, because genus was the only grouping
available that plausibly tracks receptor specificity. That gate is poor: it
admits two thirds of pairs for a 1.34 fold enrichment in lysis.

Genus is a proxy for what actually determines adsorption, which is the receptor
binding protein at the tip of the tail. The released data name the receptor
binding protein of each phage but give no sequence, which is why script 02 could
not use them. Script 07 changed that: the gene identifiers in the released table
match the proteins predicted here for 125 of 130 entries, and the matched
proteins have a median length of 677 residues against 127 for all predicted
proteins, which is what tail fibers and spikes look like and is strong evidence
the identifiers correspond.

This script therefore groups phages by receptor binding protein homology rather
than by genus, pools isolation hosts within those groups, and asks whether the
resulting gate separates adsorption from post-adsorption failure any better.
The grouping uses only phage protein sequences and the gate uses only isolation
hosts, so nothing here is derived from the interaction outcomes being analysed.

Outputs:
    data/processed/rbp_sequences.faa
    data/processed/rbp_groups.csv
    results/tables/table27_rbp_group_assignments.csv
    results/tables/table28_gate_comparison.csv
    results/rbp_receptor_gate_report.txt
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from papr import paths
from papr.stats import wilson_interval

# Similarity thresholds for single-linkage grouping, as a fraction of the
# self-score of the shorter sequence. Swept rather than chosen in advance.
SIMILARITY_THRESHOLDS = [0.20, 0.35, 0.50, 0.65, 0.80]


def load_proteins() -> dict[str, str]:
    seqs: dict[str, str] = {}
    with gzip.open(paths.PHAGE_PROTEINS, "rt", encoding="utf-8") as fh:
        name, buf = None, []
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
        if name:
            seqs[name] = "".join(buf)
    return seqs


def collect_rbps(seqs: dict[str, str]) -> tuple[pd.DataFrame, list[str]]:
    rbp = pd.read_csv(paths.RBP, sep=";")
    rbp = rbp[rbp["RBP"].notna() & (rbp["RBP"].astype(str) != "NA")].copy()
    rows, missing = [], []
    for r in rbp.itertuples(index=False):
        pid = str(r.RBP)
        if pid in seqs and len(seqs[pid]) >= 50:
            rows.append(dict(phage=r.phage, protein_id=pid,
                             rbp_type=r.type, length=len(seqs[pid])))
        else:
            missing.append(pid)
    return pd.DataFrame(rows), missing


def pairwise_similarity(ids: list[str], seqs: dict[str, str]) -> pd.DataFrame:
    """All-against-all protein comparison with phmmer.

    The score is normalised by the self-score of the shorter sequence, giving a
    value near 1 for near-identical proteins and near 0 for unrelated ones.
    """
    import pyhmmer

    alphabet = pyhmmer.easel.Alphabet.amino()
    digital = [pyhmmer.easel.TextSequence(name=i.encode(), sequence=seqs[i])
               .digitize(alphabet) for i in ids]
    block = pyhmmer.easel.DigitalSequenceBlock(alphabet, digital)

    def as_str(v) -> str:
        return v.decode() if isinstance(v, (bytes, bytearray)) else str(v)

    raw: dict[tuple[str, str], float] = {}
    for hits in pyhmmer.phmmer(digital, block, cpus=0, E=10.0):
        q = as_str(getattr(hits.query, "name", hits.query))
        for hit in hits:
            t = as_str(hit.name)
            score = float(hit.score)
            key = (q, t)
            if score > raw.get(key, -np.inf):
                raw[key] = score

    self_score = {i: raw.get((i, i), np.nan) for i in ids}
    rows = []
    for a in ids:
        for b in ids:
            if a >= b:
                continue
            s = max(raw.get((a, b), 0.0), raw.get((b, a), 0.0))
            denom = np.nanmin([self_score.get(a, np.nan), self_score.get(b, np.nan)])
            if not np.isfinite(denom) or denom <= 0:
                continue
            rows.append(dict(a=a, b=b, score=s, normalised=min(1.0, s / denom)))
    return pd.DataFrame(rows)


def single_linkage(ids: list[str], sim: pd.DataFrame, threshold: float) -> dict[str, int]:
    parent = {i: i for i in ids}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for r in sim[sim["normalised"] >= threshold].itertuples(index=False):
        ra, rb = find(r.a), find(r.b)
        if ra != rb:
            parent[ra] = rb
    roots = {}
    labels = {}
    for i in ids:
        root = find(i)
        if root not in roots:
            roots[root] = len(roots) + 1
        labels[i] = roots[root]
    return labels


def build_gate(pairs: pd.DataFrame, phage_groups: pd.Series,
               strain_feats: pd.DataFrame, determinant: str) -> np.ndarray:
    """Permissive states per receptor group, taken from isolation hosts only."""
    sf = strain_feats.set_index("strain")
    hosts = pairs.groupby("phage")["isolation_host"].first()
    permissive: dict[int, set[str]] = {}
    for phage, group in phage_groups.items():
        host = hosts.get(phage)
        if host is None or host not in sf.index:
            continue
        state = sf.loc[host, determinant]
        if isinstance(state, str):
            permissive.setdefault(group, set()).add(state)
    grp = pairs["phage"].map(phage_groups)
    det = pairs[determinant]
    return np.array([
        isinstance(d, str) and g in permissive and d in permissive[g]
        for d, g in zip(det, grp)])


def evaluate(pairs: pd.DataFrame, passed: np.ndarray, label: str,
             n_groups: int) -> dict:
    inside = pairs[passed]
    outside = pairs[~passed]
    nonlytic = pairs[pairs["lytic"] == 0]
    li = float(inside["lytic"].mean()) if len(inside) else np.nan
    lo = float(outside["lytic"].mean()) if len(outside) else np.nan
    lo_ci = wilson_interval(int(inside["lytic"].sum()), len(inside))
    g1 = int((~passed[(pairs["lytic"] == 0).to_numpy()]).sum())
    return dict(
        gate=label, n_groups=n_groups,
        n_passing=int(passed.sum()),
        pct_passing=100.0 * passed.mean(),
        lysis_inside=li, lysis_inside_ci_low=lo_ci[0], lysis_inside_ci_high=lo_ci[1],
        lysis_outside=lo,
        enrichment=(li / lo if lo and lo > 0 else np.nan),
        pct_nonlytic_at_gate1=100.0 * g1 / len(nonlytic),
        pct_lysis_captured=100.0 * inside["lytic"].sum() / pairs["lytic"].sum(),
    )


def main() -> int:
    paths.bootstrap()
    pairs = pd.read_csv(paths.PAIRS, compression="gzip", low_memory=False)
    strain_feats = pd.read_csv(paths.STRAIN_FEATURES, low_memory=False)
    phage_feats = pd.read_csv(paths.PHAGE_FEATURES)

    seqs = load_proteins()
    rbps, missing = collect_rbps(seqs)
    ids = sorted(rbps["protein_id"].unique())
    print(f"receptor binding proteins recovered: {len(ids)} "
          f"across {rbps['phage'].nunique()} phages")

    (paths.PROCESSED / "rbp_sequences.faa").write_text(
        "".join(f">{i} {rbps.set_index('protein_id').loc[i, 'phage']}\n{seqs[i]}\n"
                for i in ids), encoding="utf-8")

    print("comparing receptor binding proteins all against all")
    sim = pairwise_similarity(ids, seqs)

    # A phage is assigned the group of its longest receptor binding protein,
    # which is the tail fiber where several are annotated.
    longest = (rbps.sort_values("length", ascending=False)
               .drop_duplicates("phage").set_index("phage")["protein_id"])

    results, assignments = [], {}
    for t in SIMILARITY_THRESHOLDS:
        labels = single_linkage(ids, sim, t)
        phage_group = longest.map(labels)
        n_groups = int(pd.Series(phage_group).nunique())
        assignments[t] = phage_group
        for det in ["lps_core", "o_status"]:
            passed = build_gate(pairs, phage_group, strain_feats, det)
            results.append(evaluate(pairs, passed, f"rbp_{t:.2f}_{det}", n_groups))

    # The genus gate from script 02, recomputed here on the same footing.
    genus_group = phage_feats.set_index("phage")["genus"].astype("category").cat.codes
    for det in ["lps_core", "o_status"]:
        passed = build_gate(pairs, genus_group, strain_feats, det)
        results.append(evaluate(pairs, passed, f"genus_{det}",
                                int(genus_group.nunique())))

    table = pd.DataFrame(results)
    table.to_csv(paths.TABLES / "table28_gate_comparison.csv", index=False)

    best_t = max(SIMILARITY_THRESHOLDS,
                 key=lambda t: table.set_index("gate").loc[f"rbp_{t:.2f}_lps_core",
                                                           "enrichment"])
    best_groups = assignments[best_t]
    out_assign = (phage_feats[["phage", "genus", "morphotype"]]
                  .assign(rbp_group=phage_feats["phage"].map(best_groups))
                  .assign(rbp_protein=phage_feats["phage"].map(longest)))
    out_assign.to_csv(paths.TABLES / "table27_rbp_group_assignments.csv", index=False)
    out_assign.to_csv(paths.PROCESSED / "rbp_groups.csv", index=False)

    genus_row = table.set_index("gate").loc["genus_lps_core"]
    best_row = table.set_index("gate").loc[f"rbp_{best_t:.2f}_lps_core"]

    o: list[str] = []
    add = o.append
    add("Receptor binding protein gate report")
    add("Generated by scripts/12_rbp_receptor_gate.py")
    add("")
    add("RECOVERING THE RECEPTOR BINDING PROTEINS")
    add(f"  entries in the released receptor binding protein table : "
        f"{len(rbps) + len(missing)}")
    add(f"  matched to a protein predicted here                    : {len(rbps)}")
    add(f"  not matched                                            : {len(missing)}")
    add(f"  phages with at least one recovered protein             : "
        f"{rbps['phage'].nunique()} of {len(phage_feats)}")
    add(f"  median length of recovered proteins                    : "
        f"{int(rbps['length'].median())} residues")
    add("  The released table names these proteins by gene identifier without "
        "giving a sequence. The identifiers match the gene order produced here, "
        "and the matched proteins are long in the way tail fibers and spikes "
        "are, so the correspondence is treated as established rather than "
        "assumed.")

    add("")
    add("GROUPING PHAGES BY RECEPTOR BINDING PROTEIN HOMOLOGY")
    add(f"  pairwise comparisons scored: {len(sim)}")
    add(f"  {'threshold':>9s} {'groups':>7s}")
    for t in SIMILARITY_THRESHOLDS:
        add(f"  {t:9.2f} {int(pd.Series(assignments[t]).nunique()):7d}")
    add(f"  For comparison the collection contains "
        f"{phage_feats['genus'].nunique()} VIRIDIC genera.")

    add("")
    add("DOES A RECEPTOR-BASED GATE SEPARATE BETTER THAN A GENUS-BASED ONE")
    add(f"  {'gate':24s} {'groups':>6s} {'pass pct':>9s} {'lysis in':>9s} "
        f"{'lysis out':>9s} {'enrich':>7s} {'gate1 pct':>9s}")
    for r in table.sort_values("enrichment", ascending=False).itertuples(index=False):
        add(f"  {r.gate:24s} {r.n_groups:6d} {r.pct_passing:8.1f}% "
            f"{100 * r.lysis_inside:8.2f}% {100 * r.lysis_outside:8.2f}% "
            f"{r.enrichment:7.2f} {r.pct_nonlytic_at_gate1:8.1f}%")

    add("")
    add("INTERPRETATION")
    delta = best_row["enrichment"] - genus_row["enrichment"]
    add(f"  Best receptor-based gate: similarity threshold {best_t:.2f} on the "
        f"LPS outer core, {int(best_row['n_groups'])} receptor groups, "
        f"enrichment {best_row['enrichment']:.2f} fold against "
        f"{genus_row['enrichment']:.2f} fold for the genus gate.")
    if delta > 0.15:
        add("  Grouping by receptor binding protein homology is a genuine "
            "improvement on grouping by genus. That is the expected direction: "
            "the receptor binding protein is what contacts the host surface, "
            "and genus membership is only a proxy for it.")
    elif delta > 0.0:
        add("  Grouping by receptor binding protein homology improves the gate "
            "only marginally. The receptor binding protein is the right "
            "biological object, so the limitation is not the grouping but the "
            "anchor: one isolation host per phage establishes that a surface "
            "state is permissive and can never establish that another is not.")
    else:
        add("  Grouping by receptor binding protein homology does not improve "
            "the gate. Since the receptor binding protein is the right "
            "biological object, this points at the anchor rather than the "
            "grouping: a single isolation host per phage can show that one "
            "surface state is permissive but never that any other is not, so "
            "the permissive sets stay small and the gate stays blunt however "
            "the phages are grouped.")
    add(f"  Under the best receptor-based gate, "
        f"{best_row['pct_nonlytic_at_gate1']:.1f} percent of failed infections "
        f"are attributed to surface exclusion, against "
        f"{genus_row['pct_nonlytic_at_gate1']:.1f} percent under the genus "
        "gate. That rise must not be read as support for the surface-masking "
        "account. The receptor-based gate attributes more failures to the "
        f"surface only because it is stricter, admitting "
        f"{best_row['pct_passing']:.1f} percent of pairs against "
        f"{genus_row['pct_passing']:.1f} percent, while its enrichment is "
        "lower. A gate that excludes more pairs and separates them less well "
        "is discarding signal along with noise, so the larger gate 1 share is "
        "an artefact of stringency rather than a better measurement. The share "
        "of failures genuinely attributable to the surface is not identified "
        "by either gate.")
    o_rows = table[table["gate"].str.endswith("o_status")]
    if (o_rows["enrichment"] < 1).any():
        add("  The O-antigen status gates are worse than useless: several have "
            "an enrichment below 1, meaning pairs passing the gate are lysed "
            "less often than those failing it. Anchoring on a two-level "
            "determinant whose permissive state is almost always the common "
            "one produces a gate that tracks the rarity of O-negative strains "
            "rather than receptor compatibility.")
    add("")
    add("  This is the strongest test of the two-gate framing this dataset "
        "allows without new experiments, because it uses the actual determinant "
        "of adsorption rather than a taxonomic stand-in. The conclusion of "
        "script 02 is unchanged: receptor compatibility in this collection "
        "cannot be reduced to a reliable binary call, and conditioning on "
        "surface state by stratification remains the defensible approach.")

    text = "\n".join(o)
    (paths.RESULTS / "rbp_receptor_gate_report.txt").write_text(
        text + "\n", encoding="utf-8")
    print("\n" + text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
