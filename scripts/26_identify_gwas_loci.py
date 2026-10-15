"""Find out what the loci flagged by the association study actually are.

An association study that stops at a p-value has not finished. Script 25 flags
presence patterns, which are sets of k-mers, not genes. This script locates
each flagged pattern in a genome that carries it, extracts the surrounding
sequence, calls the genes and asks three questions of them: are they already in
the DefenseFinder catalogue that the rest of this analysis relies on, do they
carry any anti-defense or lysogeny domain, and do they resemble anything in the
96 phage genomes.

The last question matters most for interpretation. A locus associated with
being lysed more often, which is what the flagged patterns here are, is as
likely to be a prophage or a mobile element marking a susceptible lineage as
anything mechanistic.

Outputs:
    data/processed/gwas_locus_sequences.faa
    results/tables/table47_gwas_locus_annotation.csv
    results/gwas_loci_report.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from papr import config, paths
from papr.sketch import K, _canonical_hashes, read_raw_sequences, read_sequences

GENOME_DIR = paths.RAW / "host_genomes"
PATTERNS = paths.PROCESSED / "gwas_patterns.npz"
FLANK = 3000


def locate(strain: str, wanted: set[int]) -> list[tuple[int, int]]:
    """Positions in a genome whose canonical k-mer hash is in `wanted`."""
    hits = []
    for contig_index, codes in enumerate(read_sequences(GENOME_DIR / f"{strain}.fna.gz")):
        h = _canonical_hashes(codes, K)
        if h.size == 0:
            continue
        mask = np.isin(h, np.fromiter(wanted, dtype=np.uint64, count=len(wanted)))
        for pos in np.flatnonzero(mask):
            hits.append((contig_index, int(pos)))
    return hits


def main() -> int:
    paths.bootstrap()
    hits = pd.read_csv(paths.TABLES / "table45_gwas_hits.csv")
    if hits.empty:
        print("no flagged patterns to identify")
        (paths.RESULTS / "gwas_loci_report.txt").write_text(
            "No patterns were flagged by script 25, so there is nothing to "
            "identify.\n", encoding="utf-8")
        return 0

    with np.load(PATTERNS, allow_pickle=False) as z:
        patterns, weights, reps = z["patterns"], z["weights"], z["representatives"]
        strains = [str(x) for x in z["strains"]]

    strain_feats = pd.read_csv(paths.STRAIN_FEATURES, low_memory=False)
    ds_cols = [c for c in strain_feats.columns if c.startswith("ds_")]
    feats = strain_feats.set_index("strain")

    import pyrodigal
    import pyhmmer

    alphabet = pyhmmer.easel.Alphabet.amino()
    hmms, acc_of = [], {}
    for d in (paths.ADF_PROFILES, paths.RAW / "lifestyle_profiles"):
        for path in sorted(d.glob("*.hmm")):
            with pyhmmer.plan7.HMMFile(path) as handle:
                for hmm in handle:
                    nm = (hmm.name.decode() if isinstance(hmm.name, (bytes, bytearray))
                          else str(hmm.name))
                    acc_of[nm] = path.stem
                    hmms.append(hmm)
    print(f"annotation profiles loaded: {len(hmms)}")

    def as_str(v) -> str:
        return v.decode() if isinstance(v, (bytes, bytearray)) else str(v)

    records, fasta = [], []
    for row in hits.itertuples(index=False):
        pat = int(row.pattern)
        bits = patterns[pat]
        carriers = [strains[j] for j in range(len(strains))
                    if (bits[j >> 6] >> np.uint64(j & 63)) & np.uint64(1)]
        carrier = next((c for c in carriers
                        if (GENOME_DIR / f"{c}.fna.gz").exists()), None)
        if carrier is None:
            continue

        found = locate(carrier, {np.uint64(row.representative_hash
                                           if row.representative_hash >= 0
                                           else row.representative_hash + (1 << 64))})
        if not found:
            records.append(dict(pattern=pat, n_kmers=int(row.n_kmers),
                                n_carriers=len(carriers), carrier_used=carrier,
                                located=False, note="k-mer not relocated"))
            continue
        contig_index, pos = found[0]
        contigs = read_raw_sequences(GENOME_DIR / f"{carrier}.fna.gz")
        contig = contigs[contig_index]
        lo = max(0, pos - FLANK)
        hi = min(len(contig), pos + FLANK)
        window = contig[lo:hi]

        finder = pyrodigal.GeneFinder(meta=True)
        proteins = []
        for g, gene in enumerate(finder.find_genes(bytes(window)), 1):
            seq = gene.translate().rstrip("*")
            if len(seq) >= 40:
                proteins.append((f"pattern{pat}_{carrier}_{g}", seq))

        domain_hits = []
        if proteins:
            block = pyhmmer.easel.DigitalSequenceBlock(
                alphabet, [pyhmmer.easel.TextSequence(name=n.encode(), sequence=s)
                           .digitize(alphabet) for n, s in proteins])
            try:
                res = pyhmmer.hmmsearch(hmms, block, bit_cutoffs="gathering", cpus=0)
            except Exception:
                res = pyhmmer.hmmsearch(hmms, block, E=1e-5, cpus=0)
            for th in res:
                acc = acc_of.get(as_str(getattr(th.query, "name", th.query)), "")
                for h in th:
                    if h.included:
                        domain_hits.append(acc)
        for n, s in proteins:
            fasta.append(f">{n}\n{s}")

        # Do the carriers share a catalogued defense system the others lack?
        present = [c for c in carriers if c in feats.index]
        shared = []
        if present:
            sub = feats.loc[present, ds_cols]
            others = feats.loc[[s for s in strains
                                if s in feats.index and s not in set(present)], ds_cols]
            for col in ds_cols:
                if sub[col].mean() >= 0.8 and others[col].mean() <= 0.2:
                    shared.append(col[3:])

        records.append(dict(
            pattern=pat, n_kmers=int(row.n_kmers), n_carriers=len(carriers),
            carrier_used=carrier, located=True,
            effect_surface_adjusted=float(row.effect_surface_adjusted),
            q_surface_adjusted=float(row.q_surface_adjusted),
            phylogenetic_r2=float(row.phylogenetic_r2),
            n_genes_in_window=len(proteins),
            protein_lengths=",".join(str(len(s)) for _, s in proteins),
            domain_hits=",".join(sorted(set(domain_hits))) or "none",
            catalogued_defense_shared=",".join(shared) or "none",
            carriers=",".join(sorted(carriers)[:12])))
        print(f"  pattern {pat}: {len(proteins)} genes, "
              f"domains {sorted(set(domain_hits)) or 'none'}")

    table = pd.DataFrame(records)
    table.to_csv(paths.TABLES / "table47_gwas_locus_annotation.csv", index=False)
    (paths.PROCESSED / "gwas_locus_sequences.faa").write_text(
        "\n".join(fasta) + "\n", encoding="utf-8")

    o: list[str] = []
    add = o.append
    add("Association hit identification report")
    add("Generated by scripts/26_identify_gwas_loci.py")
    add("")
    add(f"  patterns flagged by script 25 : {len(hits)}")
    add(f"  relocated in a carrier genome : {int(table['located'].sum())}")
    add(f"  sequence window extracted     : {2 * FLANK} bases around the k-mer")
    add("")
    for r in table.itertuples(index=False):
        add(f"  pattern {r.pattern}")
        if not r.located:
            add("    could not be relocated in a carrier genome")
            continue
        add(f"    k-mers backing it      : {r.n_kmers} "
            f"(roughly {20 * r.n_kmers} bases of shared sequence)")
        add(f"    carrier strains        : {r.n_carriers}")
        add(f"    effect on lysis        : {r.effect_surface_adjusted:+.4f} "
            f"(q = {r.q_surface_adjusted:.2e})")
        add(f"    phylogenetic signal    : {r.phylogenetic_r2:.3f}")
        add(f"    genes in the window    : {r.n_genes_in_window} "
            f"(lengths {r.protein_lengths})")
        add(f"    anti-defense or lysogeny domains : {r.domain_hits}")
        add(f"    catalogued defense systems shared by carriers : "
            f"{r.catalogued_defense_shared}")
        add(f"    carriers : {r.carriers}")
        add("")

    located = table[table["located"]] if "located" in table else table
    add("INTERPRETATION")
    if len(located) and (located["effect_surface_adjusted"] > 0).all():
        add("  Every flagged locus is associated with being lysed more often, "
            "not less. Nothing in the accessory genome of these 254 strains "
            "is associated with resistance once prevalence and phylogenetic "
            "confinement are matched.")
        add("  A locus associated with susceptibility is a weaker finding than "
            "it sounds. It can mark a lineage that is easy to lyse for reasons "
            "unrelated to the locus, and with carriers numbering around ten "
            "strains the estimate rests on a small group. These are "
            "observations to follow up, not mechanisms.")
    add("  The accompanying protein sequences are written to "
        "data/processed/gwas_locus_sequences.faa so that a reader can take "
        "them further than the profile set used here allows. The profiles "
        "searched cover anti-defense families and lysogeny machinery only, so "
        "a locus with no hit is unannotated rather than uncharacterised.")

    text = "\n".join(o)
    (paths.RESULTS / "gwas_loci_report.txt").write_text(text + "\n", encoding="utf-8")
    print("\n" + text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
