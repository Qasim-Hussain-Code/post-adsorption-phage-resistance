"""Annotate the 96 phage genomes for anti-defense systems directly.

Script 05 had to impute anti-defense content from the prevalence of each family
among a genus's sequenced relatives, which is a lineage-level prior rather than
a measurement and attenuates any real interaction toward the null. That was the
largest single weakness in the analysis. This script removes it by annotating
each genome itself.

The annotation reproduces the AntiDefenseFinder procedure with pure-Python
implementations of the two external tools it depends on:

  gene calling      pyrodigal, a binding of the Prodigal algorithm, run in
                    metagenomic mode because phage genomes are short and
                    compositionally unlike their hosts
  profile search    pyhmmer, a binding of HMMER3, searching the 178 curated
                    anti-defense profiles shipped in defense-finder-models

Thresholds are the curated gathering cutoffs carried inside the profiles
themselves, which is what AntiDefenseFinder uses, so hit calling is not
reinvented here. System calls then apply each model's own rules from its XML
definition: the number of mandatory genes required and the minimum number of
genes required. Of the 156 anti-defense definitions, 150 require a single gene,
so for almost all families a profile hit above the gathering cutoff is the
system call. For the few that require two or more genes, the model's
inter_gene_max_space constraint is applied over Prodigal gene order, as
MacSyFinder does.

The result is validated against the published survey: per-genus prevalence
computed from these 96 genomes is compared with prevalence among RefSeq
relatives of the same genera.

Outputs:
    data/processed/phage_antidefense_hits.csv
    data/processed/phage_antidefense_systems.csv
    results/tables/table14_phage_antidefense_calls.csv
    results/tables/table15_antidefense_direct_vs_imputed.csv
    results/phage_antidefense_report.txt
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
import tarfile
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from papr import paths
from papr.stats import wilson_interval

COLI = "https://raw.githubusercontent.com/mdmparis/coli_phage_interactions_2023/main/"
MODELS_URL = ("https://github.com/mdmparis/defense-finder-models/releases/download/"
              "3.1.0/defense-finder-models-3.1.0.tar.gz")
MODELS_VERSION = "3.1.0"


def fetch_phage_genomes(phages: list[str]) -> list[Path]:
    """Download one FASTA per phage, skipping any already present."""
    paths.PHAGE_GENOMES.mkdir(parents=True, exist_ok=True)
    out = []
    for name in phages:
        dest = paths.PHAGE_GENOMES / f"{name}.fna"
        if not dest.exists() or dest.stat().st_size == 0:
            url = COLI + f"data/genomics/phages/FNA/{name}.fna"
            req = urllib.request.Request(
                url, headers={"User-Agent": "post-adsorption-analysis/1.0"})
            with urllib.request.urlopen(req, timeout=120) as resp:
                dest.write_bytes(resp.read())
        out.append(dest)
    return out


def fetch_and_extract_models(workdir: Path) -> tuple[int, int]:
    """Pull the model release and keep only the anti-defense parts.

    The archive is about 50 MB and almost all of it is bacterial defense
    profiles this analysis does not need. Only the anti-defense profiles and
    definitions are written to disk and the archive is deleted immediately
    afterwards, so the project's footprint grows by a few megabytes rather than
    by fifty.
    """
    if paths.ADF_PROFILES.exists() and any(paths.ADF_PROFILES.glob("*.hmm")):
        n_prof = len(list(paths.ADF_PROFILES.glob("*.hmm")))
        n_def = len(list(paths.ADF_DEFINITIONS.rglob("*.xml")))
        return n_prof, n_def

    workdir.mkdir(parents=True, exist_ok=True)
    tarball = workdir / "defense-finder-models.tar.gz"
    if not tarball.exists():
        # Download to a temporary name and rename only on success, so an
        # interrupted run can never leave a truncated archive that the next run
        # mistakes for a complete one.
        partial = tarball.with_suffix(tarball.suffix + ".part")
        req = urllib.request.Request(
            MODELS_URL, headers={"User-Agent": "post-adsorption-analysis/1.0"})
        with urllib.request.urlopen(req, timeout=600) as resp, partial.open("wb") as fh:
            while True:
                chunk = resp.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
        partial.replace(tarball)

    paths.ADF_PROFILES.mkdir(parents=True, exist_ok=True)
    paths.ADF_DEFINITIONS.mkdir(parents=True, exist_ok=True)

    # The archive is gzip-compressed, so random access by member name forces a
    # fresh decompression each time and turns a few hundred extractions into
    # minutes of work. Both passes below stream through the archive in order
    # instead, which costs two decompressions in total.
    wanted_genes: set[str] = set()
    n_defs = 0
    with tarfile.open(tarball, "r|gz") as tar:
        for member in tar:
            name = member.name
            if not (member.isfile() and "/ADF/" in name and name.endswith(".xml")):
                continue
            handle = tar.extractfile(member)
            if handle is None:
                continue
            text = handle.read().decode("utf-8", "replace")
            wanted_genes.update(re.findall(r'name="([^"]+)"', text))
            family = name.split("/ADF/")[1].split("/")[0]
            dest_dir = paths.ADF_DEFINITIONS / family
            dest_dir.mkdir(parents=True, exist_ok=True)
            (dest_dir / Path(name).name).write_text(text, encoding="utf-8")
            n_defs += 1

    n_prof = 0
    with tarfile.open(tarball, "r|gz") as tar:
        for member in tar:
            name = member.name
            if not (member.isfile() and name.endswith(".hmm")):
                continue
            gene = Path(name).stem
            if gene not in wanted_genes:
                continue
            handle = tar.extractfile(member)
            if handle is None:
                continue
            (paths.ADF_PROFILES / f"{gene}.hmm").write_bytes(handle.read())
            n_prof += 1

    tarball.unlink(missing_ok=True)
    return n_prof, n_defs


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records, name, chunks = [], None, []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def call_genes(genome_files: list[Path]) -> tuple[dict[str, list], int]:
    """Predict proteins per phage with Prodigal in metagenomic mode."""
    import pyrodigal

    finder = pyrodigal.GeneFinder(meta=True)
    proteins: dict[str, list] = {}
    total = 0
    for path in genome_files:
        phage = path.stem
        entries = []
        index = 0
        for contig, sequence in read_fasta(path):
            for gene in finder.find_genes(sequence.encode()):
                index += 1
                entries.append({
                    "protein_id": f"{phage}_{index:05d}",
                    "contig": contig,
                    "order": index,
                    "begin": int(gene.begin),
                    "end": int(gene.end),
                    "strand": int(gene.strand),
                    "sequence": gene.translate().rstrip("*"),
                })
        proteins[phage] = entries
        total += len(entries)
    return proteins, total


def search_profiles(proteins: dict[str, list]) -> pd.DataFrame:
    """Search every predicted protein against the anti-defense profiles.

    Hits are kept at the profile's own gathering cutoff, the curated
    per-profile bit score threshold that AntiDefenseFinder relies on, rather
    than at an arbitrary E-value.
    """
    import pyhmmer

    alphabet = pyhmmer.easel.Alphabet.amino()
    sequences = []
    owner: dict[str, str] = {}
    for phage, entries in proteins.items():
        for entry in entries:
            if not entry["sequence"]:
                continue
            name = entry["protein_id"].encode()
            owner[entry["protein_id"]] = phage
            sequences.append(
                pyhmmer.easel.TextSequence(
                    name=name, sequence=entry["sequence"]).digitize(alphabet))
    block = pyhmmer.easel.DigitalSequenceBlock(alphabet, sequences)

    hmms = []
    for path in sorted(paths.ADF_PROFILES.glob("*.hmm")):
        with pyhmmer.plan7.HMMFile(path) as handle:
            for hmm in handle:
                hmms.append(hmm)

    def as_str(value) -> str:
        """pyhmmer returns names as bytes in some versions and str in others."""
        return value.decode() if isinstance(value, (bytes, bytearray)) else str(value)

    rows = []
    for top_hits in pyhmmer.hmmsearch(hmms, block, bit_cutoffs="gathering", cpus=0):
        query = top_hits.query
        gene = as_str(getattr(query, "name", query))
        for hit in top_hits:
            if not hit.included:
                continue
            protein_id = as_str(hit.name)
            rows.append({
                "phage": owner[protein_id],
                "protein_id": protein_id,
                "gene": gene,
                "bit_score": float(hit.score),
                "evalue": float(hit.evalue),
            })
    return pd.DataFrame(rows)


def load_definitions() -> list[dict]:
    """Read each anti-defense model's gene requirements from its XML."""
    models = []
    for path in sorted(paths.ADF_DEFINITIONS.rglob("*.xml")):
        root = ET.parse(path).getroot()
        genes = {g.get("name"): g.get("presence") for g in root.findall(".//gene")}
        mandatory = [n for n, p in genes.items() if p == "mandatory"]
        models.append({
            "family": path.parent.name,
            "system": path.stem,
            "genes": set(genes),
            "mandatory": set(mandatory),
            "min_mandatory": int(root.get("min_mandatory_genes_required", 1)),
            "min_genes": int(root.get("min_genes_required", 1)),
            "max_space": int(root.get("inter_gene_max_space", 5)),
        })
    return models


def call_systems(hits: pd.DataFrame, proteins: dict[str, list],
                 models: list[dict], phages: list[str]) -> pd.DataFrame:
    """Apply each model's own rules to the profile hits."""
    order = {}
    for phage, entries in proteins.items():
        for entry in entries:
            order[entry["protein_id"]] = entry["order"]

    by_phage: dict[str, pd.DataFrame] = {
        p: g for p, g in hits.groupby("phage")} if len(hits) else {}

    rows = []
    for phage in phages:
        frame = by_phage.get(phage)
        present_genes = set(frame["gene"]) if frame is not None else set()
        for model in models:
            found = present_genes & model["genes"]
            found_mandatory = present_genes & model["mandatory"]
            ok = (len(found_mandatory) >= model["min_mandatory"]
                  and len(found) >= model["min_genes"])
            if ok and model["min_genes"] > 1 and frame is not None:
                # Co-localisation only binds when more than one gene is
                # required, which is true for six of the 156 models.
                positions = sorted(order[p] for p in
                                   frame[frame["gene"].isin(found)]["protein_id"])
                clustered, run = 1, 1
                for a, b in zip(positions, positions[1:]):
                    run = run + 1 if (b - a - 1) <= model["max_space"] else 1
                    clustered = max(clustered, run)
                ok = clustered >= model["min_genes"]
            if ok:
                rows.append({
                    "phage": phage, "family": model["family"],
                    "system": model["system"],
                    "genes_found": ",".join(sorted(found)),
                    "n_genes_found": len(found),
                })
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", default=None,
                    help="scratch directory for the model archive")
    args = ap.parse_args()

    paths.bootstrap()
    phage_feats = pd.read_csv(paths.PHAGE_FEATURES)
    phages = sorted(phage_feats["phage"].astype(str))

    print(f"fetching {len(phages)} phage genomes")
    genome_files = fetch_phage_genomes(phages)
    genome_mb = sum(p.stat().st_size for p in genome_files) / 1024 / 1024

    workdir = Path(args.workdir) if args.workdir else (paths.ROOT / ".model_cache")
    print("preparing anti-defense models")
    n_profiles, n_defs = fetch_and_extract_models(workdir)

    print("calling genes")
    proteins, n_proteins = call_genes(genome_files)

    print("searching anti-defense profiles")
    hits = search_profiles(proteins)

    models = load_definitions()
    systems = call_systems(hits, proteins, models, phages)

    # Persist the proteome so the annotation can be re-checked without rerunning.
    with gzip.open(paths.PHAGE_PROTEINS, "wt", encoding="utf-8", compresslevel=9) as fh:
        for phage in phages:
            for entry in proteins.get(phage, []):
                fh.write(f">{entry['protein_id']} {phage} "
                         f"{entry['begin']}..{entry['end']}:{entry['strand']}\n")
                fh.write(entry["sequence"] + "\n")
    hits.to_csv(paths.PHAGE_ADF_HITS, index=False)

    # Wide presence matrix, one row per phage, one column per family.
    families = sorted({m["family"] for m in models})
    wide = pd.DataFrame(0, index=pd.Index(phages, name="phage"), columns=families)
    if not systems.empty:
        for phage, family in zip(systems["phage"], systems["family"]):
            wide.loc[phage, family] = 1
    wide["n_antidefense_systems"] = (
        systems.groupby("phage").size().reindex(phages).fillna(0).astype(int).to_numpy()
        if not systems.empty else 0)
    wide["any_antidefense"] = (wide[families].sum(axis=1) > 0).astype(int)
    wide.reset_index().to_csv(paths.PHAGE_ADF_SYSTEMS, index=False)
    wide.reset_index().to_csv(
        paths.TABLES / "table14_phage_antidefense_calls.csv", index=False)

    # Validation against the published RefSeq survey, per genus.
    genus = phage_feats.set_index("phage")["genus"]
    imputed = pd.read_csv(paths.ADF_GENUS, index_col=0)
    direct_rows = []
    for g, members in genus.groupby(genus):
        names = [p for p in members.index if p in wide.index]
        sub = wide.loc[names]
        n = len(names)
        k = int(sub["any_antidefense"].sum())
        lo, hi = wilson_interval(k, n)
        row = {"genus": g, "n_phages_here": n, "n_with_any_direct": k,
               "prevalence_direct": k / n if n else np.nan,
               "direct_ci_low": lo, "direct_ci_high": hi}
        if g in imputed.index:
            row["prevalence_refseq_imputed"] = float(imputed.loc[g, "prevalence_any"])
            row["n_refseq_surveyed"] = int(imputed.loc[g, "n_surveyed"])
            for fam in ["Anti_RM", "Anti_CRISPR", "Anti_CBASS", "Anti_TA"]:
                row[f"direct_{fam}"] = (float(sub[fam].mean())
                                        if fam in sub.columns else np.nan)
                row[f"refseq_{fam}"] = (float(imputed.loc[g, fam])
                                        if fam in imputed.columns else np.nan)
        direct_rows.append(row)
    comparison = pd.DataFrame(direct_rows).sort_values("n_phages_here", ascending=False)
    comparison.to_csv(
        paths.TABLES / "table15_antidefense_direct_vs_imputed.csv", index=False)

    out: list[str] = []
    add = out.append
    add("Direct anti-defense annotation of the 96 phage genomes")
    add("Generated by scripts/07_annotate_phage_antidefense.py")
    add("")
    add("PROCEDURE")
    add(f"  genomes annotated                 : {len(genome_files)} "
        f"({genome_mb:.1f} MiB of FASTA)")
    add(f"  anti-defense profiles searched    : {n_profiles}")
    add(f"  anti-defense model definitions    : {n_defs}")
    add(f"  defense-finder-models version     : {MODELS_VERSION}")
    add(f"  proteins predicted by Prodigal    : {n_proteins}")
    add(f"  profile hits above gathering cutoff: {len(hits)}")
    add(f"  system calls after applying rules : {len(systems)}")
    add("  Hit thresholds are the curated gathering cutoffs carried in the "
        "profiles, not a chosen E-value, so the calls follow the same criterion "
        "AntiDefenseFinder applies.")

    add("")
    add("WHAT THE 96 PHAGES CARRY")
    add(f"  phages with at least one anti-defense system: "
        f"{int(wide['any_antidefense'].sum())} of {len(phages)} "
        f"({100 * wide['any_antidefense'].mean():.1f} percent)")
    add(f"  median systems per phage          : "
        f"{int(wide['n_antidefense_systems'].median())}")
    add(f"  maximum systems in one phage      : "
        f"{int(wide['n_antidefense_systems'].max())}")
    add("")
    add("  Prevalence by family across the 96 phages:")
    for fam in families:
        n = int(wide[fam].sum())
        if n == 0:
            continue
        lo, hi = wilson_interval(n, len(phages))
        add(f"    {fam:16s} {n:3d} phages ({100 * n / len(phages):5.1f} percent, "
            f"95 percent CI {100 * lo:4.1f} to {100 * hi:4.1f})")
    absent = [f for f in families if int(wide[f].sum()) == 0]
    add(f"    absent entirely: {', '.join(absent) if absent else 'none'}")

    add("")
    add("VALIDATION AGAINST THE PUBLISHED REFSEQ SURVEY")
    add(f"  {'genus':22s} {'phages':>6s} {'direct':>8s} {'RefSeq':>8s} {'surveyed':>8s}")
    for r in comparison.itertuples(index=False):
        imp = getattr(r, "prevalence_refseq_imputed", np.nan)
        ns = getattr(r, "n_refseq_surveyed", np.nan)
        imp_s = "     n/a" if not np.isfinite(imp) else f"{100 * imp:7.1f}%"
        ns_s = "     n/a" if not np.isfinite(ns) else f"{int(ns):8d}"
        add(f"  {r.genus:22s} {r.n_phages_here:6d} "
            f"{100 * r.prevalence_direct:7.1f}% {imp_s} {ns_s}")
    both = comparison.dropna(subset=["prevalence_refseq_imputed"])
    if len(both) > 3:
        rho = both[["prevalence_direct", "prevalence_refseq_imputed"]].corr(
            method="spearman").iloc[0, 1]
        mad = float((both["prevalence_direct"]
                     - both["prevalence_refseq_imputed"]).abs().mean())
        add("")
        add(f"  Spearman correlation between per-genus prevalence measured here "
            f"and imputed from RefSeq relatives: {rho:+.3f} over {len(both)} genera.")
        add(f"  Mean absolute difference in prevalence: {100 * mad:.1f} percentage points.")
        add("  Agreement at the genus level is the expected check on the "
            "annotation. Disagreement where it occurs is informative in its own "
            "right: a genus-level prior cannot capture variation between phages "
            "of the same genus, and that variation is exactly what the direct "
            "annotation adds.")

    within = []
    for g, members in genus.groupby(genus):
        names = [p for p in members.index if p in wide.index]
        if len(names) >= 3:
            v = wide.loc[names, "any_antidefense"]
            if v.nunique() > 1:
                within.append(g)
    add("")
    add(f"  Genera of three or more phages that are not uniform in carrying any "
        f"anti-defense system: {len(within)} "
        f"({', '.join(within) if within else 'none'})")
    add("  Those genera are where the direct annotation can separate phages that "
        "the genus-level prior necessarily treated as identical.")

    text = "\n".join(out)
    (paths.RESULTS / "phage_antidefense_report.txt").write_text(
        text + "\n", encoding="utf-8")
    print("\n" + text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
