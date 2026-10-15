"""Resolve the Picard strain names to public assemblies and download them.

The bacterial genomes are the one input this analysis has never had. They are
not in the primary repository, which ships only BLAST databases for the 33
isolation hosts, and the companion Zenodo record holds plaque images rather than
sequence. No project accession for the collection could be located. The strains
themselves, however, are drawn from published collections (NILS, IAI, ECOR and
the H series account for most of them), and the great majority are individually
deposited.

Resolving names to assemblies is the risky step, and the risk is silent: a
wrongly matched genome produces a plausible-looking result with no error. Three
things guard against it here. Candidates must be Escherichia and must carry the
strain name in the assembly's own strain field, not merely somewhere in the
record. Where NCBI embeds a Sequence Type or serotype in that field, which the
NILS series does, it is captured for comparison against the values recorded in
the collection metadata. And script 19 checks every downloaded genome against
the released core-genome distance matrix, which is the decisive test: a
misidentified genome cannot reproduce the phylogenetic relationships the
collection records for that strain.

Nothing downloaded here is redistributed. Accessions and checksums are recorded
so the set can be rebuilt.

Outputs:
    data/raw/host_genomes/<strain>.fna.gz
    data/processed/host_genome_provenance.csv
    results/host_genome_retrieval_report.txt
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

from papr import paths

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
GENOME_DIR = paths.RAW / "host_genomes"
PROVENANCE = paths.PROCESSED / "host_genome_provenance.csv"

# NCBI asks for no more than three requests a second without an API key.
REQUEST_PAUSE = 0.38

STATUS_RANK = {"complete genome": 4, "chromosome": 3, "scaffold": 2, "contig": 1}


def fetch(url: str, binary: bool = False, tries: int = 4):
    req = urllib.request.Request(url, headers={"User-Agent": "papr-host-fetch/1.0"})
    last = None
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=90) as resp:
                return resp.read() if binary else resp.read().decode()
        except Exception as exc:
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"failed after {tries} attempts: {url} ({last})")


def normalise(value: str) -> str:
    """Reduce a strain label to a comparable token.

    NCBI strain fields for this collection often carry annotations after the
    name, as in 'NILS13:A /STw 1137 /STp 837 / O16:H48'. Only the leading token
    identifies the strain.
    """
    head = re.split(r"[:/,]", str(value))[0]
    return re.sub(r"[^a-z0-9]", "", head.lower())


def extract_annotations(value: str) -> dict:
    """Pull any Sequence Type and serotype NCBI embedded in the strain field."""
    out: dict[str, str] = {}
    m = re.search(r"STw\s*(\d+)", str(value), re.I)
    if m:
        out["ncbi_st_warwick"] = m.group(1)
    m = re.search(r"\b(O\w+):(H\w+)\b", str(value))
    if m:
        out["ncbi_o_type"] = m.group(1)
        out["ncbi_h_type"] = m.group(2)
    return out


def search_assemblies(name: str) -> list[str]:
    term = f'"{name}"[All Fields] AND Escherichia[Organism]'
    url = EUTILS + "esearch.fcgi?db=assembly&retmode=json&retmax=20&term=" + \
        urllib.parse.quote(term)
    js = json.loads(fetch(url))
    return js.get("esearchresult", {}).get("idlist", [])


def summarise(uids: list[str]) -> list[dict]:
    url = EUTILS + "esummary.fcgi?db=assembly&retmode=json&id=" + ",".join(uids)
    js = json.loads(fetch(url))
    res = js.get("result", {})
    out = []
    for uid in res.get("uids", []):
        d = res[uid]
        strain = ""
        for item in (d.get("biosource", {}).get("infraspecieslist") or []):
            if item.get("sub_type") == "strain":
                strain = item.get("sub_value", "")
        out.append(dict(
            uid=uid,
            accession=d.get("assemblyaccession", ""),
            strain_field=strain,
            organism=d.get("speciesname", ""),
            status=str(d.get("assemblystatus", "")).lower(),
            contig_n50=int(d.get("contign50") or 0),
            ftp_refseq=d.get("ftppath_refseq", ""),
            ftp_genbank=d.get("ftppath_genbank", ""),
        ))
    return out


def choose(name: str, candidates: list[dict]) -> tuple[dict | None, str]:
    """Pick the best candidate, requiring the strain field to match the name."""
    target = normalise(name)
    scored = []
    for c in candidates:
        if "escherichia" not in c["organism"].lower():
            continue
        if not (c["ftp_refseq"] or c["ftp_genbank"]):
            continue
        field = normalise(c["strain_field"])
        if field == target:
            match = "exact_strain_field"
            base = 100
        elif target and field and (target in field or field in target):
            match = "partial_strain_field"
            base = 40
        else:
            continue
        score = (base
                 + 5 * STATUS_RANK.get(c["status"], 0)
                 + (3 if c["accession"].startswith("GCF") else 0)
                 + min(c["contig_n50"], 5_000_000) / 5_000_000)
        scored.append((score, match, c))
    if not scored:
        return None, "no_strain_field_match"
    scored.sort(key=lambda t: -t[0])
    return scored[0][2], scored[0][1]


def download_genome(record: dict, dest: Path) -> int:
    ftp = (record["ftp_refseq"] or record["ftp_genbank"]).replace("ftp://", "https://")
    base = ftp.rstrip("/").split("/")[-1]
    url = f"{ftp}/{base}_genomic.fna.gz"
    blob = fetch(url, binary=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    tmp.write_bytes(blob)
    tmp.replace(dest)
    return len(blob)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--limit", type=int, default=None,
                    help="resolve at most this many strains, for a trial run")
    args = ap.parse_args()

    paths.bootstrap()
    GENOME_DIR.mkdir(parents=True, exist_ok=True)
    pairs = pd.read_csv(paths.PAIRS, compression="gzip",
                        usecols=["strain"], low_memory=False)
    strain_feats = pd.read_csv(paths.STRAIN_FEATURES, low_memory=False)
    strains = sorted(pairs["strain"].astype(str).unique())
    if args.limit:
        strains = strains[:args.limit]

    previous: dict[str, dict] = {}
    if PROVENANCE.exists():
        previous = pd.read_csv(PROVENANCE).set_index("strain").to_dict("index")

    meta = strain_feats.set_index("strain")
    rows = []
    for i, name in enumerate(strains, 1):
        dest = GENOME_DIR / f"{name}.fna.gz"
        if dest.exists() and dest.stat().st_size > 0 and name in previous:
            rec = dict(previous[name])
            rec["strain"] = name
            rows.append(rec)
            continue

        row = dict(strain=name, accession="", match_type="", strain_field="",
                   status="", contig_n50=0, bytes_gz=0, sha256="",
                   ncbi_st_warwick="", ncbi_o_type="", ncbi_h_type="",
                   recorded_st=str(meta.loc[name, "sequence_type"])
                   if name in meta.index else "",
                   recorded_o_type=str(meta.loc[name, "o_antigen"])
                   if name in meta.index else "")
        try:
            uids = search_assemblies(name)
            time.sleep(REQUEST_PAUSE)
            if uids:
                cands = summarise(uids)
                time.sleep(REQUEST_PAUSE)
                best, match = choose(name, cands)
                row["match_type"] = match
                if best is not None:
                    row.update(accession=best["accession"],
                               strain_field=best["strain_field"],
                               status=best["status"],
                               contig_n50=best["contig_n50"])
                    row.update(extract_annotations(best["strain_field"]))
                    n = download_genome(best, dest)
                    row["bytes_gz"] = n
                    row["sha256"] = hashlib.sha256(dest.read_bytes()).hexdigest()
                    time.sleep(REQUEST_PAUSE)
            else:
                row["match_type"] = "no_search_hit"
        except Exception as exc:
            row["match_type"] = f"error: {type(exc).__name__}"
        rows.append(row)
        if i % 20 == 0 or i == len(strains):
            got = sum(1 for r in rows if r.get("accession"))
            print(f"  {i}/{len(strains)} resolved, {got} downloaded")
            pd.DataFrame(rows).to_csv(PROVENANCE, index=False)

    prov = pd.DataFrame(rows)
    prov.to_csv(PROVENANCE, index=False)

    got = prov[prov["accession"] != ""]
    total_mb = got["bytes_gz"].sum() / 1024 / 1024

    # Cross-check the annotations NCBI embedded against the collection metadata.
    st_check = got[(got["ncbi_st_warwick"] != "") & (got["recorded_st"] != "")].copy()
    st_check["st_agrees"] = (
        st_check["ncbi_st_warwick"].astype(str).str.strip()
        == st_check["recorded_st"].astype(str).str.replace(r"\.0$", "", regex=True).str.strip())
    o_check = got[(got["ncbi_o_type"] != "") & (got["recorded_o_type"] != "")].copy()
    o_check["o_agrees"] = (o_check["ncbi_o_type"].astype(str).str.strip()
                           == o_check["recorded_o_type"].astype(str).str.strip())

    out: list[str] = []
    add = out.append
    add("Host genome retrieval report")
    add("Generated by scripts/18_fetch_host_genomes.py")
    add("")
    add("RESOLUTION")
    add(f"  strains in the interaction matrix : {len(strains)}")
    add(f"  resolved and downloaded           : {len(got)} "
        f"({100 * len(got) / len(strains):.1f} percent)")
    add(f"  downloaded size                   : {total_mb:.0f} MiB")
    add("")
    add("  outcome by match type:")
    for kind, n in prov["match_type"].value_counts().items():
        add(f"    {str(kind):26s} {n}")
    add("")
    add("  assembly status of the genomes obtained:")
    for status, n in got["status"].value_counts().items():
        add(f"    {str(status):26s} {n}")

    add("")
    add("INDEPENDENT CHECKS ON IDENTITY")
    if len(st_check):
        add(f"  strains where NCBI embeds a Warwick Sequence Type : {len(st_check)}")
        add(f"    agreeing with the collection metadata           : "
            f"{int(st_check['st_agrees'].sum())} "
            f"({100 * st_check['st_agrees'].mean():.1f} percent)")
        bad = st_check[~st_check["st_agrees"]]
        if len(bad):
            add("    disagreements:")
            for r in bad.head(12).itertuples(index=False):
                add(f"      {r.strain:18s} NCBI ST {r.ncbi_st_warwick:>6s} "
                    f"vs recorded {r.recorded_st}")
    else:
        add("  no Sequence Types embedded in the NCBI strain fields")
    if len(o_check):
        add(f"  strains where NCBI embeds an O-antigen            : {len(o_check)}")
        add(f"    agreeing with the collection metadata           : "
            f"{int(o_check['o_agrees'].sum())} "
            f"({100 * o_check['o_agrees'].mean():.1f} percent)")
    add("")
    add("  These checks cover only the strains whose NCBI record happens to "
        "carry the annotation, and they confirm identity rather than "
        "establishing it. The decisive test is in script 19, which asks "
        "whether each downloaded genome occupies the position in the "
        "core-genome distance matrix that the collection records for that "
        "strain.")

    missing = prov[prov["accession"] == ""]
    add("")
    add("WHAT COULD NOT BE RESOLVED")
    add(f"  strains without a usable assembly : {len(missing)}")
    if len(missing):
        add(f"  {', '.join(missing['strain'].head(40))}"
            + (" ..." if len(missing) > 40 else ""))
    add("  Any analysis built on these genomes covers a subset of the "
        "collection and must say which strains it covers, because the strains "
        "that are deposited are unlikely to be a random sample of those that "
        "are not.")

    text = "\n".join(out)
    (paths.RESULTS / "host_genome_retrieval_report.txt").write_text(
        text + "\n", encoding="utf-8")
    print("\n" + text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
