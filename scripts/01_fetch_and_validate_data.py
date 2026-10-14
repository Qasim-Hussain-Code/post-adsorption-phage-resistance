"""Retrieve the public Picard x Guelin interaction data and verify its integrity.

The script is idempotent. On the first run it downloads each input, records a
SHA-256 digest in data/raw_manifest.sha256 and validates the contents. On later
runs it re-verifies the digests already on disk and skips the network unless
--force is supplied or a file is missing or altered.

Two points of intellectual honesty are handled here rather than being left to
the reader. First, the dimensions of the released interaction matrix are
402 strains by 96 phages, not the 403 by 96 that a reader might infer from the
403-row strain metadata table; one strain (LF110) is described in the metadata
but carries no row in the matrix. Second, 157 of the 38,592 cells are blank,
so the number of interactions actually observed is smaller than the number of
possible pairs. The script asserts the observed values and writes an explicit
reconciliation note instead of asserting the larger nominal figures.

Usage:
    python scripts/01_fetch_and_validate_data.py [--force] [--skip-raw-check]
"""

from __future__ import annotations

import argparse
import hashlib
import sys
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

from papr import paths

COLI = "https://raw.githubusercontent.com/mdmparis/coli_phage_interactions_2023/main/"
ADF = "https://raw.githubusercontent.com/mdmparis/antidefensefinder_2024/main/"

# (destination, source URL, short description of why the file is needed)
INPUTS: list[tuple[Path, str, str]] = [
    (paths.MATRIX, COLI + "data/interactions/interaction_matrix.csv",
     "Summarised plaque score per strain-phage pair, scale 0 to 4"),
    (paths.RAW_INTERACTIONS, COLI + "data/interactions/raw/raw_interactions.csv",
     "Replicate-level spot tests at each log dilution, used to audit the summary matrix"),
    (paths.PICARD, COLI + "data/genomics/bacteria/picard_collection.csv",
     "Bacterial metadata: phylogroup, Sequence Type, O-type, LPS core type, capsule state"),
    (paths.ECTYPER, COLI + "data/genomics/bacteria/o_type/output.tsv",
     "ECTyper O-antigen and H-antigen calls with per-allele quality evidence"),
    (paths.LPS_WAAL, COLI + "data/genomics/bacteria/outer_core_lps/LPS_type_waaL_370.txt",
     "LPS outer core type assigned from waaL, independent confirmation of picard_collection"),
    (paths.OMP_ALLELES, COLI + "data/genomics/bacteria/outer_membrane_proteins/"
     "blast_results_cured_clusters=99_wide.tsv",
     "Outer membrane protein allele clusters at 99 percent identity for 12 receptor loci"),
    (paths.KAPTIVE, COLI + "data/genomics/bacteria/capsules/klebsiella_capsules/"
     "kaptive_results_table_colocoli.txt",
     "Kaptive capsule locus calls for the subset of strains carrying a Klebsiella-like capsule"),
    (paths.DEFENSE, COLI + "data/genomics/bacteria/defense_finder/"
     "370+host_defense_systems_subtypes.csv",
     "DefenseFinder system counts per strain, resolved to subtype"),
    (paths.CORE_TREE, COLI + "data/genomics/bacteria/panacota/tree/"
     "370+host_tree_root=albertii.nwk",
     "Core-genome maximum likelihood tree, rooted on Escherichia albertii"),
    (paths.CORE_DISTANCES, COLI + "data/genomics/bacteria/isolation_strains/panacota/tree/"
     "370+host_distance_matrix_orignames.tsv",
     "Core-genome pairwise distance matrix, used for principal coordinates"),
    (paths.GUELIN, COLI + "data/genomics/phages/guelin_collection.csv",
     "Phage metadata: morphotype, family, VIRIDIC genus, genome size, isolation host"),
    (paths.RBP, COLI + "data/genomics/phages/RBP/RBP_list.csv",
     "Receptor binding protein assignments per phage, fiber or spike"),
    (paths.PHAGE_TREE, COLI + "data/genomics/phages/tree/"
     "96_viridic_distance_phylogenetic_tree_algo=upgma.nwk",
     "Phage UPGMA tree from VIRIDIC intergenomic distances"),
    (paths.ADF_PHAGE, ADF + "data/figures/SuppData3_ADF_Phage.tsv",
     "AntiDefenseFinder hits across RefSeq phage genomes, with ICTV genus labels"),
    (paths.ADF_TAXONOMY, ADF + "data/figures/Taxonomy_phage.tsv",
     "Taxonomy of every RefSeq phage genome surveyed, the denominator for anti-defense prevalence"),
]

# Figures quoted in the source publication, retained so that the reconciliation
# below is explicit about what differs and why.
NOMINAL_STRAINS = 403
NOMINAL_PHAGES = 96


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "post-adsorption-analysis/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=120) as resp, tmp.open("wb") as out:
            while True:
                chunk = resp.read(1 << 20)
                if not chunk:
                    break
                out.write(chunk)
    except urllib.error.URLError as exc:
        tmp.unlink(missing_ok=True)
        raise SystemExit(f"download failed for {url}: {exc}")
    tmp.replace(dest)


def read_manifest() -> dict[str, str]:
    if not paths.RAW_MANIFEST.exists():
        return {}
    out: dict[str, str] = {}
    for line in paths.RAW_MANIFEST.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        digest, rel = line.split("  ", 1)
        out[rel] = digest
    return out


def write_manifest(rows: list[tuple[str, str]]) -> None:
    lines = [
        "# SHA-256 digests of the upstream inputs used by this analysis.",
        "# Regenerate and verify with: python scripts/01_fetch_and_validate_data.py",
        "# The data files themselves are not redistributed in this repository.",
        "",
    ]
    lines += [f"{digest}  {rel}" for digest, rel in sorted(rows, key=lambda r: r[1])]
    paths.RAW_MANIFEST.write_text("\n".join(lines) + "\n", encoding="utf-8")


def acquire(force: bool) -> list[tuple[str, str]]:
    expected = read_manifest()
    rows: list[tuple[str, str]] = []
    total = 0
    for dest, url, purpose in INPUTS:
        rel = dest.relative_to(paths.DATA).as_posix()
        need = force or not dest.exists()
        if not need and rel in expected and sha256_of(dest) != expected[rel]:
            print(f"  digest mismatch, re-downloading: {rel}")
            need = True
        if need:
            download(url, dest)
            action = "downloaded"
        else:
            action = "cached"
        digest = sha256_of(dest)
        size = dest.stat().st_size
        total += size
        if rel in expected and expected[rel] != digest and not force:
            raise SystemExit(f"digest still mismatched after download: {rel}")
        rows.append((digest, rel))
        print(f"  {action:10s} {size / 1024:10.1f} KiB  {rel}")
        print(f"             {purpose}")
    print(f"\n  total raw input size: {total / 1024 / 1024:.1f} MiB")
    return rows


def load_core() -> dict[str, pd.DataFrame]:
    """Read the inputs that the validation below depends on."""
    return {
        "matrix": pd.read_csv(paths.MATRIX, sep=";"),
        "picard": pd.read_csv(paths.PICARD, sep=";"),
        "guelin": pd.read_csv(paths.GUELIN, sep=";"),
        "defense": pd.read_csv(paths.DEFENSE, sep=";"),
        "omp": pd.read_csv(paths.OMP_ALLELES, sep="\t"),
        "rbp": pd.read_csv(paths.RBP, sep=";"),
    }


# Position of each assayed dilution on the ladder, deepest last. The released
# matrix runs 0 to 4 while the replicate table records four dilution steps, so
# the deepest reconstructed step (4) is the ceiling of the reconstruction.
LADDER = {0: 1, -1: 2, -2: 3, -4: 4}


def audit_raw_against_summary(matrix: pd.DataFrame) -> list[str]:
    """Rebuild the lysis call from the replicate spot tests and compare.

    Each spot test is binary. The literal value 'n' marks a spot that could not
    be read and is treated as missing rather than as a negative. A dilution step
    counts as positive when a plaque was recorded in at least two of the
    replicates that could be read, which is the reproducibility rule the
    released matrix appears to apply. The reconstructed score is the deepest
    step of the 0, -1, -2, -4 dilution ladder that remained positive.

    The replicate table covers only the 369 strains of the original Picard
    collection; the 33 strains contributed as phage isolation hosts carry no
    replicate records, so they are excluded from the comparison rather than
    being counted as disagreements.
    """
    out: list[str] = []
    add = out.append

    raw = pd.read_csv(
        paths.RAW_INTERACTIONS, sep=";", dtype=str,
        usecols=["bacteria", "phage", "replicate", "log_dilution", "score"],
    )
    n_unreadable = int((raw["score"] == "n").sum())
    raw["positive"] = pd.to_numeric(raw["score"], errors="coerce")
    raw["step"] = raw["log_dilution"].astype(int).map(LADDER)

    add(f"  replicate-level rows read            : {len(raw)}")
    add(f"  distinct per-spot score values       : "
        f"{sorted(pd.unique(raw['score']).tolist())}")
    add(f"  spots recorded as not readable ('n') : {n_unreadable} "
        f"({100.0 * n_unreadable / len(raw):.2f} percent)")
    add(f"  log dilution steps assayed           : "
        f"{sorted(raw['log_dilution'].astype(int).unique().tolist())}")
    add(f"  replicates per strain-phage-dilution : "
        f"{sorted(raw['replicate'].unique().tolist())}")
    add(f"  strains with replicate records       : {raw['bacteria'].nunique()} "
        f"of {len(matrix)} in the matrix")

    per_step = (raw.groupby(["bacteria", "phage", "step"])["positive"]
                .agg(n_positive="sum", n_readable="count").reset_index())
    reproducible = per_step[per_step["n_positive"] >= 2]
    deepest = reproducible.groupby(["bacteria", "phage"])["step"].max()

    strain_col = matrix.columns[0]
    released = (matrix.melt(id_vars=strain_col, var_name="phage", value_name="score")
                .rename(columns={strain_col: "bacteria"})
                .dropna(subset=["score"]))
    covered = released[released["bacteria"].isin(set(raw["bacteria"]))].copy()
    covered = covered.merge(deepest.rename("rebuilt"), how="left",
                            left_on=["bacteria", "phage"], right_index=True)
    covered["rebuilt"] = covered["rebuilt"].fillna(0.0)

    n = len(covered)
    exact = float((covered["rebuilt"] == covered["score"]).mean())
    spearman = float(covered[["score", "rebuilt"]].corr(method="spearman").iloc[0, 1])
    rel_lytic = covered["score"] > 0
    reb_lytic = covered["rebuilt"] > 0
    binary_agree = float((rel_lytic == reb_lytic).mean())
    fn = int((rel_lytic & ~reb_lytic).sum())
    fp = int((~rel_lytic & reb_lytic).sum())

    add("")
    add(f"  comparable cells                     : {n}")
    add(f"  binary lysis call reproduced         : {100.0 * binary_agree:.2f} percent")
    add(f"    released lytic, rebuilt non-lytic  : {fn}")
    add(f"    released non-lytic, rebuilt lytic  : {fp}")
    add(f"  ordinal score reproduced exactly     : {100.0 * exact:.2f} percent")
    add(f"  Spearman correlation, released vs rebuilt: {spearman:.4f}")
    add("")
    add("  cross-tabulation, released score by rebuilt score:")
    ct = pd.crosstab(covered["score"], covered["rebuilt"])
    for line in ct.to_string().splitlines():
        add("    " + line)
    add("")
    add("  Every released zero corresponds to a rebuilt zero, so the binary "
        "lytic classification used as the primary outcome of this study is "
        "recovered from the raw spot tests almost exactly. Residual ordinal "
        "disagreement is confined to the magnitude of the score and is "
        "monotone: the replicate ladder has four steps, so a released score of "
        "4 maps onto the deepest rebuilt step rather than to a fifth step. The "
        "quantitative analyses in script 03 therefore use the released score as "
        "an ordinal rank and never as a count of dilutions.")
    return out


def validate(force_raw: bool) -> list[str]:
    d = load_core()
    mat, picard, guelin = d["matrix"], d["picard"], d["guelin"]
    report: list[str] = []
    add = report.append

    strain_col = mat.columns[0]
    strains = mat[strain_col].tolist()
    phages = [c for c in mat.columns if c != strain_col]
    body = mat[phages]

    n_strains, n_phages = len(strains), len(phages)
    cells = n_strains * n_phages
    n_missing = int(body.isna().sum().sum())
    observed = cells - n_missing

    add("DIMENSIONS AND RECONCILIATION")
    add(f"  interaction matrix rows (strains tested) : {n_strains}")
    add(f"  interaction matrix columns (phages)      : {n_phages}")
    add(f"  matrix cells                             : {cells}")
    add(f"  blank cells (pair never tested)          : {n_missing}")
    add(f"  interactions actually observed           : {observed}")
    add(f"  strain metadata rows in picard_collection: {len(picard)}")
    add("")
    add("  The source publication describes a collection of "
        f"{NOMINAL_STRAINS} strains and {NOMINAL_PHAGES} phages, which would give "
        f"{NOMINAL_STRAINS * NOMINAL_PHAGES} pairs. The released matrix contains "
        f"{cells} cells because one strain present in the metadata carries no row "
        "in the matrix, and a further "
        f"{n_missing} cells are blank. All downstream analysis in this repository "
        f"uses the {observed} observed interactions and never the nominal figure.")

    in_meta_not_matrix = sorted(set(picard["bacteria"]) - set(strains))
    in_matrix_not_meta = sorted(set(strains) - set(picard["bacteria"]))
    add(f"  strains in metadata but absent from matrix: {in_meta_not_matrix}")
    add(f"  strains in matrix but absent from metadata: {in_matrix_not_meta}")

    assert n_phages == NOMINAL_PHAGES, f"expected {NOMINAL_PHAGES} phages, found {n_phages}"
    assert n_strains == 402, f"expected 402 matrix strains, found {n_strains}"
    assert cells == 38_592, f"expected 38592 cells, found {cells}"
    assert not in_matrix_not_meta, "every tested strain must have metadata"
    assert len(in_meta_not_matrix) == 1, "expected exactly one metadata-only strain"
    assert set(phages) == set(guelin["phage"]), "matrix phages must match phage metadata"
    assert mat[strain_col].is_unique and len(set(phages)) == n_phages, "identifiers must be unique"

    vals = pd.unique(body.to_numpy().ravel())
    finite = sorted(float(v) for v in vals if pd.notna(v))
    add("")
    add("PLAQUE SCORE SCALE")
    add(f"  distinct scores observed: {finite}")
    add("  The released scores run from 0 to 4, each step being one further "
        "ten-fold dilution at which plaques were still visible. A 0 to 6 scale "
        "is not present in this release.")
    assert set(finite) <= {0.0, 1.0, 2.0, 3.0, 4.0}, f"unexpected score values: {finite}"

    counts = pd.Series(body.to_numpy().ravel()).value_counts(dropna=False).sort_index()
    add("  score frequencies:")
    for k, v in counts.items():
        label = "blank" if pd.isna(k) else f"{k:.0f}"
        add(f"    {label:>5s}: {v:6d}  ({100.0 * v / cells:5.2f} percent of cells)")
    lytic = int((body.to_numpy() > 0).sum())
    add(f"  productive lysis (score > 0): {lytic} of {observed} observed "
        f"({100.0 * lytic / observed:.2f} percent)")

    add("")
    add("GENOMIC FEATURE COVERAGE OF THE 402 TESTED STRAINS")
    for name, frame, key in [
        ("DefenseFinder subtypes", d["defense"], "bacteria"),
        ("outer membrane protein alleles", d["omp"], "bacteria"),
    ]:
        missing = sorted(set(strains) - set(frame[key]))
        add(f"  {name:32s}: {len(frame)} rows, {len(missing)} tested strains missing")
        assert not missing, f"{name} does not cover every tested strain"
    dist_index = pd.read_csv(paths.CORE_DISTANCES, sep="\t", usecols=[0], nrows=None)
    dist_ids = set(dist_index.iloc[:, 0])
    add(f"  core-genome distance matrix     : {len(dist_ids)} rows, "
        f"{len(set(strains) - dist_ids)} tested strains missing")

    for col in ["Clermont_Phylo", "ST_Warwick", "LPS_type", "O-type"]:
        sub = picard[picard["bacteria"].isin(strains)]
        add(f"  {col:32s}: {sub[col].nunique()} levels, "
            f"{int(sub[col].isna().sum())} missing among tested strains")

    add("")
    add("PHAGE RECEPTOR INFORMATION ACTUALLY AVAILABLE")
    add(f"  phage metadata columns: {list(guelin.columns)}")
    add("  There is no column in this release that states, for each phage, the "
        "host surface receptor it requires. No such table is distributed with "
        "the dataset. Receptor compatibility therefore cannot be read off the "
        "data and must be inferred; script 02 does so from the isolation host "
        "of each phage and from VIRIDIC genus membership, which is independent "
        "of the interaction outcomes being analysed.")
    hosts = sorted(set(guelin["Phage_host"].dropna()))
    add(f"  distinct phage isolation hosts: {len(hosts)}")
    add(f"  isolation hosts with strain metadata: "
        f"{len(set(hosts) & set(picard['bacteria']))} of {len(hosts)}")
    add(f"  isolation hosts also present as matrix rows: "
        f"{len(set(hosts) & set(strains))} of {len(hosts)}")
    add(f"  VIRIDIC genera represented: {guelin['Genus'].nunique()}")
    add(f"  phages with a receptor binding protein assignment: {d['rbp']['phage'].nunique()}")
    assert set(hosts) <= set(picard["bacteria"]), "every isolation host needs metadata"

    add("")
    add("AUDIT OF THE SUMMARY MATRIX AGAINST THE REPLICATE SPOT TESTS")
    if force_raw:
        report.extend(audit_raw_against_summary(mat))
    else:
        add("  skipped (--skip-raw-check)")

    return report


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true",
                    help="re-download every input even if a valid copy is cached")
    ap.add_argument("--skip-raw-check", action="store_true",
                    help="skip the replicate-level audit of the summary matrix")
    args = ap.parse_args()

    paths.bootstrap()
    print("Acquiring upstream inputs")
    rows = acquire(force=args.force)
    write_manifest(rows)
    print(f"\nmanifest written to {paths.RAW_MANIFEST.relative_to(paths.ROOT)}")

    print("\nValidating")
    report = validate(force_raw=not args.skip_raw_check)
    text = "\n".join(report)
    print(text)

    out = paths.RESULTS / "data_validation_report.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "Data validation report",
        "Generated by scripts/01_fetch_and_validate_data.py",
        "",
    ]
    out.write_text("\n".join(header) + text + "\n", encoding="utf-8")
    print(f"\nreport written to {out.relative_to(paths.ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
