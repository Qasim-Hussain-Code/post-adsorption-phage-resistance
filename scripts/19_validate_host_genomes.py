"""Check that each downloaded genome really is the strain it is labelled with.

Matching a strain name to a public assembly is the step where this part of the
analysis could go silently wrong. A wrong genome produces results that look
entirely normal. The collection, however, ships a core-genome distance matrix
covering all 402 strains, and that matrix is a fingerprint: it records how far
each strain sits from every other one. A genome that has been matched to the
wrong strain cannot reproduce those relationships.

Each downloaded genome is sketched, all pairwise k-mer distances are computed,
and for every strain the vector of its distances to the others is compared with
the corresponding row of the released matrix. A correctly identified genome
gives a strong positive correlation. A misidentified one does not, and is
dropped rather than carried forward.

Two further checks are applied: genome length has to be plausible for
Escherichia, and the strain the released matrix says is closest has to be among
the closest by k-mer distance.

Outputs:
    data/processed/host_genome_sketches.npz
    data/processed/host_genome_validation.csv
    results/tables/table36_host_genome_validation.csv
    results/host_genome_validation_report.txt
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
from scipy import stats

from papr import paths
from papr.sketch import K, SCALED, jaccard, mash_distance, read_sequences, sketch_file

GENOME_DIR = paths.RAW / "host_genomes"
SKETCH_FILE = paths.PROCESSED / "host_genome_sketches.npz"
VALIDATION = paths.PROCESSED / "host_genome_validation.csv"

# A genome whose distance profile correlates with the released matrix below
# this is not the strain it claims to be.
MIN_PROFILE_CORRELATION = 0.80
# Escherichia genomes outside this range indicate a contaminated or wrong
# assembly.
MIN_GENOME_MB, MAX_GENOME_MB = 4.0, 6.5
# How far down the k-mer ranking the released nearest neighbour may fall.
NEAREST_RANK_TOLERANCE = 10


def build_sketches(strains: list[str]) -> tuple[dict[str, np.ndarray], dict[str, float]]:
    cached: dict[str, np.ndarray] = {}
    lengths: dict[str, float] = {}
    if SKETCH_FILE.exists():
        with np.load(SKETCH_FILE, allow_pickle=False) as z:
            for name in z.files:
                if name.endswith("__len"):
                    lengths[name[:-5]] = float(z[name])
                else:
                    cached[name] = z[name]
    todo = [s for s in strains if s not in cached]
    for i, s in enumerate(todo, 1):
        path = GENOME_DIR / f"{s}.fna.gz"
        if not path.exists():
            continue
        cached[s] = sketch_file(path)
        lengths[s] = sum(c.size for c in read_sequences(path)) / 1e6
        if i % 25 == 0 or i == len(todo):
            print(f"  sketched {i}/{len(todo)}")
    if todo:
        payload = dict(cached)
        payload.update({f"{k}__len": np.array(v) for k, v in lengths.items()})
        np.savez_compressed(SKETCH_FILE, **payload)
    return cached, lengths


def main() -> int:
    paths.bootstrap()
    prov = pd.read_csv(paths.PROCESSED / "host_genome_provenance.csv")
    have = prov[prov["accession"].astype(str) != ""].copy()
    strains = sorted(s for s in have["strain"].astype(str)
                     if (GENOME_DIR / f"{s}.fna.gz").exists())
    print(f"genomes on disk: {len(strains)}")
    if len(strains) < 10:
        raise SystemExit("too few genomes to validate; run script 18 first")

    sketches, lengths = build_sketches(strains)
    strains = [s for s in strains if s in sketches and sketches[s].size > 0]
    print(f"sketched: {len(strains)}")

    core = pd.read_csv(paths.CORE_DISTANCES, sep="\t", index_col=0)
    core.index = core.index.astype(str)
    core.columns = core.columns.astype(str)
    usable = [s for s in strains if s in core.index]
    print(f"with a row in the released distance matrix: {len(usable)}")

    n = len(usable)
    kd = np.zeros((n, n))
    for i in range(n):
        si = sketches[usable[i]]
        for j in range(i + 1, n):
            d = mash_distance(si, sketches[usable[j]], K)
            kd[i, j] = kd[j, i] = d
        if (i + 1) % 25 == 0 or i == n - 1:
            print(f"  distances {i + 1}/{n}")

    cd = core.loc[usable, usable].to_numpy(dtype=float)
    cd = 0.5 * (cd + cd.T)
    np.fill_diagonal(cd, 0.0)

    iu = np.triu_indices(n, 1)
    overall_spearman = float(stats.spearmanr(kd[iu], cd[iu]).statistic)
    overall_pearson = float(np.corrcoef(kd[iu], cd[iu])[0, 1])

    rows = []
    for i, s in enumerate(usable):
        others = [j for j in range(n) if j != i]
        rho = float(stats.spearmanr(kd[i, others], cd[i, others]).statistic)
        nearest_core = others[int(np.argmin(cd[i, others]))]
        order = sorted(others, key=lambda j: kd[i, j])
        rank = order.index(nearest_core) + 1
        mb = lengths.get(s, np.nan)
        ok_len = MIN_GENOME_MB <= mb <= MAX_GENOME_MB if np.isfinite(mb) else False
        ok_rho = np.isfinite(rho) and rho >= MIN_PROFILE_CORRELATION
        ok_rank = rank <= NEAREST_RANK_TOLERANCE
        rows.append(dict(
            strain=s, genome_mb=mb, n_hashes=int(sketches[s].size),
            profile_correlation=rho,
            nearest_core_neighbour=usable[nearest_core],
            kmer_rank_of_core_nearest=rank,
            kmer_distance_to_core_nearest=float(kd[i, nearest_core]),
            length_ok=ok_len, correlation_ok=ok_rho, nearest_ok=ok_rank,
            validated=bool(ok_len and ok_rho and ok_rank)))
    val = pd.DataFrame(rows).sort_values("profile_correlation")
    val = val.merge(have[["strain", "accession", "status", "match_type"]],
                    on="strain", how="left")
    val.to_csv(VALIDATION, index=False)
    val.to_csv(paths.TABLES / "table36_host_genome_validation.csv", index=False)

    good = val[val["validated"]]

    out: list[str] = []
    add = out.append
    add("Host genome validation report")
    add("Generated by scripts/19_validate_host_genomes.py")
    add("")
    add("WHAT WAS CHECKED")
    add(f"  genomes downloaded and sketched : {len(strains)}")
    add(f"  with a row in the released core-genome matrix : {len(usable)}")
    add(f"  sketch parameters : k = {K}, scaled = {SCALED}")
    add("")
    add("AGREEMENT WITH THE RELEASED DISTANCE MATRIX")
    add(f"  Spearman over all {len(iu[0])} strain pairs : {overall_spearman:+.4f}")
    add(f"  Pearson over the same pairs               : {overall_pearson:+.4f}")
    add("  The released matrix comes from a core-genome alignment and the "
        "distances here come from whole-genome k-mer containment; they are "
        "different measurements of the same phylogeny, so a strong positive "
        "correlation is what correct identification looks like.")

    add("")
    add("PER-STRAIN OUTCOME")
    add(f"  passed all three checks : {int(val['validated'].sum())} of {len(val)} "
        f"({100 * val['validated'].mean():.1f} percent)")
    add(f"    failed on distance profile correlation : "
        f"{int((~val['correlation_ok']).sum())}")
    add(f"    failed on genome length                : {int((~val['length_ok']).sum())}")
    add(f"    failed on nearest neighbour rank       : {int((~val['nearest_ok']).sum())}")
    add(f"  median profile correlation among those that passed : "
        f"{good['profile_correlation'].median():.3f}")

    worst = val.head(12)
    add("")
    add("  Weakest profile correlations:")
    add(f"    {'strain':18s} {'rho':>7s} {'Mb':>6s} {'rank':>5s} {'verdict':>10s}")
    for r in worst.itertuples(index=False):
        add(f"    {r.strain:18s} {r.profile_correlation:7.3f} {r.genome_mb:6.2f} "
            f"{r.kmer_rank_of_core_nearest:5d} "
            f"{'kept' if r.validated else 'dropped':>10s}")

    add("")
    add("INTERPRETATION")
    if overall_spearman > 0.8:
        add("  The downloaded genomes reproduce the phylogenetic structure the "
            "collection records. Name-based matching to public assemblies has "
            "therefore worked for the great majority of strains, and the few "
            "that fail are identifiable and removed rather than left to "
            "contaminate the analysis.")
    else:
        add("  The downloaded genomes do not reproduce the recorded "
            "phylogenetic structure well. Name-based matching cannot be "
            "trusted at this level of agreement and nothing downstream should "
            "rely on these genomes.")
    add(f"  {len(good)} genomes are carried forward. That is "
        f"{100 * len(good) / 402:.0f} percent of the 402 strains in the "
        "interaction matrix, so any analysis built on them describes a subset "
        "of the collection and is reported as such.")
    add("  Strains deposited in public archives are not a random sample of the "
        "collection: reference and pathogenic isolates are deposited more "
        "often than unremarkable commensals. Associations estimated on this "
        "subset inherit that selection.")

    text = "\n".join(out)
    (paths.RESULTS / "host_genome_validation_report.txt").write_text(
        text + "\n", encoding="utf-8")
    print("\n" + text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
