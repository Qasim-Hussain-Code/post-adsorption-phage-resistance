# Data provenance

## Sources

All inputs are publicly released data from two repositories of the Microbial
Evolutionary Genomics unit, Institut Pasteur.

**Primary.** `mdmparis/coli_phage_interactions_2023`, the open release
accompanying:

> Gaborieau, B. et al. Prediction of strain level phage-host interactions across
> the Escherichia genus using only genomic information. *Nature Microbiology*
> (2024).

**Anti-defense.** `mdmparis/antidefensefinder_2024`, the open release
accompanying the AntiDefenseFinder survey of anti-defense systems across RefSeq
phage genomes.

## Redistribution

No upstream data file is committed to this repository. `data/raw/` is excluded
by `.gitignore`. What is tracked is `data/raw_manifest.sha256`, which records a
SHA-256 digest, size and source URL for each of the 15 retrieved files, and
`scripts/01_fetch_and_validate_data.py`, which re-downloads them and verifies
the digests. A reader reproduces the inputs rather than receiving a copy of
them, and the upstream data remains under the license of its originating
repository.

## Files retrieved

| Local path | Purpose |
|---|---|
| `interactions/interaction_matrix.csv` | summarised plaque score per pair, 0 to 4 |
| `interactions/raw_interactions.csv` | replicate-level spot tests at each dilution |
| `bacteria/picard_collection.csv` | phylogroup, Sequence Type, O-type, LPS core, capsule |
| `bacteria/ectyper_o_type.tsv` | ECTyper O-antigen and H-antigen calls with evidence |
| `bacteria/lps_type_waaL.txt` | LPS outer core assigned from waaL |
| `bacteria/omp_allele_clusters.tsv` | outer membrane protein alleles at 12 loci |
| `bacteria/kaptive_capsule_results.txt` | Kaptive capsule locus calls |
| `bacteria/defense_systems_subtypes.csv` | DefenseFinder systems per strain, by subtype |
| `bacteria/core_genome_tree.nwk` | core-genome tree rooted on E. albertii |
| `bacteria/core_genome_distance_matrix.tsv` | pairwise core-genome distances |
| `phages/guelin_collection.csv` | morphotype, family, VIRIDIC genus, isolation host |
| `phages/rbp_list.csv` | receptor binding protein assignments |
| `phages/phage_viridic_tree.nwk` | UPGMA tree from VIRIDIC distances |
| `antidefense/adf_hits_refseq_phage.tsv` | AntiDefenseFinder hits, with ICTV genus |
| `antidefense/phage_taxonomy_refseq.tsv` | taxonomy of every surveyed genome |

Total download for the core pipeline is approximately 21 MB.

## Additional inputs for the direct phage annotation

Steps 07 to 10 fetch two further sets of inputs, handled by script 07 rather
than script 01 so that the core pipeline stays light.

| Local path | Source | Size |
|---|---|---|
| `data/raw/phage_genomes/*.fna` | the 96 phage genome FASTA files from the same primary repository | 6.9 MB |
| `data/raw/adf_models/profiles/*.hmm` | 178 anti-defense profile hidden Markov models | about 14 MB |
| `data/raw/adf_models/definitions/**/*.xml` | 156 anti-defense system definitions | small |

The profiles and definitions come from release 3.1.0 of
`mdmparis/defense-finder-models`. That release is distributed as a single 49 MB
archive containing mostly bacterial defense profiles the analysis does not need.
Script 07 streams through the archive twice, writing out only the anti-defense
parts, and deletes the archive immediately afterwards, so the project grows by
about 14 MB rather than 49 MB. The archive is streamed rather than accessed
member by member because random access into a gzip stream forces a fresh
decompression each time.

With these added, the repository and all downloaded data occupy about 51 MB,
against the 3.5 GB budget this analysis was held to.

## Discrepancies between the release and the nominal description

These were found during validation and are recorded here because downstream
numbers depend on them. `scripts/01_fetch_and_validate_data.py` asserts the
observed values and writes the reconciliation to
`results/data_validation_report.txt`.

1. **402 strains, not 403.** The strain metadata table has 403 rows. The
   interaction matrix has 402. The strain `LF110` is described in the metadata,
   is the isolation host of four phages, and carries no row in the matrix. It is
   excluded from all pair-level analysis and retained where it serves as an
   isolation host anchor.

2. **38,592 cells, of which 157 are blank.** The product 403 by 96 would give
   38,688. The matrix as released contains 38,592 cells and 157 are empty, so
   **38,435** interactions were observed. Every denominator in this repository is
   38,435 or a subset of it.

3. **Plaque score runs 0 to 4.** Each step is one further ten-fold dilution at
   which plaques remained visible. There is no 0 to 6 scale in this release.

4. **Replicate coverage is partial.** The replicate-level table covers the 369
   strains of the original Picard collection. The 33 strains contributed as
   phage isolation hosts have no replicate records, so the audit of the summary
   matrix is restricted to the 35,268 comparable cells.

5. **Capsule calls are missing for the isolation-host strains.** The four capsule
   columns are empty for the same 33 strains, which is why capsule state is used
   only in sensitivity analysis and never in the primary receptor gate.

6. **Several files are semicolon-delimited** despite a `.csv` extension
   (`interaction_matrix.csv`, `picard_collection.csv`, `guelin_collection.csv`,
   `RBP_list.csv`, the DefenseFinder subtype table). Reading them with a comma
   separator yields a single-column frame that silently parses.

7. **The per-spot score column is not numeric.** It takes values `0`, `1` and
   `n`, the last marking a reading that could not be made. It is treated as
   missing rather than as a negative.

## Bacterial genomes, retrieved individually

The 370 Picard genomes are not in the primary repository, which ships only
BLAST databases for the 33 isolation hosts. The companion Zenodo record
(10.5281/zenodo.10202713) contains a 9.9 GB archive of plaque assay images and
no sequence. No ENA or BioProject accession for the collection could be located
from the paper, the preprint or a search of ENA and NCBI.

Script 18 therefore resolves strain names against NCBI Assembly one at a time,
requiring the candidate to be Escherichia and to carry the strain name in its
own strain field. Each downloaded genome, its accession, assembly status and
SHA-256 digest are recorded in `data/processed/host_genome_provenance.csv`,
which is tracked; the genomes themselves live under `data/raw/host_genomes/`
and are not redistributed.

Coverage is partial and the reasons are worth separating. Some strains are
simply not deposited, mostly the internally numbered isolates. Others are
deposited under names shared with unrelated isolates and are rejected by the
validation in script 19. And a full pass is throttled by NCBI, which is why the
run recorded here stops short of all 402 names; the script is resumable and
re-running retries only what failed.

| Local path | Source | Size |
|---|---|---|
| `data/raw/host_genomes/<strain>.fna.gz` | NCBI Assembly, per strain | about 237 MB for 166 genomes |
| `data/raw/lifestyle_profiles/*.hmm` | Pfam via InterPro, 7 lysogeny families | small |

With these, the repository and all downloaded data occupy roughly 300 MB
against the 3.5 GB budget.

## Integrity checking

`python scripts/01_fetch_and_validate_data.py` is idempotent. On first run it
downloads and records digests. On later runs it re-verifies the cached copies
and re-downloads only what is missing or altered. Pass `--force` to refetch
everything. The validation step additionally asserts that the matrix has the
expected shape, that every tested strain has metadata, that the phage set
matches the phage metadata exactly, that score values lie in the expected set,
and that the DefenseFinder, outer membrane protein and distance matrix tables
cover all 402 tested strains.
