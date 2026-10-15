# Limitations

These are the reasons a reader should hold the conclusions in this repository
loosely, ordered roughly by how much they matter.

## A claim in an earlier version of this analysis was wrong

The first version of the variance decomposition reported that defense carriage
contributes more to the variation in lysis than host surface does. That was an
artefact of two defects: the blocks were entered in one fixed order, and the
surface block held six terms against thirty-three for defense. With a surface
block containing everything the collection types about the cell surface and an
order-independent Shapley allocation, host surface accounts for about three
times what defense does. The earlier statement is retracted in
`docs/results.md`, and the table that produced it is retained in script 03 only
to show the nested fits, with a note directing the reader to script 11.

The lesson generalises beyond this one number: any comparison between blocks of
unequal size, entered in a chosen order, can be made to favour either block.

## Receptor compatibility is inferred and the inference is weak

The central methodological problem is that this dataset contains no statement of
which surface receptor each phage requires. The gate built here uses the only
outcome-independent receptor evidence available, the strain each phage was
isolated on, pooled across a VIRIDIC genus. Measured against the data, it is a
poor discriminator: 1.34 fold enrichment in lysis, and only 35.9 percent of
failed infections attributed to the surface.

Three specific weaknesses:

- Gate stringency tracks genus size. A genus represented by one phage yields one
  permissive state and a strict gate; a genus with eighteen phages yields four
  permissive states and a permissive gate. This has nothing to do with receptor
  biology.
- One isolation host establishes that a surface state is permissive. It cannot
  establish that any other state is not, so the gate has no power to exclude
  correctly and will pass many genuinely incompatible pairs.
- Enrichment is real in some genera and absent in the two largest, so a single
  species-wide gate is not well defined.

Consequently the two-gate decomposition should be read as one defensible
estimate under a stated definition, not as a measurement of how often adsorption
fails. The primary analysis avoids depending on it by stratifying on surface
state instead.

Script 12 tested whether the weakness lies in the grouping by rebuilding the
gate from receptor binding protein homology, which is the actual determinant of
adsorption rather than a taxonomic proxy. It does not help: 1.26 fold enrichment
against 1.34 for the genus gate. The weakness is therefore in the **anchor**.
One isolation host per phage establishes that a surface state is permissive and
can never establish that another is not, so permissive sets remain small however
the phages are grouped. Fixing this would require adsorption assays rather than
plaque assays, which no reanalysis of this dataset can supply.

## Anti-defense content: this limitation has since been closed

An earlier version of this analysis imputed anti-defense content from each
phage's genus, which attenuates a real interaction toward the null and was the
main reason to doubt the negative counter-defense result. The 96 genomes have
since been annotated directly (script 07), and the rescue test repeated on the
measured calls (script 08). The answer did not change: one significant
interaction either way, direction agreeing for 11 of 12 systems, estimates
correlating at Spearman +0.993. The negative result does not rest on imputation.

What remains is narrower. Profile hits are accepted at the curated gathering
cutoffs and the models' gene-count rules are applied, but MacSyFinder's full
co-localisation logic is only exercised for the six definitions that require two
or more genes, so a handful of multi-gene systems could in principle be called
where the real genes are scattered. Anti_CRISPR, Anti_BREX, Anti_Gabija and
Anti_Hachiman are absent from all 96 genomes, so CRISPR-Cas, BREX, Gabija and
Hachiman cannot be tested for rescue in this collection at all. Their absence is
a real feature of these phage lineages, but it leaves the rescue model untested
for several of the systems most discussed in the literature.

## Resident prophages were not annotated

Superinfection exclusion by a resident prophage is a well-documented route to
post-adsorption resistance and is a plausible contributor to the 75.8 percent
residual failure rate. Identifying prophages would require running a prophage
caller over 402 bacterial genomes, which are not part of the released data used
here and which exceed the storage budget this analysis was held to. This is a
gap in the analysis, not a question the data answer negatively.

## Receptor binding protein sequences were not compared

The release names the receptor binding protein of each phage but does not give
its sequence, and the phage genome files carry nucleotides without gene
coordinates. Grouping phages by receptor binding protein similarity, which would
give a sharper receptor proxy than genus membership, was therefore not possible
without re-annotating the genomes.

## Statistical caveats

- The Mantel-Haenszel intervals reported in step 3 treat strain-phage pairs as
  independent. They are too narrow, visibly so: 25 of 33 systems reach
  significance there against 6 of 33 under the same contrast fitted with
  clustered standard errors. The clustered estimates are the inferential ones
  and the stratified estimates should be read as effect sizes, not as evidence
  strength.
- Nested-model pseudo R2 increments are descriptive. The correlated structure of
  the data inflates the accompanying likelihood ratio statistics.
- The clonal cluster threshold (0.0005, giving 224 clusters) is a judgement
  call. A sweep from 1e-5 to 1e-2 is reported so its influence can be judged,
  but the main results were not re-derived at every threshold.
- Adjusting for phylogroup and eight principal coordinates may remove genuine
  defense signal where a system is confined to one lineage. The internal check
  (all biologically impossible associations vanish, both surviving associations
  are protective) argues against gross overcorrection, but it cannot rule out
  partial overcorrection for a lineage-restricted system.
- The 20-carrier threshold excludes 104 of the 137 annotated systems. Those are
  untested, not shown to be inert.
- O-antigen serotype has 94 levels across 402 strains, so serotype-level
  stratification is sparse by construction and its estimates are correspondingly
  fragile.

## The two surviving associations are candidates, not findings

The robustness analyses in script 09 qualify the headline substantially and
should be read alongside it.

- **RM_Type_IIG loses its effect among near-identical genomes.** In the
  within-cluster design, which holds the genomic background fixed, its odds
  ratio is 0.968 (95 percent CI 0.487 to 1.923) across 7 discordant clusters
  spanning 36 strains. That is moderate power, so this is evidence against the
  association rather than mere absence of evidence.
- **Zorya_TypeII cannot be tested that way at all.** Only 2 clusters spanning
  11 strains are informative, so its interval runs from 0.060 to 1.978. Its
  point estimate is preserved (0.344 against 0.444 adjusted), but preserved
  point estimates with uninformative intervals are weak evidence.
- Zorya_TypeII fails the permutation test at both levels; RM_Type_IIG passes it
  within phylogroup but not within Sequence Type.
- No system passes all four stricter designs. Conversely, **RM_Type_III, which
  did not survive the main lineage-adjusted regression, is significant and
  protective in two of them** and may be the better candidate.

Each stricter design discards most of the data in order to remove a different
source of confounding, so disagreement between them is expected. The
conservative reading is that this dataset yields a short candidate list, not a
demonstration.

## Most of the nulls are uninformative, not negative

Across the 31 systems reported as null after clonal control, the protective odds
ratio that each test could have detected at 80 percent power ranges from 0.43 to
0.70. A system whose figure is near 0.70 can exclude only very large protective
effects. Reporting those systems as showing no association is accurate but
easily misread: for most of them the analysis simply had no power to see a
moderate effect, and they should not be cited as evidence that the system does
nothing.

## Scope

This is an observational re-analysis of a single collection of 402 Escherichia
strains and 96 phages, assayed by spot test on a four-step dilution ladder. Spot
tests measure plaque formation, not adsorption, so no result here distinguishes
failure to adsorb from failure after adsorption except through the inferred gate.

Nothing in this analysis establishes that Zorya_TypeII or RM_Type_IIG cause
resistance. It establishes that, among the systems and covariates available,
their association with reduced lysis is not explained by phage identity, host
surface type, phylogroup or core-genome position. Demonstrating causation would
require knockouts in an isogenic background.
