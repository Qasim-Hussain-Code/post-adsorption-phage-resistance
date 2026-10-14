# Limitations

These are the reasons a reader should hold the conclusions in this repository
loosely, ordered roughly by how much they matter.

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

## Anti-defense content is imputed, not observed

The 96 phages were not annotated for anti-defense genes. Each genus carries the
prevalence of each anti-defense family among its sequenced RefSeq relatives.
Imputation of this kind attenuates a real interaction toward the null, so the
largely negative counter-defense result is some way short of evidence of
absence. Annotating the 96 genomes directly is the single change most likely to
alter that conclusion.

Survey depth also varies from 2 to 207 genomes per genus. Genera with fewer than
10 were excluded, which removes four of the nineteen, including Lambdavirus,
whose lysis rate of 0.75 percent makes it uninformative in any case.

Anti_CRISPR, Anti_Gabija and Anti_Hachiman are absent from every adequately
surveyed genus here, so CRISPR-Cas, Gabija and Hachiman could not be tested for
rescue at all. Their absence from this collection's phage lineages is itself
worth noting, but it means the rescue model is untested for three of the
systems most discussed in the literature.

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
