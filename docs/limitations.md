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

## Resident prophages: tested, and untestable here

An earlier version listed this as the largest open gap, on the grounds that the
bacterial genomes were unavailable. They were retrieved (script 18), validated
against the released core-genome distance matrix (script 19) and screened for
prophages related to each challenging phage (script 20). The gap is closed in
the sense that the measurement was made; it is not closed in the sense of
producing an answer.

Across 24,382 analysable pairs covering 254 strains, only 8 have a host
carrying a prophage closely related to the challenging phage, and they involve
just 2 of the 96 phages. Those two plaque on 0.98 percent of strains, so 0.08
lysis events are expected among the carriers under the null and 0 are observed.
That is uninformative rather than negative, and nearly doubling genome coverage
added no further pairs.

The obstacle is the panel. Superinfection exclusion requires a temperate phage
that nonetheless forms plaques, and the Guelin collection was assembled for
lytic activity: 2 of 96 phages carry an integrase and both barely plaque.
Testing the mechanism needs a collection built for that question.

Two narrower limits remain. Exact k-mer containment detects prophages closely
related to a collection phage and misses diverged ones, which is the right
behaviour for repressor-mediated immunity but means this is not a measure of
prophage content in general. And the recombinase counts reported alongside are
an index of mobile element content, not a prophage count: the same profiles
match XerC, XerD and integron recombinases, and the median of 28 per genome is
far above the number of prophages an Escherichia genome actually carries.

## Receptor binding proteins: compared, and it did not help

Also superseded. The released table names each phage's receptor binding protein
without a sequence, but those identifiers match the proteins predicted in
script 07 for 125 of 130 entries, and the matched proteins have a median length
of 677 residues against 127 for all predicted proteins. Script 12 recovered
them and rebuilt the receptor gate from receptor binding protein homology.

It is no better than the genus-based gate: 1.26 fold lysis enrichment against
1.34. Since the binding protein is the right biological object, the weakness
lies in the anchor rather than the grouping, and no reanalysis fixes that.
## The shrunken intervals are narrower than they should be

The joint estimation in script 13 is empirical Bayes: the pooled mean and the
between-system variance are estimated from the same data and then treated as
known when each system's posterior interval is formed. Those intervals are
therefore anticonservative, most of all for the rare systems that depend most on
the pooled quantities. A fully Bayesian treatment would widen them and would
probably remove the weakest of the nine systems whose intervals currently
exclude an odds ratio of 1. The standard error of the pooled mean alone is not
negligible against the median posterior standard deviation.

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

## Lysogeny is only partly ruled out

Script 17 shows that temperate phages cannot explain much of the unexplained
failure in this collection, because only two of ninety-six carry a tyrosine
integrase. That argument rests on a genomic marker, not on an assay. Carrying an
integrase shows a phage can integrate, not that it did so in any particular
well, and the absence of one does not mean a phage always lyses: poor adsorption
kinetics, abortive infection without plaque formation, or a host growth
condition all produce the same empty spot. The check bounds how much lysogeny
could contribute; it does not measure it.

## The host genomes cover two thirds of the collection, non-randomly

Scripts 18 to 20 rest on 254 genomes validated out of 287 retrieved, which is
63 percent of the 402 strains. Strains deposited in public archives are not a
random sample: reference and pathogenic isolates are deposited more often than
unremarkable commensals. Identity was validated against the released core-genome distance matrix rather
than assumed from the name match, and 33 genomes were dropped. That check
cannot catch a wrong genome whose true strain happens to sit in the same part
of the tree, and it is relative to the genomes retrieved: ED1a and BL21 both
failed on a smaller set of 162 and both pass at rank 1 on the full 287, because
their nearest relatives had not been retrieved yet. A genome is only as
checkable as its neighbourhood is complete.

## The accessory genome scan has a specific blind spot

Script 25 found no accessory locus associated with resistance, and that
statement needs its boundaries attached. The null deliberately removes
lineage-confined signal, so a locus carried by exactly one clade is invisible by
construction, and in Escherichia a great deal of accessory content is exactly
that. For the same reason the empirical null is only defined for 3,144 of the
10,351 loci: the rest sit in bins holding too few comparable loci, or in bins
where clade markers collapse the spread, and those 7,207 were left untested
rather than found null. The scan covers 254 of 402 strains, tests presence
rather than sequence so
point mutations are invisible, and uses each strain's average across 96 phages
as the outcome, which dilutes anything acting on one phage lineage. It also
samples one k-mer in twenty, so a locus shorter than about a hundred bases could
fall through.

The two loci it did flag are associated with increased lysis, rest on ten
carrier strains each, and carry no recognised domain. They are observations for
follow-up, not mechanisms.

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
