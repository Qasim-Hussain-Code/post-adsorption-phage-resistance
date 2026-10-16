# Methods

## Data and phenotype definition

The analysis uses the open release accompanying Gaborieau et al. (Nature
Microbiology, 2024). The retrieved inputs are the summarised interaction matrix,
the replicate-level spot tests, bacterial metadata (Clermont phylogroup, Warwick
Sequence Type, ECTyper O-antigen and H-antigen calls, LPS outer core type
assigned from waaL, capsule state, Kaptive capsule locus calls), outer membrane
protein allele clusters at twelve loci, DefenseFinder system calls resolved to
subtype, the core-genome maximum likelihood tree and pairwise distance matrix,
phage metadata including VIRIDIC genus and isolation host, receptor binding
protein assignments, and the published AntiDefenseFinder survey of RefSeq phage
genomes. Each file is verified against a SHA-256 digest recorded in
`data/raw_manifest.sha256`.

The interaction matrix contains 402 strains by 96 phages. One strain present in
the metadata, `LF110`, has no row. A further 157 of the 38,592 cells are blank.
All analyses use the 38,435 observed interactions. The released plaque score
runs from 0 to 4, each step being one further ten-fold dilution at which plaques
remained visible. The primary outcome is binary: an interaction is lytic when
the score exceeds 0.

The summary matrix was audited against the 318,816 replicate-level spot tests.
Each spot is binary, with the literal value `n` marking a reading that could not
be made, treated as missing rather than negative. A dilution step was counted as
positive when plaques were recorded in at least two readable replicates, and the
reconstructed score taken as the deepest positive step on the 0, -1, -2, -4
ladder. This reproduces the binary lytic call for 99.69 percent of the 35,268
comparable cells, with all 28,124 released zeros recovered as zeros, and
correlates with the released ordinal score at Spearman 0.991. The replicate
table covers only the 369 strains of the original collection, so the 33 strains
contributed as phage isolation hosts are excluded from that comparison rather
than counted as disagreements.

## Receptor compatibility

No table in the release states which host surface receptor each phage requires.
Compatibility is therefore inferred from the one piece of receptor evidence
independent of the interaction outcomes under analysis: the strain each phage
was isolated on. A strain that yielded a phage necessarily presented a surface
that phage could bind. Isolation hosts are pooled across the members of a
VIRIDIC genus, whose receptor binding proteins are homologous, giving a
permissive surface state set per genus. A strain-phage pair passes the gate when
the strain's surface state appears in the permissive set of that phage's genus.

Three gate tiers were defined before any association was tested: LPS outer core
alone, LPS outer core with O-antigen status, and LPS outer core with full
O-antigen serotype. A serotype recorded as `ONT` is not typeable rather than
absent and is left missing. Capsule state is absent for the 33 isolation-host
strains and is therefore not used in the gate.

Gate performance was then measured, and reported whether or not it was
favourable. The LPS gate admits 65.7 percent of pairs and raises the lysis rate
from 16.98 to 22.72 percent, an enrichment of 1.34 fold, attributing 35.92
percent of failed infections to surface exclusion. Enrichment varies from 3.16
fold (Felixounavirus) to 1.05 fold (Vectrevirus). Gate stringency also tracks
the number of phages in a genus, because a genus represented by one phage yields
one permissive state, which is a structural artefact unrelated to receptor
biology.

Because the gate is weak, the primary analysis conditions on surface state by
stratification rather than by filtering. The gate is retained as a secondary
analysis and its weakness reported as a finding.

## Association testing

A defense system is tested when at least 20 of the 402 strains carry it, which
leaves 33 of the 137 annotated systems. The threshold was fixed before testing;
below it, most phage strata contain no carrier and a stratified estimate rests
on a handful of tables.

Effects are estimated along a conditioning ladder:

| Level | Estimator | Stratification |
|---|---|---|
| L0 | Fisher exact odds ratio, Woolf interval | none |
| L1 | Cochran-Mantel-Haenszel | phage genus |
| L2 | Cochran-Mantel-Haenszel | phage identity |
| L3 | Cochran-Mantel-Haenszel | phage identity and LPS outer core |
| L4 | Cochran-Mantel-Haenszel | L3 plus O-antigen status |
| L5 | Cochran-Mantel-Haenszel | phage identity, inside the receptor gate |

L3 is the primary receptor-conditioned estimate. Mantel-Haenszel standard errors
use the Robins, Breslow and Greenland estimator, which is consistent both for
few large strata and many sparse ones. Strata in which either exposure or
outcome does not vary carry no information about the common odds ratio and are
dropped; the number retained is reported with every estimate. Benjamini-Hochberg
correction is applied within each conditioning level, each being a separate
family of 33 simultaneous tests.

Quantitative shifts in the ordinal plaque score are tested with the
Mann-Whitney U test and the Brunner-Munzel test. Brunner-Munzel is reported as
the primary rank test because it does not assume equal variances or equal
distribution shapes, which matters when defense carriage compresses the upper
tail without shifting the median. The reported `p_hat` is the probability that a
randomly chosen non-carrier pair shows a higher score than a randomly chosen
carrier pair, ties counted as half. Two scopes are given: all pairs within the
gate, and pairs that produced plaques. The second isolates the magnitude of a
successful infection but conditions on the outcome and can therefore carry
selection bias; it is labelled accordingly.

## Variance attribution

Nested logistic regressions measure how much of the lysis signal each block of
variables explains: phage identity, host surface (LPS outer core and O-antigen
status), the 33 defense systems, and phylogroup. McFadden pseudo R2 increments
and likelihood ratio statistics are reported. These are descriptive. Pairs
sharing a bacterial genome are correlated, so the likelihood ratio statistics are
anticonservative and are used only to order the blocks by contribution.

## Population structure control

Core-genome principal coordinates are obtained by classical multidimensional
scaling of the published pairwise distance matrix: double centring of the matrix
of squared distances followed by eigendecomposition. This is a deterministic
linear projection with no fitted parameters, chosen in preference to a nonlinear
embedding so that nothing in the pipeline depends on an optimisation. Eight axes
carry 89.6 percent of the variation in core-genome distance.

Clonal clusters are single-linkage groups of the same distance matrix cut at
0.0005, giving 224 clusters among the 402 strains, against 163 Sequence Types
and 15 phylogroups. A sweep across thresholds from 1e-5 to 1e-2 is reported so
the choice can be judged.

Each system is then refitted as a binomial regression with phage identity and
LPS outer core as fixed effects, under four adjustments: none, phylogroup,
principal coordinates, and both. Standard errors are clustered on clonal
cluster throughout, because a strain contributes up to 96 correlated
observations and strains within a clonal cluster are near-replicates. A
generalized linear mixed model with a random intercept per strain was
considered; the cluster-robust fit is reported instead because it assumes
nothing about the shape of the between-lineage distribution and is stable at
this number of fixed effects. Models are fitted by Newton iteration on the
logistic likelihood rather than by iteratively reweighted least squares, which
exhausts the available LAPACK workspace at this sample size.

A joint model containing all 33 systems simultaneously, with the same lineage
covariates, is reported alongside the per-system fits.

## Phage counter-defense

The 96 phages were not annotated individually for anti-defense content; doing so
would require re-annotating 96 genomes with HMMER and MacSyFinder, which the
platform and storage budget did not allow. Instead the published
AntiDefenseFinder survey of RefSeq phage genomes supplies, for each ICTV genus,
the prevalence of each anti-defense family among its sequenced relatives. All 19
genera in this collection appear in that survey, with 2 to 207 surveyed genomes
each; genera with fewer than 10 are excluded from the rescue analysis. Every
resulting statement is about lineages, not about which genes a given phage
carries.

A defense system is paired with the anti-defense family that counteracts it
(restriction-modification subtypes with Anti_RM, CRISPR-Cas with Anti_CRISPR,
CBASS with Anti_CBASS, toxin-antitoxin systems with Anti_TA, and so on). Genera
are split at the median prevalence of that family among adequately surveyed
genera. Systems with no published counterpart are not tested; an earlier design
fell back to total anti-defense burden, which was abandoned because prevalence of
any anti-defense family saturates at 100 percent in most genera and leaves a
median split with nothing on one side.

The test is an interaction, not a comparison of lysis between phage groups.
Anti-defense status is a property of the phage, so comparing lysis across
anti-defense-rich and anti-defense-poor phages would mostly measure host range.
The model instead estimates the protective effect of the host system separately
against each group, with phage identity absorbed as a fixed effect. Host defense
carriage varies within each phage, so the interaction is identified while phage
host range is fully controlled. Host surface type, phylogroup and principal
coordinates are included, and standard errors are clustered on clonal cluster.

## Direct annotation of the phage genomes

The 96 phage genomes were annotated for anti-defense content rather than having
it imputed from genus. Protein coding genes were predicted with pyrodigal, a
binding of the Prodigal algorithm, run in metagenomic mode because phage genomes
are short and compositionally unlike their hosts. The predicted proteins were
searched with pyhmmer, a binding of HMMER3, against the 178 curated anti-defense
profiles distributed in defense-finder-models 3.1.0.

Hits were accepted at each profile's own gathering cutoff, the curated
per-profile bit score threshold AntiDefenseFinder relies on, so no E-value was
chosen by hand. System calls then applied each model's rules as given in its XML
definition: the number of mandatory genes required and the minimum number of
genes required. Of the 156 anti-defense definitions, 150 require a single gene,
so for almost all families a hit above the gathering cutoff is the system call.
For the six requiring two or more genes, the model's inter_gene_max_space
constraint was applied over Prodigal gene order, as MacSyFinder does.

The annotation was validated against the published RefSeq survey by comparing
per-genus prevalence computed from these 96 genomes with prevalence among
sequenced relatives of the same genera.

The rescue test was then repeated with the measured calls in place of the genus
prior, using an identical estimator, so the two are directly comparable.

## Robustness analyses

Four designs that do not depend on the covariate adjustment were used to stress
the surviving associations.

**Within clonal cluster.** For each system, the data were restricted to clusters
containing the system both present and absent, since no other cluster carries
information about it, and a binomial regression was fitted with phage and
cluster fixed effects and standard errors clustered on cluster. The cluster
fixed effect absorbs everything the shared genomic background contributes. A
Mantel-Haenszel estimate over the same strata is reported beside it for
comparison; it treats pairs as independent and is miscalibrated, which the
internal check below makes visible.

**Permutation inside a lineage.** Carrier labels were shuffled among strains of
the same phylogroup, and separately of the same Sequence Type, holding the
number of carriers fixed, with 10,000 permutations. The statistic is the
difference in mean per-strain lysis rate between carriers and non-carriers;
every strain was challenged with the same 96 phages, so per-strain rates are
directly comparable. Because shuffling happens inside the lineage, a system that
merely marks a lineage cannot produce a low p-value.

**Minimum detectable effect.** For each system, the protective odds ratio its
own fitted standard error would have detected at 80 percent power and a
two-sided alpha of 0.05, computed as exp(-(1.96 + 0.84) times the standard
error). This uses the standard error the design actually produced, including the
loss from clustering, rather than an idealised calculation.

**Sensitivity and unmeasured confounding.** The headline count was re-derived at
nine clonal distance thresholds. E-values of VanderWeele and Ding were computed
for each estimate, converting the odds ratio to an approximate risk ratio first
because the outcome is common at about 21 percent. Effect heterogeneity across
phage genera was summarised with Cochran's Q and I squared over inverse-variance
weighted per-genus log odds ratios.

**Internal calibration check.** Throughout the analysis, the number of systems
significantly associated with *increased* lysis is used as a calibration
statistic. No defense mechanism can make a host easier to lyse, so any such
result is confounding or miscalibration. The lineage-adjusted regression returns
none; the within-cluster regression returns none; the Mantel-Haenszel version of
the within-cluster design returns nine, which is why its p-values are reported
but not interpreted.

## Partitioning explained variation

Blocks of variables are compared by the exact Shapley value of the McFadden
pseudo R2: each block's share is the average of the increase it produces when
added to a model, taken over every ordering of the blocks. With four blocks this
requires fitting all sixteen subsets. The allocation is order-independent and
sums to the full-model pseudo R2, which an ordered sequence of increments does
not. Results are reported alongside each block fitted alone and each block added
last, which bracket the Shapley figure, and per degree of freedom, because a
block with more terms can capture more variation by flexibility alone.

The surface block contains the LPS outer core, O-antigen serotype, O-antigen
status, capsule state and the allele called at twelve outer membrane protein
loci. Levels of a categorical variable carried by fewer than five strains are
pooled into one "other" level so that a 94-level serotype does not spend ninety
degrees of freedom on singletons. Capsule columns, which are missing for the
isolation-host strains, enter as a presence indicator and a missingness
indicator rather than causing those rows to be dropped.

## Receptor binding protein grouping

The released receptor binding protein table names each protein by gene
identifier without a sequence. Those identifiers match the proteins predicted in
script 07 for 125 of 130 entries, and the matched proteins have a median length
of 677 residues against 127 for all predicted proteins, with 85 percent above
400 residues against 12 percent overall. The correspondence is therefore treated
as established.

Recovered proteins are compared all against all with phmmer. The similarity of
a pair is its bit score normalised by the self-score of the shorter sequence,
giving a value near 1 for near-identical proteins. Phages are grouped by single
linkage at five thresholds from 0.20 to 0.80, each phage taking the group of its
longest annotated receptor binding protein. Gates are then rebuilt exactly as in
script 02, pooling isolation hosts within receptor group instead of within
genus, and compared on the same footing. Grouping uses only phage protein
sequences and the gate uses only isolation hosts, so the construction remains
independent of the interaction outcomes.

## Joint estimation across systems

Every system carried by at least five strains is first fitted separately with
the lineage adjustment of script 04. The resulting log odds ratios are combined
in a random effects model, with the between-system variance estimated by the
Paule-Mandel iteration and DerSimonian-Laird reported alongside it. Each
system's estimate is shrunk toward the pooled mean in proportion to its own
precision, and the posterior interval is formed from the shrunken variance.
Estimates with a standard error above 5 on the log odds scale are treated as
degenerate fits and excluded rather than shrunk.

The between-system variance is reported as a quantity of interest in its own
right: it measures how much genuine variation in effect exists between defense
systems once sampling error is removed.

## Combinations

Pairs of systems are fitted with both main effects, their interaction, phage
identity, host surface type and the clonal background covariates, with standard
errors clustered on clonal cluster. A pair is tested only when all four carriage
combinations occur in at least ten strains. Systems are additionally grouped
into seven mechanistic classes by name prefix, and the number of distinct
classes a strain carries is compared with the raw system count as a predictor of
its lysis rate.

## Consolidation

The union of systems nominated in the protective direction by the lineage
adjusted regression, the within-cluster design, either permutation scheme or the
shrinkage model is assembled, and the full battery is applied to all of them,
including rare systems that had not previously been tested. Ranking is by the
number of independent checks passed rather than by any single p-value, because
the checks disagree often enough that any one of them alone would mislead.

## Host genome retrieval and identity validation

The bacterial genomes are not part of the released data. Strain names were
resolved individually against NCBI, requiring the candidate to be Escherichia
and to carry the strain name in the assembly's own strain field rather than
anywhere in the record, preferring RefSeq and more complete assemblies.

Identity is then verified rather than assumed. Every genome is sketched with a
FracMinHash sketch (k = 21, one hash kept in a hundred, implemented directly
with a SplitMix64 finaliser over canonical 2-bit encoded k-mers) and all
pairwise Mash distances computed. A genome is accepted when its length is
plausible for Escherichia and when the strain the released core-genome matrix
calls nearest also falls within the ten nearest by k-mer distance, which has
probability 6.2 percent by chance among 162 genomes.

A distance-profile rank correlation was tried first as the gate and rejected:
among genomes whose nearest neighbour matches exactly, it still ranges from
below zero to 0.95, because most pairs of Escherichia strains sit in a narrow
band of core-genome distance. It is reported as a diagnostic.

## Prophage detection

Two measures. The specific one is the fraction of each phage's k-mers present
in each host genome, from the same sketches. Exact k-mer containment detects
prophages closely related to a collection phage and falls away quickly with
divergence, which is the right behaviour for superinfection exclusion because
repressor-mediated immunity requires close relatedness. It varies within a
strain across phages, so the contrast is identified with the host's genetic
background held fixed.

The general one counts tyrosine recombinase domains, with genes called by
Prodigal trained on each genome rather than in metagenomic mode, which is both
appropriate for whole genomes and sixteen times faster. That count is an index
of mobile element content and not a prophage count: the same profiles match
XerC, XerD and integron recombinases.

## Phylogenetically matched nulls

For each system, the share of carriage variance explained by the eight
core-genome principal coordinates is measured, then null traits are simulated
with that same share and the same prevalence: a random direction through the
coordinates supplies the structured part, white noise the rest, mixed to the
measured ratio and thresholded at the observed carrier count. The observed
association is compared with 10,000 such traits.

The test conditions on ancestry and on nothing else, which is not sufficient:
host surface type explains 15.3 percent of per-strain lysis rate and is only
partly collinear with ancestry. The same null is therefore also run on lysis
rates with LPS outer core and O-antigen status regressed out.

## Full Bayesian pooling

The between-system standard deviation is given a half-Cauchy prior of scale
0.5 and the pooled mean a flat prior. The pooled mean is integrated out
analytically at each value of tau, leaving a one-dimensional integral evaluated
on a 2000-point grid; each system's posterior is a mixture over that grid, so
uncertainty in both the pooled mean and tau is carried into every interval. No
Markov chain is used and nothing is approximated by sampling.

## Ordinal thresholds

Each system is fitted at all four plaque-score thresholds with identical
covariates and clustered standard errors. The thresholds are nested and reuse
the same observations, so no combination of them is treated as a test; a
weighted average is reported as an effect size and the per-threshold estimates
carry the honest intervals. Reporting them separately also tests the
proportional odds assumption that a single ordinal model would impose.

## Accessory genome association

Every variable k-mer in the validated genomes is tested, with no catalogue.
K-mers are sampled at one in twenty with the same FracMinHash scheme used
elsewhere, and those carried by fewer than ten or more than all-but-ten strains
are dropped as uninformative. K-mers occurring in exactly the same strains are
collapsed into one presence pattern, which is both what makes the problem
tractable and the right treatment of linkage; patterns backed by fewer than
five k-mers are discarded as too short to be loci.

Significance is assessed against an empirical null rather than a simulated one.
Within bins of carrier count and phylogenetic signal, the spread of the test
statistic across loci, which cannot nearly all be causal, is taken as the null
and summarised by its median and scaled median absolute deviation so that a few
true associations cannot inflate it. A bin whose spread collapses, which happens
where near-identical clade markers dominate, carries no usable null and is
skipped.

Fixed bins leave most loci untested, so the same null is also computed on an
adaptive neighbourhood: each locus is matched to its 400 nearest neighbours in
the two covariates, ranked to quantiles so that neither dominates the metric.
The reference distribution is not assumed normal, which the data reject for the
loci fixed bins exclude. Instead the 98th percentile of the standardised
statistic within the neighbourhood is matched to the corresponding Student-t
quantile over a grid of degrees of freedom from 2 to 50, and the locus is read
against the t distribution so obtained. Two diagnostics accompany every locus:
the local spread as a fraction of the spread across all loci, and the fraction
of the neighbourhood lying within a tenth of the global spread of the local
centre. A locus flagged in a neighbourhood that is collapsed and clumped on
those measures is reported with its diagnostics and not interpreted, since such
a neighbourhood estimates a null from one set of near-identical clade markers. The simulated null of script 21 is computed alongside and reported, not
interpreted: it is about 2.8 times too liberal here, because traits built from
eight principal coordinates are less clade-confined than real accessory loci.

Flagged patterns are then relocated in a carrier genome by rescanning it for the
pattern's k-mer hash, a window of three kilobases either side is extracted, genes
are called and searched against the anti-defense and lysogeny profiles, and the
carriers are checked for catalogued defense systems that non-carriers lack.

## Phage-resolved effects

Each defense system is tested against each phage separately, as a Rao score
test for adding the system to a logistic model of lysis fitted for that phage
with an intercept and eight core-genome principal coordinates. The covariate
model is fitted once per phage and every candidate trait is then scored against
it in closed form, which is what makes a simulation-based null affordable
across thousands of traits.

The score form is also the only one of the three classical tests that stays
defined here. In 317 of the 3,102 system-by-phage tests, a tenth of them and
touching 27 of the 33 systems, every carrier shares the same outcome. A
logistic fit then has no finite maximum likelihood estimate, so Wald and
likelihood-ratio statistics are unusable, while the score statistic is
evaluated under the null where the coefficient is zero and remains well
behaved.

Phages are collapsed into groups by single linkage on patristic distance from
the VIRIDIC tree at a threshold of 0.1, and the unit counted is the group
rather than the phage, because near-identical phages are not independent tests.
The headline count is reported at 0.05, 0.1 and 0.2 so that it does not rest on
one cut.

Two nulls are generated per system, each with 1000 traits. The first matches prevalence and the share of carriage variance
explained by the principal coordinates, bisecting on the share the simulated
binary traits actually realise rather than on the latent mixture behind them,
because thresholding a latent mixture loses enough structure to miss the
target by half for the most lineage-confined systems. The second permutes the number of carriers in each clonal cluster
among clusters of the same size, which preserves prevalence and the whole
within-cluster carriage pattern while randomising which lineages carry the
trait. The second exists because the first does not reproduce the clonal
coherence of real defense systems, measured as the fraction of same-cluster
strain pairs agreeing on carriage; that measure is used rather than a
regression on cluster indicators, which would overfit at 224 clusters for 402
strains. Neither null dominates, so every reported p-value is the larger of the
two.

Both tails are tested, not only the protective one, since an association with
increased lysis cannot be a defense mechanism and that pattern recurs in this
collection. Everything below is applied to whichever direction a system
survives in.

A system that survives is also re-tested with the collection's typing of the
host surface added to the covariates, since the per-phage model carries only
core-genome position and surface explains about three times the variation
defense systems do. For a system acting after adsorption the surface is a
confounder rather than a mediator, so conditioning on it is the appropriate
question rather than an over-adjustment.

A system that survives is re-tested with the four systems most enriched among
its carriers added to the covariates, since defense systems travel together on
islands and a survivor may be carrying the credit for a neighbour. It is then
refitted with each clonal cluster holding one of its carriers dropped in turn.
Because losing carriers costs power whatever they are, the most damaging of
those drops is compared against 40 drops of the same number of strains and the
same number of carriers chosen at random: a result that collapses only when a
particular cluster goes is resting on that group of relatives, while one that
collapses equally under random drops is merely thin.

The clonal null is checked for degeneracy, since a null that reproduced the
clonal pattern exactly might have no freedom left to move it. Agreement between
each null trait and the real carriage vector is compared with the agreement
that the trait's prevalence forces on its own.

Permutation blocks draw from separate generators seeded from the project seed,
so that adding an analysis does not shift the p-values of the ones beside it.

The procedure was measured against known answers rather than assumed correct.
Applying it to traits built to have exactly the clonal structure of a real
defense system but no relation to lysis, and judging each against fresh
independent draws, it rejects at 1.5 percent where 5 percent is nominal, and at
7.0 and 16.0 percent where 10 and 20 are nominal: it is conservative
throughout. Its median p-value is 1, which is a property of the statistic
rather than a fault, since a count of phage groups is zero for most null traits
and a p-value built from a count is 1 whenever the count cannot be exceeded.

Power was measured the same way, by suppressing lysis among carriers for twelve
phages drawn from distinct groups with lysis rates between 0.41 and 0.65. A
suppression producing a difference of 0.267 on those phages is detected in 77
percent of replicates, one producing 0.355 in all of them, and one producing
0.470 in all of them. The one candidate this analysis reports sits at 0.222,
just below the point where power is 77 percent, so the test is working in a
range where it can see effects of the size claimed; effects much smaller than
that on a subset of phages would not be reliably detected, which is the honest
bound on the negative results beside it.

## Software

Python 3.13.11 with numpy 2.3.5, pandas 2.3.3, scipy 1.16.3, statsmodels 0.14.6,
matplotlib 3.10.8, seaborn 0.13.2 and biopython 1.86, plus pyrodigal 3.7.1 and
pyhmmer 0.12.3 for the genome annotation. No predictive model is fitted anywhere
in the pipeline. Every estimate comes from a closed-form contingency table test,
a rank test, a generalized linear model fitted by maximum likelihood, a
permutation null, or an eigendecomposition of a precomputed distance matrix.
Gene prediction and profile hidden Markov model search are deterministic
algorithms with fixed published parameters, not fitted models.
