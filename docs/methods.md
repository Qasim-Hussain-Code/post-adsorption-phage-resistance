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

## Software

Python 3.13.11 with numpy 2.3.5, pandas 2.3.3, scipy 1.16.3, statsmodels 0.14.6,
matplotlib 3.10.8, seaborn 0.13.2 and biopython 1.86, plus pyrodigal 3.7.1 and
pyhmmer 0.12.3 for the genome annotation. No predictive model is fitted anywhere
in the pipeline. Every estimate comes from a closed-form contingency table test,
a rank test, a generalized linear model fitted by maximum likelihood, a
permutation null, or an eigendecomposition of a precomputed distance matrix.
Gene prediction and profile hidden Markov model search are deterministic
algorithms with fixed published parameters, not fitted models.
