# Results

All numbers below are produced by the committed scripts and appear in
`results/` and `results/tables/`. Negative and null findings are reported with
the same prominence as positive ones.

## 1. The interaction landscape

Of 38,435 observed strain-phage interactions, 7,976 (20.75 percent) produced
plaques and 30,459 (79.25 percent) did not. The failure rate is close to the 80
to 85 percent often quoted for surface-level failure across a species, but that
coincidence is not evidence about where the failures occur, and the gate
analysis below shows most of them are not at the surface.

Host surface chemistry plainly matters. Lysis rate varies with LPS outer core
type by a factor of 2.4, from 10.76 percent on R2 to 26.09 percent on R1.
Phage lineage matters more: Myoviridae lyse 41.30 percent of the pairs they are
tested in, Podoviridae 11.20 percent and Siphoviridae 9.14 percent, and across
VIRIDIC genera the lysis rate spans from 0.75 percent (Lambdavirus) to 56.61
percent (Mosigvirus).

## 2. Surface exclusion explains a minority of failures

The isolation-host receptor gate admits 65.7 percent of pairs. Inside the gate
22.72 percent of pairs are lysed (95 percent CI 22.21 to 23.24); outside it
16.98 percent (16.35 to 17.63). The enrichment is 1.34 fold, and the gate
captures 71.9 percent of all lysis events.

On this definition, **35.92 percent** of failed infections are attributable to
surface incompatibility and **64.08 percent** occur inside a host whose surface
state is compatible. Requiring the full O-antigen serotype to match as well
raises the apparent surface share to 88.70 percent, but only by discarding 87.7
percent of the data, and enrichment barely moves (1.39 fold). That is the
signature of a filter removing signal and noise together rather than sharpening
a contrast.

Gate performance is heterogeneous across genera. It separates cleanly in
Felixounavirus (57.77 percent inside against 18.27 outside, 3.16 fold) and
Phapecoctavirus (1.79 fold), and not at all in the two largest genera,
Vectrevirus (1.05 fold) and Kayfunavirus (1.20 fold). Two genera show
enrichment below 1. A gate that works for a minority of genera cannot carry a
species-wide claim, which is why the primary analysis stratifies on surface
state rather than filtering on it.

## 3. Conditioning on surface state widens defense effects

Median absolute log2 odds ratio across the 33 testable systems rises
monotonically with conditioning depth:

| Level | Conditioning | Median abs log2 OR | Significant at q <= 0.05 | Mean strata |
|---|---|---|---|---|
| L0 | none | 0.282 | 26 of 33 | 1 |
| L1 | phage genus | 0.330 | 26 | 19 |
| L2 | phage identity | 0.343 | 26 | 96 |
| L3 | phage and LPS outer core | 0.469 | 25 | 449 |
| L4 | plus O-antigen status | 0.517 | 24 | 560 |
| L5 | phage identity inside the gate | 0.567 | 25 | 95 |

Conditioning on phage identity and surface type moves the median from 0.282 to
0.469, a factor of 1.66. Effects grow rather than shrink, the direction the
surface-masking hypothesis predicts: pooling across phages and host surface
types dilutes defense associations toward the null. The magnitude is modest, so
surface variation masks these associations partially rather than dominantly.
Four of 33 systems reverse direction between the unstratified and
receptor-conditioned estimate.

At L3, 11 systems are associated with reduced lysis (Zorya_TypeII 0.348,
RM_Type_IIG 0.432, CAS_Class1-Subtype-I-E 0.543, Hachiman 0.571, Kiwa 0.582,
Shango 0.612, Mokosh_Type_I_A 0.656, PD-Lambda-4 0.729, Mok_Hok_Sok 0.786,
Gao_Tmn 0.786, RM_Type_III 0.861) and **14 with increased lysis** (up to SanaTA
1.973 and PsyrTA 1.771). No defense mechanism can make a host easier to lyse, so
the second group is the clearest available warning that these marginal
associations do not measure defense function on their own.

## 4. What actually explains lysis

Nested logistic decomposition of the binary outcome:

| Model | Parameters | Pseudo R2 | Increment |
|---|---|---|---|
| intercept only | 1 | 0 | |
| plus phage identity | 96 | 0.1757 | 0.1757 |
| plus host surface | 102 | 0.2065 | 0.0308 |
| plus 33 defense systems | 135 | 0.2497 | 0.0432 |
| plus phylogroup | 149 | 0.2612 | 0.0114 |

Which phage is tested dominates. Host surface type and defense carriage
contribute comparably, with defense slightly the larger of the two, and neither
approaches the phage term. The premise that surface incompatibility accounts for
the overwhelming majority of variation in this kind of screen is not supported
here: a large part of the variation is phage-specific host range that does not
reduce to the surface determinants typed in this collection. These are
descriptive quantities; the correlated structure of the data makes the
likelihood ratio statistics anticonservative.

## 5. Clonal background accounts for almost all of it

Refitted as regressions with standard errors clustered on clonal cluster, 6 of
33 systems are significant before any lineage covariate. Adding phylogroup drops
this to 2; adding principal coordinates alone to 1; adding both to 2. Median
absolute log2 odds ratio falls from 0.440 to 0.225.

Two systems survive with direction preserved:

| System | CMH OR (L3) | Cluster-robust, unadjusted | Fully adjusted | 95 percent CI | q |
|---|---|---|---|---|---|
| Zorya_TypeII | 0.348 | 0.401 | 0.444 | 0.282 to 0.700 | 0.0077 |
| RM_Type_IIG | 0.432 | 0.471 | 0.549 | 0.398 to 0.757 | 0.0077 |

Four systems significant before adjustment are lost: CAS_Class1-Subtype-I-E,
Kiwa, PsyrTA and SanaTA. Median retained fraction of the log2 effect among
systems keeping their direction is 0.56.

The adjustment passes an internal check. Fourteen systems were significantly
associated with increased lysis under the stratified estimates and two under the
unadjusted regression. After adjusting for phylogroup and core-genome position
there are **none**. The covariates remove exactly the kind of signal they were
introduced to remove, which makes the attenuation of the protective associations
harder to attribute to overcorrection.

In a joint model with all 33 systems mutually adjusted and the same lineage
covariates, 3 are significant: RM_Type_IIG (0.457), Shango (0.490) and SanaTA
(1.911). Defense systems co-occur, so this is a stricter and more collinear
contrast; RM_Type_IIG is the only system significant in both framings.

## 6. Defense burden: a clean null

Lysis rate by quartile of systems carried is 19.4, 22.3, 21.4 and 19.9 percent.
Across the full range of 1 to 16 systems there is no monotone trend. The
strain-level Spearman correlation between burden and per-strain lysis rate is
+0.053 (p = 0.29). The pair-level correlation is +0.022 (p = 1.4e-05), nominally
significant only because 402 strains contribute 38,435 correlated observations,
and is reported alongside the strain-level figure for that reason.

Counting systems treats a type I restriction enzyme and a Thoeris operon as
interchangeable units of protection. The absence of any trend suggests they do
not behave that way.

## 7. Quantitative plaque score

Within the gate, 25 of 33 systems shift the score distribution at q <= 0.05 by
Brunner-Munzel; restricted to pairs that produced plaques, 8 of 33 do. Direction
agrees with the binary analysis throughout. The largest shifts toward weaker
plaquing are RM_Type_IIG (p_hat 0.570, mean score 0.189 against 0.525),
CAS_Class1-Subtype-I-E (0.561) and Zorya_TypeII (0.558). The largest shifts the
other way are Lamassu-Mrr (0.437), Gao_Qat (0.439) and Septu (0.440), the same
systems that fail to survive lineage control.

## 8. Phage counter-defense

Of 13 systems with a published anti-defense counterpart and a usable genus
split, **one** interaction is significant: RM_Type_IV is less protective against
genera whose sequenced relatives commonly carry anti-RM genes (interaction odds
ratio 1.626, 95 percent CI 1.213 to 2.180, q = 0.015). Four of the five RM
subtypes have an interaction odds ratio above 1 (RM_Type_IV 1.63, RM_Type_III
1.40, RM_Type_I 1.29, RM_Type_II 1.17, RM_Type_IIG 0.91). Anti-RM proteins such
as Ocr and ArdA are the best characterised and most widespread anti-defense
proteins known, so this is the pairing where a lineage-level prior should carry
most information. The pattern is suggestive. It rests on one significant test
out of thirteen and is not established here.

Fifteen systems have no published anti-defense counterpart and were not tested.
Five have a counterpart that gives no usable split, including both CRISPR-Cas
subtypes: Anti_CRISPR, Anti_Gabija and Anti_Hachiman are absent from every
adequately surveyed genus in this collection, so the systems they counteract
cannot be tested for rescue at all.

Across genera, lysis rate and total anti-defense prevalence correlate at
Spearman +0.613 over 19 genera. With this few genera, and with lysis rate
varying by nearly two orders of magnitude for reasons of host range, this is not
interpretable as causal in either direction.

## 9. Residual resistance

| Subset | Pairs | Failure rate |
|---|---|---|
| all observed pairs | 38,435 | 79.25 percent |
| receptor gate passed | 25,256 | 77.28 percent |
| gate passed, lacks both lineage-validated systems | 22,436 | 75.83 percent |
| gate passed, lacks all 33 tested systems | 0 | not estimable |

Applying the surface gate and removing hosts carrying either validated defense
system moves the failure rate by 3.4 percentage points. The great majority of
failed infections in this collection are not attributable to the surface
determinants typed here, nor to the intracellular defense systems DefenseFinder
detects, nor to anti-defense content imputed from phage lineage.

The last row is a real constraint rather than an oversight. Every one of the 402
strains carries at least one of the 33 testable systems, MazEF alone being
present in 338 of them, so there is no defense-free comparison group and a
question of the form "how often does infection fail in a host with no defenses"
cannot be answered in this collection.

## 10. Annotating the phages directly

The counter-defense analysis above imputed anti-defense content from a phage's
genus. The 96 genomes were subsequently annotated themselves: genes called with
Prodigal, then searched against the 178 curated anti-defense profiles at their
own gathering cutoffs, with each model's gene requirements applied.

From 10,159 predicted proteins came 191 profile hits above threshold and 149
system calls. **73 of 96 phages (76.0 percent) carry at least one anti-defense
system.** Anti_RM is in 48 of 96, Anti_TA in 18, NADP in 17, Anti_CBASS in 14.
Seven families are absent entirely, among them Anti_CRISPR, Anti_BREX,
Anti_Gabija and Anti_Hachiman.

The annotation agrees with the published RefSeq survey at the genus level
(Spearman +0.856 across 19 genera, mean absolute difference 6.0 percentage
points), which is the expected check. More usefully, four genera turn out not to
be uniform in anti-defense carriage, so the genus prior had been treating
genuinely different phages as identical.

**Re-running the rescue test on the measured calls does not change its answer.**
With a balanced 48 against 48 split for Anti_RM, one interaction of twelve is
significant: RM_Type_IV again protects less well against phages carrying anti-RM
genes (interaction odds ratio 1.538, 95 percent CI 1.168 to 2.024, q = 0.026).
Now **all five** RM subtypes have an interaction above 1, up from four of five.
Direction agrees with the imputed analysis for 11 of 12 systems and the two sets
of estimates correlate at Spearman +0.993.

This closes the largest acknowledged weakness in the earlier analysis. The
negative counter-defense result was suspected of being an artefact of
imputation, since imputation attenuates toward the null. It is not: measuring
the annotation directly yields the same number of significant interactions.

## 11. Stressing the two surviving associations

Script 04 removed clonal confounding by adjustment, which only works if the
covariates capture the structure that matters. Four designs that do not rely on
that assumption give a more complicated answer than the headline.

**Within clonal cluster.** Restricting to clusters containing a system both
present and absent, with a cluster fixed effect and errors clustered on the same
variable, 3 of 33 systems are significant and protective: RM_Type_III (0.328),
Retron_I_C (0.249) and RosmerTA (0.532). This design passes the internal check
with 0 biologically impossible significant results. A Mantel-Haenszel version of
the same stratification returns 9 such impossible results and its intervals are
not usable, which is why only the clustered regression is interpreted.

Within this design the two headline survivors diverge. **RM_Type_IIG loses its
effect entirely (odds ratio 0.968, 95 percent CI 0.487 to 1.923) across 7
discordant clusters spanning 36 strains**, which is moderate power, so this is
evidence against rather than absence of evidence. **Zorya_TypeII keeps its
effect size (0.344 against 0.444 adjusted) but only 2 clusters spanning 11
strains are informative**, so the interval runs to 1.978 and the test is
uninformative rather than negative.

**Permutation inside a lineage.** Shuffling carrier labels among strains of the
same phylogroup, RM_Type_IIG is the only system to reach significance
(q = 0.0033). Shuffling within Sequence Type, it does not (q = 0.38), while
RM_Type_III (q = 0.0066) and SanaTA (q = 0.040) do. Zorya_TypeII reaches
significance under neither (q = 0.30 and q = 0.38), with only 3 usable Sequence
Type groups.

**What the nulls exclude.** Across the 31 systems reported as null, the
protective odds ratio detectable at 80 percent power ranges from 0.43 to 0.70.
Systems at the upper end can only exclude very large effects, so most of the
nulls in this analysis are uninformative about moderate protection rather than
evidence of its absence.

**Threshold sensitivity.** The headline count is stable: RM_Type_IIG and
Zorya_TypeII are significant at all nine clonal distance thresholds tried, from
386 clusters down to 43. Coarser thresholds add Shango and then RM_Type_III, so
the chosen threshold is the conservative end of the range rather than a lucky
one.

**Unmeasured confounding.** An unmeasured confounder would need a risk ratio of
at least 3.40 with both carriage and lysis to explain Zorya_TypeII away (2.01 to
reach the confidence bound), and 2.69 for RM_Type_IIG (1.82 to the bound).
Phylogroup and core-genome position are already adjusted for, so such a
confounder would have to act beyond the lineage structure they capture.

**Consistency across phage genera.** RM_Type_IIG has an odds ratio below 1 in
**15 of 15** estimable genera (pooled 0.422, I squared 86 percent). Zorya_TypeII
is below 1 in 11 of 13 (pooled 0.516, I squared 81 percent). Both are
heterogeneous in magnitude but consistent in direction, which is hard to produce
by chance.

**Synthesis.** No system passes all four stricter tests.

| System | Adjusted OR | Within-cluster OR | Clusters | Perm q (phylogroup) | Perm q (ST) | Genera below 1 | E-value | Tests passed |
|---|---|---|---|---|---|---|---|---|
| RM_Type_IIG | 0.549 | 0.968 | 7 | 0.0033 | 0.38 | 15/15 | 2.69 | 2 of 4 |
| RM_Type_III | 0.626 | 0.328 | 12 | 0.30 | 0.0066 | n/a | 2.31 | 2 of 4 |
| Zorya_TypeII | 0.444 | 0.344 | 2 | 0.30 | 0.38 | 11/13 | 3.40 | 1 of 4 |
| Retron_I_C | 0.722 | 0.249 | 4 | 0.30 | 0.78 | n/a | 1.94 | 1 of 4 |
| RosmerTA | 1.171 | 0.532 | 11 | 0.76 | 0.90 | n/a | 1.52 | 1 of 4 |

Each design discards most of the data to remove a different source of
confounding, so none retains much power, and disagreement between them is
expected rather than alarming. Still, the pattern matters: **RM_Type_III, which
did not survive the lineage-adjusted regression, is significant and protective
in two of the stricter designs**, and is arguably the best-supported candidate
in the whole analysis. It is not a claim, but it is where a targeted experiment
would start.

## 12. Correcting the variance decomposition

Section 4 compared blocks of variables by adding them to a model in one fixed
order. That comparison was defective twice over: the increments depend on the
order, and the surface block held only six terms against thirty-three for
defense. Script 11 repeats it with a surface block containing everything the
collection types about the cell surface (LPS outer core, O-antigen serotype and
status, capsule state, twelve outer membrane protein loci) and allocates the
shared variation by exact Shapley value over all orderings, which requires
fitting all sixteen subsets of four blocks.

| Block | Terms | Shapley pseudo R2 | Share | Alone | Added last |
|---|---|---|---|---|---|
| phage identity | 95 | 0.1920 | 59.7 percent | 0.1751 | 0.2049 |
| host surface | 157 | 0.0830 | 25.8 percent | 0.1010 | 0.0685 |
| defense systems | 33 | 0.0269 | 8.4 percent | 0.0410 | 0.0179 |
| clonal lineage | 22 | 0.0199 | 6.2 percent | 0.0387 | 0.0055 |

**This reverses the conclusion of section 4.** Host surface accounts for about
**three times** the variation that all thirty-three defense systems do, and the
ordering holds whether each block is fitted alone, added last, or allocated by
Shapley value. It also agrees with the published analysis of the same dataset,
which reported adsorption-related features as the more informative ones.

One qualification cuts the other way: per term, defense is the more efficient
block (0.00081 against 0.00053), because surface spends 157 terms to defense's
33. The fair statement is that host surface explains substantially more of the
variation in total, using more parameters, and that both are dwarfed by which
phage is being tested.

## 13. A receptor gate built from the actual receptor

The weakest step in this analysis was always the definition of receptor
compatibility, which script 02 anchored on isolation hosts pooled by VIRIDIC
genus. Genus is only a proxy; the receptor binding protein is what contacts the
host surface. The released table names each phage's receptor binding protein
without giving a sequence, but those identifiers match the proteins predicted in
script 07 for 125 of 130 entries, and the matched proteins have a median length
of 677 residues against 127 for all predicted proteins, which is what tail
fibers look like.

Grouping the 87 phages with a recovered protein by receptor binding protein
homology gives 30 to 45 receptor groups depending on the similarity threshold,
against 19 genera. **It does not improve the gate.** The best receptor-based
gate reaches 1.26 fold enrichment against 1.34 fold for the genus gate.

That is a useful negative result because it localises the problem. Since the
receptor binding protein is the right biological object, the limitation is not
the grouping but the **anchor**: one isolation host per phage can establish that
a surface state is permissive and can never establish that another is not, so
the permissive sets stay small and the gate stays blunt however the phages are
grouped. The receptor-based gate attributes a larger share of failures to the
surface (60.6 against 35.9 percent) only because it is stricter, admitting 40.9
percent of pairs against 65.7 percent while separating them less well. That is a
stringency artefact, not a better measurement.

## 14. All 79 systems estimated jointly

Scripts 03 and 04 tested only the 33 systems carried by at least 20 strains.
Script 13 fits every system carried by at least five, 79 in all, and combines
them in a random effects model so that each estimate is shrunk toward the pooled
mean by an amount set by its own precision.

- Pooled odds ratio across all systems: **0.923** (95 percent CI 0.847 to 1.005)
- Between-system standard deviation: 0.269 on the log odds scale, I squared 56
  percent (Cochran Q = 178.9 on 78 df, p = 7e-10)
- Implied spread: 95 percent of true system effects lie between odds ratios of
  0.544 and 1.563
- Median shrinkage weight 0.511; 7 systems significant alone, 9 with a shrunken
  interval excluding 1

So there is real variation between systems beyond sampling error, but the
average system has essentially no association with lysis. Extending coverage
surfaces four candidates the 20-carrier threshold had hidden: **PD-T7-4**,
**AbiD**, **Lamassu-Cap4_nuclease** and **BstA** (the last associated with more
lysis, not less). Their median shrinkage weight is 0.382, so most of their
individual signal is replaced by the pooled mean, which is the correct treatment
of an estimate from a handful of strains.

## 15. Combinations do not rescue the burden null

Carrying more systems does not help (section 6). That null assumes systems act
additively, so script 14 tests pairs directly. Of 528 possible pairs, 135 have
at least ten strains in each of the four carriage combinations and can support
an interaction term. **None shows a significant interaction** after clonal
adjustment and FDR correction. Grouping systems into seven mechanistic classes
does not help either: the correlation between per-strain lysis rate and the
number of distinct classes carried (rho = +0.056, p = 0.26) is no better than
for the raw count (rho = +0.053, p = 0.29).

The burden null is therefore not hiding combination effects at the pairwise
level. Three-way and higher combinations cannot be estimated in a collection of
402 strains with repertoires this correlated.

## 16. One consolidated candidate list

Four analyses nominate different systems. Script 15 takes the union in the
protective direction, applies the same battery to all of them including the rare
ones that had never been through it, and ranks by how many independent checks
each passes.

| System | Strains | Adjusted OR | Shrunk OR | Within-cluster OR | Clusters | Genera below 1 | E-value | Checks |
|---|---|---|---|---|---|---|---|---|
| RM_Type_IIG | 52 | 0.549 | 0.632 | 0.968 | 7 | 15/15 | 2.69 | 3/5 |
| RM_Type_III | 54 | 0.626 | 0.718 | 0.328 | 12 | 10/18 | 2.31 | 3/5 |
| Shango | 37 | 0.583 | 0.705 | 0.492 | 9 | 16/17 | 2.51 | 2/5 |
| Lamassu-Cap4_nuclease | 18 | 0.451 | 0.629 | 0.417 | 5 | 12/13 | 3.34 | 2/5 |
| Zorya_TypeII | 20 | 0.444 | 0.607 | 0.344 | 2 | 11/13 | 3.40 | 2/5 |
| AbiD | 14 | 0.270 | 0.575 | 0.349 | 6 | 11/15 | 5.75 | 2/5 |
| PD-T7-4 | 7 | 0.382 | 0.571 | not estimable | 1 | 8/10 | 3.99 | 1/5 |
| Retron_I_C | 20 | 0.722 | 0.822 | 0.249 | 4 | 11/16 | 1.94 | 1/5 |
| RosmerTA | 47 | 1.171 | 1.085 | 0.532 | 11 | 6/17 | 1.52 | 1/5 |

Six systems pass at least two independent checks. **None passes all five**, and
the two headline survivors of script 04 are not clearly the best entries:
RM_Type_III matches RM_Type_IIG at three of five while having twelve informative
clonal clusters to RM_Type_IIG's seven, and AbiD has the highest E-value in the
table at 5.75 despite only fourteen carriers. Zorya_TypeII has two informative
clusters, so its within-cluster column is uninformative rather than supportive.

## 17. Is some of the failure lysogeny rather than resistance

Every analysis so far has treated a failure to plaque as host resistance. A
temperate phage entering lysogeny also leaves no plaque, and the two are
indistinguishable in a spot test. Script 17 classifies the 96 phages by
searching their predicted proteins against Pfam lysogeny profiles at the
profiles' gathering cutoffs.

Only the tyrosine integrase is used to classify. Six phages carry a serine
recombinase, but the same fold serves as a DNA invertase that switches tail
fiber genes in strictly virulent phages, so it is recorded and not used; 28
carry a repressor-like helix-turn-helix, which is common in virulent phages too.

**Two of 96 phages carry a tyrosine integrase**, one Lambdavirus and one
Wanchaivirus. They are exactly the two genera whose near-total failure to
plaque prompted the question: they form plaques on 0.9 percent of strains
against 21.2 percent for the rest (Mann-Whitney on per-phage rates p = 0.017).

So the direction supports the idea and the magnitude refutes it. Those two
phages account for 793 of 30,459 failures, 2.6 percent. Counting every one of
them as a silent lysogenisation rather than a resistance event moves the overall
failure rate only from 79.25 to 78.83 percent. The Guelin collection was
selected for lytic activity, so temperate phages are rare by construction.

Dropping them changes nothing else: the median absolute shift in log2 odds ratio
across the 33 defense systems is 0.002, the largest is 0.015 (Hachiman), and one
system of 33 changes direction.
No conclusion in this analysis rests on interactions with phages that were never
going to form a plaque.

## 18. Resident prophages: the question this collection cannot answer

Superinfection exclusion has been listed as untested throughout, because the
bacterial genomes were missing. They are not in the primary repository, the
companion Zenodo record holds plaque images, and no project accession could be
located. Scripts 18 and 19 therefore resolved strain names against NCBI
individually and validated every download.

**Validation.** 287 genomes were obtained, 71 percent of the collection, every
one matched on the assembly's own strain field. Identity was decided by whether
the strain the released core-genome matrix calls nearest is also nearest by
whole-genome k-mer distance; with 287 genomes, landing in the top ten by chance
has probability 3.5 percent. **151 genomes have the released nearest neighbour
as their single nearest k-mer neighbour**, 217 within the top three, and 254
within the top ten. Those 254 are carried forward, **63 percent of the
collection**; 33 were dropped.

The test is relative to the genomes present, and that matters. ED1a and BL21
both failed on an earlier, smaller set of 162 and both pass at rank 1 here,
because their true nearest relatives had not yet been retrieved. A genome is
only as checkable as its neighbourhood is complete.

A note on method: a distance-profile correlation was tried first as the identity
gate and abandoned. Among the 83 genomes whose nearest neighbour matches
exactly, that correlation still ranges from below zero to 0.95, because most
pairs of Escherichia strains sit in a narrow band of core-genome distance and
the rank correlation is dominated by noise among near-equidistant strains. It is
reported as a diagnostic only.

**Result.** Across **24,382** analysable strain-phage pairs, exactly **8** have
a host carrying at least 20 percent of the challenging phage's k-mers. They
involve 5 strains and just **2 phages: 411_P1 and 411_P2**. All 8 are non-lytic.

That count did not move when coverage was nearly doubled. An earlier pass over
135 validated genomes found the same 8 pairs, the same 5 strains and the same 2
phages; adding 119 more genomes added none. The rarity is a property of the
collection, not of how much of it was screened.

That looks like superinfection exclusion and is not evidence of it. Those two
phages plaque on 0.98 percent of strains overall, so the number of lysis events
expected among 8 carriers under the null is **0.08**. Observing zero is exactly
what no effect predicts. The regression is barely estimable and says so: odds
ratio 0.130 with a 95 percent interval from 0.009 to 1.844 (p = 0.13).

The limitation is the collection, not the method. The Guelin panel was assembled
for lytic activity, so only 2 of 96 phages are temperate, and those two barely
plaque. Superinfection exclusion needs a temperate phage that also forms
plaques, and this panel contains almost none.

**An internal check that does hold.** The two phages found integrated in host
chromosomes are exactly the two that script 17 classified temperate from their
own proteomes carrying a tyrosine integrase. Two analyses sharing no input
beyond the phage genomes, one searching bacterial chromosomes for phage
sequence and one searching phage proteins for integration machinery, pick out
the same pair.

**One substantive by-product.** Strains richer in site-specific recombinases
carry more defense systems (Spearman +0.266, p = 1.7e-05, n = 254). Defense systems
travel on mobile elements, so defense carriage is partly a marker of how much
mobile DNA a strain has acquired, which is one more thing a defense system can
stand in for.

## 19. A harder null, and what it does to the headline

Script 09 permuted carriage inside a phylogroup or a Sequence Type, which
requires choosing a grouping level and applies the same one to every system.
Script 21 removes the choice: for each system it measures how much of carriage
is explained by core-genome position, then simulates null traits with exactly
that phylogenetic confinement and exactly that prevalence.

**The size of the population structure problem, measured.** Median inflation of
the null standard deviation is 1.09, but it reaches **2.31 for
CAS_Class1-Subtype-I-E**, whose carriage is 59 percent explained by core-genome
position. That system goes from p = 0.001 against an unstructured null to
p = 0.079 against a matched one: a textbook lineage artifact, and the same
system script 04 lost on adjustment.

**The matched null alone is not enough.** Five systems survive it, and four
point the biologically impossible direction. Conditioning on ancestry does not
remove them because they are not lineage-confined. Host surface type is the
missing variable: it explains 15.3 percent of per-strain lysis rate and is only
partly collinear with ancestry. Gao_Qat (87 percent of carriers on a permissive
outer core, OR 4.77) and Septu (OR 2.16) are carried disproportionately by
strains that are easy to lyse whatever they carry.

Re-running the same null on surface-residualised lysis drops the impossible
survivors from 4 to 2.

**What this does to RM_Type_IIG.** It is strongly *depleted* in permissive outer
cores (OR 0.33, p = 0.0004): 37 percent of its carriers against 60 percent of
strains overall. Its carriers are drawn from strains intrinsically harder to
lyse. Removing surface type attenuates its effect by 31 percent, but it remains
significant (q = 0.0088). Adjusting for surface also brings **Zorya_TypeII**
across the threshold, which surface composition had been masking.

## 20. Proper posteriors, and the whole score

**Full Bayesian pooling** (script 22) replaces the empirical Bayes shortcut.
With a half-Cauchy prior on the between-system spread and the pooled mean
integrated out analytically on a 2000-point grid, tau has posterior median
0.281 (95 percent interval 0.205 to 0.371), comfortably containing the
empirical Bayes point estimate of 0.269. The same 9 systems exclude no effect,
so the shortcut was anticonservative in principle but harmless here. The gain is
an interpretable quantity: PD-T7-4 and RM_Type_IIG each have posterior
probability 0.999 of being protective.

**The ordinal score** (script 23) adds no power and was not expected to: the
four nested thresholds reuse the same observations, so combining them as
independent studies would be double counting, and that combination is reported
as an effect size only. What it adds is a consistency check the binary analysis
cannot perform. **No system of 33 shows heterogeneity across thresholds.** The
associations are not artifacts of where the lysed/not-lysed boundary was drawn.

## 21. Searching the whole accessory genome

Every association in this study has been chosen from a catalogue. With 254
validated genomes the restriction can be lifted: script 25 tests every variable
k-mer in the collection against resistance, named or not.

From 63.2 million sampled k-mer observations, 616,300 are variable across
strains. Collapsing k-mers that occur in exactly the same strains, and keeping
only patterns backed by at least five k-mers, leaves **10,351 loci** to test.

**Two methodological failures are worth recording, because both were fixed.**

The first was sampling density. At one k-mer in fifty, 83 percent of presence
patterns were backed by a single k-mer: k-mers from the same gene landed in
slightly different sets of strains, so the collapsing step that is supposed to
turn a locus into one test did almost nothing. At one in twenty the median
pattern is backed by seven k-mers and the tests are loci rather than noise.

The second was the null. The phylogenetically matched null of script 21 fails
badly here: median p of 0.276 and 23.7 percent of loci below 0.05, about five
times the expected rate. The reason is specific. Simulated traits are built
from eight principal coordinates, and a real accessory locus confined to one
clade has sharper phylogenetic structure than eight axes can describe, so the
matched null is too weak. That null was adequate for 33 defense systems, whose
distributions are far less clade-restricted than a typical accessory k-mer.

The replacement is an empirical null: within each prevalence and phylogenetic
signal bin, the spread of statistics among loci that cannot nearly all be
causal is itself the null, estimated by median and scaled median absolute
deviation. It calibrates correctly and slightly conservatively (median p 0.500,
1.2 and 3.9 percent below 0.05 for the two outcomes).

**The result is a clean negative.** **No locus in the accessory genome of these
254 strains is associated with resistance.** Two are associated with the
opposite, being lysed more often, each present in 10 strains and neither
lineage-confined.

Script 26 then asks what they are. Neither carries an anti-defense or lysogeny
domain. But the carriers of one of them share **SanaTA and Septu**, two of the
systems that have shown an association with increased lysis throughout this
analysis and that no defense mechanism can explain.

That convergence is the most informative thing in this section. A scan that
used no catalogue at all landed on the same strains that the catalogued
analysis kept flagging as anomalous. The anomaly is therefore not an artifact
of how defense systems were annotated. The likeliest reading is that a mobile
element carries both the locus and those systems and marks a group of strains
that are easy to lyse for some other reason, which is consistent with the
finding in section 18 that recombinase content and defense content covary.

**What the scan cannot see** is as important as what it found. It cannot see a
locus confined to one lineage, because the null removes exactly that signal and
no observational design can separate a lineage-confined cause from the lineage.
It cannot see point mutations, since this is a presence test. It covers 254 of
402 strains, so an effect confined to a handful is invisible. And the outcome
is a strain's average across 96 phages, which dilutes anything acting on a
single phage lineage.

## Summary

Surface receptor variation is real but does not dominate: it explains about a
third of failed infections on the best available non-circular definition and,
on a fair order-independent partition, about three times the variation that all
thirty-three defense systems account for. Conditioning on it modestly widens
defense associations. Almost all of those associations are nonetheless
properties of the bacterial lineage: two of thirty-three survive control for
clonal background, and the biologically impossible associations disappear
entirely, which is the strongest evidence that the control is doing its job.
Three quarters of interactions still fail after accounting for everything
measured here.

Every attempt to strengthen the case for a specific system has instead
qualified it. Annotating the phage genomes directly left the counter-defense
conclusion unchanged. Rebuilding the receptor gate from actual receptor binding
protein homology did not improve it, which localises that weakness to the
single-isolation-host anchor rather than to the grouping. Extending from 33
systems to 79 left a pooled odds ratio of 0.923. None of 135 testable system
pairs interacts. A phylogenetically matched null removes CRISPR-Cas subtype I-E,
the most lineage-confined system in the collection, and shows that conditioning
on ancestry alone is not enough because host surface composition confounds
separately. Replacing empirical Bayes with full Bayesian inference changed no
conclusion. Using the whole plaque score rather than a binary outcome added no
power but confirmed that no association depends on where the threshold was
drawn.

Two things were established rather than qualified. The first is negative and
about design: resident prophages closely related to a collection phage exist in
only 5 of 254 strains and involve only the 2 temperate phages, which plaque on
0.98 percent of strains, so this collection cannot test superinfection exclusion
at all. The second is an internal check that passed: the phages found integrated
in host chromosomes are exactly the two independently classified temperate from
their own proteomes.

An unbiased scan of the accessory genome found nothing associated with
resistance either, which bounds the search beyond the catalogue: whatever
accounts for the unexplained failures, it is not a single widely distributed
accessory locus that this design could see.

The defensible output is a ranked candidate list. Eleven systems are nominated
somewhere in the study; the best, RM_Type_IIG, passes five of seven independent
checks and still fails the one that holds the entire genomic background fixed.
RM_Type_III passes four with twelve informative clonal clusters to RM_Type_IIG's
seven. No system passes all seven. That list is a priority order for knockout
experiments, not a set of findings.
