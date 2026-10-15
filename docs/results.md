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

## Summary

Surface receptor variation is real but does not dominate: it explains about a
third of failed infections on the best available non-circular definition and
contributes less to the variance in lysis than defense carriage does.
Conditioning on it modestly widens defense associations. Almost all of those
associations are nonetheless properties of the bacterial lineage: two of
thirty-three survive control for clonal background, and the biologically
impossible associations disappear entirely, which is the strongest evidence that
the control is doing its job. Three quarters of interactions still fail after
accounting for everything measured here.

Stress-testing qualifies the headline rather than overturning it. The count of
two is stable across every clonal threshold, and both systems are protective in
direction across nearly every phage genus, but neither survives all four
stricter designs: RM_Type_IIG loses its effect among near-identical genomes
where the test had reasonable power, and Zorya_TypeII is simply untestable
there. Annotating the phage genomes directly, which was the largest outstanding
methodological gap, left every counter-defense conclusion intact. The honest
position is that this dataset supports a short list of candidates worth testing
experimentally, not a demonstration that any particular system confers
resistance.
