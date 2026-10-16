# Post-adsorption phage resistance in Escherichia

A re-analysis of the Picard and Guelin phage-host interaction collection asking
a narrow question: once you account for which phage is being tested, for the
host surface it has to bind, and for the clonal lineage the host belongs to,
how many intracellular defense systems are still associated with resisting
infection?

The short answer is two, out of thirty-three systems with enough carriers to
test. Most of what looks like defense-mediated resistance in this dataset is a
property of the bacterial lineage rather than of the defense system.

## The question and why it needs care

Roughly four fifths of strain-phage pairs in a species-wide screen fail to
produce plaques. A common reading is that most of those failures happen at the
cell surface, because the phage cannot adsorb, and that this surface variation
masks the contribution of intracellular defense systems. The implied remedy is
to restrict the analysis to receptor-compatible pairs and look for defense
effects there.

Two things complicate that plan in this dataset, and both are reported here
rather than worked around.

First, the released data contain no table stating which surface receptor each
phage requires. Receptor compatibility has to be inferred, and the inference is
weak (see below). A hard receptor-compatible subset presented as though it were
validated would overstate what the data support.

Second, defense systems are inherited vertically and co-segregate with clonal
lineages. Any association between carrying a system and resisting a phage is
confounded with everything else that lineage carries. The scale of that problem
is visible in the data itself: before lineage adjustment, fourteen systems are
significantly associated with *more* lysis, which no defense mechanism can
cause.

## Data

Open data from Gaborieau et al., Nature Microbiology, 2024, retrieved from
<https://github.com/mdmparis/coli_phage_interactions_2023>, with anti-defense
annotations from <https://github.com/mdmparis/antidefensefinder_2024>. Nothing
is redistributed here. `scripts/01_fetch_and_validate_data.py` downloads the
inputs (about 21 MB) and verifies them against `data/raw_manifest.sha256`.

Three details differ from the figures a reader might expect, and the pipeline
asserts the observed values rather than the nominal ones:

| Expected | Observed in the release |
|---|---|
| 403 strains by 96 phages, 38,688 pairs | 402 by 96, **38,592 cells**. One strain (`LF110`) appears in the metadata but has no row in the matrix |
| every pair measured | **157 cells are blank**, so **38,435 interactions** were actually observed |
| plaque score on a 0 to 6 scale | the released scale runs **0 to 4** |

Of the 38,435 observed interactions, 7,976 (20.75 percent) produced plaques and
30,459 (79.25 percent) did not.

The summary matrix was checked against the 318,816 replicate-level spot tests
that accompany it. Treating a dilution step as positive when plaques appeared in
at least two readable replicates reproduces the **binary lytic call for 99.69
percent** of comparable cells, with every released zero recovered as a zero, and
correlates with the ordinal score at Spearman 0.991. The phenotype used
throughout is therefore well supported by the underlying measurements.

## What the analysis does

**Receptor gating, and why it is only a secondary analysis.** Each phage has a
known isolation host, and a strain that yielded a phage must present a surface
that phage can bind. Pooling isolation hosts across a VIRIDIC genus, whose
members share receptor binding proteins, gives a permissive surface state set
per genus that uses no information from the interaction outcomes. That gate is
honest but weak: it admits 65.7 percent of pairs and lifts the lysis rate only
from 17.0 to 22.7 percent, an enrichment of 1.34 fold. It assigns **35.9
percent** of failed infections to surface exclusion, not the 80 to 85 percent a
strict surface-masking account predicts. Adding O-antigen serotype cuts coverage
to 12.3 percent of pairs while enrichment stays near 1.4 fold, which is what
happens when a filter removes signal and noise together. Enrichment is strong in
some genera (Felixounavirus 3.2 fold) and absent in the two largest
(Vectrevirus 1.05 fold, Kayfunavirus 1.20 fold).

**So the primary analysis stratifies on surface state instead of filtering on
it.** Comparing carriers with non-carriers inside the same phage and the same
LPS outer core type removes the surface contribution without discarding data and
without needing to know which state is permissive. Effects are estimated along a
conditioning ladder, from unstratified through to stratified on phage identity
and surface type, then refitted with standard errors clustered on clonal
background and covariates for phylogroup and eight core-genome principal
coordinates.

## Results

**Surface conditioning widens defense effects rather than removing them.** The
median absolute log2 odds ratio across 33 systems rises monotonically with
conditioning depth: 0.282 unstratified, 0.343 stratified by phage, 0.469 adding
LPS outer core, 0.517 adding O-antigen status. This is the direction the masking
hypothesis predicts, but the magnitude is modest. Surface variation masks these
associations partially, not dominantly.

**Which phage is being tested matters far more than anything about the host,
and host surface outweighs defense.** Allocating explained variation by exact
Shapley value over all orderings of the blocks, with a surface block containing
everything the collection types about the cell surface:

| Block | Terms | Shapley pseudo R2 | Share |
|---|---|---|---|
| phage identity | 95 | 0.1920 | 59.7 percent |
| host surface | 157 | 0.0830 | 25.8 percent |
| defense systems | 33 | 0.0269 | 8.4 percent |
| clonal lineage | 22 | 0.0199 | 6.2 percent |

Host surface accounts for about **three times** what all 33 defense systems do,
and that ordering holds whether each block is fitted alone, added last, or
allocated by Shapley value. An earlier version of this analysis reported the
opposite; that comparison gave surface six terms against defense's thirty-three
and depended on the order the blocks were entered. Per term the two are closer,
with defense marginally ahead, since surface spends 157 terms to defense's 33.

**Almost nothing survives clonal control.** Fitted as a regression with
clustered standard errors, 6 of 33 systems are significant before lineage
adjustment. After adjusting for phylogroup and core-genome position, **2
remain**:

| System | Odds ratio | 95 percent CI | q |
|---|---|---|---|
| Zorya_TypeII | 0.444 | 0.282 to 0.700 | 0.0077 |
| RM_Type_IIG | 0.549 | 0.398 to 0.757 | 0.0077 |

The adjustment passes an internal check. All fourteen associations with
increased lysis, which are biologically impossible as defense effects,
disappear once lineage is controlled. The covariates remove the kind of signal
they were introduced to remove, which makes the attenuation of the protective
associations harder to dismiss as overcorrection.

**Defense burden is flatly null.** The number of systems a strain carries has no
monotone relationship with how often it is lysed (lysis rate 19.4, 22.3, 21.4
and 19.9 percent across burden quartiles; strain-level Spearman rho +0.053,
p = 0.29). The pair-level correlation is nominally significant only because 402
strains contribute 38,435 correlated observations.

**Counter-defense: one suggestive signal, not a result.** The 96 phage genomes
were annotated directly (genes called with Prodigal, searched against the 178
curated anti-defense profiles at their own gathering cutoffs). 73 of 96 carry at
least one anti-defense system and 48 carry Anti_RM. Of 12 systems with a
published counterpart, one interaction is significant: RM_Type_IV protects less
well against phages carrying anti-RM genes (interaction odds ratio 1.54, 95
percent CI 1.17 to 2.02, q = 0.026), and all five RM subtypes point that way.
This rests on one significant test out of twelve. Anti-CRISPR, Anti-BREX,
Anti-Gabija and Anti-Hachiman are absent from all 96 genomes, so the systems
they counteract cannot be tested for rescue here at all.

Earlier versions of this analysis imputed anti-defense content from the phage's
genus, which attenuates toward the null and was the main reason to doubt the
negative result. Direct annotation gives the same answer (one significant
interaction either way, direction agreeing for 11 of 12 systems, estimates
correlating at Spearman +0.993), so the conclusion does not rest on imputation.

**The two survivors do not survive everything.** Four designs that make no use
of covariate adjustment give a more complicated picture:

| System | Adjusted OR | Within-cluster OR | Informative clusters | Perm q (phylogroup) | Genera with OR below 1 | E-value | Strict tests passed |
|---|---|---|---|---|---|---|---|
| RM_Type_IIG | 0.549 | 0.968 | 7 | 0.0033 | 15 of 15 | 2.69 | 2 of 4 |
| RM_Type_III | 0.626 | 0.328 | 12 | 0.30 | n/a | 2.31 | 2 of 4 |
| Zorya_TypeII | 0.444 | 0.344 | 2 | 0.30 | 11 of 13 | 3.40 | 1 of 4 |

RM_Type_IIG loses its effect among near-identical genomes, where the test had
reasonable power, which is evidence against rather than absence of evidence.
Zorya_TypeII keeps its point estimate there but only two clusters are
informative, so the test is uninformative. RM_Type_III, which the main analysis
discarded, is significant and protective in two of the stricter designs and is
arguably the best-supported candidate in the whole study. The count of two is
nonetheless stable across all nine clonal distance thresholds tried, and both
survivors are protective in direction across nearly every phage genus.

**Most nulls are uninformative rather than negative.** Across the 31 systems
reported as null, the protective odds ratio detectable at 80 percent power
ranges only from 0.43 to 0.70. For most of them the analysis had no power to
see a moderate effect, and they should not be read as evidence that the system
does nothing.

**All 79 systems, estimated jointly.** Fitting every system carried by at least
five strains and shrinking each toward the pooled mean extends coverage from 33
to 79. The pooled odds ratio is **0.923** (95 percent CI 0.847 to 1.005) with a
between-system standard deviation of 0.269, so systems genuinely differ from
one another but the average has essentially no association with lysis. Four
candidates below the old 20-carrier threshold surface here: PD-T7-4, AbiD,
Lamassu-Cap4_nuclease and BstA.

**A receptor gate built from the actual receptor is no better.** Receptor
binding protein sequences were recovered for 87 phages and used to group them by
homology (30 to 45 groups against 19 genera). The resulting gate reaches 1.26
fold enrichment against the genus gate's 1.34. Since the binding protein is the
right biological object, this localises the weakness to the **anchor**, one
isolation host per phage, which can show a surface state is permissive but never
that another is not.

**Combinations do not rescue the burden null.** Of 135 system pairs with enough
strains in all four carriage combinations, **none** shows a significant
interaction. Counting systems fails to predict lysis, and so does pairing them.

**Lysogeny is a real but marginal explanation for failure.** A temperate phage
entering lysogeny leaves no plaque, which a spot test cannot distinguish from
resistance. Searching the phage proteomes against Pfam lysogeny profiles finds a
tyrosine integrase in just **2 of 96 phages**, one Lambdavirus and one
Wanchaivirus, and those two are exactly the ones that plaque on under 1 percent
of strains (against 21.2 percent for the rest). The direction fits lysogeny; the
magnitude does not. They account for 2.6 percent of all failures, and removing
them shifts the overall failure rate from 79.25 to 78.83 percent and the median
defense odds ratio by 0.002 in log2.

**Resident prophages: this collection cannot test them.** The host genomes are
not in the released data, so 287 were resolved against NCBI individually and
validated against the released core-genome distance matrix; **254 passed, 63
percent of the collection**. Across 24,382 analysable pairs, exactly **8** have
a host carrying a prophage closely related to the challenging phage, involving 5
strains and only 2 phages, and nearly doubling coverage added none. All 8 are
non-lytic, but those phages plaque on 0.98 percent of strains, so **0.08** lysis
events were expected under the null (odds ratio 0.13, interval 0.01 to 1.84). The
limitation is the panel: it was assembled for lytic activity, so only 2 of 96
phages are temperate. One check does hold: the two phages found integrated in
host chromosomes are exactly the two independently classified temperate from
their own proteomes.

**A phylogenetically matched null, and what it exposes.** Simulating null traits
with each system's own phylogenetic confinement measures the population
structure problem directly: median inflation 1.13, but **4.17 for
CAS_Class1-Subtype-I-E**, whose carriage is 59 percent explained by core-genome
position and which drops from p = 0.0005 to p = 0.50. Conditioning on ancestry
alone is not enough, though. Host surface type confounds separately, and
RM_Type_IIG turns out to be strongly depleted in permissive outer cores
(OR 0.33, p = 0.0004): its carriers are drawn from strains that are harder to
lyse whatever they carry. Removing surface type attenuates its effect by 31
percent; it remains significant.

**Proper posteriors and the whole score.** Full Bayesian pooling replaces the
empirical Bayes shortcut; tau's posterior (0.205 to 0.371) contains the earlier
point estimate and no conclusion changes. Using all four plaque-score
thresholds adds no power, because they reuse the same observations, but
confirms that **no system of 33** shows heterogeneity across thresholds.

**The whole accessory genome, searched without a catalogue.** With 254
validated genomes, script 25 tests every variable k-mer rather than a
pre-defined feature list: 616,300 variable k-mers collapsing into **10,351
loci**, every one of them tested. **No locus outside the lineage-confined blind
spot is associated with resistance.** Two are associated with
being lysed more often, and the carriers of one share **SanaTA and Septu**,
two systems that have shown that same impossible direction throughout. A scan
using no catalogue landed on the same strains the catalogued analysis kept
flagging, so the anomaly is not an artifact of annotation.

Getting there required fixing three failures, all documented: sampling k-mers
at one in fifty left 83 percent of patterns backed by a single k-mer; the
phylogenetically matched null is still 2.8 times too liberal here, even
after correcting its generator, because eight
principal coordinates cannot describe a clade-confined accessory locus; and the
binned empirical null that replaced it was well calibrated but defined for only
3,144 loci. Matching each locus to its nearest neighbours and reading it
against the tail actually present there covers all 10,351. That last fix
explained the earlier restriction rather than merely lifting it: the loci the
fixed bins dropped are sharply non-normal (fitted degrees of freedom around 8,
against 50 for the loci they kept), so the bin threshold had been selecting the
loci where a normal reference holds. The three loci the full scan does flag are
lineage-confined clade markers whose neighbourhood spread has collapsed, and
the report says so with the numbers rather than presenting them as findings.

**The panel resolved into its phages.** Every association above is an average
over 96 phages, and defense systems are narrow in spectrum, so script 28 tests
each system against each phage separately. Two corrections decide the answer:
near-identical phages are collapsed into 54 groups, since a system associated
with one member gets the rest for free, and a second null is built because the
matched null of section 18 is not clonal enough (real systems agree across
same-cluster strain pairs 0.911 of the time against 0.810 for that null, with
25 of 33 systems exceeding it). The second null permutes carrier counts among
clonal clusters of the same size and reproduces 0.911 exactly. Reporting the
worse of the two p-values is the difference between six findings and one.

**The direction that should not exist comes close, and does not arrive.**
Applying the same two nulls to the opposite tail, five systems cluster just
above the line and **none clears the correction**: BREX_I (q = 0.053),
CAS_Class1-Subtype-I-F and MazEF (0.059), Gao_Qat (0.073), SanaTA (0.086). The
nearest was re-tested on its own against 50,000 null traits of each kind and
gives p = 0.0018 against the 0.0015 required, so the verdict is the data's and
not the generator's. An earlier version of this analysis reported three of them
as clearing while Retron_I_C failed; both halves were artifacts of running the
null at a thousand draws, where the smallest attainable p sits almost exactly
on the correction threshold.

**One system clears both nulls, and then fails a stability check.** Retron_I_C
reduces lysis by **0.222** on the 12 phages it acts on against **0.016** across
the whole panel, a fourteen-fold gap that explains why every other test here
missed it, and retrons act by abortive infection, which is post-adsorption and
exactly what this project set out to find. It clears at q = 0.015 and keeps its
signal when its island neighbours are conditioned on. But dropping each clonal
cluster of carriers in turn, one cluster holding 6 of its 20 carriers takes the
result to **p = 0.160**, while 40 random drops of the same size leave a median
p of **0.009** and only 8 percent do as much damage. Six near-identical strains
are carrying it. And conditioning on host surface type, which the per-phage
model does not carry, takes it from 11 phage groups to **4** and its p-value to
**0.016**: 55 percent of its carriers sit on a K12 outer core against 8 percent
of other strains, and K12 strains are lysed 0.138 of the time against 0.262 for
the common R1 core. **That makes it the study's best candidate and not a
finding.**

### The candidate list

Nine systems are nominated somewhere in the study. Putting all of them through
the same battery:

Eleven systems are nominated somewhere in the study. All are put through the
same seven checks: lineage-adjusted regression, within clonal cluster, two
permutation schemes, shrinkage, a surface-adjusted phylogenetic null, and the
Bayesian posterior.

| System | Strains | Adjusted OR | Within-cluster OR | Clusters | Genera below 1 | E-value | P(protective) | Checks |
|---|---|---|---|---|---|---|---|---|
| RM_Type_IIG | 52 | 0.549 | 0.968 | 7 | 15/15 | 2.69 | 0.999 | 5/7 |
| Zorya_TypeII | 20 | 0.444 | 0.344 | 2 | 11/13 | 3.40 | 0.998 | 4/7 |
| RM_Type_III | 54 | 0.626 | 0.328 | 12 | 10/18 | 2.31 | 0.981 | 4/7 |
| Shango | 37 | 0.583 | 0.492 | 9 | 16/17 | 2.51 | 0.978 | 3/7 |
| Lamassu-Cap4_nuclease | 18 | 0.451 | 0.417 | 5 | 12/13 | 3.34 | 0.994 | 3/7 |
| AbiD | 14 | 0.270 | 0.349 | 6 | 11/15 | 5.75 | 0.995 | 3/7 |

**No system passes all seven.** RM_Type_IIG leads at five but still fails the
one test that holds the entire genomic background fixed, where it had moderate
power. RM_Type_III passes four with twelve informative clonal clusters to
RM_Type_IIG's seven. AbiD has the highest E-value in the study on only fourteen
carriers.

**Most resistance in this collection remains unexplained.** Among pairs passing
the surface gate whose host carries neither lineage-validated system, 75.8
percent of infections still fail, against 79.25 percent across the whole matrix.
Applying the surface gate and removing validated defense barely moves the
failure rate. A cleaner comparison is impossible: every one of the 402 strains
carries at least one of the 33 testable systems, so there is no defense-free
group in this collection.

## Figures

Each is written as 300 dpi PNG and vector PDF. Every figure is drawn from the
tables in `results/tables`, so a figure cannot disagree with a reported number.

### Figure 1. Surface exclusion accounts for a minority of failed infections

![Outcome of all 38,435 strain-phage pairs split into productive lysis, failure
in a surface-compatible host, and surface exclusion, shown overall and broken
down by phage genus and host phylogroup](figures/fig1_two_gate_filtering.png)

### Figure 2. Conditioning on surface state widens defense effects

![Paired forest plot of log2 odds ratios for 33 defense systems, comparing the
unstratified estimate with one stratified on phage identity and LPS outer core
type](figures/fig2_global_vs_conditioned_forest.png)

### Figure 3. Almost every association is a property of the clonal lineage

![Effect sizes before and after controlling for phylogroup and core-genome
position, with a scatter showing attenuation toward the null](figures/fig3_lineage_controlled_survival.png)

### Figure 4. Plaque score tracks phage morphotype, not defense burden

![Plaque score distributions by defense burden quartile and by phage
morphotype, with per-system shifts in mean score](figures/fig4_quantitative_titer_shifts.png)

### Figure 5. Annotating the phages directly leaves the answer unchanged

![Anti-defense families found in the 96 phage genomes, validation against the
published RefSeq survey, and rescue interactions from measured annotation
against the genus prior](figures/fig5_direct_antidefense_annotation.png)

### Figure 6. Neither surviving association is robust to every design

![Each system's estimate beside the protective odds ratio its own standard
error could have detected, sensitivity to the clonal distance threshold, and
per-genus consistency](figures/fig6_robustness_of_surviving_effects.png)

### Figure 7. Host surface outweighs defense, and no candidate survives every check

![Shapley partition of explained variation, 79 systems shrunk jointly toward
the pooled mean, and the consolidated candidate list](figures/fig7_consolidated_evidence.png)

### Figure 8. Harder nulls and a candidate list no single test supports

![The population structure problem measured per system, the phylogenetically
matched null before and after removing host surface type, and the final
verdict](figures/fig8_inference_and_verdict.png)

### Figure 9. Why the accessory genome scan needed its null rebuilt twice

![Four nulls against the uniform diagonal, their behaviour in the tail that
FDR control depends on, the heavier tails among loci the fixed bins dropped,
and the collapsed neighbourhood that disqualifies every locus the full scan
flags](figures/fig9_association_calibration.png)

### Figure 10. What averaging over a panel of 96 phages costs, and which direction survives

![The panel-average effect against the effect on the phages each system acts
on, the one protective candidate across the whole panel with the clonal cluster
that carries it, the clonal structure the coordinate null fails to reproduce,
what each null admits, the three systems clearing both nulls in the
increased-lysis direction against the one in the protective direction, and the
leave-one-cluster-out check that separates
them](figures/fig10_phage_resolved_spectrum.png)

## Reproducing

Requires Python 3.13 and the pinned versions in `requirements.txt`.

```bash
pip install -r requirements.txt
python scripts/run_all.py                 # everything
python scripts/run_all.py --skip-genomes  # omit the steps needing host genomes
python scripts/run_all.py --core-only     # steps 01 to 06 only
```

`--skip-genomes` omits steps 18 to 20, which download several hundred bacterial
genomes from NCBI, together with steps 25 to 27, which need those genomes. The
phage-resolved analysis of steps 28 and 29 does not, so it still runs.

Individual stages can be run in order instead; each writes its own report to
`results/` and its tables to `results/tables/`.

```bash
python tests/test_statistics.py   # the estimators, against independent versions
```

That compares every estimator in `src/papr` against statsmodels, scipy or a
brute-force calculation and exits non-zero if any disagrees. It is worth
running before trusting anything else here: a defect in one of those functions
would be invisible to a check that compares the written results against the
computed ones, since the documentation would faithfully report the wrong
number. That is how the one defect found this way came to light, a factor of
one half in the Mantel-Haenszel variance that made its intervals about a tenth
too narrow.

The timings below are measured, from one `--skip-fetch --skip-genomes` run on
the 16 GB machine this was developed on, which completed every step in
**6906 s**. Steps 02 to 06, the core analysis, take **7 minutes**, almost all of
it the cluster-robust regression fitting of step 04 (337 s). Steps 07 to 17 add
the 96 phage genomes (6.9 MB) and the DefenseFinder model release (49 MB, from
which only the anti-defense profiles are kept before the archive is deleted),
and take **61 minutes**, dominated by the permutation tests of step 09 (1577 s)
and the sixteen nested fits of step 11 (1013 s). Steps 21 to 24 add **47
minutes**, almost entirely the four-threshold ordinal fits of step 23 (2782 s).
Step 28 postdates that run and was timed separately at roughly 20 minutes, for
its 66,000 scored traits and the refits behind the stability checks. Peak
memory stays below 1 GB.

The repository plus all downloaded data occupies about 642 MB: 414 MB of
bacterial genomes, 100 MB of derived k-mer sketches (regenerable and
untracked), 72 MB of git history and the rest the public input tables. Run with
`--skip-genomes` it stays near 110 MB, and the genomes can be deleted after
step 20 without affecting any other result.

Every output is deterministic. Rerunning the whole pipeline reproduces all 121
tables, reports and figures byte for byte, which is checked by regenerating
them and comparing checksums rather than assumed. One caveat is worth stating
precisely: the coefficients that come from iterative model fits are written at
full float precision, and their last digits depend on the order the linear
algebra library accumulates its sums. They are stable on one machine, and one
such coefficient has been seen to differ in the eighth significant figure
between environments. Nothing rounds differently at any precision the analysis
reports.

## Layout

```text
scripts/01 .. 06   core pipeline: fetch, stratify, audit, lineage control, figures
scripts/07 .. 10   direct phage annotation, counter-defense, robustness, figures
scripts/11 .. 17   variance partition, receptor gate, shrinkage, combinations,
                   consolidated candidates, summary figure, phage lifestyle
scripts/18 .. 20   host genome retrieval, identity validation, prophages
scripts/21 .. 24   phylogenetic null, Bayesian pooling, ordinal score, figure
scripts/25 .. 27   accessory genome scan, locus identification, calibration figure
scripts/28 .. 29   phage-resolved defense effects, spectrum figure
scripts/run_all.py runs them in order; --core-only stops after 06 and
                   --skip-genomes omits 18 to 20 and 25 to 27
src/papr/          shared estimators, paths, configuration, figure style
tests/             the statistical core checked against statsmodels and scipy
data/raw/          downloaded inputs, not tracked; verified against data/raw_manifest.sha256
data/processed/    derived tables
results/           per-stage reports in plain text
results/tables/    every number behind the figures, as CSV
figures/           manuscript figures, PNG and PDF
docs/              methods, results, limitations, data provenance
```

## Limitations

The main ones, stated in full in `docs/limitations.md`:

- Receptor compatibility is inferred, not measured, and the inference is weak.
- Neither surviving association passes all four stricter designs, and most of
  the 31 nulls are uninformative rather than negative.
- Resident prophages were not annotated, so superinfection exclusion remains an
  untested explanation for part of the residual resistance. This is now the
  largest open gap.
- Anti_CRISPR, Anti_BREX, Anti_Gabija and Anti_Hachiman are absent from all 96
  phages, so CRISPR-Cas, BREX, Gabija and Hachiman cannot be tested for rescue
  in this collection.
- Mantel-Haenszel intervals in step 3 treat correlated pairs as independent and
  are too narrow; the clustered estimates in steps 4 and 9 are the inferential
  ones.
- This is an observational analysis of one collection. Nothing here establishes
  that any system causes resistance, only that some associations are not
  explained by the covariates available.

The anti-defense imputation that earlier versions relied on is no longer a
limitation: the genomes are now annotated directly and the conclusion is
unchanged.

## License

MIT, see `LICENSE`. The upstream data remains under the license of its
originating repository and is not redistributed; see `docs/data_provenance.md`.
