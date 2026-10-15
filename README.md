# Post-adsorption phage resistance in Escherichia

A re-analysis of the Picard and Guelin phage-host interaction collection asking
a narrow question: once you account for which phage is being tested, for the
host surface it has to bind, and for the clonal lineage the host belongs to,
how many intracellular defense systems are still associated with resisting
infection?

The short answer is two, out of thirty-three systems with enough carriers to
test. Most of what looks like defense-mediated resistance in this dataset is a
property of the bacterial lineage rather than of the defense system.

Stress-testing those two against designs that do not rely on covariate
adjustment qualifies even that. Neither survives every test, and a third system
that the main analysis discarded holds up better in the strictest ones. The
defensible output of this work is a short candidate list for experiment, not a
demonstration that any system confers resistance.

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

### The candidate list

Nine systems are nominated somewhere in the study. Putting all of them through
the same battery:

| System | Strains | Adjusted OR | Within-cluster OR | Clusters | Genera below 1 | E-value | Checks |
|---|---|---|---|---|---|---|---|
| RM_Type_IIG | 52 | 0.549 | 0.968 | 7 | 15/15 | 2.69 | 3/5 |
| RM_Type_III | 54 | 0.626 | 0.328 | 12 | 10/18 | 2.31 | 3/5 |
| Shango | 37 | 0.583 | 0.492 | 9 | 16/17 | 2.51 | 2/5 |
| Lamassu-Cap4_nuclease | 18 | 0.451 | 0.417 | 5 | 12/13 | 3.34 | 2/5 |
| Zorya_TypeII | 20 | 0.444 | 0.344 | 2 | 11/13 | 3.40 | 2/5 |
| AbiD | 14 | 0.270 | 0.349 | 6 | 11/15 | 5.75 | 2/5 |

**No system passes all five checks**, and the two the main analysis singled out
are not clearly its best entries. RM_Type_III matches RM_Type_IIG while having
twelve informative clonal clusters to its seven; AbiD has the highest E-value in
the table on only fourteen carriers.

**Most resistance in this collection remains unexplained.** Among pairs passing
the surface gate whose host carries neither lineage-validated system, 75.8
percent of infections still fail, against 79.25 percent across the whole matrix.
Applying the surface gate and removing validated defense barely moves the
failure rate. A cleaner comparison is impossible: every one of the 402 strains
carries at least one of the 33 testable systems, so there is no defense-free
group in this collection.

## Figures

| Figure | Content |
|---|---|
| `figures/fig1_two_gate_filtering` | Outcome of all 38,435 pairs split into lysis, post-adsorption failure and surface exclusion, overall and by phage genus and host phylogroup |
| `figures/fig2_global_vs_conditioned_forest` | Paired forest plot of log2 odds ratios, unstratified against receptor-conditioned |
| `figures/fig3_lineage_controlled_survival` | Effect sizes before and after clonal control, and the attenuation toward the null |
| `figures/fig4_quantitative_titer_shifts` | Plaque score distributions by defense burden and phage morphotype, and per-system score shifts |
| `figures/fig5_direct_antidefense_annotation` | Anti-defense families found in the 96 genomes, validation against the RefSeq survey, and rescue interactions measured against imputed |
| `figures/fig6_robustness_of_surviving_effects` | Each estimate beside what it could have detected, threshold sensitivity, and per-genus consistency |
| `figures/fig7_consolidated_evidence` | Shapley partition of explained variation, 79 systems shrunk jointly, and the consolidated candidate list |

Each is written as 300 dpi PNG and vector PDF. Every figure is drawn from the
tables in `results/tables`, so a figure cannot disagree with a reported number.

## Reproducing

Requires Python 3.13 and the pinned versions in `requirements.txt`.

```bash
pip install -r requirements.txt
python scripts/run_all.py              # everything, 12 to 25 minutes
python scripts/run_all.py --core-only  # steps 01 to 06 only, about 2 minutes
```

Individual stages can be run in order instead; each writes its own report to
`results/` and its tables to `results/tables/`.

Steps 01 to 06 are the core analysis and take about two minutes, most of it in
the regression fitting of step 4. Steps 07 to 16 additionally download the 96
phage genomes (6.9 MB) and the DefenseFinder model release (49 MB, from which
only the anti-defense profiles are kept and the archive is then deleted). They
take 12 to 25 minutes depending on machine load, dominated by the permutation
tests in step 9, the sixteen nested fits in step 11 and the 135 pair fits in
step 14. Peak memory
stays below 1 GB. The repository plus all downloaded data occupies about 56 MB.

## Layout

```text
scripts/01 .. 06   core pipeline: fetch, stratify, audit, lineage control, figures
scripts/07 .. 10   direct phage annotation, counter-defense, robustness, figures
scripts/11 .. 17   variance partition, receptor gate, shrinkage, combinations,
                   consolidated candidates, summary figure, phage lifestyle
scripts/run_all.py runs them in order, --core-only stops after 06
src/papr/          shared estimators, paths, configuration, figure style
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

## Note on tooling

The analysis code in this repository was developed with AI assistance. All
results were produced by executing the committed scripts against the public
data, and every number quoted above is reproducible by running the pipeline.
