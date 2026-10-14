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

**Which phage is being tested matters far more than anything about the host.**
In a nested logistic decomposition, phage identity alone reaches a McFadden
pseudo R2 of 0.176. Host surface type adds 0.031, all 33 defense systems
together add 0.043, and phylogroup adds 0.011. Defense carriage contributes
slightly more than the surface determinants typed in this collection, and
neither approaches the phage term.

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

**Counter-defense: one suggestive signal, not a result.** Anti-defense content
was imputed at the level of the phage genus from the published survey of RefSeq
phage genomes, so these are statements about lineages, not about which genes a
phage carries. Of 13 systems with a published anti-defense counterpart and a
usable split, one interaction is significant: RM_Type_IV protects less well
against lineages whose relatives commonly carry anti-RM genes (interaction odds
ratio 1.63, q = 0.015). Four of the five RM subtypes point the same way, which
is the pairing where a lineage-level prior should carry the most information,
but this rests on one significant test out of thirteen. Anti-CRISPR proteins are
absent from every adequately surveyed genus here, so CRISPR-Cas systems could
not be tested for rescue at all.

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

Each is written as 300 dpi PNG and vector PDF. Every figure is drawn from the
tables in `results/tables`, so a figure cannot disagree with a reported number.

## Reproducing

Requires Python 3.13 and the pinned versions in `requirements.txt`.

```
pip install -r requirements.txt
python scripts/run_all.py
```

Individual stages can be run in order instead; each writes its own report to
`results/` and its tables to `results/tables/`. The full pipeline takes a few
minutes, most of it in the regression fitting of step 4. Peak memory stays below
1 GB and the repository plus downloaded data stays near 25 MB.

## Layout

```
scripts/     numbered pipeline stages, run in order
src/papr/    shared estimators, paths, configuration, figure style
data/raw/    downloaded inputs, not tracked; verified against data/raw_manifest.sha256
data/processed/  derived tables
results/     per-stage reports in plain text
results/tables/  every number behind the figures, as CSV
figures/     manuscript figures, PNG and PDF
docs/        methods, results, limitations, data provenance
```

## Limitations

The main ones, stated in full in `docs/limitations.md`:

- Receptor compatibility is inferred, not measured, and the inference is weak.
- Anti-defense content is imputed per genus from sequenced relatives, which
  attenuates any real interaction toward the null.
- Resident prophages were not annotated, so superinfection exclusion remains an
  untested explanation for part of the residual resistance.
- Mantel-Haenszel intervals in step 3 treat correlated pairs as independent and
  are too narrow; the clustered estimates in step 4 are the inferential ones.
- This is an observational analysis of one collection. Nothing here establishes
  that Zorya_TypeII or RM_Type_IIG cause resistance, only that their association
  with it is not explained by the covariates available.

## License

MIT, see `LICENSE`. The upstream data remains under the license of its
originating repository and is not redistributed; see `docs/data_provenance.md`.

## Note on tooling

The analysis code in this repository was developed with AI assistance. All
results were produced by executing the committed scripts against the public
data, and every number quoted above is reproducible by running the pipeline.
