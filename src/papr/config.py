"""Analysis constants, declared once so that every script shares them.

Thresholds were fixed before the association tests were run and are not tuned
to the results. Where a threshold is a judgement call the reasoning is given.
"""

from __future__ import annotations

# A defense system needs enough carrier strains for a stratified contingency
# test to be informative. Twenty carriers across 402 strains is the point below
# which most phage strata contain no carrier at all and the Mantel-Haenszel
# estimate rests on a handful of tables.
MIN_CARRIER_STRAINS = 20

# False discovery rate for the Benjamini-Hochberg correction applied across
# defense systems within each family of tests.
FDR_ALPHA = 0.05

# Bacterial surface determinants, in the order they gate adsorption. The LPS
# outer core is complete for all tested strains; the capsule columns are absent
# for the 33 strains contributed as phage isolation hosts, which is why capsule
# state is used only in a sensitivity analysis.
LPS_COLUMN = "LPS_type"
O_TYPE_COLUMN = "O-type"
CAPSULE_COLUMNS = [
    "Capsule_ABC",
    "Capsule_GroupIV_e",
    "Capsule_GroupIV_s",
    "Capsule_Wzy_stricte",
]

# Outer membrane protein loci typed by BLAST at 99 percent identity. Several
# are established phage receptors in Escherichia coli (LamB, OmpC, OmpA, FhuA,
# BtuB, TolC, Tsx, FadL, NfrA).
OMP_LOCI = ["BTUB", "FADL", "FHUA", "LAMB", "LPTD", "NFRA",
            "OMPA", "OMPC", "OMPF", "TOLC", "TSX", "YNCD"]

# Clonal background covariates.
PHYLOGROUP_COLUMN = "Clermont_Phylo"
ST_COLUMN = "ST_Warwick"

# Number of principal coordinates retained from the core-genome distance
# matrix. Eight axes carry the between-phylogroup structure; adding more mainly
# splits individual Sequence Types, which the clonal cluster term already
# absorbs.
N_PCOA_AXES = 8

# Patristic distance at which two strains are treated as the same clonal
# background when forming clusters for the robust standard errors.
CLONAL_DISTANCE_THRESHOLD = 0.0005

# Resampling settings.
RANDOM_SEED = 20260101
N_BOOTSTRAP = 2000
N_PERMUTATIONS = 2000

# Figure output.
FIGURE_DPI = 300
