"""Canonical filesystem layout, resolved relative to the repository root.

Scripts call `bootstrap()` so that `import papr` works regardless of the
working directory the user launched them from.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

SRC = ROOT / "src"
SCRIPTS = ROOT / "scripts"
DATA = ROOT / "data"
RAW = DATA / "raw"
PROCESSED = DATA / "processed"
RESULTS = ROOT / "results"
TABLES = RESULTS / "tables"
FIGURES = ROOT / "figures"
DOCS = ROOT / "docs"

RAW_MANIFEST = DATA / "raw_manifest.sha256"

# Raw inputs, named once here so that no script hard-codes an upstream path.
MATRIX = RAW / "interactions" / "interaction_matrix.csv"
RAW_INTERACTIONS = RAW / "interactions" / "raw_interactions.csv"
PICARD = RAW / "bacteria" / "picard_collection.csv"
ECTYPER = RAW / "bacteria" / "ectyper_o_type.tsv"
LPS_WAAL = RAW / "bacteria" / "lps_type_waaL.txt"
OMP_ALLELES = RAW / "bacteria" / "omp_allele_clusters.tsv"
KAPTIVE = RAW / "bacteria" / "kaptive_capsule_results.txt"
DEFENSE = RAW / "bacteria" / "defense_systems_subtypes.csv"
CORE_TREE = RAW / "bacteria" / "core_genome_tree.nwk"
CORE_DISTANCES = RAW / "bacteria" / "core_genome_distance_matrix.tsv"
GUELIN = RAW / "phages" / "guelin_collection.csv"
RBP = RAW / "phages" / "rbp_list.csv"
PHAGE_TREE = RAW / "phages" / "phage_viridic_tree.nwk"
ADF_PHAGE = RAW / "antidefense" / "adf_hits_refseq_phage.tsv"
ADF_TAXONOMY = RAW / "antidefense" / "phage_taxonomy_refseq.tsv"

# Inputs for the direct annotation of the 96 phage genomes (script 07).
PHAGE_GENOMES = RAW / "phage_genomes"
ADF_MODELS = RAW / "adf_models"
ADF_PROFILES = ADF_MODELS / "profiles"
ADF_DEFINITIONS = ADF_MODELS / "definitions"

# Derived tables written by the pipeline.
PAIRS = PROCESSED / "pair_level_table.csv.gz"
STRAIN_FEATURES = PROCESSED / "strain_features.csv"
PHAGE_FEATURES = PROCESSED / "phage_features.csv"
GATE_DEFINITIONS = PROCESSED / "receptor_gate_definitions.csv"
PCOA = PROCESSED / "core_genome_pcoa.csv"
ADF_GENUS = PROCESSED / "antidefense_genus_prevalence.csv"
PHAGE_PROTEINS = PROCESSED / "phage_proteins.faa.gz"
PHAGE_ADF_HITS = PROCESSED / "phage_antidefense_hits.csv"
PHAGE_ADF_SYSTEMS = PROCESSED / "phage_antidefense_systems.csv"


def bootstrap() -> Path:
    """Put `src` on the import path and create the output directories."""
    if str(SRC) not in sys.path:
        sys.path.insert(0, str(SRC))
    for d in (RAW, PROCESSED, TABLES, FIGURES, DOCS):
        d.mkdir(parents=True, exist_ok=True)
    return ROOT
