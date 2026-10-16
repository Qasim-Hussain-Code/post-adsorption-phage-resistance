"""Run the whole pipeline in order and stop at the first failure.

Usage:
    python scripts/run_all.py [--skip-fetch]
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

STEPS = [
    ("01_fetch_and_validate_data.py", "retrieve and verify the public inputs"),
    ("02_stratify_receptor_haplotypes.py", "harmonise data and build the receptor gate"),
    ("03_audit_defense_systems.py", "defense associations along the conditioning ladder"),
    ("04_control_population_structure.py", "control for clonal background"),
    ("05_evaluate_phage_counterdefense.py", "counter-defense from genus-level priors"),
    ("06_generate_figures.py", "main figures"),
    ("07_annotate_phage_antidefense.py", "annotate the 96 phage genomes directly"),
    ("08_counterdefense_direct.py", "counter-defense from the direct annotation"),
    ("09_robustness_and_sensitivity.py", "stress the surviving associations"),
    ("10_generate_extended_figures.py", "extended figures"),
    ("11_variance_partitioning.py", "order-independent partition of explained variation"),
    ("12_rbp_receptor_gate.py", "receptor gate from binding protein homology"),
    ("13_hierarchical_shrinkage.py", "all systems estimated jointly"),
    ("14_defense_combinations.py", "pairwise and class-level combinations"),
    ("15_consolidated_candidates.py", "one ranked candidate list"),
    ("16_generate_summary_figure.py", "summary figure"),
    ("17_phage_lifestyle.py", "temperate against virulent phages"),
    ("18_fetch_host_genomes.py", "resolve and download host genomes"),
    ("19_validate_host_genomes.py", "check each genome is the strain claimed"),
    ("20_prophage_analysis.py", "resident prophages and superinfection exclusion"),
    ("21_phylogenetic_null.py", "null traits matched on phylogeny and surface"),
    ("22_full_bayesian_pooling.py", "posterior intervals without the shortcut"),
    ("23_ordinal_thresholds.py", "the whole plaque score, not just lysis"),
    ("24_generate_inference_figure.py", "inference figure"),
    ("25_accessory_genome_gwas.py", "search the whole accessory genome"),
    ("26_identify_gwas_loci.py", "identify what the flagged loci are"),
    ("27_generate_calibration_figure.py", "association calibration figure"),
    ("28_phage_resolved_effects.py", "resolve the panel into its phages"),
    ("29_generate_spectrum_figure.py", "phage spectrum figure"),
    ("30_generate_workflow_diagram.py", "workflow diagram"),
]

# Steps 07 onward download phage genomes and a model archive and are much
# slower than the core pipeline, so they can be skipped independently.
EXTENDED_PREFIXES = ("07_", "08_", "09_", "10_", "11_", "12_", "13_", "14_",
                     "15_", "16_", "17_", "18_", "19_", "20_", "21_", "22_",
                     "23_", "24_", "25_", "26_", "27_", "28_", "29_", "30_")

# Steps 18 to 20 download several hundred bacterial genomes from NCBI, which
# throttles long runs. They are separable so that a network problem cannot
# stop the rest, and step 18 is resumable: completed downloads are skipped.
GENOME_PREFIXES = ("18_", "19_", "20_", "25_", "26_", "27_")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--skip-fetch", action="store_true",
                    help="assume the raw inputs are already present and verified")
    ap.add_argument("--core-only", action="store_true",
                    help="run only steps 01 to 06, skipping the phage annotation")
    ap.add_argument("--skip-genomes", action="store_true",
                    help="skip steps 18 to 20, which download host genomes")
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    started = time.time()
    for script, description in STEPS:
        if args.skip_fetch and script.startswith("01_"):
            print(f"[skip] {script}")
            continue
        if args.core_only and script.startswith(EXTENDED_PREFIXES):
            print(f"[skip] {script}  ({description})")
            continue
        if args.skip_genomes and script.startswith(GENOME_PREFIXES):
            print(f"[skip] {script}  ({description})")
            continue
        print(f"\n[run ] {script}  ({description})")
        t0 = time.time()
        result = subprocess.run([sys.executable, str(here / script)])
        if result.returncode != 0:
            print(f"[fail] {script} exited with {result.returncode}")
            return result.returncode
        print(f"[done] {script} in {time.time() - t0:.1f} s")
    print(f"\npipeline completed in {time.time() - started:.1f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
