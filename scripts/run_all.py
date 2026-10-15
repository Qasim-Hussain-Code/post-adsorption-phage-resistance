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
]

# Steps 07 onward download phage genomes and a model archive and are much
# slower than the core pipeline, so they can be skipped independently.
EXTENDED_PREFIXES = ("07_", "08_", "09_", "10_", "11_", "12_", "13_", "14_",
                     "15_", "16_")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--skip-fetch", action="store_true",
                    help="assume the raw inputs are already present and verified")
    ap.add_argument("--core-only", action="store_true",
                    help="run only steps 01 to 06, skipping the phage annotation")
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
