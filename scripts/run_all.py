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
    ("05_evaluate_phage_counterdefense.py", "phage counter-defense and residual resistance"),
    ("06_generate_figures.py", "manuscript figures"),
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--skip-fetch", action="store_true",
                    help="assume the raw inputs are already present and verified")
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    started = time.time()
    for script, description in STEPS:
        if args.skip_fetch and script.startswith("01_"):
            print(f"[skip] {script}")
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
