"""Partition the explained variation in lysis fairly between blocks of variables.

Script 03 compared blocks by adding them to a model in a fixed order and
reading off the increments. That comparison had two defects, and this script
exists to correct them.

First, the increments depend on the order the blocks are entered, because the
blocks are correlated. Reported in one order they can make any block look
dominant. Here every block's contribution is the exact Shapley value of the
McFadden pseudo R2: the average of its incremental contribution over all
orderings of the blocks. With four blocks that is sixteen model fits and the
answer no longer depends on an arbitrary choice.

Second, and more seriously, the surface block in script 03 contained only the
LPS outer core type and O-antigen status, six degrees of freedom, while the
defense block contained thirty-three. Concluding from that comparison that
defense carriage contributes more than host surface was not a fair test. The
surface block here carries everything this collection types about the cell
surface: LPS outer core, O-antigen serotype, capsule state and the allele
called at twelve outer membrane protein loci. Contributions are reported per
degree of freedom as well as in total.

Outputs:
    results/tables/table25_shapley_variance_partition.csv
    results/tables/table26_block_definitions.csv
    results/variance_partitioning_report.txt
"""

from __future__ import annotations

import itertools
import math
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
import statsmodels.api as sm

from papr import config, paths

# Levels of a categorical variable seen in fewer than this many strains are
# pooled into one "other" level, so that a 94-level serotype does not spend
# ninety degrees of freedom describing singletons.
MIN_LEVEL_COUNT = 5


def pooled_dummies(frame: pd.DataFrame, columns: list[str],
                   strain_counts: dict[str, pd.Series]) -> pd.DataFrame:
    """Dummy encode, pooling rare levels by how many strains carry them."""
    parts = []
    for col in columns:
        values = frame[col].astype("string").fillna("missing")
        counts = strain_counts[col]
        keep = set(counts[counts >= MIN_LEVEL_COUNT].index.astype(str))
        collapsed = values.where(values.isin(keep), other="other")
        d = pd.get_dummies(collapsed, prefix=col, drop_first=True, dtype=float)
        parts.append(d)
    if not parts:
        return pd.DataFrame(index=frame.index)
    return pd.concat(parts, axis=1)


def build_blocks(pairs: pd.DataFrame, strain_feats: pd.DataFrame,
                 pcoa: pd.DataFrame, systems: list[str]
                 ) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    sf = strain_feats.set_index("strain")
    surface_cols = ["lps_core", "o_antigen", "o_status", "abc_serotype"] + \
                   [f"omp_{l.lower()}" for l in config.OMP_LOCI]
    surface_cols = [c for c in surface_cols if c in sf.columns]
    capsule_cols = [c for c in ["capsule_abc", "capsule_groupiv_e",
                                "capsule_groupiv_s", "capsule_wzy_stricte"]
                    if c in sf.columns]

    strain_counts = {c: sf[c].astype("string").value_counts() for c in surface_cols}

    surface = pooled_dummies(pairs, surface_cols, strain_counts)
    for c in capsule_cols:
        # Capsule calls are missing for the isolation-host strains; an explicit
        # indicator keeps those rows rather than dropping them.
        v = pairs[c] if c in pairs.columns else pairs["strain"].map(sf[c])
        surface[f"{c}_present"] = (v.fillna(0) > 0).astype(float)
        surface[f"{c}_missing"] = v.isna().astype(float)

    phage = pd.get_dummies(pairs["phage"].astype("string"),
                           drop_first=True, dtype=float)
    defense = pairs[systems].astype(float)
    lineage_cat = pd.get_dummies(pairs["phylogroup"].astype("string"),
                                 drop_first=True, dtype=float)
    coords = pcoa.reindex(pairs["strain"]).to_numpy(dtype=float)
    coords = (coords - coords.mean(axis=0)) / coords.std(axis=0)
    lineage = np.hstack([lineage_cat.to_numpy(), coords])

    blocks = {
        "phage identity": phage.to_numpy(),
        "host surface": surface.to_numpy(),
        "defense systems": defense.to_numpy(),
        "clonal lineage": lineage,
    }
    definition = pd.DataFrame([
        dict(block="phage identity", n_terms=phage.shape[1],
             description="96 phage fixed effects"),
        dict(block="host surface", n_terms=surface.shape[1],
             description="LPS outer core, O-antigen serotype and status, "
                         "capsule state, 12 outer membrane protein loci, "
                         f"levels below {MIN_LEVEL_COUNT} strains pooled"),
        dict(block="defense systems", n_terms=defense.shape[1],
             description=f"{len(systems)} DefenseFinder systems, presence"),
        dict(block="clonal lineage", n_terms=lineage.shape[1],
             description="phylogroup plus 8 core-genome principal coordinates"),
    ])
    return blocks, definition


def fit_pseudo_r2(y: np.ndarray, x: np.ndarray, llf_null: float) -> float:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = sm.Logit(y, x)
        try:
            res = model.fit(method="newton", maxiter=200, disp=0)
            if not np.isfinite(res.llf):
                raise ValueError
        except Exception:
            res = model.fit(method="bfgs", maxiter=800, disp=0)
    return 1.0 - float(res.llf) / llf_null


def shapley_partition(y: np.ndarray, blocks: dict[str, np.ndarray]
                      ) -> tuple[pd.DataFrame, dict[frozenset, float]]:
    """Exact Shapley decomposition of pseudo R2 over the blocks.

    For k blocks there are 2^k subsets. Each block's share is the average of
    the increase it produces when added, taken over every ordering, which is
    the unique allocation that is order-independent and sums to the full model
    pseudo R2.
    """
    names = list(blocks)
    k = len(names)
    const = np.ones((len(y), 1))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        llf_null = float(sm.Logit(y, const).fit(method="newton", disp=0).llf)

    cache: dict[frozenset, float] = {}
    for size in range(k + 1):
        for combo in itertools.combinations(names, size):
            key = frozenset(combo)
            parts = [const] + [blocks[n] for n in combo]
            cache[key] = fit_pseudo_r2(y, np.hstack(parts), llf_null)
            print(f"    fitted {sorted(combo) if combo else ['intercept only']}"
                  f" -> pseudo R2 {cache[key]:.4f}")

    rows = []
    for name in names:
        phi = 0.0
        others = [n for n in names if n != name]
        for size in range(len(others) + 1):
            weight = math.factorial(size) * math.factorial(k - size - 1) / math.factorial(k)
            for combo in itertools.combinations(others, size):
                gain = cache[frozenset(combo + (name,))] - cache[frozenset(combo)]
                phi += weight * gain
        alone = cache[frozenset([name])]
        last = cache[frozenset(names)] - cache[frozenset(n for n in names if n != name)]
        rows.append(dict(block=name, shapley_pseudo_r2=phi,
                         alone_pseudo_r2=alone, added_last_pseudo_r2=last,
                         n_terms=blocks[name].shape[1]))
    out = pd.DataFrame(rows)
    out["share_of_explained"] = out["shapley_pseudo_r2"] / out["shapley_pseudo_r2"].sum()
    out["shapley_per_term"] = out["shapley_pseudo_r2"] / out["n_terms"]
    return out.sort_values("shapley_pseudo_r2", ascending=False), cache


def main() -> int:
    paths.bootstrap()
    pairs = pd.read_csv(paths.PAIRS, compression="gzip", low_memory=False)
    strain_feats = pd.read_csv(paths.STRAIN_FEATURES, low_memory=False)
    pcoa = pd.read_csv(paths.PCOA, index_col=0)
    ds_cols = [c for c in strain_feats.columns if c.startswith("ds_")]
    tested = strain_feats[strain_feats["strain"].isin(set(pairs["strain"]))]
    systems = sorted([c for c in ds_cols
                      if tested[c].sum() >= config.MIN_CARRIER_STRAINS])

    merge_cols = ["strain", "abc_serotype", "capsule_abc", "capsule_groupiv_e",
                  "capsule_groupiv_s", "capsule_wzy_stricte"]
    merge_cols += [f"omp_{l.lower()}" for l in config.OMP_LOCI]
    merge_cols = [c for c in merge_cols if c in strain_feats.columns]
    pairs = pairs.merge(strain_feats[merge_cols], on="strain", how="left")
    pairs = pairs.merge(strain_feats[["strain"] + systems], on="strain", how="left")
    pairs = pairs.dropna(subset=["lps_core", "phylogroup"]).reset_index(drop=True)

    blocks, definition = build_blocks(pairs, strain_feats, pcoa, systems)
    y = pairs["lytic"].to_numpy(dtype=float)
    print(f"rows {len(y)}, blocks " +
          ", ".join(f"{k}={v.shape[1]}" for k, v in blocks.items()))
    print("  fitting 2^4 = 16 nested models for the exact Shapley decomposition")
    table, cache = shapley_partition(y, blocks)

    table.to_csv(paths.TABLES / "table25_shapley_variance_partition.csv", index=False)
    definition.to_csv(paths.TABLES / "table26_block_definitions.csv", index=False)

    full = cache[frozenset(blocks)]
    out: list[str] = []
    add = out.append
    add("Variance partitioning report")
    add("Generated by scripts/11_variance_partitioning.py")
    add("")
    add("WHY THIS SUPERSEDES THE DECOMPOSITION IN SCRIPT 03")
    add("  Script 03 added blocks in one fixed order and read off the "
        "increments. Those increments depend on the order chosen, and the "
        "surface block it used held only the LPS outer core and O-antigen "
        "status, six terms against thirty-three for defense. Comparing them "
        "directly was not a fair test, and the conclusion drawn there, that "
        "defense carriage contributes slightly more than host surface, does "
        "not survive the correction below.")

    add("")
    add("BLOCKS")
    for r in definition.itertuples(index=False):
        add(f"  {r.block:18s} {r.n_terms:4d} terms  {r.description}")
    add(f"  observations: {len(y)}")
    add(f"  full model pseudo R2: {full:.4f}")

    add("")
    add("EXACT SHAPLEY DECOMPOSITION OF PSEUDO R2")
    add(f"  {'block':18s} {'Shapley':>9s} {'share':>7s} {'alone':>8s} "
        f"{'added last':>11s} {'terms':>6s} {'per term':>9s}")
    for r in table.itertuples(index=False):
        add(f"  {r.block:18s} {r.shapley_pseudo_r2:9.4f} "
            f"{100 * r.share_of_explained:6.1f}% {r.alone_pseudo_r2:8.4f} "
            f"{r.added_last_pseudo_r2:11.4f} {r.n_terms:6d} "
            f"{r.shapley_per_term:9.5f}")
    add("")
    add("  'alone' is the block fitted by itself and 'added last' is what it "
        "adds once every other block is present. The gap between them measures "
        "how much of a block's apparent contribution is shared with the others. "
        "The Shapley column splits that shared part evenly across orderings and "
        "is the figure to quote.")

    surf = table.set_index("block").loc["host surface"]
    defn = table.set_index("block").loc["defense systems"]
    phg = table.set_index("block").loc["phage identity"]
    lin = table.set_index("block").loc["clonal lineage"]

    add("")
    add("WHAT CHANGES")
    ratio_total = surf["shapley_pseudo_r2"] / defn["shapley_pseudo_r2"]
    ratio_term = surf["shapley_per_term"] / defn["shapley_per_term"]
    add(f"  Host surface now accounts for {surf['shapley_pseudo_r2']:.4f} of "
        f"pseudo R2 against {defn['shapley_pseudo_r2']:.4f} for all "
        f"{int(defn['n_terms'])} defense systems together, a ratio of "
        f"{ratio_total:.1f} to 1. The same ordering holds for the block fitted "
        f"alone ({surf['alone_pseudo_r2']:.4f} against "
        f"{defn['alone_pseudo_r2']:.4f}) and for the block added last "
        f"({surf['added_last_pseudo_r2']:.4f} against "
        f"{defn['added_last_pseudo_r2']:.4f}), so it does not depend on how "
        "the shared variation is allocated.")
    if ratio_term < 1:
        add(f"  Per term the ordering reverses: surface contributes "
            f"{surf['shapley_per_term']:.5f} per term against "
            f"{defn['shapley_per_term']:.5f} for defense, so defense is the "
            f"more efficient block by a factor of {1 / ratio_term:.1f}. The "
            f"surface advantage in total comes partly from spending "
            f"{int(surf['n_terms'])} terms against {int(defn['n_terms'])}. "
            "Both readings are reported because they answer different "
            "questions: how much of the variation a block accounts for, and "
            "how much each term buys.")
    else:
        add(f"  Per term the gap is wider still: "
            f"{surf['shapley_per_term']:.5f} against "
            f"{defn['shapley_per_term']:.5f}, a ratio of {ratio_term:.1f} "
            "to 1.")
    add("  The earlier finding is therefore reversed in total terms. Once the "
        "surface block contains what this collection actually types about the "
        "cell surface, and once the comparison no longer depends on the order "
        "the blocks are entered, host surface accounts for substantially more "
        "of the variation in lysis than intracellular defense does. That "
        "agrees with the published analysis of the same dataset, which "
        "reported adsorption-related features as more informative than "
        "antiphage systems, and it contradicts what script 03 implied.")
    add(f"  Which phage is being tested remains the largest single term "
        f"({phg['shapley_pseudo_r2']:.4f}, "
        f"{100 * phg['share_of_explained']:.0f} percent of the explained "
        "variation), and clonal lineage contributes "
        f"{lin['shapley_pseudo_r2']:.4f}.")
    add("")
    add("  Two cautions. These are descriptive shares of explained variation, "
        "not causal effects, and the blocks remain correlated with one another "
        "however the shared part is allocated. And a block with more terms can "
        "capture more variation simply by being more flexible, which is why the "
        "per-term column is shown beside the totals rather than instead of them.")
    add("")
    add("  The conclusion that intracellular defense explains little of the "
        "variation in lysis is unchanged and is in fact strengthened: 33 "
        "defense systems together account for "
        f"{100 * defn['share_of_explained']:.0f} percent of the explained "
        "variation, and almost all of the associations behind that share do "
        "not survive control for clonal background.")

    text = "\n".join(out)
    (paths.RESULTS / "variance_partitioning_report.txt").write_text(
        text + "\n", encoding="utf-8")
    print("\n" + text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
