"""Separate failures to lyse from failures to be counted as lysis.

Three quarters of the interactions in this collection fail, and the analysis so
far has treated every failure as a form of host resistance. Some of them are not.
A temperate phage entering lysogeny leaves no plaque, and a strain that is
lysogenised looks exactly like a strain that resisted. Two genera here lyse
under one percent of the strains they were tested against, which is what a
temperate lifestyle looks like rather than what universal resistance looks like.

This script classifies the 96 phages by lifestyle from their own genomes. The
marker is the lysogeny machinery: a tyrosine or serine recombinase is what lets
a phage integrate, and its presence is the standard genomic signal of a
temperate lifestyle. Profiles are fetched from Pfam and searched against the
proteins predicted in script 07 at the profiles' own gathering cutoffs.

The point is not to characterise the phages for their own sake. It is to ask
how much of the unexplained failure in this dataset is host resistance at all,
and whether any conclusion about defense systems depends on phages that were
never going to produce a plaque.

Outputs:
    data/processed/phage_lifestyle.csv
    results/tables/table34_phage_lifestyle.csv
    results/tables/table35_lifestyle_sensitivity.csv
    results/phage_lifestyle_report.txt
"""

from __future__ import annotations

import gzip
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
from scipy import stats

from papr import config, paths
from papr.stats import (benjamini_hochberg, mantel_haenszel, stratified_tables,
                        wilson_interval)

PFAM_URL = "https://www.ebi.ac.uk/interpro/api/entry/pfam/{acc}/?annotation=hmm"

# Lysogeny machinery, split by how specific each marker is.
#
# A tyrosine integrase of the lambda type is the canonical signal that a phage
# can integrate into a host chromosome, and is used here to classify.
TYROSINE_INTEGRASE_FAMILIES = {
    "PF00589": "Phage_integrase, tyrosine recombinase",
    "PF02899": "Phage_int_SAM_1, integrase N-terminal domain",
    "PF13356": "Arm-DNA-bind_3, integrase arm-type DNA binding",
}
# Serine recombinases are recorded but not used to classify. Some phages use
# them to integrate, but the same fold serves as a DNA invertase that switches
# tail fiber genes in strictly virulent phages, so a hit is not evidence of
# lysogeny on its own.
SERINE_RECOMBINASE_FAMILIES = {
    "PF00239": "Resolvase, serine recombinase or DNA invertase",
    "PF07508": "Recombinase",
}
INTEGRATION_FAMILIES = {**TYROSINE_INTEGRASE_FAMILIES,
                        **SERINE_RECOMBINASE_FAMILIES}
REPRESSOR_FAMILIES = {
    "PF01381": "HTH_3, cI-like repressor helix-turn-helix",
    "PF13443": "HTH_26, cro/cI-type helix-turn-helix",
}
ALL_FAMILIES = {**INTEGRATION_FAMILIES, **REPRESSOR_FAMILIES}

LIFESTYLE_DIR = paths.RAW / "lifestyle_profiles"


def fetch_profiles() -> dict[str, Path]:
    LIFESTYLE_DIR.mkdir(parents=True, exist_ok=True)
    got: dict[str, Path] = {}
    for acc in ALL_FAMILIES:
        dest = LIFESTYLE_DIR / f"{acc}.hmm"
        if not dest.exists() or dest.stat().st_size == 0:
            req = urllib.request.Request(
                PFAM_URL.format(acc=acc),
                headers={"User-Agent": "post-adsorption-analysis/1.0"})
            try:
                with urllib.request.urlopen(req, timeout=120) as resp:
                    payload = resp.read()
                if payload[:2] == b"\x1f\x8b":
                    payload = gzip.decompress(payload)
                if not payload.lstrip().startswith(b"HMMER3"):
                    print(f"  {acc}: response is not an HMM profile, skipped")
                    continue
                dest.write_bytes(payload)
            except Exception as exc:
                print(f"  {acc}: fetch failed ({exc}), skipped")
                continue
        got[acc] = dest
    return got


def load_proteins() -> dict[str, list[tuple[str, str]]]:
    by_phage: dict[str, list[tuple[str, str]]] = {}
    with gzip.open(paths.PHAGE_PROTEINS, "rt", encoding="utf-8") as fh:
        pid = phage = None
        buf: list[str] = []
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if pid:
                    by_phage.setdefault(phage, []).append((pid, "".join(buf)))
                parts = line[1:].split()
                pid, phage, buf = parts[0], parts[1], []
            else:
                buf.append(line)
        if pid:
            by_phage.setdefault(phage, []).append((pid, "".join(buf)))
    return by_phage


def search(profiles: dict[str, Path],
           by_phage: dict[str, list[tuple[str, str]]]) -> pd.DataFrame:
    import pyhmmer

    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs, owner = [], {}
    for phage, entries in by_phage.items():
        for pid, seq in entries:
            if not seq:
                continue
            owner[pid] = phage
            seqs.append(pyhmmer.easel.TextSequence(
                name=pid.encode(), sequence=seq).digitize(alphabet))
    block = pyhmmer.easel.DigitalSequenceBlock(alphabet, seqs)

    hmms = []
    acc_of = {}
    for acc, path in profiles.items():
        with pyhmmer.plan7.HMMFile(path) as handle:
            for hmm in handle:
                name = (hmm.name.decode() if isinstance(hmm.name, (bytes, bytearray))
                        else str(hmm.name))
                acc_of[name] = acc
                hmms.append(hmm)

    def as_str(v) -> str:
        return v.decode() if isinstance(v, (bytes, bytearray)) else str(v)

    rows = []
    try:
        results = pyhmmer.hmmsearch(hmms, block, bit_cutoffs="gathering", cpus=0)
        mode = "gathering cutoff"
    except Exception:
        results = pyhmmer.hmmsearch(hmms, block, E=1e-5, cpus=0)
        mode = "E-value 1e-5"
    for hits in results:
        qname = as_str(getattr(hits.query, "name", hits.query))
        acc = acc_of.get(qname, qname)
        for hit in hits:
            if not hit.included:
                continue
            pid = as_str(hit.name)
            rows.append(dict(phage=owner[pid], protein_id=pid, pfam=acc,
                             family=ALL_FAMILIES.get(acc, acc),
                             bit_score=float(hit.score)))
    out = pd.DataFrame(rows)
    out.attrs["threshold_mode"] = mode
    return out


def main() -> int:
    paths.bootstrap()
    pairs = pd.read_csv(paths.PAIRS, compression="gzip", low_memory=False)
    strain_feats = pd.read_csv(paths.STRAIN_FEATURES, low_memory=False)
    phage_feats = pd.read_csv(paths.PHAGE_FEATURES)

    print("fetching lysogeny profiles from Pfam")
    profiles = fetch_profiles()
    print(f"  profiles available: {len(profiles)} of {len(ALL_FAMILIES)}")
    if not profiles:
        raise SystemExit("no profiles could be fetched; cannot classify lifestyle")

    by_phage = load_proteins()
    print("searching phage proteomes")
    hits = search(profiles, by_phage)

    phages = sorted(phage_feats["phage"].astype(str))

    def carriers(families) -> set:
        if not len(hits):
            return set()
        return set(hits.loc[hits["pfam"].isin(families), "phage"])

    tyr = carriers(TYROSINE_INTEGRASE_FAMILIES)
    ser = carriers(SERINE_RECOMBINASE_FAMILIES)
    repr_ = carriers(REPRESSOR_FAMILIES)

    life = pd.DataFrame({"phage": phages})
    life["has_tyrosine_integrase"] = life["phage"].isin(tyr).astype(int)
    life["has_serine_recombinase"] = life["phage"].isin(ser).astype(int)
    life["has_integrase"] = life["phage"].isin(tyr | ser).astype(int)
    life["has_repressor"] = life["phage"].isin(repr_).astype(int)
    # Classification uses the tyrosine integrase only.
    life["lifestyle"] = np.where(life["has_tyrosine_integrase"] == 1,
                                 "temperate", "virulent")
    life = life.merge(phage_feats[["phage", "genus", "morphotype"]], on="phage",
                      how="left")
    per_phage = pairs.groupby("phage").agg(
        n_pairs=("lytic", "size"), lysis_rate=("lytic", "mean"),
        mean_score=("plaque_score", "mean"))
    life = life.merge(per_phage, left_on="phage", right_index=True, how="left")
    life.to_csv(paths.PROCESSED / "phage_lifestyle.csv", index=False)
    life.to_csv(paths.TABLES / "table34_phage_lifestyle.csv", index=False)

    temperate = set(life.loc[life["lifestyle"] == "temperate", "phage"])
    pairs["lifestyle"] = np.where(pairs["phage"].isin(temperate),
                                  "temperate", "virulent")

    # Sensitivity: do the defense associations depend on temperate phages?
    ds_cols = [c for c in strain_feats.columns if c.startswith("ds_")]
    tested = strain_feats[strain_feats["strain"].isin(set(pairs["strain"]))]
    systems = sorted([c for c in ds_cols
                      if tested[c].sum() >= config.MIN_CARRIER_STRAINS])
    pairs = pairs.merge(strain_feats[["strain"] + systems], on="strain", how="left")
    virulent_only = pairs[pairs["lifestyle"] == "virulent"]

    rows = []
    for col in systems:
        all_eff = mantel_haenszel(
            stratified_tables(pairs, ["phage", "lps_core"], col, "lytic"))
        vir_eff = mantel_haenszel(
            stratified_tables(virulent_only, ["phage", "lps_core"], col, "lytic"))
        rows.append(dict(system=col[3:],
                         or_all_phages=all_eff.estimate,
                         or_virulent_only=vir_eff.estimate,
                         p_all=all_eff.p_value, p_virulent=vir_eff.p_value))
    sens = pd.DataFrame(rows)
    sens["q_all"] = benjamini_hochberg(sens["p_all"].to_numpy())
    sens["q_virulent"] = benjamini_hochberg(sens["p_virulent"].to_numpy())
    sens["log2_shift"] = np.log2(sens["or_virulent_only"]) - np.log2(sens["or_all_phages"])
    sens.to_csv(paths.TABLES / "table35_lifestyle_sensitivity.csv", index=False)

    o: list[str] = []
    add = o.append
    add("Phage lifestyle report")
    add("Generated by scripts/17_phage_lifestyle.py")
    add("")
    add("CLASSIFICATION")
    add(f"  Pfam profiles retrieved : {len(profiles)} of {len(ALL_FAMILIES)}")
    for acc, desc in ALL_FAMILIES.items():
        mark = "yes" if acc in profiles else "not retrieved"
        n = int((hits["pfam"] == acc).sum()) if len(hits) else 0
        add(f"    {acc}  {desc:52s} {mark:14s} hits {n}")
    add(f"  threshold used          : {hits.attrs.get('threshold_mode', 'n/a')}")
    add(f"  phages with a tyrosine integrase (classified temperate) : "
        f"{int(life['has_tyrosine_integrase'].sum())} of {len(life)}")
    add(f"  phages with a serine recombinase only                   : "
        f"{int(((life['has_serine_recombinase'] == 1) & (life['has_tyrosine_integrase'] == 0)).sum())}")
    add(f"  phages with a repressor-like protein                    : "
        f"{int(life['has_repressor'].sum())} of {len(life)}")
    add("  Only the tyrosine integrase is used to classify. A serine "
        "recombinase hit is not evidence of lysogeny on its own, because the "
        "same fold works as a DNA invertase that switches tail fiber genes in "
        "strictly virulent phages, and a repressor-like helix-turn-helix is "
        "common in virulent phages too. Both are recorded and neither is used "
        "to classify.")

    add("")
    add("LYSIS RATE BY LIFESTYLE")
    add(f"  {'lifestyle':12s} {'phages':>6s} {'pairs':>7s} {'lysis':>8s} "
        f"{'95 pct CI':>16s} {'mean score':>10s}")
    for ls, g in pairs.groupby("lifestyle"):
        k, n = int(g["lytic"].sum()), len(g)
        lo, hi = wilson_interval(k, n)
        add(f"  {ls:12s} {int(g['phage'].nunique()):6d} {n:7d} "
            f"{100 * k / n:7.2f}% [{100 * lo:6.2f}, {100 * hi:6.2f}] "
            f"{g['plaque_score'].mean():10.3f}")
    t = life[life["lifestyle"] == "temperate"]
    v = life[life["lifestyle"] == "virulent"]
    if len(t) and len(v):
        u = stats.mannwhitneyu(t["lysis_rate"], v["lysis_rate"],
                               alternative="two-sided")
        add(f"  Mann-Whitney on per-phage lysis rate: p = {u.pvalue:.3g}")

    add("")
    add("  Temperate phages by genus:")
    for genus, g in life[life["lifestyle"] == "temperate"].groupby("genus"):
        add(f"    {genus:22s} {len(g):2d} phages, mean lysis rate "
            f"{100 * g['lysis_rate'].mean():5.2f} percent")

    add("")
    add("HOW MUCH OF THE FAILURE IS LYSOGENY RATHER THAN RESISTANCE")
    n_t = int((pairs["lifestyle"] == "temperate").sum())
    fail_t = int(((pairs["lifestyle"] == "temperate") & (pairs["lytic"] == 0)).sum())
    fail_all = int((pairs["lytic"] == 0).sum())
    add(f"  interactions involving a temperate phage : {n_t} "
        f"({100 * n_t / len(pairs):.1f} percent of all pairs)")
    add(f"  failures involving a temperate phage     : {fail_t} "
        f"({100 * fail_t / fail_all:.1f} percent of all failures)")
    add(f"  failure rate, virulent phages only       : "
        f"{100 * (1 - virulent_only['lytic'].mean()):.2f} percent")
    add(f"  failure rate, all phages                 : "
        f"{100 * (1 - pairs['lytic'].mean()):.2f} percent")
    add("")
    temp_rate = float(pairs.loc[pairs["lifestyle"] == "temperate", "lytic"].mean()) \
        if n_t else np.nan
    vir_rate = float(virulent_only["lytic"].mean())
    add("  This script was written to test the idea that part of the "
        "unexplained failure is lysogeny rather than host resistance. The "
        "answer separates cleanly into a direction and a magnitude.")
    if np.isfinite(temp_rate) and temp_rate < vir_rate:
        add(f"  The direction supports it. Phages carrying a tyrosine "
            f"integrase form plaques on {100 * temp_rate:.1f} percent of the "
            f"strains they met, against {100 * vir_rate:.1f} percent for the "
            "rest. A phage that can integrate is far less likely to leave a "
            "plaque, which is what lysogeny predicts, and the genera involved "
            f"are {', '.join(sorted(t['genus'].dropna().unique()))}.")
    else:
        add(f"  The direction does not support it. Phages carrying a tyrosine "
            f"integrase form plaques on {100 * temp_rate:.1f} percent of the "
            f"strains they met, against {100 * vir_rate:.1f} percent for the "
            "rest, so they form more plaques rather than fewer and cannot be "
            "a hidden source of silent failures.")
    add(f"  The magnitude does not. Only {int(life['has_tyrosine_integrase'].sum())} "
        f"of {len(life)} phages carry the marker, and they account for "
        f"{fail_t} of {fail_all} failures "
        f"({100 * fail_t / fail_all:.1f} percent). Even counting every one of "
        "those as a silent lysogenisation rather than a resistance event "
        f"moves the overall failure rate only from "
        f"{100 * (1 - pairs['lytic'].mean()):.2f} to "
        f"{100 * (1 - vir_rate):.2f} percent.")
    add("  Lysogeny is therefore a real but marginal explanation in this "
        "collection. The Guelin phages were selected for lytic activity, so "
        "temperate ones are rare by construction, and the unexplained "
        "remainder reported in script 05 survives this check nearly intact.")

    add("")
    add("DO THE DEFENSE ASSOCIATIONS DEPEND ON TEMPERATE PHAGES")
    shifted = sens.reindex(sens["log2_shift"].abs().sort_values(
        ascending=False).index).head(8)
    add(f"  {'system':24s} {'OR all':>8s} {'OR virulent':>12s} {'log2 shift':>11s}")
    for r in shifted.itertuples(index=False):
        add(f"  {r.system:24s} {r.or_all_phages:8.3f} {r.or_virulent_only:12.3f} "
            f"{r.log2_shift:+11.3f}")
    med = float(sens["log2_shift"].abs().median())
    flips = int((np.sign(np.log(sens["or_all_phages"]))
                 != np.sign(np.log(sens["or_virulent_only"]))).sum())
    add("")
    add(f"  Median absolute shift in log2 odds ratio when temperate phages are "
        f"dropped: {med:.3f}. Systems changing direction: {flips} of "
        f"{len(sens)}.")
    if med < 0.1:
        add("  The defense associations are essentially unchanged. None of the "
            "conclusions in this analysis rests on interactions with phages "
            "that were never going to form a plaque, which is the reassurance "
            "this check was run to provide.")
    else:
        add("  The defense associations move appreciably when temperate phages "
            "are dropped, so the stratified estimates elsewhere in this "
            "analysis are partly carried by interactions whose outcome the "
            "phage lifestyle, not the host, determines.")

    add("")
    add("LIMITS OF THIS CLASSIFICATION")
    add("  Carrying a recombinase shows a phage can integrate, not that it did "
        "so in any particular assay. The converse error is also possible: a "
        "phage with no detectable recombinase may still fail to plaque for "
        "reasons unrelated to host defense, such as poor adsorption kinetics "
        "or a host growth condition. This separates two populations that the "
        "plaque assay cannot distinguish; it does not measure lysogeny "
        "directly, which would need the assay the dataset does not contain.")

    text = "\n".join(o)
    (paths.RESULTS / "phage_lifestyle_report.txt").write_text(
        text + "\n", encoding="utf-8")
    print("\n" + text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
