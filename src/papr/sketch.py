"""FracMinHash sketching and containment, implemented directly.

Used to compare downloaded genomes with one another and to measure how much of
a phage genome is present inside a bacterial genome. A scaled sketch keeps a
deterministic pseudorandom fraction of a sequence's k-mers, so two sketches can
be intersected to estimate both Jaccard similarity and containment without
holding the full k-mer sets. The estimator is unbiased, which matters here
because a prophage occupies under one percent of a host genome and containment
has to be measured against a very large background.

There is no external dependency: the hash is a 64-bit mix applied to the
2-bit encoding of each canonical k-mer.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np

K = 21
SCALED = 100
_MAX_HASH = (1 << 64) - 1
_THRESHOLD = _MAX_HASH // SCALED

_CODE = np.full(256, 255, dtype=np.uint8)
for _b, _v in zip(b"ACGT", range(4)):
    _CODE[_b] = _v
for _b, _v in zip(b"acgt", range(4)):
    _CODE[_b] = _v
_COMP = np.array([3, 2, 1, 0], dtype=np.uint8)


def _mix64(x: np.ndarray) -> np.ndarray:
    """SplitMix64 finaliser: a cheap, well-distributed 64-bit hash."""
    x = x.astype(np.uint64, copy=True)
    x ^= x >> np.uint64(30)
    x *= np.uint64(0xBF58476D1CE4E5B9)
    x ^= x >> np.uint64(27)
    x *= np.uint64(0x94D049BB133111EB)
    x ^= x >> np.uint64(31)
    return x


def read_sequences(path: Path) -> list[np.ndarray]:
    """Read a FASTA, gzipped or not, as arrays of 2-bit codes with 255 for N."""
    opener = gzip.open if str(path).endswith(".gz") else open
    seqs, buf = [], []
    with opener(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if buf:
                    seqs.append(b"".join(buf))
                    buf = []
            else:
                buf.append(line.strip())
    if buf:
        seqs.append(b"".join(buf))
    return [_CODE[np.frombuffer(s, dtype=np.uint8)] for s in seqs if s]


def _canonical_hashes(codes: np.ndarray, k: int = K) -> np.ndarray:
    """Hashes of every valid canonical k-mer in one sequence.

    Built with k vectorised passes over the sequence rather than one pass per
    base. A bacterial genome holds about five million k-mers, so a per-base
    Python loop would dominate the whole analysis; this runs in well under a
    second per genome.
    """
    n = codes.size
    m = n - k + 1
    if m <= 0:
        return np.empty(0, dtype=np.uint64)

    vals = codes.astype(np.uint64)
    comp = np.where(codes < 4, _COMP[np.minimum(codes, 3)], 0).astype(np.uint64)

    two = np.uint64(2)
    fwd = np.zeros(m, dtype=np.uint64)
    rev = np.zeros(m, dtype=np.uint64)
    for j in range(k):
        # Forward: most significant base first.
        fwd = (fwd << two) | vals[j:j + m]
        # Reverse complement: complement of the same base, least significant.
        rev |= comp[j:j + m] << np.uint64(2 * j)

    # A k-mer is usable only if none of its bases was ambiguous.
    bad = (codes >= 4).astype(np.int32)
    cs = np.concatenate(([0], np.cumsum(bad)))
    ok = (cs[k:] - cs[:-k]) == 0

    canon = np.minimum(fwd, rev)[ok]
    if canon.size == 0:
        return np.empty(0, dtype=np.uint64)
    return _mix64(canon)


def sketch_file(path: Path, k: int = K, scaled: int = SCALED) -> np.ndarray:
    """Scaled sketch of a whole FASTA: the hashes below the keep threshold."""
    threshold = _MAX_HASH // scaled
    keep = []
    for codes in read_sequences(path):
        h = _canonical_hashes(codes, k)
        if h.size:
            keep.append(h[h <= threshold])
    if not keep:
        return np.empty(0, dtype=np.uint64)
    return np.unique(np.concatenate(keep))


def containment(small: np.ndarray, large: np.ndarray) -> float:
    """Fraction of the smaller sketch's hashes that appear in the larger."""
    if small.size == 0:
        return float("nan")
    return float(np.intersect1d(small, large, assume_unique=True).size / small.size)


def jaccard(a: np.ndarray, b: np.ndarray) -> float:
    if a.size == 0 and b.size == 0:
        return float("nan")
    inter = np.intersect1d(a, b, assume_unique=True).size
    union = a.size + b.size - inter
    return float(inter / union) if union else float("nan")


def mash_distance(a: np.ndarray, b: np.ndarray, k: int = K) -> float:
    """Convert Jaccard similarity to an approximate substitution distance."""
    j = jaccard(a, b)
    if not np.isfinite(j) or j <= 0:
        return 1.0
    if j >= 1:
        return 0.0
    return float(-np.log(2 * j / (1 + j)) / k)
