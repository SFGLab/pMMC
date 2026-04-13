#!/usr/bin/env python3
# Compare pMMC / FLAMINGO / LorDG HLA PDBs on *shape* metrics
# that are rotation/scale invariant:
#   - N beads
#   - Bounding-box diameter (after unit rescaling)
#   - Radius of gyration (Rg) after per-structure rescale to mean-nn = 1
#   - Mean / median nearest-neighbour (consecutive-bead) distance
#   - Dist(|i-j|) distribution decay (short-range vs long-range ratio)
#   - Pearson corr of pairwise distances vs. the input Hi-C
# Writes a summary to stdout + CSV.

import os
import numpy as np

PDBS = {
    "pMMC":     r"D:/git/pMMC/output_hla_comparison/pmmc/hla_pmmc.pdb",
    "FLAMINGO": r"D:/git/pMMC/output_hla_comparison/flamingo/hla_flamingo.pdb",
    "LorDG":    r"D:/git/pMMC/output_hla_comparison/lordg/hla_lordg.pdb",
}
INPUT_MATRIX = r"D:/git/pMMC/data/gm12878_chr6_hla_25kb.txt"
OUT_CSV      = r"D:/git/pMMC/output_hla_comparison/hla_shape_metrics.csv"


def read_ca(path, anchor_occupancy=None):
    xs = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if anchor_occupancy is not None:
                try:
                    occ = float(line[54:60])
                    if abs(occ - anchor_occupancy) > 0.01:
                        continue
                except ValueError:
                    pass
            # Try strict PDB columns first (31-38, 39-46, 47-54)
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                tok = line.split()
                # Fallback: find three consecutive float tokens
                floats = []
                for t in tok:
                    try:
                        floats.append(float(t))
                    except ValueError:
                        if floats and len(floats) < 3:
                            floats = []
                if len(floats) >= 3:
                    x, y, z = floats[:3]
                else:
                    continue
            xs.append((x, y, z))
    return np.array(xs, dtype=float)


def metrics(name, coords, hic):
    n = len(coords)
    # Consecutive-bead distances
    d_nn = np.linalg.norm(np.diff(coords, axis=0), axis=1)
    mean_nn = float(np.mean(d_nn))
    # Normalise so mean_nn == 1 (unit-free shape)
    if mean_nn > 0:
        c = coords / mean_nn
    else:
        c = coords
    # Bounding box diameter (normalised units)
    bbox = c.max(axis=0) - c.min(axis=0)
    bbox_diag = float(np.linalg.norm(bbox))
    # Radius of gyration (normalised units)
    center = c.mean(axis=0)
    rg = float(np.sqrt(np.mean(np.sum((c - center) ** 2, axis=1))))
    # End-to-end distance (normalised)
    e2e = float(np.linalg.norm(c[-1] - c[0]))
    # Short- vs long-range: mean distance for |i-j|<=5 vs |i-j|>=n/4
    i, j = np.triu_indices(n, k=1)
    dij = np.linalg.norm(c[i] - c[j], axis=1)
    sep = j - i
    short = dij[sep <= 5]
    longr = dij[sep >= n // 4]
    sr_mean = float(short.mean()) if len(short) else float("nan")
    lr_mean = float(longr.mean()) if len(longr) else float("nan")
    # Correlation with Hi-C: d ~ 1 / IF^alpha ⇒ corr(dist, IF) should be negative
    if hic is not None and hic.shape[0] == n:
        hij = hic[i, j]
        mask = (hij > 0) & (dij > 1e-6)
        if mask.sum() > 10:
            pear = float(np.corrcoef(dij[mask], np.log1p(hij[mask]))[0, 1])
            pear_loglog = float(np.corrcoef(np.log(dij[mask]),
                                             np.log1p(hij[mask]))[0, 1])
            from scipy.stats import spearmanr
            spearman = float(spearmanr(dij[mask], hij[mask])[0])
        else:
            pear = float("nan")
            pear_loglog = float("nan")
            spearman = float("nan")
    else:
        pear = float("nan")
        pear_loglog = float("nan")
        spearman = float("nan")
    return {
        "tool":          name,
        "n":             n,
        "mean_nn":       mean_nn,
        "bbox_diag/nn":  bbox_diag,
        "Rg/nn":         rg,
        "e2e/nn":        e2e,
        "short_mean":    sr_mean,
        "long_mean":     lr_mean,
        "long/short":    lr_mean / sr_mean if sr_mean else float("nan"),
        "pearson_d_vs_logIF": pear,
        "pearson_logd_vs_logIF": pear_loglog,
        "spearman_d_vs_IF":  spearman,
    }


hic = None
if os.path.exists(INPUT_MATRIX):
    hic = np.loadtxt(INPUT_MATRIX)

rows = []
for name, path in PDBS.items():
    if not os.path.exists(path):
        print(f"  MISSING: {path}")
        continue
    # pMMC's PDB contains anchor-level (occ=4) and sub-anchor (occ=5).
    # Filter to anchor level for parity with FLAMINGO/LorDG.
    occ_filter = 4.0 if name == "pMMC" else None
    coords = read_ca(path, anchor_occupancy=occ_filter)
    if len(coords) == 0:
        print(f"  EMPTY: {path}")
        continue
    m = metrics(name, coords, hic)
    rows.append(m)

if not rows:
    raise SystemExit("no PDBs to compare")

cols = list(rows[0].keys())
def _fmt(v):
    return f"{v:.4g}" if isinstance(v, float) else str(v)
widths = {c: max(len(c), max(len(_fmt(r[c])) for r in rows)) for c in cols}
print(" | ".join(c.ljust(widths[c]) for c in cols))
print("-+-".join("-" * widths[c] for c in cols))
for r in rows:
    cells = []
    for c in cols:
        v = r[c]
        s = f"{v:.4g}" if isinstance(v, float) else str(v)
        cells.append(s.ljust(widths[c]))
    print(" | ".join(cells))

with open(OUT_CSV, "w") as fh:
    fh.write(",".join(cols) + "\n")
    for r in rows:
        fh.write(",".join(f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
print(f"\nwrote {OUT_CSV}")
