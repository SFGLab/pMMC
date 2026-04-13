#!/usr/bin/env python3
# Run pMMC N times with different seeds, keep the structure with the
# best Pearson(log(d), log(IF)) vs. the input Hi-C matrix.
#
# Usage: python pmmc_best_of_n.py <N> [outdir]

import os
import sys
import shutil
import subprocess
import numpy as np

N = int(sys.argv[1]) if len(sys.argv) > 1 else 8
ROOT = r"D:/git/pMMC"
INI  = r"./output_hla_comparison/pmmc_hla.ini"
TUNED_INI = os.path.join(ROOT, "output_hla_comparison", "pmmc_hla_seeded.ini")
HIC  = r"D:/git/pMMC/data/gm12878_chr6_hla_25kb.txt"
PMMC_DIR = r"D:/git/pMMC/output_hla_comparison/pmmc"
BEST_PDB = r"D:/git/pMMC/output_hla_comparison/pmmc/hla_pmmc.pdb"

WSL = ["wsl", "-d", "Ubuntu", "--", "bash", "-c"]


def read_anchor_ca(path):
    xs = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            try:
                occ = float(line[54:60])
            except ValueError:
                continue
            if abs(occ - 4.0) > 0.01:
                continue
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            xs.append((x, y, z))
    return np.array(xs, dtype=float)


def log_log_pearson(coords, hic):
    n = len(coords)
    i, j = np.triu_indices(n, k=1)
    d = np.linalg.norm(coords[i] - coords[j], axis=1)
    h = hic[i, j]
    m = (h > 0) & (d > 1e-9)
    return float(np.corrcoef(np.log(d[m]), np.log1p(h[m]))[0, 1])


def run_once(seed, run_tag):
    # edit seed in ini
    with open(os.path.join(ROOT, "output_hla_comparison", "pmmc_hla.ini")) as f:
        txt = f.read()
    new = []
    for line in txt.splitlines():
        if line.strip().startswith("seed"):
            new.append(f"seed = {seed}")
        else:
            new.append(line)
    with open(TUNED_INI, "w") as f:
        f.write("\n".join(new) + "\n")

    # nuke caches
    for fn in os.listdir(PMMC_DIR):
        if fn.endswith((".dat", ".heat")) or fn.startswith("loops_hla_pmmc"):
            try: os.remove(os.path.join(PMMC_DIR, fn))
            except OSError: pass

    cmd = (
        f"cd /mnt/d/git/pMMC && "
        f"./build_wsl/pMMC-reconstruct -s ./output_hla_comparison/pmmc_hla_seeded.ini "
        f"> ./output_hla_comparison/pmmc/seed_{seed}.log 2>&1"
    )
    subprocess.run(WSL + [cmd], check=False)
    return BEST_PDB


hic = np.loadtxt(HIC)
best_seed = None
best_score = 1.0  # we want the most negative
results = []

for s in range(42, 42 + N):
    pdb = run_once(s, f"seed{s}")
    if not os.path.exists(pdb):
        print(f"seed {s}: NO PDB")
        continue
    coords = read_anchor_ca(pdb)
    if len(coords) == 0:
        print(f"seed {s}: empty PDB")
        continue
    p = log_log_pearson(coords, hic)
    results.append((s, p))
    print(f"seed {s:3d}: log-log Pearson = {p:+.4f}")
    if p < best_score:
        best_score = p
        best_seed = s
        shutil.copy(pdb, pdb + f".best")

print()
print(f"best seed = {best_seed}  pearson_loglog = {best_score:+.4f}")
# Put the best structure back in the canonical location
if best_seed is not None:
    shutil.copy(BEST_PDB + ".best", BEST_PDB)
