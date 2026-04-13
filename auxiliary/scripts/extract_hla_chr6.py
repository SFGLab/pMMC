#!/usr/bin/env python3
# Extract chr6:25,000,000-34,000,000 @ 5kb contacts from the
# GM12878 mcool and emit three input files:
#   1) NxN dense matrix        -> data/gm12878_chr6_hla_5kb.txt      (FLAMINGO)
#   2) pMMC "pairs" format     -> data/gm12878_chr6_hla_5kb_pairs.txt
#      (chr  pos  chr  pos  count)
#   3) LorDG IF list (1-based) -> data/gm12878_chr6_hla_5kb_IFList.txt
#      (idx1  idx2  count)

import numpy as np
import cooler

MCOOL  = r"D:/git/DATASET/hic_gm12878/4DNFIXP4QG5B.mcool"
OUTDIR = r"D:/git/pMMC/data"
CHR    = "chr6"
START  = 25_000_000
END    = 34_000_000
RES    = 25_000
REGION = f"{CHR}:{START}-{END}"

c = cooler.Cooler(f"{MCOOL}::/resolutions/{RES}")
print(f"Loaded {MCOOL} @ {RES}bp; chroms: {len(c.chromnames)}")

mat = c.matrix(balance=False, sparse=False).fetch(REGION)
mat = np.nan_to_num(mat, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
N = mat.shape[0]
print(f"Matrix: {mat.shape}  nnz={int((mat > 0).sum())}  sum={float(mat.sum()):.0f}")

# (1) dense NxN for FLAMINGO
p_dense = f"{OUTDIR}/gm12878_chr6_hla_25kb.txt"
np.savetxt(p_dense, mat, fmt="%.6g")
print(f"wrote {p_dense}")

# (2) pMMC pairs: only upper triangle, non-zero
p_pairs = f"{OUTDIR}/gm12878_chr6_hla_25kb_pairs.txt"
iu, ju = np.nonzero(np.triu(mat))
with open(p_pairs, "w") as fh:
    for i, j in zip(iu, ju):
        p1 = START + int(i) * RES
        p2 = START + int(j) * RES
        fh.write(f"{CHR}\t{p1}\t{CHR}\t{p2}\t{mat[i, j]:.6g}\n")
print(f"wrote {p_pairs}  ({len(iu)} pairs)")

# (3) LorDG IF list (1-based indices, strict upper triangle, non-zero only)
p_lor = f"{OUTDIR}/gm12878_chr6_hla_25kb_IFList.txt"
iu2, ju2 = np.nonzero(np.triu(mat, k=1))
with open(p_lor, "w") as fh:
    for i, j in zip(iu2, ju2):
        fh.write(f"{int(i) + 1}\t{int(j) + 1}\t{mat[i, j]:.6g}\n")
print(f"wrote {p_lor}  ({len(iu2)} pairs)")
