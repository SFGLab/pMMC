#!/usr/bin/env Rscript
# Run FLAMINGOrLite::flamingo_basic() on a dense NxN contact matrix
# and write (1) raw coordinates TSV, (2) a minimal PDB.
#
# Usage: Rscript run_flamingo_basic.R <matrix.txt> <out_prefix>

suppressPackageStartupMessages(library(FLAMINGOrLite))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript run_flamingo_basic.R <matrix.txt> <out_prefix>\n")
  quit(status = 1)
}
input_file <- args[1]
out_prefix <- args[2]

cat("Reading matrix:", input_file, "\n")
mat <- as.matrix(read.table(input_file))
n <- nrow(mat)
cat("Matrix:", n, "x", ncol(mat), "  sum=", sum(mat), "\n")

# Symmetrise (defensive — cooler output is already symmetric)
mat <- (mat + t(mat)) / 2

# FLAMINGOrLite expects a sparseMatrix
suppressPackageStartupMessages(library(Matrix))
sm <- as(mat, "sparseMatrix")

set.seed(42)
cat("Running flamingo_basic()...\n")
res <- flamingo_basic(
  input_if     = sm,
  sample_rate  = 0.75,
  lambda       = 10,
  r            = 1,
  max_dist     = 0.01,
  error_threshold = 1e-3,
  max_iter     = 500,
  alpha        = -0.25,
  inf_dist     = 4
)

# res is typically a "flamingo_prediction" S4 object with @coordinates (n x 3)
coords <- tryCatch(res@coordinates, error = function(e) NULL)
if (is.null(coords)) coords <- tryCatch(res$coordinates, error = function(e) NULL)
if (is.null(coords)) stop("Could not extract coordinates from flamingo_basic() result")
if (ncol(coords) != 3 && nrow(coords) == 3) coords <- t(coords)

cat("Coords:", nrow(coords), "x", ncol(coords), "\n")

tsv <- paste0(out_prefix, ".tsv")
write.table(coords, tsv, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("wrote", tsv, "\n")

# Write a minimal PDB (CA atoms, chain A, resSeq = bead index)
# Rescale so RMSD between tools is not dominated by arbitrary units.
# FLAMINGO coords are typically in [-1, 1]; scale so max extent ≈ 100 Å for viewing.
span <- max(apply(coords, 2, function(x) diff(range(x))))
if (span > 0) coords <- coords * (100 / span)

pdb <- paste0(out_prefix, ".pdb")
fh <- file(pdb, "w")
writeLines(sprintf("HEADER    FLAMINGO HLA chr6:25-34Mb  n=%d", nrow(coords)), fh)
for (i in seq_len(nrow(coords))) {
  line <- sprintf(
    "ATOM  %5d  CA  GLY A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C",
    i, i, coords[i, 1], coords[i, 2], coords[i, 3]
  )
  writeLines(line, fh)
}
writeLines("END", fh)
close(fh)
cat("wrote", pdb, "\n")
