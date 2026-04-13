#!/usr/bin/env Rscript
# Run FLAMINGO on HLA region (chr6:25M-34M) from ChIA-PET-derived contact matrix
#
# Uses FLAMINGOr's flamingo_basic() which takes a raw interaction frequency
# matrix directly — no .hic/.mcool parsing needed.
#
# Usage: Rscript run_flamingo_hla.R <contact_matrix.txt> <output_coords.txt>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript run_flamingo_hla.R <contact_matrix.txt> <output_coords.txt>\n")
  quit(status = 1)
}

input_file <- args[1]
output_file <- args[2]

cat("Loading FLAMINGO...\n")

# Try loading FLAMINGOrLite, fall back to FLAMINGOr, fall back to manual
pkg_loaded <- FALSE

if (requireNamespace("FLAMINGOrLite", quietly = TRUE)) {
  library(FLAMINGOrLite)
  cat("Using FLAMINGOrLite\n")
  pkg_loaded <- TRUE
} else if (requireNamespace("FLAMINGOr", quietly = TRUE)) {
  library(FLAMINGOr)
  cat("Using FLAMINGOr\n")
  pkg_loaded <- TRUE
}

# Load contact matrix
cat("Loading contact matrix:", input_file, "\n")
mat <- as.matrix(read.table(input_file))
n <- nrow(mat)
cat("Matrix size:", n, "x", n, "\n")

# Symmetrize
mat <- (mat + t(mat)) / 2

if (pkg_loaded) {
  # Use FLAMINGO's core function
  cat("Running flamingo_basic()...\n")
  res <- flamingo_basic(mat)
  coords <- res$coordinates  # Expected: 3 x N matrix
  if (nrow(coords) == 3) {
    coords <- t(coords)  # Transpose to N x 3
  }
} else {
  # Manual FLAMINGO implementation: low-rank matrix completion + MDS
  # Based on the FLAMINGO algorithm (Wang et al., 2022)
  cat("FLAMINGOr not installed, using manual low-rank completion + MDS\n")

  library(Matrix)

  # Step 1: Convert IF to distance matrix
  alpha <- -0.25
  pd <- mat
  pd[pd > 0] <- pd[pd > 0]^alpha
  pd[mat == 0] <- 0  # missing entries

  # Step 2: Low-rank matrix completion via soft-thresholded SVD
  # (Simplified FLAMINGO core: iterative nuclear norm minimization)
  cat("Low-rank matrix completion...\n")

  observed <- mat > 0
  # Initialize missing entries with column/row means
  for (i in 1:n) {
    missing <- !observed[i, ]
    if (any(missing) && any(observed[i, ])) {
      pd[i, missing] <- mean(pd[i, observed[i, ]])
    }
  }
  # Symmetrize
  pd <- (pd + t(pd)) / 2

  # Iterative soft-thresholded SVD (low-rank completion)
  max_iter <- 50
  rank_k <- min(50, n - 1)  # target rank
  lambda <- 0.01

  for (iter in 1:max_iter) {
    # SVD
    s <- svd(pd, nu = rank_k, nv = rank_k)
    # Soft threshold singular values
    d_thresh <- pmax(s$d[1:rank_k] - lambda, 0)
    # Reconstruct
    pd_new <- s$u %*% diag(d_thresh, nrow = rank_k) %*% t(s$v)
    # Only update missing entries
    pd_new[observed] <- pd[observed]
    pd_new <- (pd_new + t(pd_new)) / 2

    # Check convergence
    diff <- norm(pd_new - pd, "F") / (norm(pd, "F") + 1e-10)
    if (iter %% 10 == 0) cat("  Iter", iter, "change:", diff, "\n")
    pd <- pd_new
    if (diff < 1e-4) break
  }

  # Step 3: Classical MDS on completed distance matrix
  cat("Classical MDS...\n")
  pd[pd < 0] <- 0
  diag(pd) <- 0

  # Double centering
  D2 <- pd^2
  H <- diag(n) - 1/n * matrix(1, n, n)
  B <- -0.5 * H %*% D2 %*% H

  # Eigendecomposition
  e <- eigen(B, symmetric = TRUE)
  # Take top 3 positive eigenvalues
  pos_idx <- which(e$values > 0)
  if (length(pos_idx) < 3) {
    cat("Warning: fewer than 3 positive eigenvalues\n")
    pos_idx <- 1:3
    e$values[e$values <= 0] <- 1e-6
  }

  coords <- e$vectors[, pos_idx[1:3]] %*% diag(sqrt(abs(e$values[pos_idx[1:3]])))
}

cat("Output coordinates:", nrow(coords), "x", ncol(coords), "\n")

# Write coordinates
write.table(coords, file = output_file, row.names = FALSE, col.names = FALSE, sep = "\t")
cat("Coordinates written to", output_file, "\n")
