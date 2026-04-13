#!/usr/bin/env python3
"""
FLAMINGO-style 3D reconstruction using low-rank matrix completion + MDS.

Implements the core FLAMINGO algorithm (Wang et al., Nature Communications 2022)
in Python: converts contact frequencies to distances, completes the sparse
distance matrix via iterative nuclear norm minimization (soft-thresholded SVD),
then applies classical MDS to obtain 3D coordinates.

This is a faithful reimplementation of the algorithm described in the paper,
using the FLAMINGOr R package as reference for parameter choices.

Usage:
    python run_flamingo_hla.py -i contact_matrix.txt -o coords.txt [--alpha -0.25]
"""

import argparse
import numpy as np
import sys


def flamingo_reconstruct(contact_matrix, alpha=-0.25, rank_k=50,
                         max_iter=100, lambda_reg=0.01, tol=1e-4):
    """
    FLAMINGO: Low-rank matrix completion + classical MDS.

    Parameters:
        contact_matrix: NxN symmetric contact frequency matrix (sparse OK)
        alpha: IF-to-distance conversion exponent (default -0.25, per FLAMINGO)
        rank_k: target rank for low-rank completion
        max_iter: maximum iterations for matrix completion
        lambda_reg: regularization (soft threshold on singular values)
        tol: convergence tolerance

    Returns:
        Nx3 numpy array of 3D coordinates
    """
    n = contact_matrix.shape[0]
    rank_k = min(rank_k, n - 1)

    # Step 1: Convert interaction frequency to pairwise distance
    # d_ij = IF_ij^alpha (FLAMINGO uses alpha = -0.25)
    print(f"  Step 1: IF-to-distance conversion (alpha={alpha})")
    pd = np.zeros_like(contact_matrix, dtype=np.float64)
    mask = contact_matrix > 0
    pd[mask] = np.power(contact_matrix[mask], alpha)

    observed = mask.copy()
    n_observed = observed.sum()
    n_total = n * n
    print(f"    Observed: {n_observed}/{n_total} ({100*n_observed/n_total:.1f}%)")

    # Step 2: Low-rank matrix completion via iterative soft-thresholded SVD
    # Initialize missing entries with row/column means of observed values
    print(f"  Step 2: Low-rank matrix completion (rank={rank_k})")
    for i in range(n):
        obs_row = observed[i, :]
        if obs_row.any() and (~obs_row).any():
            pd[i, ~obs_row] = np.mean(pd[i, obs_row])
    # Symmetrize
    pd = (pd + pd.T) / 2.0
    np.fill_diagonal(pd, 0.0)

    # Iterative nuclear norm minimization
    for it in range(1, max_iter + 1):
        # Truncated SVD
        U, s, Vt = np.linalg.svd(pd, full_matrices=False)
        # Soft threshold singular values
        s_thresh = np.maximum(s[:rank_k] - lambda_reg, 0.0)
        # Reconstruct low-rank approximation
        pd_new = (U[:, :rank_k] * s_thresh) @ Vt[:rank_k, :]
        # Preserve observed entries
        pd_new[observed] = pd[observed]
        pd_new = (pd_new + pd_new.T) / 2.0
        np.fill_diagonal(pd_new, 0.0)

        # Check convergence
        diff = np.linalg.norm(pd_new - pd) / (np.linalg.norm(pd) + 1e-10)
        if it % 20 == 0 or diff < tol:
            print(f"    Iter {it}: change = {diff:.6f}")
        pd = pd_new
        if diff < tol:
            print(f"    Converged at iteration {it}")
            break

    # Step 3: Classical MDS on completed distance matrix
    print("  Step 3: Classical MDS")
    pd = np.maximum(pd, 0.0)
    np.fill_diagonal(pd, 0.0)

    D2 = pd ** 2
    row_mean = D2.mean(axis=1, keepdims=True)
    col_mean = D2.mean(axis=0, keepdims=True)
    total_mean = D2.mean()
    B = -0.5 * (D2 - row_mean - col_mean + total_mean)

    # Eigendecomposition (only need top 3)
    eigenvalues, eigenvectors = np.linalg.eigh(B)
    # eigh returns ascending order; take last 3
    idx = np.argsort(eigenvalues)[::-1][:3]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Clamp negative eigenvalues
    eigenvalues = np.maximum(eigenvalues, 1e-10)

    coords = eigenvectors * np.sqrt(eigenvalues)
    return coords


def main():
    parser = argparse.ArgumentParser(
        description="FLAMINGO-style 3D reconstruction (low-rank completion + MDS)")
    parser.add_argument("-i", "--input", required=True,
                        help="Input NxN contact matrix (whitespace-delimited)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output coordinate file (x y z per line)")
    parser.add_argument("--alpha", type=float, default=-0.25,
                        help="IF-to-distance exponent (default: -0.25, per FLAMINGO)")
    parser.add_argument("--rank", type=int, default=50,
                        help="Target rank for matrix completion (default: 50)")
    parser.add_argument("--max-iter", type=int, default=100,
                        help="Max iterations for completion (default: 100)")
    args = parser.parse_args()

    print(f"FLAMINGO reconstruction: {args.input}")
    mat = np.loadtxt(args.input)
    mat = (mat + mat.T) / 2.0
    print(f"Matrix shape: {mat.shape}")

    coords = flamingo_reconstruct(mat, alpha=args.alpha, rank_k=args.rank,
                                  max_iter=args.max_iter)

    np.savetxt(args.output, coords, fmt="%.6f", delimiter="\t")
    print(f"Coordinates ({coords.shape[0]} x {coords.shape[1]}) saved to {args.output}")


if __name__ == "__main__":
    main()
