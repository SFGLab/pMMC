#!/usr/bin/env python3
"""
ShRec3D — Shortest-path Reconstruction in 3D

A Python implementation of the ShRec3D algorithm for 3D chromatin structure
reconstruction from Hi-C contact matrices.

Reference:
    Lesne A, Riposo J, Roger P, Cournac A, Mozziconacci J.
    "3D genome reconstruction from chromosomal contacts."
    Nature Methods, 11(11):1141-1143, 2014.

Algorithm:
    1. Convert contact frequencies to distances: d_ij = 1 / f_ij^alpha
    2. Complete missing distances via shortest-path (Floyd-Warshall or Dijkstra)
    3. Apply classical MDS (Torgerson) to obtain 3D coordinates

Usage:
    python shrec3d.py -i contact_matrix.txt -o coords.txt [--alpha 1.0]
"""

import argparse
import numpy as np
import sys


def contact_to_distance(contact_matrix, alpha=1.0):
    """Convert contact frequency matrix to distance matrix.

    d_ij = 1 / f_ij^alpha  for f_ij > 0
    d_ij = inf              for f_ij = 0
    """
    n = contact_matrix.shape[0]
    dist = np.full((n, n), np.inf)
    np.fill_diagonal(dist, 0.0)

    mask = contact_matrix > 0
    dist[mask] = 1.0 / np.power(contact_matrix[mask], alpha)

    return dist


def shortest_path_completion(dist_matrix):
    """Complete distance matrix using Floyd-Warshall shortest paths.

    This fills in missing pairwise distances by finding the shortest
    path through known contacts, following the ShRec3D approach.
    """
    n = dist_matrix.shape[0]
    D = dist_matrix.copy()

    # Floyd-Warshall with progress reporting
    print(f"  Floyd-Warshall on {n}x{n} matrix...", end="", flush=True)
    for k in range(n):
        if k % 200 == 0 and k > 0:
            print(f" {k}/{n}", end="", flush=True)
        # Vectorized update for all (i,j) pairs through intermediate k
        new_dist = D[:, k:k+1] + D[k:k+1, :]
        np.minimum(D, new_dist, out=D)
    print(" done")

    return D


def classical_mds(dist_matrix, ndim=3):
    """Classical (Torgerson) multidimensional scaling.

    Given a complete distance matrix, recover ndim-dimensional coordinates
    by eigendecomposition of the double-centered squared distance matrix.
    """
    n = dist_matrix.shape[0]
    D_sq = dist_matrix ** 2

    # Double centering: B = -0.5 * J * D^2 * J, where J = I - 1/n * 11^T
    row_mean = D_sq.mean(axis=1, keepdims=True)
    col_mean = D_sq.mean(axis=0, keepdims=True)
    total_mean = D_sq.mean()
    B = -0.5 * (D_sq - row_mean - col_mean + total_mean)

    # Eigendecomposition (only need top ndim eigenvalues)
    print(f"  Eigendecomposition ({n}x{n})...", end="", flush=True)
    eigenvalues, eigenvectors = np.linalg.eigh(B)
    print(" done")

    # Take the ndim largest eigenvalues (eigh returns ascending order)
    idx = np.argsort(eigenvalues)[::-1][:ndim]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Clamp negative eigenvalues to small positive (numerical stability)
    eigenvalues = np.maximum(eigenvalues, 1e-10)

    # Coordinates = V * sqrt(Lambda)
    coords = eigenvectors * np.sqrt(eigenvalues)

    return coords


def shrec3d(contact_matrix, alpha=1.0, ndim=3):
    """Full ShRec3D pipeline: contact -> distance -> shortest-path -> MDS.

    Parameters:
        contact_matrix: NxN numpy array of contact frequencies
        alpha: exponent for frequency-to-distance conversion (default 1.0)
        ndim: number of dimensions for output (default 3)

    Returns:
        Nx3 numpy array of 3D coordinates
    """
    n = contact_matrix.shape[0]
    print(f"ShRec3D: {n} bins, alpha={alpha}")

    # Step 1: Contact to distance
    print("  Step 1: Contact-to-distance conversion...")
    dist = contact_to_distance(contact_matrix, alpha)
    n_finite = np.isfinite(dist).sum() - n  # exclude diagonal
    n_total = n * (n - 1)
    print(f"    {n_finite}/{n_total} finite distances "
          f"({100*n_finite/n_total:.1f}%)")

    # Step 2: Shortest-path distance completion
    print("  Step 2: Shortest-path distance completion...")
    dist_complete = shortest_path_completion(dist)

    # Check for disconnected components
    n_inf = np.isinf(dist_complete).sum()
    if n_inf > 0:
        print(f"    WARNING: {n_inf} unreachable pairs remain")
        # Replace remaining inf with max finite distance
        max_d = dist_complete[np.isfinite(dist_complete)].max()
        dist_complete[np.isinf(dist_complete)] = max_d * 1.5

    # Step 3: Classical MDS
    print("  Step 3: Classical MDS (Torgerson)...")
    coords = classical_mds(dist_complete, ndim)

    print(f"  Output: {coords.shape[0]} x {coords.shape[1]} coordinates")
    return coords


def main():
    parser = argparse.ArgumentParser(
        description="ShRec3D: 3D reconstruction from Hi-C contact matrix")
    parser.add_argument("-i", "--input", required=True,
                        help="Input NxN contact matrix (whitespace-delimited)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output coordinate file (x y z per line)")
    parser.add_argument("--alpha", type=float, default=1.0,
                        help="Frequency-to-distance exponent (default: 1.0)")
    parser.add_argument("--ndim", type=int, default=3,
                        help="Output dimensions (default: 3)")
    args = parser.parse_args()

    print(f"Loading contact matrix: {args.input}")
    mat = np.loadtxt(args.input)
    print(f"Matrix shape: {mat.shape}")

    # Ensure symmetric
    mat = (mat + mat.T) / 2.0

    coords = shrec3d(mat, alpha=args.alpha, ndim=args.ndim)

    np.savetxt(args.output, coords, fmt="%.6f", delimiter="\t")
    print(f"Coordinates saved to {args.output}")


if __name__ == "__main__":
    main()
