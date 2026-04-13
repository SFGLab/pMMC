/**
 * @file MdsInitializer.hpp
 * @brief Classical Multidimensional Scaling (MDS) initializer for 3D bead positions.
 *
 * Implements Torgerson's classical MDS (cMDS) algorithm to convert a pairwise
 * distance matrix (such as the heatmap_dist from LooperSolver) into a set of
 * initial 3D Cartesian coordinates.  The result can seed Monte Carlo
 * simulations with a geometrically meaningful starting structure rather than
 * a random walk, reducing the number of MC steps required to converge.
 *
 * Algorithm overview:
 *  1. Symmetrize and square the input distance matrix D to form D2.
 *  2. Double-centre D2 to produce the Gram matrix B.
 *  3. Compute the top 3 eigenpairs of B via the Jacobi eigenvalue algorithm.
 *  4. Return coordinates X[i][k] = sqrt(max(0, lambda_k)) * V[i][k].
 *
 * @see LooperSolver::heatmap_dist for the expected-distance heatmap.
 */

#ifndef MDSINITIALIZER_H_
#define MDSINITIALIZER_H_

#include <vector>

#include <Heatmap.hpp>
#include <common.h>   // provides vector3 via mtxlib.h

/**
 * @class MdsInitializer
 * @brief Converts a pairwise distance matrix into initial 3D bead coordinates
 *        using classical (Torgerson) Multidimensional Scaling.
 *
 * All methods are static: the class acts as a namespaced collection of
 * functions and requires no instance state.
 *
 * Example usage:
 * @code
 *   // From a Heatmap distance matrix
 *   std::vector<vector3> coords =
 *       MdsInitializer::computeCoordinates(solver.heatmap_dist);
 *   if (!coords.empty()) {
 *       // Use coords as initial positions for MC
 *   }
 * @endcode
 */
class MdsInitializer {
public:
    /**
     * @brief Compute 3D coordinates from a Heatmap distance matrix.
     *
     * Reads the N×N matrix stored in @p distanceMatrix and applies classical
     * MDS to produce N bead positions in 3D space.
     *
     * @param distanceMatrix  Square symmetric matrix of expected spatial
     *                        distances between beads (v[i][j] = d_ij).
     * @return Vector of N vector3 positions, one per bead, or an empty vector
     *         if MDS fails (N < 3, all-zero input, or eigendecomposition did
     *         not converge).
     */
    static std::vector<vector3> computeCoordinates(const Heatmap& distanceMatrix);

    /**
     * @brief Compute 3D coordinates from a raw float** distance matrix.
     *
     * Convenience overload for callers that already hold the distance data as
     * a plain 2D C array rather than a Heatmap object.
     *
     * @param distMatrix  Row-pointer array: distMatrix[i][j] = d_ij.
     *                    Must have at least @p n rows, each of length @p n.
     * @param n           Matrix dimension (number of beads).
     * @return Vector of N vector3 positions, or empty on failure.
     */
    static std::vector<vector3> computeCoordinates(float** distMatrix, int n);

private:
    /**
     * @brief Jacobi eigenvalue decomposition for a real symmetric matrix.
     *
     * Iteratively applies Givens (Jacobi) rotations until all off-diagonal
     * elements are below a convergence threshold.  The algorithm is simple,
     * robust, and free of external dependencies — appropriate for the moderate
     * matrix sizes (~10–1000) encountered in this application.
     *
     * @param[in,out] matrix      Flat row-major representation of the n×n
     *                            symmetric input matrix.  The matrix is
     *                            destroyed (overwritten) during iteration.
     * @param[in]     n           Matrix dimension.
     * @param[out]    eigenvalues Vector of n eigenvalues, sorted in descending
     *                            order on successful return.
     * @param[out]    eigenvectors Flat column-major representation of the n×n
     *                            eigenvector matrix: column k holds the
     *                            eigenvector for eigenvalue k.
     * @param[in]     maxIter     Maximum number of Jacobi sweep iterations
     *                            (default 100; each sweep visits all n*(n-1)/2
     *                            off-diagonal pairs once).
     * @return @c true if the algorithm converged within @p maxIter sweeps;
     *         @c false otherwise (output may be partially correct).
     */
    static bool jacobiEigen(std::vector<double>& matrix, int n,
                            std::vector<double>& eigenvalues,
                            std::vector<double>& eigenvectors,
                            int maxIter = 100);
};


// ============================================================================
// Implementation
// ============================================================================

/**
 * @file MdsInitializer.cpp
 * @brief Implementation of classical MDS initializer using the Jacobi
 *        eigenvalue algorithm.
 *
 * Implements Torgerson's classical Multidimensional Scaling (cMDS):
 *
 *   Given an N×N distance matrix D:
 *     1. Symmetrize:     D_sym[i][j] = 0.5*(D[i][j] + D[j][i])
 *     2. Square:         D2[i][j]    = D_sym[i][j]^2
 *     3. Double-centre:  B[i][j]     = -0.5*(D2[i][j]
 *                                            - row_mean[i]
 *                                            - col_mean[j]
 *                                            + grand_mean)
 *     4. Eigenpairs of B (top 3 by magnitude):
 *                        B = V * Lambda * V^T
 *     5. Coordinates:    X[i][k] = sqrt(max(0, lambda_k)) * V[i][k]
 *
 * All internal arithmetic uses double precision to limit accumulation errors.
 *
 * Output coordinates are normalised so that the mean pairwise Euclidean
 * distance across all bead pairs equals 1.0 (unit-distance scale).  If the
 * scale is zero (all beads are coincident) normalisation is skipped.
 */


#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <vector>

/* =========================================================================
 * Internal helpers — file-local, not part of the public interface
 * ========================================================================= */

namespace {

/**
 * @brief Inline accessor for a flat row-major matrix stored in a std::vector.
 * @param M   Flat storage (length n*n).
 * @param n   Matrix dimension.
 * @param r   Row index.
 * @param c   Column index.
 * @return    Reference to M[r][c].
 */
static inline double& mat(std::vector<double>& M, int n, int r, int c) {
    return M[r * n + c];
}

static inline double mat(const std::vector<double>& M, int n, int r, int c) {
    return M[r * n + c];
}

/**
 * @brief Compute the Frobenius norm of all off-diagonal elements of a
 *        symmetric matrix stored flat row-major.
 */
static double offDiagonalNorm(const std::vector<double>& M, int n) {
    double s = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j)
                s += M[i * n + j] * M[i * n + j];
    return std::sqrt(s);
}

} // anonymous namespace

/* =========================================================================
 * MdsInitializer::jacobiEigen
 * ========================================================================= */

/**
 * Classical Jacobi algorithm for symmetric eigenvalue problems.
 *
 * Each sweep visits every off-diagonal pair (p,q) with p < q in row-major
 * order.  For each pair it computes the 2×2 Jacobi rotation angle theta that
 * annihilates M[p][q], then applies the rotation to the full n×n matrix and
 * accumulates it in the eigenvector matrix V.
 *
 * Convergence criterion: off-diagonal Frobenius norm < n * n * 1e-12.
 *
 * After convergence eigenvalues and eigenvectors are sorted in descending
 * eigenvalue order.
 */
inline bool MdsInitializer::jacobiEigen(std::vector<double>& M, int n,
                                 std::vector<double>& eigenvalues,
                                 std::vector<double>& eigenvectors,
                                 int maxIter)
{
    // Initialise eigenvector matrix to the identity.
    eigenvectors.assign(static_cast<size_t>(n * n), 0.0);
    for (int i = 0; i < n; ++i)
        mat(eigenvectors, n, i, i) = 1.0;

    const double tol = static_cast<double>(n) * static_cast<double>(n) * 1e-12;

    for (int iter = 0; iter < maxIter; ++iter) {
        if (offDiagonalNorm(M, n) < tol)
            break;

        // One sweep: visit all (p,q) pairs with p < q
        for (int p = 0; p < n - 1; ++p) {
            for (int q = p + 1; q < n; ++q) {
                double Mpq = mat(M, n, p, q);
                if (std::fabs(Mpq) < 1e-15)
                    continue;   // already zero — skip

                // Compute rotation angle
                double Mpp = mat(M, n, p, p);
                double Mqq = mat(M, n, q, q);
                double theta = 0.5 * std::atan2(2.0 * Mpq, Mpp - Mqq);

                double c = std::cos(theta);
                double s = std::sin(theta);

                // Apply Givens rotation: M <- G^T * M * G
                // Update the p-th and q-th rows/columns of M.
                // First pass: update columns p and q for all rows r.
                for (int r = 0; r < n; ++r) {
                    double Mrp = mat(M, n, r, p);
                    double Mrq = mat(M, n, r, q);
                    mat(M, n, r, p) =  c * Mrp + s * Mrq;
                    mat(M, n, r, q) = -s * Mrp + c * Mrq;
                }

                // Second pass: update rows p and q for all columns c_idx.
                for (int c_idx = 0; c_idx < n; ++c_idx) {
                    double Mpc = mat(M, n, p, c_idx);
                    double Mqc = mat(M, n, q, c_idx);
                    mat(M, n, p, c_idx) =  c * Mpc + s * Mqc;
                    mat(M, n, q, c_idx) = -s * Mpc + c * Mqc;
                }

                // Accumulate rotation in eigenvector matrix (column updates)
                for (int r = 0; r < n; ++r) {
                    double Vrp = mat(eigenvectors, n, r, p);
                    double Vrq = mat(eigenvectors, n, r, q);
                    mat(eigenvectors, n, r, p) =  c * Vrp + s * Vrq;
                    mat(eigenvectors, n, r, q) = -s * Vrp + c * Vrq;
                }
            }
        }
    }

    // Check final convergence
    bool converged = (offDiagonalNorm(M, n) < tol * 1e3);

    // Extract eigenvalues from the diagonal of the (now nearly diagonal) M
    eigenvalues.resize(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
        eigenvalues[static_cast<size_t>(i)] = mat(M, n, i, i);

    // Sort eigenvalues and eigenvectors in descending order
    // Build an index array and sort by eigenvalue descending
    std::vector<int> idx(static_cast<size_t>(n));
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) {
        return eigenvalues[static_cast<size_t>(a)] >
               eigenvalues[static_cast<size_t>(b)];
    });

    std::vector<double> sortedEvals(static_cast<size_t>(n));
    std::vector<double> sortedEvecs(static_cast<size_t>(n * n));

    for (int k = 0; k < n; ++k) {
        sortedEvals[static_cast<size_t>(k)] =
            eigenvalues[static_cast<size_t>(idx[static_cast<size_t>(k)])];
        for (int r = 0; r < n; ++r) {
            mat(sortedEvecs, n, r, k) =
                mat(eigenvectors, n, r, idx[static_cast<size_t>(k)]);
        }
    }

    eigenvalues  = std::move(sortedEvals);
    eigenvectors = std::move(sortedEvecs);

    return converged;
}

/* =========================================================================
 * MdsInitializer::computeCoordinates (raw float** overload)
 * ========================================================================= */

inline std::vector<vector3> MdsInitializer::computeCoordinates(float** distMatrix,
                                                         int n)
{
    // -----------------------------------------------------------------------
    // Edge-case: need at least 3 beads for a 3D embedding
    // -----------------------------------------------------------------------
    if (n < 3 || distMatrix == nullptr) {
        fprintf(stderr, "[MdsInitializer] N=%d < 3 or null matrix — "
                        "skipping MDS.\n", n);
        return {};
    }

    // -----------------------------------------------------------------------
    // Step 1: Symmetrize and square the distance matrix
    //         D2[i][j] = (0.5*(D[i][j] + D[j][i]))^2
    // -----------------------------------------------------------------------
    std::vector<double> D2(static_cast<size_t>(n * n), 0.0);

    bool allZero = true;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double dij = 0.5 * (static_cast<double>(distMatrix[i][j]) +
                                static_cast<double>(distMatrix[j][i]));
            double d2  = dij * dij;
            mat(D2, n, i, j) = d2;
            if (d2 > 1e-12)
                allZero = false;
        }
    }

    if (allZero) {
        fprintf(stderr, "[MdsInitializer] Distance matrix is all-zero — "
                        "skipping MDS.\n");
        return {};
    }

    // -----------------------------------------------------------------------
    // Step 2: Double-centre to obtain the Gram matrix B
    //
    //   row_mean[i]  = (1/N) * sum_j D2[i][j]
    //   col_mean[j]  = (1/N) * sum_i D2[i][j]  (= row_mean[j] for symmetric)
    //   grand_mean   = (1/N^2) * sum_ij D2[i][j]
    //
    //   B[i][j] = -0.5 * (D2[i][j] - row_mean[i] - col_mean[j] + grand_mean)
    // -----------------------------------------------------------------------
    std::vector<double> rowMean(static_cast<size_t>(n), 0.0);
    std::vector<double> colMean(static_cast<size_t>(n), 0.0);
    double grandMean = 0.0;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double v = mat(D2, n, i, j);
            rowMean[static_cast<size_t>(i)] += v;
            colMean[static_cast<size_t>(j)] += v;
            grandMean += v;
        }
    }

    double inv_n = 1.0 / static_cast<double>(n);
    for (int i = 0; i < n; ++i) {
        rowMean[static_cast<size_t>(i)] *= inv_n;
        colMean[static_cast<size_t>(i)] *= inv_n;
    }
    grandMean *= inv_n * inv_n;

    // Reuse D2 storage for B (avoids an extra allocation)
    std::vector<double> B(static_cast<size_t>(n * n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            mat(B, n, i, j) = -0.5 * (mat(D2, n, i, j)
                                       - rowMean[static_cast<size_t>(i)]
                                       - colMean[static_cast<size_t>(j)]
                                       + grandMean);
        }
    }
    D2.clear();

    // -----------------------------------------------------------------------
    // Step 3: Jacobi eigendecomposition of B
    //         We need only the top 3 eigenpairs but the Jacobi algorithm
    //         computes all n — for the sizes used here (up to ~1000) this is
    //         still fast enough.
    // -----------------------------------------------------------------------
    std::vector<double> eigenvalues;
    std::vector<double> eigenvectors;

    // jacobiEigen destroys B (uses it as working storage)
    bool ok = jacobiEigen(B, n, eigenvalues, eigenvectors, 200);
    if (!ok) {
        fprintf(stderr, "[MdsInitializer] Jacobi eigendecomposition did not "
                        "fully converge (N=%d). Using best available result.\n",
                n);
        // Continue with the unconverged result — it is often still usable.
    }

    // -----------------------------------------------------------------------
    // Step 4: Compute 3D coordinates
    //         X[i][k] = sqrt(max(0, lambda_k)) * V[i][k],  k = 0,1,2
    //
    // Negative eigenvalues arise from noise / non-Euclidean distances.
    // They are clamped to zero, producing zero contribution in that dimension.
    // -----------------------------------------------------------------------
    // Number of usable positive eigenvalues (capped at 3)
    int dims = std::min(3, n);

    std::vector<double> scale(static_cast<size_t>(dims));
    for (int k = 0; k < dims; ++k) {
        double lam = eigenvalues[static_cast<size_t>(k)];
        scale[static_cast<size_t>(k)] = (lam > 0.0) ? std::sqrt(lam) : 0.0;
    }

    std::vector<vector3> coords(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        double cx = (dims > 0) ? scale[0] * mat(eigenvectors, n, i, 0) : 0.0;
        double cy = (dims > 1) ? scale[1] * mat(eigenvectors, n, i, 1) : 0.0;
        double cz = (dims > 2) ? scale[2] * mat(eigenvectors, n, i, 2) : 0.0;
        coords[static_cast<size_t>(i)] = vector3(
            static_cast<float>(cx),
            static_cast<float>(cy),
            static_cast<float>(cz));
    }

    // -----------------------------------------------------------------------
    // Step 5: Normalise so that the mean pairwise distance equals 1.0
    //
    // For large N computing all N*(N-1)/2 pairs is expensive.  We subsample
    // up to 10 000 pairs chosen by striding through the upper triangle.
    // -----------------------------------------------------------------------
    if (n >= 2) {
        double distSum  = 0.0;
        long long count = 0;

        // Stride chosen so we sample at most ~10000 pairs
        int stride = std::max(1, n / 150);

        for (int i = 0; i < n; i += stride) {
            for (int j = i + 1; j < n; j += stride) {
                float dx = coords[static_cast<size_t>(i)].x
                           - coords[static_cast<size_t>(j)].x;
                float dy = coords[static_cast<size_t>(i)].y
                           - coords[static_cast<size_t>(j)].y;
                float dz = coords[static_cast<size_t>(i)].z
                           - coords[static_cast<size_t>(j)].z;
                distSum += std::sqrt(static_cast<double>(dx*dx + dy*dy + dz*dz));
                ++count;
            }
        }

        if (count > 0 && distSum > 1e-12) {
            double meanDist = distSum / static_cast<double>(count);
            float  invScale = static_cast<float>(1.0 / meanDist);
            for (auto& pt : coords) {
                pt.x *= invScale;
                pt.y *= invScale;
                pt.z *= invScale;
            }
        }
        // If distSum == 0 all beads are coincident — leave as-is.
    }

    return coords;
}

/* =========================================================================
 * MdsInitializer::computeCoordinates (Heatmap overload)
 * ========================================================================= */

inline std::vector<vector3> MdsInitializer::computeCoordinates(
    const Heatmap& distanceMatrix)
{
    int n = static_cast<int>(distanceMatrix.size);
    return computeCoordinates(distanceMatrix.v, n);
}

#endif /* MDSINITIALIZER_H_ */
