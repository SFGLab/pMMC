/**
 * @file MdsInitializer.h
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

#include <Heatmap.h>
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

#endif /* MDSINITIALIZER_H_ */
