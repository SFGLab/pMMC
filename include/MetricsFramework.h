/**
 * @file MetricsFramework.h
 * @brief Structural comparison metrics for 3D genome models.
 *
 * Provides:
 *   - Pearson distance correlation between structure pairs
 *   - Contact decay curve P(s) as a function of genomic separation
 *   - Loop enrichment analysis (contact frequency at loop loci vs background)
 *   - Structural similarity scoring (SCC-like)
 *   - RMSD after Procrustes alignment
 *
 * All metrics are designed to work with Chromosome and Heatmap objects
 * from the cudaMMC framework.
 *
 * @see Chromosome for 3D structure representation.
 * @see Heatmap for contact/distance matrices.
 */

#ifndef METRICSFRAMEWORK_H_
#define METRICSFRAMEWORK_H_

#include <string>
#include <utility>
#include <vector>

#include <Chromosome.h>
#include <Heatmap.h>

/**
 * @struct ContactDecayPoint
 * @brief A single point on the contact decay curve P(s).
 */
struct ContactDecayPoint {
  int genomic_separation;   /**< Genomic distance in beads. */
  float mean_contact;       /**< Average contact frequency at this separation. */
  int count;                /**< Number of pairs contributing to this average. */
};

/**
 * @struct LoopEnrichmentResult
 * @brief Loop enrichment analysis result for a single loop.
 */
struct LoopEnrichmentResult {
  int anchor_a;        /**< First anchor bead index. */
  int anchor_b;        /**< Second anchor bead index. */
  float loop_contact;  /**< Contact frequency at the loop locus. */
  float bg_contact;    /**< Background contact frequency at same genomic separation. */
  float enrichment;    /**< Enrichment ratio: loop_contact / bg_contact. */
};

/**
 * @class MetricsFramework
 * @brief Computes structural comparison and quality metrics.
 *
 * Usage:
 * @code
 *   MetricsFramework metrics;
 *
 *   // Contact decay curve
 *   auto decay = metrics.computeContactDecay(heatmap, max_sep, bin_size);
 *   metrics.writeContactDecay(decay, "output/contact_decay.csv");
 *
 *   // Loop enrichment
 *   auto enrichment = metrics.computeLoopEnrichment(heatmap, loop_pairs);
 *   metrics.writeLoopEnrichment(enrichment, "output/loop_enrichment.csv");
 *
 *   // Distance correlation
 *   float corr = metrics.distanceCorrelation(chr_a, chr_b);
 *
 *   // Structural similarity
 *   float scc = metrics.structuralSimilarity(heatmap_a, heatmap_b, h_param);
 * @endcode
 */
class MetricsFramework {
public:
  MetricsFramework();

  /**
   * @brief Compute Pearson correlation of pairwise distance vectors.
   * @param a First polymer structure.
   * @param b Second polymer structure (same number of beads).
   * @param subsample_step Step size for subsampling (1 = all pairs).
   * @return Pearson correlation coefficient in [-1, 1].
   */
  float distanceCorrelation(const Chromosome &a, const Chromosome &b,
                            int subsample_step = 1) const;

  /**
   * @brief Compute RMSD after Procrustes alignment.
   * @param a First structure (will be aligned in-place).
   * @param b Second structure (reference).
   * @return RMSD value.
   */
  float alignedRMSD(Chromosome &a, const Chromosome &b) const;

  /**
   * @brief Compute the contact decay curve P(s).
   *
   * For each genomic separation s (in beads), computes the average
   * contact frequency across all pairs at that separation.
   *
   * @param heatmap Contact frequency matrix.
   * @param max_separation Maximum genomic separation to compute.
   * @param bin_size Number of bead separations to bin together.
   * @return Vector of ContactDecayPoint entries.
   */
  std::vector<ContactDecayPoint>
  computeContactDecay(const Heatmap &heatmap, int max_separation,
                      int bin_size = 1) const;

  /**
   * @brief Write contact decay curve to CSV file.
   * @param decay Decay curve data.
   * @param output_path Output file path.
   */
  void writeContactDecay(const std::vector<ContactDecayPoint> &decay,
                         const std::string &output_path) const;

  /**
   * @brief Compute loop enrichment for a set of loop pairs.
   *
   * For each loop pair (a, b), computes the contact frequency at
   * the loop locus and compares to the average background frequency
   * at the same genomic separation.
   *
   * @param heatmap Contact frequency matrix.
   * @param loop_pairs Vector of (anchor_a, anchor_b) bead index pairs.
   * @return Vector of LoopEnrichmentResult entries.
   */
  std::vector<LoopEnrichmentResult> computeLoopEnrichment(
      const Heatmap &heatmap,
      const std::vector<std::pair<int, int>> &loop_pairs) const;

  /**
   * @brief Write loop enrichment results to CSV file.
   * @param results Enrichment data.
   * @param output_path Output file path.
   */
  void
  writeLoopEnrichment(const std::vector<LoopEnrichmentResult> &results,
                      const std::string &output_path) const;

  /**
   * @brief Compute stratum-adjusted correlation coefficient (SCC).
   *
   * Compares two heatmaps by computing per-diagonal Pearson correlations
   * and averaging them with distance-weighted contributions.
   *
   * @param a First heatmap.
   * @param b Second heatmap (same size as a).
   * @param h_param Smoothing parameter for distance weighting.
   * @return SCC value in [-1, 1].
   */
  float structuralSimilarity(const Heatmap &a, const Heatmap &b,
                             int h_param = 5) const;

  /**
   * @brief Compute per-diagonal Pearson correlations between two heatmaps.
   * @param a First heatmap.
   * @param b Second heatmap.
   * @return Vector of (diagonal_offset, correlation) pairs.
   */
  std::vector<std::pair<int, float>>
  perDiagonalCorrelation(const Heatmap &a, const Heatmap &b) const;

  /**
   * @brief Write a full metrics report comparing two structures.
   * @param chr_a First structure.
   * @param chr_b Second structure.
   * @param heatmap_a Contact map from first structure.
   * @param heatmap_b Contact map from second structure.
   * @param output_path Output file path.
   */
  void writeMetricsReport(const Chromosome &chr_a, const Chromosome &chr_b,
                          const Heatmap &heatmap_a, const Heatmap &heatmap_b,
                          const std::string &output_path) const;

private:
  /**
   * @brief Compute Pearson correlation between two float vectors.
   * @param a First vector.
   * @param b Second vector (same length).
   * @return Pearson correlation coefficient.
   */
  float pearsonCorrelation(const std::vector<float> &a,
                           const std::vector<float> &b) const;
};

#endif /* METRICSFRAMEWORK_H_ */
