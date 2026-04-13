/**
 * @file MetricsFramework.hpp
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

#include <Chromosome.hpp>
#include <Heatmap.hpp>

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

  /**
   * @brief Compute Spearman rank correlation between two float vectors.
   * @param a First vector.
   * @param b Second vector (same length).
   * @return Spearman rank correlation coefficient.
   */
  float spearmanCorrelation(const std::vector<float> &a,
                            const std::vector<float> &b) const;

  /**
   * @brief Compute intra-ensemble structural similarity distribution.
   *
   * For each pair of structures in the ensemble, computes Spearman
   * correlation of pairwise distance vectors. Returns all pairwise
   * correlation values.
   *
   * @param structures Vector of Chromosome structures (ensemble members).
   * @param subsample_step Step size for subsampling distance pairs.
   * @return Vector of pairwise Spearman correlation values.
   */
  std::vector<float> ensembleSimilarityDistribution(
      const std::vector<Chromosome> &structures,
      int subsample_step = 1) const;

  /**
   * @brief Write ensemble similarity distribution to CSV file.
   * @param similarities Pairwise similarity values.
   * @param output_path Output file path.
   */
  void writeEnsembleSimilarity(const std::vector<float> &similarities,
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


// ============================================================================
// Implementation
// ============================================================================


#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>

#ifdef _OPENMP
  #include <omp.h>
#endif

inline MetricsFramework::MetricsFramework() {}

inline float MetricsFramework::pearsonCorrelation(const std::vector<float> &a,
                                           const std::vector<float> &b) const {
  if (a.size() != b.size() || a.empty())
    return 0.0f;

  int n = static_cast<int>(a.size());
  double sum_a = 0.0, sum_b = 0.0;
  for (int i = 0; i < n; ++i) {
    sum_a += a[i];
    sum_b += b[i];
  }
  double mean_a = sum_a / n;
  double mean_b = sum_b / n;

  double cov = 0.0, var_a = 0.0, var_b = 0.0;
  for (int i = 0; i < n; ++i) {
    double da = a[i] - mean_a;
    double db = b[i] - mean_b;
    cov += da * db;
    var_a += da * da;
    var_b += db * db;
  }

  double denom = std::sqrt(var_a * var_b);
  if (denom < 1e-12)
    return 0.0f;

  return static_cast<float>(cov / denom);
}

inline float MetricsFramework::distanceCorrelation(const Chromosome &a,
                                            const Chromosome &b,
                                            int subsample_step) const {
  int n = std::min(a.size, b.size);
  if (n < 2)
    return 0.0f;

  int step = std::max(1, subsample_step);

  std::vector<float> dist_a, dist_b;
  for (int i = 0; i < n; i += step) {
    for (int j = i + step; j < n; j += step) {
      dist_a.push_back((a.points[i] - a.points[j]).length());
      dist_b.push_back((b.points[i] - b.points[j]).length());
    }
  }

  return pearsonCorrelation(dist_a, dist_b);
}

inline float MetricsFramework::alignedRMSD(Chromosome &a,
                                     const Chromosome &b) const {
  a.align(b);
  return a.calcRMSD(b);
}

inline std::vector<ContactDecayPoint>
MetricsFramework::computeContactDecay(const Heatmap &heatmap,
                                      int max_separation,
                                      int bin_size) const {
  int n = static_cast<int>(heatmap.size);
  if (max_separation <= 0 || max_separation > n)
    max_separation = n;
  if (bin_size < 1)
    bin_size = 1;

  int num_bins = (max_separation + bin_size - 1) / bin_size;
  std::vector<double> sum(num_bins, 0.0);
  std::vector<int> count(num_bins, 0);

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      int sep = j - i;
      if (sep > max_separation)
        break;
      int bin = (sep - 1) / bin_size;
      if (bin >= 0 && bin < num_bins) {
        sum[bin] += heatmap.v[i][j];
        count[bin]++;
      }
    }
  }

  std::vector<ContactDecayPoint> result;
  for (int b = 0; b < num_bins; ++b) {
    ContactDecayPoint pt;
    pt.genomic_separation = (b * bin_size) + bin_size / 2 + 1;
    pt.mean_contact = (count[b] > 0) ? static_cast<float>(sum[b] / count[b]) : 0.0f;
    pt.count = count[b];
    result.push_back(pt);
  }

  return result;
}

inline void MetricsFramework::writeContactDecay(
    const std::vector<ContactDecayPoint> &decay,
    const std::string &output_path) const {
  FILE *f = fopen(output_path.c_str(), "w");
  if (!f) {
    printf("Error: cannot write %s\n", output_path.c_str());
    return;
  }

  fprintf(f, "genomic_separation,mean_contact,pair_count\n");
  for (const auto &pt : decay) {
    fprintf(f, "%d,%.8f,%d\n", pt.genomic_separation, pt.mean_contact,
            pt.count);
  }
  fclose(f);
  printf("Contact decay curve (%d points) written to %s\n",
         static_cast<int>(decay.size()), output_path.c_str());
}

inline std::vector<LoopEnrichmentResult> MetricsFramework::computeLoopEnrichment(
    const Heatmap &heatmap,
    const std::vector<std::pair<int, int>> &loop_pairs) const {

  int n = static_cast<int>(heatmap.size);

  // Precompute average contact per genomic separation
  std::vector<double> sep_sum(n, 0.0);
  std::vector<int> sep_count(n, 0);

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      int sep = j - i;
      sep_sum[sep] += heatmap.v[i][j];
      sep_count[sep]++;
    }
  }

  std::vector<float> bg_avg(n, 0.0f);
  for (int s = 1; s < n; ++s) {
    if (sep_count[s] > 0)
      bg_avg[s] = static_cast<float>(sep_sum[s] / sep_count[s]);
  }

  std::vector<LoopEnrichmentResult> results;
  for (const auto &lp : loop_pairs) {
    int a = lp.first;
    int b = lp.second;
    if (a < 0 || b < 0 || a >= n || b >= n)
      continue;
    if (a > b)
      std::swap(a, b);

    int sep = b - a;

    LoopEnrichmentResult res;
    res.anchor_a = a;
    res.anchor_b = b;
    res.loop_contact = heatmap.v[a][b];
    res.bg_contact = (sep < n) ? bg_avg[sep] : 0.0f;
    res.enrichment =
        (res.bg_contact > 1e-10f) ? res.loop_contact / res.bg_contact : 0.0f;
    results.push_back(res);
  }

  return results;
}

inline void MetricsFramework::writeLoopEnrichment(
    const std::vector<LoopEnrichmentResult> &results,
    const std::string &output_path) const {
  FILE *f = fopen(output_path.c_str(), "w");
  if (!f) {
    printf("Error: cannot write %s\n", output_path.c_str());
    return;
  }

  fprintf(f, "anchor_a,anchor_b,loop_contact,bg_contact,enrichment\n");
  for (const auto &r : results) {
    fprintf(f, "%d,%d,%.8f,%.8f,%.4f\n", r.anchor_a, r.anchor_b,
            r.loop_contact, r.bg_contact, r.enrichment);
  }
  fclose(f);
  printf("Loop enrichment (%d loops) written to %s\n",
         static_cast<int>(results.size()), output_path.c_str());
}

inline float MetricsFramework::structuralSimilarity(const Heatmap &a,
                                             const Heatmap &b,
                                             int h_param) const {
  int n = static_cast<int>(std::min(a.size, b.size));
  if (n < 2)
    return 0.0f;

  auto diag_corrs = perDiagonalCorrelation(a, b);
  if (diag_corrs.empty())
    return 0.0f;

  // Weight by number of pairs at each diagonal and distance kernel
  double weighted_sum = 0.0;
  double weight_total = 0.0;

  for (const auto &dc : diag_corrs) {
    int k = dc.first;
    float corr = dc.second;

    if (std::isnan(corr))
      continue;

    int num_pairs = n - k;
    if (num_pairs < 2)
      continue;

    // Gaussian distance weighting
    double w = num_pairs * std::exp(-static_cast<double>(k * k) /
                                    (2.0 * h_param * h_param));
    weighted_sum += w * corr;
    weight_total += w;
  }

  if (weight_total < 1e-12)
    return 0.0f;

  return static_cast<float>(weighted_sum / weight_total);
}

inline std::vector<std::pair<int, float>>
MetricsFramework::perDiagonalCorrelation(const Heatmap &a,
                                         const Heatmap &b) const {
  int n = static_cast<int>(std::min(a.size, b.size));
  std::vector<std::pair<int, float>> result;

  for (int k = 1; k < n; ++k) {
    std::vector<float> vals_a, vals_b;
    for (int i = 0; i + k < n; ++i) {
      vals_a.push_back(a.v[i][i + k]);
      vals_b.push_back(b.v[i][i + k]);
    }

    if (vals_a.size() < 2)
      continue;

    float corr = pearsonCorrelation(vals_a, vals_b);
    result.push_back({k, corr});
  }

  return result;
}

inline float MetricsFramework::spearmanCorrelation(const std::vector<float> &a,
                                            const std::vector<float> &b) const {
  if (a.size() != b.size() || a.empty())
    return 0.0f;

  int n = static_cast<int>(a.size());

  // Compute ranks for a
  std::vector<int> idx_a(n), idx_b(n);
  std::iota(idx_a.begin(), idx_a.end(), 0);
  std::iota(idx_b.begin(), idx_b.end(), 0);

  std::sort(idx_a.begin(), idx_a.end(),
            [&](int i, int j) { return a[i] < a[j]; });
  std::sort(idx_b.begin(), idx_b.end(),
            [&](int i, int j) { return b[i] < b[j]; });

  std::vector<float> rank_a(n), rank_b(n);

  // Assign ranks with tie averaging
  auto assign_ranks = [&](const std::vector<int> &idx,
                          const std::vector<float> &vals,
                          std::vector<float> &ranks) {
    int i = 0;
    while (i < n) {
      int j = i;
      while (j < n && vals[idx[j]] == vals[idx[i]])
        ++j;
      float avg_rank = 0.5f * (i + j - 1);
      for (int k = i; k < j; ++k)
        ranks[idx[k]] = avg_rank;
      i = j;
    }
  };

  assign_ranks(idx_a, a, rank_a);
  assign_ranks(idx_b, b, rank_b);

  return pearsonCorrelation(rank_a, rank_b);
}

inline std::vector<float> MetricsFramework::ensembleSimilarityDistribution(
    const std::vector<Chromosome> &structures, int subsample_step) const {
  int n = static_cast<int>(structures.size());
  if (n < 2)
    return {};

  int step = std::max(1, subsample_step);
  int n_pairs = n * (n - 1) / 2;
  std::vector<float> similarities(n_pairs, 0.0f);

#ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
  for (int k = 0; k < n_pairs; ++k) {
    // Map linear index to upper triangle (i, j) with i < j
    int i = 0, cumul = n - 1;
    while (k >= cumul) { i++; cumul += n - 1 - i; }
    int j = k - (cumul - (n - 1 - i)) + i + 1;

    int sz = std::min(structures[i].size, structures[j].size);
    if (sz < 2) continue;

    std::vector<float> dist_a, dist_b;
    for (int a = 0; a < sz; a += step) {
      for (int b = a + step; b < sz; b += step) {
        dist_a.push_back(
            (structures[i].points[a] - structures[i].points[b]).length());
        dist_b.push_back(
            (structures[j].points[a] - structures[j].points[b]).length());
      }
    }

    if (!dist_a.empty())
      similarities[k] = spearmanCorrelation(dist_a, dist_b);
  }

  return similarities;
}

inline void MetricsFramework::writeEnsembleSimilarity(
    const std::vector<float> &similarities,
    const std::string &output_path) const {
  FILE *f = fopen(output_path.c_str(), "w");
  if (!f) {
    printf("Error: cannot write %s\n", output_path.c_str());
    return;
  }

  fprintf(f, "pair_index,spearman_correlation\n");
  for (size_t i = 0; i < similarities.size(); ++i) {
    fprintf(f, "%d,%.8f\n", static_cast<int>(i), similarities[i]);
  }
  fclose(f);

  // Compute summary statistics
  if (!similarities.empty()) {
    float sum = 0.0f;
    float min_val = similarities[0], max_val = similarities[0];
    for (float v : similarities) {
      sum += v;
      if (v < min_val) min_val = v;
      if (v > max_val) max_val = v;
    }
    float mean = sum / similarities.size();
    printf("Ensemble similarity (%d pairs): mean=%.4f min=%.4f max=%.4f -> %s\n",
           static_cast<int>(similarities.size()), mean, min_val, max_val,
           output_path.c_str());
  }
}

inline void MetricsFramework::writeMetricsReport(const Chromosome &chr_a,
                                          const Chromosome &chr_b,
                                          const Heatmap &heatmap_a,
                                          const Heatmap &heatmap_b,
                                          const std::string &output_path) const {
  FILE *f = fopen(output_path.c_str(), "w");
  if (!f) {
    printf("Error: cannot write %s\n", output_path.c_str());
    return;
  }

  // Distance correlation
  float dist_corr = distanceCorrelation(chr_a, chr_b,
                                        std::max(1, chr_a.size / 500));

  // RMSD (make copies for alignment)
  Chromosome a_copy = chr_a;
  Chromosome b_copy = chr_b;
  a_copy.align(b_copy);
  float rmsd = a_copy.calcRMSD(b_copy);

  // Structural similarity
  float scc = structuralSimilarity(heatmap_a, heatmap_b);

  fprintf(f, "# Structural Comparison Metrics\n");
  fprintf(f, "distance_correlation = %.6f\n", dist_corr);
  fprintf(f, "rmsd = %.6f\n", rmsd);
  fprintf(f, "structural_similarity_scc = %.6f\n", scc);
  fprintf(f, "structure_a_beads = %d\n", chr_a.size);
  fprintf(f, "structure_b_beads = %d\n", chr_b.size);
  fprintf(f, "heatmap_a_size = %d\n", static_cast<int>(heatmap_a.size));
  fprintf(f, "heatmap_b_size = %d\n", static_cast<int>(heatmap_b.size));

  fclose(f);
  printf("Metrics report written to %s\n", output_path.c_str());
}

#endif /* METRICSFRAMEWORK_H_ */
