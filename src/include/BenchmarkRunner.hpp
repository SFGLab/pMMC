/**
 * @file BenchmarkRunner.hpp
 * @brief Automated benchmark pipeline: generate synthetic data, reconstruct, compare.
 *
 * Runs the full cycle for multiple loop counts (10, 100, 1000) on a
 * single chromosome, measuring distance correlation, RMSD, and wall-clock time.
 *
 * @see SyntheticGenerator for the data generation step.
 * @see LooperSolver for the reconstruction engine.
 */

#ifndef BENCHMARKRUNNER_H_
#define BENCHMARKRUNNER_H_

#include <string>
#include <vector>

/**
 * @class BenchmarkRunner
 * @brief Generate-reconstruct-compare benchmark for cudaMMC.
 *
 * Usage:
 * @code
 *   BenchmarkRunner bench;
 *   bench.setChromosome("chr22", 51304566);
 *   bench.setResolution(25000);
 *   bench.setEnsembleSize(20);
 *   bench.setOutputDir("./benchmark/");
 *   bench.run();
 * @endcode
 */
class BenchmarkRunner {
public:
  /** @brief Default constructor (chr22, 25kb resolution, 20 ensemble). */
  BenchmarkRunner();

  /**
   * @brief Set the target chromosome and its length.
   * @param chr Chromosome name (e.g. "chr22").
   * @param length_bp Chromosome length in base pairs.
   */
  void setChromosome(const std::string &chr, int length_bp);

  /**
   * @brief Set the reconstruction resolution.
   * @param bp_per_bead Base pairs per bead (e.g. 25000 for 25kb).
   */
  void setResolution(int bp_per_bead);

  /**
   * @brief Set the number of ensemble members per reconstruction.
   * @param n Ensemble size (default 20).
   */
  void setEnsembleSize(int n);

  /**
   * @brief Set the output directory for benchmark results.
   * @param dir Path to the output directory (created if missing).
   */
  void setOutputDir(const std::string &dir);

  /**
   * @brief Set a label prefix for output filenames.
   * @param lbl Label string (e.g. "bench").
   */
  void setLabel(const std::string &lbl);

  /**
   * @brief Run the full benchmark for loop counts {10, 100, 1000}.
   *
   * For each loop count: generates synthetic polymer + ChIA-PET data,
   * runs cudaMMC reconstruction, compares to ground truth, and reports
   * distance correlation, RMSD, and elapsed time.
   */
  void run();

private:
  std::string chr_name;    /**< Target chromosome name. */
  int chr_length_bp;       /**< Chromosome length in bp. */
  int resolution_bp;       /**< Base pairs per bead. */
  int ensemble_size;       /**< Number of ensemble members. */
  std::string output_dir;  /**< Output directory path. */
  std::string label;       /**< Label prefix for filenames. */

  /** @brief Results from a single benchmark run. */
  struct BenchmarkResult {
    int num_loops;           /**< Number of synthetic loops used. */
    float dist_correlation;  /**< Pearson correlation of pairwise distances (ground truth vs reconstructed). */
    float rmsd;              /**< Root-mean-square deviation after Procrustes alignment. */
    double elapsed_seconds;  /**< Wall-clock time for generate + reconstruct. */
  };

  /**
   * @brief Run a single benchmark with the given loop count.
   * @param num_loops Number of synthetic loops to generate.
   * @param subdir Subdirectory for this run's output.
   * @return BenchmarkResult with correlation, RMSD, and timing.
   */
  BenchmarkResult runSingle(int num_loops, const std::string &subdir);

  /**
   * @brief Compute Pearson correlation between two distance vectors.
   * @param a First vector of pairwise distances.
   * @param b Second vector of pairwise distances (same length as a).
   * @return Pearson correlation coefficient in [-1, 1].
   */
  float calcDistanceCorrelation(const std::vector<float> &a,
                                const std::vector<float> &b);

  /**
   * @brief Print a formatted results table to stdout.
   * @param results Vector of BenchmarkResult structs.
   */
  void printResults(const std::vector<BenchmarkResult> &results);

  /**
   * @brief Write results to a CSV file in the output directory.
   * @param results Vector of BenchmarkResult structs.
   */
  void writeResultsCSV(const std::vector<BenchmarkResult> &results);
};


// ============================================================================
// Implementation
// ============================================================================

#include <SvgChartGenerator.hpp>
#include <SyntheticGenerator.hpp>

#include <BedRegion.hpp>
#include <Chromosome.hpp>
#include <Heatmap.hpp>
#include <HierarchicalChromosome.h>
#include <LooperSolver.hpp>
#include <Settings.hpp>
#include <common.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <string>
#include <vector>

#include <platform.h>

inline BenchmarkRunner::BenchmarkRunner() {
  chr_name = "chr22";
  chr_length_bp = 51304566;
  resolution_bp = 25000;
  ensemble_size = 20;
  output_dir = "./benchmark/";
  label = "bench";
}

inline void BenchmarkRunner::setChromosome(const std::string &chr, int length_bp) {
  chr_name = chr;
  chr_length_bp = length_bp;
}

inline void BenchmarkRunner::setResolution(int bp_per_bead) {
  resolution_bp = bp_per_bead;
}

inline void BenchmarkRunner::setEnsembleSize(int n) { ensemble_size = n; }

inline void BenchmarkRunner::setOutputDir(const std::string &dir) {
  output_dir = dir;
}

inline void BenchmarkRunner::setLabel(const std::string &lbl) { label = lbl; }

inline void BenchmarkRunner::run() {
  portable_mkdir(output_dir.c_str());

  printf("=== Benchmark: %s (%d bp, %d bp/bead, ensemble=%d) ===\n",
         chr_name.c_str(), chr_length_bp, resolution_bp, ensemble_size);

  std::vector<int> loop_counts = {10, 100, 1000};
  std::vector<BenchmarkResult> results;

  for (int nl : loop_counts) {
    std::string subdir =
        ftext("%sloops_%d/", output_dir.c_str(), nl);
    portable_mkdir(subdir.c_str());

    BenchmarkResult res = runSingle(nl, subdir);
    results.push_back(res);

    printf("\n--- loops=%d: dist_corr=%.3f, rmsd=%.2f, time=%.1fs ---\n\n",
           res.num_loops, res.dist_correlation, res.rmsd, res.elapsed_seconds);
  }

  printResults(results);
  writeResultsCSV(results);
}

inline BenchmarkRunner::BenchmarkResult
BenchmarkRunner::runSingle(int num_loops, const std::string &subdir) {
  BenchmarkResult result;
  result.num_loops = num_loops;
  result.dist_correlation = 0.0f;
  result.rmsd = 0.0f;
  result.elapsed_seconds = 0.0;

  auto t_start = std::chrono::high_resolution_clock::now();

  // --- Step 1: Generate synthetic data ---
  printf("\n========== BENCHMARK: %d loops ==========\n", num_loops);
  printf("Step 1: Generate synthetic data\n");

  SyntheticGenerator gen;
  gen.setChromosome(chr_name, chr_length_bp);
  gen.setResolution(resolution_bp);
  gen.setNumLoops(num_loops);
  gen.setEnsembleSize(ensemble_size);
  gen.setOutputDir(subdir);
  gen.setLabel(ftext("%s_%d", label.c_str(), num_loops));
  gen.generate();

  // Save ground truth reference for comparison
  const Chromosome &ground_truth = gen.getGroundTruth();

  // --- Step 2: Load settings and run reconstruction ---
  printf("Step 2: Reconstruct from generated data\n");

  std::string settings_path = subdir + "settings.ini";
  Settings stg;
  stg.loadFromINI(settings_path);

  // Prepare chromosome list (single chr)
  std::vector<std::string> chrs;
  chrs.push_back(chr_name);
  BedRegion region_of_interest;

  // Build file paths from settings (following runLooper pattern)
  std::string path_data = Settings::dataDirectory;
  std::string path_anchors = Settings::dataAnchors;
  std::vector<std::string> path_pets = split(Settings::dataPetClusters);
  std::vector<std::string> path_singletons = split(Settings::dataSingletons);
  std::vector<std::string> path_singletons_inter =
      split(Settings::dataSingletonsInter);
  std::vector<std::string> factors = split(Settings::dataFactors);

  std::string anchors =
      ftext("%s%s", path_data.c_str(), path_anchors.c_str());

  std::vector<std::string> arcs_clusters;
  std::vector<std::string> arcs_singletons;
  std::vector<std::string> arcs_singletons_inter;

  for (size_t i = 0; i < path_pets.size(); ++i)
    arcs_clusters.push_back(
        ftext("%s%s", path_data.c_str(), path_pets[i].c_str()));

  for (size_t i = 0; i < path_singletons.size(); ++i)
    arcs_singletons.push_back(
        ftext("%s%s", path_data.c_str(), path_singletons[i].c_str()));

  for (size_t i = 0; i < path_singletons_inter.size(); ++i)
    arcs_singletons_inter.push_back(
        ftext("%s%s", path_data.c_str(), path_singletons_inter[i].c_str()));

  std::string recon_outdir = subdir + "output/";
  portable_mkdir(recon_outdir.c_str());

  std::string recon_label =
      ftext("%s_%d_recon", label.c_str(), num_loops);

  LooperSolver lsm(recon_label, recon_outdir);
  lsm.setContactData(chrs, region_of_interest, anchors, factors, arcs_clusters,
                     arcs_singletons, arcs_singletons_inter);

  printf("  Creating tree...\n");
  lsm.createTreeGenome();

  printf("  Reconstructing heatmap levels...\n");
  lsm.reconstructClustersHeatmap();

  printf("  Reconstructing arc distances...\n");
  lsm.reconstructClustersArcsDistances();

  printf("  Extracting model...\n");
  HierarchicalChromosome hc = lsm.getModel();
  hc.toFile(ftext("%sresult_%d.hcm", recon_outdir.c_str(), num_loops));

  auto t_end = std::chrono::high_resolution_clock::now();
  result.elapsed_seconds =
      std::chrono::duration<double>(t_end - t_start).count();

  // --- Step 3: Compare reconstructed to ground truth ---
  printf("Step 3: Compare to ground truth\n");

  // Extract flat chromosome from hierarchical model at subanchor level
  hc.setLevel(LVL_SUBANCHOR);
  Chromosome reconstructed = hc.createEqidistantModel(resolution_bp, chr_name);

  if (reconstructed.size < 2 || ground_truth.size < 2) {
    printf("  Warning: insufficient points for comparison (recon=%d, gt=%d)\n",
           reconstructed.size, (int)ground_truth.points.size());
    return result;
  }

  // Align reconstructed to ground truth
  Chromosome gt_copy = ground_truth;
  reconstructed.align(gt_copy);
  result.rmsd = reconstructed.calcRMSD(gt_copy);

  // Distance correlation: compare pairwise distance matrices
  // Use a subsample for large chromosomes
  int n_gt = (int)gt_copy.points.size();
  int n_recon = reconstructed.size;
  int n_compare = std::min(n_gt, n_recon);
  int step = 1;
  if (n_compare > 500)
    step = n_compare / 500;

  std::vector<float> gt_dists, recon_dists;
  for (int i = 0; i < n_compare; i += step) {
    for (int j = i + step; j < n_compare; j += step) {
      float d_gt = (gt_copy.points[i] - gt_copy.points[j]).length();
      float d_re =
          (reconstructed.points[i] - reconstructed.points[j]).length();
      gt_dists.push_back(d_gt);
      recon_dists.push_back(d_re);
    }
  }

  result.dist_correlation = calcDistanceCorrelation(gt_dists, recon_dists);

  printf("  Points compared: %d (step=%d)\n", n_compare, step);
  printf("  Distance pairs: %d\n", (int)gt_dists.size());
  printf("  Distance correlation: %.4f\n", result.dist_correlation);
  printf("  RMSD: %.4f\n", result.rmsd);

  return result;
}

inline float BenchmarkRunner::calcDistanceCorrelation(const std::vector<float> &a,
                                               const std::vector<float> &b) {
  if (a.size() != b.size() || a.empty())
    return 0.0f;

  int n = (int)a.size();
  double sum_a = 0, sum_b = 0;
  for (int i = 0; i < n; ++i) {
    sum_a += a[i];
    sum_b += b[i];
  }
  double mean_a = sum_a / n;
  double mean_b = sum_b / n;

  double cov = 0, var_a = 0, var_b = 0;
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

  return (float)(cov / denom);
}

inline void BenchmarkRunner::printResults(
    const std::vector<BenchmarkResult> &results) {
  printf("\n");
  printf("========================================\n");
  printf("        BENCHMARK RESULTS\n");
  printf("========================================\n");
  printf("Loops | Dist.Corr | RMSD    | Time(s)\n");
  printf("------|-----------|---------|--------\n");
  for (const auto &r : results) {
    printf("%-5d | %9.3f | %7.2f | %6.1f\n", r.num_loops, r.dist_correlation,
           r.rmsd, r.elapsed_seconds);
  }
  printf("========================================\n");
}

inline void BenchmarkRunner::writeResultsCSV(
    const std::vector<BenchmarkResult> &results) {
  std::string path = output_dir + "benchmark_results.csv";
  FILE *f = fopen(path.c_str(), "w");
  if (!f) {
    printf("Error: cannot write %s\n", path.c_str());
    return;
  }

  fprintf(f, "loops,dist_correlation,rmsd,elapsed_seconds\n");
  for (const auto &r : results) {
    fprintf(f, "%d,%.6f,%.6f,%.3f\n", r.num_loops, r.dist_correlation, r.rmsd,
            r.elapsed_seconds);
  }
  fclose(f);
  printf("Results saved to %s\n", path.c_str());

  // Generate SVG chart
  std::string svg_path = output_dir + "benchmark_results.svg";
  SvgChartGenerator::benchmarkResultsChart(path, svg_path);
}

#endif /* BENCHMARKRUNNER_H_ */
