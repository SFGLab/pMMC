/**
 * @file BenchmarkRunner.h
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

#endif /* BENCHMARKRUNNER_H_ */
