/**
 * @file SyntheticGenerator.h
 * @brief Generate synthetic polymer + ChIA-PET data for benchmarking.
 *
 * Creates a random walk polymer as ground truth, places CTCF-like anchors,
 * generates loop pairs, and writes all files needed by cudaMMC (anchors,
 * PET clusters, singletons, segment splits, centromeres, settings.ini).
 *
 * @see BenchmarkRunner for the full generate-reconstruct-compare pipeline.
 */

#ifndef SYNTHETICGENERATOR_H_
#define SYNTHETICGENERATOR_H_

#include <string>
#include <utility>
#include <vector>

#include <Chromosome.h>
#include <Heatmap.h>

/**
 * @class SyntheticGenerator
 * @brief Generates synthetic chromatin data for controlled benchmarking.
 *
 * Usage:
 * @code
 *   SyntheticGenerator gen;
 *   gen.setChromosome("chr22", 51304566);
 *   gen.setResolution(25000);
 *   gen.setNumLoops(100);
 *   gen.setOutputDir("./synth/");
 *   gen.generate();
 *   const Chromosome& gt = gen.getGroundTruth();
 * @endcode
 *
 * The generated polymer uses fractal-globule-like random walk with
 * harmonic loop constraints enforced by energy relaxation.
 */
class SyntheticGenerator {
public:
  /** @brief Default constructor (chr22, 25kb, 100 loops, 20 ensemble). */
  SyntheticGenerator();

  /**
   * @brief Set the target chromosome and its length.
   * @param chr Chromosome name (e.g. "chr22").
   * @param length_bp Chromosome length in base pairs.
   */
  void setChromosome(const std::string &chr, int length_bp);

  /**
   * @brief Set the reconstruction resolution.
   * @param bp_per_bead Base pairs per bead.
   */
  void setResolution(int bp_per_bead);

  /**
   * @brief Set the number of chromatin loops to generate.
   * @param n Number of loop pairs.
   */
  void setNumLoops(int n);

  /**
   * @brief Set the ensemble size written to the settings file.
   * @param n Number of ensemble members.
   */
  void setEnsembleSize(int n);

  /**
   * @brief Set the segment size for predefined splits.
   * @param bp Segment size in base pairs (default ~2Mb).
   */
  void setSegmentSize(int bp);

  /**
   * @brief Set the output directory for all generated files.
   * @param dir Directory path (created if missing).
   */
  void setOutputDir(const std::string &dir);

  /**
   * @brief Set a label prefix for filenames.
   * @param lbl Label string.
   */
  void setLabel(const std::string &lbl);

  /**
   * @brief Run the full synthetic data generation pipeline.
   *
   * Steps: place anchors, generate loops, create polymer ensemble,
   * write anchors, clusters, singletons, segment split, centromeres,
   * config INI, and ground truth files.
   */
  void generate();

  /**
   * @brief Get the ground truth polymer structure.
   * @return Const reference to the generated Chromosome.
   */
  const Chromosome &getGroundTruth() const;

  /**
   * @brief Get the generated contact map (ensemble-averaged).
   * @return Const reference to the contact Heatmap.
   */
  const Heatmap &getContactMap() const;

private:
  /** @name Configuration */
  ///@{
  std::string chr_name;     /**< Chromosome name. */
  int chr_length_bp;        /**< Chromosome length (bp). */
  int resolution_bp;        /**< Base pairs per bead. */
  int num_loops;            /**< Number of loop pairs. */
  int ensemble_size;        /**< Ensemble size for settings. */
  int segment_size_bp;      /**< Segment size (bp). */
  int loop_min_beads;       /**< Minimum loop span in beads. */
  int loop_max_beads;       /**< Maximum loop span in beads (capped by segment size). */
  int anchor_size_bp;       /**< Anchor width in bp. */
  float step_size;          /**< Random walk step size. */
  float contact_alpha;      /**< Distance-to-contact exponent: freq ~ dist^(-alpha). */
  std::string output_dir;   /**< Output directory. */
  std::string label;        /**< Filename label prefix. */
  ///@}

  /** @name Generated Data */
  ///@{
  std::vector<int> anchor_bead_indices;           /**< Bead indices where anchors are placed. */
  std::vector<std::pair<int, int>> loop_pairs;    /**< Pairs of anchor indices forming loops. */
  std::vector<std::pair<int, int>> loop_bead_pairs;  /**< Pairs of bead indices forming loops. */
  Chromosome ground_truth;                        /**< Ground truth polymer structure. */
  Heatmap contact_map;                            /**< Ensemble-averaged contact frequency map. */
  ///@}

  /** @name Pipeline Steps */
  ///@{
  /**
   * @brief Place anchors at random non-overlapping positions.
   * @param num_beads Total number of beads in the polymer.
   */
  void placeAnchors(int num_beads);

  /**
   * @brief Generate random loop pairs between anchors within the same segment.
   * @param num_beads Total number of beads.
   */
  void generateLoops(int num_beads);

  /**
   * @brief Generate a polymer ensemble and average the contact maps.
   * @param num_beads Total number of beads.
   */
  void generatePolymerEnsemble(int num_beads);

  /**
   * @brief Generate a single random walk polymer.
   * @param num_beads Number of beads.
   * @return Vector of 3D bead positions.
   */
  std::vector<vector3> generateRandomPolymer(int num_beads);

  /**
   * @brief Relax the polymer with harmonic loop constraints.
   * @param positions Bead positions (modified in place).
   * @param steps Number of relaxation steps.
   * @param dt Time step for gradient descent.
   */
  void relaxWithLoops(std::vector<vector3> &positions, int steps, float dt);
  ///@}

  /** @name Output Writers */
  ///@{
  void writeAnchorsFile();                  /**< Write anchors.txt (BED + orientation). */
  void writeClustersFile();                 /**< Write PET clusters file. */
  void writeSingletonsFile(int num_beads);  /**< Write distance-decay singletons. */
  void writeSegmentSplit();                 /**< Write segment split BED file. */
  void writeCentromeres();                  /**< Write centromere BED file. */
  void writeConfigINI();                    /**< Write settings.ini for reconstruction. */
  void writeGroundTruth();                  /**< Write ground truth structure file. */
  ///@}

  /** @name Utility */
  ///@{
  /**
   * @brief Sample from a Poisson distribution.
   * @param lambda Expected value (mean).
   * @return Random integer sample.
   */
  static int poissonSample(float lambda);

  /**
   * @brief Check if two beads are in the same segment.
   * @param bead_a First bead index.
   * @param bead_b Second bead index.
   * @return True if both beads belong to the same segment.
   */
  bool sameSegment(int bead_a, int bead_b) const;
  ///@}
};

#endif /* SYNTHETICGENERATOR_H_ */
