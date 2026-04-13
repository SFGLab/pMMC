/**
 * @file SyntheticGenerator.hpp
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

#include <Chromosome.hpp>
#include <Heatmap.hpp>

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


// ============================================================================
// Implementation
// ============================================================================

#include <common.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include <platform.h>

inline SyntheticGenerator::SyntheticGenerator() {
  chr_name = "chr22";
  chr_length_bp = 51304566;
  resolution_bp = 25000;
  num_loops = 100;
  ensemble_size = 50;
  segment_size_bp = 2000000;
  loop_min_beads = 10;
  loop_max_beads = 400;
  anchor_size_bp = 5000;
  step_size = 1.0f;
  contact_alpha = 2.0f;
  output_dir = "./synthetic/";
  label = "synthetic";
}

inline void SyntheticGenerator::setChromosome(const std::string &chr, int length_bp) {
  chr_name = chr;
  chr_length_bp = length_bp;
}

inline void SyntheticGenerator::setResolution(int bp_per_bead) {
  resolution_bp = bp_per_bead;
}

inline void SyntheticGenerator::setNumLoops(int n) { num_loops = n; }

inline void SyntheticGenerator::setEnsembleSize(int n) { ensemble_size = n; }

inline void SyntheticGenerator::setSegmentSize(int bp) { segment_size_bp = bp; }

inline void SyntheticGenerator::setOutputDir(const std::string &dir) {
  output_dir = dir;
}

inline void SyntheticGenerator::setLabel(const std::string &lbl) { label = lbl; }

inline const Chromosome &SyntheticGenerator::getGroundTruth() const {
  return ground_truth;
}

inline const Heatmap &SyntheticGenerator::getContactMap() const {
  return contact_map;
}

inline bool SyntheticGenerator::sameSegment(int bead_a, int bead_b) const {
  int pos_a = bead_a * resolution_bp;
  int pos_b = bead_b * resolution_bp;
  if (pos_a > pos_b)
    std::swap(pos_a, pos_b);
  for (int boundary = segment_size_bp; boundary < chr_length_bp;
       boundary += segment_size_bp) {
    if (pos_a <= boundary && pos_b >= boundary)
      return false;
  }
  return true;
}

inline void SyntheticGenerator::generate() {
  portable_mkdir(output_dir.c_str());

  int num_beads = chr_length_bp / resolution_bp;
  printf("Synthetic data generation:\n");
  printf("  Chromosome: %s (%d bp)\n", chr_name.c_str(), chr_length_bp);
  printf("  Resolution: %d bp/bead -> %d beads\n", resolution_bp, num_beads);
  printf("  Loops: %d\n", num_loops);
  printf("  Ensemble size: %d\n", ensemble_size);

  // Cap loop_max_beads to segment size
  int beads_per_segment = segment_size_bp / resolution_bp;
  if (loop_max_beads > beads_per_segment - 2) {
    printf("  Capping loop_max_beads from %d to %d (segment size)\n",
           loop_max_beads, beads_per_segment - 2);
    loop_max_beads = beads_per_segment - 2;
  }

  placeAnchors(num_beads);
  generateLoops(num_beads);
  printf("  Generated %d anchors, %d loop pairs\n",
         (int)anchor_bead_indices.size(), (int)loop_pairs.size());

  generatePolymerEnsemble(num_beads);
  printf("  Ensemble contact map built (%dx%d)\n", (int)contact_map.size,
         (int)contact_map.size);

  writeAnchorsFile();
  writeClustersFile();
  writeSingletonsFile(num_beads);
  writeSegmentSplit();
  writeCentromeres();
  writeConfigINI();
  writeGroundTruth();

  printf("  All files written to %s\n", output_dir.c_str());
}

inline void SyntheticGenerator::placeAnchors(int num_beads) {
  anchor_bead_indices.clear();

  int num_anchors = std::max(num_loops * 2, num_loops + 20);
  num_anchors = std::min(num_anchors, num_beads / 2);

  // Generate candidate positions and shuffle
  std::vector<int> candidates;
  for (int i = 2; i < num_beads - 2; ++i)
    candidates.push_back(i);

  for (int i = (int)candidates.size() - 1; i > 0; --i) {
    int j = std::rand() % (i + 1);
    std::swap(candidates[i], candidates[j]);
  }

  // Sort a subset and pick with minimum spacing
  int subset = std::min(num_anchors * 3, (int)candidates.size());
  std::sort(candidates.begin(), candidates.begin() + subset);

  int min_spacing = 3;
  int last_used = -min_spacing - 1;
  for (int i = 0;
       i < (int)candidates.size() && (int)anchor_bead_indices.size() < num_anchors;
       ++i) {
    if (candidates[i] - last_used >= min_spacing) {
      anchor_bead_indices.push_back(candidates[i]);
      last_used = candidates[i];
    }
  }

  std::sort(anchor_bead_indices.begin(), anchor_bead_indices.end());
}

inline void SyntheticGenerator::generateLoops(int /*num_beads*/) {
  loop_pairs.clear();
  loop_bead_pairs.clear();

  int max_attempts = num_loops * 200;
  int attempts = 0;
  int rejected_cross_seg = 0;

  while ((int)loop_pairs.size() < num_loops && attempts < max_attempts) {
    attempts++;
    int a = std::rand() % (int)anchor_bead_indices.size();
    int b = std::rand() % (int)anchor_bead_indices.size();
    if (a == b)
      continue;
    if (a > b)
      std::swap(a, b);

    int bead_dist = anchor_bead_indices[b] - anchor_bead_indices[a];
    if (bead_dist < loop_min_beads || bead_dist > loop_max_beads)
      continue;

    // Ensure both anchors are within the same segment
    if (!sameSegment(anchor_bead_indices[a], anchor_bead_indices[b])) {
      rejected_cross_seg++;
      continue;
    }

    // Check for duplicate
    bool dup = false;
    for (const auto &lp : loop_pairs) {
      if (lp.first == a && lp.second == b) {
        dup = true;
        break;
      }
    }
    if (dup)
      continue;

    loop_pairs.push_back({a, b});
    loop_bead_pairs.push_back(
        {anchor_bead_indices[a], anchor_bead_indices[b]});
  }

  if (rejected_cross_seg > 0)
    printf("  Rejected %d cross-segment loop candidates\n", rejected_cross_seg);

  if ((int)loop_pairs.size() < num_loops)
    printf("  Warning: only generated %d/%d loops\n", (int)loop_pairs.size(),
           num_loops);
}

inline std::vector<vector3> SyntheticGenerator::generateRandomPolymer(int num_beads) {
  std::vector<vector3> positions(num_beads);

  positions[0] = vector3(0.0f, 0.0f, 0.0f);
  for (int i = 1; i < num_beads; ++i) {
    vector3 step = random_vector(step_size, false);
    float len = step.length();
    if (len > 1e-8f)
      step = step * (step_size / len);
    positions[i] = positions[i - 1] + step;
  }

  // Relax with loop constraints
  if (!loop_bead_pairs.empty())
    relaxWithLoops(positions, 150, 0.002f);

  return positions;
}

inline void SyntheticGenerator::relaxWithLoops(std::vector<vector3> &positions,
                                        int steps, float dt) {
  int n = (int)positions.size();
  float loop_target = step_size * 2.0f;
  float excluded_r = step_size * 0.3f;
  float excluded_r2 = excluded_r * excluded_r;

  for (int step = 0; step < steps; ++step) {
    std::vector<vector3> forces(n, vector3(0, 0, 0));

    // Chain connectivity springs
    for (int i = 0; i < n - 1; ++i) {
      vector3 diff = positions[i + 1] - positions[i];
      float dist = diff.length();
      if (dist > 1e-8f) {
        float stretch = dist - step_size;
        vector3 f = diff * (stretch / dist);
        forces[i] = forces[i] + f;
        forces[i + 1] = forces[i + 1] - f;
      }
    }

    // Loop constraint springs
    for (const auto &lp : loop_bead_pairs) {
      vector3 diff = positions[lp.second] - positions[lp.first];
      float dist = diff.length();
      if (dist > 1e-8f) {
        float stretch = dist - loop_target;
        float k = 2.0f; // stronger spring for loops
        vector3 f = diff * (k * stretch / dist);
        forces[lp.first] = forces[lp.first] + f;
        forces[lp.second] = forces[lp.second] - f;
      }
    }

    // Excluded volume (sparse check for performance)
    int check_step = std::max(1, n / 500);
    for (int i = 0; i < n; i += check_step) {
      for (int j = i + 2; j < n && j < i + 50; ++j) {
        vector3 diff = positions[j] - positions[i];
        float d2 = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
        if (d2 < excluded_r2 && d2 > 1e-12f) {
          float d = std::sqrt(d2);
          float repulsion = (excluded_r - d) / d;
          vector3 f = diff * repulsion;
          forces[i] = forces[i] - f;
          forces[j] = forces[j] + f;
        }
      }
    }

    // Apply forces
    for (int i = 0; i < n; ++i)
      positions[i] = positions[i] + forces[i] * dt;
  }
}

inline void SyntheticGenerator::generatePolymerEnsemble(int num_beads) {
  int beads_per_bin = segment_size_bp / resolution_bp;
  if (beads_per_bin < 1)
    beads_per_bin = 1;
  int num_bins = (num_beads + beads_per_bin - 1) / beads_per_bin;

  contact_map.init(num_bins);
  contact_map.diagonal_size = 1;

  std::vector<std::vector<double>> accum_dist(
      num_bins, std::vector<double>(num_bins, 0.0));
  std::vector<std::vector<int>> accum_count(num_bins,
                                            std::vector<int>(num_bins, 0));

  for (int e = 0; e < ensemble_size; ++e) {
    if ((e + 1) % 10 == 0 || e == 0)
      printf("  Generating polymer %d/%d\n", e + 1, ensemble_size);

    std::vector<vector3> positions = generateRandomPolymer(num_beads);

    // Save first realization as ground truth
    if (e == 0) {
      ground_truth.points = positions;
      ground_truth.size = num_beads;
      ground_truth.genomic_position.resize(num_beads);
      for (int i = 0; i < num_beads; ++i)
        ground_truth.genomic_position[i] = i * resolution_bp;
    }

    // Binned contact frequencies
    for (int bi = 0; bi < num_bins; ++bi) {
      int ci = std::min(bi * beads_per_bin + beads_per_bin / 2, num_beads - 1);
      for (int bj = bi + 1; bj < num_bins; ++bj) {
        int cj =
            std::min(bj * beads_per_bin + beads_per_bin / 2, num_beads - 1);
        float dist = (positions[ci] - positions[cj]).length();
        if (dist > 1e-6f) {
          accum_dist[bi][bj] += 1.0 / std::pow(dist, contact_alpha);
          accum_count[bi][bj]++;
        }
      }
    }
  }

  // Average and fill symmetric
  for (int i = 0; i < num_bins; ++i) {
    for (int j = i + 1; j < num_bins; ++j) {
      float freq = 0.0f;
      if (accum_count[i][j] > 0)
        freq = (float)(accum_dist[i][j] / accum_count[i][j]);
      contact_map.v[i][j] = freq;
      contact_map.v[j][i] = freq;
    }
  }

  // Normalize so near-diagonal average ~ 1.0
  float diag_sum = 0.0f;
  int diag_count = 0;
  for (int i = 0; i < num_bins - 1; ++i) {
    diag_sum += contact_map.v[i][i + 1];
    diag_count++;
  }
  float diag_avg = (diag_count > 0) ? diag_sum / diag_count : 1.0f;
  if (diag_avg > 1e-10f) {
    float s = 1.0f / diag_avg;
    for (int i = 0; i < num_bins; ++i)
      for (int j = 0; j < num_bins; ++j)
        contact_map.v[i][j] *= s;
  }
}

inline void SyntheticGenerator::writeAnchorsFile() {
  std::string path = output_dir + "anchors.txt";
  FILE *f = fopen(path.c_str(), "w");
  if (!f) {
    printf("Error: cannot write %s\n", path.c_str());
    return;
  }

  for (int idx : anchor_bead_indices) {
    int genomic_pos = idx * resolution_bp;
    int start = genomic_pos - anchor_size_bp / 2;
    int end = genomic_pos + anchor_size_bp / 2;
    if (start < 0)
      start = 0;
    if (end > chr_length_bp)
      end = chr_length_bp;

    char orient = 'N';
    int r = std::rand() % 3;
    if (r == 0)
      orient = 'L';
    else if (r == 1)
      orient = 'R';

    fprintf(f, "%s\t%d\t%d\t%c\n", chr_name.c_str(), start, end, orient);
  }
  fclose(f);
}

inline void SyntheticGenerator::writeClustersFile() {
  std::string path = output_dir + "clusters.txt";
  FILE *f = fopen(path.c_str(), "w");
  if (!f) {
    printf("Error: cannot write %s\n", path.c_str());
    return;
  }

  for (const auto &lp : loop_pairs) {
    int anchor_a = anchor_bead_indices[lp.first];
    int anchor_b = anchor_bead_indices[lp.second];

    int start_a = anchor_a * resolution_bp - anchor_size_bp / 2;
    int end_a = anchor_a * resolution_bp + anchor_size_bp / 2;
    int start_b = anchor_b * resolution_bp - anchor_size_bp / 2;
    int end_b = anchor_b * resolution_bp + anchor_size_bp / 2;

    if (start_a < 0)
      start_a = 0;
    if (start_b < 0)
      start_b = 0;
    if (end_a > chr_length_bp)
      end_a = chr_length_bp;
    if (end_b > chr_length_bp)
      end_b = chr_length_bp;

    // PET count biased by genomic distance
    int genomic_dist = std::abs(anchor_b - anchor_a) * resolution_bp;
    int pet_count = 3 + std::rand() % 15;
    if (genomic_dist < 500000)
      pet_count += 5 + std::rand() % 10;

    fprintf(f, "%s\t%d\t%d\t%s\t%d\t%d\t%d\n", chr_name.c_str(), start_a,
            end_a, chr_name.c_str(), start_b, end_b, pet_count);
  }
  fclose(f);
}

inline void SyntheticGenerator::writeSingletonsFile(int num_beads) {
  std::string path = output_dir + "singletons.txt";
  FILE *f = fopen(path.c_str(), "w");
  if (!f) {
    printf("Error: cannot write %s\n", path.c_str());
    return;
  }

  int beads_per_bin = segment_size_bp / resolution_bp;
  if (beads_per_bin < 1)
    beads_per_bin = 1;
  int num_bins = (int)contact_map.size;
  // Root-cause fix: LooperSolver::normalizeHeatmap computes row-wise sums and
  // divides by them (expected_sum / row_sum).  If a segment row is all-zero
  // (no singletons connect that segment to any other), the division yields Inf,
  // and 0 * Inf = NaN which propagates to diagonal_avg and breaks the MC.
  //
  // For loops_100/1000 the last LooperSolver segment spans chr22 tail 24-51Mb
  // (~27Mb with no anchors).  The contact_map entries for long-range bin pairs
  // crossing that boundary are essentially 0 → all skipped by the freq<0.001
  // guard → entire segment row stays zero → NaN.
  //
  // Fix: guarantee at least 1 singleton per off-diagonal bin pair regardless
  // of contact frequency.  This ensures every segment row is non-zero while
  // adding only O(num_bins^2) minimal extra singletons.
  float singleton_scale = 100.0f;
  int total = 0;

  for (int bi = 0; bi < num_bins; ++bi) {
    for (int bj = bi + 1; bj < num_bins; ++bj) {
      float freq = contact_map.v[bi][bj];

      // Always write at least 1 singleton per pair so that every LooperSolver
      // segment row has non-zero connectivity (prevents NaN in normalizeHeatmap).
      int count = std::max(1, poissonSample(freq * singleton_scale));

      int start_i = bi * beads_per_bin * resolution_bp;
      int end_i =
          std::min((bi + 1) * beads_per_bin * resolution_bp, chr_length_bp);
      int start_j = bj * beads_per_bin * resolution_bp;
      int end_j =
          std::min((bj + 1) * beads_per_bin * resolution_bp, chr_length_bp);

      for (int r = 0; r < count; ++r) {
        int pos1 = start_i + std::rand() % std::max(1, end_i - start_i);
        int pos2 = start_j + std::rand() % std::max(1, end_j - start_j);
        int read_len = 50 + std::rand() % 100;

        fprintf(f, "%s\t%d\t%d\t%s\t%d\t%d\t1\n", chr_name.c_str(), pos1,
                pos1 + read_len, chr_name.c_str(), pos2, pos2 + read_len);
        total++;
      }
    }
  }
  fclose(f);
  printf("  Generated %d singleton interactions\n", total);
}

inline void SyntheticGenerator::writeSegmentSplit() {
  std::string path = output_dir + "segments.bed";
  FILE *f = fopen(path.c_str(), "w");
  if (!f)
    return;

  for (int pos = segment_size_bp; pos < chr_length_bp;
       pos += segment_size_bp) {
    fprintf(f, "%s %d %d\n", chr_name.c_str(), pos, pos);
  }
  fclose(f);
}

inline void SyntheticGenerator::writeCentromeres() {
  std::string path = output_dir + "centromeres.bed";
  FILE *f = fopen(path.c_str(), "w");
  if (!f)
    return;
  fprintf(f, "%s\t1\t100\n", chr_name.c_str());
  fclose(f);
}

inline void SyntheticGenerator::writeConfigINI() {
  std::string path = output_dir + "settings.ini";
  FILE *f = fopen(path.c_str(), "w");
  if (!f)
    return;

  fprintf(f,
          "[main]\n"
          "output_level = 2\n"
          "random_walk = no\n"
          "loop_density = 5\n"
          "steps_lvl1 = 2\n"
          "steps_lvl2 = 2\n"
          "steps_arcs = 3\n"
          "steps_smooth = 3\n"
          "noise_lvl1 = 0.5\n"
          "noise_lvl2 = 0.5\n"
          "noise_arcs = 0.01\n"
          "noise_smooth = 5.0\n"
          "max_pet_length = 1000000\n"
          "long_pet_power = 2.0\n"
          "long_pet_scale = 1.0\n"
          "\n"
          "[cuda]\n"
          "blocks_multiplier = 4\n"
          "num_threads = 512\n"
          "milestone_fails = 3\n"
          "\n"
          "[data]\n"
          "data_dir = %s\n"
          "anchors = anchors.txt\n"
          "clusters = clusters.txt\n"
          "factors = synthetic\n"
          "singletons = singletons.txt\n"
          "split_singleton_files_by_chr = no\n"
          "singletons_inter = \n"
          "segment_split = %ssegments.bed\n"
          "centromeres = %scentromeres.bed\n"
          "\n"
          "[distance]\n"
          "freq_dist_scale = 25.0\n"
          "freq_dist_power = -0.6\n"
          "freq_dist_scale_inter = 120.0\n"
          "freq_dist_power_inter = -1.0\n"
          "count_dist_a = 0.2\n"
          "count_dist_scale = 1.8\n"
          "count_dist_shift = 8\n"
          "count_dist_base_level = 0.2\n"
          "genomic_dist_power = 0.5\n"
          "genomic_dist_scale = 1.0\n"
          "genomic_dist_base = 0.0\n"
          "\n"
          "[template]\n"
          "template_scale = 7.0\n"
          "dist_heatmap_scale = 15.0\n"
          "\n"
          "[motif_orientation]\n"
          "use_motif_orientation = no\n"
          "weight = 50.0\n"
          "\n"
          "[anchor_heatmap]\n"
          "use_anchor_heatmap = no\n"
          "heatmap_influence = 0.5\n"
          "\n"
          "[subanchor_heatmap]\n"
          "use_subanchor_heatmap = no\n"
          "estimate_distances_steps = 4\n"
          "estimate_distances_replicates = 4\n"
          "heatmap_influence = 0.1\n"
          "heatmap_dist_weight = 0.01\n"
          "\n"
          "[heatmaps]\n"
          "inter_scaling = 1.0\n"
          "distance_heatmap_stretching = 2.5\n"
          "\n"
          "[springs]\n"
          "stretch_constant = 0.1\n"
          "squeeze_constant = 0.1\n"
          "angular_constant = 0.1\n"
          "stretch_constant_arcs = 1.0\n"
          "squeeze_constant_arcs = 1.0\n"
          "\n"
          "[simulation_heatmap]\n"
          "max_temp_heatmap = 5.0\n"
          "delta_temp_heatmap = 0.9999\n"
          "jump_temp_scale_heatmap = 50.0\n"
          "jump_temp_coef_heatmap = 20.0\n"
          "stop_condition_improvement_threshold_heatmap = 0.99\n"
          "stop_condition_successes_threshold_heatmap = 10\n"
          "stop_condition_steps_heatmap = 20000\n"
          "\n"
          "[simulation_arcs]\n"
          "max_temp = 5.0\n"
          "delta_temp = 0.9999\n"
          "jump_temp_scale = 50.0\n"
          "jump_temp_coef = 20.0\n"
          "stop_condition_improvement_threshold = 0.975\n"
          "stop_condition_successes_threshold = 100\n"
          "stop_condition_steps = 20000\n"
          "\n"
          "[simulation_arcs_smooth]\n"
          "dist_weight = 1.0\n"
          "angle_weight = 1.0\n"
          "max_temp = 5.0\n"
          "delta_temp = 0.9999\n"
          "jump_temp_scale = 50.0\n"
          "jump_temp_coef = 20.0\n"
          "stop_condition_improvement_threshold = 0.99\n"
          "stop_condition_successes_threshold = 50\n"
          "stop_condition_steps = 20000\n",
          output_dir.c_str(), output_dir.c_str(), output_dir.c_str());

  fclose(f);
}

inline void SyntheticGenerator::writeGroundTruth() {
  std::string path = output_dir + "ground_truth.txt";
  FILE *f = fopen(path.c_str(), "w");
  if (!f)
    return;

  fprintf(f, "%d %d\n", (int)ground_truth.points.size(), resolution_bp);
  for (size_t i = 0; i < ground_truth.points.size(); ++i) {
    int gp = (int)i * resolution_bp;
    fprintf(f, "%d %f %f %f\n", gp, ground_truth.points[i].x,
            ground_truth.points[i].y, ground_truth.points[i].z);
  }
  fclose(f);
}

inline int SyntheticGenerator::poissonSample(float lambda) {
  if (lambda <= 0.0f)
    return 0;
  if (lambda > 30.0f) {
    float u1 = (float)std::rand() / RAND_MAX;
    float u2 = (float)std::rand() / RAND_MAX;
    float z = std::sqrt(-2.0f * std::log(u1 + 1e-10f)) *
              std::cos(6.2832f * u2);
    int val = (int)(lambda + std::sqrt(lambda) * z + 0.5f);
    return std::max(0, val);
  }
  double L = std::exp(-(double)lambda);
  int k = 0;
  double p = 1.0;
  do {
    k++;
    p *= (double)std::rand() / RAND_MAX;
  } while (p > L);
  return k - 1;
}

#endif /* SYNTHETICGENERATOR_H_ */
