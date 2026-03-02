/**
 * @file MultiscaleEnergy.h
 * @brief Explicit two-scale energy decomposition for multiscale Monte Carlo.
 *
 * Provides a structured representation of the energy landscape across
 * hierarchy levels: domain scale (chromosome/segment heatmap-driven)
 * and loop scale (anchor/subanchor arc-driven). Each scale has its own
 * energy terms, acceptance statistics, and move history.
 *
 * @see LooperSolver for the MC engine that populates these terms.
 */

#ifndef MULTISCALEENERGY_H_
#define MULTISCALEENERGY_H_

#include <cstdio>
#include <string>
#include <vector>

/**
 * @class EnergyTerm
 * @brief A single named energy contribution with tracking.
 */
class EnergyTerm {
public:
  std::string name;
  double value;
  double previous_value;
  double weight;

  EnergyTerm();
  EnergyTerm(const std::string &name, double weight);

  /** @brief Update value and save previous. */
  void update(double new_value);

  /** @brief Revert to previous value (on MC rejection). */
  void revert();

  /** @brief Weighted contribution: value * weight. */
  double weighted() const;

  /** @brief Change in weighted energy from last update. */
  double delta() const;
};

/**
 * @class ScaleEnergy
 * @brief Energy decomposition at a single hierarchy scale.
 *
 * Aggregates multiple EnergyTerm components (e.g., heatmap fit,
 * spring stretch, orientation) and tracks acceptance statistics.
 */
class ScaleEnergy {
public:
  std::string scale_name;
  std::vector<EnergyTerm> terms;

  int moves_proposed;
  int moves_accepted;
  int moves_rejected;

  ScaleEnergy();
  explicit ScaleEnergy(const std::string &name);

  /** @brief Add a new energy term to this scale. */
  void addTerm(const std::string &name, double weight);

  /** @brief Get a term by name (returns nullptr if not found). */
  EnergyTerm *getTerm(const std::string &name);

  /** @brief Total weighted energy across all terms. */
  double totalEnergy() const;

  /** @brief Total change in energy from last update across all terms. */
  double totalDelta() const;

  /** @brief Record a proposed move. */
  void recordProposal();

  /** @brief Record acceptance of the last proposal. */
  void recordAcceptance();

  /** @brief Record rejection of the last proposal. */
  void recordRejection();

  /** @brief Revert all terms to their previous values. */
  void revertAll();

  /** @brief Acceptance ratio = accepted / proposed. */
  double acceptanceRatio() const;

  /** @brief Reset all counters to zero. */
  void resetCounters();
};

/**
 * @class MultiscaleEnergy
 * @brief Two-scale energy model: domain scale + loop scale.
 *
 * Domain scale encompasses chromosome-level and segment-level
 * heatmap-driven MC. Loop scale encompasses anchor-level and
 * subanchor-level arc-distance and smoothing MC.
 *
 * Provides:
 *   - Explicit energy decomposition per scale and per term
 *   - Separate acceptance statistics per scale
 *   - Formatted output for logging and diagnostics
 */
class MultiscaleEnergy {
public:
  ScaleEnergy domain_scale;
  ScaleEnergy loop_scale;

  MultiscaleEnergy();

  /** @brief Initialize with standard energy terms for both scales. */
  void initStandardTerms();

  /** @brief Total system energy (domain + loop). */
  double totalEnergy() const;

  /** @brief Print a formatted decomposition table to stdout. */
  void printDecomposition() const;

  /**
   * @brief Write decomposition to a file.
   * @param path Output file path.
   * @param append If true, append to existing file.
   */
  void writeDecomposition(const std::string &path, bool append = false) const;

  /**
   * @brief Write a single-line summary suitable for per-step logging.
   * @param f Open FILE pointer.
   * @param step Current MC step number.
   * @param temperature Current MC temperature.
   */
  void writeStepLog(FILE *f, int step, double temperature) const;

  /** @brief Write CSV header for step logging. */
  static void writeStepLogHeader(FILE *f);
};

#endif /* MULTISCALEENERGY_H_ */
