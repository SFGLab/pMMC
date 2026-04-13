/**
 * @file MultiscaleEnergy.hpp
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


// ============================================================================
// Implementation
// ============================================================================


#include <cmath>
#include <cstdio>

// --- EnergyTerm ---

inline EnergyTerm::EnergyTerm() : name(""), value(0.0), previous_value(0.0), weight(1.0) {}

inline EnergyTerm::EnergyTerm(const std::string &name, double weight)
    : name(name), value(0.0), previous_value(0.0), weight(weight) {}

inline void EnergyTerm::update(double new_value) {
  previous_value = value;
  value = new_value;
}

inline void EnergyTerm::revert() { value = previous_value; }

inline double EnergyTerm::weighted() const { return value * weight; }

inline double EnergyTerm::delta() const {
  return (value - previous_value) * weight;
}

// --- ScaleEnergy ---

inline ScaleEnergy::ScaleEnergy()
    : scale_name("unnamed"), moves_proposed(0), moves_accepted(0),
      moves_rejected(0) {}

inline ScaleEnergy::ScaleEnergy(const std::string &name)
    : scale_name(name), moves_proposed(0), moves_accepted(0),
      moves_rejected(0) {}

inline void ScaleEnergy::addTerm(const std::string &name, double weight) {
  terms.push_back(EnergyTerm(name, weight));
}

inline EnergyTerm *ScaleEnergy::getTerm(const std::string &name) {
  for (auto &t : terms) {
    if (t.name == name)
      return &t;
  }
  return nullptr;
}

inline double ScaleEnergy::totalEnergy() const {
  double sum = 0.0;
  for (const auto &t : terms)
    sum += t.weighted();
  return sum;
}

inline double ScaleEnergy::totalDelta() const {
  double sum = 0.0;
  for (const auto &t : terms)
    sum += t.delta();
  return sum;
}

inline void ScaleEnergy::recordProposal() { moves_proposed++; }

inline void ScaleEnergy::recordAcceptance() { moves_accepted++; }

inline void ScaleEnergy::recordRejection() { moves_rejected++; }

inline void ScaleEnergy::revertAll() {
  for (auto &t : terms)
    t.revert();
}

inline double ScaleEnergy::acceptanceRatio() const {
  if (moves_proposed == 0)
    return 0.0;
  return static_cast<double>(moves_accepted) / moves_proposed;
}

inline void ScaleEnergy::resetCounters() {
  moves_proposed = 0;
  moves_accepted = 0;
  moves_rejected = 0;
}

// --- MultiscaleEnergy ---

inline MultiscaleEnergy::MultiscaleEnergy()
    : domain_scale("domain"), loop_scale("loop") {}

inline void MultiscaleEnergy::initStandardTerms() {
  // Domain scale: heatmap-driven energy terms
  domain_scale.addTerm("heatmap_fit", 1.0);
  domain_scale.addTerm("inter_chr_heatmap", 1.0);
  domain_scale.addTerm("density_fit", 1.0);

  // Loop scale: arc and smoothing energy terms
  loop_scale.addTerm("arc_distance", 1.0);
  loop_scale.addTerm("linker_stretch", 1.0);
  loop_scale.addTerm("linker_squeeze", 1.0);
  loop_scale.addTerm("angular_bending", 1.0);
  loop_scale.addTerm("ctcf_orientation", 1.0);
  loop_scale.addTerm("subanchor_heatmap", 1.0);
}

inline double MultiscaleEnergy::totalEnergy() const {
  return domain_scale.totalEnergy() + loop_scale.totalEnergy();
}

inline void MultiscaleEnergy::printDecomposition() const {
  printf("\n========================================\n");
  printf("  MULTISCALE ENERGY DECOMPOSITION\n");
  printf("========================================\n");

  printf("\n  DOMAIN SCALE (%s):\n", domain_scale.scale_name.c_str());
  printf("  %-25s %12s %12s %12s\n", "Term", "Raw", "Weight", "Weighted");
  printf("  %-25s %12s %12s %12s\n", "----", "---", "------", "--------");
  for (const auto &t : domain_scale.terms) {
    if (std::fabs(t.value) > 1e-12)
      printf("  %-25s %12.4f %12.4f %12.4f\n", t.name.c_str(), t.value,
             t.weight, t.weighted());
  }
  printf("  %-25s %12s %12s %12.4f\n", "SUBTOTAL", "", "",
         domain_scale.totalEnergy());
  printf("  Acceptance: %d/%d (%.3f)\n", domain_scale.moves_accepted,
         domain_scale.moves_proposed, domain_scale.acceptanceRatio());

  printf("\n  LOOP SCALE (%s):\n", loop_scale.scale_name.c_str());
  printf("  %-25s %12s %12s %12s\n", "Term", "Raw", "Weight", "Weighted");
  printf("  %-25s %12s %12s %12s\n", "----", "---", "------", "--------");
  for (const auto &t : loop_scale.terms) {
    if (std::fabs(t.value) > 1e-12)
      printf("  %-25s %12.4f %12.4f %12.4f\n", t.name.c_str(), t.value,
             t.weight, t.weighted());
  }
  printf("  %-25s %12s %12s %12.4f\n", "SUBTOTAL", "", "",
         loop_scale.totalEnergy());
  printf("  Acceptance: %d/%d (%.3f)\n", loop_scale.moves_accepted,
         loop_scale.moves_proposed, loop_scale.acceptanceRatio());

  printf("\n  TOTAL ENERGY: %.4f\n", totalEnergy());
  printf("========================================\n\n");
}

inline void MultiscaleEnergy::writeDecomposition(const std::string &path,
                                          bool append) const {
  FILE *f = fopen(path.c_str(), append ? "a" : "w");
  if (!f)
    return;

  if (!append) {
    fprintf(f, "# Multiscale Energy Decomposition\n");
    fprintf(f, "scale,term,raw_value,weight,weighted_value\n");
  }

  for (const auto &t : domain_scale.terms) {
    fprintf(f, "%s,%s,%.6f,%.6f,%.6f\n", domain_scale.scale_name.c_str(),
            t.name.c_str(), t.value, t.weight, t.weighted());
  }
  for (const auto &t : loop_scale.terms) {
    fprintf(f, "%s,%s,%.6f,%.6f,%.6f\n", loop_scale.scale_name.c_str(),
            t.name.c_str(), t.value, t.weight, t.weighted());
  }

  fclose(f);
}

inline void MultiscaleEnergy::writeStepLogHeader(FILE *f) {
  fprintf(f, "step,temperature,domain_energy,loop_energy,total_energy,"
             "domain_accept_ratio,loop_accept_ratio\n");
}

inline void MultiscaleEnergy::writeStepLog(FILE *f, int step,
                                    double temperature) const {
  fprintf(f, "%d,%.6f,%.6f,%.6f,%.6f,%.4f,%.4f\n", step, temperature,
          domain_scale.totalEnergy(), loop_scale.totalEnergy(), totalEnergy(),
          domain_scale.acceptanceRatio(), loop_scale.acceptanceRatio());
}

#endif /* MULTISCALEENERGY_H_ */
