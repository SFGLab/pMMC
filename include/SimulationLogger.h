/**
 * @file SimulationLogger.h
 * @brief RAII-based simulation logging for reproducibility and diagnostics.
 *
 * Tracks and persists:
 *   - Random seed used for the run
 *   - All simulation parameters
 *   - Per-step energy values and acceptance ratios
 *   - Runtime statistics
 *   - Final summary
 *
 * Automatically flushes and closes on destruction (RAII).
 */

#ifndef SIMULATIONLOGGER_H_
#define SIMULATIONLOGGER_H_

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

/**
 * @class SimulationLogger
 * @brief Logs simulation state to files for reproducibility.
 *
 * Usage:
 * @code
 *   SimulationLogger logger("./output/", "run1", 42);
 *   logger.logParameter("temperature", "5.0");
 *   logger.logEnergy(step, temp, domain_e, loop_e, total_e, domain_ar, loop_ar);
 *   // ... destructor closes files automatically
 * @endcode
 */
class SimulationLogger {
public:
  /**
   * @brief Construct and open log files.
   * @param output_dir Directory for log files.
   * @param label Run label for filenames.
   * @param seed Random seed used.
   */
  SimulationLogger(const std::string &output_dir, const std::string &label,
                   unsigned int seed);

  /** @brief Flush and close all open log files. */
  ~SimulationLogger();

  // Non-copyable, movable
  SimulationLogger(const SimulationLogger &) = delete;
  SimulationLogger &operator=(const SimulationLogger &) = delete;
  SimulationLogger(SimulationLogger &&other) noexcept;
  SimulationLogger &operator=(SimulationLogger &&other) noexcept;

  /**
   * @brief Log a parameter name-value pair.
   * @param name Parameter name.
   * @param value Parameter value as string.
   */
  void logParameter(const std::string &name, const std::string &value);

  /**
   * @brief Log a per-step energy entry.
   * @param step MC step number.
   * @param temperature Current temperature.
   * @param domain_energy Domain scale energy.
   * @param loop_energy Loop scale energy.
   * @param total_energy Total energy.
   * @param domain_accept_ratio Domain acceptance ratio.
   * @param loop_accept_ratio Loop acceptance ratio.
   */
  void logEnergy(int step, double temperature, double domain_energy,
                 double loop_energy, double total_energy,
                 double domain_accept_ratio, double loop_accept_ratio);

  /**
   * @brief Log a milestone event (e.g., level transition).
   * @param message Descriptive message.
   */
  void logMilestone(const std::string &message);

  /**
   * @brief Log acceptance ratio for a completed MC phase.
   * @param phase_name Phase identifier (e.g., "heatmap_chr", "arcs_anchor").
   * @param proposed Number of proposed moves.
   * @param accepted Number of accepted moves.
   */
  void logAcceptanceRatio(const std::string &phase_name, int proposed,
                          int accepted);

  /** @brief Write the final summary (runtime, total steps, etc.). */
  void writeSummary();

  /** @brief Override the total step count (e.g., from LooperSolver). */
  void setTotalSteps(int steps);

  /** @brief Get the seed used for this run. */
  unsigned int getSeed() const;

  /** @brief Get elapsed time since construction in seconds. */
  double getElapsedSeconds() const;

  /**
   * @brief Log current memory usage (CPU RSS and GPU if available).
   * @param phase_label Optional label for the measurement point.
   */
  void logMemoryUsage(const std::string &phase_label = "");

  /**
   * @brief Get current process RSS (Resident Set Size) in MB.
   * @return RSS in megabytes, or -1.0 if unavailable.
   */
  static double getCurrentRSS_MB();

private:
  void closeFiles();

  std::string output_dir;
  std::string label;
  unsigned int seed;
  int total_steps;

  FILE *energy_file;
  FILE *params_file;
  FILE *acceptance_file;
  FILE *summary_file;

  std::chrono::high_resolution_clock::time_point start_time;
  bool energy_header_written;
};

#endif /* SIMULATIONLOGGER_H_ */
