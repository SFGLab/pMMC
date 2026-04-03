/**
 * @file SimulationLogger.hpp
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


// ============================================================================
// Implementation
// ============================================================================

#include <SvgChartGenerator.hpp>

#include <cstring>
#include <utility>

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#elif defined(__linux__)
#include <fstream>
#include <sstream>
#elif defined(__APPLE__)
#include <mach/mach.h>
#endif

// Optional CUDA header for GPU memory query
#ifdef __CUDACC__
#include <cuda_runtime.h>
#endif

inline SimulationLogger::SimulationLogger(const std::string &output_dir,
                                   const std::string &label, unsigned int seed)
    : output_dir(output_dir), label(label), seed(seed), total_steps(0),
      energy_file(nullptr), params_file(nullptr), acceptance_file(nullptr),
      summary_file(nullptr), energy_header_written(false) {

  start_time = std::chrono::high_resolution_clock::now();

  std::string energy_path = output_dir + label + "_energy.csv";
  energy_file = fopen(energy_path.c_str(), "w");

  std::string params_path = output_dir + label + "_params.log";
  params_file = fopen(params_path.c_str(), "w");

  std::string accept_path = output_dir + label + "_acceptance.csv";
  acceptance_file = fopen(accept_path.c_str(), "w");

  if (params_file) {
    fprintf(params_file, "# Simulation Parameters\n");
    fprintf(params_file, "seed = %u\n", seed);
  }

  if (acceptance_file) {
    fprintf(acceptance_file, "phase,proposed,accepted,ratio\n");
  }
}

inline SimulationLogger::~SimulationLogger() { closeFiles(); }

inline SimulationLogger::SimulationLogger(SimulationLogger &&other) noexcept
    : output_dir(std::move(other.output_dir)),
      label(std::move(other.label)), seed(other.seed),
      total_steps(other.total_steps), energy_file(other.energy_file),
      params_file(other.params_file),
      acceptance_file(other.acceptance_file),
      summary_file(other.summary_file), start_time(other.start_time),
      energy_header_written(other.energy_header_written) {
  other.energy_file = nullptr;
  other.params_file = nullptr;
  other.acceptance_file = nullptr;
  other.summary_file = nullptr;
}

inline SimulationLogger &SimulationLogger::operator=(SimulationLogger &&other) noexcept {
  if (this != &other) {
    closeFiles();
    output_dir = std::move(other.output_dir);
    label = std::move(other.label);
    seed = other.seed;
    total_steps = other.total_steps;
    energy_file = other.energy_file;
    params_file = other.params_file;
    acceptance_file = other.acceptance_file;
    summary_file = other.summary_file;
    start_time = other.start_time;
    energy_header_written = other.energy_header_written;
    other.energy_file = nullptr;
    other.params_file = nullptr;
    other.acceptance_file = nullptr;
    other.summary_file = nullptr;
  }
  return *this;
}

inline void SimulationLogger::closeFiles() {
  if (energy_file) {
    fclose(energy_file);
    energy_file = nullptr;
  }
  if (params_file) {
    fclose(params_file);
    params_file = nullptr;
  }
  if (acceptance_file) {
    fclose(acceptance_file);
    acceptance_file = nullptr;
  }
  if (summary_file) {
    fclose(summary_file);
    summary_file = nullptr;
  }
}

inline void SimulationLogger::logParameter(const std::string &name,
                                    const std::string &value) {
  if (params_file)
    fprintf(params_file, "%s = %s\n", name.c_str(), value.c_str());
}

inline void SimulationLogger::logEnergy(int step, double temperature,
                                 double domain_energy, double loop_energy,
                                 double total_energy,
                                 double domain_accept_ratio,
                                 double loop_accept_ratio) {
  if (!energy_file)
    return;

  if (!energy_header_written) {
    fprintf(energy_file,
            "step,temperature,domain_energy,loop_energy,total_energy,"
            "domain_accept_ratio,loop_accept_ratio\n");
    energy_header_written = true;
  }

  fprintf(energy_file, "%d,%.6f,%.6f,%.6f,%.6f,%.4f,%.4f\n", step,
          temperature, domain_energy, loop_energy, total_energy,
          domain_accept_ratio, loop_accept_ratio);

  total_steps = step;
}

inline void SimulationLogger::logMilestone(const std::string &message) {
  if (params_file) {
    double elapsed = getElapsedSeconds();
    fprintf(params_file, "# [%.1fs] %s\n", elapsed, message.c_str());
    fflush(params_file);
  }
}

inline void SimulationLogger::logAcceptanceRatio(const std::string &phase_name,
                                          int proposed, int accepted) {
  if (!acceptance_file)
    return;

  double ratio = (proposed > 0) ? static_cast<double>(accepted) / proposed : 0.0;
  fprintf(acceptance_file, "%s,%d,%d,%.6f\n", phase_name.c_str(), proposed,
          accepted, ratio);
  fflush(acceptance_file);
}

inline void SimulationLogger::writeSummary() {
  std::string summary_path = output_dir + label + "_summary.log";
  summary_file = fopen(summary_path.c_str(), "w");
  if (!summary_file)
    return;

  double elapsed = getElapsedSeconds();

  double rss_mb = getCurrentRSS_MB();

  fprintf(summary_file, "# Simulation Summary\n");
  fprintf(summary_file, "label = %s\n", label.c_str());
  fprintf(summary_file, "seed = %u\n", seed);
  fprintf(summary_file, "total_steps = %d\n", total_steps);
  fprintf(summary_file, "runtime_seconds = %.3f\n", elapsed);
  fprintf(summary_file, "runtime_minutes = %.2f\n", elapsed / 60.0);
  if (rss_mb >= 0.0) {
    fprintf(summary_file, "peak_rss_mb = %.1f\n", rss_mb);
  }

  fclose(summary_file);
  summary_file = nullptr;

  printf("Simulation summary written to %s\n", summary_path.c_str());

  // Close energy file so it's fully flushed before reading
  if (energy_file) {
    fclose(energy_file);
    energy_file = nullptr;
  }

  // Generate energy convergence SVG chart if energy trace was enabled
  std::string energy_csv = output_dir + label + "_energy.csv";
  FILE *echeck = fopen(energy_csv.c_str(), "r");
  if (echeck) {
    // Check it has more than just the header
    char buf[256];
    if (fgets(buf, sizeof(buf), echeck) && fgets(buf, sizeof(buf), echeck)) {
      fclose(echeck);
      std::string energy_svg = output_dir + label + "_energy.svg";
      SvgChartGenerator::energyConvergenceChart(energy_csv, energy_svg);
    } else {
      fclose(echeck);
    }
  }
}

inline void SimulationLogger::setTotalSteps(int steps) { total_steps = steps; }

inline unsigned int SimulationLogger::getSeed() const { return seed; }

inline double SimulationLogger::getElapsedSeconds() const {
  auto now = std::chrono::high_resolution_clock::now();
  return std::chrono::duration<double>(now - start_time).count();
}

inline double SimulationLogger::getCurrentRSS_MB() {
#ifdef _WIN32
  PROCESS_MEMORY_COUNTERS pmc;
  if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
    return static_cast<double>(pmc.WorkingSetSize) / (1024.0 * 1024.0);
  }
  return -1.0;
#elif defined(__linux__)
  std::ifstream status("/proc/self/status");
  std::string line;
  while (std::getline(status, line)) {
    if (line.compare(0, 6, "VmRSS:") == 0) {
      std::istringstream iss(line.substr(6));
      long rss_kb = 0;
      iss >> rss_kb;
      return static_cast<double>(rss_kb) / 1024.0;
    }
  }
  return -1.0;
#elif defined(__APPLE__)
  struct mach_task_basic_info info;
  mach_msg_type_number_t size = MACH_TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                (task_info_t)&info, &size) == KERN_SUCCESS) {
    return static_cast<double>(info.resident_size) / (1024.0 * 1024.0);
  }
  return -1.0;
#else
  return -1.0;
#endif
}

inline void SimulationLogger::logMemoryUsage(const std::string &phase_label) {
  double rss_mb = getCurrentRSS_MB();
  double elapsed = getElapsedSeconds();

  if (params_file) {
    if (rss_mb >= 0.0) {
      fprintf(params_file, "# [%.1fs] memory_rss_mb = %.1f  (%s)\n",
              elapsed, rss_mb, phase_label.c_str());
    } else {
      fprintf(params_file, "# [%.1fs] memory_rss_mb = unavailable  (%s)\n",
              elapsed, phase_label.c_str());
    }
    fflush(params_file);
  }

  // Attempt GPU memory query if CUDA is available at runtime
  // (compiled separately, so we use dlsym/weak linking approach is not
  //  portable; instead we guard with a simple function-exists check)
  printf("[%.1fs] CPU RSS: %.1f MB", elapsed, rss_mb >= 0 ? rss_mb : 0.0);
  if (!phase_label.empty())
    printf("  (%s)", phase_label.c_str());
  printf("\n");
}

#endif /* SIMULATIONLOGGER_H_ */
