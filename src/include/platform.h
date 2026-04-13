/*
 * platform.h - Cross-platform compatibility layer
 *
 * Provides portable mkdir, memory query, and platform-specific includes
 * for both Linux/POSIX and Windows (MSVC) builds.
 */

#ifndef PLATFORM_H_
#define PLATFORM_H_

#include <cerrno>
#include <cstdio>
#include <atomic>

#ifdef _WIN32
  #ifndef NOMINMAX
    #define NOMINMAX
  #endif
  #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
  #endif
  #include <direct.h>
  #include <windows.h>
  #include <psapi.h>
  #define portable_mkdir(path) _mkdir(path)
  typedef unsigned int uint;
#else
  #include <sys/stat.h>
  #include <sys/types.h>
  #define portable_mkdir(path) mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
  #include <cstdio>
  #include <cstring>
#endif

/* ------------------------------------------------------------------ */
/*  Global cancellation flag                                          */
/* ------------------------------------------------------------------ */
extern std::atomic<bool> g_cancel_requested;

/* ------------------------------------------------------------------ */
/*  Memory query utilities                                            */
/* ------------------------------------------------------------------ */

/** Return current resident set size (working set) in bytes.  0 on failure. */
inline size_t getCurrentRSS() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
        return pmc.WorkingSetSize;
    return 0;
#else
    FILE *f = fopen("/proc/self/status", "r");
    if (!f) return 0;
    char line[256];
    size_t rss = 0;
    while (fgets(line, sizeof(line), f)) {
        if (strncmp(line, "VmRSS:", 6) == 0) {
            long kb = 0;
            sscanf(line + 6, "%ld", &kb);
            rss = (size_t)kb * 1024;
            break;
        }
    }
    fclose(f);
    return rss;
#endif
}

/** Return peak RSS (peak working set) in bytes.  0 on failure. */
inline size_t getPeakRSS() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
        return pmc.PeakWorkingSetSize;
    return 0;
#else
    FILE *f = fopen("/proc/self/status", "r");
    if (!f) return 0;
    char line[256];
    size_t rss = 0;
    while (fgets(line, sizeof(line), f)) {
        if (strncmp(line, "VmHWM:", 6) == 0) {
            long kb = 0;
            sscanf(line + 6, "%ld", &kb);
            rss = (size_t)kb * 1024;
            break;
        }
    }
    fclose(f);
    return rss;
#endif
}

/** Print current / peak RSS to stdout with a label. */
inline void logMemoryUsage(const char *label) {
    size_t curr = getCurrentRSS();
    size_t peak = getPeakRSS();
    printf("[MEM] %-40s  current=%zuMB  peak=%zuMB\n",
           label, curr / (1024 * 1024), peak / (1024 * 1024));
}

/* ------------------------------------------------------------------ */
/*  GPU seed propagation (defined in ParallelMonteCarloHeatmap.cu)     */
/* ------------------------------------------------------------------ */
/** Set the GPU RNG seed for deterministic runs (called from main). */
void setGpuSeed(unsigned int seed);

/** Check memory budget.  Returns true if budget exceeded. budgetMB==0 means unlimited. */
inline bool checkMemoryBudget(size_t budgetMB, const char *context) {
    if (budgetMB == 0) return false;
    size_t curr = getCurrentRSS();
    size_t limitBytes = budgetMB * 1024ULL * 1024ULL;
    if (curr > limitBytes) {
        printf("\n[ERROR] Memory budget exceeded at '%s': using %zuMB, limit %zuMB\n",
               context, curr / (1024 * 1024), budgetMB);
        printf("[ERROR] Reduce dataset size or increase --max-memory-mb\n");
        return true;
    }
    return false;
}

#endif /* PLATFORM_H_ */
