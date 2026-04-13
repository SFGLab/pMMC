#pragma once
// AppCommon.hpp — Shared CLI helpers for all pMMC apps.
// Include this exactly once (in the main_*.cpp) to get the
// arg-parsing utilities and common globals.

#include <algorithm>
#include <chrono>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdio.h>
#include <string>
#include <time.h>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

#include <platform.h>

// Global cancellation flag (REQ-5.1)
std::atomic<bool> g_cancel_requested(false);

#ifdef _WIN32
  #include "getopt.h"
#else
  #include <unistd.h>
#endif

#include <common.h>

// ── CLI argument storage ────────────────────────────────────────────
std::map<char, std::string> args;
std::string g_output_format;  // Saved from -F flag before args.clear()

// ── Helpers ─────────────────────────────────────────────────────────
inline void error_msg(std::string msg) {
  printf("%s\n", msg.c_str());
  args.clear();
  exit(0);
}

inline bool flag_set(char c) { return args.find(c) != args.end(); }

inline std::string get_arg(char c, std::string _default = "") {
  if (args.find(c) == args.end())
    return _default;
  return args[c];
}

// example formats: "chr1", "chr2,chr3", "chr8,chr10-chr15", "genome"
inline std::vector<std::string> parseChromosomeDescription(std::string desc) {
  std::vector<std::string> ret;

  if (desc == "genome") {
    for (int i = 1; i <= 22; i++)
      ret.push_back(ftext("chr%d", i));
    ret.push_back("chrX");
    return ret;
  }

  std::vector<std::string> tmp = split(desc);
  for (std::string token : tmp) {
    if (token.find('-') != std::string::npos) {
      std::vector<std::string> tmp2 = split(desc, '-');
      if (tmp2.size() == 2) {
        int st = atoi(tmp2[0].substr(3).c_str());
        int end = atoi(tmp2[1].substr(3).c_str());
        for (int i = st; i <= end; i++)
          ret.push_back(ftext("chr%d", i));
      }
    } else
      ret.push_back(token);
  }

  return ret;
}

// ── Common CLI parsing + seed setup ─────────────────────────────────
// Call this from main() after getopt loop.
// Returns the seed used.
inline unsigned int setupSeed() {
  unsigned int seed;
  if (flag_set('j')) {
    seed = (unsigned int)atoi(get_arg('j').c_str());
    printf("Using fixed random seed: %u\n", seed);
  } else {
    auto tmp = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    seed = (unsigned int)tmp;
    printf("Using time-based random seed: %u\n", seed);
  }
  srand(seed);
  setGpuSeed(seed);
  return seed;
}

// ── Common getopt loop ──────────────────────────────────────────────
inline void parseArgs(int argc, char **argv) {
  opterr = 0;
  int opt = 0;

  while ((opt = getopt(argc, argv,
                       "s:a:o:c:n:l:u:t:d:p:b:e:m:i:v:zf:hr:g:x:j:M:F:EL:I:B:P:A:")) !=
         EOF) {
    if (optarg == NULL)
      optarg = (char *)"";
    printf("opt = [%c %s]\n", opt, optarg);
    args[opt] = optarg;
  }
}
