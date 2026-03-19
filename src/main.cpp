#include <algorithm>
#include <chrono>
#include <map>
#include <sstream>
#include <stdio.h>
#include <string>
#include <time.h>
#include <vector>

#include <platform.h>

// Global cancellation flag (REQ-5.1)
std::atomic<bool> g_cancel_requested(false);

#ifdef _WIN32
  #include "getopt.h"
#else
  #include <unistd.h>
#endif

#include <BedRegion.h>
#include <BenchmarkRunner.h>
#include <CifWriter.h>
#include <DistanceMapGenerator.h>
#include <GenomeReconstructor.h>
#include <HierarchicalChromosome.h>
#include <LooperSolver.h>
#include <IbedFileIO.h>
#include <MergedFileIO.h>
#include <MetricsFramework.h>
#include <MultiscaleEnergy.h>
#include <PdbWriter.h>
#include <SimulationLogger.h>
#include <SyntheticGenerator.h>
#include <HeatmapImageWriter.h>
#include <common.h>

std::map<char, std::string> args;
std::string g_output_format;  // Saved from -F flag before args.clear()

void usage(const char *act = "", bool quit = true) {
  printf("Usage: pMMC -a <action> [options]\n\n");
  printf("Actions:\n");
  printf("  create (c)     Reconstruct 3D structure (default)\n");
  printf("  generate (g)   Generate synthetic polymer + ChIA-PET data\n");
  printf("  benchmark (k)  Run benchmark (generate + reconstruct + compare)\n");
  printf("  distmap        Compute distance/contact maps from .hcm file\n");
  printf("  metrics        Compute structural comparison metrics\n");
  printf("  ensemble (e)   Ensemble analysis of multiple structures\n");
  printf("  extract (r)    Extract fragment of structure\n");
  printf("  position (p)   Get 3D positions for BED regions\n");
  printf("  smooth (s)     Create equidistant smoothed model\n");
  printf("  distance (d)   Compute distance matrix from .hcm file\n");
  printf("  flatten (f)    Flatten hierarchical structure to text\n");
  printf("  rewiring (w)   Rewiring analysis\n");
  printf("  merge (M)      Merge BED + BEDPE into a single merged file\n");
  printf("  simulate (S)   Reconstruct from merged file\n");
  printf("  ibed (I)       Reconstruct from .ibed file\n");
  printf("  convert (X)    Convert HCM file to PDB/CIF format\n");
  printf("\nCommon options:\n");
  printf("  -a  action to perform\n");
  printf("  -s  path to settings file\n");
  printf("  -n  label for output file names\n");
  printf("  -c  chromosomes/region (e.g. 'chr14:1:2500000', default: genome)\n");
  printf("  -o  output directory (with trailing '/', default: './')\n");
  printf("  -j  random seed for reproducibility (default: time-based)\n");
  printf("  -l  loop count (generate) or length limit (create)\n");
  printf("  -m  ensemble size\n");
  printf("  -i  input file/directory\n");
  printf("\nOutput format:\n");
  printf("  -F  additional output format: cif, mmcif, pdb, or both (HCM always written)\n");
  printf("  -E  enable per-step energy trace CSV\n");
  printf("  -L  energy log interval in MC steps (default: 1000)\n");
  printf("\nMemory:\n");
  printf("  -M  max memory budget in MB (0 = unlimited, default: 0)\n");
  printf("\nInitialization:\n");
  printf("  -I  initialization method: random (default), mds\n");
  printf("\nMerge/Simulate options:\n");
  printf("  -B  BED file path (for merge action)\n");
  printf("  -P  BEDPE file path (for merge action)\n");
  printf("\nDebugging:\n");
  printf("  -u  max chromosomes to reconstruct\n");

  if (quit) {
    args.clear();
    exit(0);
  }
}

void error_msg(string msg) {
  printf("%s\n", msg.c_str());
  args.clear();
  exit(0);
}

bool flag_set(char c) { return args.find(c) != args.end(); }

string get_arg(char c, string _default = "") {
  if (args.find(c) == args.end())
    return _default;
  return args[c];
}

// Save model in the requested output format(s).
// Always saves HCM (backward compat). Additionally saves CIF/PDB if -F flag requests it.
void saveModelWithFormat(HierarchicalChromosome &hc, const string &outdir,
                         const string &label, bool use_new_file_format,
                         const string &suffix = "",
                         bool useCurrentLevel = false) {
  // Always save HCM
  string hcm_path;
  if (use_new_file_format)
    hcm_path = ftext("%sloops_new_%s%s.hcm", outdir.c_str(), label.c_str(), suffix.c_str());
  else
    hcm_path = ftext("%sloops_%s%s.hcm", outdir.c_str(), label.c_str(), suffix.c_str());

  if (use_new_file_format)
    hc.toFile(hcm_path);
  else
    hc.toFilePreviousFormat(hcm_path);

  // Write additional formats (CIF/PDB) based on -F flag saved before args.clear()
  if (!g_output_format.empty()) {
    string base = ftext("%s%s%s", outdir.c_str(), label.c_str(), suffix.c_str());
    if (g_output_format == "cif" || g_output_format == "mmcif") {
      CifWriter::write(hc, base + ".cif", "cudaMMC", useCurrentLevel);
    } else if (g_output_format == "pdb") {
      // PDB always uses lowest level for maximum atom detail
      PdbWriter::write(hc, base + ".pdb", false);
    } else if (g_output_format == "both") {
      CifWriter::write(hc, base + ".cif", "cudaMMC", useCurrentLevel);
      PdbWriter::write(hc, base + ".pdb", false);
    } else {
      printf("[WARN] Unknown output format '%s'. Use: cif, pdb, or both\n", g_output_format.c_str());
    }
  }
}

void get3DPositions(string input_structure, BedRegions &regions,
                    string output_file, bool add_header = true) {
  // read 3D model
  HierarchicalChromosome chr;
  chr.fromFile(input_structure);
  chr.setLevel(3);

  printv(chr.chrs, "chromosomes found");

  FILE *f = open(output_file, "w");
  if (f == NULL)
    error("could not open the output file");
  if (add_header)
    fprintf(f, "%zu\n", regions.regions.size());

  for (uint i = 0; i < regions.regions.size(); i++) {
    int p = (regions.regions[i].start + regions.regions[i].end) / 2;
    vector3 pos = chr.find3DPosition(regions.regions[i].chr, p);
    printf("%d -> %f %f %f\n", p, pos.x, pos.y, pos.z);
    fprintf(f, "%f %f %f\n", pos.x, pos.y, pos.z);
  }

  fclose(f);
}

void preparePosition() {
  if (!flag_set('i'))
    error_msg("No input structure specified [-i]\n");
  if (!flag_set('o'))
    error_msg("No output file specified [-o]\n");
  if (!flag_set('b'))
    error_msg("No bed regions specified [-b]\n");

  string input = get_arg('i');
  string output = get_arg('o');
  string regions_path = get_arg('b');
  bool add_header = flag_set('h');

  // read regions to extract the position for
  BedRegions regions;
  regions.fromFile(regions_path);
  printf("%d regions loaded\n", (int)regions.regions.size());

  get3DPositions(input, regions, output, add_header);
}

void prepareSmooth() {
  if (!flag_set('i'))
    error_msg("No input structure specified [-i]\n");
  if (!flag_set('r'))
    error_msg("No resolution specified [-r]\n");

  string input = get_arg('i');
  int res = atoi(get_arg('r').c_str());

  string output = flag_set('o') ? get_arg('o') : input + ".smooth.txt";
  string regions_path = get_arg('b');
  int level = flag_set('l') ? atoi(get_arg('l').c_str()) : 3;
  string chr = flag_set('c') ? get_arg('c') : "";

  HierarchicalChromosome hc;
  hc.fromFile(input);
  hc.setLevel(level);
  // hc.createCurrentLevelStructure();

  Chromosome ch = hc.createEqidistantModel(res, chr);
  ch.toFile(output);

  // create artificial regions
  // BedRegions regions;
  // regions.addNewIntervals("chr", 1000, 2000, 50);
  // regions.print();
  // get3DPositions(input, regions, output, true);
}

// calculate structural distance matrix (i.e. matrix with distances between
// every pair of beads) input:
//    - hc: hierarchical chromosome model
//    - chr: chromosome to calculate the distances for (ignored if level=0)
//    - level: level of the structure (0-chr, 1-segment, 2-subanchor)
// output: filename where the output matrix is to be saved
void runCalcDistanceMatrix(HierarchicalChromosome hc, string chr, int level,
                           string output) {
  hc.setLevel(level);
  Heatmap h = hc.createStructuralHeatmap(chr, level);
  h.toFile(output);
}

void prepareCalcDistanceMatrix() {
  if (!flag_set('i'))
    error_msg("No input structure specified [-i]\n");
  if (!flag_set('l'))
    error_msg("No level specified [-l]\n");
  if (!flag_set('c'))
    error_msg("No chromosome specified (may be anything if level=0) [-c]\n");

  string input = get_arg('i');
  string output =
      flag_set('o') ? get_arg('o') : ftext("%s.dist.heat", input.c_str());
  string chr = get_arg('c');
  int level = atoi(get_arg('l').c_str());

  HierarchicalChromosome hc;
  hc.fromFile(input);
  runCalcDistanceMatrix(hc, chr, level, output);
}

// given a path to .hcm file create a flat representation of the specified
// chromosome on a given level if chr=="" then save all the chromosomes save the
// result to txt file output_path
void flatten(string path, int level, string chr = "", string output_path = "") {
  HierarchicalChromosome hc;
  hc.fromFile(path);
  hc.setLevel(level);
  hc.createCurrentLevelStructure();

  string out_path = output_path == "" ? path : output_path;
  if (chr == "") {
    for (auto const &item : hc.chr) {
      printf("do [%s]\n", item.first.c_str());
      if (hc.chr.size() > 1)
        hc.chr[item.first].toFile(
            ftext("%s.%s.txt", path.c_str(), item.first.c_str()));
      else
        hc.chr[item.first].toFile(out_path);
    }
  } else {
    if (hc.chr.find(chr) == hc.chr.end()) {
      printf("Could not find chromosome [%s]\nExisting chromosomes:",
             chr.c_str());
      for (auto const &item : hc.chr)
        printf(" %s", item.first.c_str());
      return;
    }
    hc.chr[chr].toFile(ftext("%s.txt", path.c_str()));
  }
}

void prepareFlatten() {
  if (!flag_set('i'))
    error_msg("No input structure specified [-i]\n");
  if (!flag_set('l'))
    error_msg("No level specified [-l]\n");

  string input = get_arg('i');
  int level = atoi(get_arg('l').c_str());

  string chr = flag_set('c') ? get_arg('c') : "";
  string output = flag_set('o') ? get_arg('o') : input + ".txt";
  flatten(input, level, chr, output);
}

void extract(string input_path, string output_path, int start, int end) {
  printf("extract fragment: %d %d\n", start, end);
  HierarchicalChromosome h;
  h.fromFile(input_path);
  printf("size = %d\n", (int)h.clusters.size());
  if (h.clusters.size() > 0) {
    HierarchicalChromosome extr = h.extractFragment(start, end);
    extr.toFile(output_path);
  }
}

void prepareExtract() {
  if (!flag_set('i'))
    error_msg("No input specified [-i]\n");
  if (!flag_set('o'))
    error_msg("No input specified [-o]\n");
  if (!flag_set('s'))
    error_msg("No start position specified [-s]\n");
  if (!flag_set('e'))
    error_msg("No end position specified [-e]\n");

  string input = get_arg('i');
  string output = get_arg('o');

  int start = atoi(get_arg('s').c_str());
  int end = atoi(get_arg('e').c_str());

  extract(input, output, start, end);
}

void ensembleAnalysis(string input_dir, string pattern, string output_dir) {
  // read all structures (hierarchical), and for each level (chr, segment,
  // subanchor) get the chr then calculate the matrix of pairwise distances
  // between structures
  printf("ensemble analysis:\n");

  // first read all hierarchical structures to the memory
  vector<HierarchicalChromosome> hchr; // keep all structures loaded
  string repl = "{N}";

  if (pattern.find(repl) == string::npos)
    error_msg(ftext("pattern does not contain replacement token (\"%s\")",
                    repl.c_str()));

  int p = 0;
  while (1) {
    string str = pattern;
    str.replace(str.find(repl), repl.size(), ftext("%d", p));
    string path = input_dir + str;
    if (!file_exists(path))
      break;
    printf("[%s]\n", path.c_str());
    p++;

    HierarchicalChromosome h;
    h.fromFile(path);
    hchr.push_back(h);

    if (p > 10000)
      break; // safety check to prevent infinite loop
  }

  int N = static_cast<int>(hchr.size()); // number of structures in the ensemble
  if (N < 2) {
    printf("Too few structures (%d) to conduct an analysis\n", N);
    return;
  }

  // now calculate the distances matrices
  vector<Heatmap>
      heat; // keep structural heatmaps for all structures on a given level
  Heatmap heat_pairs(N); // heatmap with pair-wise distances

  vector<string> chrs = hchr[0].chrs; // save chromosomes

  // first chr level
  if (chrs.size() < 2)
    printf("Only %d chromosome available, skip chr level\n", (int)chrs.size());
  else {
    printf("chr level\n");

    printf("create structural heatmaps\n");

    for (int j = 0; j < N; ++j) { // for all structures
      heat.push_back(hchr[j].createStructuralHeatmap("", 0));
    }

    printf("create distance heatmap\n");
    for (int i = 0; i < N; ++i) { // for all structures
      for (int j = i + 1; j < N; ++j) {
        heat_pairs.v[i][j] = heat_pairs.v[j][i] = heat[i].calcDistance(heat[j]);
      }
    }

    printf("save\n");
    heat_pairs.toFile(input_dir + "structural_distances_lvl0.heat", false);
    heat_pairs.zero();

    heat.clear();
  }

  // segment and subanchor levels
  int i, j;
  for (j = 0; j < N; ++j)
    hchr[j].levelDown(); // move to the segment level

  float d;
  for (int lvl = 1; lvl <= 2; ++lvl) {
    // i is current level (1-segment, 2-subanchor)
    printf("do level %d\n", lvl);

    // go through each chromosome one by one
    for (string chr : chrs) {
      for (j = 0; j < N; ++j) { // for each structure in the ensemble
        hchr[j].setLevel(lvl);
        heat.push_back(hchr[j].createStructuralHeatmap(chr, lvl));
      }

      for (i = 0; i < N; ++i) { // for all structures
        for (j = i + 1; j < N; ++j) {
          d = heat[i].calcDistance(heat[j]);
          heat_pairs.v[i][j] += d;
          heat_pairs.v[j][i] += d;
        }
      }
    }

    heat_pairs.toFile(input_dir + ftext("structural_distances_lvl%d.heat", lvl),
                      false);

    heat_pairs.zero();
    heat.clear();
  }
}

void prepareEnsembleAnalysis() {
  if (!flag_set('i'))
    error_msg("No input directory specified [-i]\n");
  if (!flag_set('p'))
    error_msg("No pattern specified [-p]\n");

  string input = get_arg('i');
  string pattern = get_arg('p');
  string output = flag_set('o') ? get_arg('o') : input;

  ensembleAnalysis(input, pattern, output);
}

void rewiringAnalysis(string input_dir, string pattern, int level = 1,
                      int resolution = 1000000, string output_dir = "") {
  // create structural heatmaps for all structures (hierarchical)
  printf("rewiring analysis:\n");

  string repl = "{N}";
  if (pattern.find(repl) == string::npos)
    error_msg(ftext("pattern does not contain the replacement token (\"%s\")",
                    repl.c_str()));

  int p = 0;
  while (1) {
    string str = pattern;
    str.replace(str.find(repl), repl.size(), ftext("%d", p));
    string path = input_dir + str;
    if (!file_exists(path)) {
      if (p == 0)
        error(ftext("No file (%s) found\n", path.c_str()));
      break;
    }
    printf("[%s]\n", path.c_str());
    p++;

    HierarchicalChromosome h;
    h.fromFile(path);
    h.setLevel(level);
    // hc.createCurrentLevelStructure();

    // go through each chromosome one by one
    for (string chr : h.chrs) {

      printf("   chr %s\n", chr.c_str());

      Chromosome ch = h.createEqidistantModel(resolution, chr);
      Heatmap h = ch.createInverseHeatmap();
      printf("save to %s\n",
             ftext("%s_%s_%d.heat", path.c_str(), chr.c_str(), p).c_str());
      h.toFile(ftext("%s_%s_%d.heat", path.c_str(), chr.c_str(), p));
    }

    if (p > 10000)
      break; // safety check to prevent infinite loop
  }
}

void prepareRewiring() {
  if (!flag_set('i'))
    error_msg("No input directory specified [-i]\n");
  if (!flag_set('p'))
    error_msg("No pattern specified [-p]\n");

  string input = get_arg('i');
  string pattern = get_arg('p');
  int level = flag_set('l') ? atoi(get_arg('l').c_str()) : 1;
  int resolution = flag_set('r') ? atoi(get_arg('r').c_str()) : 1000000;
  string output = flag_set('o') ? get_arg('o') : input;

  rewiringAnalysis(input, pattern, level, resolution, output);
}

// example formats: "chr1", "chr2,chr3", "chr8,chr10-chr15", "genome"
std::vector<std::string> parseChromosomeDescription(string desc) {
  std::vector<string> ret;

  if (desc == "genome") {
    for (int i = 1; i <= 22; i++)
      ret.push_back(ftext("chr%d", i));
    ret.push_back("chrX");
    return ret;
  }

  std::vector<string> tmp = split(desc);
  for (string token : tmp) {
    if (token.find('-') != string::npos) {
      std::vector<string> tmp2 = split(desc, '-');
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

// ---------------------------------------------------------------------------
// Output validation: create ERROR_WARNING.txt only when serious issues exist
// ---------------------------------------------------------------------------

// Return file size in bytes, or -1 if the file cannot be opened.
static long getFileSize(const std::string &path) {
  FILE *f = fopen(path.c_str(), "rb");
  if (!f) return -1;
  fseek(f, 0, SEEK_END);
  long sz = ftell(f);
  fclose(f);
  return sz;
}

// Validate expected output files and write ERROR_WARNING.txt only if serious
// issues are found (missing or empty output files).  If everything is fine,
// remove any stale ERROR_WARNING.txt left over from a previous run.
void validateOutputsAndWriteErrorWarning(const std::string &outdir,
                                         const std::string &label,
                                         bool use_new_file_format,
                                         int ensemble_size) {
  std::string ew_path = outdir + "ERROR_WARNING.txt";
  std::vector<std::string> missing_files;
  std::vector<std::string> empty_files;

  // Helper: check one candidate file.  If it exists and is empty, record it.
  // If it does not exist, record it.
  auto checkFile = [&](const std::string &path, const std::string &desc) {
    long sz = getFileSize(path);
    if (sz < 0) {
      missing_files.push_back(desc + " (" + path + ")");
    } else if (sz == 0) {
      empty_files.push_back(desc + " (" + path + ")");
    }
  };

  // Helper: check a file that has two possible names (old vs new format).
  // At least one of the two should exist and be non-empty.
  auto checkFileEither = [&](const std::string &pathA,
                             const std::string &pathB,
                             const std::string &desc) {
    long szA = getFileSize(pathA);
    long szB = getFileSize(pathB);
    bool a_ok = (szA > 0);
    bool b_ok = (szB > 0);
    if (a_ok || b_ok) return;  // at least one is fine
    // Both are problematic.  Report the one that exists-but-empty, or missing.
    if (szA == 0) {
      empty_files.push_back(desc + " (" + pathA + ")");
    } else if (szB == 0) {
      empty_files.push_back(desc + " (" + pathB + ")");
    } else {
      missing_files.push_back(desc + " (" + pathA + " or " + pathB + ")");
    }
  };

  if (ensemble_size <= 1) {
    // --- Single-run expected outputs ---
    // Final HCM
    {
      std::string hcmA = ftext("%sloops_%s.hcm", outdir.c_str(), label.c_str());
      std::string hcmB = ftext("%sloops_new_%s.hcm", outdir.c_str(), label.c_str());
      checkFileEither(hcmA, hcmB, "Final HCM");
    }
    // Final PDB
    {
      std::string pdb = ftext("%s%s.pdb", outdir.c_str(), label.c_str());
      checkFile(pdb, "Final PDB");
    }
    // Initial HCM
    {
      std::string hcmA = ftext("%sloops_%s_initial.hcm", outdir.c_str(), label.c_str());
      std::string hcmB = ftext("%sloops_new_%s_initial.hcm", outdir.c_str(), label.c_str());
      checkFileEither(hcmA, hcmB, "Initial HCM");
    }
    // Initial PDB
    {
      std::string pdb = ftext("%s%s_initial.pdb", outdir.c_str(), label.c_str());
      checkFile(pdb, "Initial PDB");
    }
  } else {
    // --- Ensemble expected outputs ---
    for (int i = 0; i < ensemble_size; ++i) {
      std::string idx = ftext("_%d", i);
      // Final HCM
      {
        std::string hcmA = ftext("%sloops_%s%s.hcm", outdir.c_str(), label.c_str(), idx.c_str());
        std::string hcmB = ftext("%sloops_new_%s%s.hcm", outdir.c_str(), label.c_str(), idx.c_str());
        checkFileEither(hcmA, hcmB, ftext("Ensemble %d final HCM", i));
      }
      // Final PDB
      {
        std::string pdb = ftext("%s%s%s.pdb", outdir.c_str(), label.c_str(), idx.c_str());
        checkFile(pdb, ftext("Ensemble %d final PDB", i));
      }
      // Initial HCM
      {
        std::string ini_suffix = ftext("_%d_initial", i);
        std::string hcmA = ftext("%sloops_%s%s.hcm", outdir.c_str(), label.c_str(), ini_suffix.c_str());
        std::string hcmB = ftext("%sloops_new_%s%s.hcm", outdir.c_str(), label.c_str(), ini_suffix.c_str());
        checkFileEither(hcmA, hcmB, ftext("Ensemble %d initial HCM", i));
      }
      // Initial PDB
      {
        std::string pdb = ftext("%s%s_%d_initial.pdb", outdir.c_str(), label.c_str(), i);
        checkFile(pdb, ftext("Ensemble %d initial PDB", i));
      }
    }
  }

  // Check .heat files: look for the known suffixes
  {
    std::vector<std::string> heat_suffixes = {
      "singletons_chromosomes_", "singletons_segment_",
      "singletons_chromosomes_norm_", "singletons_segment_norm_",
      "heat_dist_chromosomes_", "heat_dist_segment_"
    };
    for (const auto &suffix : heat_suffixes) {
      std::string heat_path = ftext("%s%s%s.heat", outdir.c_str(), suffix.c_str(), label.c_str());
      long sz = getFileSize(heat_path);
      // Only flag if the file exists but is empty.
      // Heat files are optional — not all runs produce all of them —
      // so a missing heat file is NOT an error.
      if (sz == 0) {
        empty_files.push_back("Heatmap file (" + heat_path + ")");
      }
    }
  }

  // --- Decision: create ERROR_WARNING.txt only if there are genuine issues ---
  bool has_issues = !missing_files.empty() || !empty_files.empty();

  if (!has_issues) {
    // No problems — remove stale ERROR_WARNING.txt from a previous run
    ::remove(ew_path.c_str());
    return;
  }

  // Write ERROR_WARNING.txt
  FILE *f = fopen(ew_path.c_str(), "w");
  if (!f) {
    printf("[WARN] Could not create %s\n", ew_path.c_str());
    return;
  }

  fprintf(f, "ERROR / WARNING REPORT\n");
  fprintf(f, "======================\n\n");
  fprintf(f, "One or more expected output files are missing or empty.\n");
  fprintf(f, "The 3D reconstruction output may be incomplete or unusable.\n\n");

  if (!missing_files.empty()) {
    fprintf(f, "MISSING FILES (expected but not found):\n");
    for (const auto &m : missing_files) {
      fprintf(f, "  - %s\n", m.c_str());
    }
    fprintf(f, "\n  Possible causes: the reconstruction failed or was cancelled\n");
    fprintf(f, "  before these files could be written, the output directory is\n");
    fprintf(f, "  not writable, or insufficient disk space.\n\n");
  }

  if (!empty_files.empty()) {
    fprintf(f, "EMPTY FILES (exist but contain 0 bytes):\n");
    for (const auto &e : empty_files) {
      fprintf(f, "  - %s\n", e.c_str());
    }
    fprintf(f, "\n  Possible causes: the file was created but the write failed\n");
    fprintf(f, "  (disk full, I/O error), or the input data produced a\n");
    fprintf(f, "  degenerate structure with no content to write.\n\n");
  }

  fprintf(f, "These issues seriously affect output quality.  Please check\n");
  fprintf(f, "the log output above for additional error messages.\n");
  fclose(f);

  printf("[WARN] Wrote %s — output validation found issues.\n", ew_path.c_str());
}

// ensemble_mode - whether to create 1 structure or a whole ensemble. accepted
// values:
//	   -1 - single structure
//		0 - generate N independent structures
//		1 - generate structures in hierarchical way (start with N
//structures for chr level, then for each of these create N 			structures using the
//parent as a template, and so on.)
// ensemble_multiplicity - controls how many structures to create for an
// ensemble, corresponds to 'N' selected factors - if provided, uses only the
// cluster and singleton files related to specified clusters (order in stg file
// matters). If used, the number of factors must be equal to
//    		the number of cluster and singleton files. If empty, all input files
//    are used
void runLooper(std::vector<string> chrs, BedRegion region_of_interest,
               string label, string outdir, bool use_new_file_format,
               int max_level, int chr_number_limit = -1, int length_limit = -1,
               int ensemble_size = 1,
               const std::vector<string> selected_factors = vector<string>(),
               SimulationLogger *logger = nullptr) {

  printf("run [%s]\n", label.c_str());

  int res = portable_mkdir(outdir.c_str());
  if (res != 0 && errno != EEXIST)
    error("Couldn't create the output directory");

  string path_data = Settings::dataDirectory;
  string path_anchors = Settings::dataAnchors;
  std::vector<string> path_pets = split(Settings::dataPetClusters);
  std::vector<string> path_singletons = split(Settings::dataSingletons);
  std::vector<string> path_singletons_inter =
      split(Settings::dataSingletonsInter);
  std::vector<string> factors = split(Settings::dataFactors);

  if (selected_factors.size() > 0) {
    // check if there is a propoer number of cluster/singletons files
    if (path_pets.size() != path_singletons.size() ||
        factors.size() != path_pets.size()) {
      error("Invalid number of input files! There should be exactly one "
            "cluster and singleton file per factor when -g option is used");
      return;
    }

    // remove unwanted files
    for (uint i = static_cast<uint>(factors.size()) - 1; i >= 1; --i) {
      bool found = false;
      for (uint j = 0; j < selected_factors.size(); ++j) {
        if (factors[i] == selected_factors[j])
          found = true;
      }
      if (!found) {
        printf("ignore factor: %s\n", factors[i].c_str());
        path_pets.erase(path_pets.begin() + i);
        path_singletons.erase(path_singletons.begin() + i);
        factors.erase(factors.begin() + i);
      }
    }
  }

  string anchors = ftext("%s%s", path_data.c_str(), path_anchors.c_str());

  std::vector<string> arcs_clusters;
  std::vector<string> arcs_singletons;
  std::vector<string> arcs_singletons_inter;

  for (uint i = 0; i < path_pets.size(); ++i) {
    arcs_clusters.push_back(
        ftext("%s%s", path_data.c_str(), path_pets[i].c_str()));
  }

  for (uint i = 0; i < path_singletons.size(); ++i) {
    arcs_singletons.push_back(
        ftext("%s%s", path_data.c_str(), path_singletons[i].c_str()));
  }

  for (uint i = 0; i < path_singletons_inter.size(); ++i) {
    arcs_singletons_inter.push_back(
        ftext("%s%s", path_data.c_str(), path_singletons_inter[i].c_str()));
  }

  LooperSolver lsm(label, outdir);
  lsm.setDebuggingOptions(chr_number_limit, length_limit);
  if (logger)
    lsm.setLogger(logger);

  lsm.setContactData(chrs, region_of_interest, anchors, factors, arcs_clusters,
                     arcs_singletons, arcs_singletons_inter);

  if (Settings::useDensity)
    lsm.initDensityData();

  // print_vector(lsm.densityCoordToStructure(vector3(17, 42, 28)), "telomer1");
  // print_vector(lsm.densityCoordToStructure(vector3(13, 62, 20)), "telomer2");

  logMemoryUsage("before createTreeGenome");
  printf("create tree\n");
  lsm.createTreeGenome();
  logMemoryUsage("after createTreeGenome");

  // Load explicit loop constraints if provided in config
  if (!Settings::loopConstraintsFile.empty()) {
    printf("loading loop constraints from %s\n",
           Settings::loopConstraintsFile.c_str());
    lsm.loadLoopConstraints(Settings::loopConstraintsFile);
  }

  if (ensemble_size == 1) {
    // proceed the normal way
    logMemoryUsage("before reconstructClustersHeatmap");
    printf("\nreconstruct heatmap\n");
    auto t_heatmap_start = chrono::high_resolution_clock::now();
    lsm.reconstructClustersHeatmap();
    auto t_heatmap_end = chrono::high_resolution_clock::now();
    double heatmap_ms = chrono::duration<double, std::milli>(t_heatmap_end - t_heatmap_start).count();
    printf("[TIMING] Heatmap phase: %.1f ms\n", heatmap_ms);
    logMemoryUsage("after reconstructClustersHeatmap");

    // Always save initial conformation (HCM + PDB) for reproducibility.
    {
        printf("save initial conformation\n");
        HierarchicalChromosome hc_initial = lsm.getModel();
        // HCM: save at current level (preserves segment-level view)
        saveModelWithFormat(hc_initial, outdir, label, use_new_file_format, "_initial", true);
        // PDB: always write using lowest level for maximum atom detail
        if (g_output_format != "pdb" && g_output_format != "both") {
            std::string pdb_path = ftext("%s%s_initial.pdb", outdir.c_str(), label.c_str());
            PdbWriter::write(hc_initial, pdb_path, false);
        }
    }

    if (max_level >= LVL_INTERACTION_BLOCK) {
      printf("reconstruct arcs\n");
      auto t_arcs_start = chrono::high_resolution_clock::now();
      lsm.reconstructClustersArcsDistances();
      auto t_arcs_end = chrono::high_resolution_clock::now();
      double arcs_ms = chrono::duration<double, std::milli>(t_arcs_end - t_arcs_start).count();
      printf("[TIMING] Arcs+Smooth phase: %.1f ms\n", arcs_ms);
      logMemoryUsage("after reconstructClustersArcsDistances");
    }

    // Print final multiscale energy decomposition
    lsm.energy.printDecomposition();

    printf("save\n");
    HierarchicalChromosome hc = lsm.getModel();

    // hc.center(false); // be careful! when using density this will shift the
    // structure from the density map

    lsm.showCoordinatesRange(1);

    saveModelWithFormat(hc, outdir, label, use_new_file_format);

    // Always ensure PDB is written for final conformation (R2.1)
    if (g_output_format != "pdb" && g_output_format != "both") {
        std::string pdb_path = ftext("%s%s.pdb", outdir.c_str(), label.c_str());
        PdbWriter::write(hc, pdb_path);
    }
  } else if (ensemble_size > 1) {
    // generate N independent structures
    printf("\nreconstruct %d structures\n", ensemble_size);

    // note: subanchor beads are added during
    // reconstructClustersArcsDistances(), so we need to remove them between
    // runs otherwise all the subanchor beads from one run will be treated as
    // anchor beads in the next one

    for (int i = 0; i < ensemble_size; ++i) {
      printf("structure %d/%d\n", i + 1, ensemble_size);

      // REQ-5.1: cancellation check between ensemble members
      if (g_cancel_requested.load(std::memory_order_relaxed)) {
        printf("\n[CANCEL] Ensemble interrupted after %d/%d structures\n", i, ensemble_size);
        break;
      }
      // REQ-4.1: memory budget check between ensemble members
      if (checkMemoryBudget(Settings::maxMemoryMB, "ensemble loop")) {
        printf("[ERROR] Aborting ensemble at member %d due to memory budget.\n", i);
        break;
      }
      logMemoryUsage(ftext("ensemble member %d/%d start", i + 1, ensemble_size).c_str());

      printf("reconstruct heatmap\n");
      lsm.reconstructClustersHeatmap();

      // Always save initial conformation (HCM + PDB) before arc reconstruction.
      {
          printf("save initial conformation %d/%d\n", i + 1, ensemble_size);
          HierarchicalChromosome hc_initial = lsm.getModel();
          std::string ens_suffix = ftext("_%d_initial", i);
          saveModelWithFormat(hc_initial, outdir, label, use_new_file_format,
                              ens_suffix, true);
          // Ensure PDB is always written for initial conformation, even if -F
          // flag was not set to "pdb" or "both".
          if (g_output_format != "pdb" && g_output_format != "both") {
              std::string pdb_path = ftext("%s%s%s.pdb", outdir.c_str(), label.c_str(), ens_suffix.c_str());
              PdbWriter::write(hc_initial, pdb_path, false);
          }
      }

      if (max_level >= LVL_INTERACTION_BLOCK) {
        printf("reconstruct arcs\n");
        lsm.reconstructClustersArcsDistances();
      }

      HierarchicalChromosome hc = lsm.getModel();
      // hc.center(false);

      saveModelWithFormat(hc, outdir, label, use_new_file_format,
                          ftext("_%d", i));

      // Always ensure PDB is written for ensemble final conformation (R2.1)
      if (g_output_format != "pdb" && g_output_format != "both") {
          std::string pdb_path = ftext("%s%s_%d.pdb", outdir.c_str(), label.c_str(), i);
          PdbWriter::write(hc, pdb_path);
      }

      // remove subanchor beads
      lsm.removeSubanchorBeads();
      lsm.reset();
    }
  } else
    error("ensemble size not correct");

  // Generate heatmap images from .heat files (R4.2, R4.3)
  {
    std::vector<std::string> heat_suffixes = {
        "singletons_chromosomes_", "singletons_segment_",
        "singletons_chromosomes_norm_", "singletons_segment_norm_",
        "heat_dist_chromosomes_", "heat_dist_segment_"
    };
    for (const auto& suffix : heat_suffixes) {
      std::string heat_path = ftext("%s%s%s.heat", outdir.c_str(), suffix.c_str(), label.c_str());
      if (file_exists(heat_path)) {
        Heatmap h;
        h.fromFile(heat_path);
        if (h.size > 0) {
          std::string base = heat_path.substr(0, heat_path.size() - 5); // remove .heat
          std::string title = suffix + label;
          HeatmapImageWriter::writePNG(h, base + ".png", title);
          HeatmapImageWriter::writeSVG(h, base + ".svg", title);
        }
      }
    }
  }

  // Transfer accumulated MC step count to the logger
  if (logger)
    logger->setTotalSteps(lsm.total_mc_steps);
  printf("Total MC steps accumulated: %d\n", lsm.total_mc_steps);

  // Validate outputs and create ERROR_WARNING.txt only if serious issues exist
  validateOutputsAndWriteErrorWarning(outdir, label, use_new_file_format,
                                      ensemble_size);
}

void prepareLooper(SimulationLogger *logger = nullptr) {
  printf("prepare\n");

  // we need to know what regions we should reconstruct
  // we can either provide a number of chromosomes, or a single chromosome
  // region in the first case we need to put all chromosomes id to 'chrs' in the
  // second we put region info in 'region_of_interest', but we also need to
  // update 'chrs' to contain the selected chromosome (we need it because all
  // tree manipulations are based on 'chrs')

  string chromosomes =
      "genome"; // which chromosomes (eg. "chr2", "chr1-chr22,chrX")
  BedRegion region_of_interest;
  std::vector<string> chrs;
  if (flag_set('c')) {
    chromosomes = get_arg('c');
    if (BedRegion::tryParse(
            chromosomes)) { // check if the value is a BED region
      // it seems like it is, so parse the description and update 'chrs'
      region_of_interest.parse(chromosomes);
      chrs.push_back(region_of_interest.chr);
    } else
      chrs = parseChromosomeDescription(chromosomes);
  } else {
    // if there was no -c parameter use default parameters
    chrs = parseChromosomeDescription(chromosomes);
  }

  string name = "";
  if (flag_set('n'))
    name = get_arg('n');
  else {
    name = chromosomes;
    std::replace(name.begin(), name.end(), ':', '_');
    std::replace(name.begin(), name.end(), '-', '_');
  }

  if (flag_set('t') && flag_set('d'))
    error("-t and -d options are exclusive, you may only provide one of them");

  std::vector<string> selected_factors = vector<string>();
  if (flag_set('g'))
    selected_factors = split(get_arg('g'));

  // if there is a structural template provided
  if (flag_set('t')) {
    Settings::templateSegment = get_arg('t');
    if (flag_set('p'))
      Settings::templateScale = static_cast<float>(atof(get_arg('p').c_str()));
  }

  // if there is a distance heatmap provided
  if (flag_set('d')) {
    Settings::distHeatmap = get_arg('d');
    if (flag_set('b'))
      Settings::distHeatmapScale = static_cast<float>(atof(get_arg('b').c_str()));
  }
  if (flag_set('f')) {
    Settings::dataSegmentHeatmap = get_arg('f');
    if (flag_set('b'))
      Settings::distHeatmapScale = static_cast<float>(atof(get_arg('b').c_str()));
  }

  // -x allows to override the values in the settings file
  if (flag_set('x')) {
    vector<string> arr = split(get_arg('x'), ',');
    for (std::size_t i = 0; i < arr.size(); ++i) {
      vector<string> arr_par = split(arr[i], ':');

      if (arr_par[0] == "ib_random_walk_jumps")
        Settings::ibRandomWalkJumps = static_cast<float>(atof(arr_par[1].c_str()));
      else if (arr_par[0] == "freq_dist_scale")
        Settings::freqToDistHeatmapScale = static_cast<float>(atof(arr_par[1].c_str()));
      else if (arr_par[0] == "freq_dist_power")
        Settings::freqToDistHeatmapPower = static_cast<float>(atof(arr_par[1].c_str()));
      else if (arr_par[0] == "freq_dist_scale_inter")
        Settings::freqToDistHeatmapScaleInter = static_cast<float>(atof(arr_par[1].c_str()));
      else if (arr_par[0] == "freq_dist_power_inter")
        Settings::freqToDistHeatmapPowerInter = static_cast<float>(atof(arr_par[1].c_str()));
      else if (arr_par[0] == "count_dist_a")
        Settings::countToDistA = static_cast<float>(atof(arr_par[1].c_str()));
      else if (arr_par[0] == "count_dist_scale")
        Settings::countToDistScale = static_cast<float>(atof(arr_par[1].c_str()));
      else if (arr_par[0] == "count_dist_shift")
        Settings::countToDistShift = static_cast<float>(atof(arr_par[1].c_str()));
      else if (arr_par[0] == "count_dist_base_level")
        Settings::countToDistBaseLevel = static_cast<float>(atof(arr_par[1].c_str()));
      else if (arr_par[0] == "genomic_dist_power")
        Settings::genomicLengthToDistPower = static_cast<float>(atof(arr_par[1].c_str()));
      else if (arr_par[0] == "genomic_dist_scale")
        Settings::genomicLengthToDistScale = static_cast<float>(atof(arr_par[1].c_str()));
      else if (arr_par[0] == "genomic_dist_base")
        Settings::genomicLengthToDistBase = static_cast<float>(atof(arr_par[1].c_str()));
    }
  }

  bool use_new_file_format = flag_set('z');
  string outdir = flag_set('o') ? get_arg('o') : "./";
  int length_limit = flag_set('l') ? atoi(get_arg('l').c_str()) : -1;
  int chr_number_limit = flag_set('u') ? atoi(get_arg('u').c_str()) : -1;
  int max_level = flag_set('v') ? atoi(get_arg('v').c_str()) : 5;

  int ensemble_size = flag_set('m') ? atoi(get_arg('m').c_str()) : 1;

  bool ignore_cache = flag_set('h');
  if (ignore_cache)
    Settings::useInputCache = false;

  // Save -F flag before args.clear() so saveModelWithFormat can use it
  g_output_format = flag_set('F') ? get_arg('F') : "";

  args.clear(); // clean here, otherwise heap corruption occurs when program
                // abnormally stops (i.e. exit(0))

  runLooper(chrs, region_of_interest, name, outdir, use_new_file_format,
            max_level, chr_number_limit, length_limit, ensemble_size,
            selected_factors, logger);
}

void prepareGenerate() {
  std::string outdir = flag_set('o') ? get_arg('o') : "./synthetic/";
  std::string chr = flag_set('c') ? get_arg('c') : "chr22";
  std::string name = flag_set('n') ? get_arg('n') : "synthetic";
  int num_loops = flag_set('l') ? atoi(get_arg('l').c_str()) : 100;
  int ens = flag_set('m') ? atoi(get_arg('m').c_str()) : 20;

  // chr22 length for hg19
  int chr_length = 51304566;

  SyntheticGenerator gen;
  gen.setChromosome(chr, chr_length);
  gen.setResolution(25000);
  gen.setNumLoops(num_loops);
  gen.setEnsembleSize(ens);
  gen.setOutputDir(outdir);
  gen.setLabel(name);
  gen.generate();
}

void prepareBenchmark() {
  std::string outdir = flag_set('o') ? get_arg('o') : "./benchmark/";
  std::string chr = flag_set('c') ? get_arg('c') : "chr22";
  std::string name = flag_set('n') ? get_arg('n') : "bench";
  int ens = flag_set('m') ? atoi(get_arg('m').c_str()) : 20;

  int chr_length = 51304566;

  BenchmarkRunner bench;
  bench.setChromosome(chr, chr_length);
  bench.setResolution(25000);
  bench.setEnsembleSize(ens);
  bench.setOutputDir(outdir);
  bench.setLabel(name);
  bench.run();
}

void prepareDistMap() {
  if (!flag_set('i'))
    error_msg("No input structure specified [-i]\n");

  std::string input = get_arg('i');
  std::string outdir = flag_set('o') ? get_arg('o') : "./";
  std::string chr = flag_set('c') ? get_arg('c') : "";
  int level = flag_set('l') ? atoi(get_arg('l').c_str()) : 3;
  int resolution = flag_set('r') ? atoi(get_arg('r').c_str()) : 25000;

  portable_mkdir(outdir.c_str());

  HierarchicalChromosome hc;
  hc.fromFile(input);
  hc.setLevel(level);

  Chromosome flat = hc.createEqidistantModel(resolution, chr);
  if (flat.size < 2) {
    printf("Error: extracted chromosome has too few beads (%d)\n", flat.size);
    return;
  }

  DistanceMapGenerator gen;
  gen.setContactThreshold(2.0f);
  gen.setContactExponent(2.0f);

  std::string prefix = outdir + (flag_set('n') ? get_arg('n') : "structure");

  gen.writeDistanceMap(flat, prefix + "_distmap.heat");
  gen.writeContactMap(flat, prefix + "_contactmap.heat");

  Heatmap freq = gen.computeContactFrequencyMap(flat);
  freq.toFile(prefix + "_freqmap.heat", false);
  printf("Contact frequency map written to %s\n",
         (prefix + "_freqmap.heat").c_str());

  // Contact decay curve
  MetricsFramework metrics;
  auto decay = metrics.computeContactDecay(freq, flat.size, 1);
  metrics.writeContactDecay(decay, prefix + "_contact_decay.csv");
}

void prepareMetrics() {
  if (!flag_set('i'))
    error_msg("No input structures specified [-i] (comma-separated pair)\n");

  std::string input_str = get_arg('i');
  std::vector<std::string> inputs = split(input_str);
  if (inputs.size() < 2)
    error_msg("Need two input files separated by comma [-i file1,file2]\n");

  std::string outdir = flag_set('o') ? get_arg('o') : "./";
  std::string chr = flag_set('c') ? get_arg('c') : "";
  int level = flag_set('l') ? atoi(get_arg('l').c_str()) : 3;
  int resolution = flag_set('r') ? atoi(get_arg('r').c_str()) : 25000;

  portable_mkdir(outdir.c_str());

  HierarchicalChromosome hc_a, hc_b;
  hc_a.fromFile(inputs[0]);
  hc_b.fromFile(inputs[1]);

  hc_a.setLevel(level);
  hc_b.setLevel(level);

  Chromosome chr_a = hc_a.createEqidistantModel(resolution, chr);
  Chromosome chr_b = hc_b.createEqidistantModel(resolution, chr);

  if (chr_a.size < 2 || chr_b.size < 2) {
    printf("Error: insufficient beads for comparison (a=%d, b=%d)\n",
           chr_a.size, chr_b.size);
    return;
  }

  DistanceMapGenerator dmg;
  dmg.setContactExponent(2.0f);

  Heatmap heat_a = dmg.computeContactFrequencyMap(chr_a);
  Heatmap heat_b = dmg.computeContactFrequencyMap(chr_b);

  MetricsFramework metrics;

  std::string prefix = outdir + (flag_set('n') ? get_arg('n') : "comparison");

  metrics.writeMetricsReport(chr_a, chr_b, heat_a, heat_b,
                             prefix + "_metrics.log");

  // Per-diagonal correlation
  auto diag_corr = metrics.perDiagonalCorrelation(heat_a, heat_b);
  std::string diag_path = prefix + "_diagonal_corr.csv";
  FILE *f = fopen(diag_path.c_str(), "w");
  if (f) {
    fprintf(f, "diagonal_offset,correlation\n");
    for (const auto &dc : diag_corr)
      fprintf(f, "%d,%.6f\n", dc.first, dc.second);
    fclose(f);
    printf("Per-diagonal correlations written to %s\n", diag_path.c_str());
  }

  // SCC
  float scc = metrics.structuralSimilarity(heat_a, heat_b);
  printf("Structural similarity (SCC): %.4f\n", scc);
}

void prepareMerge() {
  if (!flag_set('B'))
    error_msg("No BED file specified [-B]\n");
  if (!flag_set('P'))
    error_msg("No BEDPE file specified [-P]\n");
  if (!flag_set('o'))
    error_msg("No output file specified [-o]\n");

  std::string bed_path = get_arg('B');
  std::string bedpe_path = get_arg('P');
  std::string output_path = get_arg('o');

  printf("Merging BED [%s] + BEDPE [%s] -> [%s]\n",
         bed_path.c_str(), bedpe_path.c_str(), output_path.c_str());

  if (!MergedFileIO::writeMerged(bed_path, bedpe_path, output_path))
    error_msg("Failed to create merged file\n");
}

void prepareSimulate() {
  if (!flag_set('i'))
    error_msg("No merged input file specified [-i]\n");
  if (!flag_set('s'))
    error_msg("No settings file specified [-s]\n");

  std::string merged_path = get_arg('i');
  std::string stg_path = get_arg('s');
  std::string outdir = flag_set('o') ? get_arg('o') : "./";

  // Load settings
  Settings stg;
  printf("Load settings from [%s]\n", stg_path.c_str());
  stg.loadFromINI(stg_path);

  // Apply CLI overrides
  if (flag_set('M'))
    Settings::maxMemoryMB = (size_t)atoi(get_arg('M').c_str());
  if (flag_set('E'))
    Settings::energyTraceEnabled = true;
  if (flag_set('L'))
    Settings::energyTraceInterval = atoi(get_arg('L').c_str());
  if (flag_set('I'))
    Settings::initMethod = get_arg('I');

  // Read merged file
  std::vector<MergedLoop> loops;
  std::vector<MergedAnchor> anchors;
  if (!MergedFileIO::readMerged(merged_path, loops, anchors))
    error_msg("Failed to read merged file\n");

  // Parse chromosomes from -c flag
  std::vector<std::string> chrs;
  std::string chromosomes = "genome";
  BedRegion region_of_interest;

  if (flag_set('c')) {
    chromosomes = get_arg('c');
    if (BedRegion::tryParse(chromosomes)) {
      region_of_interest.parse(chromosomes);
      chrs.push_back(region_of_interest.chr);
    } else {
      chrs = parseChromosomeDescription(chromosomes);
    }
    // Filter merged data by chromosome
    MergedFileIO::filterByChromosomes(loops, anchors, chrs);
  } else {
    chrs = parseChromosomeDescription(chromosomes);
  }

  // Create temp directory for extracted files
  std::string tmpdir = outdir + "simulate_tmp/";
  portable_mkdir(outdir.c_str());
  portable_mkdir(tmpdir.c_str());

  // Write temp BED and BEDPE
  std::string tmp_bed = tmpdir + "anchors.bed";
  std::string tmp_bedpe = tmpdir + "clusters.bedpe";

  MergedFileIO::writeBed(anchors, tmp_bed);
  MergedFileIO::writeBedpe(loops, tmp_bedpe);

  // Override Settings data paths to point to extracted files
  Settings::dataDirectory = "";
  Settings::dataAnchors = tmp_bed;
  Settings::dataPetClusters = tmp_bedpe;
  // Clear singletons since we only have loop data from merged file
  Settings::dataSingletons = "";
  Settings::dataSingletonsInter = "";
  Settings::dataFactors = "merged";

  std::string name = flag_set('n') ? get_arg('n') : "simulate";
  bool use_new_file_format = flag_set('z');
  int length_limit = flag_set('l') ? atoi(get_arg('l').c_str()) : -1;
  int chr_number_limit = flag_set('u') ? atoi(get_arg('u').c_str()) : -1;
  int max_level = flag_set('v') ? atoi(get_arg('v').c_str()) : 5;
  int ensemble_size = flag_set('m') ? atoi(get_arg('m').c_str()) : 1;

  // Save -F flag before args.clear()
  g_output_format = flag_set('F') ? get_arg('F') : "";

  // Create simulation logger
  portable_mkdir(outdir.c_str());
  unsigned int seed = 0;
  if (flag_set('j'))
    seed = (unsigned int)atoi(get_arg('j').c_str());
  Settings::gpuSeed = seed;

  SimulationLogger logger(outdir, name, seed);
  logger.logParameter("merged_input", merged_path);
  logger.logParameter("chromosomes", chromosomes);
  logger.logMilestone("Starting simulate reconstruction");

  args.clear();

  runLooper(chrs, region_of_interest, name, outdir, use_new_file_format,
            max_level, chr_number_limit, length_limit, ensemble_size,
            std::vector<std::string>(), &logger);

  // Run territory analysis on the result
  {
    std::string hcm_path = ftext("%sloops_%s.hcm", outdir.c_str(), name.c_str());
    if (file_exists(hcm_path)) {
      HierarchicalChromosome hc;
      hc.fromFile(hcm_path);
      GenomeReconstructor gr;
      gr.analyze(hc, 1);
      gr.writeReport(outdir, name);
    }
  }

  // Final output validation (also called inside runLooper, but re-check after
  // territory analysis in case post-processing introduced issues)
  validateOutputsAndWriteErrorWarning(outdir, name, use_new_file_format,
                                      ensemble_size);

  logger.logMilestone("Simulate reconstruction complete");
  logger.writeSummary();
}

void prepareIbed() {
  if (!flag_set('i'))
    error_msg("No ibed input file specified [-i]\n");
  if (!flag_set('s'))
    error_msg("No settings file specified [-s]\n");

  std::string ibed_path = get_arg('i');
  std::string stg_path = get_arg('s');
  std::string outdir = flag_set('o') ? get_arg('o') : "./";

  // Load settings
  Settings stg;
  printf("Load settings from [%s]\n", stg_path.c_str());
  stg.loadFromINI(stg_path);

  // Apply CLI overrides
  if (flag_set('M'))
    Settings::maxMemoryMB = (size_t)atoi(get_arg('M').c_str());
  if (flag_set('E'))
    Settings::energyTraceEnabled = true;
  if (flag_set('L'))
    Settings::energyTraceInterval = atoi(get_arg('L').c_str());
  if (flag_set('I'))
    Settings::initMethod = get_arg('I');

  // Read ibed file
  std::vector<MergedLoop> loops;
  std::vector<MergedAnchor> anchors;

  // Parse chromosomes from -c flag
  std::vector<std::string> chrs;
  std::string chromosomes = "genome";
  BedRegion region_of_interest;

  if (flag_set('c')) {
    chromosomes = get_arg('c');
    if (BedRegion::tryParse(chromosomes)) {
      region_of_interest.parse(chromosomes);
      chrs.push_back(region_of_interest.chr);
    } else {
      chrs = parseChromosomeDescription(chromosomes);
    }
    // Read and filter
    if (!IbedFileIO::readIbedFiltered(ibed_path, loops, anchors, chrs))
      error_msg("Failed to read ibed file\n");
  } else {
    if (!IbedFileIO::readIbed(ibed_path, loops, anchors))
      error_msg("Failed to read ibed file\n");
    chrs = parseChromosomeDescription(chromosomes);
  }

  printf("[ibed] Loaded %d loops, %d anchors from ibed\n",
         (int)loops.size(), (int)anchors.size());

  // Create temp directory for extracted files
  std::string tmpdir = outdir + "ibed_tmp/";
  portable_mkdir(outdir.c_str());
  portable_mkdir(tmpdir.c_str());

  // Write temp BED and BEDPE
  std::string tmp_bed = tmpdir + "anchors.bed";
  std::string tmp_bedpe = tmpdir + "clusters.bedpe";

  MergedFileIO::writeBed(anchors, tmp_bed);
  MergedFileIO::writeBedpe(loops, tmp_bedpe);

  // Override Settings data paths
  Settings::dataDirectory = "";
  Settings::dataAnchors = tmp_bed;
  Settings::dataPetClusters = tmp_bedpe;
  Settings::dataSingletons = "";
  Settings::dataSingletonsInter = "";
  Settings::dataFactors = "ibed";

  std::string name = flag_set('n') ? get_arg('n') : "ibed";
  bool use_new_file_format = flag_set('z');
  int length_limit = flag_set('l') ? atoi(get_arg('l').c_str()) : -1;
  int chr_number_limit = flag_set('u') ? atoi(get_arg('u').c_str()) : -1;
  int max_level = flag_set('v') ? atoi(get_arg('v').c_str()) : 5;
  int ensemble_size = flag_set('m') ? atoi(get_arg('m').c_str()) : 1;

  g_output_format = flag_set('F') ? get_arg('F') : "";

  // Create simulation logger
  portable_mkdir(outdir.c_str());
  unsigned int seed = 0;
  if (flag_set('j'))
    seed = (unsigned int)atoi(get_arg('j').c_str());
  Settings::gpuSeed = seed;

  SimulationLogger logger(outdir, name, seed);
  logger.logParameter("ibed_input", ibed_path);
  logger.logParameter("chromosomes", chromosomes);
  logger.logMilestone("Starting ibed reconstruction");

  args.clear();

  runLooper(chrs, region_of_interest, name, outdir, use_new_file_format,
            max_level, chr_number_limit, length_limit, ensemble_size,
            std::vector<std::string>(), &logger);

  // Run territory analysis
  {
    std::string hcm_path = ftext("%sloops_%s.hcm", outdir.c_str(), name.c_str());
    if (file_exists(hcm_path)) {
      HierarchicalChromosome hc;
      hc.fromFile(hcm_path);
      GenomeReconstructor gr;
      gr.analyze(hc, 1);
      gr.writeReport(outdir, name);
    }
  }

  // Final output validation (also called inside runLooper, but re-check after
  // territory analysis in case post-processing introduced issues)
  validateOutputsAndWriteErrorWarning(outdir, name, use_new_file_format,
                                      ensemble_size);

  logger.logMilestone("ibed reconstruction complete");
  logger.writeSummary();
}

// Convert HCM file to PDB and/or CIF format.
// Usage: pMMC -a convert -i <input.hcm> -o <output_dir/> [-F cif|pdb|both]
void prepareConvert() {
  if (!flag_set('i'))
    error_msg("No input HCM file specified [-i]\n");

  string input = get_arg('i');
  string outdir = flag_set('o') ? get_arg('o') : "./";
  string format = flag_set('F') ? get_arg('F') : "both";

  printf("[convert] Reading HCM: %s\n", input.c_str());
  HierarchicalChromosome hc;
  hc.fromFile(input);

  // Ensure current_level is populated (fromFilePreviousFormat may not call useTopLevel)
  hc.useTopLevel();
  hc.useLowestLevel();

  // Derive base name from input filename (strip path and .hcm extension)
  string base = input;
  size_t slash = base.find_last_of("/\\");
  if (slash != string::npos)
    base = base.substr(slash + 1);
  if (base.size() > 4 && base.substr(base.size() - 4) == ".hcm")
    base = base.substr(0, base.size() - 4);

  string outbase = outdir + base;

  if (format == "cif" || format == "both") {
    string cif_path = outbase + ".cif";
    CifWriter::write(hc, cif_path, "cudaMMC", false);
    printf("[convert] Wrote CIF: %s\n", cif_path.c_str());
  }
  if (format == "pdb" || format == "both") {
    string pdb_path = outbase + ".pdb";
    PdbWriter::write(hc, pdb_path, false);
    printf("[convert] Wrote PDB: %s\n", pdb_path.c_str());
  }
}

int main(int argc, char **argv) {
  setbuf(stdout, NULL);

  opterr = 0;
  int opt = 0;

  while ((opt = getopt(argc, argv,
                       "s:a:o:c:n:l:u:t:d:p:b:e:m:i:v:zf:hr:g:x:j:M:F:EL:I:B:P:")) !=
         EOF) {
    if (optarg == NULL)
      optarg = (char *)"";
    printf("opt = [%c %s]\n", opt, optarg);
    args[opt] = optarg;
  }

  // Seed handling: -j for deterministic seed, otherwise time-based
  unsigned int seed;
  if (flag_set('j')) {
    seed = (unsigned int)atoi(get_arg('j').c_str());
    printf("Using fixed random seed: %u\n", seed);
  } else {
    auto tmp = chrono::high_resolution_clock::now().time_since_epoch().count();
    seed = (unsigned int)tmp;
    printf("Using time-based random seed: %u\n", seed);
  }
  srand(seed);
  setGpuSeed(seed); // propagate CLI seed to GPU curand initialization

  if (!flag_set('a') || args['a'] == "c" || args['a'] == "create") {

    if (!flag_set('s')) {
      printf("no settings file specified!\n");
      usage();
    }

    Settings stg;
    string stg_path = get_arg('s');
    printf("Load settings from [%s]\n", stg_path.c_str());
    stg.loadFromINI(stg_path);

    // Apply CLI overrides
    if (flag_set('M'))
      Settings::maxMemoryMB = (size_t)atoi(get_arg('M').c_str());
    if (flag_set('E'))
      Settings::energyTraceEnabled = true;
    if (flag_set('L'))
      Settings::energyTraceInterval = atoi(get_arg('L').c_str());
    if (flag_set('I'))
      Settings::initMethod = get_arg('I');
    Settings::gpuSeed = seed;
    if (Settings::maxMemoryMB > 0)
      printf("Memory budget: %zuMB\n", Settings::maxMemoryMB);

    logMemoryUsage("after settings load");

    // Create simulation logger for reproducibility
    string outdir = flag_set('o') ? get_arg('o') : "./";
    string label = flag_set('n') ? get_arg('n') : "run";
    portable_mkdir(outdir.c_str());
    SimulationLogger logger(outdir, label, seed);
    logger.logParameter("settings_file", stg_path);
    logger.logParameter("chromosomes", flag_set('c') ? get_arg('c') : "genome");
    logger.logMilestone("Starting reconstruction");

    // Capture label/outdir before prepareLooper clears args
    string ter_label = label;
    string ter_outdir = outdir;

    prepareLooper(&logger);

    // Run territory analysis on the reconstructed model
    {
      string hcm_path = ftext("%sloops_%s.hcm", ter_outdir.c_str(), ter_label.c_str());
      if (file_exists(hcm_path)) {
        HierarchicalChromosome hc_ter;
        hc_ter.fromFile(hcm_path);
        GenomeReconstructor gr;
        gr.analyze(hc_ter, 1);
        gr.writeReport(ter_outdir, ter_label);
      }
    }

    logMemoryUsage("after reconstruction");
    logger.logMilestone("Reconstruction complete");
    logger.writeSummary();
  }

  else if (args['a'] == "e" || args['a'] == "ensemble")
    prepareEnsembleAnalysis();
  else if (args['a'] == "r" || args['a'] == "extract")
    prepareExtract();
  else if (args['a'] == "p" || args['a'] == "position")
    preparePosition();
  else if (args['a'] == "s" || args['a'] == "smooth")
    prepareSmooth();
  else if (args['a'] == "d" || args['a'] == "distance")
    prepareCalcDistanceMatrix();
  else if (args['a'] == "f" || args['a'] == "flatten")
    prepareFlatten();
  else if (args['a'] == "w" || args['a'] == "rewiring")
    prepareRewiring();
  else if (args['a'] == "g" || args['a'] == "generate")
    prepareGenerate();
  else if (args['a'] == "k" || args['a'] == "benchmark")
    prepareBenchmark();
  else if (args['a'] == "distmap")
    prepareDistMap();
  else if (args['a'] == "metrics")
    prepareMetrics();
  else if (args['a'] == "M" || args['a'] == "merge")
    prepareMerge();
  else if (args['a'] == "S" || args['a'] == "simulate")
    prepareSimulate();
  else if (args['a'] == "I" || args['a'] == "ibed")
    prepareIbed();
  else if (args['a'] == "X" || args['a'] == "convert")
    prepareConvert();
  else {
    printf("No matching action specified [%s]!\n", args['a'].c_str());
    usage();
  }

  args.clear();
  printf("end\n");
  return 0;
}
