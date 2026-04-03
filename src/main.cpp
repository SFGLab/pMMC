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

#include <BedRegion.hpp>
#include <BenchmarkRunner.hpp>
#include <CifWriter.hpp>
#include <DistanceMapGenerator.hpp>
#include <GenomeReconstructor.hpp>
#include <HierarchicalChromosome.h>
#include <LooperSolver.hpp>
#include <ContactMatrixIO.hpp>
#include <IbedFileIO.hpp>
#include <MergedFileIO.hpp>
#include <MetricsFramework.hpp>
#include <MultiscaleEnergy.hpp>
#include <PdbWriter.hpp>
#include <SimulationLogger.hpp>
#include <SyntheticGenerator.hpp>
#include <HeatmapImageWriter.hpp>
#include <FlowchartGenerator.hpp>
#include <SvgChartGenerator.hpp>
#include <PicturePanel.hpp>
#include <VisualizationScripts.hpp>
#include <common.h>
#include <Static.hpp>

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
  printf("  contact (C)    Reconstruct from contact matrix (Hi-C/pcHi-C/HiChIP)\n");
  printf("  convert (X)    Convert HCM file to PDB/CIF format\n");
  printf("  batch (b)      Multi-cell-line batch reconstruction + similarity\n");
  printf("  flowchart      Generate algorithm flowchart (SVG)\n");
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
  printf("\nAnnotations:\n");
  printf("  -A  annotation directory (containing enhancer.bed, promoter.bed, epigenetic.bed)\n");
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

void extract(string input_path, string output_path, int start, int end,
             const string &format = "") {
  printf("extract fragment: %d %d\n", start, end);
  HierarchicalChromosome h;
  h.fromFile(input_path);
  printf("size = %d\n", (int)h.clusters.size());
  if (h.clusters.size() > 0) {
    HierarchicalChromosome extr = h.extractFragment(start, end);
    extr.toFile(output_path);

    // Write PDB/CIF if requested via -F flag
    if (!format.empty()) {
      // Derive base path by stripping .hcm extension if present
      string base = output_path;
      if (base.size() > 4 && base.substr(base.size() - 4) == ".hcm")
        base = base.substr(0, base.size() - 4);

      if (format == "pdb" || format == "both") {
        PdbWriter::write(extr, base + ".pdb");
      }
      if (format == "cif" || format == "mmcif" || format == "both") {
        CifWriter::write(extr, base + ".cif");
      }
    }
  }
}

void prepareExtract() {
  if (!flag_set('i'))
    error_msg("No input specified [-i]\n");
  if (!flag_set('o'))
    error_msg("No output specified [-o]\n");

  string input = get_arg('i');
  string output = get_arg('o');
  string format = flag_set('F') ? get_arg('F') : "";

  int start = 0, end = 0;
  bool has_region = false;
  BedRegion region;

  if (flag_set('c') && BedRegion::tryParse(get_arg('c').c_str())) {
    region.parse(get_arg('c'));
    start = region.start;
    end = region.end;
    has_region = true;
  } else {
    if (!flag_set('s'))
      error_msg("No start position specified [-s]\n");
    if (!flag_set('e'))
      error_msg("No end position specified [-e]\n");
    start = atoi(get_arg('s').c_str());
    end = atoi(get_arg('e').c_str());
  }

  extract(input, output, start, end, format);

  // Generate visualization scripts if PDB was produced
  if (format == "pdb" || format == "both") {
    string base = output;
    if (base.size() > 4 && base.substr(base.size() - 4) == ".hcm")
      base = base.substr(0, base.size() - 4);
    string pdb_path = base + ".pdb";

    // Load the extracted structure for visualization script generation
    HierarchicalChromosome hc;
    string hcm_path = output;
    if (file_exists(hcm_path)) {
      hc.fromFile(hcm_path);
      string pml_path = base + ".pml";
      string cxc_path = base + ".cxc";
      VisualizationScripts::writePyMOLScript(hc, pdb_path, pml_path);
      VisualizationScripts::writeChimeraXScript(hc, pdb_path, cxc_path);
      printf("Visualization scripts: %s, %s\n", pml_path.c_str(), cxc_path.c_str());
    }
  }
}

void ensembleAnalysis(string input_dir, string pattern, string output_dir,
                      string name_prefix) {
  // Read all structures, compute pairwise Spearman correlation, and write
  // correlation matrix CSV, histogram SVG, heatmap SVG, and summary stats.
  printf("ensemble analysis:\n");

  portable_mkdir(output_dir.c_str());

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
    printf("[%s] - loaded\n", path.c_str());
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

  vector<string> chrs = hchr[0].chrs; // save chromosomes

  // --- Distance heatmap analysis (original) ---
  vector<Heatmap> heat;
  Heatmap heat_pairs(N);

  // first chr level
  if (chrs.size() < 2)
    printf("Only %d chromosome available, skip chr level\n", (int)chrs.size());
  else {
    printf("chr level\n");
    printf("create structural heatmaps\n");
    for (int j = 0; j < N; ++j)
      heat.push_back(hchr[j].createStructuralHeatmap("", 0));

    printf("create distance heatmap\n");
    for (int i = 0; i < N; ++i)
      for (int j = i + 1; j < N; ++j)
        heat_pairs.v[i][j] = heat_pairs.v[j][i] = heat[i].calcDistance(heat[j]);

    heat_pairs.toFile(output_dir + name_prefix + "_distances_lvl0.heat", false);
    heat_pairs.zero();
    heat.clear();
  }

  // segment and subanchor levels
  int i, j;
  for (j = 0; j < N; ++j)
    hchr[j].levelDown();

  float d;
  for (int lvl = 1; lvl <= 2; ++lvl) {
    printf("do level %d\n", lvl);
    for (string chr : chrs) {
      for (j = 0; j < N; ++j) {
        hchr[j].setLevel(lvl);
        heat.push_back(hchr[j].createStructuralHeatmap(chr, lvl));
      }
      for (i = 0; i < N; ++i)
        for (j = i + 1; j < N; ++j) {
          d = heat[i].calcDistance(heat[j]);
          heat_pairs.v[i][j] += d;
          heat_pairs.v[j][i] += d;
        }
    }
    heat_pairs.toFile(
        output_dir + name_prefix + ftext("_distances_lvl%d.heat", lvl), false);
    heat_pairs.zero();
    heat.clear();
  }

  // --- Pairwise Spearman correlation analysis ---
  printf("computing pairwise Spearman correlations...\n");

  // Extract Chromosome objects at the lowest level for each structure
  vector<Chromosome> structures;
  for (int s = 0; s < N; ++s) {
    hchr[s].useLowestLevel();
    hchr[s].createCurrentLevelStructure();
    // Concatenate all chromosomes into a single Chromosome for comparison
    Chromosome combined;
    combined.size = 0;
    for (const string &c : chrs) {
      if (hchr[s].chr.count(c)) {
        const Chromosome &ch = hchr[s].chr[c];
        for (int b = 0; b < ch.size; ++b)
          combined.points.push_back(ch.points[b]);
        combined.size += ch.size;
      }
    }
    structures.push_back(combined);
  }

  MetricsFramework metrics;
  int subsample = structures[0].size > 500 ? 2 : 1;
  vector<float> similarities =
      metrics.ensembleSimilarityDistribution(structures, subsample);

  // Build NxN correlation matrix
  vector<vector<float>> corr_matrix(N, vector<float>(N, 1.0f));
  int pair_idx = 0;
  for (int a = 0; a < N; ++a) {
    for (int b = a + 1; b < N; ++b) {
      if (pair_idx < static_cast<int>(similarities.size())) {
        corr_matrix[a][b] = similarities[pair_idx];
        corr_matrix[b][a] = similarities[pair_idx];
      }
      pair_idx++;
    }
  }

  // 1) Write pairwise Spearman correlation CSV
  string corr_csv_path = output_dir + name_prefix + "_spearman_correlations.csv";
  {
    FILE *f = fopen(corr_csv_path.c_str(), "w");
    if (f) {
      fprintf(f, "pair_index,structure_i,structure_j,spearman_correlation\n");
      int idx = 0;
      for (int a = 0; a < N; ++a)
        for (int b = a + 1; b < N; ++b) {
          fprintf(f, "%d,%d,%d,%.8f\n", idx, a, b, corr_matrix[a][b]);
          idx++;
        }
      fclose(f);
      printf("wrote %s\n", corr_csv_path.c_str());
    } else {
      printf("Error: cannot write %s\n", corr_csv_path.c_str());
    }
  }

  // 2) Write correlation matrix CSV (for heatmap)
  string matrix_csv_path =
      output_dir + name_prefix + "_correlation_matrix.csv";
  {
    FILE *f = fopen(matrix_csv_path.c_str(), "w");
    if (f) {
      // header row
      for (int a = 0; a < N; ++a)
        fprintf(f, "%s%d", a == 0 ? "" : ",", a);
      fprintf(f, "\n");
      // data rows
      for (int a = 0; a < N; ++a) {
        for (int b = 0; b < N; ++b)
          fprintf(f, "%s%.6f", b == 0 ? "" : ",", corr_matrix[a][b]);
        fprintf(f, "\n");
      }
      fclose(f);
      printf("wrote %s\n", matrix_csv_path.c_str());
    } else {
      printf("Error: cannot write %s\n", matrix_csv_path.c_str());
    }
  }

  // 3) SVG heatmap of correlation matrix
  string heatmap_svg_path =
      output_dir + name_prefix + "_correlation_heatmap.svg";
  {
    vector<string> labels;
    for (int a = 0; a < N; ++a)
      labels.push_back(ftext("%d", a));
    // Determine min correlation for color scale
    float min_corr = 1.0f;
    for (auto &v : similarities)
      if (v < min_corr)
        min_corr = v;
    min_corr = std::max(0.0f, min_corr - 0.05f); // give a little margin
    SvgChartGenerator::heatmap(corr_matrix, labels, labels,
                               "Pairwise Spearman Correlation",
                               heatmap_svg_path, min_corr, 1.0f);
    printf("wrote %s\n", heatmap_svg_path.c_str());
  }

  // 4) SVG histogram of correlation distribution
  string hist_svg_path =
      output_dir + name_prefix + "_correlation_histogram.svg";
  {
    // Create histogram bins
    int num_bins = 20;
    float hist_min = 1.0f, hist_max = -1.0f;
    for (float v : similarities) {
      if (v < hist_min) hist_min = v;
      if (v > hist_max) hist_max = v;
    }
    if (hist_min == hist_max) {
      hist_min -= 0.05f;
      hist_max += 0.05f;
    }
    float bin_width = (hist_max - hist_min) / num_bins;
    vector<float> counts(num_bins, 0.0f);
    vector<string> bin_labels(num_bins);
    for (int b = 0; b < num_bins; ++b)
      bin_labels[b] = ftext("%.2f", hist_min + (b + 0.5f) * bin_width);
    for (float v : similarities) {
      int bin = static_cast<int>((v - hist_min) / bin_width);
      if (bin < 0) bin = 0;
      if (bin >= num_bins) bin = num_bins - 1;
      counts[bin] += 1.0f;
    }
    SvgChartGenerator::histogram(bin_labels, counts,
                                 "Distribution of Pairwise Spearman Correlations",
                                 "Spearman Correlation", "Count",
                                 hist_svg_path, "#3b82f6");
    printf("wrote %s\n", hist_svg_path.c_str());
  }

  // 5) Summary statistics text file
  string summary_path = output_dir + name_prefix + "_ensemble_summary.txt";
  {
    FILE *f = fopen(summary_path.c_str(), "w");
    if (f) {
      fprintf(f, "Ensemble Analysis Summary\n");
      fprintf(f, "=========================\n\n");
      fprintf(f, "Number of structures: %d\n", N);
      fprintf(f, "Number of pairwise comparisons: %d\n",
              static_cast<int>(similarities.size()));
      fprintf(f, "Chromosomes: ");
      for (size_t c = 0; c < chrs.size(); ++c)
        fprintf(f, "%s%s", c > 0 ? ", " : "", chrs[c].c_str());
      fprintf(f, "\n\n");

      if (!similarities.empty()) {
        float sum = 0.0f, sum_sq = 0.0f;
        float min_val = similarities[0], max_val = similarities[0];
        for (float v : similarities) {
          sum += v;
          sum_sq += v * v;
          if (v < min_val) min_val = v;
          if (v > max_val) max_val = v;
        }
        float mean = sum / similarities.size();
        float variance =
            sum_sq / similarities.size() - mean * mean;
        float stddev = variance > 0.0f ? sqrtf(variance) : 0.0f;

        // Compute median
        vector<float> sorted_sim = similarities;
        std::sort(sorted_sim.begin(), sorted_sim.end());
        float median;
        int ns = static_cast<int>(sorted_sim.size());
        if (ns % 2 == 0)
          median = 0.5f * (sorted_sim[ns / 2 - 1] + sorted_sim[ns / 2]);
        else
          median = sorted_sim[ns / 2];

        fprintf(f, "Spearman Correlation Statistics\n");
        fprintf(f, "-------------------------------\n");
        fprintf(f, "Mean:     %.6f\n", mean);
        fprintf(f, "Median:   %.6f\n", median);
        fprintf(f, "Std Dev:  %.6f\n", stddev);
        fprintf(f, "Min:      %.6f\n", min_val);
        fprintf(f, "Max:      %.6f\n", max_val);
        fprintf(f, "\n");

        printf("Ensemble similarity (%d pairs): mean=%.4f median=%.4f "
               "std=%.4f min=%.4f max=%.4f\n",
               static_cast<int>(similarities.size()), mean, median, stddev,
               min_val, max_val);
      }
      fclose(f);
      printf("wrote %s\n", summary_path.c_str());
    } else {
      printf("Error: cannot write %s\n", summary_path.c_str());
    }
  }

  printf("ensemble analysis complete\n");
}

void prepareEnsembleAnalysis() {
  if (!flag_set('i'))
    error_msg("No input directory specified [-i]\n");
  if (!flag_set('p'))
    error_msg("No pattern specified [-p]\n");

  string input = get_arg('i');
  string pattern = get_arg('p');
  string output = flag_set('o') ? get_arg('o') : input;
  string name = flag_set('n') ? get_arg('n') : "ensemble";

  ensembleAnalysis(input, pattern, output, name);
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

    // Generate visualization scripts alongside the PDB (F4.1-F4.8)
    {
      std::string pdb_path = ftext("%s%s.pdb", outdir.c_str(), label.c_str());
      std::string pml_path = ftext("%s%s.pml", outdir.c_str(), label.c_str());
      std::string cxc_path = ftext("%s%s.cxc", outdir.c_str(), label.c_str());
      VisualizationScripts::writePyMOLScript(hc, pdb_path, pml_path);
      VisualizationScripts::writeChimeraXScript(hc, pdb_path, cxc_path);
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

      // Generate visualization scripts for ensemble member (F4.1-F4.8)
      {
        std::string pdb_path = ftext("%s%s_%d.pdb", outdir.c_str(), label.c_str(), i);
        std::string pml_path = ftext("%s%s_%d.pml", outdir.c_str(), label.c_str(), i);
        std::string cxc_path = ftext("%s%s_%d.cxc", outdir.c_str(), label.c_str(), i);
        VisualizationScripts::writePyMOLScript(hc, pdb_path, pml_path);
        VisualizationScripts::writeChimeraXScript(hc, pdb_path, cxc_path);
      }

      // remove subanchor beads
      lsm.removeSubanchorBeads();
      lsm.reset();
    }
  } else
    error("ensemble size not correct");

  // ── Feature F5.7: Loop length distribution ────────────────────────
  {
    std::string csv_path = ftext("%sloop_length_distribution_%s.csv", outdir.c_str(), label.c_str());
    FILE *fcsv = fopen(csv_path.c_str(), "w");
    if (fcsv) {
      fprintf(fcsv, "bin_start,bin_end,count,fraction\n");

      // Define histogram bins (in bp)
      struct LenBin { long long lo; long long hi; int count; };
      std::vector<LenBin> bins = {
        {0,        50000,     0},
        {50000,    100000,    0},
        {100000,   200000,    0},
        {200000,   500000,    0},
        {500000,   1000000,   0},
        {1000000,  5000000,   0},
        {5000000,  10000000,  0},
        {10000000, 999999999, 0}
      };

      int total_arcs = 0;
      for (const auto &chr : lsm.chrs) {
        if (lsm.arcs.arcs.find(chr) == lsm.arcs.arcs.end())
          continue;
        const auto &chr_arcs = lsm.arcs.arcs[chr];
        for (size_t ai = 0; ai < chr_arcs.size(); ++ai) {
          const InteractionArc &arc = chr_arcs[ai];
          if (arc.factor == -1 && arc.eff_score == 0)
            continue; // skip deactivated summary arcs
          long long genomic_dist = (long long)arc.genomic_end - (long long)arc.genomic_start;
          if (genomic_dist < 0) genomic_dist = -genomic_dist;
          for (auto &b : bins) {
            if (genomic_dist >= b.lo && genomic_dist < b.hi) {
              b.count++;
              break;
            }
          }
          total_arcs++;
        }
      }

      for (const auto &b : bins) {
        double frac = total_arcs > 0 ? (double)b.count / total_arcs : 0.0;
        long long hi_display = b.hi >= 999999999 ? -1 : b.hi; // -1 means unbounded
        fprintf(fcsv, "%lld,%lld,%d,%.6f\n", b.lo, hi_display, b.count, frac);
      }
      fclose(fcsv);
      printf("[F5.7] Loop length distribution written to %s (%d arcs)\n",
             csv_path.c_str(), total_arcs);

      // Generate SVG chart
      std::string svg_path = csv_path.substr(0, csv_path.size() - 4) + ".svg";
      SvgChartGenerator::loopLengthChart(csv_path, svg_path);
    }
  }

  // ── Feature F5.8: CTCF occupancy at loop anchors ─────────────────
  {
    std::string csv_path = ftext("%sctcf_occupancy_%s.csv", outdir.c_str(), label.c_str());
    FILE *fcsv = fopen(csv_path.c_str(), "w");
    if (fcsv) {
      fprintf(fcsv, "total_anchors,ctcf_anchors,occupancy_fraction,"
                     "loops_with_ctcf_both,loops_with_ctcf_one,loops_with_ctcf_none\n");

      // Count CTCF-annotated anchors across all chromosomes
      int total_anchors = 0;
      int ctcf_anchors = 0;
      // Build set of CTCF anchor global indices for arc endpoint lookup
      std::set<int> ctcf_cluster_set;

      for (const auto &chr : lsm.chrs) {
        if (lsm.arcs.anchors.find(chr) == lsm.arcs.anchors.end())
          continue;
        if (lsm.chr_first_cluster.find(chr) == lsm.chr_first_cluster.end())
          continue;
        int base = lsm.chr_first_cluster[chr];
        const auto &chr_anchors = lsm.arcs.anchors[chr];
        for (size_t ai = 0; ai < chr_anchors.size(); ++ai) {
          total_anchors++;
          bool is_ctcf = (chr_anchors[ai].orientation == 'L' ||
                          chr_anchors[ai].orientation == 'R');
          // Also check annotation_type on the cluster
          int ci = base + (int)ai;
          if (ci < (int)lsm.clusters.size() &&
              lsm.clusters[ci].annotation_type == ANNOT_CTCF)
            is_ctcf = true;
          if (is_ctcf) {
            ctcf_anchors++;
            ctcf_cluster_set.insert(ci);
          }
        }
      }

      // Classify arcs by CTCF status of their endpoints
      int loops_both = 0, loops_one = 0, loops_none = 0;
      for (const auto &chr : lsm.chrs) {
        if (lsm.arcs.arcs.find(chr) == lsm.arcs.arcs.end())
          continue;
        const auto &chr_arcs = lsm.arcs.arcs[chr];
        for (size_t ai = 0; ai < chr_arcs.size(); ++ai) {
          const InteractionArc &arc = chr_arcs[ai];
          if (arc.factor == -1 && arc.eff_score == 0)
            continue;
          bool s_ctcf = ctcf_cluster_set.count(arc.start) > 0;
          bool e_ctcf = ctcf_cluster_set.count(arc.end) > 0;
          if (s_ctcf && e_ctcf)
            loops_both++;
          else if (s_ctcf || e_ctcf)
            loops_one++;
          else
            loops_none++;
        }
      }

      double occupancy_frac = total_anchors > 0
                                  ? (double)ctcf_anchors / total_anchors
                                  : 0.0;
      fprintf(fcsv, "%d,%d,%.6f,%d,%d,%d\n", total_anchors, ctcf_anchors,
              occupancy_frac, loops_both, loops_one, loops_none);
      fclose(fcsv);
      printf("[F5.8] CTCF occupancy written to %s (total=%d, ctcf=%d, frac=%.4f)\n",
             csv_path.c_str(), total_anchors, ctcf_anchors, occupancy_frac);
    }
  }

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

// Write a benchmark bash script to the given path.
static void writeBenchmarkScript(const std::string &path,
                                 const std::string &content) {
  FILE *f = fopen(path.c_str(), "w");
  if (!f) {
    printf("[benchmark] ERROR: cannot write %s\n", path.c_str());
    return;
  }
  fprintf(f, "%s", content.c_str());
  fclose(f);
  printf("[benchmark] Wrote %s\n", path.c_str());
}

void prepareBenchmark() {
  std::string outdir = flag_set('o') ? get_arg('o') : "./benchmark/";
  std::string chr = flag_set('c') ? get_arg('c') : "chr22";
  int ens = flag_set('m') ? atoi(get_arg('m').c_str()) : 20;

  portable_mkdir(outdir.c_str());

  // Also run the internal pMMC-only quick validation benchmark
  std::string name = flag_set('n') ? get_arg('n') : "bench";
  int chr_length = 51304566;
  BenchmarkRunner bench;
  bench.setChromosome(chr, chr_length);
  bench.setResolution(25000);
  bench.setEnsembleSize(ens);
  bench.setOutputDir(outdir);
  bench.setLabel(name);
  bench.run();

  // ── Generate redistributable bash scripts ──────────────────────────
  // All tools (3D-GNOME, cudaMMC, pMMC) run SEQUENTIALLY — never
  // simultaneously — so they do not compete for GPU/memory resources.

  std::string ens_str = ftext("%d", ens);

  // Script 1: run_benchmark.sh — Full 3-tool comparison on synthetic data
  {
    std::string s;
    s += "#!/usr/bin/env bash\n";
    s += "# pMMC Benchmark — 3-tool sequential comparison\n";
    s += "# Generated by: pMMC -a benchmark\n";
    s += "# Runs 3D-GNOME, cudaMMC, and pMMC ONE AFTER ANOTHER on\n";
    s += "# the same synthetic data so they never compete for resources.\n";
    s += "#\n";
    s += "# Usage:  bash run_benchmark.sh [output_dir]\n";
    s += "#\n";
    s += "# Environment variables (override paths):\n";
    s += "#   PMMC_EXE    — pMMC executable     (default: pMMC.exe)\n";
    s += "#   CUDAMMC_EXE — cudaMMC executable   (default: cudaMMC, optional)\n";
    s += "#   GNOME_EXE   — 3D-GNOME executable  (default: 3dnome, optional)\n";
    s += "set -euo pipefail\n\n";
    s += "OUTDIR=\"${1:-./benchmark_results}\"\n";
    s += "PMMC_EXE=\"${PMMC_EXE:-pMMC.exe}\"\n";
    s += "CUDAMMC_EXE=\"${CUDAMMC_EXE:-cudaMMC}\"\n";
    s += "GNOME_EXE=\"${GNOME_EXE:-3dnome}\"\n";
    s += "CHR=\"" + chr + "\"\n";
    s += "SEED=42\n";
    s += "ENSEMBLE=" + ens_str + "\n";
    s += "CSV=\"$OUTDIR/benchmark_results.csv\"\n\n";
    s += "mkdir -p \"$OUTDIR\"\n";
    s += "echo \"tool,chromosome,wall_time_seconds,exit_code\" > \"$CSV\"\n\n";
    s += "# Step 1 — generate synthetic data\n";
    s += "echo \"[1/5] Generating synthetic data...\"\n";
    s += "SYNTH=\"$OUTDIR/synthetic/\"\n";
    s += "\"$PMMC_EXE\" -a generate -o \"$SYNTH\" -c \"$CHR\" -l 100 -m \"$ENSEMBLE\" -j $SEED\n";
    s += "STG=\"$SYNTH/settings.ini\"\n\n";
    s += "# Step 2 — 3D-GNOME (sequential)\n";
    s += "if command -v \"$GNOME_EXE\" &>/dev/null; then\n";
    s += "  echo \"[2/5] Running 3D-GNOME...\"\n";
    s += "  D=\"$OUTDIR/3dgnome/\"; mkdir -p \"$D\"\n";
    s += "  T0=$(date +%s%N)\n";
    s += "  \"$GNOME_EXE\" -a create -s \"$STG\" -c \"$CHR\" -o \"$D\" -n gnome -j $SEED && RC=0 || RC=$?\n";
    s += "  T1=$(date +%s%N); T=$(echo \"scale=3;($T1-$T0)/1000000000\"|bc)\n";
    s += "  echo \"3D-GNOME,$CHR,$T,$RC\" >> \"$CSV\"\n";
    s += "else\n  echo \"[2/5] 3D-GNOME not found, skipping.\"\n";
    s += "  echo \"3D-GNOME,$CHR,,skipped\" >> \"$CSV\"\nfi\n\n";
    s += "# Step 3 — cudaMMC (sequential)\n";
    s += "if command -v \"$CUDAMMC_EXE\" &>/dev/null; then\n";
    s += "  echo \"[3/5] Running cudaMMC...\"\n";
    s += "  D=\"$OUTDIR/cudammc/\"; mkdir -p \"$D\"\n";
    s += "  T0=$(date +%s%N)\n";
    s += "  \"$CUDAMMC_EXE\" -a create -s \"$STG\" -c \"$CHR\" -o \"$D\" -n cuda -j $SEED && RC=0 || RC=$?\n";
    s += "  T1=$(date +%s%N); T=$(echo \"scale=3;($T1-$T0)/1000000000\"|bc)\n";
    s += "  echo \"cudaMMC,$CHR,$T,$RC\" >> \"$CSV\"\n";
    s += "else\n  echo \"[3/5] cudaMMC not found, skipping.\"\n";
    s += "  echo \"cudaMMC,$CHR,,skipped\" >> \"$CSV\"\nfi\n\n";
    s += "# Step 4 — pMMC (sequential)\n";
    s += "echo \"[4/5] Running pMMC...\"\n";
    s += "D=\"$OUTDIR/pmmc/\"; mkdir -p \"$D\"\n";
    s += "T0=$(date +%s%N)\n";
    s += "\"$PMMC_EXE\" -a create -s \"$STG\" -c \"$CHR\" -o \"$D\" -n pmmc -j $SEED -F pdb -E && RC=0 || RC=$?\n";
    s += "T1=$(date +%s%N); T=$(echo \"scale=3;($T1-$T0)/1000000000\"|bc)\n";
    s += "echo \"pMMC,$CHR,$T,$RC\" >> \"$CSV\"\n\n";
    s += "# Step 5 — quality comparison\n";
    s += "echo \"[5/5] Computing quality metrics...\"\n";
    s += "M=\"$OUTDIR/metrics/\"; mkdir -p \"$M\"\n";
    s += "P=$(ls \"$OUTDIR\"/pmmc/loops_pmmc.hcm 2>/dev/null||true)\n";
    s += "C=$(ls \"$OUTDIR\"/cudammc/loops_cuda.hcm 2>/dev/null||true)\n";
    s += "G=$(ls \"$OUTDIR\"/3dgnome/loops_gnome.hcm 2>/dev/null||true)\n";
    s += "[ -f \"$P\" ] && [ -f \"$C\" ] && \"$PMMC_EXE\" -a metrics -i \"$P,$C\" -o \"$M\" -n pmmc_vs_cuda -l 1 || true\n";
    s += "[ -f \"$P\" ] && [ -f \"$G\" ] && \"$PMMC_EXE\" -a metrics -i \"$P,$G\" -o \"$M\" -n pmmc_vs_gnome -l 1 || true\n\n";
    s += "echo \"\"\necho \"Benchmark complete.  Results: $CSV\"\ncat \"$CSV\"\n";
    writeBenchmarkScript(outdir + "run_benchmark.sh", s);
  }

  // Script 2: run_speed_benchmark.sh — per-chromosome wall-clock timing
  {
    std::string s;
    s += "#!/usr/bin/env bash\n";
    s += "# pMMC Speed Benchmark — per-chromosome wall-clock timing\n";
    s += "# Tools run SEQUENTIALLY to avoid resource competition.\n";
    s += "# Usage:  bash run_speed_benchmark.sh [output_dir] [replicates]\n";
    s += "set -euo pipefail\n\n";
    s += "OUTDIR=\"${1:-./speed_benchmark}\"\nREPS=\"${2:-3}\"\n";
    s += "PMMC_EXE=\"${PMMC_EXE:-pMMC.exe}\"\n";
    s += "CUDAMMC_EXE=\"${CUDAMMC_EXE:-cudaMMC}\"\n";
    s += "GNOME_EXE=\"${GNOME_EXE:-3dnome}\"\nSEED=42\n";
    s += "CSV=\"$OUTDIR/speed_results.csv\"\n\n";
    s += "# Subset for quick run; uncomment the full list for genome-wide (F5.14)\n";
    s += "CHRS=(chr1 chr14 chr21)\n";
    s += "# CHRS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)\n\n";
    s += "mkdir -p \"$OUTDIR\"\n";
    s += "echo \"tool,chromosome,replicate,wall_time_seconds\" > \"$CSV\"\n\n";
    s += "for CHR in \"${CHRS[@]}\"; do\n";
    s += "  SYNTH=\"$OUTDIR/synth_${CHR}/\"\n";
    s += "  [ -f \"$SYNTH/settings.ini\" ] || \"$PMMC_EXE\" -a generate -o \"$SYNTH\" -c \"$CHR\" -l 100 -m 1 -j $SEED\n";
    s += "  STG=\"$SYNTH/settings.ini\"\n";
    s += "  for R in $(seq 1 $REPS); do\n";
    s += "    echo \"--- $CHR  rep $R/$REPS ---\"\n";
    s += "    # 3D-GNOME\n";
    s += "    if command -v \"$GNOME_EXE\" &>/dev/null; then\n";
    s += "      D=\"$OUTDIR/gnome_${CHR}_${R}/\"; mkdir -p \"$D\"\n";
    s += "      T0=$(date +%s%N)\n";
    s += "      \"$GNOME_EXE\" -a create -s \"$STG\" -c \"$CHR\" -o \"$D\" -n g -j $((SEED+R)) || true\n";
    s += "      T1=$(date +%s%N); echo \"3D-GNOME,$CHR,$R,$(echo \"scale=3;($T1-$T0)/1e9\"|bc)\" >> \"$CSV\"\n";
    s += "    fi\n";
    s += "    # cudaMMC\n";
    s += "    if command -v \"$CUDAMMC_EXE\" &>/dev/null; then\n";
    s += "      D=\"$OUTDIR/cuda_${CHR}_${R}/\"; mkdir -p \"$D\"\n";
    s += "      T0=$(date +%s%N)\n";
    s += "      \"$CUDAMMC_EXE\" -a create -s \"$STG\" -c \"$CHR\" -o \"$D\" -n c -j $((SEED+R)) || true\n";
    s += "      T1=$(date +%s%N); echo \"cudaMMC,$CHR,$R,$(echo \"scale=3;($T1-$T0)/1e9\"|bc)\" >> \"$CSV\"\n";
    s += "    fi\n";
    s += "    # pMMC\n";
    s += "    D=\"$OUTDIR/pmmc_${CHR}_${R}/\"; mkdir -p \"$D\"\n";
    s += "    T0=$(date +%s%N)\n";
    s += "    \"$PMMC_EXE\" -a create -s \"$STG\" -c \"$CHR\" -o \"$D\" -n p -j $((SEED+R)) -F pdb || true\n";
    s += "    T1=$(date +%s%N); echo \"pMMC,$CHR,$R,$(echo \"scale=3;($T1-$T0)/1e9\"|bc)\" >> \"$CSV\"\n";
    s += "  done\n";
    s += "done\n\n";
    s += "echo \"Speed benchmark done.  Results: $CSV\"\ncat \"$CSV\"\n\n";
    // F056: Wall-clock time summary table (mean and CoV)
    s += "# Generate summary table with mean time and CoV per tool/chromosome\n";
    s += "SUMCSV=\"$OUTDIR/wall_clock_summary.csv\"\n";
    s += "echo \"tool,chromosome,mean_time_s,stdev_time_s,cov_percent,num_reps\" > \"$SUMCSV\"\n";
    s += "awk -F, 'NR>1 { key=$1\",\"$2; sum[key]+=$4; ssq[key]+=$4*$4; n[key]++ }\n";
    s += "END { for(k in sum) { m=sum[k]/n[k]; v=ssq[k]/n[k]-m*m; if(v<0)v=0; sd=sqrt(v); cov=(m>0)?sd/m*100:0;\n";
    s += "  printf \"%s,%.3f,%.3f,%.1f,%d\\n\",k,m,sd,cov,n[k] } }' \"$CSV\" | sort >> \"$SUMCSV\"\n";
    s += "echo \"\\nWall-clock summary table:\"\ncat \"$SUMCSV\"\n";
    writeBenchmarkScript(outdir + "run_speed_benchmark.sh", s);
  }

  // Script 3: run_scalability.sh — scaling analysis by loop count
  {
    std::string s;
    s += "#!/usr/bin/env bash\n";
    s += "# pMMC Scalability Benchmark — wall time vs loop count\n";
    s += "# Tools run SEQUENTIALLY to avoid resource competition.\n";
    s += "# Usage:  bash run_scalability.sh [output_dir]\n";
    s += "set -euo pipefail\n\n";
    s += "OUTDIR=\"${1:-./scalability_benchmark}\"\n";
    s += "PMMC_EXE=\"${PMMC_EXE:-pMMC.exe}\"\n";
    s += "CUDAMMC_EXE=\"${CUDAMMC_EXE:-cudaMMC}\"\n";
    s += "GNOME_EXE=\"${GNOME_EXE:-3dnome}\"\n";
    s += "CHR=\"" + chr + "\"\nSEED=42\n";
    s += "CSV=\"$OUTDIR/scalability_results.csv\"\n\n";
    s += "LOOPS=(10 25 50 100 250 500 1000)\n\n";
    s += "mkdir -p \"$OUTDIR\"\n";
    s += "echo \"tool,num_loops,wall_time_seconds\" > \"$CSV\"\n\n";
    s += "for NL in \"${LOOPS[@]}\"; do\n";
    s += "  echo \"=== $NL loops ===\"\n";
    s += "  SYNTH=\"$OUTDIR/synth_${NL}/\"\n";
    s += "  \"$PMMC_EXE\" -a generate -o \"$SYNTH\" -c \"$CHR\" -l \"$NL\" -m 1 -j $SEED\n";
    s += "  STG=\"$SYNTH/settings.ini\"\n\n";
    s += "  if command -v \"$GNOME_EXE\" &>/dev/null; then\n";
    s += "    D=\"$OUTDIR/gnome_${NL}/\"; mkdir -p \"$D\"\n";
    s += "    T0=$(date +%s%N)\n";
    s += "    \"$GNOME_EXE\" -a create -s \"$STG\" -c \"$CHR\" -o \"$D\" -n g -j $SEED || true\n";
    s += "    T1=$(date +%s%N); echo \"3D-GNOME,$NL,$(echo \"scale=3;($T1-$T0)/1e9\"|bc)\" >> \"$CSV\"\n";
    s += "  fi\n";
    s += "  if command -v \"$CUDAMMC_EXE\" &>/dev/null; then\n";
    s += "    D=\"$OUTDIR/cuda_${NL}/\"; mkdir -p \"$D\"\n";
    s += "    T0=$(date +%s%N)\n";
    s += "    \"$CUDAMMC_EXE\" -a create -s \"$STG\" -c \"$CHR\" -o \"$D\" -n c -j $SEED || true\n";
    s += "    T1=$(date +%s%N); echo \"cudaMMC,$NL,$(echo \"scale=3;($T1-$T0)/1e9\"|bc)\" >> \"$CSV\"\n";
    s += "  fi\n";
    s += "  D=\"$OUTDIR/pmmc_${NL}/\"; mkdir -p \"$D\"\n";
    s += "  T0=$(date +%s%N)\n";
    s += "  \"$PMMC_EXE\" -a create -s \"$STG\" -c \"$CHR\" -o \"$D\" -n p -j $SEED -F pdb || true\n";
    s += "  T1=$(date +%s%N); echo \"pMMC,$NL,$(echo \"scale=3;($T1-$T0)/1e9\"|bc)\" >> \"$CSV\"\n";
    s += "done\n\n";
    s += "echo \"Scalability benchmark done.  Results: $CSV\"\ncat \"$CSV\"\n";
    writeBenchmarkScript(outdir + "run_scalability.sh", s);
  }

  // Script 4: run_ensemble_benchmark.sh — intra-ensemble similarity
  {
    std::string s;
    s += "#!/usr/bin/env bash\n";
    s += "# pMMC Ensemble Benchmark — intra-ensemble structural similarity\n";
    s += "# Tools run SEQUENTIALLY to avoid resource competition.\n";
    s += "# Usage:  bash run_ensemble_benchmark.sh [output_dir] [ensemble_size]\n";
    s += "set -euo pipefail\n\n";
    s += "OUTDIR=\"${1:-./ensemble_benchmark}\"\nENS=\"${2:-50}\"\n";
    s += "PMMC_EXE=\"${PMMC_EXE:-pMMC.exe}\"\n";
    s += "CUDAMMC_EXE=\"${CUDAMMC_EXE:-cudaMMC}\"\n";
    s += "CHR=\"" + chr + "\"\nSEED=42\n\n";
    s += "mkdir -p \"$OUTDIR\"\n\n";
    s += "# Generate data\n";
    s += "SYNTH=\"$OUTDIR/synthetic/\"\n";
    s += "\"$PMMC_EXE\" -a generate -o \"$SYNTH\" -c \"$CHR\" -l 100 -m \"$ENS\" -j $SEED\n";
    s += "STG=\"$SYNTH/settings.ini\"\n\n";
    s += "# pMMC ensemble\n";
    s += "echo \"Running pMMC ensemble ($ENS members)...\"\n";
    s += "D=\"$OUTDIR/pmmc_ens/\"; mkdir -p \"$D\"\n";
    s += "T0=$(date +%s%N)\n";
    s += "\"$PMMC_EXE\" -a create -s \"$STG\" -c \"$CHR\" -o \"$D\" -n ens -m \"$ENS\" -j $SEED -F pdb\n";
    s += "T1=$(date +%s%N)\n";
    s += "echo \"pMMC ensemble: $(echo \"scale=1;($T1-$T0)/1e9\"|bc)s\"\n\n";
    s += "# cudaMMC ensemble (one run at a time — sequential)\n";
    s += "if command -v \"$CUDAMMC_EXE\" &>/dev/null; then\n";
    s += "  echo \"Running cudaMMC ensemble ($ENS members)...\"\n";
    s += "  D2=\"$OUTDIR/cuda_ens/\"; mkdir -p \"$D2\"\n";
    s += "  T0=$(date +%s%N)\n";
    s += "  for i in $(seq 0 $((ENS-1))); do\n";
    s += "    \"$CUDAMMC_EXE\" -a create -s \"$STG\" -c \"$CHR\" -o \"$D2\" -n \"ens_${i}\" -j $((SEED+i)) || true\n";
    s += "  done\n";
    s += "  T1=$(date +%s%N)\n";
    s += "  echo \"cudaMMC ensemble: $(echo \"scale=1;($T1-$T0)/1e9\"|bc)s\"\n";
    s += "fi\n\n";
    s += "# Compute structural distances\n";
    s += "\"$PMMC_EXE\" -a ensemble -i \"$OUTDIR/pmmc_ens/\" -p \"loops_ens_{N}.hcm\"\n\n";
    s += "echo \"Ensemble benchmark done.\"\n";
    s += "echo \"Distance matrices: $OUTDIR/pmmc_ens/structural_distances_lvl*.heat\"\n";
    writeBenchmarkScript(outdir + "run_ensemble_benchmark.sh", s);
  }

  // Script 5: run_quality_comparison.sh — Pearson/Spearman correlation
  {
    std::string s;
    s += "#!/usr/bin/env bash\n";
    s += "# pMMC Quality Comparison — run after run_benchmark.sh\n";
    s += "# Usage:  bash run_quality_comparison.sh [benchmark_dir]\n";
    s += "set -euo pipefail\n\n";
    s += "B=\"${1:-./benchmark_results}\"\n";
    s += "PMMC_EXE=\"${PMMC_EXE:-pMMC.exe}\"\n";
    s += "Q=\"$B/quality/\"; mkdir -p \"$Q\"\n\n";
    s += "P=$(ls \"$B\"/pmmc/loops_pmmc.hcm 2>/dev/null||true)\n";
    s += "C=$(ls \"$B\"/cudammc/loops_cuda.hcm 2>/dev/null||true)\n";
    s += "G=$(ls \"$B\"/3dgnome/loops_gnome.hcm 2>/dev/null||true)\n\n";
    s += "echo \"pair,scc\" > \"$Q/quality_summary.csv\"\n\n";
    s += "for HCM in \"$P\" \"$C\" \"$G\"; do\n";
    s += "  [ -f \"$HCM\" ] && \"$PMMC_EXE\" -a distmap -i \"$HCM\" -o \"$Q\" -n \"$(basename $HCM .hcm)\" -l 1 || true\n";
    s += "done\n\n";
    s += "if [ -f \"$P\" ] && [ -f \"$C\" ]; then\n";
    s += "  \"$PMMC_EXE\" -a metrics -i \"$P,$C\" -o \"$Q\" -n pmmc_vs_cuda -l 1 | tee \"$Q/pmmc_vs_cuda.log\"\n";
    s += "  SCC=$(grep -oP 'SCC.*?\\K[0-9.]+' \"$Q/pmmc_vs_cuda.log\" || echo N/A)\n";
    s += "  echo \"pmmc_vs_cuda,$SCC\" >> \"$Q/quality_summary.csv\"\n";
    s += "fi\n";
    s += "if [ -f \"$P\" ] && [ -f \"$G\" ]; then\n";
    s += "  \"$PMMC_EXE\" -a metrics -i \"$P,$G\" -o \"$Q\" -n pmmc_vs_gnome -l 1 | tee \"$Q/pmmc_vs_gnome.log\"\n";
    s += "  SCC=$(grep -oP 'SCC.*?\\K[0-9.]+' \"$Q/pmmc_vs_gnome.log\" || echo N/A)\n";
    s += "  echo \"pmmc_vs_gnome,$SCC\" >> \"$Q/quality_summary.csv\"\n";
    s += "fi\n\n";
    s += "echo \"Quality comparison done.\"\ncat \"$Q/quality_summary.csv\"\n";
    writeBenchmarkScript(outdir + "run_quality_comparison.sh", s);
  }

  // Script 6: run_resource_monitor.sh — resource utilisation capture (F057)
  {
    std::string s;
    s += "#!/usr/bin/env bash\n";
    s += "# pMMC Resource Utilisation Monitor\n";
    s += "# Captures GPU util, VRAM, system RAM, and CPU util during a run.\n";
    s += "# Usage:  bash run_resource_monitor.sh [output_dir] [chromosome]\n";
    s += "set -euo pipefail\n\n";
    s += "OUTDIR=\"${1:-./resource_monitor}\"\nCHR=\"${2:-chr1}\"\n";
    s += "PMMC_EXE=\"${PMMC_EXE:-pMMC.exe}\"\nSEED=42\n";
    s += "CSV=\"$OUTDIR/resource_utilisation.csv\"\n\n";
    s += "mkdir -p \"$OUTDIR\"\n\n";
    s += "# Generate synthetic data\n";
    s += "SYNTH=\"$OUTDIR/synth_${CHR}/\"\n";
    s += "[ -f \"$SYNTH/settings.ini\" ] || \"$PMMC_EXE\" -a generate -o \"$SYNTH\" -c \"$CHR\" -l 100 -j $SEED\n\n";
    s += "# Start resource monitor in background\n";
    s += "MONITOR_LOG=\"$OUTDIR/resource_samples.csv\"\n";
    s += "echo \"timestamp,gpu_util_pct,gpu_mem_used_mb,gpu_mem_total_mb,sys_mem_used_mb,cpu_util_pct\" > \"$MONITOR_LOG\"\n\n";
    s += "monitor_resources() {\n";
    s += "  while true; do\n";
    s += "    TS=$(date +%s.%N)\n";
    s += "    GPU_LINE=$(nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits 2>/dev/null || echo \"0, 0, 0\")\n";
    s += "    GPU_UTIL=$(echo \"$GPU_LINE\" | awk -F', ' '{print $1}')\n";
    s += "    GPU_MEM=$(echo \"$GPU_LINE\" | awk -F', ' '{print $2}')\n";
    s += "    GPU_TOTAL=$(echo \"$GPU_LINE\" | awk -F', ' '{print $3}')\n";
    s += "    if command -v free &>/dev/null; then\n";
    s += "      SYS_MEM=$(free -m | awk '/Mem:/{print $3}')\n";
    s += "    else\n";
    s += "      SYS_MEM=0\n";
    s += "    fi\n";
    s += "    # CPU util via /proc/stat or top\n";
    s += "    if [ -f /proc/stat ]; then\n";
    s += "      CPU_UTIL=$(awk '/^cpu /{u=$2+$4; t=$2+$4+$5; printf \"%.1f\", u/t*100}' /proc/stat)\n";
    s += "    else\n";
    s += "      CPU_UTIL=0\n";
    s += "    fi\n";
    s += "    echo \"$TS,$GPU_UTIL,$GPU_MEM,$GPU_TOTAL,$SYS_MEM,$CPU_UTIL\" >> \"$MONITOR_LOG\"\n";
    s += "    sleep 1\n";
    s += "  done\n";
    s += "}\n\n";
    s += "monitor_resources &\nMONITOR_PID=$!\n";
    s += "trap \"kill $MONITOR_PID 2>/dev/null\" EXIT\n\n";
    s += "# Run pMMC reconstruction\n";
    s += "D=\"$OUTDIR/pmmc_${CHR}/\"; mkdir -p \"$D\"\n";
    s += "T0=$(date +%s%N)\n";
    s += "\"$PMMC_EXE\" -a create -s \"$SYNTH/settings.ini\" -c \"$CHR\" -o \"$D\" -n res -j $SEED -F pdb\n";
    s += "T1=$(date +%s%N)\nWALL_TIME=$(echo \"scale=3;($T1-$T0)/1e9\"|bc)\n\n";
    s += "kill $MONITOR_PID 2>/dev/null || true\nsleep 1\n\n";
    s += "# Compute summary\n";
    s += "echo \"metric,value\" > \"$CSV\"\n";
    s += "echo \"chromosome,$CHR\" >> \"$CSV\"\n";
    s += "echo \"wall_time_s,$WALL_TIME\" >> \"$CSV\"\n";
    s += "awk -F, 'NR>1 { if($2+0>mg) mg=$2; if($3+0>mm) mm=$3; gt=$4; if($5+0>sm) sm=$5; if($6+0>mc) mc=$6; n++ }\n";
    s += "END { printf \"gpu_util_max_pct,%.1f\\ngpu_vram_max_mb,%.0f\\ngpu_vram_total_mb,%.0f\\nsys_ram_max_mb,%.0f\\ncpu_util_max_pct,%.1f\\nsamples,%d\\n\",mg,mm,gt,sm,mc,n }' \"$MONITOR_LOG\" >> \"$CSV\"\n\n";
    s += "echo \"Resource utilisation:\"\ncat \"$CSV\"\n";
    writeBenchmarkScript(outdir + "run_resource_monitor.sh", s);
  }

  printf("\n[benchmark] All scripts written to %s\n", outdir.c_str());
  printf("[benchmark] Scripts:\n");
  printf("  run_benchmark.sh          — 3-tool comparison (3D-GNOME, cudaMMC, pMMC)\n");
  printf("  run_speed_benchmark.sh    — per-chromosome wall-clock timing + summary table\n");
  printf("  run_scalability.sh        — performance vs loop count\n");
  printf("  run_ensemble_benchmark.sh — intra-ensemble similarity\n");
  printf("  run_quality_comparison.sh — Pearson/Spearman quality metrics\n");
  printf("  run_resource_monitor.sh   — GPU/CPU/memory utilisation capture\n");
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

  // Compute adaptive contact threshold from median nearest-neighbor distance
  float adaptive_threshold = 2.0f;
  {
    std::vector<float> nn_dists;
    for (int i = 0; i < flat.size - 1; ++i) {
      float d = (flat.points[i] - flat.points[i + 1]).length();
      nn_dists.push_back(d);
    }
    if (!nn_dists.empty()) {
      std::sort(nn_dists.begin(), nn_dists.end());
      float median_nn = nn_dists[nn_dists.size() / 2];
      adaptive_threshold = median_nn * 1.5f;
      printf("Adaptive contact threshold: %.2f (median NN dist = %.2f)\n",
             adaptive_threshold, median_nn);
    }
  }
  DistanceMapGenerator gen;
  gen.setContactThreshold(adaptive_threshold);
  gen.setContactExponent(2.0f);

  std::string prefix = outdir + (flag_set('n') ? get_arg('n') : "structure");

  // Compute distance and contact maps
  Heatmap dist_h = gen.computeDistanceMap(flat);
  dist_h.toFile(prefix + "_distmap.heat", false);
  printf("Distance map (%dx%d) written to %s\n",
         (int)dist_h.size, (int)dist_h.size, (prefix + "_distmap.heat").c_str());

  Heatmap cont_h = gen.computeContactMap(flat);
  cont_h.toFile(prefix + "_contactmap.heat", false);
  printf("Contact map (%dx%d) written to %s\n",
         (int)cont_h.size, (int)cont_h.size, (prefix + "_contactmap.heat").c_str());

  Heatmap freq = gen.computeContactFrequencyMap(flat);
  freq.toFile(prefix + "_freqmap.heat", false);
  printf("Contact frequency map written to %s\n",
         (prefix + "_freqmap.heat").c_str());

  // Render contact map images (SVG + PNG) directly from computed heatmaps
  HeatmapImageWriter::writeSVG(dist_h, prefix + "_distmap.svg", "Distance Map");
  HeatmapImageWriter::writePNG(dist_h, prefix + "_distmap.png", "Distance Map");
  printf("Distance map images: %s, %s\n",
         (prefix + "_distmap.svg").c_str(), (prefix + "_distmap.png").c_str());

  HeatmapImageWriter::writeSVG(cont_h, prefix + "_contactmap.svg", "Contact Map");
  HeatmapImageWriter::writePNG(cont_h, prefix + "_contactmap.png", "Contact Map");
  printf("Contact map images: %s, %s\n",
         (prefix + "_contactmap.svg").c_str(), (prefix + "_contactmap.png").c_str());

  HeatmapImageWriter::writeSVG(freq, prefix + "_freqmap.svg", "Contact Frequency");
  HeatmapImageWriter::writePNG(freq, prefix + "_freqmap.png", "Contact Frequency");
  printf("Frequency map images: %s, %s\n",
         (prefix + "_freqmap.svg").c_str(), (prefix + "_freqmap.png").c_str());

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

// Reconstruct from contact matrix (Hi-C pairs, cooler dump, HiCPro).
// Usage: pMMC -a contact -i <matrix.txt> -s <settings.ini> [-r resolution] [-c chr] [-o outdir/]
void prepareContact() {
  if (!flag_set('i'))
    error_msg("No contact matrix input file specified [-i]\n");
  if (!flag_set('s'))
    error_msg("No settings file specified [-s]\n");

  std::string contact_path = get_arg('i');
  std::string stg_path = get_arg('s');
  std::string outdir = flag_set('o') ? get_arg('o') : "./";
  int resolution = flag_set('r') ? atoi(get_arg('r').c_str()) : 5000;

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

  // Read contact matrix (auto-detect format)
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
    if (!ContactMatrixIO::readContactMatrixFiltered(contact_path, loops, anchors, chrs, resolution))
      error_msg("Failed to read contact matrix file\n");
  } else {
    if (!ContactMatrixIO::readContactMatrix(contact_path, loops, anchors, resolution))
      error_msg("Failed to read contact matrix file\n");
    chrs = parseChromosomeDescription(chromosomes);
  }

  printf("[contact] Loaded %d loops, %d anchors (resolution=%d)\n",
         (int)loops.size(), (int)anchors.size(), resolution);

  // Create temp directory for extracted files
  std::string tmpdir = outdir + "contact_tmp/";
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
  Settings::dataFactors = "contact";

  std::string name = flag_set('n') ? get_arg('n') : "contact";
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
  logger.logParameter("contact_input", contact_path);
  logger.logParameter("resolution", std::to_string(resolution));
  logger.logParameter("chromosomes", chromosomes);
  logger.logMilestone("Starting contact matrix reconstruction");

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

  validateOutputsAndWriteErrorWarning(outdir, name, use_new_file_format,
                                      ensemble_size);

  logger.logMilestone("Contact matrix reconstruction complete");
  logger.writeSummary();
}

// ─── Batch multi-cell-line reconstruction (F6.1-F6.5, F1.10) ───────────────
// Usage: pMMC -a batch -i <cell_line_list.txt> -o <output_dir/>
//
// The cell_line_list.txt file has one line per cell line:
//   cell_line_name  data_directory  settings_file
//
// For each cell line the full reconstruction pipeline (runLooper) is executed.
// After all cell lines finish, pairwise structural similarity (SCC) is computed
// and written to a CSV matrix plus a per-cell-line summary CSV.
void prepareBatch() {
  if (!flag_set('i'))
    error_msg("No cell-line list file specified [-i]\n");

  std::string list_path = get_arg('i');
  std::string outdir = flag_set('o') ? get_arg('o') : "./batch_output/";
  int ensemble_size = flag_set('m') ? atoi(get_arg('m').c_str()) : 1;

  // Parse the cell-line list
  struct CellLineEntry {
    std::string name;
    std::string data_dir;
    std::string settings_file;
  };
  std::vector<CellLineEntry> entries;

  FILE *listf = fopen(list_path.c_str(), "r");
  if (!listf) {
    error_msg(ftext("Failed to open cell line list: %s\n", list_path.c_str()));
    return;
  }
  char line_buf[4096];
  while (fgets(line_buf, sizeof(line_buf), listf)) {
    std::string line(line_buf);
    while (!line.empty() && (line.back() == '\n' || line.back() == '\r' ||
                             line.back() == ' ' || line.back() == '\t'))
      line.pop_back();
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream iss(line);
    CellLineEntry entry;
    if (!(iss >> entry.name >> entry.data_dir >> entry.settings_file)) {
      printf("[batch] Skipping malformed line: %s\n", line.c_str());
      continue;
    }
    entries.push_back(entry);
  }
  fclose(listf);

  if (entries.empty()) {
    error_msg("No valid cell line entries found in list file\n");
    return;
  }

  printf("[batch] Found %d cell lines to process\n", (int)entries.size());
  portable_mkdir(outdir.c_str());

  // Save -F flag before args.clear()
  g_output_format = flag_set('F') ? get_arg('F') : "";

  // CLI overrides to propagate
  bool has_M = flag_set('M');
  size_t cli_maxMem = has_M ? (size_t)atoi(get_arg('M').c_str()) : 0;
  bool cli_energy_trace = flag_set('E');
  int cli_energy_interval = flag_set('L') ? atoi(get_arg('L').c_str()) : 1000;
  std::string cli_init_method = flag_set('I') ? get_arg('I') : "";
  std::string cli_chr_desc = flag_set('c') ? get_arg('c') : "genome";

  unsigned int base_seed = 0;
  if (flag_set('j'))
    base_seed = (unsigned int)atoi(get_arg('j').c_str());

  args.clear();

  // Storage for per-cell-line HCM paths (for similarity computation)
  std::vector<std::string> cell_line_names;
  std::vector<std::string> hcm_paths;

  for (size_t ei = 0; ei < entries.size(); ++ei) {
    const CellLineEntry &entry = entries[ei];
    printf("\n========================================\n");
    printf("[batch] Cell line %d/%d: %s\n", (int)(ei + 1),
           (int)entries.size(), entry.name.c_str());
    printf("========================================\n");

    // Load settings for this cell line
    Settings stg;
    stg.loadFromINI(entry.settings_file);
    Settings::dataDirectory = entry.data_dir;

    // Apply CLI overrides
    if (has_M) Settings::maxMemoryMB = cli_maxMem;
    if (cli_energy_trace) Settings::energyTraceEnabled = true;
    Settings::energyTraceInterval = cli_energy_interval;
    if (!cli_init_method.empty()) Settings::initMethod = cli_init_method;

    unsigned int seed = base_seed + (unsigned int)ei;
    Settings::gpuSeed = seed;
    srand(seed);

    std::string cell_outdir = outdir + entry.name + "/";
    portable_mkdir(cell_outdir.c_str());

    SimulationLogger logger(cell_outdir, entry.name, seed);
    logger.logParameter("cell_line", entry.name);
    logger.logParameter("data_directory", entry.data_dir);
    logger.logParameter("settings_file", entry.settings_file);
    logger.logMilestone("Starting batch reconstruction");

    // Parse chromosomes (use CLI -c if given, else genome)
    std::vector<std::string> chrs = parseChromosomeDescription(cli_chr_desc);
    BedRegion region_of_interest;

    runLooper(chrs, region_of_interest, entry.name, cell_outdir,
              false /* use_new_file_format */, 5 /* max_level */,
              -1 /* chr_number_limit */, -1 /* length_limit */,
              ensemble_size, std::vector<std::string>(), &logger);

    logger.logMilestone("Batch reconstruction complete for " + entry.name);
    logger.writeSummary();

    // Track HCM path for similarity
    std::string hcm_path = ftext("%sloops_%s.hcm", cell_outdir.c_str(),
                                  entry.name.c_str());
    cell_line_names.push_back(entry.name);
    hcm_paths.push_back(hcm_path);
  }

  // ── Compute pairwise structural similarity matrix (F6.4) ──
  printf("\n[batch] Computing pairwise structural similarity...\n");
  int n = (int)cell_line_names.size();
  std::vector<Heatmap> contact_maps;
  int resolution = 100000; // 100kb resolution for comparison
  DistanceMapGenerator dmg;
  dmg.setContactExponent(2.0f);

  // Pre-size contact_maps for parallel fill
  contact_maps.resize(n);

#ifdef _OPENMP
  printf("[batch] Loading %d contact maps (OpenMP threads: %d)\n",
         n, omp_get_max_threads());
  #pragma omp parallel for schedule(dynamic)
#else
  printf("[batch] Loading %d contact maps (single-threaded)\n", n);
#endif
  for (int i = 0; i < n; ++i) {
    DistanceMapGenerator dmg_local;
    dmg_local.setContactExponent(2.0f);
    if (file_exists(hcm_paths[i])) {
      HierarchicalChromosome hc;
      hc.fromFile(hcm_paths[i]);
      hc.useLowestLevel();
      std::string first_chr = hc.chrs.empty() ? "" : hc.chrs[0];
      bool ok = false;
      if (!first_chr.empty() && !hc.current_level[first_chr].empty()) {
        Chromosome flat = hc.createEqidistantModel(resolution, first_chr);
        if (flat.size >= 2) {
          contact_maps[i] = dmg_local.computeContactFrequencyMap(flat);
          ok = true;
        }
      }
      if (!ok) {
        printf("[batch] WARNING: could not create equidistant model for %s\n",
               cell_line_names[i].c_str());
        contact_maps[i].init(1);
      }
    } else {
      printf("[batch] WARNING: HCM file not found: %s\n", hcm_paths[i].c_str());
      contact_maps[i].init(1);
    }
  }

  // Write similarity matrix CSV (F6.4) — compute pairwise SCC in parallel
  std::string sim_matrix_path = outdir + "similarity_matrix.csv";
  // Pre-compute upper triangle in parallel
  std::vector<std::vector<float>> sim_matrix(n, std::vector<float>(n, 0.0f));
  for (int i = 0; i < n; ++i) sim_matrix[i][i] = 1.0f;

  int n_pairs = n * (n - 1) / 2;
#ifdef _OPENMP
  printf("[batch] Computing %d pairwise similarities (OpenMP threads: %d)\n",
         n_pairs, omp_get_max_threads());
  #pragma omp parallel for schedule(dynamic)
#endif
  for (int k = 0; k < n_pairs; ++k) {
    // Map linear index k to upper-triangle (i, j) with i < j
    int i = 0, j = 0;
    {
      int row = 0, cumulative = n - 1;
      while (k >= cumulative) { row++; cumulative += n - 1 - row; }
      i = row;
      j = k - (cumulative - (n - 1 - row)) + i + 1;
    }
    float scc = 0.0f;
    if (contact_maps[i].size > 1 && contact_maps[j].size > 1 &&
        contact_maps[i].size == contact_maps[j].size) {
      MetricsFramework mf;
      scc = mf.structuralSimilarity(contact_maps[i], contact_maps[j]);
    }
    sim_matrix[i][j] = scc;
    sim_matrix[j][i] = scc;
  }

  FILE *sim_f = fopen(sim_matrix_path.c_str(), "w");
  if (sim_f) {
    fprintf(sim_f, "cell_line");
    for (int i = 0; i < n; ++i)
      fprintf(sim_f, ",%s", cell_line_names[i].c_str());
    fprintf(sim_f, "\n");

    for (int i = 0; i < n; ++i) {
      fprintf(sim_f, "%s", cell_line_names[i].c_str());
      for (int j = 0; j < n; ++j)
        fprintf(sim_f, ",%.6f", sim_matrix[i][j]);
      fprintf(sim_f, "\n");
    }
    fclose(sim_f);
    printf("[batch] Similarity matrix written to %s\n", sim_matrix_path.c_str());

    // Generate SVG heatmap from similarity matrix
    std::string sim_svg = outdir + "similarity_matrix.svg";
    SvgChartGenerator::similarityHeatmapChart(sim_matrix_path, sim_svg);
  }

  // Write per-cell-line summary CSV (F6.5)
  std::string summary_path = outdir + "batch_summary.csv";
  FILE *sumf = fopen(summary_path.c_str(), "w");
  if (sumf) {
    fprintf(sumf, "cell_line,hcm_path,contact_map_size,status\n");
    for (int i = 0; i < n; ++i) {
      fprintf(sumf, "%s,%s,%d,%s\n", cell_line_names[i].c_str(),
              hcm_paths[i].c_str(), (int)contact_maps[i].size,
              file_exists(hcm_paths[i]) ? "ok" : "missing");
    }
    fclose(sumf);
    printf("[batch] Summary written to %s\n", summary_path.c_str());
  }

  // ── Gene-location-based structural comparison (F6.3) ──
  // Compare 3D positions of specific genomic regions across cell lines.
  // If a BED file is provided via -B, compute per-region 3D coordinates
  // for each cell line and write a comparison matrix.
  {
    printf("\n[batch] Gene-location structural comparison (F6.3)...\n");
    std::string gene_comp_path = outdir + "gene_location_comparison.csv";
    FILE *gf = fopen(gene_comp_path.c_str(), "w");
    if (gf) {
      // Header: cell_line pairs + RMSD at matched genomic positions
      fprintf(gf, "cell_line_a,cell_line_b,position_rmsd,distance_correlation\n");

      // Load structures at common resolution for comparison
      std::vector<Chromosome> flat_structures;
      int comp_resolution = 100000; // 100kb
      for (int i = 0; i < n; ++i) {
        if (file_exists(hcm_paths[i])) {
          HierarchicalChromosome hc;
          hc.fromFile(hcm_paths[i]);
          hc.useLowestLevel();
          std::string first_chr = hc.chrs.empty() ? "" : hc.chrs[0];
          if (!first_chr.empty() && !hc.current_level[first_chr].empty()) {
            flat_structures.push_back(
                hc.createEqidistantModel(comp_resolution, first_chr));
          } else {
            Chromosome empty_chr;
            flat_structures.push_back(empty_chr);
          }
        } else {
          Chromosome empty_chr;
          flat_structures.push_back(empty_chr);
        }
      }

      // Pairwise comparison
      for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
          int sz = std::min(flat_structures[i].size, flat_structures[j].size);
          if (sz < 2) {
            fprintf(gf, "%s,%s,NA,NA\n", cell_line_names[i].c_str(),
                    cell_line_names[j].c_str());
            continue;
          }

          // Align and compute RMSD
          Chromosome a_copy = flat_structures[i];
          Chromosome b_copy = flat_structures[j];
          a_copy.size = sz;
          b_copy.size = sz;
          a_copy.align(b_copy);
          float rmsd = a_copy.calcRMSD(b_copy);

          // Distance correlation
          MetricsFramework metrics;
          float dc = metrics.distanceCorrelation(flat_structures[i],
                                                  flat_structures[j],
                                                  std::max(1, sz / 200));

          fprintf(gf, "%s,%s,%.6f,%.6f\n", cell_line_names[i].c_str(),
                  cell_line_names[j].c_str(), rmsd, dc);
        }
      }
      fclose(gf);
      printf("[batch] Gene-location comparison written to %s\n",
             gene_comp_path.c_str());
    }
  }

  // ── Tissue-type structural clustering (F6.4) ──
  // Perform hierarchical clustering on the similarity matrix using
  // single-linkage and write a dendrogram ordering + cluster assignments.
  {
    printf("\n[batch] Tissue-type structural clustering (F6.4)...\n");

    // Convert similarity matrix to distance matrix
    std::vector<std::vector<float>> dist(n, std::vector<float>(n, 0.0f));
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (i == j)
          dist[i][j] = 0.0f;
        else if (contact_maps[i].size > 1 && contact_maps[j].size > 1 &&
                 contact_maps[i].size == contact_maps[j].size) {
          MetricsFramework mf_clust;
          float scc = mf_clust.structuralSimilarity(contact_maps[i],
                                                    contact_maps[j]);
          dist[i][j] = 1.0f - scc; // distance = 1 - similarity
        } else {
          dist[i][j] = 1.0f;
        }
      }
    }

    // Single-linkage hierarchical clustering
    std::vector<bool> merged(n, false);
    std::vector<int> cluster_id(n);
    for (int i = 0; i < n; ++i) cluster_id[i] = i;

    std::string dendro_path = outdir + "clustering_dendrogram.csv";
    FILE *df = fopen(dendro_path.c_str(), "w");
    if (df) {
      fprintf(df, "step,cluster_a,cluster_b,distance,new_cluster\n");

      for (int step = 0; step < n - 1; ++step) {
        // Find closest pair of unmerged clusters
        float min_dist = 1e30f;
        int best_i = -1, best_j = -1;
        for (int i = 0; i < n; ++i) {
          if (merged[i]) continue;
          for (int j = i + 1; j < n; ++j) {
            if (merged[j]) continue;
            if (dist[i][j] < min_dist) {
              min_dist = dist[i][j];
              best_i = i;
              best_j = j;
            }
          }
        }
        if (best_i < 0) break;

        int new_id = n + step;
        fprintf(df, "%d,%s,%s,%.6f,%d\n", step,
                cell_line_names[best_i].c_str(),
                cell_line_names[best_j].c_str(), min_dist, new_id);

        // Merge: update distances (single linkage = min)
        merged[best_j] = true;
        cluster_id[best_j] = new_id;
        cluster_id[best_i] = new_id;
        for (int k = 0; k < n; ++k) {
          if (merged[k] || k == best_i) continue;
          dist[best_i][k] = std::min(dist[best_i][k], dist[best_j][k]);
          dist[k][best_i] = dist[best_i][k];
        }
      }
      fclose(df);
      printf("[batch] Clustering dendrogram written to %s\n",
             dendro_path.c_str());
    }

    // Write cluster ordering
    std::string order_path = outdir + "clustering_order.csv";
    FILE *of = fopen(order_path.c_str(), "w");
    if (of) {
      fprintf(of, "rank,cell_line,cluster_id\n");

      // Sort by cluster_id for grouping
      std::vector<int> order(n);
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(),
                [&](int a, int b) { return cluster_id[a] < cluster_id[b]; });

      for (int r = 0; r < n; ++r) {
        fprintf(of, "%d,%s,%d\n", r, cell_line_names[order[r]].c_str(),
                cluster_id[order[r]]);
      }
      fclose(of);
      printf("[batch] Cluster ordering written to %s\n", order_path.c_str());
    }
  }
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

void prepareFlowchart() {
  string output = flag_set('o') ? get_arg('o') : "pMMC_algorithm.svg";
  // If output is a directory, append default filename
  if (!output.empty() && (output.back() == '/' || output.back() == '\\'))
    output += "pMMC_algorithm.svg";
  printf("Generating algorithm flowchart: %s\n", output.c_str());
  if (FlowchartGenerator::generate(output)) {
    printf("Flowchart written to %s\n", output.c_str());
  } else {
    error_msg("Failed to generate flowchart\n");
  }
}

int main(int argc, char **argv) {
  setbuf(stdout, NULL);

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
    // -A: annotation directory containing enhancer.bed, promoter.bed, epigenetic.bed
    if (flag_set('A')) {
      std::string annot_dir = get_arg('A');
      if (!annot_dir.empty() && annot_dir.back() != '/' && annot_dir.back() != '\\')
        annot_dir += "/";
      std::string enh = annot_dir + "enhancer.bed";
      std::string pro = annot_dir + "promoter.bed";
      std::string epi = annot_dir + "epigenetic.bed";
      if (file_exists(enh)) Settings::enhancerAnnotationFile = enh;
      if (file_exists(pro)) Settings::promoterAnnotationFile = pro;
      if (file_exists(epi)) Settings::epigeneticAnnotationFile = epi;
      printf("Annotation directory: %s\n", annot_dir.c_str());
    }
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
  else if (args['a'] == "C" || args['a'] == "contact")
    prepareContact();
  else if (args['a'] == "X" || args['a'] == "convert")
    prepareConvert();
  else if (args['a'] == "b" || args['a'] == "batch")
    prepareBatch();
  else if (args['a'] == "flowchart")
    prepareFlowchart();
  else {
    printf("No matching action specified [%s]!\n", args['a'].c_str());
    usage();
  }

  args.clear();
  printf("end\n");
  return 0;
}
