// main_analyze.cpp — "analyze" app for pMMC
// Commands: distance, distmap, ensemble, metrics, rewiring, position

#include <AppCommon.hpp>
#include <Static.hpp>

#include <BedRegion.hpp>
#include <BedRegions.hpp>
#include <Chromosome.hpp>
#include <DistanceMapGenerator.hpp>
#include <HeatmapImageWriter.hpp>
#include <HierarchicalChromosome.h>
#include <Heatmap.hpp>
#include <MetricsFramework.hpp>
#include <SvgChartGenerator.hpp>

// ── Extracted functions ────────────────────────────────────────────────

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
  string input = Settings::inputFile;
  string output = Settings::outputDir;
  string regions_path = Settings::bedRegionsFile;
  bool add_header = Settings::addHeader;

  if (input.empty()) error_msg("No input_file specified in settings\n");
  if (output.empty()) error_msg("No output_dir specified in settings\n");
  if (regions_path.empty()) error_msg("No bed_regions_file specified in settings\n");

  // read regions to extract the position for
  BedRegions regions;
  regions.fromFile(regions_path);
  printf("%d regions loaded\n", (int)regions.regions.size());

  get3DPositions(input, regions, output, add_header);
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
  string input = Settings::inputFile;
  string chr = Settings::chromosomes;
  int level = Settings::level;
  string output = Settings::outputDir.empty() ? ftext("%s.dist.heat", input.c_str()) : Settings::outputDir;

  if (input.empty()) error_msg("No input_file specified in settings\n");
  if (chr.empty()) error_msg("No chromosomes specified in settings\n");

  HierarchicalChromosome hc;
  hc.fromFile(input);
  runCalcDistanceMatrix(hc, chr, level, output);
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
  string input = Settings::inputFile;
  string pattern = Settings::pattern;
  string output = Settings::outputDir.empty() ? input : Settings::outputDir;
  string name = Settings::label.empty() ? "ensemble" : Settings::label;

  if (input.empty()) error_msg("No input_file specified in settings\n");
  if (pattern.empty()) error_msg("No pattern specified in settings\n");

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
  string input = Settings::inputFile;
  string pattern = Settings::pattern;
  int level = Settings::level;
  int resolution = Settings::resolution;
  string output = Settings::outputDir.empty() ? input : Settings::outputDir;

  if (input.empty()) error_msg("No input_file specified in settings\n");
  if (pattern.empty()) error_msg("No pattern specified in settings\n");

  rewiringAnalysis(input, pattern, level, resolution, output);
}

void prepareDistMap() {
  std::string input = Settings::inputFile;
  std::string outdir = Settings::outputDir.empty() ? "./" : Settings::outputDir;
  std::string chr = Settings::chromosomes;
  int level = Settings::level;
  int resolution = Settings::resolution;

  if (input.empty()) error_msg("No input_file specified in settings\n");

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
      adaptive_threshold = median_nn * 1.5f;  // 1.5x median neighbor distance
      printf("Adaptive contact threshold: %.2f (median NN dist = %.2f)\n",
             adaptive_threshold, median_nn);
    }
  }
  DistanceMapGenerator gen;
  gen.setContactThreshold(adaptive_threshold);
  gen.setContactExponent(2.0f);

  std::string prefix = outdir + (Settings::label.empty() ? "structure" : Settings::label);

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
  std::string input_a = Settings::inputFile;
  std::string input_b = Settings::inputFile2;

  if (input_a.empty()) error_msg("No input_file specified in settings\n");
  if (input_b.empty()) error_msg("No input_file_2 specified in settings\n");

  std::vector<std::string> inputs = {input_a, input_b};

  std::string outdir = Settings::outputDir.empty() ? "./" : Settings::outputDir;
  std::string chr = Settings::chromosomes;
  int level = Settings::level;
  int resolution = Settings::resolution;

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

  std::string prefix = outdir + (Settings::label.empty() ? "comparison" : Settings::label);

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

// ── Usage ──────────────────────────────────────────────────────────────

void usage() {
  printf("pMMC-analyze — analysis commands\n\n");
  printf("Usage: pMMC-analyze -s <settings.ini>\n\n");
  printf("All parameters are controlled by the settings file.\n\n");
  printf("Key settings.ini fields in [main]:\n");
  printf("  action            Command: distance, distmap, ensemble, metrics, rewiring, position\n");
  printf("  input_file        Input HCM file or directory\n");
  printf("  input_file_2      Second input file (for metrics)\n");
  printf("  output_dir        Output directory\n");
  printf("  label             Name prefix for output files\n");
  printf("  chromosomes       Chromosome name\n");
  printf("  level             Hierarchy level (0=chr, 1=segment, 2=subanchor, 3=bead)\n");
  printf("  resolution        Resolution in bp\n");
  printf("  pattern           File pattern with {N} (for ensemble/rewiring)\n");
  printf("  bed_regions_file  BED regions file (for position)\n");
  printf("  add_header        Add header line (for position)\n");
  printf("  seed              Random seed (0=time-based)\n");
}

// ── Main ───────────────────────────────────────────────────────────────

int main(int argc, char **argv) {
  setbuf(stdout, NULL);

  if (argc < 3 || std::string(argv[1]) != "-s") {
    usage();
    return 0;
  }

  std::string stg_path = argv[2];
  Settings stg;
  printf("Load settings from [%s]\n", stg_path.c_str());
  if (!stg.loadFromINI(stg_path)) {
    printf("Failed to load settings file!\n");
    return 1;
  }

  // Seed setup
  unsigned int seed;
  if (Settings::seed != 0) {
    seed = (unsigned int)Settings::seed;
    printf("Using fixed random seed: %u\n", seed);
  } else {
    auto tmp = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    seed = (unsigned int)tmp;
    printf("Using time-based random seed: %u\n", seed);
  }
  srand(seed);

  string action = Settings::action;
  if (action.empty()) {
    printf("No action specified in settings\n");
    usage();
    return 1;
  }

  if (action == "distance")
    prepareCalcDistanceMatrix();
  else if (action == "distmap")
    prepareDistMap();
  else if (action == "ensemble")
    prepareEnsembleAnalysis();
  else if (action == "metrics")
    prepareMetrics();
  else if (action == "rewiring")
    prepareRewiring();
  else if (action == "position")
    preparePosition();
  else {
    printf("Unknown action: %s\n\n", action.c_str());
    usage();
    return 1;
  }

  printf("end\n");
  return 0;
}
