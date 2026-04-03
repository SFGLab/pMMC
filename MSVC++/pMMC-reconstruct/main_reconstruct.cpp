// main_reconstruct.cpp — pMMC "reconstruct" app
// Commands: create (default), simulate, ibed, contact, batch

#include <AppCommon.hpp>
#include <SaveModelWithFormat.hpp>
#include <BedRegion.hpp>
#include <BedRegions.hpp>
#include <ContactMatrixIO.hpp>
#include <CifWriter.hpp>
#include <DistanceMapGenerator.hpp>
#include <GenomeReconstructor.hpp>
#include <HierarchicalChromosome.h>
#include <HeatmapImageWriter.hpp>
#include <IbedFileIO.hpp>
#include <LooperSolver.hpp>
#include <MergedFileIO.hpp>
#include <MetricsFramework.hpp>
#include <MultiscaleEnergy.hpp>
#include <PdbWriter.hpp>
#include <Settings.hpp>
#include <SimulationLogger.hpp>
#include <SvgChartGenerator.hpp>
#include <VisualizationScripts.hpp>
#include <Static.hpp>

// Parsed input data for non-chiapet formats (held in memory, no temp files)
struct ParsedInputData {
  bool active = false;
  std::string format;
  std::vector<MergedLoop> loops;
  std::vector<MergedAnchor> anchors;
};
static ParsedInputData g_parsed_input;

// ---------------------------------------------------------------------------
// Usage
// ---------------------------------------------------------------------------
void usage() {
  printf("pMMC-reconstruct — 3D genome reconstruction\n\n");
  printf("Usage: pMMC-reconstruct -s <settings.ini>\n\n");
  printf("All parameters are controlled by the settings file.\n");
  printf("See config.ini for a complete reference of all fields.\n\n");
  printf("Key settings.ini sections:\n");
  printf("  [main]  chromosomes, output_dir, label, seed, ensemble_size,\n");
  printf("          output_format, max_level, annotation_dir\n");
  printf("  [data]  input_format (chiapet|merged|ibed|contact),\n");
  printf("          input_file, cell_line_list (enables batch mode)\n");
  exit(0);
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

  if (g_parsed_input.active) {
    // Feed pre-parsed data directly — no temp files
    lsm.setContactDataFromMemory(chrs, region_of_interest,
                                  g_parsed_input.anchors,
                                  g_parsed_input.loops,
                                  g_parsed_input.format);
    // Free parsed data now that it's been consumed
    g_parsed_input = ParsedInputData();
  } else {
    lsm.setContactData(chrs, region_of_interest, anchors, factors, arcs_clusters,
                       arcs_singletons, arcs_singletons_inter);
  }

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
          // Skip rendering if all values are zero (no data)
          float hmin, hmax;
          h.getRange(hmin, hmax);
          if (hmax > 1e-12f) {
            std::string base = heat_path.substr(0, heat_path.size() - 5);
            std::string title = suffix + label;
            HeatmapImageWriter::writePNG(h, base + ".png", title);
            HeatmapImageWriter::writeSVG(h, base + ".svg", title);
          }
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

// Forward declaration
bool loadInputData(const std::vector<std::string> &chrs);

void prepareLooper(SimulationLogger *logger = nullptr) {
  printf("prepare\n");

  // Parse chromosomes from settings
  string chromosomes = Settings::chromosomes;
  BedRegion region_of_interest;
  std::vector<string> chrs;

  if (BedRegion::tryParse(chromosomes)) {
    region_of_interest.parse(chromosomes);
    chrs.push_back(region_of_interest.chr);
  } else {
    chrs = parseChromosomeDescription(chromosomes);
  }

  string name = Settings::label;
  if (name.empty()) {
    name = chromosomes;
    std::replace(name.begin(), name.end(), ':', '_');
    std::replace(name.begin(), name.end(), '-', '_');
  }

  if (!Settings::templateSegment.empty() && !Settings::distHeatmap.empty())
    error("template_segment and dist_heatmap are exclusive, you may only provide one of them");

  std::vector<string> selected_factors = vector<string>();
  if (!Settings::selectedFactors.empty())
    selected_factors = split(Settings::selectedFactors);

  string outdir = Settings::outputDir;
  int length_limit = Settings::lengthLimit;
  int chr_number_limit = Settings::chrNumberLimit;
  int max_level = Settings::maxLevel;
  int ensemble_size = Settings::ensembleSize;
  bool use_new_file_format = Settings::useNewFileFormat;

  g_output_format = Settings::outputFormat;

  // Parse non-chiapet input if configured in settings
  loadInputData(chrs);

  runLooper(chrs, region_of_interest, name, outdir, use_new_file_format,
            max_level, chr_number_limit, length_limit, ensemble_size,
            selected_factors, logger);
}

// ---------------------------------------------------------------------------
// loadInputData — Parse input file into memory based on Settings::inputFormat.
// No temp files are written. Data is stored in g_parsed_input for runLooper.
// ---------------------------------------------------------------------------
bool loadInputData(const std::vector<std::string> &chrs) {
  std::string fmt = Settings::inputFormat;
  if (fmt.empty() || fmt == "chiapet")
    return false;  // use existing ChIA-PET file-based pipeline

  std::string input_path = Settings::inputFile;
  if (input_path.empty())
    error_msg("input_format is set but input_file is empty in settings\n");

  g_parsed_input = ParsedInputData();
  g_parsed_input.format = fmt;

  if (fmt == "merged") {
    if (!MergedFileIO::readMerged(input_path, g_parsed_input.loops, g_parsed_input.anchors))
      error_msg("Failed to read merged file\n");
    if (!chrs.empty())
      MergedFileIO::filterByChromosomes(g_parsed_input.loops, g_parsed_input.anchors, chrs);
  } else if (fmt == "ibed") {
    if (!chrs.empty()) {
      if (!IbedFileIO::readIbedFiltered(input_path, g_parsed_input.loops, g_parsed_input.anchors, chrs))
        error_msg("Failed to read ibed file\n");
    } else {
      if (!IbedFileIO::readIbed(input_path, g_parsed_input.loops, g_parsed_input.anchors))
        error_msg("Failed to read ibed file\n");
    }
  } else if (fmt == "contact") {
    int res = Settings::contactResolution;
    if (!chrs.empty()) {
      if (!ContactMatrixIO::readContactMatrixFiltered(input_path, g_parsed_input.loops, g_parsed_input.anchors, chrs, res))
        error_msg("Failed to read contact matrix\n");
    } else {
      if (!ContactMatrixIO::readContactMatrix(input_path, g_parsed_input.loops, g_parsed_input.anchors, res))
        error_msg("Failed to read contact matrix\n");
    }
  } else {
    error_msg(ftext("Unknown input_format: %s\n", fmt.c_str()));
  }

  g_parsed_input.active = true;
  printf("[input] Loaded %d loops, %d anchors (format=%s) — held in memory\n",
         (int)g_parsed_input.loops.size(), (int)g_parsed_input.anchors.size(), fmt.c_str());
  return true;
}

// ── DELETED: prepareSimulate, prepareIbed, prepareContact ──
// Input format is now driven by settings.ini [data] input_format field.
// The loadInputData() function above handles all non-chiapet formats.

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

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char **argv) {
  setbuf(stdout, NULL);

  // Only flag accepted: -s <settings.ini>
  if (argc < 3 || std::string(argv[1]) != "-s") {
    usage();
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
  setGpuSeed(seed);
  Settings::gpuSeed = seed;

  // Apply annotation directory
  if (!Settings::annotationDir.empty()) {
    std::string annot_dir = Settings::annotationDir;
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

  // Dispatch: batch mode if cell_line_list is set, otherwise single reconstruction
  if (!Settings::cellLineList.empty()) {
    // Inject cell_line_list path into args for prepareBatch (it reads -i flag)
    args['i'] = Settings::cellLineList;
    args['o'] = Settings::outputDir;
    if (Settings::ensembleSize > 1)
      args['m'] = std::to_string(Settings::ensembleSize);
    if (!Settings::outputFormat.empty())
      args['F'] = Settings::outputFormat;
    if (!Settings::chromosomes.empty() && Settings::chromosomes != "genome")
      args['c'] = Settings::chromosomes;
    if (Settings::seed != 0)
      args['j'] = std::to_string(Settings::seed);
    g_output_format = Settings::outputFormat;

    prepareBatch();
  } else {
    // Single reconstruction
    string outdir = Settings::outputDir;
    string label = Settings::label;
    portable_mkdir(outdir.c_str());

    SimulationLogger logger(outdir, label, seed);
    logger.logParameter("settings_file", stg_path);
    logger.logParameter("input_format", Settings::inputFormat.empty() ? "chiapet" : Settings::inputFormat);
    logger.logParameter("chromosomes", Settings::chromosomes);
    logger.logMilestone("Starting reconstruction");

    prepareLooper(&logger);

    // Run territory analysis on the reconstructed model
    {
      string hcm_path = ftext("%sloops_%s.hcm", outdir.c_str(), label.c_str());
      if (file_exists(hcm_path)) {
        HierarchicalChromosome hc_ter;
        hc_ter.fromFile(hcm_path);
        GenomeReconstructor gr;
        gr.analyze(hc_ter, 1);
        gr.writeReport(outdir, label);
      }
    }

    logMemoryUsage("after reconstruction");
    logger.logMilestone("Reconstruction complete");
    logger.writeSummary();
  }

  printf("end\n");
  return 0;
}
