// main_convert.cpp — "convert" app for pMMC
// Commands: convert, extract, flatten, smooth, merge

#include <AppCommon.hpp>
#include <Static.hpp>
#include <BedRegion.hpp>
#include <CifWriter.hpp>
#include <HierarchicalChromosome.h>
#include <Heatmap.hpp>
#include <MergedFileIO.hpp>
#include <PdbWriter.hpp>
#include <Settings.hpp>
#include <VisualizationScripts.hpp>

using namespace std;

// ── Internal helpers ──────────────────────────────────────────────

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

void extract(string input_path, string output_path, int start, int end,
             const string &format = "") {
  printf("extract fragment: %d %d\n", start, end);
  HierarchicalChromosome h;
  h.fromFile(input_path);
  printf("size = %d\n", (int)h.clusters.size());
  if (h.clusters.size() > 0) {
    HierarchicalChromosome extr = h.extractFragment(start, end);
    extr.toFile(output_path);

    if (!format.empty()) {
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

// ── Command implementations (settings-driven) ────────────────────

void prepareSmooth() {
  string input = Settings::inputFile;
  int res = Settings::resolution;
  int level = Settings::level;
  string chr = Settings::chromosomes;
  string output = Settings::outputDir.empty() ? input + ".smooth.txt" : Settings::outputDir;

  if (input.empty()) error_msg("No input_file specified in settings\n");
  if (res <= 0) error_msg("No resolution specified in settings\n");

  HierarchicalChromosome hc;
  hc.fromFile(input);
  hc.setLevel(level);

  Chromosome ch = hc.createEqidistantModel(res, chr);
  ch.toFile(output);
}

void prepareFlatten() {
  string input = Settings::inputFile;
  int level = Settings::level;
  string chr = Settings::chromosomes;
  string output = Settings::outputDir.empty() ? input + ".txt" : Settings::outputDir;

  if (input.empty()) error_msg("No input_file specified in settings\n");

  flatten(input, level, chr, output);
}

void prepareExtract() {
  string input = Settings::inputFile;
  string output = Settings::outputDir;
  string format = Settings::outputFormat;

  if (input.empty()) error_msg("No input_file specified in settings\n");
  if (output.empty()) error_msg("No output_dir specified in settings\n");

  int start = 0, end = 0;

  // Try chromosomes field as a BED region (chr:start:end)
  string region_str = Settings::chromosomes;
  if (!region_str.empty() && BedRegion::tryParse(region_str.c_str())) {
    BedRegion region;
    region.parse(region_str);
    start = region.start;
    end = region.end;
  } else {
    start = Settings::extractStart;
    end = Settings::extractEnd;
    if (start == 0 && end == 0)
      error_msg("No extract region specified. Use chromosomes=chr:start:end or extract_start/extract_end\n");
  }

  extract(input, output, start, end, format);

  // Generate visualization scripts if PDB was produced
  if (format == "pdb" || format == "both") {
    string base = output;
    if (base.size() > 4 && base.substr(base.size() - 4) == ".hcm")
      base = base.substr(0, base.size() - 4);
    string pdb_path = base + ".pdb";

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

void prepareMerge() {
  std::string bed_path = Settings::bedFilePath;
  std::string bedpe_path = Settings::bedpeFilePath;
  std::string output_path = Settings::outputDir;

  if (bed_path.empty()) error_msg("No bed_file specified in settings [data] section\n");
  if (bedpe_path.empty()) error_msg("No bedpe_file specified in settings [data] section\n");
  if (output_path.empty()) error_msg("No output_dir specified in settings\n");

  printf("Merging BED [%s] + BEDPE [%s] -> [%s]\n",
         bed_path.c_str(), bedpe_path.c_str(), output_path.c_str());

  if (!MergedFileIO::writeMerged(bed_path, bedpe_path, output_path))
    error_msg("Failed to create merged file\n");
}

void prepareConvert() {
  string input = Settings::inputFile;
  string outdir = Settings::outputDir.empty() ? "./" : Settings::outputDir;
  string format = Settings::outputFormat.empty() ? "both" : Settings::outputFormat;

  if (input.empty()) error_msg("No input_file specified in settings\n");

  printf("[convert] Reading HCM: %s\n", input.c_str());
  HierarchicalChromosome hc;
  hc.fromFile(input);

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

// ── Usage ──────────────────────────────────────────────────────────

void usage() {
  printf("pMMC-convert — file conversion and manipulation\n\n");
  printf("Usage: pMMC-convert -s <settings.ini>\n\n");
  printf("All parameters are controlled by the settings file.\n\n");
  printf("Key settings.ini fields in [main]:\n");
  printf("  action          Command: convert, extract, flatten, smooth, merge\n");
  printf("  output_dir      Output file or directory\n");
  printf("  output_format   Format: pdb, cif, both (for convert/extract)\n");
  printf("  chromosomes     Chromosome or region chr:start:end (for extract/flatten/smooth)\n");
  printf("  level           Hierarchy level (for flatten/smooth)\n");
  printf("  resolution      Resolution in bp (for smooth)\n");
  printf("  extract_start   Start position (for extract, alternative to chromosomes)\n");
  printf("  extract_end     End position (for extract, alternative to chromosomes)\n");
  printf("\nKey settings.ini fields in [data]:\n");
  printf("  input_file      Input HCM file\n");
  printf("  bed_file        BED file path (for merge)\n");
  printf("  bedpe_file      BEDPE file path (for merge)\n");
}

// ── Main ───────────────────────────────────────────────────────────

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

  string action = Settings::action;
  if (action.empty()) {
    printf("No action specified in settings\n");
    usage();
    return 1;
  }

  if (action == "convert")
    prepareConvert();
  else if (action == "extract")
    prepareExtract();
  else if (action == "flatten")
    prepareFlatten();
  else if (action == "smooth")
    prepareSmooth();
  else if (action == "merge")
    prepareMerge();
  else {
    printf("Unknown action: %s\n\n", action.c_str());
    usage();
    return 1;
  }

  printf("end\n");
  return 0;
}
