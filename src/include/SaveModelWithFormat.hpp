#pragma once
// SaveModelWithFormat.hpp — Shared output format helper for apps that
// write HCM/PDB/CIF files.

#include <string>
#include <HierarchicalChromosome.h>
#include <CifWriter.hpp>
#include <PdbWriter.hpp>
#include <common.h>

// Requires: extern std::string g_output_format;  (from AppCommon.hpp)
extern std::string g_output_format;

inline void saveModelWithFormat(HierarchicalChromosome &hc, const std::string &outdir,
                                const std::string &label, bool use_new_file_format,
                                const std::string &suffix = "",
                                bool useCurrentLevel = false) {
  // Always save HCM
  std::string hcm_path;
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
    std::string base = ftext("%s%s%s", outdir.c_str(), label.c_str(), suffix.c_str());
    if (g_output_format == "cif" || g_output_format == "mmcif") {
      CifWriter::write(hc, base + ".cif", "cudaMMC", useCurrentLevel);
    } else if (g_output_format == "pdb") {
      PdbWriter::write(hc, base + ".pdb", false);
    } else if (g_output_format == "both") {
      CifWriter::write(hc, base + ".cif", "cudaMMC", useCurrentLevel);
      PdbWriter::write(hc, base + ".pdb", false);
    } else {
      printf("[WARN] Unknown output format '%s'. Use: cif, pdb, or both\n", g_output_format.c_str());
    }
  }
}
