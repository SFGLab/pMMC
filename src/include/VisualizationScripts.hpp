/**
 * @file VisualizationScripts.hpp
 * @brief Generates PyMOL (.pml) and ChimeraX (.cxc) visualization scripts
 *        from HierarchicalChromosome models.
 *
 * These scripts load the corresponding PDB file and apply coloring,
 * representation, and annotation commands to visualize genomic features
 * such as CTCF anchors, TAD boundaries, and rainbow coloring by genomic
 * position.
 */

#ifndef VISUALIZATIONSCRIPTS_H_
#define VISUALIZATIONSCRIPTS_H_

#include <string>

class HierarchicalChromosome;

/**
 * @class VisualizationScripts
 * @brief Static methods to generate molecular visualization scripts.
 */
class VisualizationScripts {
public:
  /**
   * @brief Generate a PyMOL .pml script for the given structure and PDB.
   *
   * The script:
   *   - Loads the PDB file (F4.1)
   *   - Sets ribbon/tube representation (F4.2)
   *   - Applies rainbow coloring by residue sequence (F4.3, genomic position)
   *   - Colors CTCF anchors in red, other beads by chain color (F4.4)
   *   - Sets TAD/segment boundaries as dashed markers (F4.6)
   *   - Mouse rotation/zoom is standard in PyMOL (F4.7)
   *   - Includes commands to export PNG (F4.8)
   *
   * @param hc The hierarchical chromosome model.
   * @param pdb_path Path to the PDB file to load.
   * @param output_script_path Path for the output .pml script.
   * @return True on success.
   */
  static bool writePyMOLScript(HierarchicalChromosome &hc,
                                const std::string &pdb_path,
                                const std::string &output_script_path);

  /**
   * @brief Generate a ChimeraX .cxc script for the given structure and PDB.
   *
   * The script:
   *   - Opens the structure
   *   - Sets ribbon display
   *   - Colors by rainbow (bfactor contains genomic position)
   *   - Highlights CTCF anchors in red
   *
   * @param hc The hierarchical chromosome model.
   * @param pdb_path Path to the PDB file to open.
   * @param output_script_path Path for the output .cxc script.
   * @return True on success.
   */
  static bool writeChimeraXScript(HierarchicalChromosome &hc,
                                   const std::string &pdb_path,
                                   const std::string &output_script_path);
};


// ============================================================================
// Implementation
// ============================================================================

/**
 * @file VisualizationScripts.cpp
 * @brief Implementation of PyMOL and ChimeraX script generators.
 */

#include <HierarchicalChromosome.h>
#include <Cluster.hpp>
#include <cstdio>
#include <string>
#include <vector>
#include <set>

// ─── PyMOL Script (.pml) ────────────────────────────────────────────────────

inline bool VisualizationScripts::writePyMOLScript(
    HierarchicalChromosome &hc, const std::string &pdb_path,
    const std::string &output_script_path) {

  FILE *f = fopen(output_script_path.c_str(), "w");
  if (!f) {
    printf("[VisualizationScripts] Failed to open %s for writing\n",
           output_script_path.c_str());
    return false;
  }

  // Ensure lowest level is active for full detail
  hc.useLowestLevel();

  // ── Collect CTCF anchor residue indices and TAD boundary residue indices ──
  // CTCF anchors: clusters at level 2 (LVL_INTERACTION_BLOCK) that have arcs
  // We traverse current_level (lowest) and check which beads have arcs
  // (indicating they are CTCF anchors).
  std::set<int> ctcf_residues;       // 1-based residue sequence numbers
  std::set<int> tad_boundary_resids; // 1-based residue sequence numbers at segment boundaries

  int globalResSeq = 0;
  for (size_t ci = 0; ci < hc.chrs.size(); ++ci) {
    const std::string &chr = hc.chrs[ci];
    if (hc.current_level.find(chr) == hc.current_level.end())
      continue;
    const std::vector<int> &level = hc.current_level[chr];

    // Find segment boundaries: first and last bead index of each segment
    // Walk up from each bead to find its segment parent
    std::set<int> segment_starts_genomic;
    for (size_t j = 0; j < level.size(); ++j) {
      int idx = level[j];
      if (idx < 0 || idx >= (int)hc.clusters.size())
        continue;
      const Cluster &cl = hc.clusters[idx];
      // Walk up to find segment (level 1) parent
      int parent_idx = cl.parent;
      while (parent_idx >= 0 && parent_idx < (int)hc.clusters.size() &&
             hc.clusters[parent_idx].level > 1) {
        parent_idx = hc.clusters[parent_idx].parent;
      }
      if (parent_idx >= 0 && parent_idx < (int)hc.clusters.size() &&
          hc.clusters[parent_idx].level == 1) {
        segment_starts_genomic.insert(hc.clusters[parent_idx].start);
      }
    }

    for (size_t j = 0; j < level.size(); ++j) {
      int idx = level[j];
      if (idx < 0 || idx >= (int)hc.clusters.size())
        continue;
      globalResSeq++;
      const Cluster &cl = hc.clusters[idx];

      // Check if this bead is a CTCF anchor (has interaction arcs)
      if (!cl.arcs.empty()) {
        ctcf_residues.insert(globalResSeq);
      }

      // Check if this bead is at a segment boundary
      if (segment_starts_genomic.count(cl.start) > 0 ||
          segment_starts_genomic.count(cl.genomic_pos) > 0) {
        tad_boundary_resids.insert(globalResSeq);
      }
    }
  }

  // ── Write PyMOL script ──

  fprintf(f, "# PyMOL visualization script generated by pMMC\n");
  fprintf(f, "# Structure: %s\n\n", pdb_path.c_str());

  // F4.1: Load the PDB file
  fprintf(f, "# Load structure\n");
  fprintf(f, "load %s, genome\n\n", pdb_path.c_str());

  // F4.2: Set ribbon/tube representation
  fprintf(f, "# Set ribbon/tube representation\n");
  fprintf(f, "hide everything, genome\n");
  fprintf(f, "show ribbon, genome\n");
  fprintf(f, "set ribbon_width, 3.0\n");
  fprintf(f, "set ribbon_sampling, 10\n");
  fprintf(f, "set ribbon_smooth, 5\n\n");

  // F4.3: Rainbow coloring by residue sequence (maps to genomic position)
  fprintf(f, "# Rainbow coloring by residue sequence (genomic position)\n");
  fprintf(f, "spectrum count, rainbow, genome\n\n");

  // F4.4: Color CTCF anchors in red
  if (!ctcf_residues.empty()) {
    fprintf(f, "# Highlight CTCF anchors in red\n");
    fprintf(f, "select ctcf_anchors, ");
    bool first = true;
    // Build selection using residue ranges for efficiency
    std::vector<int> sorted_ctcf(ctcf_residues.begin(), ctcf_residues.end());
    for (size_t i = 0; i < sorted_ctcf.size(); ++i) {
      if (!first)
        fprintf(f, " or ");
      // Find contiguous range
      int range_start = sorted_ctcf[i];
      while (i + 1 < sorted_ctcf.size() &&
             sorted_ctcf[i + 1] == sorted_ctcf[i] + 1) {
        ++i;
      }
      int range_end = sorted_ctcf[i];
      if (range_start == range_end)
        fprintf(f, "resi %d", range_start);
      else
        fprintf(f, "resi %d-%d", range_start, range_end);
      first = false;
    }
    fprintf(f, "\n");
    fprintf(f, "color red, ctcf_anchors\n");
    fprintf(f, "show spheres, ctcf_anchors\n");
    fprintf(f, "set sphere_scale, 0.8, ctcf_anchors\n\n");
  }

  // F4.5: Color enhancers, promoters, and epigenetic marks by B-factor
  fprintf(f, "# Highlight enhancers (green) - B-factor 20.0\n");
  fprintf(f, "select enhancers, b > 15 and b < 25\n");
  fprintf(f, "color green, enhancers\n");
  fprintf(f, "show spheres, enhancers\n");
  fprintf(f, "set sphere_scale, 0.8, enhancers\n\n");

  fprintf(f, "# Highlight promoters (blue) - B-factor 30.0\n");
  fprintf(f, "select promoters, b > 25 and b < 35\n");
  fprintf(f, "color blue, promoters\n");
  fprintf(f, "show spheres, promoters\n");
  fprintf(f, "set sphere_scale, 0.8, promoters\n\n");

  fprintf(f, "# Highlight epigenetic marks (magenta) - B-factor 40.0\n");
  fprintf(f, "select epigenetic, b > 35 and b < 45\n");
  fprintf(f, "color magenta, epigenetic\n");
  fprintf(f, "show spheres, epigenetic\n");
  fprintf(f, "set sphere_scale, 0.8, epigenetic\n\n");

  // F4.6: TAD/segment boundaries as dashed pseudoatom markers
  if (!tad_boundary_resids.empty()) {
    fprintf(f, "# TAD/segment boundary markers\n");
    int marker_idx = 0;
    for (int resid : tad_boundary_resids) {
      fprintf(f, "select tad_bnd_%d, resi %d\n", marker_idx, resid);
      fprintf(f, "color yellow, tad_bnd_%d\n", marker_idx);
      fprintf(f, "show spheres, tad_bnd_%d\n", marker_idx);
      fprintf(f, "set sphere_scale, 1.2, tad_bnd_%d\n", marker_idx);
      marker_idx++;
    }
    fprintf(f, "\n");
  }

  // F4.7: Mouse rotation/zoom is standard in PyMOL (no special commands needed)
  fprintf(f, "# Viewing setup\n");
  fprintf(f, "zoom genome\n");
  fprintf(f, "orient\n");
  fprintf(f, "set depth_cue, 1\n");
  fprintf(f, "set ray_trace_fog, 0\n\n");

  // F4.8: Export PNG commands
  fprintf(f, "# Export to PNG (uncomment to use)\n");
  fprintf(f, "# ray 2400, 2400\n");
  fprintf(f, "# png %s.png, dpi=300\n\n",
          pdb_path.substr(0, pdb_path.size() - 4).c_str());

  fprintf(f, "# Deselect all\n");
  fprintf(f, "deselect\n");

  fclose(f);
  printf("[VisualizationScripts] Wrote PyMOL script: %s\n",
         output_script_path.c_str());
  return true;
}

// ─── ChimeraX Script (.cxc) ────────────────────────────────────────────────

inline bool VisualizationScripts::writeChimeraXScript(
    HierarchicalChromosome &hc, const std::string &pdb_path,
    const std::string &output_script_path) {

  FILE *f = fopen(output_script_path.c_str(), "w");
  if (!f) {
    printf("[VisualizationScripts] Failed to open %s for writing\n",
           output_script_path.c_str());
    return false;
  }

  // Ensure lowest level is active for full detail
  hc.useLowestLevel();

  // ── Collect CTCF anchor residue indices ──
  std::set<int> ctcf_residues;
  int globalResSeq = 0;
  for (size_t ci = 0; ci < hc.chrs.size(); ++ci) {
    const std::string &chr = hc.chrs[ci];
    if (hc.current_level.find(chr) == hc.current_level.end())
      continue;
    const std::vector<int> &level = hc.current_level[chr];
    for (size_t j = 0; j < level.size(); ++j) {
      int idx = level[j];
      if (idx < 0 || idx >= (int)hc.clusters.size())
        continue;
      globalResSeq++;
      const Cluster &cl = hc.clusters[idx];
      if (!cl.arcs.empty()) {
        ctcf_residues.insert(globalResSeq);
      }
    }
  }

  // ── Write ChimeraX script ──

  fprintf(f, "# ChimeraX visualization script generated by pMMC\n");
  fprintf(f, "# Structure: %s\n\n", pdb_path.c_str());

  // Open the structure
  fprintf(f, "# Open structure\n");
  fprintf(f, "open %s\n\n", pdb_path.c_str());

  // Set ribbon display
  fprintf(f, "# Set ribbon display\n");
  fprintf(f, "hide atoms\n");
  fprintf(f, "show ribbons\n");
  fprintf(f, "style ribbon width 3\n\n");

  // Color by rainbow using B-factor (genomic position)
  fprintf(f, "# Rainbow coloring by B-factor (genomic position in Mb)\n");
  fprintf(f, "color bfactor palette rainbow\n\n");

  // Highlight CTCF anchors in red
  if (!ctcf_residues.empty()) {
    fprintf(f, "# Highlight CTCF anchors in red\n");
    // Build residue selection string
    std::vector<int> sorted_ctcf(ctcf_residues.begin(), ctcf_residues.end());
    fprintf(f, "select ");
    bool first = true;
    for (size_t i = 0; i < sorted_ctcf.size(); ++i) {
      if (!first)
        fprintf(f, ",");
      int range_start = sorted_ctcf[i];
      while (i + 1 < sorted_ctcf.size() &&
             sorted_ctcf[i + 1] == sorted_ctcf[i] + 1) {
        ++i;
      }
      int range_end = sorted_ctcf[i];
      if (range_start == range_end)
        fprintf(f, ":%d", range_start);
      else
        fprintf(f, ":%d-%d", range_start, range_end);
      first = false;
    }
    fprintf(f, "\n");
    fprintf(f, "color sel red\n");
    fprintf(f, "show sel atoms\n");
    fprintf(f, "size sel atomRadius 1.5\n");
    fprintf(f, "~select\n\n");
  }

  // Highlight enhancers, promoters, and epigenetic marks by B-factor
  fprintf(f, "# Highlight enhancers (green) - B-factor 20.0\n");
  fprintf(f, "select ::bfactor>15 & ::bfactor<25\n");
  fprintf(f, "color sel green\n");
  fprintf(f, "show sel atoms\n");
  fprintf(f, "size sel atomRadius 1.5\n");
  fprintf(f, "~select\n\n");

  fprintf(f, "# Highlight promoters (blue) - B-factor 30.0\n");
  fprintf(f, "select ::bfactor>25 & ::bfactor<35\n");
  fprintf(f, "color sel blue\n");
  fprintf(f, "show sel atoms\n");
  fprintf(f, "size sel atomRadius 1.5\n");
  fprintf(f, "~select\n\n");

  fprintf(f, "# Highlight epigenetic marks (magenta) - B-factor 40.0\n");
  fprintf(f, "select ::bfactor>35 & ::bfactor<45\n");
  fprintf(f, "color sel magenta\n");
  fprintf(f, "show sel atoms\n");
  fprintf(f, "size sel atomRadius 1.5\n");
  fprintf(f, "~select\n\n");

  // Viewing setup
  fprintf(f, "# Viewing setup\n");
  fprintf(f, "view\n");
  fprintf(f, "lighting soft\n\n");

  // Export command (commented)
  fprintf(f, "# Export to PNG (uncomment to use)\n");
  fprintf(f, "# save %s.png width 2400 height 2400 supersample 4\n",
          pdb_path.substr(0, pdb_path.size() - 4).c_str());

  fclose(f);
  printf("[VisualizationScripts] Wrote ChimeraX script: %s\n",
         output_script_path.c_str());
  return true;
}

#endif /* VISUALIZATIONSCRIPTS_H_ */
