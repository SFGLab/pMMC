/**
 * @file PicturePanel.hpp
 * @brief Packages region visualization output for a genomic BED region.
 *
 * Extracts a fragment from a HierarchicalChromosome, writes the PDB file,
 * and returns metadata describing the exported region.
 */

#ifndef PICTUREPANEL_H_
#define PICTUREPANEL_H_

#include <string>

class HierarchicalChromosome;
class BedRegion;

/**
 * @class PicturePanel
 * @brief Metadata and PDB export for a single genomic region visualization.
 *
 * The generate() factory method extracts the sub-structure for a BED region,
 * writes a PDB file, and returns a PicturePanel populated with all metadata.
 */
class PicturePanel {
public:
  std::string pdbFilePath;  /**< Path to the exported PDB file. */
  std::string regionLabel;  /**< Human-readable region label (e.g. "chr1:1000000-2000000"). */
  int beadCount;            /**< Number of beads (atoms) in the extracted fragment. */
  std::string chromosome;   /**< Chromosome name (e.g. "chr1"). */
  int startBp;              /**< Start position in base pairs. */
  int endBp;                /**< End position in base pairs. */

  /**
   * @brief Extract a region from a HierarchicalChromosome, write PDB, and return metadata.
   * @param hc The hierarchical chromosome structure to extract from.
   * @param region The BED region specifying the genomic interval.
   * @param outputDir Directory where the PDB file will be written.
   * @return PicturePanel with all metadata populated.
   */
  static PicturePanel generate(HierarchicalChromosome &hc,
                               const BedRegion &region,
                               const std::string &outputDir);

  /**
   * @brief Return a formatted summary string describing this panel.
   * @return Multi-line summary with region, bead count, and file path.
   */
  std::string toString() const;
};


// ============================================================================
// Implementation
// ============================================================================

/**
 * @file PicturePanel.cpp
 * @brief PicturePanel implementation: region extraction and PDB export.
 */

#include <BedRegion.hpp>
#include <HierarchicalChromosome.h>
#include <PdbWriter.hpp>
#include <common.h>
#include <cstdio>

inline PicturePanel PicturePanel::generate(HierarchicalChromosome &hc,
                                    const BedRegion &region,
                                    const std::string &outputDir) {
  PicturePanel panel;
  panel.chromosome = region.chr;
  panel.startBp = region.start;
  panel.endBp = region.end;
  panel.regionLabel = ftext("%s:%d-%d", region.chr.c_str(),
                            region.start, region.end);

  // Build output path: outputDir/<chr>_<start>_<end>.pdb
  std::string filename = ftext("%s%s_%d_%d.pdb", outputDir.c_str(),
                               region.chr.c_str(), region.start, region.end);
  panel.pdbFilePath = filename;

  // Extract the sub-fragment for this genomic range
  HierarchicalChromosome fragment = hc.extractFragment(region.start, region.end);

  // Count beads at the lowest level
  fragment.useLowestLevel();
  int count = 0;
  for (const auto &entry : fragment.current_level) {
    count += (int)entry.second.size();
  }
  panel.beadCount = count;

  // Write PDB file for the extracted fragment
  PdbWriter::write(fragment, filename);

  printf("[PicturePanel] Generated %s (%d beads)\n",
         panel.regionLabel.c_str(), panel.beadCount);

  return panel;
}

inline std::string PicturePanel::toString() const {
  return ftext("PicturePanel: %s\n  Chromosome: %s\n  Range: %d - %d bp\n"
               "  Beads: %d\n  PDB: %s",
               regionLabel.c_str(), chromosome.c_str(),
               startBp, endBp, beadCount, pdbFilePath.c_str());
}

#endif /* PICTUREPANEL_H_ */
