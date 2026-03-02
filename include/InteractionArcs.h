/**
 * @file InteractionArcs.h
 * @brief Manages ChIA-PET anchors and interaction arcs for all chromosomes.
 *
 * Loads anchor positions and PET cluster data from files, maps genomic
 * coordinates to anchor indices (via markArcs/parallelMarkArcs), filters
 * long arcs, and provides per-chromosome lookup tables.
 *
 * @see InteractionArc for a single arc.
 * @see Anchor for a single anchor.
 */

#ifndef INTERACTIONARCS_H_
#define INTERACTIONARCS_H_

#include <algorithm>
#include <iostream>
#include <list>
#include <map>
#include <numeric>
#include <set>
#include <stdio.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <common.h>
#include <Anchor.h>
#include <BedRegion.h>
#include <BedRegions.h>
#include <InteractionArc.h>
#include <Settings.h>

/**
 * @class InteractionArcs
 * @brief Container for all ChIA-PET anchors and interaction arcs across chromosomes.
 *
 * Workflow:
 * 1. loadAnchorsData() — read anchor BED file
 * 2. loadPetClustersData() — read PET cluster files (one per factor)
 * 3. markArcs() or parallelMarkArcs() — map genomic coords to anchor indices
 * 4. removeEmptyAnchors() — drop anchors with no arcs
 *
 * Data is organized per-chromosome in unordered_maps keyed by chromosome name.
 */
class InteractionArcs {
public:
  /** @brief Default constructor. */
  InteractionArcs();

  /**
   * @brief Read arcs and anchors from a binary cache file.
   * @param filename Input file path.
   * @return True on success.
   */
  bool fromFile(string filename);

  /**
   * @brief Write arcs and anchors to a binary cache file.
   * @param filename Output file path.
   */
  void toFile(string filename);

  /** @brief Clear all stored data. */
  void clear();

  /**
   * @brief Restrict reconstruction to a specific genomic region.
   *
   * Only anchors and arcs within the region will be loaded.
   * @param region BedRegion to restrict to.
   */
  void selectRegion(BedRegion region);

  /**
   * @brief Load anchor positions from a tab-delimited file.
   *
   * File format: chr start end orientation (one anchor per line).
   * @param anchors_path Path to the anchors file.
   */
  void loadAnchorsData(string anchors_path);

  /**
   * @brief Load PET cluster data and store as raw arcs.
   *
   * Reads interaction data, filters by segment boundaries, and stores
   * the results in raw_arcs. Long arcs (> maxPETClusterLength) are
   * stored separately in long_arcs.
   *
   * @param pet_clusters_path Path to the PET clusters file.
   * @param factor_name Name of the ChIA-PET factor (e.g. "CTCF").
   * @param predefined_segments Segment boundaries for cross-segment filtering.
   */
  void loadPetClustersData(string pet_clusters_path, string factor_name,
                           BedRegions predefined_segments);

  /**
   * @brief Print a summary of anchors and arcs per chromosome.
   * @param display_limit Maximum number of arcs to display per chromosome.
   */
  void print(int display_limit = 10);

  /**
   * @brief Map raw arcs (genomic positions) to anchor indices.
   *
   * Input: anchors and raw_arcs. Output: arcs (with anchor-index start/end).
   * @param ignore_missing If true, silently skip arcs with unmapped anchors.
   */
  void markArcs(bool ignore_missing);

  /**
   * @brief GPU-accelerated version of markArcs() using CUDA.
   * @param ignore_missing If true, silently skip arcs with unmapped anchors.
   * @see ParallelMarkArcs.cu
   */
  void parallelMarkArcs(bool ignore_missing);

  /**
   * @brief Remove anchors that have no associated arcs.
   *
   * Should be called after markArcs(), and then markArcs() should be
   * called again to update arc indices.
   */
  void removeEmptyAnchors();

  /**
   * @brief Rewire arcs based on CTCF motif orientation rules.
   *
   * Implements the convergent CTCF loop extrusion model.
   */
  void rewire();

  std::vector<string> chrs;     /**< List of chromosome names with data. */
  std::vector<string> factors;  /**< List of factor names (e.g. "CTCF", "RNAPII"). */

  std::unordered_map<string, int> anchors_cnt;  /**< Number of anchors per chromosome. */
  std::unordered_map<string, int> arcs_cnt;      /**< Number of arcs per chromosome. */

  std::unordered_map<string, std::vector<Anchor>> anchors;  /**< Anchors per chromosome. */

  /** @brief Arcs with anchor-index endpoints (after markArcs). */
  std::unordered_map<string, std::vector<InteractionArc>> arcs;

  /** @brief Raw arcs with genomic-position endpoints (before markArcs). */
  std::unordered_map<string, std::vector<InteractionArc>> raw_arcs;

  /** @brief Long-range arcs (> maxPETClusterLength), used to refine segment-level heatmap. */
  std::unordered_map<string, std::vector<InteractionArc>> long_arcs;

  std::vector<float> expected_distances;  /**< Expected spatial distances per arc (from frequency-to-distance mapping). */

  BedRegion selected_region;  /**< Active region-of-interest filter (empty = whole genome). */
};

#endif /* INTERACTIONARCS_H_ */
