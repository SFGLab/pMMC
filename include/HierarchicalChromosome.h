/**
 * @file HierarchicalChromosome.h
 * @brief Multi-level hierarchical chromosome 3D structure.
 *
 * Stores the full 4-level tree (Chromosome > Segment > Anchor > Subanchor)
 * produced by LooperSolver. Supports navigation between levels, I/O in HCM
 * format, spline smoothing, and extraction of flat equidistant models.
 *
 * @see LooperSolver for the reconstruction engine that builds this structure.
 * @see Cluster for individual tree nodes.
 */

#ifndef HIERARCHICALCHROMOSOME_H_
#define HIERARCHICALCHROMOSOME_H_

#include <map>
#include <stdio.h>
#include <vector>

#include <Chromosome.hpp>
#include <Cluster.hpp>
#include <InteractionArcs.hpp>
#include <optional>
#ifdef __CUDACC__
#include <vector_types.h>
#else
struct float3 {};
#endif

/**
 * @class HierarchicalChromosome
 * @brief Container for the complete multi-level 3D chromosome reconstruction.
 *
 * The structure is a tree of Cluster nodes. Each chromosome has a root
 * cluster whose children are segments, which contain anchors, which
 * contain subanchors. Navigation is done by setting the active level
 * and reading current_level indices.
 */
class HierarchicalChromosome {
public:
  /** @brief Default constructor. */
  HierarchicalChromosome();

  /**
   * @brief Print the structure summary at a given verbosity level.
   * @param level Verbosity (1 = summary, higher = more detail).
   */
  void print(int level = 1);

  /** @brief Print the full regions tree for all chromosomes. */
  void printRegionsTree();

  /**
   * @brief Print the regions tree for a specific chromosome.
   * @param chr Chromosome name (e.g. "chr22").
   */
  void printRegionsTree(string chr);

  /**
   * @brief Write the structure to an HCM binary file.
   * @param filename Output file path.
   */
  void toFile(string filename);

  /**
   * @brief Write to an already-open file handle (HCM format).
   * @param file Open FILE pointer for writing.
   */
  void toFile(FILE *file);

  /**
   * @brief Write in the legacy HCM format.
   * @param filename Output file path.
   */
  void toFilePreviousFormat(string filename);

  /**
   * @brief Read from an HCM binary file.
   * @param filename Input file path.
   */
  void fromFile(string filename);

  /**
   * @brief Read from the legacy HCM format.
   * @param file Open FILE pointer for reading.
   */
  void fromFilePreviousFormat(FILE *file);

  /**
   * @brief Read from an already-open file handle (HCM format).
   * @param file Open FILE pointer for reading.
   */
  void fromFile(FILE *file);

  /**
   * @brief Read from a text-based structure file.
   * @param filename Input file path.
   * @return True on success.
   */
  bool fromTxt(string filename);

  /**
   * @brief Import a structure from HiC-Evo hierarchical chromosome format.
   * @param filename Input file path.
   * @return True on success.
   */
  bool fromHiCEvo(string filename);

  /** @brief Set active level to the topmost (chromosome) level. */
  void useTopLevel();

  /** @brief Set active level to the lowest (subanchor) level. */
  void useLowestLevel();

  /**
   * @brief Set the active hierarchy level.
   * @param level Level index: 0=chromosome, 1=segment, 2=anchor, 3=subanchor.
   */
  void setLevel(int level);

  /** @brief Move the active level one step down in the hierarchy. */
  void levelDown();

  /**
   * @brief Expand a genomic region at the current level.
   * @param start Start genomic position (bp).
   * @param end End genomic position (bp).
   * @param include_external If true, include clusters partially overlapping the region.
   */
  void expandRegion(int start, int end, bool include_external = true);

  /**
   * @brief Compute the center of mass of the structure.
   * @param current If true, use current-level clusters; if false, use all.
   * @return 3D center of mass vector.
   */
  vector3 getCenter(bool current = true);

  /**
   * @brief Translate the structure so the center of mass is at the origin.
   * @param current If true, center current-level clusters only.
   */
  void center(bool current = true);

  /**
   * @brief Get average spatial distance between anchors per chromosome.
   * @return Map from chromosome name to average anchor distance.
   */
  map<string, float> getAvgAnchorDistance();

  /**
   * @brief Find the 3D position for a given genomic coordinate.
   * @param chr Chromosome name.
   * @param pos Genomic position (bp).
   * @return Interpolated 3D position.
   */
  vector3 find3DPosition(string chr, int pos);

  /**
   * @brief Find the 3D position using the spline-smoothed structure.
   * @param chr Chromosome name.
   * @param pos Genomic position (bp).
   * @return Smoothed 3D position.
   */
  vector3 find3DSmoothPosition(string chr, int pos);

  /**
   * @brief Scale all 3D coordinates by a factor.
   * @param factor Scale multiplier.
   */
  void scale(float factor);

  /** @brief Build Chromosome objects from the current-level cluster positions. */
  void createCurrentLevelStructure();

  /**
   * @brief Create a pairwise-distance heatmap for a chromosome at a given level.
   * @param chr Chromosome name.
   * @param level Hierarchy level to extract distances from.
   * @return Heatmap where h[i][j] = 3D distance between clusters i and j.
   */
  Heatmap createStructuralHeatmap(string chr, int level);

  /**
   * @brief Compute structural distance to another HierarchicalChromosome.
   * @param hc Other structure.
   * @param level Level at which to compare.
   * @return Distance metric (sum of pairwise position differences).
   */
  float calcDistance(HierarchicalChromosome hc, int level);

  /**
   * @brief Apply spline smoothing with n interpolation points per segment.
   * @param n Number of interpolation points between consecutive clusters.
   */
  void smoothSpline(int n);

  /**
   * @brief Map a cluster index to its position in the smoothed array.
   * @param ind Cluster index.
   * @return Index in the smoothed points array.
   */
  int clusterIndexToSmoothedIndex(int ind);

  /**
   * @brief Create a flat equidistant bead model at fixed resolution.
   *
   * Produces a Chromosome with beads spaced exactly resolution_bp apart,
   * interpolated from the hierarchical structure.
   *
   * @param resolution_bp Base pairs per output bead.
   * @param chr Chromosome name (empty = use first chromosome).
   * @return Flat Chromosome with evenly spaced beads.
   */
  Chromosome createEqidistantModel(int resolution_bp, string chr = "");

  /**
   * @brief Extract a sub-fragment of the structure by genomic range.
   * @param start Start genomic position (bp).
   * @param end End genomic position (bp).
   * @return New HierarchicalChromosome containing only the fragment.
   */
  HierarchicalChromosome extractFragment(int start, int end);

  /**
   * @brief Get a heatmap of pairwise 3D distances between current-level beads.
   * @return Distance heatmap.
   */
  Heatmap getDistancesHeatmap();

  /** @brief Print spatial distribution statistics (diameter, distances) per chromosome. */
  void getSpatialDistributionStats();

  /**
   * @brief Find anchor indices flanking a genomic region.
   * @param chr Chromosome name.
   * @param genomic_position_start Region start (bp).
   * @param genomic_position_end Region end (bp).
   * @return Vector of flanking anchor indices.
   */
  vector<int> findFlankingAnchors(string chr, int genomic_position_start,
                                  int genomic_position_end);

  std::map<std::string, Chromosome> chr;         /**< Per-chromosome flat structure at current level. */
  std::map<std::string, Chromosome> chr_smooth;   /**< Per-chromosome spline-smoothed structure. */
  std::vector<string> chrs;                       /**< List of chromosome names in this structure. */
  std::vector<Cluster> clusters;                  /**< All tree nodes (beads at all levels). */
  std::map<std::string, std::vector<int>> current_level;  /**< Cluster indices at the active level, per chromosome. */
  std::map<std::string, int> chr_root;            /**< Root cluster index per chromosome. */
  InteractionArcs arcs;                           /**< Associated interaction arcs data. */

  /**
   * @brief GPU-accelerated Monte Carlo helper for heatmap reconstruction.
   * @param param1 Number of iterations.
   * @param param2 Block size.
   * @param param3 Thread count.
   * @param chr Chromosome name.
   * @return Optional vector of float3 positions (empty on failure).
   */
  std::optional<std::vector<float3>> gpuHelper(const int, const int, const int,
                                               std::string);

private:
  /**
   * @brief Recursively print the tree structure for a sub-tree.
   * @param region_index_curr Root cluster index.
   * @param margin Indentation depth.
   */
  void printRegionsTreeRecursive(int region_index_curr, int margin = 0);

  /**
   * @brief Recursively expand a genomic region within a sub-tree.
   * @param region_ind Sub-tree root.
   * @param start Start genomic position.
   * @param end End genomic position.
   * @param include_external Include partially overlapping clusters.
   * @return Vector of cluster indices in the expanded region.
   */
  std::vector<int> expandRegion(int region_ind, int start, int end,
                                bool include_external);

  int smooth_factor;  /**< Spline interpolation factor (set by smoothSpline). */
};

#endif /* HIERARCHICALCHROMOSOME_H_ */
