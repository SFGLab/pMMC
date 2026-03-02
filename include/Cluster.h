/**
 * @file Cluster.h
 * @brief A node (bead) in the hierarchical chromosome tree.
 *
 * The reconstruction uses a 4-level tree:
 *   Level 0 (Chromosome) > Level 1 (Segment ~2Mb) >
 *   Level 2 (Interaction Block / Anchor) > Level 3 (Subanchor ~10kb).
 *
 * Each Cluster holds its 3D position, parent/children links, and
 * the set of interaction arcs that touch it.
 */

#ifndef CLUSTER_H_
#define CLUSTER_H_

#include <stdio.h>
#include <vector>

#include <common.h>
#include <InteractionArc.h>

/**
 * @class Cluster
 * @brief A single bead / node in the hierarchical chromosome structure.
 *
 * Clusters form a tree: each cluster has a parent index and a list of
 * child indices. At the anchor level (LVL_INTERACTION_BLOCK), clusters
 * correspond to ChIA-PET anchors. At lower levels they subdivide
 * the chromatin between anchors.
 */
class Cluster {
public:
  /** @brief Default constructor (calls init()). */
  Cluster();

  /**
   * @brief Construct a cluster spanning [start, end] in genomic coordinates.
   * @param start Genomic start position (bp).
   * @param end Genomic end position (bp).
   */
  Cluster(int start, int end);

  /** @brief Reset all fields to default values. */
  void init();

  /** @brief Print cluster info to stdout. */
  void print();

  /**
   * @brief Write cluster data to a binary file (HCM format).
   * @param file Open file handle for writing.
   */
  void toFile(FILE *file);

  /**
   * @brief Write cluster in the legacy HCM format.
   * @param file Open file handle for writing.
   */
  void toFilePreviousFormat(FILE *file);

  /**
   * @brief Read cluster data from a binary file (HCM format).
   * @param file Open file handle for reading.
   */
  void fromFile(FILE *file);

  /**
   * @brief Check if a genomic position is contained in this cluster.
   * @param genomic_pos Position in base pairs.
   * @return True if start <= genomic_pos <= end.
   */
  bool contains(int genomic_pos);

  vector3 pos;       /**< 3D spatial position of this bead. */
  int genomic_pos;   /**< Representative genomic position (bp). */
  int start, end;    /**< Genomic start and end positions (bp). */
  char orientation;  /**< CTCF motif orientation ('L', 'R', 'N'). */

  int parent;        /**< Index of parent cluster in the global clusters array. */
  int level;         /**< Hierarchy level: 0=chromosome, 1=segment, 2=anchor, 3=subanchor. */

  int base_start;    /**< Start index of base-level beads covered by this cluster. */
  int base_end;      /**< End index of base-level beads covered by this cluster. */

  std::vector<int> arcs;      /**< Indices of interaction arcs touching this cluster. */
  std::vector<int> siblings;  /**< Indices of sibling clusters (same parent). */
  std::vector<int> children;  /**< Indices of child clusters in the hierarchy. */

  bool is_fixed;        /**< If true, position is fixed during Monte Carlo simulation. */
  double dist_to_next;  /**< Expected spatial distance to the next sequential cluster. */
};

#endif /* CLUSTER_H_ */
