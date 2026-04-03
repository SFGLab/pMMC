/**
 * @file Cluster.hpp
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
#include <InteractionArc.hpp>

/**
 * @enum AnnotationType
 * @brief Epigenomic annotation type for a cluster/anchor.
 */
enum AnnotationType {
  ANNOT_NONE = 0,       /**< No annotation. */
  ANNOT_CTCF = 1,       /**< CTCF binding site. */
  ANNOT_ENHANCER = 2,   /**< Enhancer region. */
  ANNOT_PROMOTER = 3,   /**< Promoter region. */
  ANNOT_EPIGENETIC = 4  /**< Other epigenetic mark. */
};

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

  AnnotationType annotation_type; /**< Epigenomic annotation type for this cluster. */
};


// ============================================================================
// Implementation
// ============================================================================


inline Cluster::Cluster() { init(); }

inline Cluster::Cluster(int start, int end) {
  init();
  this->start = start;
  this->end = end;
  genomic_pos =
      (start + end) / 2; // Pozycja genomiczna jako położenie środka klastra
}

inline void Cluster::init() {
  start = 0;
  end = 0;
  genomic_pos = 0;
  orientation = 'N';
  base_start = 0;
  base_end = 0;
  parent = -1;
  is_fixed = false;
  level = 0;
  dist_to_next = 0.0;
  annotation_type = ANNOT_NONE;
  pos.set(0.0f, 0.0f, 0.0f);
}

inline void Cluster::print() {
  printf("pos=%d (%d-%d, %d-%d), lvl=%d, par=%d, next=%lf, pts=(%f %f %f)",
         genomic_pos, start, end, base_start, base_end, level, parent,
         dist_to_next, pos.x, pos.y, pos.z);
  if (orientation != 'N')
    printf(" %c", orientation);

  if (is_fixed)
    printf(" fixed");

  if (siblings.size() > 0) {
    printf(", sibl=[");
    for (unsigned int i = 0; i < siblings.size(); ++i)
      printf("%d ", siblings[i]);
    printf("]");
  }

  if (arcs.size() > 0) {
    printf(", arcs=[");
    for (unsigned int i = 0; i < arcs.size(); ++i)
      printf("%d ", arcs[i]);
    printf("]");
  }

  if (children.size() > 0) {
    printf(", children=[");
    printv(children, false);
    printf("]");
  }

  printf("\n");
}

inline bool Cluster::contains(int genomic_pos) {
  return genomic_pos >= start && genomic_pos <= end;
}

inline void Cluster::toFile(FILE *file) {
  fprintf(file, "%d %d %c %f %f %f %zu", start, end, orientation, pos.x, pos.y,
          pos.z, children.size());
  for (unsigned int i = 0; i < children.size(); ++i)
    fprintf(file, " %d", children[i]);
  fprintf(file, "\n");
}

inline void Cluster::toFilePreviousFormat(FILE *file) {
  fprintf(file, "%d %d %d %f %f %f %zu", (start + end) / 2, start, end, pos.x,
          pos.y, pos.z, children.size());
  for (unsigned int i = 0; i < children.size(); ++i)
    fprintf(file, " %d", children[i]);
  fprintf(file, "\n");
}

inline void Cluster::fromFile(FILE *file) {
  int st, end;
  float x, y, z;
  int children_cnt, tmp;
  char c;

  fscanf(file, "%d %d %c %f %f %f %d", &st, &end, &c, &x, &y, &z,
         &children_cnt);

  for (int i = 0; i < children_cnt; ++i) {
    fscanf(file, "%d", &tmp);
    children.push_back(tmp);
  }

  this->start = st;
  this->end = end;
  this->genomic_pos = (st + end) / 2;
  this->orientation = c;
  pos.set(x, y, z);
}

#endif /* CLUSTER_H_ */
