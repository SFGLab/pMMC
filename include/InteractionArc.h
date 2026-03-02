/**
 * @file InteractionArc.h
 * @brief An interaction (arc) between two genomic regions from ChIA-PET data.
 *
 * Arcs represent chromatin contacts detected by PET clusters. Each arc
 * connects two anchors (identified by index) with a PET count (score).
 *
 * Initially, start/end refer to local anchor indices on a specific chromosome.
 * After LooperSolver::createTree(), they are shifted to global cluster indices.
 *
 * When multiple arcs connect the same two regions with different factors,
 * a summary arc (factor=-1) is created with eff_score = sum of all scores.
 * The individual arcs then have eff_score=0, preserving the original score
 * for display purposes.
 */

#ifndef INTERACTIONARC_H_
#define INTERACTIONARC_H_

#include <stdio.h>
#include <string>

/**
 * @class InteractionArc
 * @brief A single chromatin interaction arc between two genomic loci.
 *
 * Represents a PET cluster linking two anchors. The score field holds
 * the raw PET count, while eff_score is used during reconstruction
 * (may be aggregated across factors).
 */
class InteractionArc {
public:
  /** @brief Default constructor (calls init()). */
  InteractionArc();

  /**
   * @brief Construct an arc between two anchor indices.
   * @param start Index of the first anchor (or cluster).
   * @param end Index of the second anchor (or cluster).
   * @param score PET count (interaction strength).
   * @param factor Factor index (-1 for summary arcs, >=0 for specific factors).
   */
  InteractionArc(int start, int end, int score, int factor = -1);

  /**
   * @brief Write arc data to a binary file.
   * @param file Open file handle for writing.
   */
  void toFile(FILE *file);

  /**
   * @brief Read arc data from a binary file.
   * @param file Open file handle for reading.
   */
  void fromFile(FILE *file);

  /** @brief Reset all fields to default values. */
  void init();

  /** @brief Print arc info to stdout. */
  void print();

  /**
   * @brief Genomic span of the arc (genomic_end - genomic_start).
   * @return Arc length in base pairs.
   */
  int length();

  /**
   * @brief Comparison operator for sorting arcs by (start, end, factor).
   * @param arc The other arc to compare against.
   * @return True if this arc sorts before the other.
   */
  bool operator<(const InteractionArc &arc) const {
    if (start < arc.start)
      return true;
    if (start == arc.start) {
      if (end < arc.end)
        return true;
      if (end == arc.end)
        return factor < arc.factor;
    }
    return false;
  }

  int start;         /**< First anchor/cluster index (local then global). */
  int end;           /**< Second anchor/cluster index (local then global). */
  int genomic_start; /**< Genomic position of the first anchor (bp). */
  int genomic_end;   /**< Genomic position of the second anchor (bp). */
  int score;         /**< Raw PET count (interaction strength). */
  int eff_score;     /**< Effective score used during reconstruction (may be aggregated). */
  int factor;        /**< Factor index (-1 = summary arc, >=0 = specific factor). */
};

#endif /* INTERACTIONARC_H_ */
