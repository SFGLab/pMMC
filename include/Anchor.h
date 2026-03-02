/**
 * @file Anchor.h
 * @brief ChIA-PET anchor representing a protein binding site on the genome.
 *
 * An anchor is a genomic interval where a ChIA-PET factor (e.g. CTCF, RNAPII)
 * binds. Anchors are the fundamental units that define interaction block
 * boundaries in the 4-level hierarchy: Chromosome > Segment > Anchor > Subanchor.
 */

#ifndef ANCHOR_H_
#define ANCHOR_H_

#include <InteractionArc.h>
#include <string>

/**
 * @class Anchor
 * @brief A single ChIA-PET anchor (protein binding site) on a chromosome.
 *
 * Anchors are derived from merged PET cluster endpoints. Each anchor has a
 * genomic interval [start, end) and a CTCF motif orientation ('L', 'R', or 'N').
 * The center is computed as (start + end) / 2.
 */
class Anchor {
public:
  /** @brief Default constructor. */
  Anchor();

  /**
   * @brief Construct an anchor at the given genomic interval.
   * @param chr Chromosome name (e.g. "chr22").
   * @param start Start position in base pairs.
   * @param end End position in base pairs.
   * @param orientation CTCF motif orientation: 'L' (left/forward),
   *        'R' (right/reverse), or 'N' (not available).
   */
  Anchor(std::string chr, int start, int end, char orientation);

  /**
   * @brief Get the length of the anchor interval.
   * @return end - start, in base pairs.
   */
  int length();

  /**
   * @brief Check whether a genomic position falls within this anchor.
   * @param pos Genomic position in base pairs.
   * @return True if start <= pos <= end.
   */
  bool contains(int pos);

  std::string chr;    /**< Chromosome name (e.g. "chr22"). */
  int start;          /**< Start position in base pairs. */
  int end;            /**< End position in base pairs. */
  int center;         /**< Center position: (start + end) / 2. */
  char orientation;   /**< CTCF motif orientation: 'L', 'R', or 'N' (not available). */
};

#endif /* ANCHOR_H_ */
