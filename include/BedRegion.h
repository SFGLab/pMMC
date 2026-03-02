/**
 * @file BedRegion.h
 * @brief A single genomic region in BED format (chr:start-end).
 *
 * Used to represent centromeres, avoid regions, segment boundaries,
 * and user-specified regions of interest for selective reconstruction.
 */

#ifndef BEDREGION_H_
#define BEDREGION_H_

#include <stdio.h>
#include <string>

/**
 * @class BedRegion
 * @brief Genomic interval in BED format: chromosome, start, end.
 *
 * Supports parsing from "chr:start-end" strings and containment queries.
 * An empty BedRegion (default-constructed) has chr="" and start=end=0.
 */
class BedRegion {
public:
  /**
   * @brief Construct a region with given coordinates.
   * @param _chr Chromosome name (e.g. "chr22").
   * @param _start Start position in base pairs.
   * @param _end End position in base pairs.
   */
  BedRegion(std::string _chr, int _start, int _end);

  /** @brief Default constructor (empty region). */
  BedRegion();

  /**
   * @brief Check if a string can be parsed as a BED region.
   * @param str Input string in "chr:start-end" format.
   * @return True if the string is a valid BED region.
   */
  static bool tryParse(std::string str);

  /**
   * @brief Parse a "chr:start-end" string into this region.
   * @param str Input string.
   * @return True on successful parse.
   */
  bool parse(std::string str);

  /** @brief Print this region to stdout. */
  void print();

  /**
   * @brief Check whether a genomic position falls within this region.
   * @param pos Position in base pairs.
   * @return True if start <= pos <= end.
   */
  bool contains(int pos);

  std::string chr;  /**< Chromosome name. */
  int start;        /**< Start position (bp). */
  int end;          /**< End position (bp). */
};

#endif /* BEDREGION_H_ */
