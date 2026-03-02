/**
 * @file BedRegions.h
 * @brief Collection of BedRegion intervals loaded from a BED file.
 *
 * Used to represent centromere regions, predefined segment splits,
 * and other sets of genomic intervals (e.g. avoid regions).
 */

#ifndef BEDREGIONS_H_
#define BEDREGIONS_H_

#include <stdio.h>
#include <string>
#include <vector>

#include <BedRegion.h>

using namespace std;

/**
 * @class BedRegions
 * @brief A list of BedRegion genomic intervals.
 *
 * Regions can be loaded from a tab-delimited BED file or generated
 * programmatically by splitting a chromosome range into fixed-size bins.
 */
class BedRegions {
public:
  /** @brief Default constructor (empty list). */
  BedRegions();

  /**
   * @brief Load regions from a BED file (tab-delimited: chr start end).
   * @param path File path to the BED file.
   */
  void fromFile(std::string path);

  /** @brief Print all regions to stdout. */
  void print();

  /**
   * @brief Generate evenly spaced intervals and append them to the list.
   * @param chr Chromosome name.
   * @param start Start position (bp).
   * @param end End position (bp).
   * @param step Interval size (bp).
   */
  void addNewIntervals(std::string chr, int start, int end, int step);

  std::vector<BedRegion> regions;  /**< The stored list of genomic regions. */
};

#endif /* BEDREGIONS_H_ */
