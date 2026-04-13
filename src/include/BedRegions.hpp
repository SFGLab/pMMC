/**
 * @file BedRegions.hpp
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

#include <BedRegion.hpp>

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


// ============================================================================
// Implementation
// ============================================================================


inline BedRegions::BedRegions() {}

inline void BedRegions::fromFile(std::string filename) {
  FILE *f;
  f = fopen(filename.c_str(), "r");
  if (f == NULL) {
    printf("Error opening file [%s]!\n", filename.c_str());
    return;
  }

  char chr[16], line[100];
  int start, end;

  while (!feof(f)) {
    if (fscanf(f, "%s %d %d", chr, &start, &end) != 3)
      continue;
    // printf("%s %d %d\n", chr, start, end);
    fgets(line, 100, f); // read to the end of line

    BedRegion b((std::string)chr, start, end);
    regions.push_back(b);
  }
  fclose(f);
}

inline void BedRegions::print() {
  printf("regions: %d\n", (int)regions.size());
  for (unsigned int i = 0; i < regions.size(); ++i) {
    printf("[%s] %d %d\n", regions[i].chr.c_str(), regions[i].start,
           regions[i].end);
  }
}

inline void BedRegions::addNewIntervals(std::string chr, int start, int end,
                                 int step) {

  if (step < 1) {
    printf("Step must be at least 1!\n");
    return;
  }

  while (start < end) {
    BedRegion reg(chr, start, start);
    regions.push_back(reg);
    start += step;
  }
}

#endif /* BEDREGIONS_H_ */
