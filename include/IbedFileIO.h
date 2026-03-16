#ifndef IBEDFILEIO_H_
#define IBEDFILEIO_H_

#include <MergedFileIO.h>
#include <string>
#include <vector>

/**
 * @brief Parser for .ibed (interaction BED) files.
 *
 * ibed format: tab-separated columns where the first two columns contain
 * comma-separated genomic coordinates: "chrN,start,end"
 *
 * Header: bait_frag\tother_frag\tN_reads\tscore\t...
 * Data:   chr5,687730,690639\tchr5,789001,795437\t564\t79.66\t...
 *
 * Extracts loops (pairwise interactions) and anchors (unique endpoints).
 */
class IbedFileIO {
public:
    /**
     * @brief Read an ibed file and extract loops and anchors.
     * @param path Path to .ibed file.
     * @param loops Output vector of MergedLoop (pairwise interactions).
     * @param anchors Output vector of MergedAnchor (unique anchor regions).
     * @return true on success, false on failure.
     */
    static bool readIbed(const std::string& path,
                         std::vector<MergedLoop>& loops,
                         std::vector<MergedAnchor>& anchors);

    /**
     * @brief Read ibed and filter by chromosomes.
     * @param path Path to .ibed file.
     * @param loops Output loops.
     * @param anchors Output anchors.
     * @param chromosomes List of chromosome names to keep (e.g., {"chr1","chr8"}).
     * @return true on success.
     */
    static bool readIbedFiltered(const std::string& path,
                                  std::vector<MergedLoop>& loops,
                                  std::vector<MergedAnchor>& anchors,
                                  const std::vector<std::string>& chromosomes);
};

#endif /* IBEDFILEIO_H_ */
