#ifndef MERGED_FILE_IO_H_
#define MERGED_FILE_IO_H_

#include <string>
#include <vector>

struct MergedLoop {
    std::string chr1, chr2;
    int start1, end1, start2, end2;
    int count;
    double score;
};

struct MergedAnchor {
    std::string chr;
    int start, end;
    int frequency;
};

class MergedFileIO {
public:
    // Write merged file from BED + BEDPE
    static bool writeMerged(const std::string& bed_path,
                            const std::string& bedpe_path,
                            const std::string& output_path);

    // Read merged file
    static bool readMerged(const std::string& path,
                           std::vector<MergedLoop>& loops,
                           std::vector<MergedAnchor>& anchors);

    // Filter by chromosomes
    static void filterByChromosomes(std::vector<MergedLoop>& loops,
                                    std::vector<MergedAnchor>& anchors,
                                    const std::vector<std::string>& chromosomes);

    // Write BED file from anchors
    static bool writeBed(const std::vector<MergedAnchor>& anchors,
                         const std::string& path);

    // Write BEDPE file from loops
    static bool writeBedpe(const std::vector<MergedLoop>& loops,
                           const std::string& path);
};

#endif /* MERGED_FILE_IO_H_ */
