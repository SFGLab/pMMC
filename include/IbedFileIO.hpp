#ifndef IBEDFILEIO_H_
#define IBEDFILEIO_H_

#include <MergedFileIO.hpp>
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


// ============================================================================
// Implementation
// ============================================================================

#include <stdio.h>
#include <string.h>
#include <map>
#include <set>
#include <algorithm>

// Parse a comma-separated genomic coordinate: "chrN,start,end"
// Returns true if parsed successfully.
static bool parseCoord(const char* token, std::string& chr, int& start, int& end) {
    // Find first comma
    const char* c1 = strchr(token, ',');
    if (!c1) return false;
    // Find second comma
    const char* c2 = strchr(c1 + 1, ',');
    if (!c2) return false;

    chr.assign(token, c1 - token);
    start = atoi(c1 + 1);
    end = atoi(c2 + 1);
    return !chr.empty() && end > start;
}

inline bool IbedFileIO::readIbed(const std::string& path,
                           std::vector<MergedLoop>& loops,
                           std::vector<MergedAnchor>& anchors) {
    loops.clear();
    anchors.clear();

    FILE* f = fopen(path.c_str(), "r");
    if (!f) {
        printf("[IbedFileIO] Could not open ibed file: %s\n", path.c_str());
        return false;
    }

    char line[8192];
    int line_num = 0;
    int parsed = 0;
    int skipped = 0;

    // Track unique anchors: key = "chr:start:end"
    std::map<std::string, int> anchor_index;

    while (fgets(line, sizeof(line), f)) {
        line_num++;
        // Trim trailing whitespace
        size_t len = strlen(line);
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r' || line[len - 1] == ' '))
            line[--len] = '\0';

        if (len == 0) continue;

        // Skip header line (starts with "bait_frag" or any non-chr token)
        if (line_num == 1 && strncmp(line, "chr", 3) != 0 && strncmp(line, "Chr", 3) != 0) {
            continue;
        }

        // Split on tabs to get bait_frag, other_frag, N_reads, score
        char bait[1024] = {0};
        char other[1024] = {0};
        int n_reads = 0;
        double score = 0.0;

        // Parse first 4 tab-separated fields
        const char* p = line;
        // Field 1: bait_frag
        const char* tab1 = strchr(p, '\t');
        if (!tab1) { skipped++; continue; }
        size_t f1_len = tab1 - p;
        if (f1_len >= sizeof(bait)) f1_len = sizeof(bait) - 1;
        memcpy(bait, p, f1_len);
        bait[f1_len] = '\0';

        // Field 2: other_frag
        p = tab1 + 1;
        const char* tab2 = strchr(p, '\t');
        if (!tab2) { skipped++; continue; }
        size_t f2_len = tab2 - p;
        if (f2_len >= sizeof(other)) f2_len = sizeof(other) - 1;
        memcpy(other, p, f2_len);
        other[f2_len] = '\0';

        // Field 3: N_reads
        p = tab2 + 1;
        n_reads = atoi(p);

        // Field 4: score (optional)
        const char* tab3 = strchr(p, '\t');
        if (tab3) {
            score = atof(tab3 + 1);
        }

        // Parse coordinates from bait and other fragments
        std::string chr1, chr2;
        int s1, e1, s2, e2;
        if (!parseCoord(bait, chr1, s1, e1)) { skipped++; continue; }
        if (!parseCoord(other, chr2, s2, e2)) { skipped++; continue; }

        // Create loop
        MergedLoop loop;
        loop.chr1 = chr1;
        loop.start1 = s1;
        loop.end1 = e1;
        loop.chr2 = chr2;
        loop.start2 = s2;
        loop.end2 = e2;
        loop.count = n_reads;
        loop.score = score;
        loops.push_back(loop);
        parsed++;

        // Track anchors
        char key1[512], key2[512];
        snprintf(key1, sizeof(key1), "%s:%d:%d", chr1.c_str(), s1, e1);
        snprintf(key2, sizeof(key2), "%s:%d:%d", chr2.c_str(), s2, e2);

        if (anchor_index.find(key1) == anchor_index.end()) {
            anchor_index[key1] = (int)anchors.size();
            MergedAnchor a;
            a.chr = chr1;
            a.start = s1;
            a.end = e1;
            a.frequency = 0;
            anchors.push_back(a);
        }
        anchors[anchor_index[key1]].frequency++;

        if (anchor_index.find(key2) == anchor_index.end()) {
            anchor_index[key2] = (int)anchors.size();
            MergedAnchor a;
            a.chr = chr2;
            a.start = s2;
            a.end = e2;
            a.frequency = 0;
            anchors.push_back(a);
        }
        anchors[anchor_index[key2]].frequency++;
    }

    fclose(f);
    printf("[IbedFileIO] Read ibed file: %s (%d loops, %d anchors, %d skipped lines)\n",
           path.c_str(), parsed, (int)anchors.size(), skipped);
    return parsed > 0;
}

inline bool IbedFileIO::readIbedFiltered(const std::string& path,
                                    std::vector<MergedLoop>& loops,
                                    std::vector<MergedAnchor>& anchors,
                                    const std::vector<std::string>& chromosomes) {
    if (!readIbed(path, loops, anchors))
        return false;

    if (chromosomes.empty())
        return true;

    MergedFileIO::filterByChromosomes(loops, anchors, chromosomes);
    return true;
}

#endif /* IBEDFILEIO_H_ */
