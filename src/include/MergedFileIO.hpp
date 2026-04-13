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


// ============================================================================
// Implementation
// ============================================================================

#include <stdio.h>
#include <string.h>
#include <map>
#include <algorithm>

// Helper: trim trailing newline/carriage return
static void trim_newline(char* s) {
    size_t len = strlen(s);
    while (len > 0 && (s[len - 1] == '\n' || s[len - 1] == '\r'))
        s[--len] = '\0';
}

inline bool MergedFileIO::writeMerged(const std::string& bed_path,
                                const std::string& bedpe_path,
                                const std::string& output_path) {
    // First pass: read BED anchors
    FILE* bed_f = fopen(bed_path.c_str(), "r");
    if (!bed_f) {
        printf("[MergedFileIO] Could not open BED file: %s\n", bed_path.c_str());
        return false;
    }

    std::vector<MergedAnchor> anchors;
    char line[4096];
    while (fgets(line, sizeof(line), bed_f)) {
        trim_newline(line);
        if (line[0] == '#' || line[0] == '\0')
            continue;

        char chr[256];
        int start, end;
        if (sscanf(line, "%255s %d %d", chr, &start, &end) >= 3) {
            MergedAnchor a;
            a.chr = chr;
            a.start = start;
            a.end = end;
            a.frequency = 0;
            anchors.push_back(a);
        }
    }
    fclose(bed_f);

    // Build a lookup for anchors: key = chr:start:end -> index
    std::map<std::string, int> anchor_lookup;
    for (int i = 0; i < (int)anchors.size(); i++) {
        char key[512];
        snprintf(key, sizeof(key), "%s:%d:%d",
                 anchors[i].chr.c_str(), anchors[i].start, anchors[i].end);
        anchor_lookup[key] = i;
    }

    // Read BEDPE loops
    FILE* bedpe_f = fopen(bedpe_path.c_str(), "r");
    if (!bedpe_f) {
        printf("[MergedFileIO] Could not open BEDPE file: %s\n", bedpe_path.c_str());
        return false;
    }

    std::vector<MergedLoop> loops;
    while (fgets(line, sizeof(line), bedpe_f)) {
        trim_newline(line);
        if (line[0] == '#' || line[0] == '\0')
            continue;

        char chr1[256], chr2[256];
        int s1, e1, s2, e2, count;
        double score = 0.0;
        int nread = sscanf(line, "%255s %d %d %255s %d %d %d %lf",
                           chr1, &s1, &e1, chr2, &s2, &e2, &count, &score);
        if (nread >= 7) {
            MergedLoop loop;
            loop.chr1 = chr1;
            loop.start1 = s1;
            loop.end1 = e1;
            loop.chr2 = chr2;
            loop.start2 = s2;
            loop.end2 = e2;
            loop.count = count;
            loop.score = (nread >= 8) ? score : 0.0;
            loops.push_back(loop);

            // Update anchor frequencies
            char key1[512], key2[512];
            snprintf(key1, sizeof(key1), "%s:%d:%d", chr1, s1, e1);
            snprintf(key2, sizeof(key2), "%s:%d:%d", chr2, s2, e2);
            if (anchor_lookup.find(key1) != anchor_lookup.end())
                anchors[anchor_lookup[key1]].frequency++;
            if (anchor_lookup.find(key2) != anchor_lookup.end())
                anchors[anchor_lookup[key2]].frequency++;
        }
    }
    fclose(bedpe_f);

    // Write merged output
    FILE* out_f = fopen(output_path.c_str(), "w");
    if (!out_f) {
        printf("[MergedFileIO] Could not open output file: %s\n", output_path.c_str());
        return false;
    }

    fprintf(out_f, "#MERGED_FORMAT v1\n");
    fprintf(out_f, "#BED_SOURCE %s\n", bed_path.c_str());
    fprintf(out_f, "#BEDPE_SOURCE %s\n", bedpe_path.c_str());

    for (const MergedLoop& l : loops) {
        fprintf(out_f, "LOOP\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%.6f\n",
                l.chr1.c_str(), l.start1, l.end1,
                l.chr2.c_str(), l.start2, l.end2,
                l.count, l.score);
    }

    for (const MergedAnchor& a : anchors) {
        fprintf(out_f, "ANCHOR\t%s\t%d\t%d\t%d\n",
                a.chr.c_str(), a.start, a.end, a.frequency);
    }

    fclose(out_f);
    printf("[MergedFileIO] Wrote merged file: %s (%d loops, %d anchors)\n",
           output_path.c_str(), (int)loops.size(), (int)anchors.size());
    return true;
}

inline bool MergedFileIO::readMerged(const std::string& path,
                               std::vector<MergedLoop>& loops,
                               std::vector<MergedAnchor>& anchors) {
    loops.clear();
    anchors.clear();

    FILE* f = fopen(path.c_str(), "r");
    if (!f) {
        printf("[MergedFileIO] Could not open merged file: %s\n", path.c_str());
        return false;
    }

    char line[4096];
    while (fgets(line, sizeof(line), f)) {
        trim_newline(line);
        if (line[0] == '#' || line[0] == '\0')
            continue;

        if (strncmp(line, "LOOP\t", 5) == 0) {
            char chr1[256], chr2[256];
            int s1, e1, s2, e2, count;
            double score = 0.0;
            int nread = sscanf(line + 5, "%255s %d %d %255s %d %d %d %lf",
                               chr1, &s1, &e1, chr2, &s2, &e2, &count, &score);
            if (nread >= 7) {
                MergedLoop loop;
                loop.chr1 = chr1;
                loop.start1 = s1;
                loop.end1 = e1;
                loop.chr2 = chr2;
                loop.start2 = s2;
                loop.end2 = e2;
                loop.count = count;
                loop.score = (nread >= 8) ? score : 0.0;
                loops.push_back(loop);
            }
        } else if (strncmp(line, "ANCHOR\t", 7) == 0) {
            char chr[256];
            int start, end, freq;
            if (sscanf(line + 7, "%255s %d %d %d", chr, &start, &end, &freq) >= 4) {
                MergedAnchor a;
                a.chr = chr;
                a.start = start;
                a.end = end;
                a.frequency = freq;
                anchors.push_back(a);
            }
        }
    }

    fclose(f);
    printf("[MergedFileIO] Read merged file: %s (%d loops, %d anchors)\n",
           path.c_str(), (int)loops.size(), (int)anchors.size());
    return true;
}

inline void MergedFileIO::filterByChromosomes(std::vector<MergedLoop>& loops,
                                        std::vector<MergedAnchor>& anchors,
                                        const std::vector<std::string>& chromosomes) {
    if (chromosomes.empty())
        return;

    // Build set for fast lookup
    std::map<std::string, bool> chr_set;
    for (const std::string& c : chromosomes)
        chr_set[c] = true;

    // Filter loops: keep only those where both endpoints are in the set
    std::vector<MergedLoop> filtered_loops;
    for (const MergedLoop& l : loops) {
        if (chr_set.find(l.chr1) != chr_set.end() &&
            chr_set.find(l.chr2) != chr_set.end()) {
            filtered_loops.push_back(l);
        }
    }
    loops = filtered_loops;

    // Filter anchors
    std::vector<MergedAnchor> filtered_anchors;
    for (const MergedAnchor& a : anchors) {
        if (chr_set.find(a.chr) != chr_set.end())
            filtered_anchors.push_back(a);
    }
    anchors = filtered_anchors;

    printf("[MergedFileIO] After chromosome filter: %d loops, %d anchors\n",
           (int)loops.size(), (int)anchors.size());
}

inline bool MergedFileIO::writeBed(const std::vector<MergedAnchor>& anchors,
                             const std::string& path) {
    FILE* f = fopen(path.c_str(), "w");
    if (!f) {
        printf("[MergedFileIO] Could not open BED output: %s\n", path.c_str());
        return false;
    }

    for (const MergedAnchor& a : anchors) {
        fprintf(f, "%s\t%d\t%d\n", a.chr.c_str(), a.start, a.end);
    }

    fclose(f);
    return true;
}

inline bool MergedFileIO::writeBedpe(const std::vector<MergedLoop>& loops,
                               const std::string& path) {
    FILE* f = fopen(path.c_str(), "w");
    if (!f) {
        printf("[MergedFileIO] Could not open BEDPE output: %s\n", path.c_str());
        return false;
    }

    for (const MergedLoop& l : loops) {
        fprintf(f, "%s\t%d\t%d\t%s\t%d\t%d\t%d\n",
                l.chr1.c_str(), l.start1, l.end1,
                l.chr2.c_str(), l.start2, l.end2,
                l.count);
    }

    fclose(f);
    return true;
}

#endif /* MERGED_FILE_IO_H_ */
