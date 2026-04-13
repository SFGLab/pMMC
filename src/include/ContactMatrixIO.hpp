#ifndef CONTACTMATRIXIO_H_
#define CONTACTMATRIXIO_H_

#include <MergedFileIO.hpp>
#include <string>
#include <vector>
#include <set>

/**
 * @brief Parser for Hi-C, pcHi-C, and HiChIP contact matrix data.
 *
 * Supports three tab-delimited text formats plus binary mcool.
 * All readers support streaming chromosome filtering so that only
 * matching rows are stored in memory -- critical for files that can
 * exceed 90 GB.
 */
class ContactMatrixIO {
public:
    enum class Format {
        UNKNOWN,
        PAIRS,      // chr1 pos1 chr2 pos2 count
        COOL_DUMP,  // bin1_id chr1 start1 end1 bin2_id chr2 start2 end2 count
        HICPRO,     // chr start end chr start end count
        MCOOL_BINARY // Binary .mcool/.cool file (auto-dump via cooler)
    };

    static Format detectFormat(const std::string& path);

    static bool readContactMatrix(const std::string& path,
                                  std::vector<MergedLoop>& loops,
                                  std::vector<MergedAnchor>& anchors,
                                  int resolution = 5000);

    static bool readPairs(const std::string& path,
                          std::vector<MergedLoop>& loops,
                          std::vector<MergedAnchor>& anchors,
                          int resolution,
                          const std::set<std::string>& chr_filter = {});

    static bool readCoolDump(const std::string& path,
                             std::vector<MergedLoop>& loops,
                             std::vector<MergedAnchor>& anchors,
                             const std::set<std::string>& chr_filter = {});

    static bool readHiCPro(const std::string& path,
                           std::vector<MergedLoop>& loops,
                           std::vector<MergedAnchor>& anchors,
                           const std::set<std::string>& chr_filter = {});

    /**
     * @brief Streaming filtered read -- chromosome filter is applied during
     *        parsing so non-matching rows are never stored.
     */
    static bool readContactMatrixFiltered(const std::string& path,
                                          std::vector<MergedLoop>& loops,
                                          std::vector<MergedAnchor>& anchors,
                                          const std::vector<std::string>& chromosomes,
                                          int resolution = 5000);

    static const char* formatName(Format fmt);
};


// ============================================================================
// Implementation
// ============================================================================

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <map>
#include <set>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Helper: trim trailing whitespace in-place
// ---------------------------------------------------------------------------
static void trimLine(char* line) {
    size_t len = strlen(line);
    while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r' ||
                       line[len - 1] == ' ' || line[len - 1] == '\t'))
        line[--len] = '\0';
}

// ---------------------------------------------------------------------------
// Helper: count tab-separated fields in a line
// ---------------------------------------------------------------------------
static int countFields(const char* line) {
    if (!line || line[0] == '\0') return 0;
    int n = 1;
    for (const char* p = line; *p; ++p)
        if (*p == '\t') ++n;
    return n;
}

// ---------------------------------------------------------------------------
// Helper: check if a string looks like an integer (digits, possibly negative)
// ---------------------------------------------------------------------------
static bool looksLikeInt(const char* s) {
    if (!s || *s == '\0') return false;
    if (*s == '-') ++s;
    if (*s == '\0') return false;
    while (*s) {
        if (*s < '0' || *s > '9') return false;
        ++s;
    }
    return true;
}

// ---------------------------------------------------------------------------
// Helper: check if a token looks like a chromosome name
// ---------------------------------------------------------------------------
static bool looksLikeChr(const char* s) {
    if (!s || *s == '\0') return false;
    if (strncmp(s, "chr", 3) == 0 || strncmp(s, "Chr", 3) == 0) return true;
    size_t len = strlen(s);
    if (len <= 2 && ((s[0] >= '1' && s[0] <= '9') || s[0] == 'X' || s[0] == 'Y'))
        return true;
    if (len == 2 && s[0] == 'M' && s[1] == 'T') return true;
    return false;
}

// ---------------------------------------------------------------------------
// Helper: register an anchor in the anchor map, return its index
// ---------------------------------------------------------------------------
static int registerAnchor(const std::string& chr, int start, int end,
                          std::map<std::string, int>& anchor_index,
                          std::vector<MergedAnchor>& anchors) {
    char key[512];
    snprintf(key, sizeof(key), "%s:%d:%d", chr.c_str(), start, end);
    std::string skey(key);
    auto it = anchor_index.find(skey);
    if (it == anchor_index.end()) {
        int idx = (int)anchors.size();
        anchor_index[skey] = idx;
        MergedAnchor a;
        a.chr = chr;
        a.start = start;
        a.end = end;
        a.frequency = 0;
        anchors.push_back(a);
        return idx;
    }
    return it->second;
}

// ---------------------------------------------------------------------------
// Helper: split a line into tokens by tabs. Returns count of tokens stored.
// ---------------------------------------------------------------------------
static int splitTabs(char* line, char** tokens, int max_tokens) {
    int n = 0;
    char* p = line;
    while (n < max_tokens) {
        tokens[n++] = p;
        char* tab = strchr(p, '\t');
        if (!tab) break;
        *tab = '\0';
        p = tab + 1;
    }
    return n;
}

// ---------------------------------------------------------------------------
// Helper: check whether a chromosome passes the filter.
// Empty filter means accept everything.
// ---------------------------------------------------------------------------
static inline bool chrAccepted(const std::string& chr,
                               const std::set<std::string>& filter) {
    return filter.empty() || filter.count(chr) > 0;
}

// ---------------------------------------------------------------------------
// Progress reporting every 10M lines
// ---------------------------------------------------------------------------
static const long long PROGRESS_INTERVAL = 10000000LL;

static inline void reportProgress(long long line_num, int parsed,
                                  long long skipped) {
    if (line_num > 0 && (line_num % PROGRESS_INTERVAL) == 0) {
        printf("[ContactMatrixIO]   ... %lld M lines scanned, %d contacts kept, "
               "%lld filtered/skipped\n",
               line_num / 1000000LL, parsed, skipped);
        fflush(stdout);
    }
}

// ============================================================================
// Format detection
// ============================================================================

inline ContactMatrixIO::Format ContactMatrixIO::detectFormat(const std::string& path) {
    // Check for binary HDF5/mcool/cool format by reading magic bytes
    {
        FILE* bf = fopen(path.c_str(), "rb");
        if (bf) {
            unsigned char magic[8] = {0};
            size_t nr = fread(magic, 1, 8, bf);
            fclose(bf);
            if (nr >= 8 && magic[0] == 0x89 && magic[1] == 'H' &&
                magic[2] == 'D' && magic[3] == 'F') {
                return Format::MCOOL_BINARY;
            }
        }
    }

    FILE* f = fopen(path.c_str(), "r");
    if (!f) return Format::UNKNOWN;

    char line[8192];
    Format detected = Format::UNKNOWN;

    while (fgets(line, sizeof(line), f)) {
        trimLine(line);
        if (line[0] == '\0') continue;
        if (line[0] == '#') continue;

        int nfields = countFields(line);

        char copy[8192];
        strncpy(copy, line, sizeof(copy) - 1);
        copy[sizeof(copy) - 1] = '\0';

        char* tokens[16];
        int ntok = splitTabs(copy, tokens, 16);

        if (ntok >= 5) {
            if (nfields >= 9 && ntok >= 9 &&
                looksLikeInt(tokens[0]) && looksLikeChr(tokens[1]) &&
                looksLikeInt(tokens[2]) && looksLikeInt(tokens[3]) &&
                looksLikeInt(tokens[4]) && looksLikeChr(tokens[5])) {
                detected = Format::COOL_DUMP;
                break;
            }

            if (nfields >= 7 && ntok >= 7 &&
                looksLikeChr(tokens[0]) && looksLikeInt(tokens[1]) &&
                looksLikeInt(tokens[2]) && looksLikeChr(tokens[3]) &&
                looksLikeInt(tokens[4]) && looksLikeInt(tokens[5]) &&
                looksLikeInt(tokens[6])) {
                detected = Format::HICPRO;
                break;
            }

            if (nfields >= 5 && ntok >= 5 &&
                looksLikeChr(tokens[0]) && looksLikeInt(tokens[1]) &&
                looksLikeChr(tokens[2]) && looksLikeInt(tokens[3]) &&
                looksLikeInt(tokens[4])) {
                detected = Format::PAIRS;
                break;
            }
        }
    }

    fclose(f);
    return detected;
}

inline const char* ContactMatrixIO::formatName(Format fmt) {
    switch (fmt) {
        case Format::PAIRS:        return "pairs";
        case Format::COOL_DUMP:    return "cool_dump";
        case Format::HICPRO:       return "hicpro";
        case Format::MCOOL_BINARY: return "mcool_binary";
        default:                   return "unknown";
    }
}

// ============================================================================
// Pairs reader  (streaming chromosome filter)
// ============================================================================

inline bool ContactMatrixIO::readPairs(const std::string& path,
                                std::vector<MergedLoop>& loops,
                                std::vector<MergedAnchor>& anchors,
                                int resolution,
                                const std::set<std::string>& chr_filter) {
    loops.clear();
    anchors.clear();

    FILE* f = fopen(path.c_str(), "r");
    if (!f) {
        printf("[ContactMatrixIO] Could not open file: %s\n", path.c_str());
        return false;
    }

    if (resolution <= 0) resolution = 5000;

    char line[8192];
    long long line_num = 0;
    int parsed = 0;
    long long skipped = 0;

    std::map<std::string, int> anchor_index;

    printf("[ContactMatrixIO] Reading pairs file: %s (filter: %s)\n",
           path.c_str(),
           chr_filter.empty() ? "none" : "active");

    while (fgets(line, sizeof(line), f)) {
        line_num++;
        reportProgress(line_num, parsed, skipped);

        trimLine(line);
        if (line[0] == '\0' || line[0] == '#') continue;

        char* tokens[16];
        char copy[8192];
        strncpy(copy, line, sizeof(copy) - 1);
        copy[sizeof(copy) - 1] = '\0';
        int ntok = splitTabs(copy, tokens, 16);

        if (ntok < 5) { skipped++; continue; }
        if (!looksLikeChr(tokens[0])) { continue; }

        std::string chr1(tokens[0]);
        std::string chr2(tokens[2]);

        // Streaming filter: skip rows where neither chr matches
        if (!chrAccepted(chr1, chr_filter) && !chrAccepted(chr2, chr_filter)) {
            skipped++;
            continue;
        }

        int pos1 = atoi(tokens[1]);
        int pos2 = atoi(tokens[3]);
        int count = atoi(tokens[4]);

        if (count <= 0) { skipped++; continue; }

        int start1 = (pos1 / resolution) * resolution;
        int end1 = start1 + resolution;
        int start2 = (pos2 / resolution) * resolution;
        int end2 = start2 + resolution;

        MergedLoop loop;
        loop.chr1 = chr1;
        loop.start1 = start1;
        loop.end1 = end1;
        loop.chr2 = chr2;
        loop.start2 = start2;
        loop.end2 = end2;
        loop.count = count;
        loop.score = (double)count;
        loops.push_back(loop);
        parsed++;

        int idx1 = registerAnchor(chr1, start1, end1, anchor_index, anchors);
        anchors[idx1].frequency++;
        int idx2 = registerAnchor(chr2, start2, end2, anchor_index, anchors);
        anchors[idx2].frequency++;
    }

    fclose(f);
    printf("[ContactMatrixIO] Read pairs file: %s (%d loops, %d anchors, "
           "%lld skipped, %lld total lines, resolution=%d)\n",
           path.c_str(), parsed, (int)anchors.size(), skipped, line_num,
           resolution);
    return parsed > 0;
}

// ============================================================================
// Cool dump reader  (streaming chromosome filter)
// ============================================================================

inline bool ContactMatrixIO::readCoolDump(const std::string& path,
                                   std::vector<MergedLoop>& loops,
                                   std::vector<MergedAnchor>& anchors,
                                   const std::set<std::string>& chr_filter) {
    loops.clear();
    anchors.clear();

    FILE* f = fopen(path.c_str(), "r");
    if (!f) {
        printf("[ContactMatrixIO] Could not open file: %s\n", path.c_str());
        return false;
    }

    char line[8192];
    long long line_num = 0;
    int parsed = 0;
    long long skipped = 0;

    std::map<std::string, int> anchor_index;

    printf("[ContactMatrixIO] Reading cool dump: %s (filter: %s)\n",
           path.c_str(),
           chr_filter.empty() ? "none" : "active");

    while (fgets(line, sizeof(line), f)) {
        line_num++;
        reportProgress(line_num, parsed, skipped);

        trimLine(line);
        if (line[0] == '\0' || line[0] == '#') continue;

        char* tokens[16];
        char copy[8192];
        strncpy(copy, line, sizeof(copy) - 1);
        copy[sizeof(copy) - 1] = '\0';
        int ntok = splitTabs(copy, tokens, 16);

        if (ntok < 9) { skipped++; continue; }
        if (!looksLikeInt(tokens[0])) { continue; }

        // Streaming filter: check chr fields before parsing the rest
        std::string chr1(tokens[1]);
        std::string chr2(tokens[5]);

        if (!chrAccepted(chr1, chr_filter) && !chrAccepted(chr2, chr_filter)) {
            skipped++;
            continue;
        }

        int start1 = atoi(tokens[2]);
        int end1 = atoi(tokens[3]);
        int start2 = atoi(tokens[6]);
        int end2 = atoi(tokens[7]);
        int count = atoi(tokens[8]);

        if (count <= 0) { skipped++; continue; }
        if (end1 <= start1 || end2 <= start2) { skipped++; continue; }

        MergedLoop loop;
        loop.chr1 = chr1;
        loop.start1 = start1;
        loop.end1 = end1;
        loop.chr2 = chr2;
        loop.start2 = start2;
        loop.end2 = end2;
        loop.count = count;
        loop.score = (double)count;
        loops.push_back(loop);
        parsed++;

        int idx1 = registerAnchor(chr1, start1, end1, anchor_index, anchors);
        anchors[idx1].frequency++;
        int idx2 = registerAnchor(chr2, start2, end2, anchor_index, anchors);
        anchors[idx2].frequency++;
    }

    fclose(f);
    printf("[ContactMatrixIO] Read cool dump: %s (%d loops, %d anchors, "
           "%lld skipped, %lld total lines)\n",
           path.c_str(), parsed, (int)anchors.size(), skipped, line_num);
    return parsed > 0;
}

// ============================================================================
// HiCPro reader  (streaming chromosome filter)
// ============================================================================

inline bool ContactMatrixIO::readHiCPro(const std::string& path,
                                 std::vector<MergedLoop>& loops,
                                 std::vector<MergedAnchor>& anchors,
                                 const std::set<std::string>& chr_filter) {
    loops.clear();
    anchors.clear();

    FILE* f = fopen(path.c_str(), "r");
    if (!f) {
        printf("[ContactMatrixIO] Could not open file: %s\n", path.c_str());
        return false;
    }

    char line[8192];
    long long line_num = 0;
    int parsed = 0;
    long long skipped = 0;

    std::map<std::string, int> anchor_index;

    printf("[ContactMatrixIO] Reading HiCPro file: %s (filter: %s)\n",
           path.c_str(),
           chr_filter.empty() ? "none" : "active");

    while (fgets(line, sizeof(line), f)) {
        line_num++;
        reportProgress(line_num, parsed, skipped);

        trimLine(line);
        if (line[0] == '\0' || line[0] == '#') continue;

        char* tokens[16];
        char copy[8192];
        strncpy(copy, line, sizeof(copy) - 1);
        copy[sizeof(copy) - 1] = '\0';
        int ntok = splitTabs(copy, tokens, 16);

        if (ntok < 7) { skipped++; continue; }
        if (!looksLikeChr(tokens[0])) { continue; }

        // Streaming filter: check chr fields first
        std::string chr1(tokens[0]);
        std::string chr2(tokens[3]);

        if (!chrAccepted(chr1, chr_filter) && !chrAccepted(chr2, chr_filter)) {
            skipped++;
            continue;
        }

        int start1 = atoi(tokens[1]);
        int end1 = atoi(tokens[2]);
        int start2 = atoi(tokens[4]);
        int end2 = atoi(tokens[5]);
        int count = atoi(tokens[6]);

        if (count <= 0) { skipped++; continue; }
        if (end1 <= start1 || end2 <= start2) { skipped++; continue; }

        MergedLoop loop;
        loop.chr1 = chr1;
        loop.start1 = start1;
        loop.end1 = end1;
        loop.chr2 = chr2;
        loop.start2 = start2;
        loop.end2 = end2;
        loop.count = count;
        loop.score = (double)count;
        loops.push_back(loop);
        parsed++;

        int idx1 = registerAnchor(chr1, start1, end1, anchor_index, anchors);
        anchors[idx1].frequency++;
        int idx2 = registerAnchor(chr2, start2, end2, anchor_index, anchors);
        anchors[idx2].frequency++;
    }

    fclose(f);
    printf("[ContactMatrixIO] Read HiCPro file: %s (%d loops, %d anchors, "
           "%lld skipped, %lld total lines)\n",
           path.c_str(), parsed, (int)anchors.size(), skipped, line_num);
    return parsed > 0;
}

// ============================================================================
// Auto-detect and read  (unfiltered -- entire file)
// ============================================================================

static std::string dumpMcoolToText(const std::string& mcool_path, int resolution,
                                   const std::set<std::string>& chr_filter = {}) {
    std::string uri = mcool_path;
    size_t dot = mcool_path.rfind('.');
    bool is_mcool = (dot != std::string::npos &&
                     mcool_path.substr(dot) == ".mcool");
    if (is_mcool) {
        char buf[64];
        snprintf(buf, sizeof(buf), "::resolutions/%d", resolution);
        uri += buf;
    }

    // Use a unique temp file to avoid collisions with existing dumps
    std::string dump_path = mcool_path + ".tmp_dump.tsv";

    // Build region flags for cooler dump if chromosome filter is active
    std::string region_flags;
    if (chr_filter.size() == 1) {
        // Single chromosome — use -r for efficient extraction
        region_flags = " -r " + *chr_filter.begin() + " -r2 " + *chr_filter.begin();
    }

    char cmd[4096];
    snprintf(cmd, sizeof(cmd),
             "cooler dump --join \"%s\"%s > \"%s\" 2>&1",
             uri.c_str(), region_flags.c_str(), dump_path.c_str());

    printf("[ContactMatrixIO] Running: %s\n", cmd);
    int rc = system(cmd);
    if (rc != 0) {
        // Fallback without --join
        snprintf(cmd, sizeof(cmd),
                 "cooler dump \"%s\"%s > \"%s\" 2>&1",
                 uri.c_str(), region_flags.c_str(), dump_path.c_str());
        rc = system(cmd);
        if (rc != 0) {
            printf("[ContactMatrixIO] cooler dump fallback also failed\n");
            return "";
        }
    }

    FILE* df = fopen(dump_path.c_str(), "r");
    if (!df) return "";
    // Use _fseeki64/_ftelli64 on Windows to handle files > 2GB
#ifdef _WIN32
    _fseeki64(df, 0, SEEK_END);
    long long sz = _ftelli64(df);
#else
    fseeko(df, 0, SEEK_END);
    long long sz = ftello(df);
#endif
    fclose(df);
    if (sz <= 0) {
        printf("[ContactMatrixIO] cooler dump produced empty file\n");
        return "";
    }

    printf("[ContactMatrixIO] cooler dump wrote %lld bytes to %s\n", sz, dump_path.c_str());
    return dump_path;
}

// Internal dispatcher that accepts an optional chr filter set
static bool readByFormat(ContactMatrixIO::Format fmt, const std::string& path,
                         std::vector<MergedLoop>& loops,
                         std::vector<MergedAnchor>& anchors,
                         int resolution,
                         const std::set<std::string>& chr_filter) {
    switch (fmt) {
        case ContactMatrixIO::Format::PAIRS:
            return ContactMatrixIO::readPairs(path, loops, anchors, resolution, chr_filter);
        case ContactMatrixIO::Format::COOL_DUMP:
            return ContactMatrixIO::readCoolDump(path, loops, anchors, chr_filter);
        case ContactMatrixIO::Format::HICPRO:
            return ContactMatrixIO::readHiCPro(path, loops, anchors, chr_filter);
        default:
            return false;
    }
}

inline bool ContactMatrixIO::readContactMatrix(const std::string& path,
                                        std::vector<MergedLoop>& loops,
                                        std::vector<MergedAnchor>& anchors,
                                        int resolution) {
    Format fmt = detectFormat(path);
    if (fmt == Format::UNKNOWN) {
        printf("[ContactMatrixIO] Could not detect format of: %s\n", path.c_str());
        return false;
    }

    printf("[ContactMatrixIO] Detected format: %s\n", formatName(fmt));

    if (fmt == Format::MCOOL_BINARY) {
        std::string dump_path = dumpMcoolToText(path, resolution);
        if (dump_path.empty()) {
            printf("[ContactMatrixIO] Failed to dump binary mcool/cool file.\n");
            printf("[ContactMatrixIO] Ensure 'cooler' is installed: pip install cooler\n");
            return false;
        }
        Format dump_fmt = detectFormat(dump_path);
        printf("[ContactMatrixIO] Dump format: %s\n", formatName(dump_fmt));
        std::set<std::string> empty_filter;
        bool ok = readByFormat(dump_fmt, dump_path, loops, anchors,
                               resolution, empty_filter);
        remove(dump_path.c_str());
        return ok;
    }

    std::set<std::string> empty_filter;
    return readByFormat(fmt, path, loops, anchors, resolution, empty_filter);
}

// ============================================================================
// Streaming filtered read -- chromosome filter applied during parsing
// ============================================================================

inline bool ContactMatrixIO::readContactMatrixFiltered(const std::string& path,
                                                std::vector<MergedLoop>& loops,
                                                std::vector<MergedAnchor>& anchors,
                                                const std::vector<std::string>& chromosomes,
                                                int resolution) {
    // Build a set for O(1) lookup during streaming parse
    std::set<std::string> chr_filter(chromosomes.begin(), chromosomes.end());

    if (!chr_filter.empty()) {
        printf("[ContactMatrixIO] Streaming filter active for %d chromosome(s):",
               (int)chr_filter.size());
        for (auto& c : chr_filter) printf(" %s", c.c_str());
        printf("\n");
    }

    Format fmt = detectFormat(path);
    if (fmt == Format::UNKNOWN) {
        printf("[ContactMatrixIO] Could not detect format of: %s\n", path.c_str());
        return false;
    }

    printf("[ContactMatrixIO] Detected format: %s\n", formatName(fmt));

    if (fmt == Format::MCOOL_BINARY) {
        std::string dump_path = dumpMcoolToText(path, resolution, chr_filter);
        if (dump_path.empty()) {
            printf("[ContactMatrixIO] Failed to dump binary mcool/cool file.\n");
            return false;
        }
        Format dump_fmt = detectFormat(dump_path);
        bool ok = readByFormat(dump_fmt, dump_path, loops, anchors,
                               resolution, chr_filter);
        remove(dump_path.c_str());
        return ok;
    }

    return readByFormat(fmt, path, loops, anchors, resolution, chr_filter);
}

#endif /* CONTACTMATRIXIO_H_ */
