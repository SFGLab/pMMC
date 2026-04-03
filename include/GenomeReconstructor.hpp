#ifndef GENOME_RECONSTRUCTOR_H_
#define GENOME_RECONSTRUCTOR_H_

#include <string>
#include <vector>
#include <map>
#include <HierarchicalChromosome.h>

struct ChromosomeTerritory {
    std::string chr_name;
    double center_x, center_y, center_z;
    double radius_of_gyration;
    int n_beads;
    double bounding_radius;
};

struct TerritoryPair {
    std::string chr_a, chr_b;
    double min_distance;
    int overlap_count;  // beads within threshold
};

class GenomeReconstructor {
public:
    GenomeReconstructor();

    void setOverlapThreshold(double threshold);

    // Analyze a reconstructed model
    void analyze(HierarchicalChromosome& hc, int level = 1);

    // Write territory report
    void writeReport(const std::string& outdir, const std::string& label);

    // Getters
    const std::vector<ChromosomeTerritory>& getTerritories() const;
    const std::vector<TerritoryPair>& getPairs() const;

private:
    std::vector<ChromosomeTerritory> territories_;
    std::vector<TerritoryPair> pairs_;
    double overlap_threshold_;
};


// ============================================================================
// Implementation
// ============================================================================

#include <SvgChartGenerator.hpp>
#include <common.h>
#include <math.h>
#include <stdio.h>
#include <float.h>

inline GenomeReconstructor::GenomeReconstructor() : overlap_threshold_(2.0) {}

inline void GenomeReconstructor::setOverlapThreshold(double threshold) {
    overlap_threshold_ = threshold;
}

inline void GenomeReconstructor::analyze(HierarchicalChromosome& hc, int level) {
    territories_.clear();
    pairs_.clear();

    hc.setLevel(level);
    hc.createCurrentLevelStructure();

    if (hc.chrs.size() < 1) {
        printf("[GenomeReconstructor] No chromosomes found in model\n");
        return;
    }

    // Collect per-chromosome bead positions
    std::map<std::string, std::vector<vector3>> chr_positions;

    for (const std::string& chr_name : hc.chrs) {
        if (hc.current_level.find(chr_name) == hc.current_level.end())
            continue;

        const std::vector<int>& indices = hc.current_level[chr_name];
        std::vector<vector3> positions;
        for (int idx : indices) {
            if (idx >= 0 && idx < (int)hc.clusters.size()) {
                positions.push_back(hc.clusters[idx].pos);
            }
        }
        if (!positions.empty())
            chr_positions[chr_name] = positions;
    }

    // Compute territory metrics for each chromosome
    for (auto& kv : chr_positions) {
        const std::string& chr_name = kv.first;
        const std::vector<vector3>& positions = kv.second;

        ChromosomeTerritory t;
        t.chr_name = chr_name;
        t.n_beads = (int)positions.size();

        // Center of mass
        double cx = 0, cy = 0, cz = 0;
        for (const vector3& p : positions) {
            cx += p.x;
            cy += p.y;
            cz += p.z;
        }
        t.center_x = cx / t.n_beads;
        t.center_y = cy / t.n_beads;
        t.center_z = cz / t.n_beads;

        // Radius of gyration
        double rg_sum = 0;
        double max_r = 0;
        for (const vector3& p : positions) {
            double dx = p.x - t.center_x;
            double dy = p.y - t.center_y;
            double dz = p.z - t.center_z;
            double r2 = dx * dx + dy * dy + dz * dz;
            rg_sum += r2;
            double r = sqrt(r2);
            if (r > max_r)
                max_r = r;
        }
        t.radius_of_gyration = sqrt(rg_sum / t.n_beads);
        t.bounding_radius = max_r;

        territories_.push_back(t);
    }

    // Compute pairwise territory metrics
    for (size_t i = 0; i < territories_.size(); i++) {
        for (size_t j = i + 1; j < territories_.size(); j++) {
            const std::string& chr_a = territories_[i].chr_name;
            const std::string& chr_b = territories_[j].chr_name;
            const std::vector<vector3>& pos_a = chr_positions[chr_a];
            const std::vector<vector3>& pos_b = chr_positions[chr_b];

            TerritoryPair pair;
            pair.chr_a = chr_a;
            pair.chr_b = chr_b;
            pair.min_distance = DBL_MAX;
            pair.overlap_count = 0;

            for (const vector3& pa : pos_a) {
                for (const vector3& pb : pos_b) {
                    double dx = pa.x - pb.x;
                    double dy = pa.y - pb.y;
                    double dz = pa.z - pb.z;
                    double dist = sqrt(dx * dx + dy * dy + dz * dz);

                    if (dist < pair.min_distance)
                        pair.min_distance = dist;
                    if (dist < overlap_threshold_)
                        pair.overlap_count++;
                }
            }

            pairs_.push_back(pair);
        }
    }

    printf("[GenomeReconstructor] Analyzed %d chromosomes, %d pairs\n",
           (int)territories_.size(), (int)pairs_.size());
}

inline void GenomeReconstructor::writeReport(const std::string& outdir, const std::string& label) {
    std::string path = outdir + "territory_report_" + label + ".txt";
    FILE* f = fopen(path.c_str(), "w");
    if (!f) {
        printf("[GenomeReconstructor] Could not open %s for writing\n", path.c_str());
        return;
    }

    fprintf(f, "# Chromosomal Territory Report\n");
    fprintf(f, "# Label: %s\n\n", label.c_str());

    // Per-chromosome section
    fprintf(f, "# Per-chromosome territories\n");
    fprintf(f, "chr\tcenter_x\tcenter_y\tcenter_z\tradius_of_gyration\tn_beads\tbounding_radius\n");
    for (const ChromosomeTerritory& t : territories_) {
        fprintf(f, "%s\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%.4f\n",
                t.chr_name.c_str(), t.center_x, t.center_y, t.center_z,
                t.radius_of_gyration, t.n_beads, t.bounding_radius);
    }

    fprintf(f, "\n# Pairwise territory metrics (overlap_threshold=%.2f)\n", overlap_threshold_);
    fprintf(f, "chr_a\tchr_b\tmin_distance\toverlap_count\n");
    for (const TerritoryPair& p : pairs_) {
        fprintf(f, "%s\t%s\t%.4f\t%d\n",
                p.chr_a.c_str(), p.chr_b.c_str(), p.min_distance, p.overlap_count);
    }

    fclose(f);
    printf("[GenomeReconstructor] Territory report written to %s\n", path.c_str());

    // Generate Rg bar chart SVG if we have territory data
    if (!territories_.empty()) {
        std::vector<std::string> chr_names;
        std::vector<float> rg_values;
        for (const ChromosomeTerritory& t : territories_) {
            chr_names.push_back(t.chr_name);
            rg_values.push_back((float)t.radius_of_gyration);
        }
        std::string rg_svg = outdir + "rg_plot_" + label + ".svg";
        if (SvgChartGenerator::rgBarChart(chr_names, rg_values, rg_svg)) {
            printf("[SVG] Chart written: %s\n", rg_svg.c_str());
        }
    }
}

inline const std::vector<ChromosomeTerritory>& GenomeReconstructor::getTerritories() const {
    return territories_;
}

inline const std::vector<TerritoryPair>& GenomeReconstructor::getPairs() const {
    return pairs_;
}

#endif /* GENOME_RECONSTRUCTOR_H_ */
