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

#endif /* GENOME_RECONSTRUCTOR_H_ */
