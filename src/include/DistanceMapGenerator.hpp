/**
 * @file DistanceMapGenerator.hpp
 * @brief Standalone module for pairwise distance and contact probability maps.
 *
 * Generates:
 *   - Pairwise 3D distance matrices from polymer structures
 *   - Contact probability matrices with configurable threshold
 *   - Ensemble-averaged distance and contact maps
 *   - Contact frequency maps from inverse distance power law
 *
 * @see Heatmap for the underlying matrix storage.
 * @see HierarchicalChromosome for multi-level structure extraction.
 */

#ifndef DISTANCEMAPGENERATOR_H_
#define DISTANCEMAPGENERATOR_H_

#include <string>
#include <vector>

#include <Chromosome.hpp>
#include <Heatmap.hpp>

/**
 * @class DistanceMapGenerator
 * @brief Computes distance and contact maps from 3D polymer structures.
 *
 * Usage:
 * @code
 *   DistanceMapGenerator gen;
 *   gen.setContactThreshold(2.0f);
 *   gen.setContactExponent(2.0f);
 *   Heatmap dist = gen.computeDistanceMap(chromosome);
 *   Heatmap contact = gen.computeContactMap(chromosome);
 *
 *   // Ensemble averaging
 *   gen.addEnsembleMember(chr1);
 *   gen.addEnsembleMember(chr2);
 *   Heatmap avg_dist = gen.computeEnsembleAverageDistance();
 *   Heatmap avg_contact = gen.computeEnsembleAverageContact();
 * @endcode
 */
class DistanceMapGenerator {
public:
  DistanceMapGenerator();

  /**
   * @brief Set the contact distance threshold.
   * @param threshold Pairs closer than this are "in contact".
   */
  void setContactThreshold(float threshold);

  /**
   * @brief Set the power-law exponent for contact frequency.
   * @param exponent freq ~ 1 / dist^exponent.
   */
  void setContactExponent(float exponent);

  /**
   * @brief Compute pairwise 3D distance matrix.
   * @param chr Polymer structure.
   * @return NxN Heatmap of pairwise distances.
   */
  Heatmap computeDistanceMap(const Chromosome &chr) const;

  /**
   * @brief Compute binary contact map (1 if dist < threshold, else 0).
   * @param chr Polymer structure.
   * @return NxN Heatmap of contact indicators.
   */
  Heatmap computeContactMap(const Chromosome &chr) const;

  /**
   * @brief Compute contact frequency map using power-law decay.
   * @param chr Polymer structure.
   * @return NxN Heatmap of contact frequencies (1/dist^exponent).
   */
  Heatmap computeContactFrequencyMap(const Chromosome &chr) const;

  /**
   * @brief Add a structure to the ensemble for averaging.
   * @param chr Polymer structure to add.
   */
  void addEnsembleMember(const Chromosome &chr);

  /** @brief Clear all ensemble members. */
  void clearEnsemble();

  /** @brief Number of structures in the ensemble. */
  int ensembleSize() const;

  /**
   * @brief Compute ensemble-averaged pairwise distance map.
   * @return Averaged distance Heatmap.
   */
  Heatmap computeEnsembleAverageDistance() const;

  /**
   * @brief Compute ensemble-averaged contact probability map.
   * @return Heatmap where h[i][j] = fraction of ensemble members with contact.
   */
  Heatmap computeEnsembleAverageContact() const;

  /**
   * @brief Compute ensemble-averaged contact frequency map.
   * @return Averaged contact frequency Heatmap.
   */
  Heatmap computeEnsembleAverageFrequency() const;

  /**
   * @brief Write distance map to file.
   * @param chr Polymer structure.
   * @param output_path Output file path.
   */
  void writeDistanceMap(const Chromosome &chr,
                        const std::string &output_path) const;

  /**
   * @brief Write contact map to file.
   * @param chr Polymer structure.
   * @param output_path Output file path.
   */
  void writeContactMap(const Chromosome &chr,
                       const std::string &output_path) const;

private:
  float contact_threshold;
  float contact_exponent;
  std::vector<Chromosome> ensemble;
};


// ============================================================================
// Implementation
// ============================================================================


#include <algorithm>
#include <cmath>
#include <cstdio>

inline DistanceMapGenerator::DistanceMapGenerator()
    : contact_threshold(2.0f), contact_exponent(2.0f) {}

inline void DistanceMapGenerator::setContactThreshold(float threshold) {
  contact_threshold = threshold;
}

inline void DistanceMapGenerator::setContactExponent(float exponent) {
  contact_exponent = exponent;
}

inline Heatmap DistanceMapGenerator::computeDistanceMap(const Chromosome &chr) const {
  int n = chr.size;
  Heatmap h(n);

  for (int i = 0; i < n; ++i) {
    h.v[i][i] = 0.0f;
    for (int j = i + 1; j < n; ++j) {
      float dist = (chr.points[i] - chr.points[j]).length();
      h.v[i][j] = dist;
      h.v[j][i] = dist;
    }
  }

  return h;
}

inline Heatmap DistanceMapGenerator::computeContactMap(const Chromosome &chr) const {
  int n = chr.size;
  Heatmap h(n);

  float thresh2 = contact_threshold * contact_threshold;

  for (int i = 0; i < n; ++i) {
    h.v[i][i] = 1.0f;
    for (int j = i + 1; j < n; ++j) {
      vector3 diff = chr.points[i] - chr.points[j];
      float d2 = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
      float val = (d2 < thresh2) ? 1.0f : 0.0f;
      h.v[i][j] = val;
      h.v[j][i] = val;
    }
  }

  return h;
}

Heatmap
inline DistanceMapGenerator::computeContactFrequencyMap(const Chromosome &chr) const {
  int n = chr.size;
  Heatmap h(n);

  for (int i = 0; i < n; ++i) {
    h.v[i][i] = 0.0f;
    for (int j = i + 1; j < n; ++j) {
      float dist = (chr.points[i] - chr.points[j]).length();
      float freq = 0.0f;
      if (dist > 1e-6f)
        freq = 1.0f / std::pow(dist, contact_exponent);
      h.v[i][j] = freq;
      h.v[j][i] = freq;
    }
  }

  return h;
}

inline void DistanceMapGenerator::addEnsembleMember(const Chromosome &chr) {
  ensemble.push_back(chr);
}

inline void DistanceMapGenerator::clearEnsemble() { ensemble.clear(); }

inline int DistanceMapGenerator::ensembleSize() const {
  return static_cast<int>(ensemble.size());
}

inline Heatmap DistanceMapGenerator::computeEnsembleAverageDistance() const {
  if (ensemble.empty())
    return Heatmap();

  int n = ensemble[0].size;
  Heatmap avg(n);
  int count = static_cast<int>(ensemble.size());

  for (const auto &chr : ensemble) {
    int m = std::min(n, chr.size);
    for (int i = 0; i < m; ++i) {
      for (int j = i + 1; j < m; ++j) {
        float dist = (chr.points[i] - chr.points[j]).length();
        avg.v[i][j] += dist;
        avg.v[j][i] += dist;
      }
    }
  }

  float inv_count = 1.0f / count;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      avg.v[i][j] *= inv_count;

  return avg;
}

inline Heatmap DistanceMapGenerator::computeEnsembleAverageContact() const {
  if (ensemble.empty())
    return Heatmap();

  int n = ensemble[0].size;
  Heatmap avg(n);
  int count = static_cast<int>(ensemble.size());
  float thresh2 = contact_threshold * contact_threshold;

  for (const auto &chr : ensemble) {
    int m = std::min(n, chr.size);
    for (int i = 0; i < m; ++i) {
      for (int j = i + 1; j < m; ++j) {
        vector3 diff = chr.points[i] - chr.points[j];
        float d2 = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
        if (d2 < thresh2) {
          avg.v[i][j] += 1.0f;
          avg.v[j][i] += 1.0f;
        }
      }
    }
  }

  float inv_count = 1.0f / count;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      avg.v[i][j] *= inv_count;

  return avg;
}

inline Heatmap DistanceMapGenerator::computeEnsembleAverageFrequency() const {
  if (ensemble.empty())
    return Heatmap();

  int n = ensemble[0].size;
  Heatmap avg(n);
  int count = static_cast<int>(ensemble.size());

  for (const auto &chr : ensemble) {
    int m = std::min(n, chr.size);
    for (int i = 0; i < m; ++i) {
      for (int j = i + 1; j < m; ++j) {
        float dist = (chr.points[i] - chr.points[j]).length();
        if (dist > 1e-6f) {
          float freq = 1.0f / std::pow(dist, contact_exponent);
          avg.v[i][j] += freq;
          avg.v[j][i] += freq;
        }
      }
    }
  }

  float inv_count = 1.0f / count;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      avg.v[i][j] *= inv_count;

  return avg;
}

inline void DistanceMapGenerator::writeDistanceMap(
    const Chromosome &chr, const std::string &output_path) const {
  Heatmap h = computeDistanceMap(chr);
  h.toFile(output_path, false);
  printf("Distance map (%dx%d) written to %s\n", (int)h.size, (int)h.size,
         output_path.c_str());
}

inline void DistanceMapGenerator::writeContactMap(
    const Chromosome &chr, const std::string &output_path) const {
  Heatmap h = computeContactMap(chr);
  h.toFile(output_path, false);
  printf("Contact map (%dx%d) written to %s\n", (int)h.size, (int)h.size,
         output_path.c_str());
}

#endif /* DISTANCEMAPGENERATOR_H_ */
