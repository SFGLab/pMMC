/**
 * @file DistanceMapGenerator.h
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

#include <Chromosome.h>
#include <Heatmap.h>

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

#endif /* DISTANCEMAPGENERATOR_H_ */
