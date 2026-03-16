/**
 * @file CifWriter.h
 * @brief mmCIF format writer for 3D genome structures.
 *
 * Produces standards-compliant mmCIF files with _atom_site loop
 * records, suitable for visualization in molecular viewers.
 */

#ifndef CIFWRITER_H_
#define CIFWRITER_H_

#include <string>
#include <vector>
#include <Cluster.h>

class HierarchicalChromosome;

/**
 * @class CifWriter
 * @brief Writes 3D genome structures in mmCIF format.
 */
class CifWriter {
public:
  /**
   * @brief Write a HierarchicalChromosome to an mmCIF file.
   * @param hc The structure to write.
   * @param filename Output file path.
   * @param dataName Data block name (default: "cudaMMC").
   * @return True on success.
   */
  static bool write(HierarchicalChromosome &hc, const std::string &filename,
                    const std::string &dataName = "cudaMMC",
                    bool useCurrentLevel = false);

  /**
   * @brief Write a flat list of positions to mmCIF.
   * @param positions Vector of 3D positions.
   * @param genomic_positions Corresponding genomic positions.
   * @param chrName Chromosome name for chain ID.
   * @param filename Output file path.
   * @return True on success.
   */
  static bool writeFlat(const std::vector<vector3> &positions,
                        const std::vector<int> &genomic_positions,
                        const std::string &chrName,
                        const std::string &filename);
};

#endif /* CIFWRITER_H_ */
