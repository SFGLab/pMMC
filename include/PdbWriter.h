/**
 * @file PdbWriter.h
 * @brief PDB format writer for 3D genome structures.
 *
 * Produces PDB files with ATOM records. Limited to 99,999 atoms
 * (PDB format constraint); warns if exceeded.
 */

#ifndef PDBWRITER_H_
#define PDBWRITER_H_

#include <string>
#include <vector>
#include <Cluster.h>

class HierarchicalChromosome;

/**
 * @class PdbWriter
 * @brief Writes 3D genome structures in PDB format.
 */
class PdbWriter {
public:
  /**
   * @brief Write a HierarchicalChromosome to a PDB file.
   * @param hc The structure to write.
   * @param filename Output file path.
   * @return True on success.
   */
  static bool write(HierarchicalChromosome &hc, const std::string &filename);

  /**
   * @brief Write a flat list of positions to PDB.
   * @param positions Vector of 3D positions.
   * @param chrName Chromosome name for chain ID.
   * @param filename Output file path.
   * @return True on success.
   */
  static bool writeFlat(const std::vector<vector3> &positions,
                        const std::string &chrName,
                        const std::string &filename);
};

#endif /* PDBWRITER_H_ */
