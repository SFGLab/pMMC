/**
 * @file ChromosomesSet.h
 * @brief An ensemble of chromosome 3D reconstructions.
 *
 * Stores multiple maps of {chromosome_name -> Chromosome} representing
 * an ensemble of independently reconstructed 3D structures.
 * Used by LooperSolver to collect results across ensemble members.
 */

#ifndef CHROMOSOMESSET_H_
#define CHROMOSOMESSET_H_

#include <Chromosome.h>
#include <string.h>
#include <vector>

/**
 * @class ChromosomesSet
 * @brief Collection of chromosome structure ensembles with optional descriptions.
 *
 * Each element is a map from chromosome name to its 3D structure (Chromosome).
 * The ensemble can be serialized to/from binary files.
 */
class ChromosomesSet {
public:
  /** @brief Default constructor. */
  ChromosomesSet();

  /** @brief Print summary of all stored ensembles to stdout. */
  void print();

  /**
   * @brief Add a chromosome map (one ensemble member) with empty description.
   * @param chr Map from chromosome name to Chromosome structure.
   */
  void add(std::map<std::string, Chromosome> chr);

  /**
   * @brief Add a chromosome map with a description label.
   * @param chr Map from chromosome name to Chromosome structure.
   * @param desc Description string for this ensemble member.
   */
  void add(std::map<std::string, Chromosome> chr, string desc);

  /**
   * @brief Write the entire ensemble to a binary file.
   * @param filename Output file path.
   */
  void toFile(string filename);

  /**
   * @brief Write the ensemble to an already-open file handle.
   * @param file Open FILE pointer for writing.
   */
  void toFile(FILE *file);

  /**
   * @brief Read the ensemble from a binary file.
   * @param filename Input file path.
   */
  void fromFile(string filename);

  /**
   * @brief Read the ensemble from an already-open file handle.
   * @param file Open FILE pointer for reading.
   */
  void fromFile(FILE *file);

  std::vector<std::map<std::string, Chromosome>> chromosome;  /**< Ensemble: each element maps chr name to its 3D structure. */
  std::vector<string> desc;  /**< Optional description for each ensemble member. */
};

#endif /* CHROMOSOMESSET_H_ */
