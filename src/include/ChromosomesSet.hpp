/**
 * @file ChromosomesSet.hpp
 * @brief An ensemble of chromosome 3D reconstructions.
 *
 * Stores multiple maps of {chromosome_name -> Chromosome} representing
 * an ensemble of independently reconstructed 3D structures.
 * Used by LooperSolver to collect results across ensemble members.
 */

#ifndef CHROMOSOMESSET_H_
#define CHROMOSOMESSET_H_

#include <Chromosome.hpp>
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


// ============================================================================
// Implementation
// ============================================================================


inline ChromosomesSet::ChromosomesSet() {
  // TODO Auto-generated constructor stub
}

inline void ChromosomesSet::print() {
  printf("Set size = %zu\n", chromosome.size());

  for (size_t i = 0; i < chromosome.size(); ++i) {
    //	printf("%d %s\n", chromosome[i].points.size(), desc[i].c_str());
  }
}

inline void ChromosomesSet::add(std::map<std::string, Chromosome> chr) {
  add(chr, "<no_desc>");
}

inline void ChromosomesSet::add(std::map<std::string, Chromosome> chr, string desc) {
  chromosome.push_back(chr);
  if (desc == "")
    desc = "none";
  this->desc.push_back(desc);
}

inline void ChromosomesSet::toFile(string filename) {
  FILE *f;
  f = fopen(filename.c_str(), "w");
  if (f == NULL) {
    printf("Error opening file [%s]!\n", filename.c_str());
    return;
  }

  toFile(f);

  fclose(f);
}

inline void ChromosomesSet::toFile(FILE *file) {
  // printf("set = %d\n", chromosome.size());
  fprintf(file, "%zu\n", chromosome.size());
  for (size_t i = 0; i < chromosome.size(); i++) {
    // fprintf(file, "%d\n%s\n", chromosome[i].points.size(), desc[i].c_str());
    fprintf(file, "%zu %s\n", chromosome[i].size(), desc[i].c_str());
    for (auto el : chromosome[i]) {
      fprintf(file, "%s %zu\n", el.first.c_str(), el.second.points.size());
      el.second.toFile(file);
      // chromosome[i].toFile(file);
    }
  }
}

inline void ChromosomesSet::fromFile(string filename) {
  FILE *f = open(filename, "r");
  if (f == NULL)
    return;
  fromFile(f);
  fclose(f);
}

inline void ChromosomesSet::fromFile(FILE *file) {
  chromosome.clear();

  char de[100], str_chr[10];
  int n, n_chr;
  fscanf(file, "%d", &n);

  int pts;
  for (int i = 0; i < n; i++) {
    fscanf(file, "%d %100s", &n_chr, de);
    string s(de);
    desc.push_back(s);

    std::map<std::string, Chromosome> map_chr;
    for (int j = 0; j < n_chr; ++j) {
      fscanf(file, "%s %d", str_chr, &pts);
      Chromosome chr;
      chr.fromFile(file, pts);
      // chromosome.push_back(chr);
      map_chr[str_chr] = chr;
    }

    chromosome.push_back(map_chr);
  }
}

#endif /* CHROMOSOMESSET_H_ */
