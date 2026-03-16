/**
 * @file PdbWriter.cpp
 * @brief PDB format writer implementation.
 *
 * PDB ATOM record format (columns are 1-indexed):
 *   1-6:   Record name "ATOM  "
 *   7-11:  Atom serial (right-justified)
 *  13-16:  Atom name
 *  17:     Alternate location indicator
 *  18-20:  Residue name
 *  22:     Chain ID
 *  23-26:  Residue sequence number (right-justified)
 *  27:     Code for insertion of residues
 *  31-38:  X coordinate (8.3 format)
 *  39-46:  Y coordinate (8.3 format)
 *  47-54:  Z coordinate (8.3 format)
 *  55-60:  Occupancy (6.2 format)
 *  61-66:  Temperature factor (6.2 format)
 *  77-78:  Element symbol
 */

#include <PdbWriter.h>
#include <HierarchicalChromosome.h>
#include <cstdio>
#include <cstring>
#include <algorithm>

static const int PDB_MAX_ATOMS = 99999;
static const int PDB_MAX_RESID = 9999;

static char chrToChainChar(int chrIndex) {
  if (chrIndex < 26)
    return 'A' + (char)chrIndex;
  if (chrIndex < 52)
    return 'a' + (char)(chrIndex - 26);
  if (chrIndex < 62)
    return '0' + (char)(chrIndex - 52);
  return 'Z'; // fallback
}

bool PdbWriter::write(HierarchicalChromosome &hc, const std::string &filename,
                      bool useCurrentLevel) {
  FILE *f = fopen(filename.c_str(), "w");
  if (!f) {
    printf("[PDB] Failed to open %s for writing\n", filename.c_str());
    return false;
  }

  // Header
  fprintf(f, "HEADER    CHROMOSOME STRUCTURE                        cudaMMC\n");
  fprintf(f, "TITLE     3D GENOME STRUCTURE FROM CUDA-MMC\n");

  // Count total atoms first
  if (!useCurrentLevel)
    hc.useLowestLevel();
  int totalAtoms = 0;
  for (const auto &chr : hc.chrs) {
    if (hc.current_level.find(chr) != hc.current_level.end())
      totalAtoms += (int)hc.current_level[chr].size();
  }

  if (totalAtoms > PDB_MAX_ATOMS) {
    printf("[PDB] WARNING: %d atoms exceeds PDB limit of %d. "
           "Consider using mmCIF format (-F cif) for large models.\n",
           totalAtoms, PDB_MAX_ATOMS);
  }

  int atomSerial = 1;

  for (size_t ci = 0; ci < hc.chrs.size(); ++ci) {
    std::string chr = hc.chrs[ci];
    char chainId = chrToChainChar((int)ci);

    if (hc.current_level.find(chr) == hc.current_level.end())
      continue;

    const std::vector<int> &level = hc.current_level[chr];
    int resSeq = 1;

    for (size_t j = 0; j < level.size(); ++j) {
      int idx = level[j];
      if (idx < 0 || idx >= (int)hc.clusters.size())
        continue;

      const Cluster &cl = hc.clusters[idx];

      int serial = atomSerial;
      if (serial > PDB_MAX_ATOMS)
        serial = PDB_MAX_ATOMS; // clamp for formatting

      int res = resSeq;
      if (res > PDB_MAX_RESID)
        res = resSeq % PDB_MAX_RESID + 1;

      // PDB ATOM format: strict column-based
      fprintf(f,
              "ATOM  %5d  CA  ALA %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          C \n",
              serial, chainId, res, cl.pos.x, cl.pos.y, cl.pos.z, 1.00f,
              std::min((float)cl.genomic_pos / 1000.0f, 99.99f));

      atomSerial++;
      resSeq++;
    }

    // TER record between chains
    fprintf(f, "TER   %5d      ALA %c%4d\n",
            std::min(atomSerial, PDB_MAX_ATOMS), chainId,
            std::min(resSeq - 1, PDB_MAX_RESID));
    atomSerial++;
  }

  fprintf(f, "END\n");
  fclose(f);

  printf("[PDB] Wrote %d atoms to %s\n", atomSerial - 1 - (int)hc.chrs.size(),
         filename.c_str());
  return true;
}

bool PdbWriter::writeFlat(const std::vector<vector3> &positions,
                          const std::string &chrName,
                          const std::string &filename) {
  FILE *f = fopen(filename.c_str(), "w");
  if (!f) {
    printf("[PDB] Failed to open %s for writing\n", filename.c_str());
    return false;
  }

  fprintf(f, "HEADER    CHROMOSOME STRUCTURE                        cudaMMC\n");
  fprintf(f, "TITLE     3D GENOME STRUCTURE FROM CUDA-MMC\n");

  if ((int)positions.size() > PDB_MAX_ATOMS) {
    printf("[PDB] WARNING: %d atoms exceeds PDB limit of %d. "
           "Consider using mmCIF format (-F cif) for large models.\n",
           (int)positions.size(), PDB_MAX_ATOMS);
  }

  for (size_t i = 0; i < positions.size(); ++i) {
    int serial = std::min((int)(i + 1), PDB_MAX_ATOMS);
    int resSeq = (int)(i + 1);
    if (resSeq > PDB_MAX_RESID)
      resSeq = resSeq % PDB_MAX_RESID + 1;

    fprintf(f,
            "ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          C \n",
            serial, resSeq, positions[i].x, positions[i].y, positions[i].z,
            1.00f, 0.00f);
  }

  fprintf(f, "TER   %5d      ALA A%4d\n",
          std::min((int)(positions.size() + 1), PDB_MAX_ATOMS),
          std::min((int)positions.size(), PDB_MAX_RESID));
  fprintf(f, "END\n");
  fclose(f);

  printf("[PDB] Wrote %d atoms to %s\n", (int)positions.size(),
         filename.c_str());
  return true;
}
