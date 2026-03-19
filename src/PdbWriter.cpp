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
#include <cmath>
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

  // Genomic PDB marker — enables 3DGenomeViewer to detect this as a genomic
  // structure and render it identically to HCM files.
  fprintf(f, "REMARK 999 GENOMIC_PDB CONVERTED_FROM_HCM\n");

  // Factor records (one per ChIA-PET factor)
  for (size_t fi = 0; fi < hc.arcs.factors.size(); ++fi) {
    fprintf(f, "REMARK 999 FACTOR %d %s\n", (int)fi,
            hc.arcs.factors[fi].c_str());
  }

  // Chromosome records
  for (const auto &chr : hc.chrs) {
    fprintf(f, "REMARK 999 CHROMOSOME %s\n", chr.c_str());
  }

  // Use lowest level for maximum detail
  if (!useCurrentLevel)
    hc.useLowestLevel();

  // Collect all beads across all chromosomes in order, using a global
  // sequential resSeq (matching 3DGenomeViewer's HCM→PDB converter format).
  struct BeadInfo {
    int clusterIdx;
    int chrIdx;
    std::string chrName;
    int factorIdx;
  };
  std::vector<BeadInfo> beads;

  for (size_t ci = 0; ci < hc.chrs.size(); ++ci) {
    std::string chr = hc.chrs[ci];
    if (hc.current_level.find(chr) == hc.current_level.end())
      continue;
    const std::vector<int> &level = hc.current_level[chr];
    for (size_t j = 0; j < level.size(); ++j) {
      int idx = level[j];
      if (idx < 0 || idx >= (int)hc.clusters.size())
        continue;
      const Cluster &cl = hc.clusters[idx];

      // Determine factor index from arcs
      int factorIdx = -1;
      if (hc.arcs.arcs.find(chr) != hc.arcs.arcs.end()) {
        const auto &chrArcVec = hc.arcs.arcs[chr];
        for (int ai : cl.arcs) {
          if (ai >= 0 && ai < (int)chrArcVec.size() &&
              chrArcVec[ai].factor >= 0) {
            factorIdx = chrArcVec[ai].factor;
            break;
          }
        }
      }
      beads.push_back({idx, (int)ci, chr, factorIdx});
    }
  }

  if ((int)beads.size() > PDB_MAX_ATOMS) {
    printf("[PDB] WARNING: %d atoms exceeds PDB limit of %d. "
           "Consider using mmCIF format (-F cif) for large models.\n",
           (int)beads.size(), PDB_MAX_ATOMS);
  }

  // Write per-residue REMARK 999 RESIDUE records with lossless genomic
  // positions, hierarchy levels, and factor indices.
  // Format matches what 3DGenomeViewer's CoordParser expects:
  //   REMARK 999 RESIDUE <resSeq> GENOMICPOS <bp> LEVEL <lvl> FACTOR <fIdx>
  for (size_t i = 0; i < beads.size(); ++i) {
    int resSeq = (int)(i + 1);
    if (resSeq > PDB_MAX_RESID)
      resSeq = resSeq % PDB_MAX_RESID + 1;
    const Cluster &cl = hc.clusters[beads[i].clusterIdx];
    fprintf(f, "REMARK 999 RESIDUE %d GENOMICPOS %d LEVEL %d FACTOR %d\n",
            resSeq, cl.genomic_pos, cl.level, beads[i].factorIdx);
  }

  // Write HELIX records spanning each chain (so Chimera renders wide ribbon).
  // We track which global resSeq range belongs to each chain.
  {
    int helixSerial = 1;
    int globalIdx = 0;
    for (size_t ci = 0; ci < hc.chrs.size(); ++ci) {
      std::string chr = hc.chrs[ci];
      char chainId = chrToChainChar((int)ci);
      if (hc.current_level.find(chr) == hc.current_level.end())
        continue;
      int chainLen = (int)hc.current_level[chr].size();
      if (chainLen < 1) {
        globalIdx += chainLen;
        continue;
      }

      int firstResSeq = globalIdx + 1;
      int lastResSeq = globalIdx + chainLen;
      // Apply same wrapping as ATOM records
      int firstRes = firstResSeq;
      if (firstRes > PDB_MAX_RESID)
        firstRes = firstResSeq % PDB_MAX_RESID + 1;
      int lastRes = lastResSeq;
      if (lastRes > PDB_MAX_RESID)
        lastRes = lastResSeq % PDB_MAX_RESID + 1;

      fprintf(f,
              "HELIX  %3d %3d GLY %c %4d  GLY %c %4d  1\n",
              helixSerial, helixSerial, chainId, firstRes, chainId, lastRes);
      globalIdx += chainLen;
      helixSerial++;
    }
  }

  // ── Normalization: center + scale to ±50 range ─────────────────────
  // This ensures:
  //   1. PDB coordinates fit in the %8.3f column format (max ±9999.999)
  //   2. The structure renders at a visible scale in Chimera
  //   3. Pixel-identical output across HCM, PDB, and CIF rendering paths
  float cx = 0, cy = 0, cz = 0;
  int validCount = 0;
  for (size_t i = 0; i < beads.size(); ++i) {
    const Cluster &cl = hc.clusters[beads[i].clusterIdx];
    cx += cl.pos.x; cy += cl.pos.y; cz += cl.pos.z;
    validCount++;
  }
  if (validCount > 0) { cx /= validCount; cy /= validCount; cz /= validCount; }

  float maxDist = 1.0f;
  for (size_t i = 0; i < beads.size(); ++i) {
    const Cluster &cl = hc.clusters[beads[i].clusterIdx];
    float dx = cl.pos.x - cx, dy = cl.pos.y - cy, dz = cl.pos.z - cz;
    float d = sqrtf(dx*dx + dy*dy + dz*dz);
    if (d > maxDist) maxDist = d;
  }
  float scale = 50.0f / maxDist;

  // Write ATOM and TER records
  int atomSerial = 1;
  int prevChrIdx = -1;

  for (size_t i = 0; i < beads.size(); ++i) {
    char chainId = chrToChainChar(beads[i].chrIdx);

    // Write TER before switching to a new chain (except first chain)
    if (beads[i].chrIdx != prevChrIdx && prevChrIdx >= 0) {
      char prevChainId = chrToChainChar(prevChrIdx);
      int prevResSeq = (int)i; // previous bead's resSeq
      if (prevResSeq > PDB_MAX_RESID)
        prevResSeq = prevResSeq % PDB_MAX_RESID + 1;
      fprintf(f, "TER   %5d      GLY %c%4d\n",
              std::min(atomSerial, PDB_MAX_ATOMS), prevChainId, prevResSeq);
      atomSerial++;
    }
    prevChrIdx = beads[i].chrIdx;

    const Cluster &cl = hc.clusters[beads[i].clusterIdx];
    int serial = std::min(atomSerial, PDB_MAX_ATOMS);
    int resSeq = (int)(i + 1);
    if (resSeq > PDB_MAX_RESID)
      resSeq = resSeq % PDB_MAX_RESID + 1;

    // Normalized coordinates (centroid + scale to ±50 range)
    float nx = (cl.pos.x - cx) * scale;
    float ny = (cl.pos.y - cy) * scale;
    float nz = (cl.pos.z - cz) * scale;

    // B-factor carries genomic_pos / 1e6 (Mb scale, fits up to 999.99 Mb)
    // Occupancy encodes hierarchy level (as float)
    fprintf(f,
            "ATOM  %5d  CA  GLY %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          C \n",
            serial, chainId, resSeq, nx, ny, nz,
            (float)cl.level, std::min((float)cl.genomic_pos / 1e6f, 999.99f));

    atomSerial++;
  }

  // Final TER
  if (!beads.empty()) {
    char lastChainId = chrToChainChar(beads.back().chrIdx);
    int lastResSeq = (int)beads.size();
    if (lastResSeq > PDB_MAX_RESID)
      lastResSeq = lastResSeq % PDB_MAX_RESID + 1;
    fprintf(f, "TER   %5d      GLY %c%4d\n",
            std::min(atomSerial, PDB_MAX_ATOMS), lastChainId, lastResSeq);
  }

  fprintf(f, "END\n");
  fclose(f);

  printf("[PDB] Wrote %d atoms to %s\n", (int)beads.size(), filename.c_str());
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
