/**
 * @file CifWriter.cpp
 * @brief mmCIF format writer implementation.
 */

#include <CifWriter.h>
#include <HierarchicalChromosome.h>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>

static void writeHeader(FILE *f, const std::string &dataName) {
  fprintf(f, "data_%s\n", dataName.c_str());
  fprintf(f, "#\n");
  fprintf(f, "_entry.id %s\n", dataName.c_str());
  fprintf(f, "#\n");
  fprintf(f, "_audit_conform.dict_name       mmcif_pdbx.dic\n");
  fprintf(f, "_audit_conform.dict_version    5.296\n");
  fprintf(f, "_audit_conform.dict_location   "
             "http://mmcif.pdb.org/dictionaries/ascii/mmcif_pdbx.dic\n");
  fprintf(f, "#\n");
}

static void writeGenomicMetadata(FILE *f,
                                 const std::vector<std::string> &factors,
                                 const std::vector<std::string> &chrs) {
  // Genomic CIF marker — enables 3DGenomeViewer to detect this as a genomic
  // structure and render it identically to HCM files.
  fprintf(f, "# GENOMIC_CIF\n");
  for (size_t i = 0; i < factors.size(); ++i) {
    fprintf(f, "# FACTOR %d %s\n", (int)i, factors[i].c_str());
  }
  for (const auto &chr : chrs) {
    fprintf(f, "# CHROMOSOME %s\n", chr.c_str());
  }
  fprintf(f, "#\n");
}

static void writeAtomSiteHeader(FILE *f) {
  fprintf(f, "loop_\n");
  fprintf(f, "_atom_site.group_PDB\n");
  fprintf(f, "_atom_site.id\n");
  fprintf(f, "_atom_site.type_symbol\n");
  fprintf(f, "_atom_site.label_atom_id\n");
  fprintf(f, "_atom_site.label_alt_id\n");
  fprintf(f, "_atom_site.label_comp_id\n");
  fprintf(f, "_atom_site.label_asym_id\n");
  fprintf(f, "_atom_site.label_entity_id\n");
  fprintf(f, "_atom_site.label_seq_id\n");
  fprintf(f, "_atom_site.pdbx_PDB_ins_code\n");
  fprintf(f, "_atom_site.Cartn_x\n");
  fprintf(f, "_atom_site.Cartn_y\n");
  fprintf(f, "_atom_site.Cartn_z\n");
  fprintf(f, "_atom_site.occupancy\n");
  fprintf(f, "_atom_site.B_iso_or_equiv\n");
  fprintf(f, "_atom_site.pdbx_auth_seq_id\n");
  fprintf(f, "_atom_site.auth_asym_id\n");
}

// Map chromosome name to a short chain ID (A-Z, then AA, AB, ...)
static std::string chrToChainId(const std::string &chr, int chrIndex) {
  if (chrIndex < 26) {
    char c = 'A' + (char)chrIndex;
    return std::string(1, c);
  }
  // For >26 chromosomes, use two-letter codes
  char c1 = 'A' + (char)(chrIndex / 26 - 1);
  char c2 = 'A' + (char)(chrIndex % 26);
  return std::string(1, c1) + std::string(1, c2);
}

bool CifWriter::write(HierarchicalChromosome &hc, const std::string &filename,
                      const std::string &dataName, bool useCurrentLevel) {
  FILE *f = fopen(filename.c_str(), "w");
  if (!f) {
    printf("[CIF] Failed to open %s for writing\n", filename.c_str());
    return false;
  }

  writeHeader(f, dataName);
  writeGenomicMetadata(f, hc.arcs.factors, hc.chrs);

  // Use lowest level for maximum detail, unless caller wants current level
  if (!useCurrentLevel)
    hc.useLowestLevel();

  // ── Collect all beads across chromosomes ──────────────────────────
  struct BeadRef { int clusterIdx; int chrIdx; std::string chainId; };
  std::vector<BeadRef> allBeads;
  // Track per-chain bead counts for HELIX records
  struct ChainInfo { std::string chainId; int startSeq; int endSeq; };
  std::vector<ChainInfo> chains;

  for (size_t ci = 0; ci < hc.chrs.size(); ++ci) {
    std::string chr = hc.chrs[ci];
    std::string cid = chrToChainId(chr, (int)ci);
    if (hc.current_level.find(chr) == hc.current_level.end()) continue;
    const std::vector<int> &level = hc.current_level[chr];
    int chainStart = (int)allBeads.size() + 1;
    for (size_t j = 0; j < level.size(); ++j) {
      int idx = level[j];
      if (idx >= 0 && idx < (int)hc.clusters.size())
        allBeads.push_back({idx, (int)ci, cid});
    }
    int chainEnd = (int)allBeads.size();
    if (chainEnd >= chainStart)
      chains.push_back({cid, 1, chainEnd - chainStart + 1});
  }

  // ── Write _struct_conf (HELIX records) matching PDB HELIX ─────────
  // Full mmCIF _struct_conf format matching real PDB archive CIF files.
  // Chimera v1.19 reads beg/end_label_asym_id + beg/end_label_seq_id
  // to assign isHelix to residues.
  fprintf(f, "loop_\n");
  fprintf(f, "_struct_conf.conf_type_id\n");
  fprintf(f, "_struct_conf.id\n");
  fprintf(f, "_struct_conf.pdbx_PDB_helix_id\n");
  fprintf(f, "_struct_conf.beg_label_comp_id\n");
  fprintf(f, "_struct_conf.beg_label_asym_id\n");
  fprintf(f, "_struct_conf.beg_label_seq_id\n");
  fprintf(f, "_struct_conf.pdbx_beg_PDB_ins_code\n");
  fprintf(f, "_struct_conf.end_label_comp_id\n");
  fprintf(f, "_struct_conf.end_label_asym_id\n");
  fprintf(f, "_struct_conf.end_label_seq_id\n");
  fprintf(f, "_struct_conf.pdbx_end_PDB_ins_code\n");
  fprintf(f, "_struct_conf.beg_auth_comp_id\n");
  fprintf(f, "_struct_conf.beg_auth_asym_id\n");
  fprintf(f, "_struct_conf.beg_auth_seq_id\n");
  fprintf(f, "_struct_conf.end_auth_comp_id\n");
  fprintf(f, "_struct_conf.end_auth_asym_id\n");
  fprintf(f, "_struct_conf.end_auth_seq_id\n");
  fprintf(f, "_struct_conf.pdbx_PDB_helix_class\n");
  fprintf(f, "_struct_conf.details\n");
  fprintf(f, "_struct_conf.pdbx_PDB_helix_length\n");
  for (size_t hi = 0; hi < chains.size(); ++hi) {
    int len = chains[hi].endSeq - chains[hi].startSeq + 1;
    fprintf(f, "HELX_P HELX_P%d H%d GLY %s %d ? GLY %s %d ? GLY %s %d GLY %s %d 1 ? %d\n",
            (int)(hi + 1), (int)(hi + 1),
            chains[hi].chainId.c_str(), chains[hi].startSeq,
            chains[hi].chainId.c_str(), chains[hi].endSeq,
            chains[hi].chainId.c_str(), chains[hi].startSeq,
            chains[hi].chainId.c_str(), chains[hi].endSeq,
            len);
  }
  fprintf(f, "#\n");

  // ── Write _entity_poly_seq (polymer sequence) ─────────────────────
  // Chimera's mmCIF reader uses this table to establish polymer
  // connectivity (C→N backbone bonds) and build the residueAfter chain
  // required for _struct_conf HELIX assignment.
  fprintf(f, "loop_\n");
  fprintf(f, "_entity_poly_seq.entity_id\n");
  fprintf(f, "_entity_poly_seq.num\n");
  fprintf(f, "_entity_poly_seq.mon_id\n");
  fprintf(f, "_entity_poly_seq.hetero\n");
  {
    int entId = 1;
    for (size_t hi = 0; hi < chains.size(); ++hi) {
      for (int s = chains[hi].startSeq; s <= chains[hi].endSeq; ++s) {
        fprintf(f, "%d %d GLY n\n", entId, s);
      }
      entId++;
    }
  }
  fprintf(f, "#\n");

  // ── Normalization: center + scale to ±50 range ────────────────────
  // Same normalization as PdbWriter to ensure bit-identical coordinates.
  float cx = 0, cy = 0, cz = 0;
  for (const auto &b : allBeads) {
    const Cluster &cl = hc.clusters[b.clusterIdx];
    cx += cl.pos.x; cy += cl.pos.y; cz += cl.pos.z;
  }
  if (!allBeads.empty()) {
    cx /= allBeads.size(); cy /= allBeads.size(); cz /= allBeads.size();
  }
  float maxDist = 1.0f;
  for (const auto &b : allBeads) {
    const Cluster &cl = hc.clusters[b.clusterIdx];
    float dx = cl.pos.x - cx, dy = cl.pos.y - cy, dz = cl.pos.z - cz;
    float d = sqrtf(dx*dx + dy*dy + dz*dz);
    if (d > maxDist) maxDist = d;
  }
  float scale = 50.0f / maxDist;

  // ── Write atom_site records ───────────────────────────────────────
  writeAtomSiteHeader(f);

  int atomId = 1;
  int entityId = 1;
  int prevChrIdx = -1;
  int seqId = 1;
  for (const auto &b : allBeads) {
    if (b.chrIdx != prevChrIdx) {
      seqId = 1;
      if (prevChrIdx >= 0)
        entityId++;
      prevChrIdx = b.chrIdx;
    }
    const Cluster &cl = hc.clusters[b.clusterIdx];

    // Normalized coordinates
    float nx = (cl.pos.x - cx) * scale;
    float ny = (cl.pos.y - cy) * scale;
    float nz = (cl.pos.z - cz) * scale;
    // Round through PDB %8.3f format for bit-identical coords with PDB writer
    char bx[9]={}, by[9]={}, bz[9]={};
    snprintf(bx, 9, "%8.3f", (double)nx);
    snprintf(by, 9, "%8.3f", (double)ny);
    snprintf(bz, 9, "%8.3f", (double)nz);
    float rx=0, ry=0, rz=0;
    sscanf(bx, "%f", &rx); sscanf(by, "%f", &ry); sscanf(bz, "%f", &rz);

    // Write N, CA, C atoms per residue.
    // Chimera's mmCIF parser needs backbone atoms (C→N bonds) to build
    // the residueAfter chain used for _struct_conf HELIX assignment.
    // N and C are placed at the same coordinates as CA (co-located).
    fprintf(f,
            "ATOM %d N N . GLY %s %d %d ? %.3f %.3f %.3f 1.00 %.2f %d %s\n",
            atomId, b.chainId.c_str(), entityId, seqId, rx, ry, rz,
            (float)cl.genomic_pos / 1e6f, cl.genomic_pos,
            hc.chrs[b.chrIdx].c_str());
    atomId++;
    fprintf(f,
            "ATOM %d C CA . GLY %s %d %d ? %.3f %.3f %.3f 1.00 %.2f %d %s\n",
            atomId, b.chainId.c_str(), entityId, seqId, rx, ry, rz,
            (float)cl.genomic_pos / 1e6f, cl.genomic_pos,
            hc.chrs[b.chrIdx].c_str());
    atomId++;
    fprintf(f,
            "ATOM %d C C . GLY %s %d %d ? %.3f %.3f %.3f 1.00 %.2f %d %s\n",
            atomId, b.chainId.c_str(), entityId, seqId, rx, ry, rz,
            (float)cl.genomic_pos / 1e6f, cl.genomic_pos,
            hc.chrs[b.chrIdx].c_str());
    atomId++;
    seqId++;
  }

  fprintf(f, "#\n");
  fclose(f);

  printf("[CIF] Wrote %d atoms to %s\n", atomId - 1, filename.c_str());
  return true;
}

bool CifWriter::writeFlat(const std::vector<vector3> &positions,
                          const std::vector<int> &genomic_positions,
                          const std::string &chrName,
                          const std::string &filename) {
  FILE *f = fopen(filename.c_str(), "w");
  if (!f) {
    printf("[CIF] Failed to open %s for writing\n", filename.c_str());
    return false;
  }

  writeHeader(f, "cudaMMC");

  // Write genomic metadata for flat output too
  std::vector<std::string> emptyfactors;
  std::vector<std::string> chrs = {chrName};
  writeGenomicMetadata(f, emptyfactors, chrs);
  writeAtomSiteHeader(f);

  for (size_t i = 0; i < positions.size(); ++i) {
    int gpos = (i < genomic_positions.size()) ? genomic_positions[i] : (int)i;
    // pdbx_auth_seq_id carries the exact genomic position in bp (lossless)
    fprintf(f,
            "ATOM %d C CA . GLY A 1 %d ? %.3f %.3f %.3f 1.00 %.2f %d %s\n",
            (int)(i + 1), (int)(i + 1), positions[i].x, positions[i].y,
            positions[i].z, (float)gpos / 1e6f, gpos, chrName.c_str());
  }

  fprintf(f, "#\n");
  fclose(f);

  printf("[CIF] Wrote %d atoms to %s\n", (int)positions.size(),
         filename.c_str());
  return true;
}
