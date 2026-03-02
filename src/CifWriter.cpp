/**
 * @file CifWriter.cpp
 * @brief mmCIF format writer implementation.
 */

#include <CifWriter.h>
#include <HierarchicalChromosome.h>
#include <cstdio>
#include <cstring>

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
                      const std::string &dataName) {
  FILE *f = fopen(filename.c_str(), "w");
  if (!f) {
    printf("[CIF] Failed to open %s for writing\n", filename.c_str());
    return false;
  }

  writeHeader(f, dataName);
  writeAtomSiteHeader(f);

  int atomId = 1;
  int entityId = 1;

  // Use lowest level for maximum detail
  hc.useLowestLevel();

  for (size_t ci = 0; ci < hc.chrs.size(); ++ci) {
    std::string chr = hc.chrs[ci];
    std::string chainId = chrToChainId(chr, (int)ci);

    if (hc.current_level.find(chr) == hc.current_level.end())
      continue;

    const std::vector<int> &level = hc.current_level[chr];
    int seqId = 1;

    for (size_t j = 0; j < level.size(); ++j) {
      int idx = level[j];
      if (idx < 0 || idx >= (int)hc.clusters.size())
        continue;

      const Cluster &cl = hc.clusters[idx];
      fprintf(f,
              "ATOM %d C CA . ALA %s %d %d ? %.3f %.3f %.3f 1.00 %.2f %s\n",
              atomId, chainId.c_str(), entityId, seqId, cl.pos.x, cl.pos.y,
              cl.pos.z, (float)cl.genomic_pos / 1000.0f, chr.c_str());

      atomId++;
      seqId++;
    }
    entityId++;
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
  writeAtomSiteHeader(f);

  for (size_t i = 0; i < positions.size(); ++i) {
    int gpos = (i < genomic_positions.size()) ? genomic_positions[i] : (int)i;
    fprintf(f, "ATOM %d C CA . ALA A 1 %d ? %.3f %.3f %.3f 1.00 %.2f %s\n",
            (int)(i + 1), (int)(i + 1), positions[i].x, positions[i].y,
            positions[i].z, (float)gpos / 1000.0f, chrName.c_str());
  }

  fprintf(f, "#\n");
  fclose(f);

  printf("[CIF] Wrote %d atoms to %s\n", (int)positions.size(),
         filename.c_str());
  return true;
}
