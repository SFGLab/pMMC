/**
 * @file Chromosome.h
 * @brief 3D polymer structure representing a chromosome or fragment thereof.
 *
 * Stores a sequence of 3D points (beads), supports alignment, RMSD calculation,
 * heatmap generation, and various geometric operations. This is the primary
 * output structure at each level of the hierarchy.
 *
 * @see HierarchicalChromosome for the multi-level container.
 */

#ifndef CHROMOSOME_H_
#define CHROMOSOME_H_

#include <string.h>
#include <vector>

#include <locale>

#include <common.h>
#include <rmsd.h>
#include <Heatmap.h>
#include <Settings.h>

using namespace std;

/**
 * @class Chromosome
 * @brief A 3D bead-chain polymer model of a chromosome (or sub-region).
 *
 * The points vector stores the 3D coordinates of each bead.
 * Genomic positions can be stored in genomic_position.
 * Supports I/O (binary, PDB), random generation, alignment via
 * Procrustes rotation, and conversion to/from heatmaps.
 */
class Chromosome {

public:
  /** @brief Default constructor (empty structure). */
  Chromosome();

  /** @brief Reset all fields. */
  void init();

  /** @brief Print structure summary to stdout. */
  void print();

  /**
   * @brief Write structure to a binary file.
   * @param filename Output file path.
   */
  void toFile(string filename);

  /**
   * @brief Write structure to an already-open file handle.
   * @param file Open FILE pointer.
   */
  void toFile(FILE *file);

  /**
   * @brief Read structure from a binary file.
   * @param filename Input file path.
   */
  void fromFile(string filename);

  /**
   * @brief Read structure from an already-open file handle.
   * @param file Open FILE pointer.
   * @param pts_cnt Number of points to read (0 = read from file header).
   */
  void fromFile(FILE *file, int pts_cnt = 0);

  /**
   * @brief Read structure from a PDB file.
   * @param filename Path to PDB file.
   */
  void fromFilePDB(string filename);

  /**
   * @brief Generate a random polymer structure.
   * @param pts_cnt Number of beads.
   * @param size Step size for random walk.
   * @param walk If true, use random walk; if false, place randomly.
   */
  void createRandom(int pts_cnt, float size = 0.1f, bool walk = true);

  /**
   * @brief Create a straight-line structure.
   * @param pts_cnt Number of beads.
   * @param dist Distance between consecutive beads.
   */
  void makeLine(int pts_cnt, float dist);

  /**
   * @brief Create a helical (spiral) structure.
   * @param pts_cnt Number of beads.
   * @param r Helix radius.
   * @param angle Angular step per bead (radians).
   * @param spin_height Height increment per revolution.
   */
  void makeSpiral(int pts_cnt, float r, float angle, float spin_height);

  /**
   * @brief Create a spiral between two endpoints.
   * @param pts_cnt Number of beads.
   * @param r Helix radius.
   * @param angle Angular step per bead.
   * @param p1 Start point.
   * @param p2 End point.
   */
  void makeSpiral(int pts_cnt, float r, float angle, vector3 p1, vector3 p2);

  /**
   * @brief Assemble this chromosome from a list of sub-chromosomes.
   * @param chrs Vector of sub-chromosome structures to concatenate.
   */
  void createFromSubchromosomes(const vector<Chromosome> &chrs);

  /**
   * @brief Randomize positions of sub-chromosome fragments.
   * @param dispersion Amount of random displacement.
   * @param keep_subchr_structure If true, preserve internal sub-chromosome geometry.
   */
  void randomizeSubchromosomes(float dispersion,
                               bool keep_subchr_structure = true);

  /**
   * @brief Assign sub-chromosome membership indices to beads.
   * @param subchr_index Vector mapping each bead to its sub-chromosome index.
   */
  void setSubchromosomesIndices(const vector<int> &subchr_index);

  /**
   * @brief Create a deep copy of this chromosome.
   * @return A new Chromosome with copied data.
   */
  Chromosome clone();

  /** @brief Translate the structure so its center of mass is at the origin. */
  void center();

  /**
   * @brief Scale all coordinates by a factor.
   * @param scale Scale factor.
   * @param center If true, center the structure first.
   */
  void scale(float scale, bool center = false);

  /**
   * @brief Create a contact-frequency heatmap from pairwise 3D distances.
   * @return Heatmap where h[i][j] ~ 1/dist(i,j).
   */
  Heatmap createHeatmap();

  /**
   * @brief Create an inverse heatmap (simulated "inverse Hi-C").
   * @return Heatmap with distance-based values.
   */
  Heatmap createInverseHeatmap();

  /**
   * @brief Align this structure to a reference by Procrustes rotation.
   * @param chr Reference chromosome to align to.
   * @param resize If true, allow scaling during alignment.
   * @param max_angle Maximum rotation angle (radians, default 2*pi).
   */
  void align(const Chromosome &chr, bool resize = false,
             float max_angle = 2 * 3.14f);

  /**
   * @brief Rescale this structure to match the size of a reference.
   * @param chr Reference chromosome.
   */
  void adjustSize(const Chromosome &chr);

  /** @brief Recompute the size field from the points array. */
  void updateSize();

  /**
   * @brief Compute the sum of squared distances to another structure.
   * @param chr Other chromosome (must have same number of points).
   * @return Sum of squared distances between corresponding beads.
   */
  float getDistanceSqr(const Chromosome &chr);

  /**
   * @brief Compute the center of mass.
   * @return 3D center of mass vector.
   */
  vector3 getCenter();

  /**
   * @brief Translate all points by a vector.
   * @param v Translation vector.
   */
  void translate(const vector3 &v);

  /**
   * @brief Rotate all points by a 4x4 transformation matrix.
   * @param mat Rotation/transformation matrix.
   */
  void rotate(const matrix44 &mat);

  /**
   * @brief Compute the diameter (maximum pairwise distance).
   * @return Diameter in spatial units.
   */
  float getDiameter() const;

  /**
   * @brief Get distances between consecutive beads.
   * @return Vector of distances: dist(point[i], point[i+1]).
   */
  std::vector<float> getConsecutiveBeadDistances();

  /**
   * @brief Average distance between consecutive beads.
   * @return Mean of getConsecutiveBeadDistances().
   */
  float getAvgConsecutiveBeadDistance();

  /**
   * @brief Average Euclidean distance to another structure.
   * @param chr Other chromosome.
   * @return Mean distance across corresponding bead pairs.
   */
  float calcDistance(const Chromosome &chr);

  /**
   * @brief Root-mean-square deviation to another structure.
   * @param chr Other chromosome (must have same size).
   * @return RMSD value.
   */
  float calcRMSD(const Chromosome &chr);

  /**
   * @brief Calculate base density (number of beads within a sphere around each bead).
   * @param sphere_r Sphere radius (0 = use default from settings).
   * @return 2D vector: [bead_index][density_values].
   */
  vector<vector<float>> calcBaseDensity(float sphere_r = 0.0f);

  /**
   * @brief Trim the structure to a sub-range of beads.
   * @param start First bead index to keep.
   * @param end Last bead index to keep (0 = keep to the end).
   */
  void trim(int start, int end = 0);

  int size;  /**< Number of beads (points). */

  vector<vector3> points;              /**< 3D coordinates of each bead. */
  vector<int> genomic_position;        /**< Genomic position (bp) per bead. */
  vector<vector<vector3>> points_hierarchical;  /**< Multi-level bead positions for visualization. */

  float score;  /**< Energy score from the last Monte Carlo optimization. */

private:
  /**
   * @brief Grid search for the best alignment rotation.
   * @param chr Reference chromosome.
   * @param steps Number of angular steps per axis.
   * @param angle Best rotation angles (output).
   * @param max_angle Maximum rotation angle.
   * @return Minimum distance after rotation.
   */
  float findBestAlignmentRotation(const Chromosome &chr, int steps,
                                  vector3 &angle, float max_angle);

  /**
   * @brief Find the pair of points with the largest distance.
   * @param dist Output: the maximum distance found.
   * @return Index of one of the two most distant points.
   */
  int getMostDistantPointsPair(float &dist);
};

#endif /* CHROMOSOME_H_ */
