/**
 * @file Density.h
 * @brief 3D density map for microscopy-constrained chromosome reconstruction.
 *
 * Represents a voxel grid of density values (e.g. from FISH or cryo-EM data).
 * The density map can constrain the Monte Carlo reconstruction to place
 * chromatin in regions of high experimental density.
 *
 * @see LooperSolver::MonteCarloHeatmapAndDensity for density-guided simulation.
 */

#ifndef DENSITY_H_
#define DENSITY_H_

#include <common.h>

/**
 * @class Density
 * @brief A 3D voxel grid of density values.
 *
 * Supports I/O in both raw 3D array format and segmentation format
 * (x, y, z, density per line). Can be downsampled by averaging
 * neighboring voxels.
 */
class Density {
public:
  /** @brief Default constructor. */
  Density();

  /**
   * @brief Read a 3D density array from a binary file.
   * @param filename Input file path.
   */
  void fromFile(string filename);

  /**
   * @brief Write the 3D density array to a binary file.
   * @param filename Output file path.
   */
  void toFile(string filename);

  /**
   * @brief Write as a segmentation file (one "x y z density" per line).
   * @param filename Output file path.
   */
  void toSegmentationFile(string filename);

  /**
   * @brief Read from a segmentation file (one "x y z density" per line).
   * @param filename Input file path.
   */
  void fromSegmentationFile(string filename);

  /**
   * @brief Initialize an empty density grid.
   * @param size_x Grid size in X dimension.
   * @param size_y Grid size in Y dimension.
   * @param size_z Grid size in Z dimension.
   * @param is_static True for fixed input density, false for dynamic structure density.
   */
  void init(int size_x, int size_y, int size_z, bool is_static);

  /**
   * @brief Initialize with explicit origin offsets.
   * @param size_x Grid size in X.
   * @param size_y Grid size in Y.
   * @param size_z Grid size in Z.
   * @param start_x Origin X offset.
   * @param start_y Origin Y offset.
   * @param start_z Origin Z offset.
   * @param is_static True for fixed input density.
   */
  void init(int size_x, int size_y, int size_z, int start_x, int start_y,
            int start_z, bool is_static);

  /** @brief Normalize density values to sum to 1.0. */
  void normalize();

  /** @brief Set all voxels to zero. */
  void clear();

  /** @brief Print grid dimensions and ranges to stdout. */
  void print();

  /**
   * @brief Create a downsampled copy by averaging voxels.
   * @param factor Downsampling factor (e.g. 2 = halve each dimension).
   * @return New Density at lower resolution.
   */
  Density scaleDown(int factor);

  /**
   * @brief Add a point mass to the density grid (with Gaussian spread).
   * @param x Voxel X coordinate.
   * @param y Voxel Y coordinate.
   * @param z Voxel Z coordinate.
   * @param val Density value to add.
   */
  void addPointMass(int x, int y, int z, float val);

  vector<vector<vector<float>>> t;  /**< 3D density voxel data: t[x][y][z]. */

  int size_x, size_y, size_z;  /**< Grid dimensions. */

  int range_x_start, range_x_end;  /**< X-axis data range. */
  int range_y_start, range_y_end;  /**< Y-axis data range. */
  int range_z_start, range_z_end;  /**< Z-axis data range. */
  float range_d_start, range_d_end;  /**< Density value range. */

  vector3 center;  /**< Center of mass of the density map. */
  vector3 origin;  /**< Origin offset in world coordinates. */

  float scale;  /**< Scale factor from settings (for density visualization). */

private:
  /** @brief Internal BFS helper key. */
  struct key_3d {
    int x, y, z;
    float value;
  };

  /**
   * @brief Check if voxel coordinates are within the grid bounds.
   * @param x X coordinate.
   * @param y Y coordinate.
   * @param z Z coordinate.
   * @return True if inside the grid.
   */
  bool isInside(int x, int y, int z);

  /** @brief Clear the visited-flag array (used by BFS). */
  void clearOdw();

  vector<vector<vector<bool>>> odw;  /**< Visited flags for BFS traversal. */

  bool is_static;  /**< True for input density (fixed), false for structure-derived (dynamic). */

  const int dx[6] = {-1, 1, 0, 0, 0, 0};  /**< BFS X-direction offsets. */
  const int dy[6] = {0, 0, 1, -1, 0, 0};  /**< BFS Y-direction offsets. */
  const int dz[6] = {0, 0, 0, 0, 1, -1};  /**< BFS Z-direction offsets. */
};

#endif /* DENSITY_H_ */
