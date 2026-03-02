/**
 * @file Heatmap.h
 * @brief 2D interaction frequency matrix (contact map).
 *
 * Stores a symmetric NxN matrix of interaction frequencies between
 * genomic bins. Used at both the segment level (singleton heatmap)
 * and the chromosome level during Monte Carlo reconstruction.
 *
 * @see LooperSolver::reconstructClustersHeatmap for heatmap-guided MC.
 */

#ifndef HEATMAP_H_
#define HEATMAP_H_

#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include <locale>

#include <common.h>
#include <Settings.h>

using namespace std;

/**
 * @class Heatmap
 * @brief A square symmetric matrix of interaction frequencies.
 *
 * The matrix is stored as a 2D C array (v[i][j]) dynamically allocated.
 * Supports I/O, arithmetic operations, smoothing, distance conversion,
 * and extraction of diagonal statistics.
 */
class Heatmap {

public:
  /** @brief Default constructor (empty, size=0). */
  Heatmap();

  /** @brief Destructor — frees the NxN matrix (REQ-3.7 RAII). */
  ~Heatmap();

  /** @brief Copy constructor — deep-copies the matrix. */
  Heatmap(const Heatmap &other);

  /** @brief Copy assignment — deep-copies the matrix. */
  Heatmap &operator=(const Heatmap &other);

  /**
   * @brief Construct and allocate an NxN heatmap.
   * @param size Number of bins (rows = columns).
   */
  Heatmap(int size);

  /**
   * @brief Allocate (or reallocate) the matrix to NxN.
   * @param size Number of bins.
   */
  void init(int size);

  /** @brief Set all values to zero. */
  void zero();

  /** @brief Print the matrix to stdout. */
  void print();

  /**
   * @brief Write the heatmap to a text file.
   * @param filename Output file path.
   * @param total_count If true, include a header with the total contact count.
   * @param zero_diagonal If true, set diagonal values to zero before writing.
   * @param as_integers If true, write values as integers.
   */
  void toFile(string filename, bool total_count = true,
              bool zero_diagonal = false, bool as_integers = false);

  /**
   * @brief Read a heatmap from a text file (space-delimited matrix).
   * @param filename Input file path.
   */
  void fromFile(string filename);

  /**
   * @brief Read a heatmap with optional row/column labels.
   * @param filename Input file path.
   * @param labels True if the file has a header row and label column.
   */
  void fromFile(string filename, bool labels);

  /**
   * @brief Read a heatmap from an MDS-format file.
   * @param filename Input file path.
   * @param size Expected matrix size.
   */
  void fromMDS(string filename, int size);

  /**
   * @brief Check if the heatmap is entirely zero.
   * @return True if all values are zero.
   */
  bool isEmpty();

  /**
   * @brief Get a per-column mask of empty (all-zero) columns.
   * @return Boolean vector: true if column i is all zeros.
   */
  vector<bool> getEmpty();

  /**
   * @brief Remove specified columns/rows from the heatmap.
   * @param del Boolean vector: true = remove that index.
   * @return New Heatmap with the specified rows/columns removed.
   */
  Heatmap removeColumns(vector<bool> del);

  /**
   * @brief Remove all-zero columns and rows.
   * @return New Heatmap with empty rows/columns removed.
   */
  Heatmap removeEmptyColumns();

  /**
   * @brief Multiply all values by a scalar.
   * @param scale Scale factor.
   */
  void scale(float scale);

  /**
   * @brief Element-wise add another heatmap.
   * @param heat Heatmap to add (must be same size).
   */
  void add(const Heatmap &heat);

  /**
   * @brief Add a constant to all elements.
   * @param val Value to add.
   */
  void add(float val);

  /**
   * @brief Compute element-wise difference between two heatmaps.
   * @param heat Other heatmap (same size).
   * @param abs If true, take absolute differences.
   * @return New Heatmap with difference values.
   */
  Heatmap diff(const Heatmap &heat, bool abs = false);

  /**
   * @brief Element-wise divide by another heatmap.
   * @param hmap Divisor heatmap (same size).
   */
  void divide(const Heatmap &hmap);

  /**
   * @brief Compute the distance (sum of squared differences) to another heatmap.
   * @param hmap Other heatmap.
   * @return Sum of squared element-wise differences.
   */
  float calcDistance(const Heatmap &hmap);

  /**
   * @brief Smooth the heatmap by averaging with neighbors.
   * @param threshold Minimum value threshold for smoothing.
   * @param factor Smoothing factor (blend weight).
   */
  void smooth(float threshold, float factor);

  /**
   * @brief Smooth with break points (do not smooth across segment boundaries).
   * @param threshold Minimum value threshold.
   * @param factor Smoothing factor.
   * @param breaks Set of indices where smoothing should not cross.
   */
  void smooth(float threshold, float factor, set<int> &breaks);

  /**
   * @brief Get the width of the zero-band along the diagonal.
   * @return Number of zero diagonals from the main diagonal.
   */
  int getDiagonalSize();

  /**
   * @brief Get the min and max values in the matrix.
   * @param[out] min Minimum value.
   * @param[out] max Maximum value.
   */
  void getRange(float &min, float &max);

  /**
   * @brief Compute the average of all non-diagonal values.
   * @return Average value.
   */
  float getAvg();

  /**
   * @brief Average value of the non-zero cells closest to the diagonal.
   * @return Average near-diagonal value.
   */
  float getAvgNearDiagonal();

  /**
   * @brief Compute average values for each diagonal offset k.
   *
   * For each k, avg_value[k] = mean of v[i][i+k] over all valid i.
   * Results are stored in the avg_value member.
   *
   * @param count_zeros If true, include zero values in the average.
   */
  void calcAvgValues(bool count_zeros = true);

  /**
   * @brief Set the main diagonal (and nearby diagonals) to zero.
   * @param diag Number of diagonals to clear (1 = only main diagonal).
   */
  void clearDiagonal(int diag = 1);

  /**
   * @brief Flatten the upper triangle to a 1D vector.
   * @param diag Skip this many diagonals from the main diagonal.
   * @return Vector of upper-triangle values.
   */
  vector<float> toVector(int diag = 0);

  size_t size;    /**< Matrix dimension (NxN). */
  float **v;      /**< 2D matrix data: v[row][col]. */

  int start;       /**< Genomic start position of the first bin. */
  int resolution;  /**< Genomic resolution (bp per bin). */

  int diagonal_size;  /**< Width of the zero-band on the diagonal. */

  std::vector<float> avg_value;  /**< Average value per diagonal offset k, set by calcAvgValues(). */

private:
  /** @brief Free the allocated matrix memory. */
  void clear();
};

#endif /* HEATMAP_H_ */
