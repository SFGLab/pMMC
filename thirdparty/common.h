/**
 * @file common.h
 * @brief Utility functions, constants, and helpers used throughout cudaMMC.
 *
 * Provides formatted string output, random number generation, vector math,
 * file I/O helpers, spline interpolation, and miscellaneous utilities.
 *
 * @author psz (original), extended for cudaMMC
 * @date Aug 4, 2013
 */

#ifndef COMMON_H_
#define COMMON_H_

#include <algorithm>
#include <assert.h>
#include <cstring>
#include <map>
#include <math.h>
#include <memory>
#include <sstream>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "mtxlib.h"

using namespace std;

/** @brief Small floating-point tolerance for comparisons. */
#define epsilon 1e-6

/** @name Hierarchy Level Constants
 *  Define the 4-level hierarchy used in the reconstruction tree.
 */
///@{
#define LVL_CHROMOSOME 0       /**< Level 0: whole chromosome. */
#define LVL_SEGMENT 1          /**< Level 1: segment (~2Mb). */
#define LVL_INTERACTION_BLOCK 2 /**< Level 2: interaction block / anchor. */
#define LVL_SUBANCHOR 3        /**< Level 3: subanchor (~10kb). */
///@}

/**
 * @brief printf-style formatted string (returns std::string).
 * @param fmt_str Format string with printf-style placeholders.
 * @return Formatted std::string.
 */
string ftext(const string fmt_str, ...);

/**
 * @brief Generate a random integer in [0, range).
 * @param range Upper bound (exclusive).
 * @return Random integer.
 */
int random(int range);

/**
 * @brief Generate a random float in [0, 1).
 * @return Uniform random float.
 */
float random_uniform();

/**
 * @brief Generate a random float in [0, range) or [-range, range).
 * @param range Maximum magnitude.
 * @param negative If true, allow negative values.
 * @return Random float.
 */
float random(float range, bool negative = false);

/**
 * @brief Generate a random 3D vector with components in [-max_size, max_size].
 * @param max_size Maximum component magnitude.
 * @param in2D If true, z-component is zero.
 * @return Random 3D vector.
 */
vector3 random_vector(float max_size, bool in2D = false);

/**
 * @brief Generate a random rotation matrix.
 * @param max_angle Maximum rotation angle in degrees (default 180).
 * @return Random 4x4 rotation matrix.
 */
matrix44 random_rot_matrix(float max_angle = 180.0f);

/**
 * @brief Euclidean distance between two 3D points.
 * @param v1 First point.
 * @param v2 Second point.
 * @return Distance.
 */
float dist(const vector3 &v1, const vector3 &v2);

/**
 * @brief Angle (in radians) between two 3D vectors.
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return Angle in radians.
 */
float angle(vector3 v1, vector3 v2);

/**
 * @brief Angle between two already-normalized vectors.
 * @param v1 First unit vector.
 * @param v2 Second unit vector.
 * @return Angle in radians.
 */
float angle_norm(vector3 v1, vector3 v2);

/**
 * @brief Displace a point by a random amount.
 * @param v Original point (value, not modified).
 * @param displacement Maximum displacement magnitude.
 * @param in2D If true, displace only in XY.
 * @return Displaced point.
 */
vector3 displace(vector3, float displacement, bool in2D = false);

/**
 * @brief Displace a point in place by a random amount.
 * @param v Point to modify (reference).
 * @param displacement Maximum displacement magnitude.
 */
void displace_ref(vector3 &v, float displacement);

/**
 * @brief Convert 3D spatial distance to interaction frequency.
 * @param distance Spatial distance between beads.
 * @return Estimated interaction frequency.
 */
float distanceToInteractionFrequency(float distance);

/**
 * @brief Linear interpolation between two floats.
 * @param a Start value.
 * @param b End value.
 * @param p Interpolation parameter [0, 1].
 * @return Interpolated value: a*(1-p) + b*p.
 */
float interpolate(float a, float b, float p);

/**
 * @brief Linear interpolation between two 3D points.
 * @param a Start point.
 * @param b End point.
 * @param p Interpolation parameter [0, 1].
 * @return Interpolated 3D point.
 */
vector3 interpolate(vector3 a, vector3 b, float p);

/**
 * @brief Fisher-Yates shuffle of a float array.
 * @param arr Array to shuffle.
 * @param size Array length.
 */
void random_shuffle(float *arr, int size);

/** @name Print Utilities
 *  Formatted printing of vectors and maps to stdout.
 */
///@{
void printv(const vector<int> &v, bool count = false, bool newline = false,
            const char *label = "");
void printv(const vector<float> &v, bool count = false, bool newline = false,
            const char *label = "");
void printv(const vector<string> &v, bool count = false, bool newline = false,
            const char *label = "");
void printv(vector<bool> v);
void printm(const std::map<int, int> &map, const char *label = "");
void printm(const std::map<int, vector<int>> &map, const char *label = "");
void printm(const std::map<int, vector<float>> &map, const char *label = "");
void printm(const std::map<std::string, int> &map, bool flat = false,
            bool newline = false, const char *label = "");
void printm(std::map<std::string, int> &map,
            const std::vector<std::string> &keys, bool flat = false,
            bool newline = false, const char *label = "");
void printm(const std::map<std::string, std::vector<int>> &map,
            bool flat = false, const char *label = "");
///@}

/**
 * @brief Print a 3D vector to stdout with an optional label.
 * @param v The vector to print.
 * @param label Optional label prefix.
 */
void print_vector(const vector3 &v, const char *label = "");

/**
 * @brief Return true with probability 'chance'.
 * @param chance Probability in [0, 1].
 * @return True with the given probability.
 */
bool withChance(float chance);

/** @name CSV I/O
 *  Save vectors and 2D arrays to CSV files.
 */
///@{
void saveToCSV(string path, vector<int> v, bool add_index_column = false);
void saveToCSV(string path, vector<float> v, bool add_index_column = false);
void saveToCSV(string path, vector<bool> v, bool add_index_column = false);
void saveToCSV(string path, vector<vector<float>> v,
               bool add_index_column = false);
void saveToCSV(string path, vector<vector<int>> v,
               bool add_index_column = false);
///@}

/**
 * @brief Check if a file exists.
 * @param name File path.
 * @return True if the file exists and is readable.
 */
bool file_exists(const std::string &name);

/**
 * @brief Open a file with error checking.
 * @param path File path.
 * @param mode Open mode (default "r").
 * @return FILE pointer (aborts on failure).
 */
FILE *open(string path, string mode = "r");

/**
 * @brief Count the number of whitespace-delimited words in a line.
 * @param line Input C-string.
 * @return Word count.
 */
int countWords(char *line);

/**
 * @brief Interpolate intermediate points between two positions.
 * @param posl Left position.
 * @param posr Right position.
 * @param dl Left distance.
 * @param dr Right distance.
 * @param mn Minimum spacing factor (default 1.5).
 * @return Vector of interpolated integer positions.
 */
vector<int> interpolateInterval(int posl, int posr, int dl, int dr,
                                float mn = 1.5f);

/**
 * @brief Interpolate a full set of points.
 * @param points Input points.
 * @param mn Minimum spacing factor.
 * @return Vector of interpolated positions.
 */
vector<int> interpolatePoints(vector<int> points, float mn = 1.5f);

/**
 * @brief Check if a vector contains a value.
 * @param v Vector to search.
 * @param val Value to find.
 * @return True if val is in v.
 */
bool vector_contains(const std::vector<int> &v, int val);

/**
 * @brief Insert a value into a vector only if not already present.
 * @param v Vector to modify.
 * @param val Value to insert.
 */
void vector_insert_unique(std::vector<int> &v, int val);

/**
 * @brief Print an error message and abort.
 * @param msg Error message.
 */
void error(string msg);

/** @name Spline Interpolation
 *  Catmull-Rom and centripetal spline interpolation for 3D curves.
 */
///@{
/**
 * @brief Catmull-Rom spline interpolation at parameter t.
 * @param t Parameter in [0, 1].
 * @param p1 Control point before the segment.
 * @param p2 Segment start point.
 * @param p3 Segment end point.
 * @param p4 Control point after the segment.
 * @return Interpolated 3D point.
 */
vector3 interpolateSpline(float t, const vector3 &p1, const vector3 &p2,
                          const vector3 &p3, const vector3 &p4);

/**
 * @brief Catmull-Rom spline through a list of points.
 * @param points Input control points.
 * @param splits Number of interpolated points per segment.
 * @return Vector of interpolated 3D points.
 */
std::vector<vector3> interpolateSpline(const std::vector<vector3> &points,
                                       int splits);

/**
 * @brief Centripetal Catmull-Rom spline interpolation.
 * @param points Input control points (modified for endpoint duplication).
 * @param splits Number of interpolated points per segment.
 * @return Vector of interpolated 3D points.
 */
std::vector<vector3> interpolateSplineCentripetal(std::vector<vector3> &points,
                                                  int splits);

/**
 * @brief Centripetal spline interpolation at parameter t.
 * @param p0 Control point 0.
 * @param p1 Control point 1.
 * @param p2 Control point 2.
 * @param p3 Control point 3.
 * @param time Knot parameters array.
 * @param t Interpolation parameter.
 * @return Interpolated 3D point.
 */
vector3 interpolate(vector3 &p0, vector3 &p1, vector3 &p2, vector3 &p3,
                    double time[], double t);
///@}

/**
 * @brief Mirror a point across a fixed point.
 * @param fixed The fixed (mirror) point.
 * @param pt The point to mirror.
 * @return Mirrored point: 2*fixed - pt.
 */
vector3 mirrorPoint(vector3 fixed, vector3 pt);

/** @name String Splitting
 *  Split strings by delimiter.
 */
///@{
/**
 * @brief Split a string by delimiter and append to output vector.
 * @param s Input string.
 * @param delim Delimiter character.
 * @param elems Output vector (appended to).
 */
void split(const std::string &s, char delim, std::vector<std::string> &elems);

/**
 * @brief Split a string by delimiter.
 * @param s Input string.
 * @param delim Delimiter (default ',').
 * @return Vector of tokens.
 */
std::vector<std::string> split(const std::string &s, char delim = ',');

/**
 * @brief Split a string by any whitespace.
 * @param s Input string.
 * @return Vector of tokens.
 */
std::vector<std::string> split_unknown(const std::string &s);
///@}

/**
 * @brief Split a file path into directory and filename components.
 * @param s Full file path.
 * @param[out] dir Directory part.
 * @param[out] filename Filename part.
 */
void split_file_path(const std::string &s, std::string &dir,
                     std::string &filename);

#endif /* COMMON_H_ */
