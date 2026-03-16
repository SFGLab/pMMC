/**
 * @file LooperSolver.h
 * @brief Main reconstruction engine for 3D genome structure from ChIA-PET data.
 *
 * Implements the Multiscale Monte Carlo (MMC) algorithm from:
 *   Szalaj et al., "3D-GNOME: 3D Genome Organization Modeling Engine",
 *   Nucleic Acids Research, 2016.
 *
 * The reconstruction proceeds top-down through four hierarchy levels:
 *   1. Chromosome level — heatmap-guided MC (singleton data)
 *   2. Segment level (~2Mb) — heatmap-guided MC
 *   3. Anchor level — arc-distance-guided MC (PET cluster data)
 *   4. Subanchor level (~10kb) — smoothing MC (linker lengths + angles)
 *
 * @see HierarchicalChromosome for the output structure.
 * @see Settings for all configurable parameters.
 */

#ifndef LOOPERSOLVER_H_
#define LOOPERSOLVER_H_

#include <fstream>
#include <queue>
#include <stdarg.h>
#include <string>
#include <vector>

#include <BedRegion.h>
#include <BedRegions.h>
#include <ChromosomesSet.h>
#include <Cluster.h>
#include <Density.h>
#include <Heatmap.h>
#include <HierarchicalChromosome.h>
#include <InteractionArc.h>
#include <InteractionArcs.h>
#include <MultiscaleEnergy.h>
#include <SimulationLogger.h>
#include <platform.h>

// Forward declaration
class SimulationLogger;

/**
 * @brief Persistent GPU resource cache to avoid per-block alloc/free overhead.
 *
 * All CUDA types stored as void* to avoid including cuda_runtime.h in the
 * header (which would conflict with HierarchicalChromosome.h's float3 stub).
 */
struct GpuResourceCache {
  bool initialized = false;

  // Device properties (queried once)
  int sm_count = 0;
  int dev_id = 0;
  // cudaDeviceProp stored opaquely — accessed only from .cu files
  char dev_prop_storage[4096] = {};  // sizeof(cudaDeviceProp) varies by CUDA version

  // Persistent curandState for arcs kernel (void* to avoid CUDA header dep)
  void *d_arcs_states = nullptr;
  int arcs_states_count = 0;
  bool arcs_states_seeded = false;

  // Persistent curandState for smooth kernel
  void *d_smooth_states = nullptr;
  int smooth_states_count = 0;
  bool smooth_states_seeded = false;

  // Persistent isDone flags (device booleans)
  void *d_arcs_isDone = nullptr;
  void *d_smooth_isDone = nullptr;

  // Persistent device buffers for arcs kernel (avoid per-call thrust alloc)
  void *d_arcs_positions = nullptr;   // float3[max_n]
  void *d_arcs_is_fixed = nullptr;    // int[max_n]
  void *d_arcs_exp_dist = nullptr;    // float[max_n * max_n]
  void *d_arcs_loop_pairs = nullptr;  // int2[max_loops]
  void *d_arcs_loop_params = nullptr; // float2[max_loops]
  void *d_arcs_best_score = nullptr;  // float[max_warps]
  void *d_arcs_best_pos = nullptr;    // float3[max_warps * max_n]
  int arcs_max_n = 0;
  int arcs_max_warps = 0;
  int arcs_max_loops = 0;

  // Persistent device buffers for smooth kernel
  void *d_smooth_positions = nullptr; // float3[max_n]
  void *d_smooth_fixed = nullptr;     // bool[max_n] (stored as char)
  void *d_smooth_dist = nullptr;      // float[max_n]
  int smooth_max_n = 0;

  void init();
  void initArcsBuffers(int max_n, int max_warps, int max_loops);
  void initSmoothBuffers(int max_n);
  void cleanup();
  ~GpuResourceCache() { cleanup(); }
};

/**
 * @class LooperSolver
 * @brief The core 3D genome reconstruction engine using multiscale Monte Carlo.
 *
 * Usage:
 * @code
 *   LooperSolver lsm("label", "./output/");
 *   lsm.setContactData(chrs, region, anchors, factors, clusters, singletons, singletons_inter);
 *   lsm.createTreeGenome();
 *   lsm.reconstructClustersHeatmap();      // Levels 0-1
 *   lsm.reconstructClustersArcsDistances(); // Levels 2-3
 *   HierarchicalChromosome model = lsm.getModel();
 * @endcode
 */
class LooperSolver {
public:
  /**
   * @brief Construct a solver with a label and output directory.
   * @param label Label for output filenames (e.g. "run1").
   * @param outdir Output directory path.
   */
  LooperSolver(string label, string outdir);

  /**
   * @brief Print a formatted log message at a given verbosity level.
   * @param level Minimum output level required to print.
   * @param format printf-style format string.
   */
  void output(int level, const char *format, ...);

  /** @brief Print a summary of the current state to stdout. */
  void print();

  /** @brief Print the active region clusters. */
  void printActiveRegion();

  /**
   * @brief Load all contact data (anchors, PET clusters, singletons).
   *
   * This is the first step: it loads anchors, PET clusters for each factor,
   * and singleton data for heatmap reconstruction.
   *
   * @param chrs_list List of chromosome names to reconstruct.
   * @param region_of_interest Restrict to this genomic region (empty = whole genome).
   * @param anchors Path to the anchors file.
   * @param factors List of factor names (e.g. {"CTCF"}).
   * @param arcs_clusters Paths to PET cluster files (one per factor).
   * @param arcs_singletons Paths to singleton files (intra-chromosomal).
   * @param arcs_singletons_inter Paths to singleton files (inter-chromosomal).
   */
  void setContactData(std::vector<string> chrs_list,
                      BedRegion region_of_interest, string anchors,
                      std::vector<string> factors,
                      std::vector<string> arcs_clusters,
                      std::vector<string> arcs_singletons,
                      std::vector<string> arcs_singletons_inter);

  /** @brief Initialize density data from microscopy (if useDensity is enabled). */
  void initDensityData();

  /** @brief Load loop constraints from a text file (public entry point).
   * Format: cluster_i cluster_j [stiffness] [eq_distance]
   */
  void loadLoopConstraints(const std::string &filepath);

  /**
   * @brief Build the hierarchical cluster tree for the entire genome.
   *
   * Creates the 4-level tree: chromosome roots -> segments -> anchors -> subanchors.
   * Calls createTreeChromosome() for each chromosome.
   */
  void createTreeGenome();

  /**
   * @brief Build the cluster tree for a single chromosome.
   * @param chr Chromosome name.
   * @return Index of the chromosome root cluster.
   */
  int createTreeChromosome(string chr);

  /**
   * @brief Reconstruct all clusters at the arc/distance level (levels 2-3).
   *
   * For each chromosome, runs MC simulation at the anchor level (using PET
   * cluster distances) and then at the subanchor level (smoothing).
   */
  void reconstructClustersArcsDistances();

  /**
   * @brief Reconstruct a single cluster's children using arc distances.
   * @param cluster Index of the parent cluster.
   * @param current_ind Current level index.
   * @param smooth If true, also perform subanchor smoothing.
   * @param use_subanchor_heatmap If true, use subanchor heatmap in scoring.
   */
  void reconstructClusterArcsDistances(int cluster, int current_ind,
                                       bool smooth = false,
                                       bool use_subanchor_heatmap = false);

  /**
   * @brief Reconstruct all clusters at the heatmap level (levels 0-1).
   *
   * Creates singleton heatmaps and runs MC simulation at the chromosome
   * and segment levels.
   */
  void reconstructClustersHeatmap();

  /**
   * @brief Reconstruct a single heatmap level.
   * @param heatmap_ind Index in the heatmap array (0 = chromosome, 1 = segment).
   */
  void reconstructClustersHeatmapSingleLevel(int heatmap_ind);

  /**
   * @brief Normalize a heatmap by its diagonal average values.
   * @param heat Input heatmap.
   * @return Normalized heatmap.
   */
  Heatmap normalizeHeatmap(const Heatmap &heat);

  /**
   * @brief Normalize an inter-chromosomal heatmap.
   * @param heat Heatmap to normalize (modified in place).
   * @param scale Normalization scale factor (default 2.0).
   */
  void normalizeHeatmapInter(Heatmap &heat, float scale = 2.0f);

  /**
   * @brief Normalize so that the total sum of near-diagonal cells equals val.
   * @param heat Heatmap to normalize (modified in place).
   * @param val Target total (default 1.0).
   */
  void normalizeHeatmapDiagonalTotal(Heatmap &heat, float val = 1.0f);

  /**
   * @brief Get average heatmap value between neighboring regions at a given level.
   * @param heat The heatmap.
   * @param level Hierarchy level.
   * @return Average contact frequency between adjacent regions.
   */
  float getHeatmapAvgBetweenNeighboringRegions(Heatmap &heat, int level);

  /**
   * @brief Create a distance heatmap from the current 3D structure.
   * @param heatmap_ind Index in the heatmap array.
   */
  void createDistanceHeatmap(int heatmap_ind);

  /**
   * @brief Convert a 3D-distance heatmap to interaction frequencies.
   * @param h Distance heatmap.
   * @return Interaction frequency heatmap.
   */
  Heatmap createIFHeatmapFromDistances(const Heatmap &h);

  /** @brief Set active level to the topmost (chromosome) level. */
  void setTopLevel();

  /**
   * @brief Set the active hierarchy level.
   * @param level 0=chromosome, 1=segment, 2=anchor, 3=subanchor.
   */
  void setLevel(int level);

  /** @brief Move one level down in the hierarchy. */
  void levelDown();

  /** @brief Reset state between consecutive ensemble runs. */
  void reset();

  /** @brief Remove subanchor-level beads from the cluster tree. */
  void removeSubanchorBeads();

  /**
   * @brief Extract the reconstructed model as a HierarchicalChromosome.
   * @return The complete multi-level 3D structure.
   */
  HierarchicalChromosome getModel();

  /** @brief Repair inconsistencies in the structure heatmap. */
  void repairStructureHeatmap();

  /**
   * @brief Find the length of unlooped chromatin around a cluster.
   * @param cluster Cluster index.
   * @return Genomic length (bp) of unlooped region.
   */
  int findUnloopedChromatinLength(int cluster);

  /**
   * @brief Convert interaction frequency to spatial distance (heatmap level).
   * @param freq Interaction frequency value.
   * @return Spatial distance.
   */
  float freqToDistanceHeatmap(float freq);

  /**
   * @brief Convert inter-chromosomal frequency to spatial distance.
   * @param freq Interaction frequency value.
   * @return Spatial distance.
   */
  float freqToDistanceHeatmapInter(float freq);

  /**
   * @brief Convert spatial distance to interaction frequency (inverse mapping).
   * @param freq Distance value (named freq for legacy reasons).
   * @return Interaction frequency.
   */
  float distToFreqHeatmap(float freq);

  /**
   * @brief Convert PET count to spatial distance (arc level).
   * @param freq PET count (interaction score).
   * @param memo If true, use memoization table.
   * @return Expected spatial distance.
   */
  float freqToDistance(int freq, bool memo = true);

  /** @brief Add empty clusters to split overly long linker segments. */
  void densify();

  /** @brief Print diagnostic information about the current state. */
  void diagnose();

  /**
   * @brief Save a snapshot of the current structure for debugging.
   * @param s Description label for the snapshot.
   */
  void getSnapshot(string s = "<no_desc>");

  /**
   * @brief Add an external chromosome structure as a snapshot.
   * @param chr Chromosome structure to save.
   * @param s Description label.
   */
  void addSnapshot(Chromosome chr, string s = "<no_desc>");

  /**
   * @brief Input an arbitrary heatmap for reconstruction (bypasses loading).
   * @param h The heatmap to use.
   */
  void inputArbitraryHeatmap(Heatmap h);

  /**
   * @brief Get the current-level chromosome structures.
   * @return Map from chromosome name to Chromosome.
   */
  std::map<std::string, Chromosome> getCurrentChromosomes();

  /**
   * @brief Find gaps (anchor-free regions) in a chromosome.
   * @param chr Chromosome name.
   * @return Vector of gap boundary indices.
   */
  std::vector<int> findGaps(string chr);

  /**
   * @brief Find gaps starting from a given cluster index.
   * @param start_index Starting cluster index.
   * @return Vector of gap boundary indices.
   */
  std::vector<int> findGaps(int start_index);

  /**
   * @brief Find the optimal segment split given gap positions.
   * @param gaps Vector of gap boundary indices.
   * @param size Target segment size in bp.
   * @param chr Chromosome name.
   * @return Vector of split point indices.
   */
  std::vector<int> findSplit(std::vector<int> gaps, int size, string chr);

  /**
   * @brief Set debugging limits for faster test runs.
   * @param chr_number_limit Maximum number of chromosomes (-1 = no limit).
   * @param length_limit Maximum genomic length per chromosome (-1 = no limit).
   */
  void setDebuggingOptions(int chr_number_limit = -1, int length_limit = -1);

  /**
   * @brief Print coordinate ranges at a specific hierarchy level.
   * @param level Hierarchy level.
   */
  void showCoordinatesRange(int level);

  /**
   * @brief Print coordinate ranges for the current or active region.
   * @param active_region If true, show only the active region.
   */
  void showCoordinatesRange(bool active_region = false);

  /** @brief Compute and print statistics about anchor positions and densities. */
  void calcAnchorsStats();

  /**
   * @brief Compute genomic length excluding the centromere region.
   * @param chr Chromosome name.
   * @param st Start position (bp).
   * @param end End position (bp).
   * @return Effective length in bp (minus centromere overlap).
   */
  int getGenomicLengthExcludingCentromere(string chr, int st, int end);

  /** @brief Find the cluster and segment containing a genomic position.
   * @param chr Chromosome name.
   * @param pos Genomic position (bp).
   */
  void printStructuresForGenomicPosition(string chr, int pos);

  /**
   * @brief Find the cluster index for a genomic position.
   * @param chr Chromosome name.
   * @param pos Genomic position (bp).
   * @return Cluster index, or -1 if not found.
   */
  int findClusterForGenomicPosition(string chr, int pos);

  /**
   * @brief Find the current-level index for a genomic position.
   * @param chr Chromosome name.
   * @param pos Genomic position (bp).
   * @return Index in current_level, or -1 if not found.
   */
  int findCurrentLevelForGenomicPosition(string chr, int pos);

  /**
   * @brief GPU-accelerated Monte Carlo for heatmap level (CUDA kernel).
   * @param step_size Maximum displacement per MC move.
   * @return Final energy score.
   * @see ParallelMonteCarloHeatmap.cu
   */
  float ParallelMonteCarloHeatmap(float step_size);

  /**
   * @brief GPU-accelerated Monte Carlo for subanchor smoothing (CUDA kernel).
   *
   * Mirrors MonteCarloArcsSmooth() but runs warp-level SA on the GPU.
   * CTCF orientation scoring is omitted; apply as a CPU post-step if needed.
   *
   * @param step_size  Initial maximum displacement per MC move.
   * @param seed       RNG seed (use time(NULL) for non-deterministic runs).
   * @return Final smooth energy score (calcScoreStructureSmooth, CPU-evaluated).
   * @see ParallelMonteCarloSmooth.cu
   */
  float ParallelMonteCarloSmooth(float step_size,
                                  unsigned long long seed = 0ULL);

  /**
   * @brief CPU Monte Carlo simulation for heatmap level.
   * @param step_size Maximum displacement per MC move.
   * @return Final energy score.
   */
  double MonteCarloHeatmap(float step_size);

  /**
   * @brief CPU Monte Carlo with both heatmap and density scoring.
   * @param step_size Maximum displacement per MC move.
   * @return Final energy score.
   */
  double MonteCarloHeatmapAndDensity(float step_size);

  /**
   * @brief Convert density-space coordinates to structure-space.
   * @param coord Density voxel coordinates.
   * @return Structure-space 3D position.
   */
  vector3 densityCoordToStructure(vector3 coord);

  /**
   * @brief Convert structure-space coordinates to density-space.
   * @param pos Structure-space 3D position.
   * @return Density voxel coordinates.
   */
  vector3 structureCoordToDensity(vector3 pos);

  /** @brief Set a SimulationLogger for energy tracing. */
  void setLogger(SimulationLogger *log) { sim_logger = log; }

  int total_mc_steps;                /**< Accumulated MC step count across all phases. */
  GpuResourceCache gpu_cache;        /**< Persistent GPU resources for arcs+smooth kernels. */

  MultiscaleEnergy energy;           /**< Two-scale energy decomposition tracker. */

  InteractionArcs arcs;              /**< All loaded interaction arcs data. */
  std::vector<Cluster> clusters;     /**< Global cluster array (all levels, all chromosomes). */
  std::vector<Cluster> base_clusters;/**< Backup of initial cluster state (for reset). */

  std::map<std::string, int> chr_root;  /**< Root cluster index per chromosome. */

  std::map<std::string, std::vector<int>> current_level;  /**< Active-level cluster indices per chromosome. */
  std::vector<int> active_region;       /**< Indices of clusters being reconstructed in the current MC run. */

  std::vector<Heatmap> heatmap;  /**< Singleton heatmaps: [0]=chromosome level, [1]=segment level. */
  Heatmap heatmap_dist;          /**< Expected-distance heatmap for the current level. */

  ChromosomesSet chr_set;         /**< Ensemble of reconstructed chromosome structures. */
  ChromosomesSet chromosome_set;  /**< Another ensemble storage (used during multi-run). */

  std::map<string, int> chr_first_cluster;  /**< Index of the first anchor cluster per chromosome. */

  std::map<string, int> chr_length;     /**< Chromosome lengths in bp (from anchor data). */
  std::vector<double> chr_length_v;     /**< Chromosome lengths as a vector (for sorting). */

  std::vector<string> chrs;  /**< List of chromosome names to reconstruct. */

  bool is_bed_region_provided;  /**< True if a specific region was requested with -c. */
  BedRegion selected_region;    /**< The requested BED region to reconstruct. */

private:
  /**
   * @brief Create singleton heatmap from background PET data.
   * @param diag Diagonal offset to start from.
   * @return The singleton interaction frequency heatmap.
   */
  Heatmap createSingletonHeatmap(int diag = 1);

  /**
   * @brief Create subanchor-level singleton heatmap for a chromosome.
   * @param chr Chromosome name.
   * @param anchors_gap_len Output: gap lengths between anchors.
   */
  void createSingletonSubanchorHeatmap(string chr,
                                       vector<int> &anchors_gap_len);

  /**
   * @brief Create expected-distance heatmap at subanchor level.
   * @param subanchor_avg_dist Average distance heatmap.
   */
  void createExpectedDistSubanchorHeatmap(const Heatmap &subanchor_avg_dist);

  /**
   * @brief Score the active region against the target heatmap.
   * @param moved Index of the moved cluster (-1 = score all).
   * @return Energy score (lower = better fit to heatmap).
   */
  double calcScoreHeatmapActiveRegion(int moved = -1);

  /**
   * @brief Score the current structure against the density map.
   * @return Density score.
   */
  double calcScoreDensity();

  /**
   * @brief Score arc-distance fit for the active region.
   * @return Distance score.
   */
  double calcScoreDistancesActiveRegion();

  /**
   * @brief Score arc-distance fit after moving a single cluster.
   * @param cluster_moved Index of the moved cluster.
   * @return Distance score.
   */
  double calcScoreDistancesActiveRegion(int cluster_moved);

  /**
   * @brief Score linker length consistency for the active region.
   * @return Length score.
   */
  double calcScoreStructureActiveRegionLengths();

  /**
   * @brief Score linker length after moving a single cluster.
   * @param cluster_moved Index of the moved cluster.
   * @return Length score.
   */
  double calcScoreStructureActiveRegionLengths(int cluster_moved);

  /**
   * @brief Score smoothness (linker lengths + chain angles) at subanchor level.
   * @param lengths Include length scoring.
   * @param angles Include angle scoring.
   * @return Smoothness score.
   */
  double calcScoreStructureSmooth(bool lengths, bool angles);

  /**
   * @brief Score smoothness after moving a single cluster.
   * @param cluster_moved Index of the moved cluster.
   * @param lengths Include length scoring.
   * @param angles Include angle scoring.
   * @return Smoothness score.
   */
  double calcScoreStructureSmooth(int cluster_moved, bool lengths, bool angles);

  /**
   * @brief Score CTCF motif orientation consistency.
   * @param orientation Vector of orientation unit vectors per anchor.
   * @return Orientation score.
   */
  double calcScoreOrientation(const vector<vector3> &orientation);

  /**
   * @brief Score orientation after moving a single anchor.
   * @param orientation Orientation vectors.
   * @param anchor_index Index of the moved anchor.
   * @return Orientation score.
   */
  double calcScoreOrientation(const vector<vector3> &orientation,
                              int anchor_index);

  /**
   * @brief Score subanchor heatmap fit.
   * @param cluster_moved Index of the moved cluster (-1 = score all).
   * @return Subanchor heatmap score.
   */
  double calcScoreSubanchorHeatmap(int cluster_moved = -1);

  /**
   * @brief Compute the orientation vector for a cluster.
   * @param cluster_index Index in active_region.
   * @return Unit vector representing the cluster's orientation.
   */
  vector3 calcOrientation(int cluster_index);

  /**
   * @brief GPU-accelerated MC for arc-level reconstruction.
   * @param step_size Maximum displacement.
   * @return Final energy score.
   */
  float parallelMonteCarloArcs(float step_size);

  /**
   * @brief CPU MC for anchor-level reconstruction (arcs + linker lengths).
   * @param step_size Maximum displacement.
   * @return Final energy score.
   */
  double MonteCarloArcs(float step_size);

  /**
   * @brief CPU MC for subanchor-level smoothing (lengths + angles).
   * @param step_size Maximum displacement.
   * @param use_subanchor_heatmap If true, include subanchor heatmap scoring.
   * @return Final energy score.
   */
  double MonteCarloArcsSmooth(float step_size,
                              bool use_subanchor_heatmap = false);

  /** @brief Update the density map from the current 3D structure. */
  void updateDensity();

  /** @brief Repair density at the structure boundary. */
  void repairDensityBoundary();

  /**
   * @brief Densify the active region by inserting subanchor beads.
   * @param cluster Parent cluster index.
   * @param fix If true, fix the new beads (no MC moves).
   */
  void densifyActiveRegion(int cluster, bool fix = false);

  /**
   * @brief Convert genomic distance (bp) to expected spatial distance.
   * @param dist Genomic distance in bp.
   * @return Expected spatial distance.
   */
  double genomicLengthToDistance(int dist);

  /**
   * @brief Downsample a segment-level heatmap to chromosome level.
   * @param heat Segment-level singleton heatmap.
   * @return Chromosome-level heatmap.
   */
  Heatmap downsampleSingletonHeatmap(Heatmap &heat);

  /**
   * @brief Create a heatmap of 3D distances between clusters.
   * @param clusters_ind Vector of cluster indices.
   * @return Distance heatmap.
   */
  Heatmap calculateRegionsDistances(const std::vector<int> clusters_ind);

  /**
   * @brief Find the other endpoint of an arc.
   * @param arc The interaction arc.
   * @param current_cluster One endpoint cluster index.
   * @return The other endpoint cluster index.
   */
  int otherEnd(const InteractionArc &arc, int current_cluster);

  /** @overload */
  int otherEnd(int arc, int current_cluster);

  /**
   * @brief Find the other endpoint of an arc (chromosome-specific).
   * @param chr Chromosome name.
   * @param arc Arc index.
   * @param current_cluster One endpoint.
   * @return The other endpoint.
   */
  int otherEnd(string chr, int arc, int current_cluster);

  /**
   * @brief Interpolate children positions from parent-level positions.
   * @param regions Vector of parent region indices.
   */
  void interpolateChildrenPosition(std::vector<int> &regions);

  /**
   * @brief Interpolate positions for a specific region's children.
   * @param regions Parent regions.
   * @param clusters_template Template cluster indices.
   * @param region_ind Region index.
   * @param noise If true, add random noise.
   */
  void interpolatePosition(std::vector<int> &regions,
                           std::vector<int> &clusters_template, int region_ind,
                           bool noise = false);

  /**
   * @brief Interpolate children positions using spline curves.
   * @param regions Parent region indices.
   * @param use_genomic_dist If true, use genomic distances for spacing.
   */
  void interpolateChildrenPositionSpline(std::vector<int> &regions,
                                         bool use_genomic_dist = false);

  /**
   * @brief Position interaction blocks within segments.
   * @param segments Vector of segment cluster indices.
   */
  void positionInteractionBlocks(std::vector<int> &segments);

  /**
   * @brief Compute boundary score for a point on a chromosome.
   * @param chr Chromosome name.
   * @param pt 3D point to evaluate.
   * @return Boundary penalty score.
   */
  float boundaryScore(string chr, vector3 pt);

  /**
   * @brief Map an active_region index to its chromosome name.
   * @param ind Active region index.
   * @return Chromosome name.
   */
  string activeRegionIndexToChr(int ind);

  /**
   * @brief Get the set of heatmap break indices (segment boundaries).
   * @return Set of break indices.
   */
  std::set<int> getCurrentHeatmapBreaks();

  /**
   * @brief Compute true 3D distances for a set of regions.
   * @param regions Vector of region indices.
   * @return Heatmap of true pairwise distances.
   */
  Heatmap calcTrueDistancesHeatmapForRegion(std::vector<int> &regions);

  /** @brief Compute neighbor lists for active anchors (for orientation scoring). */
  void calcActiveAnchorsNeighbors();

  /** @brief Compute expected distance heatmap at anchor level. */
  void calcAnchorExpectedDistancesHeatmap();

  /**
   * @brief Get chromosome boundaries within the heatmap.
   * @param p Position in the heatmap.
   * @param[out] start Start index of the chromosome block.
   * @param[out] end End index of the chromosome block.
   */
  void getChromosomeHeatmapBoundary(int p, int &start, int &end);

  string current_chr;  /**< Currently processed chromosome name. */
  vector3 rw_pos;      /**< Helper position for random walk initialization. */
  bool zero_distance;  /**< Flag for zero-distance edge case. */

  float freq_to_distance[101];  /**< Memoization table: PET count -> spatial distance. */
  int tmp_var;                  /**< Temporary variable for internal use. */

  std::vector<string> arcs_clusters;        /**< PET cluster file paths. */
  std::vector<string> arcs_singletons;      /**< Singleton file paths (intra-chr). */
  std::vector<string> arcs_singletons_inter;/**< Singleton file paths (inter-chr). */

  bool split_singleton_files_by_chr;  /**< Whether singleton files are per-chromosome. */

  std::map<string, int> active_region_to_chr;  /**< Maps max cluster ID to chromosome name (for chr lookup). */

  string output_dir;  /**< Output directory path. */
  string label;       /**< Label for output filenames. */

  SimulationLogger *sim_logger = nullptr;  /**< Logger for energy tracing (optional). */

  int debug_length_limit;       /**< Debug: max genomic length per chromosome. */
  int debug_chromosomes_limit;  /**< Debug: max number of chromosomes. */

  BedRegions centromeres;         /**< Centromere regions (excluded from reconstruction). */
  BedRegions interaction_blocks;  /**< Interaction block split regions. */
  BedRegions segments_predefined; /**< Predefined segment boundaries. */

  Heatmap heatmap_anchor;           /**< Singleton heatmap at anchor level. */
  Heatmap heatmap_exp_dist_anchor;  /**< Expected-distance heatmap at anchor level. */
  std::map<int, int> active_to_cluster_index;  /**< Maps active-region index to global cluster index. */

  Heatmap heatmap_subanchor;       /**< Singleton heatmap at subanchor level. */
  Heatmap heatmap_dist_subanchor;  /**< Distance heatmap at subanchor level. */

  std::map<int, vector<int>> active_anchors_neighbors;       /**< Neighbor anchor indices per active anchor. */
  std::map<int, vector<float>> active_anchors_neighbors_weight;  /**< Neighbor weights per active anchor. */

  vector<int> heatmap_chromosome_boundaries;  /**< Chromosome boundary indices in the heatmap. */

  Density density;       /**< Input 3D density map (from microscopy). */
  Density density_curr;  /**< Current structure-derived density map. */

  /** @name Loop Energy Constraints: E_loop = k * (r_ij - r0)^2 */
  ///@{
  /**
   * @brief Explicit loop constraints loaded from file.
   * Each entry: (global_cluster_i, global_cluster_j, stiffness, eq_distance).
   */
  struct LoopConstraint {
    int cluster_i;
    int cluster_j;
    float stiffness;
    float eq_distance;
  };
  std::vector<LoopConstraint> loop_constraints;

  /**
   * @brief Fast lookup: maps (active_region_i, active_region_j) -> constraint index.
   * Rebuilt each time active_region changes.
   */
  std::map<std::pair<int,int>, int> active_loop_map;

  /** @brief Rebuild active_loop_map for current active_region. */
  void rebuildActiveLoopMap();

  /** @brief Compute loop energy contribution for the full active region. */
  double calcScoreLoopEnergy() const;

  /** @brief Compute loop energy contribution for a single moved bead. */
  double calcScoreLoopEnergy(int cluster_moved) const;
  ///@}
};

#endif /* LOOPERSOLVER_H_ */
