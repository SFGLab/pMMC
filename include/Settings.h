/**
 * @file Settings.h
 * @brief Global configuration parameters loaded from an INI file.
 *
 * All members are static, making Settings a global singleton.
 * Parameters control every aspect of the reconstruction: data paths,
 * Monte Carlo temperatures, frequency-to-distance mappings, spring
 * constants, CUDA kernel configuration, and more.
 *
 * Typical usage:
 * @code
 *   Settings stg;
 *   stg.loadFromINI("settings.ini");
 *   // Access: Settings::segmentSize, Settings::maxTemp, etc.
 * @endcode
 */

#ifndef SETTINGS_H_
#define SETTINGS_H_

#include <INIReader.h>
#include <common.h>
#include <stdio.h>
#include <string>

/**
 * @class Settings
 * @brief Static global configuration class loaded from an INI file.
 *
 * All reconstruction parameters are stored as static members.
 * Call loadFromINI() once at startup to populate them.
 */
class Settings {

public:
  /** @brief Default constructor (calls init() if not yet initialized). */
  Settings();

  /** @brief Set all parameters to their default values. */
  void init();

  /**
   * @brief Print current settings to stdout.
   * @param level Verbosity level (0 = summary, higher = more detail).
   */
  void print(int level = 0);

  /**
   * @brief Load settings from an INI file.
   * @param ini_path Path to the INI configuration file.
   * @return True if the file was loaded successfully.
   */
  bool loadFromINI(std::string ini_path);

  static bool initialized;  /**< True after first initialization. */

  /** @name General Options */
  ///@{
  static bool use2D;            /**< If true, reconstruct in 2D instead of 3D. */
  static int debug;             /**< Debug verbosity level. */
  static bool randomWalk;       /**< Use random walk for initial structure. */
  static bool useInputCache;    /**< Cache loaded input data for faster restarts. */
  static size_t maxMemoryMB;    /**< Memory budget in MB (0 = unlimited). Set via --max-memory-mb or config. */
  static int outputLevel;       /**< Output verbosity (5 = milestone info only). */
  static int loopDensity;       /**< Loop density parameter for densification. */
  ///@}

  /** @name CUDA Configuration */
  ///@{
  static int milestoneFailsThreshold;  /**< Max consecutive MC failures before stopping. */
  static int cudaBlocksMultiplier;     /**< CUDA grid size multiplier. */
  static int cudaThreadsPerBlock;      /**< Threads per CUDA block (typically 128 or 256). */
  ///@}

  /** @name Data Paths */
  ///@{
  static std::string dataDirectory;         /**< Base directory for all data files. */
  static std::string dataAnchors;           /**< Anchors file path (relative to dataDirectory). */
  static std::string dataPetClusters;       /**< Comma-separated PET cluster file paths. */
  static std::string dataSingletons;        /**< Comma-separated singleton (PET 1-3) file paths. */
  static std::string dataSingletonsInter;   /**< Inter-chromosomal singleton file paths. */
  static std::string dataFactors;           /**< Comma-separated factor names (e.g. "CTCF,RNAPII"). */
  static bool dataSplitSingletonFilesByChr; /**< If true, singleton files are split per chromosome. */
  static std::string dataCentromeres;       /**< Centromere BED file path (relative to dataDirectory). */
  static std::string dataSegmentsSplit;     /**< Predefined segment split BED file path. */
  static std::string dataSegmentHeatmap;    /**< Pre-computed segment heatmap file (optional). */
  ///@}

  /** @name Structural Template */
  ///@{
  static std::string templateSegment;  /**< Path to a template segment structure (optional). */
  static float templateScale;          /**< Scale factor for the template structure. */
  ///@}

  /** @name Density Constraints (Microscopy) */
  ///@{
  static bool useDensity;          /**< Enable density-constrained reconstruction. */
  static std::string densityMapFile;  /**< Path to the 3D density map file. */
  static float densityScale;       /**< Scale factor for density voxels. */
  static float densityInfluence;   /**< Influence radius of density constraint. */
  static float densityWeight;      /**< Weight of density term in the energy function. */
  ///@}

  /** @name Telomere Positions */
  ///@{
  static bool useTelomerePositions;    /**< If true, fix telomere positions during reconstruction. */
  static vector3 telomere_1_position;  /**< 3D position of the p-arm telomere. */
  static vector3 telomere_2_position;  /**< 3D position of the q-arm telomere. */
  ///@}

  static float ibRandomWalkJumps;  /**< Number of random walk jumps for interaction block initialization. */

  /** @name Distance Matrix (Segment Level) */
  ///@{
  static std::string distHeatmap;     /**< Path to a precomputed distance heatmap (segment level). */
  static float distHeatmapScale;      /**< Scale factor for the distance heatmap. */
  ///@}

  /** @name Rewiring */
  ///@{
  static bool rewiringEnabled;           /**< Enable arc rewiring based on CTCF motif orientation. */
  static int rewiringCertaintyThreshold; /**< Minimum PET count to keep an arc during rewiring. */
  ///@}

  /** @name CTCF Motif Orientation */
  ///@{
  static bool useCTCFMotifOrientation;  /**< Use CTCF convergent orientation as a scoring term. */
  static bool motifsSymmetric;          /**< Treat L-R and R-L orientations as equivalent. */
  static float motifOrientationWeigth;  /**< Weight of the orientation term in the energy function. */
  ///@}

  /** @name Anchor-Level Heatmap */
  ///@{
  static bool useAnchorHeatmap;         /**< Use singleton heatmap at anchor level. */
  static float anchorHeatmapInfluence;  /**< Influence weight for anchor heatmap. */
  static float anchorHeatmapDistWeight; /**< Distance weight for anchor heatmap scoring. */
  ///@}

  /** @name Subanchor-Level Heatmap */
  ///@{
  static bool useSubanchorHeatmap;              /**< Use singleton heatmap at subanchor level. */
  static float subanchorHeatmapInfluence;       /**< Influence weight for subanchor heatmap. */
  static float subanchorHeatmapDistWeight;      /**< Distance weight for subanchor heatmap. */
  static int subanchorEstimateDistancesReplicates;  /**< Replicates for distance estimation at subanchor level. */
  static int subanchorEstimateDistancesSteps;       /**< MC steps for distance estimation at subanchor level. */
  static float subanchorLoopLoopInteractionsScale;  /**< Scale for loop-loop interactions at subanchor level. */
  ///@}

  /** @name PET Cluster Filtering */
  ///@{
  static int maxPETClusterLength;           /**< Maximum arc length in bp (longer arcs go to long_arcs). */
  static float longPETClustersEffectPower;  /**< Power-law exponent for long PET cluster effect. */
  static float longPETClustersEffectScale;  /**< Scale factor for long PET cluster effect. */
  ///@}

  static int segmentSize;  /**< Target segment size in bp (~2Mb). */

  /** @name Simulation Steps Per Level */
  ///@{
  static int simulationStepsLevelChr;        /**< MC steps at chromosome level. */
  static int simulationStepsLevelSegment;    /**< MC steps at segment level. */
  static int simulationStepsLevelAnchor;     /**< MC steps at anchor level. */
  static int simulationStepsLevelSubanchor;  /**< MC steps at subanchor level. */
  ///@}

  /** @name Noise Coefficients Per Level */
  ///@{
  static float noiseCoefficientLevelChr;        /**< Displacement noise at chromosome level. */
  static float noiseCoefficientLevelSegment;    /**< Displacement noise at segment level. */
  static float noiseCoefficientLevelAnchor;     /**< Displacement noise at anchor level. */
  static float noiseCoefficientLevelSubanchor;  /**< Displacement noise at subanchor level. */
  ///@}

  /** @name Heatmap Parameters */
  ///@{
  static float heatmapInterScaling;                /**< Scaling factor for inter-chromosomal heatmap. */
  static float heatmapDistanceHeatmapStretching;   /**< Max stretching factor for distance heatmap (multiples of average). */
  ///@}

  /** @name Frequency-to-Distance Mapping (Heatmap Level) */
  ///@{
  static float freqToDistHeatmapScale;       /**< Scale: d = scale * f^(-power). */
  static float freqToDistHeatmapPower;       /**< Power exponent for frequency-to-distance. */
  static float freqToDistHeatmapScaleInter;  /**< Scale for inter-chromosomal frequency-to-distance. */
  static float freqToDistHeatmapPowerInter;  /**< Power for inter-chromosomal frequency-to-distance. */
  ///@}

  /** @name Count-to-Distance Mapping (Arc Level) */
  ///@{
  static float countToDistA;          /**< 'a' parameter: d = a * (count + shift)^(-scale) + baseLevel. */
  static float countToDistScale;      /**< Scale exponent for count-to-distance. */
  static float countToDistShift;      /**< Additive shift to count before conversion. */
  static float countToDistBaseLevel;  /**< Base distance level (minimum distance). */
  ///@}

  /** @name Genomic-Distance-to-Spatial-Distance Mapping */
  ///@{
  static float genomicLengthToDistPower;  /**< Power: d = scale * genomic_length^power + base. */
  static float genomicLengthToDistScale;  /**< Scale factor. */
  static float genomicLengthToDistBase;   /**< Base distance. */
  ///@}

  /** @name Spring Constants (Linker Length Scoring) */
  ///@{
  static float springConstantStretch;      /**< Stretch penalty for linkers (segment/anchor). */
  static float springConstantSqueeze;      /**< Squeeze penalty for linkers. */
  static float springAngularConstant;      /**< Angular penalty for chain bending. */
  static float springConstantStretchArcs;  /**< Stretch penalty for arc-connected anchors. */
  static float springConstantSqueezeArcs;  /**< Squeeze penalty for arc-connected anchors. */
  ///@}

  /** @name Loop Energy Term: E_loop = k_loop * (r_ij - r0)^2 */
  ///@{
  static float loopStiffness;             /**< Spring constant for explicit loop constraints. */
  static float loopEquilibriumDistance;    /**< Equilibrium distance r0 for loop constraints. */
  static std::string loopConstraintsFile; /**< Path to external loop constraints file. */
  ///@}

  /** @name Merged Input */
  ///@{
  static std::string mergedInputFile;    /**< Path to .merged file (merged BED+BEDPE). */
  ///@}

  /** @name MC Simulated Annealing: Heatmap Level (Segment/Chromosome) */
  ///@{
  static float maxTempHeatmap;                       /**< Starting temperature. */
  static float dtTempHeatmap;                        /**< Temperature decrement per step. */
  static float tempJumpScaleHeatmap;                 /**< Scale of temperature jumps on improvement. */
  static float tempJumpCoefHeatmap;                  /**< Coefficient for adaptive temperature jumps. */
  static float MCstopConditionImprovementHeatmap;    /**< Minimum improvement to continue. */
  static int MCstopConditionMinSuccessesHeatmap;     /**< Minimum accepted moves per milestone. */
  static int MCstopConditionStepsHeatmap;            /**< Steps per milestone. */
  ///@}

  /** @name MC Simulated Annealing: Heatmap + Density */
  ///@{
  static float maxTempHeatmapDensity;                       /**< Starting temperature (with density). */
  static float dtTempHeatmapDensity;                        /**< Temperature decrement. */
  static float tempJumpScaleHeatmapDensity;                 /**< Jump scale on improvement. */
  static float tempJumpCoefHeatmapDensity;                  /**< Adaptive jump coefficient. */
  static float MCstopConditionImprovementHeatmapDensity;    /**< Min improvement threshold. */
  static int MCstopConditionMinSuccessesHeatmapDensity;     /**< Min accepted moves. */
  static int MCstopConditionStepsHeatmapDensity;            /**< Steps per milestone. */
  ///@}

  /** @name MC Simulated Annealing: Arc Level (Anchor) */
  ///@{
  static float maxTemp;                       /**< Starting temperature. */
  static float dtTemp;                        /**< Temperature decrement. */
  static float tempJumpScale;                 /**< Jump scale on improvement. */
  static float tempJumpCoef;                  /**< Adaptive jump coefficient. */
  static float MCstopConditionImprovement;    /**< Min improvement threshold. */
  static int MCstopConditionMinSuccesses;     /**< Min accepted moves per milestone. */
  static int MCstopConditionSteps;            /**< Steps per milestone. */
  ///@}

  /** @name MC Simulated Annealing: Smooth Level (Subanchor) */
  ///@{
  static float weightDistSmooth;                     /**< Weight for distance term during smoothing. */
  static float weightAngleSmooth;                    /**< Weight for angle term during smoothing. */
  static float maxTempSmooth;                        /**< Starting temperature. */
  static float dtTempSmooth;                         /**< Temperature decrement. */
  static float tempJumpScaleSmooth;                  /**< Jump scale on improvement. */
  static float tempJumpCoefSmooth;                   /**< Adaptive jump coefficient. */
  static float MCstopConditionImprovementSmooth;     /**< Min improvement threshold. */
  static int MCstopConditionMinSuccessesSmooth;      /**< Min accepted moves. */
  static int MCstopConditionStepsSmooth;             /**< Steps per milestone. */
  ///@}

  /** @name Energy Tracing */
  ///@{
  static bool energyTraceEnabled;    /**< Enable per-step energy CSV logging (-E flag). */
  static int energyTraceInterval;    /**< Log energy every N MC steps (-L flag, default 1000). */
  ///@}

  /** @name Initialization Method */
  ///@{
  static std::string initMethod;     /**< Initialization method: "random" (default) or "mds". */
  ///@}

  /** @name GPU Seed */
  ///@{
  static unsigned int gpuSeed;       /**< GPU RNG seed (from CLI -j). */
  ///@}
};

#endif /* SETTINGS_H_ */
