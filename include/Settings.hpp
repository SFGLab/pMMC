/**
 * @file Settings.hpp
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
  static size_t maxMemoryMB;    /**< Memory budget in MB (0 = unlimited). */
  static int outputLevel;       /**< Output verbosity (5 = milestone info only). */
  static int loopDensity;       /**< Loop density parameter for densification. */
  ///@}

  /** @name Run Configuration */
  ///@{
  static std::string chromosomes;    /**< Chromosomes to reconstruct (e.g. "genome", "chr14", "chr1-chr5"). */
  static std::string outputDir;      /**< Output directory (default: "./"). */
  static std::string label;          /**< Output label for file names (default: "run"). */
  static int seed;                   /**< Random seed (0 = time-based). */
  static int ensembleSize;           /**< Number of structures to generate (default: 1). */
  static std::string outputFormat;   /**< Additional output format: pdb, cif, both, or empty. */
  static std::string annotationDir;  /**< Annotation directory (enhancer.bed, promoter.bed, epigenetic.bed). */
  static int maxLevel;               /**< Maximum reconstruction level (default: 5). */
  static int lengthLimit;            /**< Length limit in bp (-1 = no limit). */
  static int chrNumberLimit;         /**< Max chromosomes to reconstruct (-1 = no limit). */
  static bool useNewFileFormat;      /**< Use new HCM file format. */
  static std::string selectedFactors; /**< Comma-separated factor subset to use (empty = all). */
  static std::string cellLineList;   /**< Path to cell-line list file for batch mode (empty = single run). */
  ///@}

  /** @name Analysis Configuration */
  ///@{
  static std::string action;         /**< Action/command to perform. */
  static std::string inputFile2;     /**< Second input file (for metrics comparison). */
  static std::string bedRegionsFile; /**< BED regions file (for position extraction). */
  static std::string pattern;        /**< File pattern with {N} placeholder (for ensemble/rewiring). */
  static int level;                  /**< Hierarchy level (0=chr, 1=segment, 2=subanchor, 3=bead). */
  static int resolution;             /**< Resolution in bp for resampling. */
  static bool addHeader;             /**< Add header line to output (for position). */
  static int extractStart;           /**< Start position for extract command. */
  static int extractEnd;             /**< End position for extract command. */
  static std::string bedFilePath;    /**< BED file path (for merge). */
  static std::string bedpeFilePath;  /**< BEDPE file path (for merge). */
  static int numLoops;               /**< Number of loops for synthetic data generation. */
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

  /** @name Input Format */
  ///@{
  static std::string inputFormat;        /**< Input format: chiapet (default), merged, ibed, contact. */
  static std::string inputFile;          /**< Path to input file (for merged/ibed/contact formats). */
  static int contactResolution;          /**< Resolution in bp for contact matrix format (default: 5000). */
  static std::string mergedInputFile;    /**< Path to .merged file (merged BED+BEDPE). [deprecated: use inputFormat/inputFile] */
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

  /** @name Spherical Boundary */
  ///@{
  static float sphere_radius;  /**< Radius of confining sphere (0 = disabled). */
  ///@}

  /** @name Annotation Files (BED) */
  ///@{
  static std::string enhancerAnnotationFile;   /**< Path to enhancer BED annotation file. */
  static std::string promoterAnnotationFile;   /**< Path to promoter BED annotation file. */
  static std::string epigeneticAnnotationFile; /**< Path to epigenetic BED annotation file. */
  ///@}
};


// ============================================================================
// Implementation
// ============================================================================


inline Settings::Settings() {
  if (!initialized)
    init();
}

inline void Settings::init() {
  initialized = true;
  debug = 0;
  outputLevel = 0;
  randomWalk = false;

  use2D = false;

  useInputCache = true;
  maxMemoryMB = 0; // 0 = unlimited

  // Run configuration
  chromosomes = "genome";
  outputDir = "./";
  label = "run";
  seed = 0;
  ensembleSize = 1;
  outputFormat = "";
  annotationDir = "";
  maxLevel = 5;
  lengthLimit = -1;
  chrNumberLimit = -1;
  useNewFileFormat = false;
  selectedFactors = "";
  cellLineList = "";

  // Analysis configuration
  action = "";
  inputFile2 = "";
  bedRegionsFile = "";
  pattern = "";
  level = 3;
  resolution = 25000;
  addHeader = false;
  extractStart = 0;
  extractEnd = 0;
  bedFilePath = "";
  bedpeFilePath = "";
  numLoops = 100;

  energyTraceEnabled = true;
  energyTraceInterval = 1000;
  initMethod = "random";
  gpuSeed = 0;

  loopDensity = 5;

  milestoneFailsThreshold = 3;
  cudaBlocksMultiplier = 4;
  cudaThreadsPerBlock = 256;

  useCTCFMotifOrientation = false;
  motifsSymmetric = true;
  motifOrientationWeigth = 1.0f;

  useSubanchorHeatmap = false;
  subanchorEstimateDistancesReplicates = 5;
  subanchorEstimateDistancesSteps = 2;
  subanchorHeatmapInfluence = 0.5f;
  subanchorHeatmapDistWeight = 1.0f;

  useAnchorHeatmap = false;
  anchorHeatmapInfluence = 0.5f;
  anchorHeatmapDistWeight = 1.0f;

  segmentSize = 2000000;
  maxPETClusterLength = 1000000;
  longPETClustersEffectPower = 2.0f;
  longPETClustersEffectScale = 10.0f;

  rewiringCertaintyThreshold = 16;

  simulationStepsLevelChr = 2;
  simulationStepsLevelSegment = 2;
  simulationStepsLevelAnchor = 5;
  simulationStepsLevelSubanchor = 5;

  templateScale = 1.0;
  distHeatmapScale = 1.0;

  useDensity = false;
  densityScale = 1.0f;
  densityInfluence = 0.95f;
  densityWeight = 1.0f;

  ibRandomWalkJumps = 10.0f;

  noiseCoefficientLevelChr = 1.0f;
  noiseCoefficientLevelSegment = 0.1f;
  noiseCoefficientLevelAnchor = 0.5f;
  noiseCoefficientLevelSubanchor = 0.5f;

  heatmapInterScaling = 1.0f;
  heatmapDistanceHeatmapStretching = 2.0f;

  freqToDistHeatmapScale = 100.0f;
  freqToDistHeatmapPower = -0.333f;
  freqToDistHeatmapScaleInter = 100.0f;
  freqToDistHeatmapPowerInter = -1.0;

  countToDistA = 0.5f;
  countToDistScale = 20.0f;
  countToDistShift = 1.0f;
  countToDistBaseLevel = 0.01f;

  genomicLengthToDistPower = 0.5f;
  genomicLengthToDistScale = 1.0f;
  genomicLengthToDistBase = 0.0f;

  maxTempHeatmap = 20.0f;
  dtTempHeatmap = 0.99995f;
  tempJumpCoefHeatmap = 20.0f;
  tempJumpScaleHeatmap = 50.0f;
  MCstopConditionImprovementHeatmap = 0.995f;
  MCstopConditionMinSuccessesHeatmap = 5;
  MCstopConditionStepsHeatmap = 10000;

  maxTempHeatmapDensity = 20.0f;
  dtTempHeatmapDensity = 0.99995f;
  tempJumpCoefHeatmapDensity = 20.0f;
  tempJumpScaleHeatmapDensity = 50.0f;
  MCstopConditionImprovementHeatmapDensity = 0.995f;
  MCstopConditionMinSuccessesHeatmapDensity = 5;
  MCstopConditionStepsHeatmapDensity = 10000;

  maxTemp = 20.0f;
  dtTemp = 0.99995f;
  tempJumpCoef = 20.0f;
  tempJumpScale = 50.0f;
  MCstopConditionImprovement = 0.995f;
  MCstopConditionMinSuccesses = 5;
  MCstopConditionSteps = 10000;

  weightAngleSmooth = 1.0f;
  weightDistSmooth = 1.0f;
  maxTempSmooth = 20.0f;
  dtTempSmooth = 0.99995f;
  tempJumpCoefSmooth = 20.0f;
  tempJumpScaleSmooth = 50.0f;
  MCstopConditionImprovementSmooth = 0.995f;
  MCstopConditionMinSuccessesSmooth = 5;
  MCstopConditionStepsSmooth = 10000;

  springConstantSqueeze = 0.1f;
  springConstantStretch = 0.1f;
  springAngularConstant = 0.1f;

  springConstantSqueezeArcs = 1.0f;
  springConstantStretchArcs = 1.0f;

  loopStiffness = 1.0f;
  loopEquilibriumDistance = 1.5f;
  loopConstraintsFile = "";
  mergedInputFile = "";
  inputFormat = "";
  inputFile = "";
  contactResolution = 5000;

  sphere_radius = 0.0f; // 0 = spherical boundary disabled

  enhancerAnnotationFile = "";
  promoterAnnotationFile = "";
  epigeneticAnnotationFile = "";
}

inline void Settings::print(int level) {

  // printf("debug: %d\n", debug);
  printf("cache input: %s\n", useInputCache ? "yes" : "no");
  printf("random walk: %s\n", randomWalk ? "yes" : "no");
  printf("2D: %s\n", use2D ? "yes" : "no");

  printf("use motif orientation: %s (symmetric: %s)\n",
         useCTCFMotifOrientation ? "yes" : "no",
         motifsSymmetric ? "yes" : "no");
  printf("use singletons on subanchor level: %s\n",
         useSubanchorHeatmap ? "yes" : "no");

  printf("loop density: %d\n", loopDensity);

  printf("use density: %s\n", useDensity ? "yes" : "no");
  if (useDensity) {
    printf("   density map file: %s\n", densityMapFile.c_str());
    printf("   density scale: %f\n", densityScale);
    printf("   density influence: %f\n", densityInfluence);
    printf("   density weight: %f\n", densityWeight);
  }

  printf("use rewiring: %s\n", rewiringEnabled ? "yes" : "no");
  if (rewiringEnabled) {
    printf("   certainty threshold: %d\n", rewiringCertaintyThreshold);
  }

  printf("use telomere positions: %s\n", useTelomerePositions ? "yes" : "no");
  if (useTelomerePositions) {
    print_vector(telomere_1_position, "   telomere 1 position");
    print_vector(telomere_2_position, "   telomere 2 position");
  }

  printf("data\n");
  printf("   data directory: %s\n", dataDirectory.c_str());
  printf("   anchors: %s\n", dataAnchors.c_str());
  printf("   PET clusters: %s\n", dataPetClusters.c_str());
  printf("   singletons: %s\n", dataSingletons.c_str());
  printf("   split singleton files by chr: %d\n", dataSplitSingletonFilesByChr);
  printf("   singletons (inter): %s\n", dataSingletonsInter.c_str());
  printf("   factors: %s\n", dataFactors.c_str());
  printf("   centromeres: %s\n", dataCentromeres.c_str());
  printf("   segment split: %s\n", dataSegmentsSplit.c_str());
  printf("   segment heatmap: %s\n", dataSegmentHeatmap.c_str());

  printf("template structure: [%s]\n", templateSegment.c_str());
  printf("template scale: %f\n", templateScale);

  // printf("segment size: %d\n", segmentSize);

  printf("simulation steps, lvl1: %d\n", simulationStepsLevelChr);
  printf("simulation steps, lvl2: %d\n", simulationStepsLevelSegment);
  printf("simulation steps, arcs: %d\n", simulationStepsLevelAnchor);
  printf("simulation steps, smooth: %d\n", simulationStepsLevelSubanchor);

  //	printf("noise coef, lvl1: %f\n", noiseCoefficientLevel1);
  //	printf("noise coef, lvl2: %f\n", noiseCoefficientLevel2);
  //	printf("noise coef, arcs: %f\n", noiseCoefficientLevelArcs);
  //	printf("noise coef, smooth: %f\n", noiseCoefficientLevelArcsSmooth);

  printf("heatmap distance stretching: %f\n", heatmapDistanceHeatmapStretching);

  printf("ib random walk jump size: %f\n", ibRandomWalkJumps);
  printf("frequency-distance scale: %f\n", freqToDistHeatmapScale);
  printf("frequency-distance power: %f\n", freqToDistHeatmapPower);
  printf("frequency-distance power (inter): %f\n", freqToDistHeatmapPowerInter);
  printf("count-distance A: %f\n", countToDistA);
  printf("count-distance scale: %f\n", countToDistScale);
  printf("count-distance shift: %f\n", countToDistShift);
  printf("count-distance min dist: %f\n", countToDistBaseLevel);
  printf("genomic length (kb) to dist power: %f\n", genomicLengthToDistPower);
  printf("genomic length (kb) to dist scale: %f\n", genomicLengthToDistScale);
  printf("genomic length (kb) to dist base: %f\n", genomicLengthToDistBase);

  printf("spring stretch: %f\n", springConstantStretch);
  printf("spring squeeze: %f\n", springConstantSqueeze);
  printf("spring angular: %f\n", springAngularConstant);
  printf("spring stretch (arcs): %f\n", springConstantStretchArcs);
  printf("spring squeeze (arcs): %f\n", springConstantSqueezeArcs);

  printf("maxTemp: %f\n", maxTemp);
  printf("dTemp: %f\n", dtTemp);
  printf("jump coef: %f\n", tempJumpCoef);
  printf("jump scale: %f\n", tempJumpScale);
  printf("stop condition steps: %d\n", MCstopConditionSteps);
  printf("stop condition improvement threshold: %f\n",
         MCstopConditionImprovement);
  printf("stop condition successes threshold: %d\n",
         MCstopConditionMinSuccesses);

  printf("maxTemp (smooth): %f\n", maxTempSmooth);
  printf("dTemp: %f\n", dtTempSmooth);
  printf("jump coef: %f\n", tempJumpCoefSmooth);
  printf("jump scale: %f\n", tempJumpScaleSmooth);
  printf("stop condition steps: %d\n", MCstopConditionStepsSmooth);
  printf("stop condition improvement threshold: %f\n",
         MCstopConditionImprovementSmooth);
  printf("stop condition successes threshold: %d\n",
         MCstopConditionMinSuccessesSmooth);

  printf("sphere radius: %f%s\n", sphere_radius,
         sphere_radius > 0.0f ? "" : " (disabled)");

  if (!enhancerAnnotationFile.empty())
    printf("enhancer annotation: %s\n", enhancerAnnotationFile.c_str());
  if (!promoterAnnotationFile.empty())
    printf("promoter annotation: %s\n", promoterAnnotationFile.c_str());
  if (!epigeneticAnnotationFile.empty())
    printf("epigenetic annotation: %s\n", epigeneticAnnotationFile.c_str());

  printf("maxTemp (heatmap): %f\n", maxTempHeatmap);
  printf("dTemp (heatmap): %f\n", dtTempHeatmap);
  printf("jump coef: %f\n", tempJumpCoefHeatmap);
  printf("jump scale: %f\n", tempJumpScaleHeatmap);
  printf("stop condition steps (heatmap): %d\n", MCstopConditionStepsHeatmap);
  printf("stop condition improvement threshold (heatmap): %f\n",
         MCstopConditionImprovementHeatmap);
  printf("stop condition successes threshold (heatmap): %d\n",
         MCstopConditionMinSuccessesHeatmap);
}

inline bool Settings::loadFromINI(std::string ini_path) {
  INIReader reader(ini_path);

  if (reader.ParseError() < 0) {
    printf("Can't load [%s]'\n", ini_path.c_str());
    return false;
  }

  cudaBlocksMultiplier =
      reader.GetInteger("cuda", "blocks_multiplier", cudaBlocksMultiplier);
  cudaThreadsPerBlock =
      reader.GetInteger("cuda", "num_threads", cudaThreadsPerBlock);
  milestoneFailsThreshold =
      reader.GetInteger("cuda", "milestone_fails", milestoneFailsThreshold);

  // debug = reader.GetInteger("main", "debug", debug);
  outputLevel = reader.GetInteger("main", "output_level", outputLevel);
  randomWalk = reader.GetBoolean("main", "random_walk", randomWalk);

  useInputCache = reader.GetBoolean("main", "cache_input", useInputCache);
  maxMemoryMB = (size_t)reader.GetInteger("main", "max_memory_mb", (long)maxMemoryMB);
  use2D = reader.GetBoolean("main", "use_2D", use2D);

  // Run configuration
  chromosomes = reader.Get("main", "chromosomes", chromosomes);
  outputDir = reader.Get("main", "output_dir", outputDir);
  label = reader.Get("main", "label", label);
  seed = reader.GetInteger("main", "seed", seed);
  ensembleSize = reader.GetInteger("main", "ensemble_size", ensembleSize);
  outputFormat = reader.Get("main", "output_format", outputFormat);
  annotationDir = reader.Get("main", "annotation_dir", annotationDir);
  maxLevel = reader.GetInteger("main", "max_level", maxLevel);
  lengthLimit = reader.GetInteger("main", "length_limit", lengthLimit);
  chrNumberLimit = reader.GetInteger("main", "chr_number_limit", chrNumberLimit);
  useNewFileFormat = reader.GetBoolean("main", "use_new_file_format", useNewFileFormat);
  selectedFactors = reader.Get("data", "selected_factors", selectedFactors);
  cellLineList = reader.Get("data", "cell_line_list", cellLineList);

  // Analysis configuration
  action = reader.Get("main", "action", action);
  inputFile2 = reader.Get("main", "input_file_2", inputFile2);
  bedRegionsFile = reader.Get("main", "bed_regions_file", bedRegionsFile);
  pattern = reader.Get("main", "pattern", pattern);
  level = reader.GetInteger("main", "level", level);
  resolution = reader.GetInteger("main", "resolution", resolution);
  addHeader = reader.GetBoolean("main", "add_header", addHeader);
  extractStart = reader.GetInteger("main", "extract_start", extractStart);
  extractEnd = reader.GetInteger("main", "extract_end", extractEnd);
  bedFilePath = reader.Get("data", "bed_file", bedFilePath);
  bedpeFilePath = reader.Get("data", "bedpe_file", bedpeFilePath);
  numLoops = reader.GetInteger("main", "num_loops", numLoops);

  loopDensity = reader.GetInteger("main", "loop_density", loopDensity);

  useCTCFMotifOrientation = reader.GetBoolean(
      "motif_orientation", "use_motif_orientation", useCTCFMotifOrientation);
  motifsSymmetric = reader.GetBoolean("motif_orientation", "symmetric_motifs",
                                      motifsSymmetric);
  motifOrientationWeigth = (float)reader.GetReal("motif_orientation", "weight",
                                                 motifOrientationWeigth);

  useSubanchorHeatmap = reader.GetBoolean(
      "subanchor_heatmap", "use_subanchor_heatmap", useSubanchorHeatmap);
  subanchorEstimateDistancesReplicates =
      reader.GetInteger("subanchor_heatmap", "estimate_distances_replicates",
                        subanchorEstimateDistancesReplicates);
  subanchorEstimateDistancesSteps =
      reader.GetInteger("subanchor_heatmap", "estimate_distances_steps",
                        subanchorEstimateDistancesSteps);
  subanchorHeatmapInfluence = (float)reader.GetReal(
      "subanchor_heatmap", "heatmap_influence", subanchorHeatmapInfluence);
  subanchorHeatmapDistWeight = (float)reader.GetReal(
      "subanchor_heatmap", "heatmap_dist_weight", subanchorHeatmapDistWeight);

  useAnchorHeatmap = reader.GetBoolean("anchor_heatmap", "use_anchor_heatmap",
                                       useAnchorHeatmap);
  anchorHeatmapInfluence = (float)reader.GetReal(
      "anchor_heatmap", "heatmap_influence", anchorHeatmapInfluence);
  // anchorHeatmapDistWeight = (float)reader.GetReal("anchor_heatmap",
  // "heatmap_dist_weight", anchorHeatmapDistWeight);

  // segmentSize = reader.GetInteger("main", "segment_size", segmentSize);

  maxPETClusterLength =
      reader.GetInteger("main", "max_pet_length", maxPETClusterLength);
  longPETClustersEffectPower = (float)reader.GetReal(
      "main", "long_pet_power", longPETClustersEffectPower);
  longPETClustersEffectScale = (float)reader.GetReal(
      "main", "long_pet_scale", longPETClustersEffectScale);

  rewiringEnabled =
      reader.GetBoolean("rewiring", "rewiring_enabled", rewiringEnabled);
  rewiringCertaintyThreshold = reader.GetInteger(
      "rewiring", "rewiring_certainty_threshold", rewiringCertaintyThreshold);

  dataDirectory = reader.Get("data", "data_dir", dataDirectory);
  dataAnchors = reader.Get("data", "anchors", dataAnchors);
  dataPetClusters = reader.Get("data", "clusters", dataPetClusters);
  dataSingletons = reader.Get("data", "singletons", dataSingletons);
  dataSingletonsInter =
      reader.Get("data", "singletons_inter", dataSingletonsInter);
  dataSplitSingletonFilesByChr = reader.GetBoolean(
      "data", "split_singleton_files_by_chr", dataSplitSingletonFilesByChr);
  dataFactors = reader.Get("data", "factors", dataFactors);

  dataCentromeres = reader.Get("data", "centromeres", dataCentromeres);
  dataSegmentsSplit = reader.Get("data", "segment_split", dataSegmentsSplit);
  mergedInputFile = reader.Get("data", "merged_input", mergedInputFile);
  inputFormat = reader.Get("data", "input_format", inputFormat);
  inputFile = reader.Get("data", "input_file", inputFile);
  contactResolution = reader.GetInteger("data", "contact_resolution", contactResolution);

  dataSegmentHeatmap =
      reader.Get("template", "segment_heatmap", dataSegmentHeatmap);

  templateSegment = reader.Get("template", "template_segment", templateSegment);
  templateScale =
      (float)reader.GetReal("template", "template_scale", templateScale);

  distHeatmap = reader.Get("template", "dist_heatmap", distHeatmap);
  distHeatmapScale =
      (float)reader.GetReal("template", "dist_heatmap_scale", distHeatmapScale);

  useDensity = reader.GetBoolean("density", "use_density", useDensity);
  densityMapFile = reader.Get("density", "density_data", densityMapFile);
  densityScale =
      (float)reader.GetReal("density", "density_scale", densityScale);
  densityInfluence =
      (float)reader.GetReal("density", "density_influence", densityInfluence);
  densityWeight =
      (float)reader.GetReal("density", "density_weight", densityWeight);

  useTelomerePositions =
      reader.GetBoolean("telomeres", "use_telomeres", useTelomerePositions);
  if (useTelomerePositions) {
    string p = reader.Get("telomeres", "telomere_1_position", "");
    sscanf(p.c_str(), "%f,%f,%f", &telomere_1_position.x,
           &telomere_1_position.y, &telomere_1_position.z);
    p = reader.Get("telomeres", "telomere_2_position", "");
    sscanf(p.c_str(), "%f,%f,%f", &telomere_2_position.x,
           &telomere_2_position.y, &telomere_2_position.z);
  }

  simulationStepsLevelChr =
      reader.GetInteger("main", "steps_lvl1", simulationStepsLevelChr);
  simulationStepsLevelSegment =
      reader.GetInteger("main", "steps_lvl2", simulationStepsLevelSegment);
  simulationStepsLevelAnchor =
      reader.GetInteger("main", "steps_arcs", simulationStepsLevelAnchor);
  simulationStepsLevelSubanchor =
      reader.GetInteger("main", "steps_smooth", simulationStepsLevelSubanchor);

  noiseCoefficientLevelChr =
      (float)reader.GetReal("main", "noise_lvl1", noiseCoefficientLevelChr);
  noiseCoefficientLevelSegment =
      (float)reader.GetReal("main", "noise_lvl2", noiseCoefficientLevelSegment);
  noiseCoefficientLevelAnchor =
      (float)reader.GetReal("main", "noise_arcs", noiseCoefficientLevelAnchor);
  noiseCoefficientLevelSubanchor = (float)reader.GetReal(
      "main", "noise_smooth", noiseCoefficientLevelSubanchor);

  heatmapInterScaling =
      (float)reader.GetReal("heatmaps", "inter_scaling", heatmapInterScaling);
  heatmapDistanceHeatmapStretching =
      (float)reader.GetReal("heatmaps", "distance_heatmap_stretching",
                            heatmapDistanceHeatmapStretching);

  ibRandomWalkJumps = (float)reader.GetReal("distance", "ib_random_walk_jumps",
                                            ibRandomWalkJumps);
  freqToDistHeatmapScale = (float)reader.GetReal("distance", "freq_dist_scale",
                                                 freqToDistHeatmapScale);
  freqToDistHeatmapPower = (float)reader.GetReal("distance", "freq_dist_power",
                                                 freqToDistHeatmapPower);
  freqToDistHeatmapScaleInter = (float)reader.GetReal(
      "distance", "freq_dist_scale_inter", freqToDistHeatmapScaleInter);
  freqToDistHeatmapPowerInter = (float)reader.GetReal(
      "distance", "freq_dist_power_inter", freqToDistHeatmapPowerInter);
  countToDistA =
      (float)reader.GetReal("distance", "count_dist_a", countToDistA);
  countToDistScale =
      (float)reader.GetReal("distance", "count_dist_scale", countToDistScale);
  countToDistShift =
      (float)reader.GetReal("distance", "count_dist_shift", countToDistShift);
  countToDistBaseLevel = (float)reader.GetReal(
      "distance", "count_dist_base_level", countToDistBaseLevel);
  genomicLengthToDistPower = (float)reader.GetReal(
      "distance", "genomic_dist_power", genomicLengthToDistPower);
  genomicLengthToDistScale = (float)reader.GetReal(
      "distance", "genomic_dist_scale", genomicLengthToDistScale);
  genomicLengthToDistBase = (float)reader.GetReal(
      "distance", "genomic_dist_base", genomicLengthToDistBase);

  springConstantStretch = (float)reader.GetReal("springs", "stretch_constant",
                                                springConstantStretch);
  springConstantSqueeze = (float)reader.GetReal("springs", "squeeze_constant",
                                                springConstantSqueeze);
  springAngularConstant = (float)reader.GetReal("springs", "angular_constant",
                                                springAngularConstant);
  springConstantStretchArcs = (float)reader.GetReal(
      "springs", "stretch_constant_arcs", springConstantStretchArcs);
  springConstantSqueezeArcs = (float)reader.GetReal(
      "springs", "squeeze_constant_arcs", springConstantSqueezeArcs);

  loopStiffness = (float)reader.GetReal(
      "loop_energy", "loop_stiffness", loopStiffness);
  loopEquilibriumDistance = (float)reader.GetReal(
      "loop_energy", "loop_equilibrium_distance", loopEquilibriumDistance);
  loopConstraintsFile = reader.Get(
      "loop_energy", "loop_constraints_file", loopConstraintsFile);

  sphere_radius = (float)reader.GetReal("simulation", "sphere_radius", sphere_radius);

  enhancerAnnotationFile = reader.Get("annotations", "enhancer_bed", enhancerAnnotationFile);
  promoterAnnotationFile = reader.Get("annotations", "promoter_bed", promoterAnnotationFile);
  epigeneticAnnotationFile = reader.Get("annotations", "epigenetic_bed", epigeneticAnnotationFile);

  maxTempHeatmap = (float)reader.GetReal("simulation_heatmap",
                                         "max_temp_heatmap", maxTempHeatmap);
  dtTempHeatmap = (float)reader.GetReal("simulation_heatmap",
                                        "delta_temp_heatmap", dtTempHeatmap);
  tempJumpCoefHeatmap = (float)reader.GetReal(
      "simulation_heatmap", "jump_temp_coef_heatmap", tempJumpCoefHeatmap);
  tempJumpScaleHeatmap = (float)reader.GetReal(
      "simulation_heatmap", "jump_temp_scale_heatmap", tempJumpScaleHeatmap);
  MCstopConditionStepsHeatmap =
      reader.GetInteger("simulation_heatmap", "stop_condition_steps_heatmap",
                        MCstopConditionStepsHeatmap);
  MCstopConditionImprovementHeatmap = (float)reader.GetReal(
      "simulation_heatmap", "stop_condition_improvement_threshold_heatmap",
      MCstopConditionImprovementHeatmap);
  MCstopConditionMinSuccessesHeatmap = reader.GetInteger(
      "simulation_heatmap", "stop_condition_successes_threshold_heatmap",
      MCstopConditionMinSuccessesHeatmap);

  maxTempHeatmapDensity = (float)reader.GetReal(
      "simulation_heatmap_density", "max_temp_heatmap", maxTempHeatmapDensity);
  dtTempHeatmapDensity = (float)reader.GetReal(
      "simulation_heatmap_density", "delta_temp_heatmap", dtTempHeatmapDensity);
  tempJumpCoefHeatmapDensity = (float)reader.GetReal(
      "simulation_heatmap_density", "jump_temp_coef_heatmap",
      tempJumpCoefHeatmapDensity);
  tempJumpScaleHeatmapDensity = (float)reader.GetReal(
      "simulation_heatmap_density", "jump_temp_scale_heatmap",
      tempJumpScaleHeatmapDensity);
  MCstopConditionStepsHeatmapDensity = reader.GetInteger(
      "simulation_heatmap_density", "stop_condition_steps_heatmap",
      MCstopConditionStepsHeatmapDensity);
  MCstopConditionImprovementHeatmapDensity =
      (float)reader.GetReal("simulation_heatmap_density",
                            "stop_condition_improvement_threshold_heatmap",
                            MCstopConditionImprovementHeatmapDensity);
  MCstopConditionMinSuccessesHeatmapDensity =
      reader.GetInteger("simulation_heatmap_density",
                        "stop_condition_successes_threshold_heatmap",
                        MCstopConditionMinSuccessesHeatmapDensity);

  maxTemp = (float)reader.GetReal("simulation_arcs", "max_temp", maxTemp);
  dtTemp = (float)reader.GetReal("simulation_arcs", "delta_temp", dtTemp);
  tempJumpCoef =
      (float)reader.GetReal("simulation_arcs", "jump_temp_coef", tempJumpCoef);
  tempJumpScale = (float)reader.GetReal("simulation_arcs", "jump_temp_scale",
                                        tempJumpScale);
  MCstopConditionSteps = reader.GetInteger(
      "simulation_arcs", "stop_condition_steps", MCstopConditionSteps);
  MCstopConditionImprovement = (float)reader.GetReal(
      "simulation_arcs", "stop_condition_improvement_threshold",
      MCstopConditionImprovement);
  MCstopConditionMinSuccesses =
      reader.GetInteger("simulation_arcs", "stop_condition_successes_threshold",
                        MCstopConditionMinSuccesses);

  weightAngleSmooth = (float)reader.GetReal("simulation_arcs_smooth",
                                            "dist_weight", weightAngleSmooth);
  weightDistSmooth = (float)reader.GetReal("simulation_arcs_smooth",
                                           "angle_weight", weightDistSmooth);
  maxTempSmooth = (float)reader.GetReal("simulation_arcs_smooth", "max_temp",
                                        maxTempSmooth);
  dtTempSmooth = (float)reader.GetReal("simulation_arcs_smooth", "delta_temp",
                                       dtTempSmooth);
  tempJumpCoefSmooth = (float)reader.GetReal(
      "simulation_arcs_smooth", "jump_temp_coef", tempJumpCoefSmooth);
  tempJumpScaleSmooth = (float)reader.GetReal(
      "simulation_arcs_smooth", "jump_temp_scale", tempJumpScaleSmooth);
  MCstopConditionStepsSmooth =
      reader.GetInteger("simulation_arcs_smooth", "stop_condition_steps",
                        MCstopConditionStepsSmooth);
  MCstopConditionImprovementSmooth = (float)reader.GetReal(
      "simulation_arcs_smooth", "stop_condition_improvement_threshold",
      MCstopConditionImprovementSmooth);
  MCstopConditionMinSuccessesSmooth = reader.GetInteger(
      "simulation_arcs_smooth", "stop_condition_successes_threshold",
      MCstopConditionMinSuccessesSmooth);
  return true;
}

#endif /* SETTINGS_H_ */
