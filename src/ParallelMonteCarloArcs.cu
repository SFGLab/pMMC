/**
 * @file ParallelMonteCarloArcs.cu
 * @brief GPU-accelerated arc-level Monte Carlo for anchor reconstruction.
 *
 * Implements LooperSolver::parallelMonteCarloArcs(float step_size), the GPU
 * counterpart of the CPU MonteCarloArcs() function.  The energy function
 * mirrors the CPU version exactly:
 *
 *   E = arc_distance_score + loop_energy
 *
 * Arc distance score (per pair connected by an arc):
 *   diff = (|pos_i - pos_j| - expected_dist) / expected_dist
 *   score += diff^2 * (diff >= 0 ? stretchConst : squeezeConst)
 *
 * Pairs with expected_dist < 0 contribute a repulsion term (1/dist).
 * Pairs with expected_dist < 1e-6 are ignored (no arc, no repulsion).
 *
 * Loop energy (per explicit loop constraint):
 *   delta = |pos_i - pos_j| - eq_distance
 *   score += stiffness * delta^2
 *
 * Kernel design (warp-level independent SA):
 *   - Each warp (32 threads) runs an independent SA trajectory.
 *   - Thread 0 in the warp selects the bead and applies the move; all 32
 *     threads help evaluate the score in parallel via warp reduction.
 *   - __shfl_down_sync is used for warp-level summation.
 *   - No inter-warp global sync is needed.
 *   - Each warp cools independently with its own temperature.
 *   - The host picks the best final position across all warps.
 *
 * Host function:
 *   1. Flatten heatmap_exp_dist_anchor (NxN) to GPU (float*).
 *   2. Flatten loop constraints to GPU arrays.
 *   3. Upload float3 positions.
 *   4. Launch kernel.
 *   5. Download the position array from the warp with the best final score.
 *   6. Fall back to CPU MonteCarloArcs() if CUDA device unavailable or N<=1.
 *
 * GPU seed:
 *   The static variable s_arcs_gpu_seed is set once (from main via
 *   LooperSolver::setArcsGpuSeed) and reused across calls so that
 *   results are reproducible when a -S seed is provided.
 */

#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include <float.h>
#include <stdio.h>
#include <time.h>

#include <LooperSolver.h>

/* ======================================================================
 * Error-check macro (mirrors ParallelMonteCarloHeatmap.cu)
 * ====================================================================== */
#define gpuErrchkArcs(ans)                                                     \
  { gpuAssertArcs((ans), __FILE__, __LINE__); }
static inline void gpuAssertArcs(cudaError_t code, const char *file, int line,
                                 bool abort_on_error = true) {
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert (arcs): %s  %s:%d\n",
            cudaGetErrorString(code), file, line);
    if (abort_on_error)
      exit(code);
  }
}

/* ======================================================================
 * Static seed — set once by LooperSolver::setArcsGpuSeed()
 * ====================================================================== */
static unsigned int s_arcs_gpu_seed = 0; /* 0 = use time(NULL) */

/* ======================================================================
 * GPU settings struct for arc-level MC
 * ====================================================================== */
struct ArcsMCSettings {
  float maxTemp;
  float dtTemp;
  float tempJumpScale;
  float tempJumpCoef;
  float MCstopConditionImprovement;
  int   MCstopConditionSteps;
  int   MCstopConditionMinSuccesses;
  float springConstantStretchArcs;
  float springConstantSqueezeArcs;
  bool  use2D;
  int   milestoneFailsThreshold; /**< Max consecutive milestone failures before stopping. */
};

/* ======================================================================
 * Device-side RNG helpers
 * ====================================================================== */

/** Return a random integer in [0, range). */
__device__ __forceinline__ int arcs_rand_int(int range, curandState *st) {
  return (int)(curand(st) % (unsigned int)range);
}

/** Return a random float in (-max_size, +max_size). */
__device__ __forceinline__ float arcs_rand_float(float max_size,
                                                 curandState *st) {
  return (2.0f * curand_uniform(st) - 1.0f) * max_size;
}

/** Fill a float3 with a random displacement vector (optionally 2D). */
__device__ __forceinline__ float3 arcs_random_vector(float max_size, bool in2D,
                                                      curandState *st) {
  float3 v;
  v.x = arcs_rand_float(max_size, st);
  v.y = arcs_rand_float(max_size, st);
  v.z = in2D ? 0.0f : arcs_rand_float(max_size, st);
  return v;
}

/* ======================================================================
 * Device: score contribution for a single bead against all others.
 *
 * The n×n expected-distance matrix is stored row-major:
 *   exp_dist[row * n + col]
 *
 * Semantics (matching CPU calcScoreDistancesActiveRegion):
 *   < 0       → repulsion (no direct arc) — we skip repulsion for GPU
 *               simplicity; the dominant term is arc-distance pairs.
 *   [0, 1e-6) → ignored (diagonal or zero arc)
 *   >= 1e-6   → arc pair, penalise deviation from expected distance
 *
 * The work is split across the warp: each thread handles a stride of
 * positions, then results are reduced.
 * ====================================================================== */
__device__ float arcs_score_bead(int bead_idx, int n, float3 bead_pos,
                                  const float3 *__restrict__ positions,
                                  const float *__restrict__ exp_dist_mat,
                                  float stretch_k, float squeeze_k,
                                  int lane) {
  float sc = 0.0f;
  /* Each lane of the warp handles every 32nd pair. */
  for (int j = lane; j < n; j += 32) {
    if (j == bead_idx)
      continue;

    float ed = exp_dist_mat[bead_idx * n + j];
    if (ed < 1e-6f)
      continue; /* no arc or zero — skip (negative values mean no arc) */

    float3 pj = positions[j];
    float dx = bead_pos.x - pj.x;
    float dy = bead_pos.y - pj.y;
    float dz = bead_pos.z - pj.z;
    float dist = sqrtf(dx * dx + dy * dy + dz * dz);
    float diff = (dist - ed) / ed;
    sc += diff * diff * (diff >= 0.0f ? stretch_k : squeeze_k);
  }
  return sc;
}

/* ======================================================================
 * Device: loop-constraint score contribution for a single bead.
 *
 * loop_pairs[k]   = {active_i, active_j}  (both stored as int2)
 * loop_params[k]  = {stiffness, eq_distance}  (stored as float2)
 * ====================================================================== */
__device__ float arcs_loop_score_bead(int bead_idx,
                                       const float3 *__restrict__ positions,
                                       const int2 *__restrict__ loop_pairs,
                                       const float2 *__restrict__ loop_params,
                                       int n_loops,
                                       int lane) {
  float sc = 0.0f;
  for (int k = lane; k < n_loops; k += 32) {
    int ai = loop_pairs[k].x;
    int aj = loop_pairs[k].y;
    if (ai != bead_idx && aj != bead_idx)
      continue;

    float3 pi = positions[ai];
    float3 pj = positions[aj];
    /* If this bead is the moved one, its current pos is in positions[] already
     * (caller updates the position array before scoring). */
    float dx = pi.x - pj.x;
    float dy = pi.y - pj.y;
    float dz = pi.z - pj.z;
    float dist = sqrtf(dx * dx + dy * dy + dz * dz);
    float delta = dist - loop_params[k].y; /* dist - eq_distance */
    sc += loop_params[k].x * delta * delta; /* stiffness * delta^2 */
  }
  return sc;
}

/* ======================================================================
 * Warp-level float reduction (sum across all 32 lanes).
 * ====================================================================== */
#define FULL_WARP_MASK 0xffffffffu

__device__ __forceinline__ float warp_reduce_sum(float v) {
#pragma unroll
  for (int offset = 16; offset > 0; offset >>= 1)
    v += __shfl_down_sync(FULL_WARP_MASK, v, offset);
  return v; /* result valid in lane 0 */
}

/* ======================================================================
 * Setup kernel — initialise curandState per thread.
 * ====================================================================== */
__global__ void arcs_setupKernel(curandState *state, unsigned int seed) {
  int tid = blockDim.x * blockIdx.x + threadIdx.x;
  curand_init(seed, tid, 0, &state[tid]);
}

/* ======================================================================
 * Device variable: total iteration count (for logging, like heatmap kernel).
 * ====================================================================== */
__device__ int d_arcs_total_iterations;

/* ======================================================================
 * Main SA kernel.
 *
 * Grid: arbitrary blocks/threads; warps are the unit of work.
 * Each warp independently runs the full SA schedule.
 * Thread 0 of each warp drives the control flow; all threads in the
 * warp cooperatively compute energy via parallel reduction.
 *
 * After the kernel finishes, d_best_score[warp_id] and the positions in
 * d_best_positions[warp_id * n ... warp_id*n + n - 1] hold the result.
 * ====================================================================== */
__global__ void MonteCarloArcsKernel(
    /* RNG */
    curandState *__restrict__ state,
    /* Cluster data */
    float3 *__restrict__ positions,  /* [n] input positions (per warp copy in smem) */
    const int *__restrict__ is_fixed, /* [n] non-zero = fixed */
    /* Arc distance matrix */
    const float *__restrict__ exp_dist_mat, /* [n*n] row-major */
    int n,                                  /* active region size */
    /* Loop constraints */
    const int2 *__restrict__ loop_pairs,    /* [n_loops] */
    const float2 *__restrict__ loop_params, /* [n_loops] {stiffness, eq_dist} */
    int n_loops,
    /* SA settings */
    ArcsMCSettings settings,
    float step_size,
    /* Output: best score and positions per warp */
    float *__restrict__ d_best_score,         /* [n_warps] */
    float3 *__restrict__ d_best_positions,    /* [n_warps * n] */
    /* Shared scratch: each warp needs n float3 positions */
    /* Allocated dynamically: n * warps_per_block * sizeof(float3) bytes */
    /* Laid out as: smem[warp_in_block * n .. warp_in_block*n + n - 1] */
    bool *__restrict__ d_isDone) {
  /* ------------------------------------------------------------------
   * Thread/warp identification
   * ------------------------------------------------------------------ */
  int tid         = blockDim.x * blockIdx.x + threadIdx.x;
  int warp_id     = tid / 32;         /* global warp index */
  int lane        = tid & 31;         /* lane within warp [0..31] */
  int warp_in_blk = threadIdx.x / 32; /* warp index within block */

  /* ------------------------------------------------------------------
   * Shared memory: each warp gets n float3 slots for local positions.
   * Layout: float3 smem_pos[warps_per_block][n]
   * ------------------------------------------------------------------ */
  extern __shared__ float3 smem_pos[];
  float3 *my_pos = smem_pos + (size_t)warp_in_blk * n;

  /* Initialise local position copy from global (all lanes help). */
  for (int i = lane; i < n; i += 32) {
    my_pos[i] = positions[i];
  }
  __syncwarp();

  /* ------------------------------------------------------------------
   * Load curand state
   * ------------------------------------------------------------------ */
  curandState localState = state[tid];

  /* ------------------------------------------------------------------
   * SA state
   * ------------------------------------------------------------------ */
  float T          = settings.maxTemp;
  float score_curr = 0.0f;
  float score_prev = 0.0f;
  int   iterations = 0;
  int   imp_misses = 0;
  int   milestone_successes = 0;
  float milestone_score = 0.0f;

  /* Compute initial total score (all lanes contribute). */
  {
    float my_sc = 0.0f;
    /* Each lane handles a subset of bead pairs to avoid double-counting.
     * We sum over all beads: for each bead i (lane i%32), score vs all j>i. */
    for (int i = lane; i < n; i += 32) {
      for (int j = i + 1; j < n; j++) {
        float ed = exp_dist_mat[i * n + j];
        if (ed < 1e-6f) continue;
        float3 pi = my_pos[i], pj = my_pos[j];
        float dx = pi.x - pj.x, dy = pi.y - pj.y, dz = pi.z - pj.z;
        float dist = sqrtf(dx*dx + dy*dy + dz*dz);
        float diff = (dist - ed) / ed;
        my_sc += diff * diff *
                 (diff >= 0.f ? settings.springConstantStretchArcs
                              : settings.springConstantSqueezeArcs);
      }
    }
    /* Sum loop energy. */
    for (int k = lane; k < n_loops; k += 32) {
      int ai = loop_pairs[k].x, aj = loop_pairs[k].y;
      float3 pi = my_pos[ai], pj = my_pos[aj];
      float dx = pi.x-pj.x, dy = pi.y-pj.y, dz = pi.z-pj.z;
      float dist = sqrtf(dx*dx+dy*dy+dz*dz);
      float delta = dist - loop_params[k].y;
      my_sc += loop_params[k].x * delta * delta;
    }
    my_sc = warp_reduce_sum(my_sc);
    /* Broadcast to all lanes. */
    score_curr = __shfl_sync(FULL_WARP_MASK, my_sc, 0);
    score_prev = score_curr;
    milestone_score = score_curr;
  }

  /* ------------------------------------------------------------------
   * SA main loop — runs until convergence or isDone flag set.
   * ------------------------------------------------------------------ */
#define ARCS_INNER_N 256

  while (!(*d_isDone)) {

#pragma unroll 4
    for (int inner = 0; inner < ARCS_INNER_N; ++inner) {

      /* Lane 0 picks a random non-fixed bead. */
      int bead_idx = 0;
      float3 disp  = {0.f, 0.f, 0.f};
      if (lane == 0) {
        /* Keep trying until we find a non-fixed bead (bounded loop). */
        for (int attempt = 0; attempt < n * 4; ++attempt) {
          int candidate = arcs_rand_int(n, &localState);
          if (!is_fixed[candidate]) {
            bead_idx = candidate;
            break;
          }
        }
        disp = arcs_random_vector(step_size, settings.use2D, &localState);
      }
      /* Broadcast bead index and displacement to all lanes. */
      bead_idx = __shfl_sync(FULL_WARP_MASK, bead_idx, 0);
      disp.x   = __shfl_sync(FULL_WARP_MASK, disp.x, 0);
      disp.y   = __shfl_sync(FULL_WARP_MASK, disp.y, 0);
      disp.z   = __shfl_sync(FULL_WARP_MASK, disp.z, 0);

      /* ------ Score BEFORE move (local contribution of bead_idx). ------ */
      float sc_before = 0.0f;
      {
        float my_sc = arcs_score_bead(bead_idx, n, my_pos[bead_idx],
                                      my_pos, exp_dist_mat,
                                      settings.springConstantStretchArcs,
                                      settings.springConstantSqueezeArcs,
                                      lane);
        my_sc += arcs_loop_score_bead(bead_idx, my_pos,
                                      loop_pairs, loop_params, n_loops, lane);
        my_sc = warp_reduce_sum(my_sc);
        sc_before = __shfl_sync(FULL_WARP_MASK, my_sc, 0);
      }

      /* ------ Apply move in shared memory (all lanes see update). ------ */
      if (lane == 0) {
        my_pos[bead_idx].x += disp.x;
        my_pos[bead_idx].y += disp.y;
        my_pos[bead_idx].z += disp.z;
      }
      __syncwarp();

      /* ------ Score AFTER move. ------ */
      float sc_after = 0.0f;
      {
        float my_sc = arcs_score_bead(bead_idx, n, my_pos[bead_idx],
                                      my_pos, exp_dist_mat,
                                      settings.springConstantStretchArcs,
                                      settings.springConstantSqueezeArcs,
                                      lane);
        my_sc += arcs_loop_score_bead(bead_idx, my_pos,
                                      loop_pairs, loop_params, n_loops, lane);
        my_sc = warp_reduce_sum(my_sc);
        sc_after = __shfl_sync(FULL_WARP_MASK, my_sc, 0);
      }

      /* ------ Update total score incrementally. ------ */
      float new_total = score_prev - sc_before + sc_after;

      /* ------ Metropolis acceptance (lane 0 decides, broadcasts). ------ */
      bool accept = false;
      if (lane == 0) {
        if (new_total <= score_prev) {
          accept = true;
        } else if (T > 0.0f && score_prev > 0.0f) {
          float tp = settings.tempJumpScale *
                     expf(-settings.tempJumpCoef *
                          (new_total / score_prev) / T);
          accept = (curand_uniform(&localState) < tp);
        }
      }
      accept = (bool)__shfl_sync(FULL_WARP_MASK, (int)accept, 0);

      if (accept) {
        score_curr = new_total;
        if (lane == 0)
          milestone_successes++;
      } else {
        /* Revert move. */
        if (lane == 0) {
          my_pos[bead_idx].x -= disp.x;
          my_pos[bead_idx].y -= disp.y;
          my_pos[bead_idx].z -= disp.z;
        }
        __syncwarp();
        score_curr = score_prev;
      }

      /* Cool temperature once per inner step. */
      if (lane == 0)
        T *= settings.dtTemp;
      T = __shfl_sync(FULL_WARP_MASK, T, 0);

      score_prev = score_curr;
    } /* inner loop */

    iterations += ARCS_INNER_N;

    /* ------ Milestone check (lane 0 only, for warp 0 sets isDone). ------ */
    if (lane == 0) {
      bool no_improvement =
          (score_curr >
               settings.MCstopConditionImprovement * milestone_score &&
           milestone_successes < settings.MCstopConditionMinSuccesses) ||
          score_curr < 1e-5f ||
          (milestone_score > 1e-10f && score_curr / milestone_score > 0.9999f);

      if (no_improvement)
        imp_misses++;
      else
        imp_misses = 0;

      if (imp_misses >= settings.milestoneFailsThreshold || score_curr < 1e-5f) {
        *d_isDone = true;
      }

      milestone_score = score_curr;
      milestone_successes = 0;
    }
    /* Sync isDone flag across the warp before checking. */
    __syncwarp();
  } /* outer while */

  /* ------------------------------------------------------------------
   * Write per-warp best score and positions to global memory.
   * ------------------------------------------------------------------ */
  if (lane == 0) {
    d_best_score[warp_id] = score_curr;
    d_arcs_total_iterations = iterations; /* rough count; last-write wins */
  }
  /* All lanes cooperate to copy local positions to output. */
  for (int i = lane; i < n; i += 32) {
    d_best_positions[(size_t)warp_id * n + i] = my_pos[i];
  }

  /* Save RNG state. */
  state[tid] = localState;
}

/* ======================================================================
 * Host function: LooperSolver::parallelMonteCarloArcs
 * ====================================================================== */
float LooperSolver::parallelMonteCarloArcs(float step_size) {

  /* ------------------------------------------------------------------
   * Basic size checks.
   * ------------------------------------------------------------------ */
  int n = static_cast<int>(active_region.size());
  if (n <= 1)
    return 0.0f;

  /* ------------------------------------------------------------------
   * Rebuild loop map and check GPU availability.
   * ------------------------------------------------------------------ */
  rebuildActiveLoopMap();

  int deviceCount = 0;
  cudaError_t cudaStatus = cudaGetDeviceCount(&deviceCount);
  if (cudaStatus != cudaSuccess || deviceCount == 0) {
    printf("[parallelMonteCarloArcs] No CUDA device — falling back to CPU.\n");
    return static_cast<float>(MonteCarloArcs(step_size));
  }

  /* ------------------------------------------------------------------
   * Determine launch configuration.
   * ------------------------------------------------------------------ */
  int threads_per_block = Settings::cudaThreadsPerBlock;
  /* Ensure multiple of 32 (warp size). */
  threads_per_block = ((threads_per_block + 31) / 32) * 32;
  if (threads_per_block < 32) threads_per_block = 32;

  int blocks = static_cast<int>(Settings::cudaBlocksMultiplier *
                                active_region.size());
  if (blocks < 1) blocks = 1;

  int total_threads = threads_per_block * blocks;
  int n_warps       = total_threads / 32;

  /* Cap curandState allocation at 256MB. */
  {
    size_t state_bytes = (size_t)total_threads * sizeof(curandState);
    size_t max_bytes   = 256ULL * 1024 * 1024;
    if (state_bytes > max_bytes) {
      int max_threads = (int)(max_bytes / sizeof(curandState));
      max_threads = (max_threads / threads_per_block) * threads_per_block;
      if (max_threads < threads_per_block) max_threads = threads_per_block;
      total_threads = max_threads;
      blocks        = max_threads / threads_per_block;
      n_warps       = total_threads / 32;
      printf("[GPU Arcs] Capping to %d blocks (curandState 256MB limit)\n",
             blocks);
    }
  }

  /* Shared memory: one float3 slot per bead per warp-in-block. */
  int warps_per_block = threads_per_block / 32;
  size_t smem_bytes   = (size_t)warps_per_block * n * sizeof(float3);

  /* Check shared memory limit and occupancy. */
  {
    int dev = 0;
    cudaDeviceProp prop;
    cudaGetDevice(&dev);
    cudaGetDeviceProperties(&prop, dev);
    if (smem_bytes > (size_t)prop.sharedMemPerBlock) {
      printf("[GPU Arcs] smem %zu bytes > limit %zu bytes — "
             "reducing warps per block\n",
             smem_bytes, (size_t)prop.sharedMemPerBlock);
      /* Reduce warps_per_block until it fits. */
      warps_per_block = (int)(prop.sharedMemPerBlock / (n * sizeof(float3)));
      if (warps_per_block < 1) {
        printf("[GPU Arcs] Active region too large for GPU shared memory "
               "(%d beads) — falling back to CPU.\n", n);
        return static_cast<float>(MonteCarloArcs(step_size));
      }
      threads_per_block = warps_per_block * 32;
      blocks            = (total_threads + threads_per_block - 1) /
                          threads_per_block;
      total_threads     = threads_per_block * blocks;
      n_warps           = total_threads / 32;
      smem_bytes        = (size_t)warps_per_block * n * sizeof(float3);
    }

    /* Verify occupancy: reduce threads_per_block until the kernel can launch.
     * This handles Debug builds where -G inflates register usage. */
    int numBlocks = 0;
    cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &numBlocks, MonteCarloArcsKernel, threads_per_block,
        smem_bytes);
    while (numBlocks == 0 && threads_per_block > 32) {
      threads_per_block /= 2;
      warps_per_block = threads_per_block / 32;
      smem_bytes      = (size_t)warps_per_block * n * sizeof(float3);
      cudaOccupancyMaxActiveBlocksPerMultiprocessor(
          &numBlocks, MonteCarloArcsKernel, threads_per_block,
          smem_bytes);
    }
    if (numBlocks == 0) {
      printf("[GPU Arcs] Cannot achieve occupancy — falling back to CPU.\n");
      return static_cast<float>(MonteCarloArcs(step_size));
    }
    /* Recompute grid dims after possible thread reduction. */
    blocks        = (n_warps * 32 + threads_per_block - 1) / threads_per_block;
    if (blocks < 1) blocks = 1;
    total_threads = threads_per_block * blocks;
    n_warps       = total_threads / 32;
  }

  printf("[GPU Arcs] n=%d, threads=%d, blocks=%d, n_warps=%d, "
         "smem=%.1fKB, loops=%d\n",
         n, threads_per_block, blocks, n_warps,
         (float)smem_bytes / 1024.0f,
         (int)active_loop_map.size());

  /* ------------------------------------------------------------------
   * Build GPU settings struct.
   * ------------------------------------------------------------------ */
  ArcsMCSettings gpu_settings;
  gpu_settings.maxTemp                    = Settings::maxTemp;
  gpu_settings.dtTemp                     = Settings::dtTemp;
  gpu_settings.tempJumpScale              = Settings::tempJumpScale;
  gpu_settings.tempJumpCoef               = Settings::tempJumpCoef;
  gpu_settings.MCstopConditionImprovement = Settings::MCstopConditionImprovement;
  gpu_settings.MCstopConditionSteps       = Settings::MCstopConditionSteps;
  gpu_settings.MCstopConditionMinSuccesses= Settings::MCstopConditionMinSuccesses;
  gpu_settings.springConstantStretchArcs  = Settings::springConstantStretchArcs;
  gpu_settings.springConstantSqueezeArcs  = Settings::springConstantSqueezeArcs;
  gpu_settings.use2D                      = Settings::use2D;
  /* Reuse milestoneFailsThreshold from heatmap settings as convergence guard. */
  gpu_settings.milestoneFailsThreshold    = Settings::milestoneFailsThreshold;

  /* ------------------------------------------------------------------
   * Build host arrays for GPU upload.
   * ------------------------------------------------------------------ */

  /* 1. Positions (float3, full precision — arc distances need it). */
  thrust::host_vector<float3> h_positions(n);
  for (int i = 0; i < n; ++i) {
    h_positions[i].x = clusters[active_region[i]].pos.x;
    h_positions[i].y = clusters[active_region[i]].pos.y;
    h_positions[i].z = clusters[active_region[i]].pos.z;
  }

  /* 2. is_fixed flags (use int to avoid thrust::device_vector<bool> bit-packing). */
  thrust::host_vector<int> h_is_fixed(n);
  for (int i = 0; i < n; ++i)
    h_is_fixed[i] = clusters[active_region[i]].is_fixed ? 1 : 0;

  /* 3. Expected-distance matrix (NxN, row-major).
   *    heatmap_exp_dist_anchor.v[i][j] is already in active-region space.
   *    Negative values = no arc / repulsion;  [0, 1e-6) = ignore;
   *    >= 1e-6 = arc pair. */
  int n_sq = n * n;
  thrust::host_vector<float> h_exp_dist(n_sq, -1.0f);
  if ((int)heatmap_exp_dist_anchor.size == n) {
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        h_exp_dist[i * n + j] = heatmap_exp_dist_anchor.v[i][j];
  } else {
    printf("[GPU Arcs] Warning: heatmap_exp_dist_anchor size mismatch "
           "(%zu vs %d) — scoring disabled.\n",
           heatmap_exp_dist_anchor.size, n);
  }

  /* 4. Loop constraints — flatten active_loop_map. */
  int n_loops = static_cast<int>(active_loop_map.size());
  thrust::host_vector<int2>   h_loop_pairs(n_loops);
  thrust::host_vector<float2> h_loop_params(n_loops);
  {
    int k = 0;
    for (const auto &entry : active_loop_map) {
      h_loop_pairs[k].x  = entry.first.first;
      h_loop_pairs[k].y  = entry.first.second;
      const LoopConstraint &lc = loop_constraints[entry.second];
      h_loop_params[k].x = lc.stiffness;
      h_loop_params[k].y = lc.eq_distance;
      ++k;
    }
  }

  /* ------------------------------------------------------------------
   * GPU allocations via thrust.
   * ------------------------------------------------------------------ */
  thrust::device_vector<float3> d_positions   = h_positions;
  thrust::device_vector<int>    d_is_fixed    = h_is_fixed;
  thrust::device_vector<float>  d_exp_dist    = h_exp_dist;

  thrust::device_vector<int2>   d_loop_pairs;
  thrust::device_vector<float2> d_loop_params;
  if (n_loops > 0) {
    d_loop_pairs  = h_loop_pairs;
    d_loop_params = h_loop_params;
  } else {
    /* Allocate dummy 1-element arrays to keep device pointers valid. */
    d_loop_pairs.resize(1);
    d_loop_params.resize(1);
  }

  /* Per-warp output: best score and corresponding positions. */
  thrust::device_vector<float>  d_best_score(n_warps, FLT_MAX);
  thrust::device_vector<float3> d_best_positions((size_t)n_warps * n);

  /* isDone flag. */
  bool *d_isDone = nullptr;
  gpuErrchkArcs(cudaMalloc(&d_isDone, sizeof(bool)));
  gpuErrchkArcs(cudaMemset(d_isDone, 0, sizeof(bool)));

  /* curandState. */
  curandState *d_states = nullptr;
  gpuErrchkArcs(
      cudaMalloc(&d_states, (size_t)total_threads * sizeof(curandState)));

  /* ------------------------------------------------------------------
   * Seed selection (reproducible if s_arcs_gpu_seed != 0).
   * ------------------------------------------------------------------ */
  unsigned int seed = (s_arcs_gpu_seed != 0)
                          ? s_arcs_gpu_seed
                          : static_cast<unsigned int>(time(NULL));

  arcs_setupKernel<<<blocks, threads_per_block>>>(d_states, seed);
  gpuErrchkArcs(cudaDeviceSynchronize());

  /* ------------------------------------------------------------------
   * Compute initial score on host (for logging).
   * ------------------------------------------------------------------ */
  double score_init =
      calcScoreDistancesActiveRegion() + calcScoreLoopEnergy();
  output(4, "[GPU Arcs] initial score = %lf, step_size = %f\n",
         score_init, step_size);

  /* ------------------------------------------------------------------
   * Launch kernel.
   * ------------------------------------------------------------------ */
  MonteCarloArcsKernel<<<blocks, threads_per_block, smem_bytes>>>(
      d_states,
      thrust::raw_pointer_cast(d_positions.data()),
      thrust::raw_pointer_cast(d_is_fixed.data()),
      thrust::raw_pointer_cast(d_exp_dist.data()),
      n,
      thrust::raw_pointer_cast(d_loop_pairs.data()),
      thrust::raw_pointer_cast(d_loop_params.data()),
      n_loops,
      gpu_settings,
      step_size,
      thrust::raw_pointer_cast(d_best_score.data()),
      thrust::raw_pointer_cast(d_best_positions.data()),
      d_isDone);

  gpuErrchkArcs(cudaPeekAtLastError());
  gpuErrchkArcs(cudaDeviceSynchronize());

  /* ------------------------------------------------------------------
   * Find the warp with the best (lowest) score.
   * ------------------------------------------------------------------ */
  thrust::host_vector<float>  h_best_score    = d_best_score;
  thrust::host_vector<float3> h_best_positions = d_best_positions;

  int best_warp   = 0;
  float best_sc   = h_best_score[0];
  for (int w = 1; w < n_warps; ++w) {
    if (h_best_score[w] < best_sc) {
      best_sc   = h_best_score[w];
      best_warp = w;
    }
  }

  /* ------------------------------------------------------------------
   * Write back best positions to CPU clusters array.
   * ------------------------------------------------------------------ */
  for (int i = 0; i < n; ++i) {
    float3 p = h_best_positions[(size_t)best_warp * n + i];
    clusters[active_region[i]].pos.x = p.x;
    clusters[active_region[i]].pos.y = p.y;
    clusters[active_region[i]].pos.z = p.z;
  }

  /* ------------------------------------------------------------------
   * Retrieve GPU iteration count for logging.
   * ------------------------------------------------------------------ */
  int gpu_iters = 0;
  cudaMemcpyFromSymbol(&gpu_iters, d_arcs_total_iterations, sizeof(int));
  total_mc_steps += gpu_iters;

  /* ------------------------------------------------------------------
   * Recompute final score on CPU for accuracy (matches return type of
   * MonteCarloArcs).
   * ------------------------------------------------------------------ */
  double score_final =
      calcScoreDistancesActiveRegion() + calcScoreLoopEnergy();

  printf("[GPU Arcs] FINAL SCORE = %f  (best warp=%d, GPU iters=%d)\n",
         (float)score_final, best_warp, gpu_iters);

  /* ------------------------------------------------------------------
   * Cleanup.
   * ------------------------------------------------------------------ */
  cudaFree(d_isDone);
  cudaFree(d_states);

  return static_cast<float>(score_final);
}

/* ======================================================================
 * Public helper — call from main() after parsing -S seed.
 * Not declared in the header (internal linkage), kept here for locality.
 * ====================================================================== */
void LooperSolver_setArcsGpuSeed(unsigned int seed) {
  s_arcs_gpu_seed = seed;
}
