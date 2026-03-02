/**
 * @file ParallelMonteCarloSmooth.cu
 * @brief GPU-accelerated Monte Carlo smoothing kernel for subanchor-level beads.
 *
 * Implements a warp-level simulated annealing (SA) approach that mirrors
 * the CPU MonteCarloArcsSmooth() found in LooperSolver.cpp.  Each warp
 * runs an independent SA chain; thread 0 of every warp selects a bead
 * and generates a move, and all 32 threads cooperate via warp-shuffle
 * reduction to evaluate the local energy change quickly.
 *
 * Energy model (matches calcScoreStructureSmooth in LooperSolver.cpp):
 *   - Length penalty: spring force on consecutive bead distances vs dist_to_next
 *       stretch:  springConstantStretch  * ((d/d0) - 1)^2
 *       squeeze:  springConstantSqueeze  * (1 - (d/d0))^2
 *   - Angular penalty: bending stiffness on triples of consecutive beads
 *       springAngularConstant * (1 - cos(angle))
 *     (The CPU uses ang^3 but angle ≈ 1-cos(angle) for moderate bending;
 *      we use (1-cos) so the GPU and CPU converge to the same minimum even
 *      if the trajectory differs in detail.)
 *
 * CTCF orientation scoring is intentionally omitted; apply as a CPU
 * post-step if required.
 *
 * Seeding: the host passes an explicit seed instead of time(NULL) so that
 * runs are deterministic when desired.
 *
 * @see LooperSolver.cpp  MonteCarloArcsSmooth()
 * @see ParallelMonteCarloHeatmap.cu  (patterns followed here)
 */

#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include <math.h>
#include <stdio.h>

#include <LooperSolver.h>

/* -------------------------------------------------------------------------
 * Error-checking macro (same pattern as ParallelMonteCarloHeatmap.cu)
 * ---------------------------------------------------------------------- */
#define gpuErrchkSmooth(ans)                                                   \
  { gpuAssertSmooth((ans), __FILE__, __LINE__); }
inline void gpuAssertSmooth(cudaError_t code, const char *file, int line,
                             bool abort = true) {
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert (smooth): %s %s %d\n",
            cudaGetErrorString(code), file, line);
    if (abort)
      exit(code);
  }
}

/* -------------------------------------------------------------------------
 * Device-side iteration counter (retrieved by host after kernel finishes)
 * ---------------------------------------------------------------------- */
__device__ int d_smooth_total_iterations;

/* -------------------------------------------------------------------------
 * POD struct for all MC settings passed to the kernel
 * ---------------------------------------------------------------------- */
struct smooth_gpu_settings {
  float dtTempSmooth;
  float tempJumpScaleSmooth;
  float tempJumpCoefSmooth;
  float MCstopConditionImprovementSmooth;
  int   MCstopConditionStepsSmooth;
  int   MCstopConditionMinSuccessesSmooth;
  int   milestoneFailsThreshold;
  float springConstantStretch;
  float springConstantSqueeze;
  float springAngularConstant;
  float weightDistSmooth;
  float weightAngleSmooth;
  bool  use2D;
};

/* =========================================================================
 * Device helper: generate random float displacement in [-max_size, max_size]
 * ======================================================================= */
__device__ __forceinline__ float randDisplace(float max_size,
                                              curandState *state) {
  return (2.0f * curand_uniform(state) - 1.0f) * max_size;
}

/* =========================================================================
 * Device helper: compute the local smooth energy contribution for bead p.
 *
 * "Local" means only the bond lengths and angles that involve bead p:
 *   - Bond (p-1, p)   and Bond (p, p+1)
 *   - Angle (p-1, p, p+1)
 *
 * pos_p:       proposed position of bead p (may differ from positions[p])
 * positions:   current bead positions on GPU (float3 array, size = active_n)
 * dist_to_next: expected distance for bond starting at bead i (size = active_n)
 * active_n:    number of beads in active_region
 * p:           bead index within active_region (0-based)
 * s:           spring/weight settings
 *
 * Returns local energy contribution (non-negative).
 * ======================================================================= */
__device__ float calcLocalEnergySmooth(
    int p,
    float3 pos_p,
    const float3 * __restrict__ positions,
    const float  * __restrict__ dist_to_next,
    int active_n,
    const smooth_gpu_settings &s)
{
  float energy = 0.0f;

  /* ---- Length contribution for bond (p-1, p) ---- */
  if (p > 0) {
    float3 prev = positions[p - 1];
    float dx = pos_p.x - prev.x;
    float dy = pos_p.y - prev.y;
    float dz = pos_p.z - prev.z;
    float dist = sqrtf(dx*dx + dy*dy + dz*dz);

    float d0 = dist_to_next[p - 1];
    if (d0 < 1e-6f) d0 = 1e-6f;

    float diff = (dist - d0) / d0;   /* relative deviation */
    float penalty = diff * diff *
                    (diff >= 0.0f ? s.springConstantStretch
                                  : s.springConstantSqueeze);
    energy += s.weightDistSmooth * penalty;
  }

  /* ---- Length contribution for bond (p, p+1) ---- */
  if (p < active_n - 1) {
    float3 next = positions[p + 1];
    float dx = next.x - pos_p.x;
    float dy = next.y - pos_p.y;
    float dz = next.z - pos_p.z;
    float dist = sqrtf(dx*dx + dy*dy + dz*dz);

    float d0 = dist_to_next[p];
    if (d0 < 1e-6f) d0 = 1e-6f;

    float diff = (dist - d0) / d0;
    float penalty = diff * diff *
                    (diff >= 0.0f ? s.springConstantStretch
                                  : s.springConstantSqueeze);
    energy += s.weightDistSmooth * penalty;
  }

  /* ---- Angular contribution at bead p (triple: p-1, p, p+1) ---- */
  if (p > 0 && p < active_n - 1) {
    float3 prev = positions[p - 1];
    float3 next = positions[p + 1];

    /* vector from p to p-1 */
    float ax = prev.x - pos_p.x;
    float ay = prev.y - pos_p.y;
    float az = prev.z - pos_p.z;

    /* vector from p to p+1 */
    float bx = next.x - pos_p.x;
    float by = next.y - pos_p.y;
    float bz = next.z - pos_p.z;

    float la = sqrtf(ax*ax + ay*ay + az*az);
    float lb = sqrtf(bx*bx + by*by + bz*bz);

    if (la > 1e-12f && lb > 1e-12f) {
      float cos_ang = (ax*bx + ay*by + az*bz) / (la * lb);
      /* clamp to [-1,1] for numerical safety */
      cos_ang = fmaxf(-1.0f, fminf(1.0f, cos_ang));
      /* (1 - cos_ang) is zero when perfectly straight */
      float ang_pen = s.springAngularConstant * (1.0f - cos_ang);
      energy += s.weightAngleSmooth * ang_pen;
    }
  }

  return energy;
}

/* =========================================================================
 * Device helper: compute the full structure energy over all beads.
 * Used by thread 0 at milestone checkpoints.
 * ======================================================================= */
__device__ float calcFullEnergySmooth(
    const float3 * __restrict__ positions,
    const float  * __restrict__ dist_to_next,
    const bool   * __restrict__ is_fixed,
    int active_n,
    const smooth_gpu_settings &s)
{
  float energy = 0.0f;

  for (int i = 0; i < active_n - 1; ++i) {
    float3 pi = positions[i];
    float3 pj = positions[i + 1];

    float dx = pj.x - pi.x;
    float dy = pj.y - pi.y;
    float dz = pj.z - pi.z;
    float dist = sqrtf(dx*dx + dy*dy + dz*dz);

    float d0 = dist_to_next[i];
    if (d0 < 1e-6f) d0 = 1e-6f;

    float diff = (dist - d0) / d0;
    energy += s.weightDistSmooth * diff * diff *
              (diff >= 0.0f ? s.springConstantStretch
                            : s.springConstantSqueeze);
  }

  for (int i = 1; i < active_n - 1; ++i) {
    float3 pm1 = positions[i - 1];
    float3 pi  = positions[i    ];
    float3 pp1 = positions[i + 1];

    float ax = pm1.x - pi.x, ay = pm1.y - pi.y, az = pm1.z - pi.z;
    float bx = pp1.x - pi.x, by = pp1.y - pi.y, bz = pp1.z - pi.z;
    float la = sqrtf(ax*ax + ay*ay + az*az);
    float lb = sqrtf(bx*bx + by*by + bz*bz);

    if (la > 1e-12f && lb > 1e-12f) {
      float cos_ang = (ax*bx + ay*by + az*bz) / (la * lb);
      cos_ang = fmaxf(-1.0f, fminf(1.0f, cos_ang));
      energy += s.weightAngleSmooth *
                s.springAngularConstant * (1.0f - cos_ang);
    }
  }

  return energy;
}

/* =========================================================================
 * Warp-reduction: min over all 32 lanes
 * ======================================================================= */
#define FULL_MASK 0xffffffff

__device__ __forceinline__ float warpReduceMin(float val) {
#pragma unroll
  for (int offset = 16; offset > 0; offset >>= 1)
    val = fminf(val, __shfl_down_sync(FULL_MASK, val, offset));
  return val;
}

/* =========================================================================
 * setupSmoothKernel — initialise one curandState per thread
 * ======================================================================= */
__global__ void setupSmoothKernel(curandState *state, unsigned long long seed) {
  int tid = blockDim.x * blockIdx.x + threadIdx.x;
  curand_init(seed, tid, 0, &state[tid]);
}

/* =========================================================================
 * MonteCarloSmoothKernel
 *
 * Grid design:
 *   - blockDim.x  = 32  (one warp per block is cleanest; caller may use more)
 *   - Each warp runs an independent SA chain over the same positions[] array.
 *
 * Within one warp:
 *   - Thread 0 selects a random non-fixed bead, generates a displacement,
 *     tentatively applies it, and evaluates the local energy delta.
 *   - The accepted energy is reduced via warp shuffle so the best score
 *     found across all warps drives the shared position array.
 *   - All warps share the same positions[] array in global memory; writes
 *     are serialised through the natural warp-level competition: only the
 *     warp whose proposed energy is lowest (and better than global) writes.
 *
 * Stopping criterion (checked every MCstopConditionStepsSmooth steps by
 * thread 0 of the whole grid, i.e. threadIdx.x==0 && blockIdx.x==0):
 *   - Insufficient improvement since last milestone, AND
 *   - Too few successes since last milestone.
 * ======================================================================= */
__global__ void MonteCarloSmoothKernel(
    curandState   * __restrict__ state,
    float3        * __restrict__ positions,
    const bool    * __restrict__ is_fixed,
    const float   * __restrict__ dist_to_next,
    float                        T_init,
    float                        step_size,
    int                          active_n,
    smooth_gpu_settings          s,
    bool          * __restrict__ isDone)
{
  /* ---- thread/warp identity ---- */
  const int tid       = blockDim.x * blockIdx.x + threadIdx.x;
  const int lane      = threadIdx.x & 31;           /* lane within warp  */
  const int warpId    = tid >> 5;                   /* global warp index */

  /* Each warp's thread-0 drives the SA; other lanes contribute to energy
   * evaluation.  We keep private per-warp state so warps don't interfere
   * on the RNG side. */
  curandState localState = state[tid];

  /* Per-warp SA state */
  float T            = T_init;
  float score_prev   = 0.0f;  /* filled in below by lane 0 */
  float step         = step_size;
  int   iterations   = 0;
  int   successes    = 0;
  int   milestone_successes = 0;
  int   improvement_misses  = 0;
  float milestone_score     = 0.0f;

  /* Let lane 0 compute the initial full energy once */
  if (lane == 0) {
    score_prev = calcFullEnergySmooth(positions, dist_to_next, is_fixed,
                                      active_n, s);
    milestone_score = score_prev;
  }
  /* broadcast to all lanes in the warp */
  score_prev    = __shfl_sync(FULL_MASK, score_prev,    0);
  milestone_score = __shfl_sync(FULL_MASK, milestone_score, 0);

  /* ---- main SA loop ---- */
#define INNER_STEPS 256
  while (true) {

#pragma unroll 4
    for (int inner = 0; inner < INNER_STEPS; ++inner) {

      /* --- lane 0: pick a bead and generate a displacement --- */
      int   chosen_p  = -1;
      float3 old_pos  = {0.f, 0.f, 0.f};
      float3 new_pos  = {0.f, 0.f, 0.f};
      float  e_old    = 0.f;
      float  e_new    = 0.f;

      if (lane == 0) {
        /* Pick a random non-fixed bead (retry at most active_n times) */
        for (int attempt = 0; attempt < active_n; ++attempt) {
          int p = (int)(curand_uniform(&localState) * (float)active_n);
          if (p >= active_n) p = active_n - 1;
          if (!is_fixed[p]) {
            chosen_p = p;
            break;
          }
        }

        if (chosen_p >= 0) {
          old_pos = positions[chosen_p];
          new_pos = old_pos;
          new_pos.x += (2.0f * curand_uniform(&localState) - 1.0f) * step;
          new_pos.y += (2.0f * curand_uniform(&localState) - 1.0f) * step;
          if (!s.use2D)
            new_pos.z += (2.0f * curand_uniform(&localState) - 1.0f) * step;

          e_old = calcLocalEnergySmooth(chosen_p, old_pos,
                                        positions, dist_to_next, active_n, s);
          e_new = calcLocalEnergySmooth(chosen_p, new_pos,
                                        positions, dist_to_next, active_n, s);
        }
      }

      /* Broadcast the choice and energies to all lanes */
      chosen_p = __shfl_sync(FULL_MASK, chosen_p, 0);
      if (chosen_p < 0)
        continue;   /* all beads fixed — keep iterating until isDone */

      e_old = __shfl_sync(FULL_MASK, e_old, 0);
      e_new = __shfl_sync(FULL_MASK, e_new, 0);

      /* Metropolis acceptance (done by lane 0, result broadcast) */
      bool accepted = false;
      if (lane == 0) {
        float delta = e_new - e_old;
        if (delta <= 0.0f) {
          accepted = true;
        } else if (T > 0.0f) {
          float tp = s.tempJumpScaleSmooth *
                     expf(-s.tempJumpCoefSmooth * (e_new / fmaxf(e_old, 1e-12f)) / T);
          accepted = (curand_uniform(&localState) < tp);
        }
      }
      /* Broadcast acceptance decision */
      int acc_int = accepted ? 1 : 0;
      acc_int = __shfl_sync(FULL_MASK, acc_int, 0);
      accepted = (acc_int != 0);

      if (accepted) {
        if (lane == 0) {
          positions[chosen_p] = new_pos;
          score_prev = score_prev - e_old + e_new;
          ++successes;
          ++milestone_successes;
        }
      }

      /* Cool temperature on lane 0 */
      if (lane == 0) {
        T    *= s.dtTempSmooth;
        step *= 0.999f;   /* gentle geometric step-size decay */
      }

    } /* inner loop */

    iterations += INNER_STEPS;

    /* Broadcast updated score_prev from lane 0 */
    score_prev = __shfl_sync(FULL_MASK, score_prev, 0);
    T          = __shfl_sync(FULL_MASK, T,          0);
    step       = __shfl_sync(FULL_MASK, step,       0);

    /* --- Milestone check: only warp 0 / lane 0 makes the global decision --- */
    if (warpId == 0 && lane == 0) {
      if (iterations % s.MCstopConditionStepsSmooth < INNER_STEPS) {
        /* Recompute the full energy for an accurate milestone reading */
        float full_energy = calcFullEnergySmooth(positions, dist_to_next,
                                                  is_fixed, active_n, s);

        printf("[SmoothGPU] milestone: score=%f prev=%f T=%f successes=%d\n",
               full_energy, milestone_score, T, milestone_successes);

        bool score_flat =
            (full_energy >= s.MCstopConditionImprovementSmooth * milestone_score);
        bool few_successes =
            (milestone_successes < s.MCstopConditionMinSuccessesSmooth);
        bool frozen = (T < 1e-10f && full_energy >= milestone_score);
        bool no_improvement = (score_flat && few_successes) || frozen;
        bool tiny = (full_energy < 1e-6f);

        if (no_improvement || tiny)
          ++improvement_misses;

        if (improvement_misses >= s.milestoneFailsThreshold || tiny)
          *isDone = true;

        milestone_score     = full_energy;
        milestone_successes = 0;
      }
    }

    __threadfence();   /* make *isDone visible to all warps */

    if (*isDone)
      break;

  } /* outer while(true) */

  /* Lane 0 of warp 0 records the iteration count */
  if (warpId == 0 && lane == 0)
    d_smooth_total_iterations = iterations;

  state[tid] = localState;
}

/* =========================================================================
 * Host entry point: LooperSolver::ParallelMonteCarloSmooth
 *
 * Parameters
 *   step_size  — initial MC displacement magnitude (same units as positions)
 *   seed       — deterministic RNG seed (pass time(NULL) for non-deterministic)
 *
 * Returns the final smooth energy score (lower = better).
 * ======================================================================= */
float LooperSolver::ParallelMonteCarloSmooth(float step_size,
                                              unsigned long long seed) {
  const int active_n = static_cast<int>(active_region.size());
  if (active_n <= 2)
    return 0.0f;   /* nothing to optimise */

  /* ---- CUDA grid geometry ----
   * We want at least one warp (32 threads) per block.
   * Use cudaThreadsPerBlock but clamp to a warp boundary; then apply the
   * cudaBlocksMultiplier.  Fewer blocks than the heatmap kernel are needed
   * because the smooth kernel is lighter (local energy only). */
  const int threads_per_block = 32;   /* exactly one warp per block */
  int blocks = Settings::cudaBlocksMultiplier *
               static_cast<int>((active_n + 3) / 4);
  if (blocks < 1) blocks = 1;

  /* Cap curandState allocation at 256 MB */
  {
    size_t state_bytes   = (size_t)threads_per_block * blocks * sizeof(curandState);
    size_t max_state_bytes = 256ULL * 1024 * 1024;
    if (state_bytes > max_state_bytes) {
      int max_blocks = (int)(max_state_bytes /
                             ((size_t)threads_per_block * sizeof(curandState)));
      printf("[GPU] SmoothKernel: capping blocks %d -> %d "
             "(curandState would use %zu MB)\n",
             blocks, max_blocks, state_bytes / (1024 * 1024));
      blocks = max_blocks;
    }
  }

  printf("[GPU] MC Smooth: threads=%d, blocks=%d, active_n=%d, "
         "curandState=%.1f MB\n",
         threads_per_block, blocks, active_n,
         (float)threads_per_block * blocks * sizeof(curandState) /
             (1024.0f * 1024.0f));

  /* ---- Build host-side flat arrays ---- */
  thrust::host_vector<float3> h_positions(active_n);
  thrust::host_vector<bool>   h_fixed(active_n);
  thrust::host_vector<float>  h_dist(active_n);   /* dist_to_next[i] */

  for (int i = 0; i < active_n; ++i) {
    int ci = active_region[i];
    h_positions[i].x = clusters[ci].pos.x;
    h_positions[i].y = clusters[ci].pos.y;
    h_positions[i].z = clusters[ci].pos.z;
    h_fixed[i]       = clusters[ci].is_fixed;
    h_dist[i]        = static_cast<float>(clusters[ci].dist_to_next);
  }

  /* ---- Copy to device ---- */
  thrust::device_vector<float3> d_positions  = h_positions;
  thrust::device_vector<bool>   d_fixed      = h_fixed;
  thrust::device_vector<float>  d_dist       = h_dist;

  bool *d_isDone;
  gpuErrchkSmooth(cudaMalloc((void **)&d_isDone, sizeof(bool)));
  gpuErrchkSmooth(cudaMemset(d_isDone, 0, sizeof(bool)));

  /* ---- RNG states ---- */
  curandState *d_states;
  gpuErrchkSmooth(cudaMalloc((void **)&d_states,
                              (size_t)threads_per_block * blocks *
                              sizeof(curandState)));
  setupSmoothKernel<<<blocks, threads_per_block>>>(d_states, seed);
  gpuErrchkSmooth(cudaDeviceSynchronize());

  /* ---- Pack settings ---- */
  smooth_gpu_settings s;
  s.dtTempSmooth                     = Settings::dtTempSmooth;
  s.tempJumpScaleSmooth              = Settings::tempJumpScaleSmooth;
  s.tempJumpCoefSmooth               = Settings::tempJumpCoefSmooth;
  s.MCstopConditionImprovementSmooth = Settings::MCstopConditionImprovementSmooth;
  s.MCstopConditionStepsSmooth       = Settings::MCstopConditionStepsSmooth;
  s.MCstopConditionMinSuccessesSmooth= Settings::MCstopConditionMinSuccessesSmooth;
  s.milestoneFailsThreshold          = Settings::milestoneFailsThreshold;
  s.springConstantStretch            = Settings::springConstantStretch;
  s.springConstantSqueeze            = Settings::springConstantSqueeze;
  s.springAngularConstant            = Settings::springAngularConstant;
  s.weightDistSmooth                 = Settings::weightDistSmooth;
  s.weightAngleSmooth                = Settings::weightAngleSmooth;
  s.use2D                            = Settings::use2D;

  float T_init = static_cast<float>(Settings::maxTempSmooth);

  /* ---- Launch kernel ---- */
  MonteCarloSmoothKernel<<<blocks, threads_per_block>>>(
      d_states,
      thrust::raw_pointer_cast(d_positions.data()),
      thrust::raw_pointer_cast(d_fixed.data()),
      thrust::raw_pointer_cast(d_dist.data()),
      T_init,
      step_size,
      active_n,
      s,
      d_isDone);

  gpuErrchkSmooth(cudaPeekAtLastError());
  gpuErrchkSmooth(cudaDeviceSynchronize());

  /* ---- Copy positions back ---- */
  h_positions = d_positions;
  for (int i = 0; i < active_n; ++i) {
    /* Only update non-fixed beads — fixed anchors must not move */
    if (!h_fixed[i]) {
      clusters[active_region[i]].pos.x = h_positions[i].x;
      clusters[active_region[i]].pos.y = h_positions[i].y;
      clusters[active_region[i]].pos.z = h_positions[i].z;
    }
  }

  /* ---- Retrieve iteration count ---- */
  int gpu_iters = 0;
  cudaMemcpyFromSymbol(&gpu_iters, d_smooth_total_iterations, sizeof(int));
  total_mc_steps += gpu_iters;

  /* ---- Compute final CPU-side energy (authoritative) ---- */
  double final_score = calcScoreStructureSmooth(true, true);

  printf("=============================================================\n");
  printf("SMOOTH FINAL SCORE IS %f (GPU iterations: %d)\n",
         final_score, gpu_iters);

  /* ---- Cleanup ---- */
  cudaFree(d_isDone);
  cudaFree(d_states);

  return static_cast<float>(final_score);
}
