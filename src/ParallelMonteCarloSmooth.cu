/**
 * @file ParallelMonteCarloSmooth.cu
 * @brief GPU-accelerated Monte Carlo smoothing kernel for subanchor-level beads.
 *
 * Multi-warp design: each warp runs an independent SA trajectory on its own
 * copy of the position array. The host picks the best result. This gives
 * N_warps parallel SA chains (typically 46+ on modern GPUs), dramatically
 * speeding up the smooth phase compared to a single sequential trajectory.
 *
 * Energy model (matches calcScoreStructureSmooth in LooperSolver.cpp):
 *   - Length penalty: spring force on consecutive bead distances vs dist_to_next
 *   - Angular penalty: bending stiffness on triples of consecutive beads
 *
 * @see LooperSolver.cpp  MonteCarloArcsSmooth()
 */

#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <vector>

#include <LooperSolver.hpp>

/* -------------------------------------------------------------------------
 * Error-checking macro
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
 * Device-side iteration counter
 * ---------------------------------------------------------------------- */
__device__ int d_smooth_total_iterations;
__device__ int d_smooth_total_accepted;

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
  float3 sphere_center;
  float sphere_radius;
};

__device__ __forceinline__ float sphereBoundaryPenalty(float3 pos, float3 center, float radius) {
    if (radius <= 0.0f) return 0.0f;
    float dx = pos.x - center.x;
    float dy = pos.y - center.y;
    float dz = pos.z - center.z;
    float dist = sqrtf(dx*dx + dy*dy + dz*dz);
    if (dist <= radius) return 0.0f;
    float overshoot = dist - radius;
    return overshoot * overshoot;
}

/* =========================================================================
 * Device helper: compute the local smooth energy for bead p.
 * Only checks bonds (p-1,p) and (p,p+1) and angles involving p.
 * Double precision. Returns sca + scb WITHOUT weights.
 * ======================================================================= */
__device__ float calcLocalEnergySmooth(
    int p,
    float3 pos_p,
    const float3 * __restrict__ positions,
    const float  * __restrict__ dist_to_next,
    int active_n,
    const smooth_gpu_settings &s)
{
  float sca = 0.0f, scb = 0.0f;

  /* Length: bonds (p-1,p) and (p,p+1) */
  if (p > 0) {
    float3 prev = positions[p - 1];
    float dx = pos_p.x - prev.x;
    float dy = pos_p.y - prev.y;
    float dz = pos_p.z - prev.z;
    float dist = sqrtf(dx*dx + dy*dy + dz*dz);
    float d0 = dist_to_next[p - 1];
    if (d0 < 1e-6f) d0 = 1e-6f;
    float diff = (dist - d0) / d0;
    sca += diff * diff * (diff >= 0.0f ? s.springConstantStretch
                                       : s.springConstantSqueeze);
  }
  if (p < active_n - 1) {
    float3 next = positions[p + 1];
    float dx = next.x - pos_p.x;
    float dy = next.y - pos_p.y;
    float dz = next.z - pos_p.z;
    float dist = sqrtf(dx*dx + dy*dy + dz*dz);
    float d0 = dist_to_next[p];
    if (d0 < 1e-6f) d0 = 1e-6f;
    float diff = (dist - d0) / d0;
    sca += diff * diff * (diff >= 0.0f ? s.springConstantStretch
                                       : s.springConstantSqueeze);
  }

  /* Angular: 3 angles around moved bead */
  {
    float prev_vx = 0.0f, prev_vy = 0.0f, prev_vz = 0.0f;
    bool have_prev = false;
    for (int i = p - 2; i < p + 2; i++) {
      if (i < 0 || i + 1 >= active_n)
        continue;
      float3 pi_pos = (i == p) ? pos_p : positions[i];
      float3 pj_pos = (i + 1 == p) ? pos_p : positions[i + 1];
      float vx = pi_pos.x - pj_pos.x;
      float vy = pi_pos.y - pj_pos.y;
      float vz = pi_pos.z - pj_pos.z;
      if (have_prev && i > p - 2 && i > 0) {
        float la = sqrtf(vx*vx + vy*vy + vz*vz);
        float lb = sqrtf(prev_vx*prev_vx + prev_vy*prev_vy + prev_vz*prev_vz);
        if (la > 1e-6f && lb > 1e-6f) {
          float dot = (vx*prev_vx + vy*prev_vy + vz*prev_vz) / (la * lb);
          if (dot < -1.0f) dot = -1.0f;
          if (dot >  1.0f) dot =  1.0f;
          float ang = (1.0f - dot) * 0.5f;
          scb += ang * ang * ang * s.springAngularConstant;
        }
      }
      prev_vx = vx; prev_vy = vy; prev_vz = vz;
      have_prev = true;
    }
  }

  /* Sphere boundary penalty for the moved bead. */
  sca += sphereBoundaryPenalty(pos_p, s.sphere_center, s.sphere_radius);

  return sca + scb;
}

/* Full energy for initial score computation (float precision for GPU perf) */
__device__ float calcFullEnergySmooth(
    const float3 * __restrict__ positions,
    const float  * __restrict__ dist_to_next,
    int active_n,
    const smooth_gpu_settings &s)
{
  float sca = 0.0f, scb = 0.0f;

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
    sca += diff * diff * (diff >= 0.0f ? s.springConstantStretch
                                       : s.springConstantSqueeze);
  }

  {
    float prev_vx = 0.0f, prev_vy = 0.0f, prev_vz = 0.0f;
    for (int i = 0; i + 1 < active_n; ++i) {
      float3 pi = positions[i], pj = positions[i + 1];
      float vx = pi.x - pj.x;
      float vy = pi.y - pj.y;
      float vz = pi.z - pj.z;
      if (i > 0) {
        float la = sqrtf(vx*vx + vy*vy + vz*vz);
        float lb = sqrtf(prev_vx*prev_vx + prev_vy*prev_vy + prev_vz*prev_vz);
        if (la > 1e-6f && lb > 1e-6f) {
          float dot = (vx*prev_vx + vy*prev_vy + vz*prev_vz) / (la * lb);
          if (dot < -1.0f) dot = -1.0f;
          if (dot >  1.0f) dot =  1.0f;
          float ang = (1.0f - dot) * 0.5f;
          scb += ang * ang * ang * s.springAngularConstant;
        }
      }
      prev_vx = vx; prev_vy = vy; prev_vz = vz;
    }
  }

  /* Sphere boundary penalty for all beads. */
  float sphere_penalty = 0.0f;
  for (int i = 0; i < active_n; ++i) {
    sphere_penalty += sphereBoundaryPenalty(positions[i], s.sphere_center,
                                            s.sphere_radius);
  }

  return sca * s.weightDistSmooth + scb * s.weightAngleSmooth + sphere_penalty;
}

/* =========================================================================
 * Setup kernel — initialise one curandState per thread
 * ======================================================================= */
#define FULL_MASK 0xffffffff

__global__ void setupSmoothKernel(curandState *state, unsigned long long seed) {
  int tid = blockDim.x * blockIdx.x + threadIdx.x;
  curand_init(seed, tid, 0, &state[tid]);
}

/* =========================================================================
 * MonteCarloSmoothKernelMultiWarp
 *
 * Each warp runs an independent SA trajectory on its own position copy.
 * Only lane 0 of each warp drives the SA (energy is local, no parallel
 * evaluation benefit). Each warp converges independently.
 *
 * per_warp_pos: [n_warps * active_n] — each warp's position array
 * d_best_score: [n_warps] — best score achieved by each warp
 * ======================================================================= */
__global__ void MonteCarloSmoothKernelMultiWarp(
    curandState   * __restrict__ state,
    float3        * __restrict__ per_warp_pos,  /* [n_warps * active_n] */
    const bool    * __restrict__ is_fixed,
    const float   * __restrict__ dist_to_next,
    float                        T_init,
    float                        step_size,
    int                          active_n,
    smooth_gpu_settings          s,
    float         * __restrict__ d_best_score)   /* [n_warps] output */
{
  const int tid    = blockDim.x * blockIdx.x + threadIdx.x;
  const int lane   = threadIdx.x & 31;
  const int warpId = tid >> 5;

  /* Only lane 0 does useful work; other lanes return early */
  if (lane != 0) return;

  /* Per-warp position array */
  float3 *my_pos = per_warp_pos + (size_t)warpId * active_n;

  curandState localState = state[tid];

  /* Per-warp SA state (float precision — 64x faster on consumer GPUs) */
  float T            = T_init;
  float score_prev   = calcFullEnergySmooth(my_pos, dist_to_next, active_n, s);
  float step         = step_size;
  int   iterations   = 0;
  int   total_accepted = 0;
  int   milestone_successes = 0;
  int   improvement_misses  = 0;
  float milestone_score     = score_prev;
  float best_score          = score_prev;

  const int max_iters = s.MCstopConditionStepsSmooth;

  while (iterations < max_iters) {

    for (int inner = 0; inner < 256; ++inner) {

      /* Pick a random non-fixed bead */
      int chosen_p = -1;
      for (int attempt = 0; attempt < active_n; ++attempt) {
        int p = (int)(curand_uniform(&localState) * (float)active_n);
        if (p >= active_n) p = active_n - 1;
        if (!is_fixed[p]) {
          chosen_p = p;
          break;
        }
      }
      if (chosen_p < 0) continue;

      float3 old_pos = my_pos[chosen_p];
      float3 new_pos = old_pos;
      new_pos.x += (2.0f * curand_uniform(&localState) - 1.0f) * step;
      new_pos.y += (2.0f * curand_uniform(&localState) - 1.0f) * step;
      if (!s.use2D)
        new_pos.z += (2.0f * curand_uniform(&localState) - 1.0f) * step;

      float e_old = calcLocalEnergySmooth(chosen_p, old_pos,
                                           my_pos, dist_to_next, active_n, s);
      float e_new = calcLocalEnergySmooth(chosen_p, new_pos,
                                           my_pos, dist_to_next, active_n, s);
      float new_total = score_prev - e_old + e_new;

      bool accepted = false;
      if (new_total < score_prev) {
        accepted = true;
      } else if (T > 0.0f) {
        float tp = s.tempJumpScaleSmooth *
                   expf(-s.tempJumpCoefSmooth *
                        (new_total / fmaxf(score_prev, 1e-12f)) / T);
        accepted = (curand_uniform(&localState) < tp);
      }

      if (accepted) {
        my_pos[chosen_p] = new_pos;
        score_prev = new_total;
        ++milestone_successes;
        ++total_accepted;
        if (score_prev < best_score) best_score = score_prev;
      }

      T *= s.dtTempSmooth;
    } /* inner loop */

    iterations += 256;

    /* Milestone check — each warp converges independently */
    if (iterations % s.MCstopConditionStepsSmooth < 256) {
      bool no_improvement =
          (score_prev >
               s.MCstopConditionImprovementSmooth * milestone_score &&
           milestone_successes < s.MCstopConditionMinSuccessesSmooth) ||
          score_prev < 1e-6f;

      if (no_improvement)
        ++improvement_misses;

      if (improvement_misses >= s.milestoneFailsThreshold || score_prev < 1e-6f)
        break;

      milestone_score     = score_prev;
      milestone_successes = 0;
    }

  } /* outer while */

  /* Write final score for this warp */
  d_best_score[warpId] = best_score;

  /* Record iteration count from warp 0 */
  if (warpId == 0) {
    d_smooth_total_iterations = iterations;
    d_smooth_total_accepted = total_accepted;
  }

  state[tid] = localState;
}

/* =========================================================================
 * Host entry point: LooperSolver::ParallelMonteCarloSmooth
 *
 * Multi-warp version: launches sm_count warps, each running an independent
 * SA trajectory. Picks the best result.
 * ======================================================================= */
float LooperSolver::ParallelMonteCarloSmooth(float step_size,
                                              unsigned long long seed) {
  const int active_n = static_cast<int>(active_region.size());
  if (active_n <= 2)
    return 0.0f;

  /* Initialize GPU cache; fall back to CPU if no CUDA device. */
  if (!gpu_cache.initialized) {
    gpu_cache.init();
    if (!gpu_cache.initialized)
      return static_cast<float>(MonteCarloArcsSmooth(step_size, false));
  }

  const int threads_per_block = 32; /* 1 warp per block */
  /* Use the memory-capped warp count from initSmoothBuffers if available,
   * otherwise fall back to sm_count. This prevents OOM on large chromosomes. */
  int blocks = (gpu_cache.smooth_n_warps > 0) ? gpu_cache.smooth_n_warps : gpu_cache.sm_count;
  if (blocks < 1) blocks = 1;
  int n_warps = blocks; /* 1 warp per block */

  /* Ensure smooth curandState buffer is large enough */
  int total_threads = threads_per_block * blocks;
  if (total_threads > gpu_cache.smooth_states_count) {
    /* Reallocate if needed */
    if (gpu_cache.d_smooth_states)
      cudaFree(gpu_cache.d_smooth_states);
    curandState *ptr = nullptr;
    cudaMalloc(&ptr, (size_t)total_threads * sizeof(curandState));
    gpu_cache.d_smooth_states = ptr;
    gpu_cache.smooth_states_count = total_threads;
    gpu_cache.smooth_states_seeded = false;
  }

  /* Build host-side flat arrays */
  std::vector<float3> h_positions(active_n);
  std::vector<char>   h_fixed(active_n);
  std::vector<float>  h_dist(active_n);

  for (int i = 0; i < active_n; ++i) {
    int ci = active_region[i];
    h_positions[i].x = clusters[ci].pos.x;
    h_positions[i].y = clusters[ci].pos.y;
    h_positions[i].z = clusters[ci].pos.z;
    h_fixed[i]       = clusters[ci].is_fixed ? 1 : 0;
    h_dist[i]        = static_cast<float>(clusters[ci].dist_to_next);
  }

  /* Upload shared data (is_fixed, dist_to_next) — also allocates warp buffers */
  gpu_cache.initSmoothBuffers(active_n);

  /* Use cached per-warp buffers (must be after initSmoothBuffers) */
  float3 *d_warp_positions = static_cast<float3 *>(gpu_cache.d_smooth_warp_positions);
  float  *d_warp_scores    = static_cast<float *>(gpu_cache.d_smooth_warp_scores);
  bool  *d_fixed_ptr = static_cast<bool *>(gpu_cache.d_smooth_fixed);
  float *d_dist_ptr  = static_cast<float *>(gpu_cache.d_smooth_dist);
  cudaMemcpy(d_fixed_ptr, h_fixed.data(), active_n * sizeof(bool), cudaMemcpyHostToDevice);
  cudaMemcpy(d_dist_ptr,  h_dist.data(),  active_n * sizeof(float), cudaMemcpyHostToDevice);

  /* Initialize all per-warp position copies from input:
   * Host→Device for warp 0, then Device→Device for remaining warps */
  cudaMemcpy(d_warp_positions, h_positions.data(),
             active_n * sizeof(float3), cudaMemcpyHostToDevice);
  for (int w = 1; w < n_warps; ++w) {
    cudaMemcpy(d_warp_positions + (size_t)w * active_n,
               d_warp_positions, active_n * sizeof(float3),
               cudaMemcpyDeviceToDevice);
  }

  /* Seed RNG only once */
  curandState *d_states = static_cast<curandState *>(gpu_cache.d_smooth_states);
  if (!gpu_cache.smooth_states_seeded) {
    setupSmoothKernel<<<blocks, threads_per_block>>>(d_states, seed);
    gpuErrchkSmooth(cudaDeviceSynchronize());
    gpu_cache.smooth_states_seeded = true;
  }

  /* Pack settings */
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
  s.sphere_center                    = make_float3(0.0f, 0.0f, 0.0f);
  s.sphere_radius                    = Settings::sphere_radius;

  float T_init = static_cast<float>(Settings::maxTempSmooth);

  /* Launch multi-warp kernel */
  MonteCarloSmoothKernelMultiWarp<<<blocks, threads_per_block>>>(
      d_states,
      d_warp_positions,
      d_fixed_ptr,
      d_dist_ptr,
      T_init,
      step_size,
      active_n,
      s,
      d_warp_scores);

  gpuErrchkSmooth(cudaPeekAtLastError());

  /* Find the best warp */
  std::vector<float> h_scores(n_warps);
  cudaMemcpy(h_scores.data(), d_warp_scores,
             n_warps * sizeof(float), cudaMemcpyDeviceToHost);

  int best_warp = 0;
  float best_sc = h_scores[0];
  for (int w = 1; w < n_warps; ++w) {
    if (h_scores[w] < best_sc) {
      best_sc = h_scores[w];
      best_warp = w;
    }
  }

  /* Copy best warp's positions back to CPU */
  cudaMemcpy(h_positions.data(),
             d_warp_positions + (size_t)best_warp * active_n,
             active_n * sizeof(float3), cudaMemcpyDeviceToHost);
  for (int i = 0; i < active_n; ++i) {
    if (h_fixed[i] == 0) {
      clusters[active_region[i]].pos.x = h_positions[i].x;
      clusters[active_region[i]].pos.y = h_positions[i].y;
      clusters[active_region[i]].pos.z = h_positions[i].z;
    }
  }

  /* Retrieve iteration and acceptance counts */
  int gpu_iters = 0, gpu_accepted = 0;
  cudaMemcpyFromSymbol(&gpu_iters, d_smooth_total_iterations, sizeof(int));
  cudaMemcpyFromSymbol(&gpu_accepted, d_smooth_total_accepted, sizeof(int));
  total_mc_steps += gpu_iters;

  if (sim_logger)
    sim_logger->logAcceptanceRatio("smooth_gpu", gpu_iters, gpu_accepted);

  /* Compute final CPU-side energy (authoritative) */
  double final_score = calcScoreStructureSmooth(true, true);

  return static_cast<float>(final_score);
}
