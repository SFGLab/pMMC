#!/usr/bin/env bash
# pMMC Scalability Benchmark — wall time vs loop count
# Tools run SEQUENTIALLY to avoid resource competition.
# Usage:  bash run_scalability.sh [output_dir]
set -euo pipefail

OUTDIR="${1:-./scalability_benchmark}"
PMMC_GENERATE="${PMMC_GENERATE:-pMMC-generate}"
PMMC_EXE="${PMMC_EXE:-pMMC-reconstruct}"
CUDAMMC_EXE="${CUDAMMC_EXE:-cudaMMC}"
GNOME_EXE="${GNOME_EXE:-3dnome}"
CHR="chr22"
SEED=42
CSV="$OUTDIR/scalability_results.csv"

write_generate_ini() {
  cat > "$1" <<GENEOF
[main]
output_dir = $2
chromosomes = $3
label = synthetic
num_loops = $4
ensemble_size = ${5:-1}
resolution = 25000
seed = $SEED
GENEOF
}

write_reconstruct_ini() {
  cat > "$1" <<INIEOF
[main]
output_level = 1
random_walk = no
loop_density = 5
chromosomes = $5
output_dir = $3
label = $4
seed = $6
ensemble_size = ${7:-1}
output_format = pdb
max_level = 5
steps_lvl1 = 2
steps_lvl2 = 2
steps_arcs = 3
steps_smooth = 3
noise_lvl1 = 0.5
noise_lvl2 = 0.5
noise_arcs = 0.01
noise_smooth = 5.0
max_pet_length = 1000000
long_pet_power = 2.0
long_pet_scale = 1.0
[cuda]
num_threads = 512
blocks_multiplier = 4
milestone_fails = 3
[data]
data_dir = $2
anchors = anchors.txt
clusters = clusters.txt
factors = synthetic
singletons = singletons.txt
split_singleton_files_by_chr = no
singletons_inter =
segment_split = ${2}segments.bed
centromeres = ${2}centromeres.bed
[distance]
genomic_dist_power = 0.5
genomic_dist_scale = 1.0
genomic_dist_base = 0.0
freq_dist_scale = 25.0
freq_dist_power = -0.6
freq_dist_scale_inter = 120.0
freq_dist_power_inter = -1.0
count_dist_a = 0.2
count_dist_scale = 1.8
count_dist_shift = 8
count_dist_base_level = 0.2
[template]
template_scale = 7.0
dist_heatmap_scale = 15.0
[motif_orientation]
use_motif_orientation = no
weight = 50.0
[springs]
stretch_constant = 0.1
squeeze_constant = 0.1
angular_constant = 0.1
stretch_constant_arcs = 1.0
squeeze_constant_arcs = 1.0
[simulation_heatmap]
max_temp_heatmap = 5.0
delta_temp_heatmap = 0.9999
jump_temp_scale_heatmap = 50.0
jump_temp_coef_heatmap = 20.0
stop_condition_improvement_threshold_heatmap = 0.99
stop_condition_successes_threshold_heatmap = 10
stop_condition_steps_heatmap = 50000
[simulation_arcs]
max_temp = 5.0
delta_temp = 0.9999
jump_temp_scale = 50.0
jump_temp_coef = 20.0
stop_condition_improvement_threshold = 0.975
stop_condition_successes_threshold = 100
stop_condition_steps = 50000
[simulation_arcs_smooth]
dist_weight = 1.0
angle_weight = 1.0
max_temp = 5.0
delta_temp = 0.9999
jump_temp_scale = 50.0
jump_temp_coef = 20.0
stop_condition_improvement_threshold = 0.99
stop_condition_successes_threshold = 50
stop_condition_steps = 50000
INIEOF
}

LOOPS=(10 25 50 100 250 500 1000)

mkdir -p "$OUTDIR"
echo "tool,num_loops,wall_time_seconds" > "$CSV"

for NL in "${LOOPS[@]}"; do
  echo "=== $NL loops ==="
  SYNTH="$OUTDIR/synth_${NL}/"; mkdir -p "$SYNTH"
  write_generate_ini "$SYNTH/generate.ini" "$SYNTH" "$CHR" "$NL" 1
  "$PMMC_GENERATE" -s "$SYNTH/generate.ini"

  if command -v "$GNOME_EXE" &>/dev/null; then
    D="$OUTDIR/gnome_${NL}/"; mkdir -p "$D"
    write_reconstruct_ini "$D/settings.ini" "$SYNTH" "$D" g "$CHR" $SEED
    T0=$(date +%s%N)
    "$GNOME_EXE" -a create -s "$D/settings.ini" -c "$CHR" -o "$D" -n g || true
    T1=$(date +%s%N); echo "3D-GNOME,$NL,$(echo "scale=3;($T1-$T0)/1000000000"|bc)" >> "$CSV"
  fi
  if command -v "$CUDAMMC_EXE" &>/dev/null; then
    D="$OUTDIR/cuda_${NL}/"; mkdir -p "$D"
    write_reconstruct_ini "$D/settings.ini" "$SYNTH" "$D" c "$CHR" $SEED
    T0=$(date +%s%N)
    "$CUDAMMC_EXE" -a create -s "$D/settings.ini" -c "$CHR" -o "$D" -n c || true
    T1=$(date +%s%N); echo "cudaMMC,$NL,$(echo "scale=3;($T1-$T0)/1000000000"|bc)" >> "$CSV"
  fi
  D="$OUTDIR/pmmc_${NL}/"; mkdir -p "$D"
  write_reconstruct_ini "$D/settings.ini" "$SYNTH" "$D" p "$CHR" $SEED
  T0=$(date +%s%N)
  "$PMMC_EXE" -s "$D/settings.ini" || true
  T1=$(date +%s%N); echo "pMMC,$NL,$(echo "scale=3;($T1-$T0)/1000000000"|bc)" >> "$CSV"
done

echo "Scalability benchmark done.  Results: $CSV"
cat "$CSV"
