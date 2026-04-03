#!/usr/bin/env bash
# pMMC Resource Utilisation Monitor
# Captures GPU util, VRAM, system RAM, and CPU util during a run.
# Usage:  bash run_resource_monitor.sh [output_dir] [chromosome]
set -euo pipefail

OUTDIR="${1:-./resource_monitor}"
CHR="${2:-chr1}"
PMMC_GENERATE="${PMMC_GENERATE:-pMMC-generate}"
PMMC_EXE="${PMMC_EXE:-pMMC-reconstruct}"
SEED=42
CSV="$OUTDIR/resource_utilisation.csv"

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

mkdir -p "$OUTDIR"

# Generate synthetic data
SYNTH="$OUTDIR/synth_${CHR}/"
if [ ! -f "$SYNTH/settings.ini" ]; then
  mkdir -p "$SYNTH"
  write_generate_ini "$SYNTH/generate.ini" "$SYNTH" "$CHR" 100 1
  "$PMMC_GENERATE" -s "$SYNTH/generate.ini"
fi

# Start resource monitor in background
MONITOR_LOG="$OUTDIR/resource_samples.csv"
echo "timestamp,gpu_util_pct,gpu_mem_used_mb,gpu_mem_total_mb,sys_mem_used_mb,cpu_util_pct" > "$MONITOR_LOG"

monitor_resources() {
  while true; do
    TS=$(date +%s.%N)
    GPU_LINE=$(nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits 2>/dev/null || echo "0, 0, 0")
    GPU_UTIL=$(echo "$GPU_LINE" | awk -F', ' '{print $1}')
    GPU_MEM=$(echo "$GPU_LINE" | awk -F', ' '{print $2}')
    GPU_TOTAL=$(echo "$GPU_LINE" | awk -F', ' '{print $3}')
    if command -v free &>/dev/null; then
      SYS_MEM=$(free -m | awk '/Mem:/{print $3}')
    else
      SYS_MEM=0
    fi
    if [ -f /proc/stat ]; then
      CPU_UTIL=$(awk '/^cpu /{u=$2+$4; t=$2+$4+$5; printf "%.1f", u/t*100}' /proc/stat)
    else
      CPU_UTIL=0
    fi
    echo "$TS,$GPU_UTIL,$GPU_MEM,$GPU_TOTAL,$SYS_MEM,$CPU_UTIL" >> "$MONITOR_LOG"
    sleep 1
  done
}

monitor_resources &
MONITOR_PID=$!
trap "kill $MONITOR_PID 2>/dev/null || true" EXIT

# Run pMMC reconstruction
D="$OUTDIR/pmmc_${CHR}/"; mkdir -p "$D"
write_reconstruct_ini "$D/settings.ini" "$SYNTH" "$D" res "$CHR" $SEED
T0=$(date +%s%N)
"$PMMC_EXE" -s "$D/settings.ini"
T1=$(date +%s%N)
WALL_TIME=$(echo "scale=3;($T1-$T0)/1000000000"|bc)

kill $MONITOR_PID 2>/dev/null || true
sleep 1

# Compute summary
echo "metric,value" > "$CSV"
echo "chromosome,$CHR" >> "$CSV"
echo "wall_time_s,$WALL_TIME" >> "$CSV"
awk -F, 'NR>1 { if($2+0>mg) mg=$2; if($3+0>mm) mm=$3; gt=$4; if($5+0>sm) sm=$5; if($6+0>mc) mc=$6; n++ }
END { printf "gpu_util_max_pct,%.1f\ngpu_vram_max_mb,%.0f\ngpu_vram_total_mb,%.0f\nsys_ram_max_mb,%.0f\ncpu_util_max_pct,%.1f\nsamples,%d\n",mg,mm,gt,sm,mc,n }' "$MONITOR_LOG" >> "$CSV"

echo "Resource utilisation:"
cat "$CSV"
