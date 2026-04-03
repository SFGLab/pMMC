#!/usr/bin/env bash
# pMMC Quality Comparison — run after run_benchmark.sh
# Usage:  bash run_quality_comparison.sh [benchmark_dir]
set -euo pipefail

B="${1:-./benchmark_results}"
PMMC_ANALYZE="${PMMC_ANALYZE:-pMMC-analyze}"
CHR="chr22"
Q="$B/quality/"; mkdir -p "$Q"

P="$B/pmmc/loops_pmmc.hcm"
C="$B/cudammc/loops_cuda.hcm"
G="$B/3dgnome/loops_gnome.hcm"

echo "pair,scc" > "$Q/quality_summary.csv"

for HCM in "$P" "$C" "$G"; do
  if [ -f "$HCM" ]; then
    NAME=$(basename "$HCM" .hcm)
    cat > "$Q/distmap_${NAME}.ini" <<DMEOF
[main]
action = distmap
output_dir = $Q
label = $NAME
chromosomes = $CHR
level = 2
resolution = 25000
[data]
input_file = $HCM
DMEOF
    "$PMMC_ANALYZE" -s "$Q/distmap_${NAME}.ini" || true
  fi
done

if [ -f "$P" ] && [ -f "$C" ]; then
  cat > "$Q/pmmc_vs_cuda.ini" <<MCEOF
[main]
action = metrics
output_dir = $Q
label = pmmc_vs_cuda
chromosomes = $CHR
level = 2
resolution = 25000
input_file_2 = $C
[data]
input_file = $P
MCEOF
  "$PMMC_ANALYZE" -s "$Q/pmmc_vs_cuda.ini" | tee "$Q/pmmc_vs_cuda.log"
  SCC=$(grep -oP 'SCC.*?\K[0-9.]+' "$Q/pmmc_vs_cuda.log" 2>/dev/null || echo N/A)
  echo "pmmc_vs_cuda,$SCC" >> "$Q/quality_summary.csv"
fi
if [ -f "$P" ] && [ -f "$G" ]; then
  cat > "$Q/pmmc_vs_gnome.ini" <<MGEOF
[main]
action = metrics
output_dir = $Q
label = pmmc_vs_gnome
chromosomes = $CHR
level = 2
resolution = 25000
input_file_2 = $G
[data]
input_file = $P
MGEOF
  "$PMMC_ANALYZE" -s "$Q/pmmc_vs_gnome.ini" | tee "$Q/pmmc_vs_gnome.log"
  SCC=$(grep -oP 'SCC.*?\K[0-9.]+' "$Q/pmmc_vs_gnome.log" 2>/dev/null || echo N/A)
  echo "pmmc_vs_gnome,$SCC" >> "$Q/quality_summary.csv"
fi

echo "Quality comparison done."
cat "$Q/quality_summary.csv"
