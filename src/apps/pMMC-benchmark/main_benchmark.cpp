#include <AppCommon.hpp>
#include <Static.hpp>
#include <Settings.hpp>
#include <BenchmarkRunner.hpp>

using namespace std;

void usage() {
  printf("pMMC-benchmark — benchmarking and script generation\n\n");
  printf("Usage: pMMC-benchmark -s <settings.ini>\n\n");
  printf("All parameters are controlled by the settings file.\n\n");
  printf("Key settings.ini fields in [main]:\n");
  printf("  output_dir      Output directory (default: ./benchmark/)\n");
  printf("  chromosomes     Chromosome (default: chr22)\n");
  printf("  label           Label (default: bench)\n");
  printf("  ensemble_size   Ensemble size (default: 20)\n");
  printf("  seed            Random seed (0=time-based)\n");
}

// Write a benchmark bash script to the given path.
static void writeBenchmarkScript(const std::string &path,
                                 const std::string &content) {
  FILE *f = fopen(path.c_str(), "w");
  if (!f) {
    printf("[benchmark] ERROR: cannot write %s\n", path.c_str());
    return;
  }
  fprintf(f, "%s", content.c_str());
  fclose(f);
  printf("[benchmark] Wrote %s\n", path.c_str());
}

void prepareBenchmark() {
  std::string outdir = Settings::outputDir.empty() ? "./benchmark/" : Settings::outputDir;
  std::string chr = Settings::chromosomes.empty() ? "chr22" : Settings::chromosomes;
  int ens = Settings::ensembleSize;

  portable_mkdir(outdir.c_str());

  std::string name = Settings::label.empty() ? "bench" : Settings::label;
  int chr_length = 51304566;
  BenchmarkRunner bench;
  bench.setChromosome(chr, chr_length);
  bench.setResolution(25000);
  bench.setEnsembleSize(ens);
  bench.setOutputDir(outdir);
  bench.setLabel(name);
  bench.run();

  // ── Generate redistributable bash scripts ──────────────────────────
  // All tools (3D-GNOME, cudaMMC, pMMC) run SEQUENTIALLY — never
  // simultaneously — so they do not compete for GPU/memory resources.
  //
  // pMMC is split into sub-apps (pMMC-generate, pMMC-reconstruct,
  // pMMC-analyze) that each take -s <settings.ini>.  cudaMMC and
  // 3D-GNOME still use the old -a/-c/-o/-n CLI.  Each script embeds
  // helper functions that write the required INI files at runtime.

  std::string ens_str = ftext("%d", ens);

  // Shared shell fragment: write_reconstruct_ini function body
  // (embedded into each script that needs it)
  auto iniBody = []() -> std::string {
    std::string b;
    b += "  cat > \"$1\" <<INIEOF\n";
    b += "[main]\n";
    b += "output_level = 1\n";
    b += "random_walk = no\n";
    b += "loop_density = 5\n";
    b += "chromosomes = $5\n";
    b += "output_dir = $3\n";
    b += "label = $4\n";
    b += "seed = $6\n";
    b += "ensemble_size = ${7:-1}\n";
    b += "output_format = pdb\n";
    b += "max_level = 5\n";
    b += "steps_lvl1 = 2\n";
    b += "steps_lvl2 = 2\n";
    b += "steps_arcs = 3\n";
    b += "steps_smooth = 3\n";
    b += "noise_lvl1 = 0.5\n";
    b += "noise_lvl2 = 0.5\n";
    b += "noise_arcs = 0.01\n";
    b += "noise_smooth = 5.0\n";
    b += "max_pet_length = 1000000\n";
    b += "long_pet_power = 2.0\n";
    b += "long_pet_scale = 1.0\n";
    b += "[cuda]\n";
    b += "num_threads = 512\n";
    b += "blocks_multiplier = 4\n";
    b += "milestone_fails = 3\n";
    b += "[data]\n";
    b += "data_dir = $2\n";
    b += "anchors = anchors.txt\n";
    b += "clusters = clusters.txt\n";
    b += "factors = synthetic\n";
    b += "singletons = singletons.txt\n";
    b += "split_singleton_files_by_chr = no\n";
    b += "singletons_inter =\n";
    b += "segment_split = ${2}segments.bed\n";
    b += "centromeres = ${2}centromeres.bed\n";
    b += "[distance]\n";
    b += "genomic_dist_power = 0.5\n";
    b += "genomic_dist_scale = 1.0\n";
    b += "genomic_dist_base = 0.0\n";
    b += "freq_dist_scale = 25.0\n";
    b += "freq_dist_power = -0.6\n";
    b += "freq_dist_scale_inter = 120.0\n";
    b += "freq_dist_power_inter = -1.0\n";
    b += "count_dist_a = 0.2\n";
    b += "count_dist_scale = 1.8\n";
    b += "count_dist_shift = 8\n";
    b += "count_dist_base_level = 0.2\n";
    b += "[template]\n";
    b += "template_scale = 7.0\n";
    b += "dist_heatmap_scale = 15.0\n";
    b += "[motif_orientation]\n";
    b += "use_motif_orientation = no\n";
    b += "weight = 50.0\n";
    b += "[springs]\n";
    b += "stretch_constant = 0.1\n";
    b += "squeeze_constant = 0.1\n";
    b += "angular_constant = 0.1\n";
    b += "stretch_constant_arcs = 1.0\n";
    b += "squeeze_constant_arcs = 1.0\n";
    b += "[simulation_heatmap]\n";
    b += "max_temp_heatmap = 5.0\n";
    b += "delta_temp_heatmap = 0.9999\n";
    b += "jump_temp_scale_heatmap = 50.0\n";
    b += "jump_temp_coef_heatmap = 20.0\n";
    b += "stop_condition_improvement_threshold_heatmap = 0.99\n";
    b += "stop_condition_successes_threshold_heatmap = 10\n";
    b += "stop_condition_steps_heatmap = 50000\n";
    b += "[simulation_arcs]\n";
    b += "max_temp = 5.0\n";
    b += "delta_temp = 0.9999\n";
    b += "jump_temp_scale = 50.0\n";
    b += "jump_temp_coef = 20.0\n";
    b += "stop_condition_improvement_threshold = 0.975\n";
    b += "stop_condition_successes_threshold = 100\n";
    b += "stop_condition_steps = 50000\n";
    b += "[simulation_arcs_smooth]\n";
    b += "dist_weight = 1.0\n";
    b += "angle_weight = 1.0\n";
    b += "max_temp = 5.0\n";
    b += "delta_temp = 0.9999\n";
    b += "jump_temp_scale = 50.0\n";
    b += "jump_temp_coef = 20.0\n";
    b += "stop_condition_improvement_threshold = 0.99\n";
    b += "stop_condition_successes_threshold = 50\n";
    b += "stop_condition_steps = 50000\n";
    b += "INIEOF\n";
    return b;
  };

  // Shared shell fragment: helper functions used across scripts
  auto helpers = [&]() -> std::string {
    std::string h;
    // write_generate_ini <ini> <outdir> <chr> <loops> <ens>
    h += "write_generate_ini() {\n";
    h += "  cat > \"$1\" <<GENEOF\n";
    h += "[main]\n";
    h += "output_dir = $2\n";
    h += "chromosomes = $3\n";
    h += "label = synthetic\n";
    h += "num_loops = $4\n";
    h += "ensemble_size = ${5:-1}\n";
    h += "resolution = 25000\n";
    h += "seed = $SEED\n";
    h += "GENEOF\n";
    h += "}\n\n";
    // write_reconstruct_ini <ini> <datadir> <outdir> <label> <chr> <seed> [ens]
    h += "write_reconstruct_ini() {\n";
    h += iniBody();
    h += "}\n\n";
    return h;
  };

  // Shared shell fragment: env-var header for all scripts
  auto envHeader = [](bool needAnalyze = false) -> std::string {
    std::string h;
    h += "PMMC_GENERATE=\"${PMMC_GENERATE:-pMMC-generate}\"\n";
    h += "PMMC_EXE=\"${PMMC_EXE:-pMMC-reconstruct}\"\n";
    if (needAnalyze)
      h += "PMMC_ANALYZE=\"${PMMC_ANALYZE:-pMMC-analyze}\"\n";
    h += "CUDAMMC_EXE=\"${CUDAMMC_EXE:-cudaMMC}\"\n";
    h += "GNOME_EXE=\"${GNOME_EXE:-3dnome}\"\n";
    return h;
  };

  // Script 1: run_benchmark.sh — Full 3-tool comparison on synthetic data
  {
    std::string s;
    s += "#!/usr/bin/env bash\n";
    s += "# pMMC Benchmark — 3-tool sequential comparison\n";
    s += "# Generated by: pMMC-benchmark\n";
    s += "# Runs 3D-GNOME, cudaMMC, and pMMC-reconstruct ONE AFTER ANOTHER on\n";
    s += "# the same synthetic data so they never compete for resources.\n";
    s += "#\n";
    s += "# Usage:  bash run_benchmark.sh [output_dir]\n";
    s += "#\n";
    s += "# Environment variables (override paths):\n";
    s += "#   PMMC_GENERATE — pMMC-generate     (default: pMMC-generate)\n";
    s += "#   PMMC_EXE      — pMMC-reconstruct  (default: pMMC-reconstruct)\n";
    s += "#   PMMC_ANALYZE  — pMMC-analyze       (default: pMMC-analyze)\n";
    s += "#   CUDAMMC_EXE   — cudaMMC            (default: cudaMMC, optional)\n";
    s += "#   GNOME_EXE     — 3D-GNOME           (default: 3dnome, optional)\n";
    s += "set -euo pipefail\n\n";
    s += "OUTDIR=\"${1:-./benchmark_results}\"\n";
    s += envHeader(true);
    s += "CHR=\"" + chr + "\"\n";
    s += "SEED=42\n";
    s += "ENSEMBLE=" + ens_str + "\n";
    s += "CSV=\"$OUTDIR/benchmark_results.csv\"\n\n";
    s += helpers();
    s += "mkdir -p \"$OUTDIR\"\n";
    s += "echo \"tool,chromosome,wall_time_seconds,exit_code\" > \"$CSV\"\n\n";
    // Step 1: generate
    s += "# Step 1 — generate synthetic data\n";
    s += "echo \"[1/5] Generating synthetic data...\"\n";
    s += "SYNTH=\"$OUTDIR/synthetic/\"; mkdir -p \"$SYNTH\"\n";
    s += "write_generate_ini \"$OUTDIR/generate.ini\" \"$SYNTH\" \"$CHR\" 100 \"$ENSEMBLE\"\n";
    s += "\"$PMMC_GENERATE\" -s \"$OUTDIR/generate.ini\"\n\n";
    // Step 2: 3D-GNOME
    s += "# Step 2 — 3D-GNOME (sequential)\n";
    s += "if command -v \"$GNOME_EXE\" &>/dev/null; then\n";
    s += "  echo \"[2/5] Running 3D-GNOME...\"\n";
    s += "  D=\"$OUTDIR/3dgnome/\"; mkdir -p \"$D\"\n";
    s += "  write_reconstruct_ini \"$D/settings.ini\" \"$SYNTH\" \"$D\" gnome \"$CHR\" $SEED\n";
    s += "  T0=$(date +%s%N)\n";
    s += "  \"$GNOME_EXE\" -a create -s \"$D/settings.ini\" -c \"$CHR\" -o \"$D\" -n gnome && RC=0 || RC=$?\n";
    s += "  T1=$(date +%s%N); T=$(echo \"scale=3;($T1-$T0)/1000000000\"|bc)\n";
    s += "  echo \"3D-GNOME,$CHR,$T,$RC\" >> \"$CSV\"\n";
    s += "else\n  echo \"[2/5] 3D-GNOME not found, skipping.\"\n";
    s += "  echo \"3D-GNOME,$CHR,,skipped\" >> \"$CSV\"\nfi\n\n";
    // Step 3: cudaMMC
    s += "# Step 3 — cudaMMC (sequential)\n";
    s += "if command -v \"$CUDAMMC_EXE\" &>/dev/null; then\n";
    s += "  echo \"[3/5] Running cudaMMC...\"\n";
    s += "  D=\"$OUTDIR/cudammc/\"; mkdir -p \"$D\"\n";
    s += "  write_reconstruct_ini \"$D/settings.ini\" \"$SYNTH\" \"$D\" cuda \"$CHR\" $SEED\n";
    s += "  T0=$(date +%s%N)\n";
    s += "  \"$CUDAMMC_EXE\" -a create -s \"$D/settings.ini\" -c \"$CHR\" -o \"$D\" -n cuda && RC=0 || RC=$?\n";
    s += "  T1=$(date +%s%N); T=$(echo \"scale=3;($T1-$T0)/1000000000\"|bc)\n";
    s += "  echo \"cudaMMC,$CHR,$T,$RC\" >> \"$CSV\"\n";
    s += "else\n  echo \"[3/5] cudaMMC not found, skipping.\"\n";
    s += "  echo \"cudaMMC,$CHR,,skipped\" >> \"$CSV\"\nfi\n\n";
    // Step 4: pMMC
    s += "# Step 4 — pMMC-reconstruct (sequential)\n";
    s += "echo \"[4/5] Running pMMC...\"\n";
    s += "D=\"$OUTDIR/pmmc/\"; mkdir -p \"$D\"\n";
    s += "write_reconstruct_ini \"$D/settings.ini\" \"$SYNTH\" \"$D\" pmmc \"$CHR\" $SEED\n";
    s += "T0=$(date +%s%N)\n";
    s += "\"$PMMC_EXE\" -s \"$D/settings.ini\" && RC=0 || RC=$?\n";
    s += "T1=$(date +%s%N); T=$(echo \"scale=3;($T1-$T0)/1000000000\"|bc)\n";
    s += "echo \"pMMC,$CHR,$T,$RC\" >> \"$CSV\"\n\n";
    // Step 5: quality
    s += "# Step 5 — quality comparison\n";
    s += "echo \"[5/5] Computing quality metrics...\"\n";
    s += "M=\"$OUTDIR/metrics/\"; mkdir -p \"$M\"\n";
    s += "P=\"$OUTDIR/pmmc/loops_pmmc.hcm\"\n";
    s += "C=\"$OUTDIR/cudammc/loops_cuda.hcm\"\n";
    s += "G=\"$OUTDIR/3dgnome/loops_gnome.hcm\"\n";
    s += "for PAIR in \"pmmc_vs_cuda:$P:$C\" \"pmmc_vs_gnome:$P:$G\"; do\n";
    s += "  NAME=\"${PAIR%%:*}\"; REST=\"${PAIR#*:}\"\n";
    s += "  FA=\"${REST%%:*}\"; FB=\"${REST#*:}\"\n";
    s += "  if [ -f \"$FA\" ] && [ -f \"$FB\" ]; then\n";
    s += "    cat > \"$M/${NAME}.ini\" <<MEOF\n";
    s += "[main]\n";
    s += "action = metrics\n";
    s += "output_dir = $M\n";
    s += "label = $NAME\n";
    s += "chromosomes = $CHR\n";
    s += "level = 2\n";
    s += "resolution = 25000\n";
    s += "input_file_2 = $FB\n";
    s += "[data]\n";
    s += "input_file = $FA\n";
    s += "MEOF\n";
    s += "    \"$PMMC_ANALYZE\" -s \"$M/${NAME}.ini\" || true\n";
    s += "  fi\n";
    s += "done\n\n";
    s += "echo \"\"\necho \"Benchmark complete.  Results: $CSV\"\ncat \"$CSV\"\n";
    writeBenchmarkScript(outdir + "run_benchmark.sh", s);
  }

  // Script 2: run_speed_benchmark.sh — per-chromosome wall-clock timing
  {
    std::string s;
    s += "#!/usr/bin/env bash\n";
    s += "# pMMC Speed Benchmark — per-chromosome wall-clock timing\n";
    s += "# Tools run SEQUENTIALLY to avoid resource competition.\n";
    s += "# Usage:  bash run_speed_benchmark.sh [output_dir] [replicates]\n";
    s += "set -euo pipefail\n\n";
    s += "OUTDIR=\"${1:-./speed_benchmark}\"\nREPS=\"${2:-3}\"\n";
    s += envHeader();
    s += "SEED=42\n";
    s += "CSV=\"$OUTDIR/speed_results.csv\"\n\n";
    s += helpers();
    s += "# Subset for quick run; uncomment the full list for genome-wide\n";
    s += "CHRS=(chr1 chr14 chr21)\n";
    s += "# CHRS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)\n\n";
    s += "mkdir -p \"$OUTDIR\"\n";
    s += "echo \"tool,chromosome,replicate,wall_time_seconds\" > \"$CSV\"\n\n";
    s += "for CHR in \"${CHRS[@]}\"; do\n";
    s += "  SYNTH=\"$OUTDIR/synth_${CHR}/\"\n";
    s += "  if [ ! -f \"$SYNTH/settings.ini\" ]; then\n";
    s += "    mkdir -p \"$SYNTH\"\n";
    s += "    write_generate_ini \"$SYNTH/generate.ini\" \"$SYNTH\" \"$CHR\" 100 1\n";
    s += "    \"$PMMC_GENERATE\" -s \"$SYNTH/generate.ini\"\n";
    s += "  fi\n";
    s += "  for R in $(seq 1 $REPS); do\n";
    s += "    echo \"--- $CHR  rep $R/$REPS ---\"\n";
    s += "    S=$((SEED+R))\n";
    s += "    # 3D-GNOME\n";
    s += "    if command -v \"$GNOME_EXE\" &>/dev/null; then\n";
    s += "      D=\"$OUTDIR/gnome_${CHR}_${R}/\"; mkdir -p \"$D\"\n";
    s += "      write_reconstruct_ini \"$D/settings.ini\" \"$SYNTH\" \"$D\" g \"$CHR\" $S\n";
    s += "      T0=$(date +%s%N)\n";
    s += "      \"$GNOME_EXE\" -a create -s \"$D/settings.ini\" -c \"$CHR\" -o \"$D\" -n g || true\n";
    s += "      T1=$(date +%s%N); echo \"3D-GNOME,$CHR,$R,$(echo \"scale=3;($T1-$T0)/1000000000\"|bc)\" >> \"$CSV\"\n";
    s += "    fi\n";
    s += "    # cudaMMC\n";
    s += "    if command -v \"$CUDAMMC_EXE\" &>/dev/null; then\n";
    s += "      D=\"$OUTDIR/cuda_${CHR}_${R}/\"; mkdir -p \"$D\"\n";
    s += "      write_reconstruct_ini \"$D/settings.ini\" \"$SYNTH\" \"$D\" c \"$CHR\" $S\n";
    s += "      T0=$(date +%s%N)\n";
    s += "      \"$CUDAMMC_EXE\" -a create -s \"$D/settings.ini\" -c \"$CHR\" -o \"$D\" -n c || true\n";
    s += "      T1=$(date +%s%N); echo \"cudaMMC,$CHR,$R,$(echo \"scale=3;($T1-$T0)/1000000000\"|bc)\" >> \"$CSV\"\n";
    s += "    fi\n";
    s += "    # pMMC\n";
    s += "    D=\"$OUTDIR/pmmc_${CHR}_${R}/\"; mkdir -p \"$D\"\n";
    s += "    write_reconstruct_ini \"$D/settings.ini\" \"$SYNTH\" \"$D\" p \"$CHR\" $S\n";
    s += "    T0=$(date +%s%N)\n";
    s += "    \"$PMMC_EXE\" -s \"$D/settings.ini\" || true\n";
    s += "    T1=$(date +%s%N); echo \"pMMC,$CHR,$R,$(echo \"scale=3;($T1-$T0)/1000000000\"|bc)\" >> \"$CSV\"\n";
    s += "  done\n";
    s += "done\n\n";
    s += "echo \"Speed benchmark done.  Results: $CSV\"\ncat \"$CSV\"\n\n";
    // Wall-clock time summary table (mean and CoV)
    s += "# Generate summary table with mean time and CoV per tool/chromosome\n";
    s += "SUMCSV=\"$OUTDIR/wall_clock_summary.csv\"\n";
    s += "echo \"tool,chromosome,mean_time_s,stdev_time_s,cov_percent,num_reps\" > \"$SUMCSV\"\n";
    s += "awk -F, 'NR>1 { key=$1\",\"$2; sum[key]+=$4; ssq[key]+=$4*$4; n[key]++ }\n";
    s += "END { for(k in sum) { m=sum[k]/n[k]; v=ssq[k]/n[k]-m*m; if(v<0)v=0; sd=sqrt(v); cov=(m>0)?sd/m*100:0;\n";
    s += "  printf \"%s,%.3f,%.3f,%.1f,%d\\n\",k,m,sd,cov,n[k] } }' \"$CSV\" | sort >> \"$SUMCSV\"\n";
    s += "echo \"\"\necho \"Wall-clock summary table:\"\ncat \"$SUMCSV\"\n";
    writeBenchmarkScript(outdir + "run_speed_benchmark.sh", s);
  }

  // Script 3: run_scalability.sh — scaling analysis by loop count
  {
    std::string s;
    s += "#!/usr/bin/env bash\n";
    s += "# pMMC Scalability Benchmark — wall time vs loop count\n";
    s += "# Tools run SEQUENTIALLY to avoid resource competition.\n";
    s += "# Usage:  bash run_scalability.sh [output_dir]\n";
    s += "set -euo pipefail\n\n";
    s += "OUTDIR=\"${1:-./scalability_benchmark}\"\n";
    s += envHeader();
    s += "CHR=\"" + chr + "\"\nSEED=42\n";
    s += "CSV=\"$OUTDIR/scalability_results.csv\"\n\n";
    s += helpers();
    s += "LOOPS=(10 25 50 100 250 500 1000)\n\n";
    s += "mkdir -p \"$OUTDIR\"\n";
    s += "echo \"tool,num_loops,wall_time_seconds\" > \"$CSV\"\n\n";
    s += "for NL in \"${LOOPS[@]}\"; do\n";
    s += "  echo \"=== $NL loops ===\"\n";
    s += "  SYNTH=\"$OUTDIR/synth_${NL}/\"; mkdir -p \"$SYNTH\"\n";
    s += "  write_generate_ini \"$SYNTH/generate.ini\" \"$SYNTH\" \"$CHR\" \"$NL\" 1\n";
    s += "  \"$PMMC_GENERATE\" -s \"$SYNTH/generate.ini\"\n\n";
    s += "  if command -v \"$GNOME_EXE\" &>/dev/null; then\n";
    s += "    D=\"$OUTDIR/gnome_${NL}/\"; mkdir -p \"$D\"\n";
    s += "    write_reconstruct_ini \"$D/settings.ini\" \"$SYNTH\" \"$D\" g \"$CHR\" $SEED\n";
    s += "    T0=$(date +%s%N)\n";
    s += "    \"$GNOME_EXE\" -a create -s \"$D/settings.ini\" -c \"$CHR\" -o \"$D\" -n g || true\n";
    s += "    T1=$(date +%s%N); echo \"3D-GNOME,$NL,$(echo \"scale=3;($T1-$T0)/1000000000\"|bc)\" >> \"$CSV\"\n";
    s += "  fi\n";
    s += "  if command -v \"$CUDAMMC_EXE\" &>/dev/null; then\n";
    s += "    D=\"$OUTDIR/cuda_${NL}/\"; mkdir -p \"$D\"\n";
    s += "    write_reconstruct_ini \"$D/settings.ini\" \"$SYNTH\" \"$D\" c \"$CHR\" $SEED\n";
    s += "    T0=$(date +%s%N)\n";
    s += "    \"$CUDAMMC_EXE\" -a create -s \"$D/settings.ini\" -c \"$CHR\" -o \"$D\" -n c || true\n";
    s += "    T1=$(date +%s%N); echo \"cudaMMC,$NL,$(echo \"scale=3;($T1-$T0)/1000000000\"|bc)\" >> \"$CSV\"\n";
    s += "  fi\n";
    s += "  D=\"$OUTDIR/pmmc_${NL}/\"; mkdir -p \"$D\"\n";
    s += "  write_reconstruct_ini \"$D/settings.ini\" \"$SYNTH\" \"$D\" p \"$CHR\" $SEED\n";
    s += "  T0=$(date +%s%N)\n";
    s += "  \"$PMMC_EXE\" -s \"$D/settings.ini\" || true\n";
    s += "  T1=$(date +%s%N); echo \"pMMC,$NL,$(echo \"scale=3;($T1-$T0)/1000000000\"|bc)\" >> \"$CSV\"\n";
    s += "done\n\n";
    s += "echo \"Scalability benchmark done.  Results: $CSV\"\ncat \"$CSV\"\n";
    writeBenchmarkScript(outdir + "run_scalability.sh", s);
  }

  // Script 4: run_ensemble_benchmark.sh — intra-ensemble similarity
  {
    std::string s;
    s += "#!/usr/bin/env bash\n";
    s += "# pMMC Ensemble Benchmark — intra-ensemble structural similarity\n";
    s += "# Tools run SEQUENTIALLY to avoid resource competition.\n";
    s += "# Usage:  bash run_ensemble_benchmark.sh [output_dir] [ensemble_size]\n";
    s += "set -euo pipefail\n\n";
    s += "OUTDIR=\"${1:-./ensemble_benchmark}\"\nENS=\"${2:-50}\"\n";
    s += envHeader(true);
    s += "CHR=\"" + chr + "\"\nSEED=42\n\n";
    s += helpers();
    s += "mkdir -p \"$OUTDIR\"\n\n";
    // Generate data
    s += "# Generate data\n";
    s += "SYNTH=\"$OUTDIR/synthetic/\"; mkdir -p \"$SYNTH\"\n";
    s += "write_generate_ini \"$SYNTH/generate.ini\" \"$SYNTH\" \"$CHR\" 100 \"$ENS\"\n";
    s += "\"$PMMC_GENERATE\" -s \"$SYNTH/generate.ini\"\n\n";
    // pMMC ensemble
    s += "# pMMC ensemble\n";
    s += "echo \"Running pMMC ensemble ($ENS members)...\"\n";
    s += "D=\"$OUTDIR/pmmc_ens/\"; mkdir -p \"$D\"\n";
    s += "write_reconstruct_ini \"$D/settings.ini\" \"$SYNTH\" \"$D\" ens \"$CHR\" $SEED \"$ENS\"\n";
    s += "T0=$(date +%s%N)\n";
    s += "\"$PMMC_EXE\" -s \"$D/settings.ini\"\n";
    s += "T1=$(date +%s%N)\n";
    s += "echo \"pMMC ensemble: $(echo \"scale=1;($T1-$T0)/1000000000\"|bc)s\"\n\n";
    // cudaMMC ensemble
    s += "# cudaMMC ensemble (one run at a time — sequential)\n";
    s += "if command -v \"$CUDAMMC_EXE\" &>/dev/null; then\n";
    s += "  echo \"Running cudaMMC ensemble ($ENS members)...\"\n";
    s += "  D2=\"$OUTDIR/cuda_ens/\"; mkdir -p \"$D2\"\n";
    s += "  T0=$(date +%s%N)\n";
    s += "  for i in $(seq 0 $((ENS-1))); do\n";
    s += "    write_reconstruct_ini \"$D2/settings_${i}.ini\" \"$SYNTH\" \"$D2\" \"ens_${i}\" \"$CHR\" $((SEED+i)) 1\n";
    s += "    \"$CUDAMMC_EXE\" -a create -s \"$D2/settings_${i}.ini\" -c \"$CHR\" -o \"$D2\" -n \"ens_${i}\" || true\n";
    s += "  done\n";
    s += "  T1=$(date +%s%N)\n";
    s += "  echo \"cudaMMC ensemble: $(echo \"scale=1;($T1-$T0)/1000000000\"|bc)s\"\n";
    s += "fi\n\n";
    // Ensemble analysis
    s += "# Compute structural distances\n";
    s += "cat > \"$OUTDIR/ensemble_analysis.ini\" <<EAEOF\n";
    s += "[main]\n";
    s += "action = ensemble\n";
    s += "output_dir = $OUTDIR/pmmc_ens/\n";
    s += "label = ensemble\n";
    s += "pattern = loops_ens_{N}.hcm\n";
    s += "seed = $SEED\n";
    s += "[data]\n";
    s += "input_file = $OUTDIR/pmmc_ens/\n";
    s += "EAEOF\n";
    s += "\"$PMMC_ANALYZE\" -s \"$OUTDIR/ensemble_analysis.ini\" || true\n\n";
    s += "echo \"Ensemble benchmark done.\"\n";
    s += "echo \"Distance matrices: $OUTDIR/pmmc_ens/*distances*.heat\"\n";
    writeBenchmarkScript(outdir + "run_ensemble_benchmark.sh", s);
  }

  // Script 5: run_quality_comparison.sh — Pearson/Spearman correlation
  {
    std::string s;
    s += "#!/usr/bin/env bash\n";
    s += "# pMMC Quality Comparison — run after run_benchmark.sh\n";
    s += "# Usage:  bash run_quality_comparison.sh [benchmark_dir]\n";
    s += "set -euo pipefail\n\n";
    s += "B=\"${1:-./benchmark_results}\"\n";
    s += "PMMC_ANALYZE=\"${PMMC_ANALYZE:-pMMC-analyze}\"\n";
    s += "CHR=\"" + chr + "\"\n";
    s += "Q=\"$B/quality/\"; mkdir -p \"$Q\"\n\n";
    s += "P=\"$B/pmmc/loops_pmmc.hcm\"\n";
    s += "C=\"$B/cudammc/loops_cuda.hcm\"\n";
    s += "G=\"$B/3dgnome/loops_gnome.hcm\"\n\n";
    s += "echo \"pair,scc\" > \"$Q/quality_summary.csv\"\n\n";
    // Distance maps
    s += "for HCM in \"$P\" \"$C\" \"$G\"; do\n";
    s += "  if [ -f \"$HCM\" ]; then\n";
    s += "    NAME=$(basename \"$HCM\" .hcm)\n";
    s += "    cat > \"$Q/distmap_${NAME}.ini\" <<DMEOF\n";
    s += "[main]\n";
    s += "action = distmap\n";
    s += "output_dir = $Q\n";
    s += "label = $NAME\n";
    s += "chromosomes = $CHR\n";
    s += "level = 2\n";
    s += "resolution = 25000\n";
    s += "[data]\n";
    s += "input_file = $HCM\n";
    s += "DMEOF\n";
    s += "    \"$PMMC_ANALYZE\" -s \"$Q/distmap_${NAME}.ini\" || true\n";
    s += "  fi\n";
    s += "done\n\n";
    // Pairwise metrics
    s += "if [ -f \"$P\" ] && [ -f \"$C\" ]; then\n";
    s += "  cat > \"$Q/pmmc_vs_cuda.ini\" <<MCEOF\n";
    s += "[main]\n";
    s += "action = metrics\noutput_dir = $Q\nlabel = pmmc_vs_cuda\nchromosomes = $CHR\nlevel = 2\nresolution = 25000\ninput_file_2 = $C\n";
    s += "[data]\ninput_file = $P\n";
    s += "MCEOF\n";
    s += "  \"$PMMC_ANALYZE\" -s \"$Q/pmmc_vs_cuda.ini\" | tee \"$Q/pmmc_vs_cuda.log\"\n";
    s += "  SCC=$(grep -oP 'SCC.*?\\K[0-9.]+' \"$Q/pmmc_vs_cuda.log\" 2>/dev/null || echo N/A)\n";
    s += "  echo \"pmmc_vs_cuda,$SCC\" >> \"$Q/quality_summary.csv\"\n";
    s += "fi\n";
    s += "if [ -f \"$P\" ] && [ -f \"$G\" ]; then\n";
    s += "  cat > \"$Q/pmmc_vs_gnome.ini\" <<MGEOF\n";
    s += "[main]\n";
    s += "action = metrics\noutput_dir = $Q\nlabel = pmmc_vs_gnome\nchromosomes = $CHR\nlevel = 2\nresolution = 25000\ninput_file_2 = $G\n";
    s += "[data]\ninput_file = $P\n";
    s += "MGEOF\n";
    s += "  \"$PMMC_ANALYZE\" -s \"$Q/pmmc_vs_gnome.ini\" | tee \"$Q/pmmc_vs_gnome.log\"\n";
    s += "  SCC=$(grep -oP 'SCC.*?\\K[0-9.]+' \"$Q/pmmc_vs_gnome.log\" 2>/dev/null || echo N/A)\n";
    s += "  echo \"pmmc_vs_gnome,$SCC\" >> \"$Q/quality_summary.csv\"\n";
    s += "fi\n\n";
    s += "echo \"Quality comparison done.\"\ncat \"$Q/quality_summary.csv\"\n";
    writeBenchmarkScript(outdir + "run_quality_comparison.sh", s);
  }

  // Script 6: run_resource_monitor.sh — resource utilisation capture
  {
    std::string s;
    s += "#!/usr/bin/env bash\n";
    s += "# pMMC Resource Utilisation Monitor\n";
    s += "# Captures GPU util, VRAM, system RAM, and CPU util during a run.\n";
    s += "# Usage:  bash run_resource_monitor.sh [output_dir] [chromosome]\n";
    s += "set -euo pipefail\n\n";
    s += "OUTDIR=\"${1:-./resource_monitor}\"\nCHR=\"${2:-chr1}\"\n";
    s += "PMMC_GENERATE=\"${PMMC_GENERATE:-pMMC-generate}\"\n";
    s += "PMMC_EXE=\"${PMMC_EXE:-pMMC-reconstruct}\"\n";
    s += "SEED=42\n";
    s += "CSV=\"$OUTDIR/resource_utilisation.csv\"\n\n";
    s += helpers();
    s += "mkdir -p \"$OUTDIR\"\n\n";
    // Generate synthetic data
    s += "# Generate synthetic data\n";
    s += "SYNTH=\"$OUTDIR/synth_${CHR}/\"\n";
    s += "if [ ! -f \"$SYNTH/settings.ini\" ]; then\n";
    s += "  mkdir -p \"$SYNTH\"\n";
    s += "  write_generate_ini \"$SYNTH/generate.ini\" \"$SYNTH\" \"$CHR\" 100 1\n";
    s += "  \"$PMMC_GENERATE\" -s \"$SYNTH/generate.ini\"\n";
    s += "fi\n\n";
    // Resource monitor
    s += "# Start resource monitor in background\n";
    s += "MONITOR_LOG=\"$OUTDIR/resource_samples.csv\"\n";
    s += "echo \"timestamp,gpu_util_pct,gpu_mem_used_mb,gpu_mem_total_mb,sys_mem_used_mb,cpu_util_pct\" > \"$MONITOR_LOG\"\n\n";
    s += "monitor_resources() {\n";
    s += "  while true; do\n";
    s += "    TS=$(date +%s.%N)\n";
    s += "    GPU_LINE=$(nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits 2>/dev/null || echo \"0, 0, 0\")\n";
    s += "    GPU_UTIL=$(echo \"$GPU_LINE\" | awk -F', ' '{print $1}')\n";
    s += "    GPU_MEM=$(echo \"$GPU_LINE\" | awk -F', ' '{print $2}')\n";
    s += "    GPU_TOTAL=$(echo \"$GPU_LINE\" | awk -F', ' '{print $3}')\n";
    s += "    if command -v free &>/dev/null; then\n";
    s += "      SYS_MEM=$(free -m | awk '/Mem:/{print $3}')\n";
    s += "    else\n";
    s += "      SYS_MEM=0\n";
    s += "    fi\n";
    s += "    if [ -f /proc/stat ]; then\n";
    s += "      CPU_UTIL=$(awk '/^cpu /{u=$2+$4; t=$2+$4+$5; printf \"%.1f\", u/t*100}' /proc/stat)\n";
    s += "    else\n";
    s += "      CPU_UTIL=0\n";
    s += "    fi\n";
    s += "    echo \"$TS,$GPU_UTIL,$GPU_MEM,$GPU_TOTAL,$SYS_MEM,$CPU_UTIL\" >> \"$MONITOR_LOG\"\n";
    s += "    sleep 1\n";
    s += "  done\n";
    s += "}\n\n";
    s += "monitor_resources &\nMONITOR_PID=$!\n";
    s += "trap \"kill $MONITOR_PID 2>/dev/null || true\" EXIT\n\n";
    // Run pMMC reconstruction
    s += "# Run pMMC reconstruction\n";
    s += "D=\"$OUTDIR/pmmc_${CHR}/\"; mkdir -p \"$D\"\n";
    s += "write_reconstruct_ini \"$D/settings.ini\" \"$SYNTH\" \"$D\" res \"$CHR\" $SEED\n";
    s += "T0=$(date +%s%N)\n";
    s += "\"$PMMC_EXE\" -s \"$D/settings.ini\"\n";
    s += "T1=$(date +%s%N)\nWALL_TIME=$(echo \"scale=3;($T1-$T0)/1000000000\"|bc)\n\n";
    s += "kill $MONITOR_PID 2>/dev/null || true\nsleep 1\n\n";
    // Summary
    s += "# Compute summary\n";
    s += "echo \"metric,value\" > \"$CSV\"\n";
    s += "echo \"chromosome,$CHR\" >> \"$CSV\"\n";
    s += "echo \"wall_time_s,$WALL_TIME\" >> \"$CSV\"\n";
    s += "awk -F, 'NR>1 { if($2+0>mg) mg=$2; if($3+0>mm) mm=$3; gt=$4; if($5+0>sm) sm=$5; if($6+0>mc) mc=$6; n++ }\n";
    s += "END { printf \"gpu_util_max_pct,%.1f\\ngpu_vram_max_mb,%.0f\\ngpu_vram_total_mb,%.0f\\nsys_ram_max_mb,%.0f\\ncpu_util_max_pct,%.1f\\nsamples,%d\\n\",mg,mm,gt,sm,mc,n }' \"$MONITOR_LOG\" >> \"$CSV\"\n\n";
    s += "echo \"Resource utilisation:\"\ncat \"$CSV\"\n";
    writeBenchmarkScript(outdir + "run_resource_monitor.sh", s);
  }

  printf("\n[benchmark] All scripts written to %s\n", outdir.c_str());
  printf("[benchmark] Scripts:\n");
  printf("  run_benchmark.sh          — 3-tool comparison (3D-GNOME, cudaMMC, pMMC)\n");
  printf("  run_speed_benchmark.sh    — per-chromosome wall-clock timing + summary table\n");
  printf("  run_scalability.sh        — performance vs loop count\n");
  printf("  run_ensemble_benchmark.sh — intra-ensemble similarity\n");
  printf("  run_quality_comparison.sh — Pearson/Spearman quality metrics\n");
  printf("  run_resource_monitor.sh   — GPU/CPU/memory utilisation capture\n");
}

int main(int argc, char **argv) {
  setbuf(stdout, NULL);

  if (argc < 3 || std::string(argv[1]) != "-s") {
    usage();
    return 0;
  }

  Settings stg;
  printf("Load settings from [%s]\n", argv[2]);
  if (!stg.loadFromINI(argv[2])) {
    printf("Failed to load settings file!\n");
    return 1;
  }

  // Seed setup
  unsigned int seed;
  if (Settings::seed != 0) {
    seed = (unsigned int)Settings::seed;
    printf("Using fixed random seed: %u\n", seed);
  } else {
    auto tmp = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    seed = (unsigned int)tmp;
    printf("Using time-based random seed: %u\n", seed);
  }
  srand(seed);

  prepareBenchmark();

  printf("end\n");
  return 0;
}
