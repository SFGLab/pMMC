# pMMC

**Parallel Multiscale Monte Carlo** approach to 3D chromatin structure reconstruction from ChIA-PET data.

pMMC (also known as **3D-GNOME**) reconstructs high-resolution 3D genome organization by combining hierarchical modeling with GPU-accelerated Monte Carlo simulations. It takes ChIA-PET interaction data (anchors, PET clusters, singletons) and produces 3D polymer models of chromatin at multiple scales.

> Based on: Szalaj et al., *"3D-GNOME: 3D Genome Organization Modeling Engine"*, Nucleic Acids Research, 2016.

## Features

- **Multiscale reconstruction** across 4 hierarchical levels: chromosome, segment (~2 Mb), anchor/interaction block, and subanchor (~10 kb)
- **CUDA-accelerated** Monte Carlo simulations with simulated annealing
- **Ensemble modeling** to generate multiple independent structures
- **Synthetic data generation** and automated benchmarking pipeline
- **Structural comparison metrics** including RMSD, Pearson correlation, SCC, and contact decay curves
- **Multiple output formats**: HCM (native binary), mmCIF, and PDB
- **Cross-platform**: Linux, Windows (MSVC), and Docker
- **Deterministic seeding** for reproducible results
- **Memory budget enforcement** with graceful shutdown

## Requirements

- NVIDIA GPU with CUDA Compute Capability >= 6.0
- CUDA Toolkit 11.8+
- CMake 3.13+
- C++17 compiler (GCC, Clang, or MSVC 2017+)

## Building

### Linux

```bash
mkdir build && cd build
cmake .. -DCUDA_ARCH="80" -GNinja
ninja
```

Set `CUDA_ARCH` to match your GPU architecture:

| GPU Generation | Architecture |
|----------------|-------------|
| Pascal (GTX 10xx) | 60 |
| Volta (V100) | 70 |
| Turing (RTX 20xx) | 75 |
| Ampere (RTX 30xx, A100) | 80, 86 |

Multiple architectures can be specified: `-DCUDA_ARCH="70;75;80"`

### Windows

Open the Visual Studio project in `MSVC++/` or use CMake:

```cmd
mkdir build && cd build
cmake .. -DCUDA_ARCH="80"
cmake --build . --config Release
```

### Docker

```bash
docker build -t cudammc:11.8 .
docker run --gpus all -it cudammc:11.8
```

The Docker image builds for architectures 60, 70, 75, 80, and 86.

## Usage

```
cudaMMC -a <action> [options]
```

### Actions

| Action | Alias | Description |
|--------|-------|-------------|
| `create` | `c` | Reconstruct 3D structure from ChIA-PET data (default) |
| `generate` | `g` | Generate synthetic polymer and ChIA-PET data |
| `benchmark` | `k` | Run full benchmark: generate, reconstruct, and compare |
| `distmap` | | Compute distance/contact/frequency maps from structure |
| `metrics` | | Compute structural comparison metrics between two structures |
| `ensemble` | `e` | Pairwise distance analysis across an ensemble of structures |
| `extract` | `r` | Extract a genomic fragment from a structure |
| `position` | `p` | Get 3D coordinates for BED regions |
| `smooth` | `s` | Create equidistant smoothed model |
| `distance` | `d` | Compute structural distance matrix |
| `flatten` | `f` | Flatten hierarchical structure to text coordinates |
| `rewiring` | `w` | Loop rewiring analysis across an ensemble |

### Common Options

| Option | Description |
|--------|-------------|
| `-s FILE` | Settings INI file (required for `create`) |
| `-c REGION` | Chromosomes or region, e.g. `genome`, `chr14`, `chr1-chr5`, `chr14:1:2500000` |
| `-n LABEL` | Label for output file names |
| `-o DIR/` | Output directory (include trailing `/`) |
| `-m N` | Ensemble size (number of structures to generate) |
| `-j SEED` | Fixed random seed for reproducibility |
| `-i FILE` | Input file or directory |
| `-F FORMAT` | Output format: `hcm` (default), `cif`, `pdb`, or `both` |
| `-M MB` | Memory budget in MB (0 = unlimited) |
| `-I METHOD` | Initialization method: `random` (default) or `mds` |
| `-E` | Enable per-step energy trace CSV |
| `-L N` | Energy logging interval in MC steps (default: 1000) |

### Examples

**Reconstruct a single chromosome:**

```bash
cudaMMC -a create -s config.ini -c chr14 -n my_run -o ./output/
```

**Reconstruct entire genome as an ensemble of 10 structures:**

```bash
cudaMMC -a create -s config.ini -c genome -m 10 -o ./ensemble/ -F both
```

**Reconstruct a specific genomic region with a fixed seed:**

```bash
cudaMMC -a create -s config.ini -c chr14:1:2500000 -j 42 -o ./region/
```

**Generate synthetic data and benchmark:**

```bash
cudaMMC -a generate -c chr22 -l 100 -m 20 -o ./synthetic/
cudaMMC -a benchmark -c chr22 -m 20 -o ./benchmark/
```

**Compare two structures:**

```bash
cudaMMC -a metrics -i structure_a.hcm,structure_b.hcm -o ./comparison/
```

**Extract distance and contact maps:**

```bash
cudaMMC -a distmap -i model.hcm -c chr14 -r 25000 -o ./maps/
```

**Ensemble pairwise analysis:**

```bash
cudaMMC -a ensemble -i ./models/ -p "model_{N}.hcm"
```

## Configuration

Reconstruction is driven by an INI settings file (see `config.ini` for a complete example with inline documentation). Key sections:

### `[data]` - Input files

```ini
data_dir = ../data_GM12878/
anchors = GM12878_anchors.bed           # Anchor BED file (e.g. CTCF binding sites)
clusters = GM12878_clusters.bedpe       # PET cluster files (comma-separated)
factors = CTCF                          # Factor names
singletons = GM12878_singletons.bedpe   # Singleton interaction files
segment_split = GM12878_segments.bed    # Segment boundary definitions
centromeres = hg38.bed                  # Centromere regions
```

### `[cuda]` - GPU parameters

```ini
num_threads = 512         # CUDA threads per block
blocks_multiplier = 16    # Grid size multiplier
milestone_fails = 3       # Max consecutive failures before stopping
```

### `[simulation_heatmap]` / `[simulation_arcs]` - Monte Carlo parameters

```ini
max_temp = 5.0                                # Starting temperature
delta_temp = 0.9999                           # Cooling rate
stop_condition_improvement_threshold = 0.975   # Convergence threshold
stop_condition_steps = 50000                   # Steps per milestone
```

### `[distance]` - Interaction-to-distance mapping

```ini
freq_dist_scale = 25.0     # Singleton frequency to 3D distance scale
freq_dist_power = -0.6     # Singleton frequency to 3D distance exponent
count_dist_a = 0.2         # PET count to distance parameter
```

### `[springs]` - Polymer constraints

```ini
stretch_constant = 0.1     # Linker stretch penalty
squeeze_constant = 0.1     # Linker compression penalty
angular_constant = 0.1     # Bending angle penalty
```

See `config.ini` for the full list of ~120 configurable parameters.

## Input Data

pMMC expects ChIA-PET data organized as:

- **Anchors** (`.bed`): Genomic coordinates of protein binding sites (e.g. CTCF peaks)
- **PET clusters** (`.bedpe`): High-confidence chromatin interaction pairs
- **Singletons** (`.bedpe`): Lower-frequency interaction pairs used to build contact heatmaps
- **Segment boundaries** (`.bed`): Defines ~2 Mb segments for the hierarchical decomposition
- **Centromeres** (`.bed`): Centromere positions for chromosome arm separation

## Output Formats

- **HCM** (`.hcm`): Native hierarchical binary format preserving the full 4-level tree structure
- **mmCIF** (`.cif`): Macromolecular Crystallographic Information File for molecular viewers
- **PDB** (`.pdb`): Protein Data Bank format for visualization in tools like PyMOL, Chimera, VMD
- **Text** (`.txt`): Flat 3D coordinate files from `flatten` and `smooth` actions
- **Heatmap** (`.heat`): Distance/contact matrices from `distmap` and `distance` actions
- **CSV**: Energy traces, contact decay curves, and metric reports

## Algorithm Overview

pMMC reconstructs 3D structure through a top-down multiscale approach:

1. **Tree construction**: Build a 4-level hierarchical tree from input data (chromosome > segment > anchor > subanchor)
2. **Heatmap-guided MC** (levels 0-1): Position chromosome territories and segments using singleton contact frequency heatmaps
3. **Arc-distance-guided MC** (level 2): Refine anchor positions using PET cluster interaction distances
4. **Smoothing MC** (level 3): Generate fine-grained subanchor structure with polymer constraints

Each level uses parallel CUDA Monte Carlo with Metropolis acceptance and simulated annealing. The energy function combines:

- Contact frequency matching (heatmap scale)
- Interaction distance matching (arc/loop scale)
- Polymer spring constraints (stretch, squeeze, angular)
- Optional microscopy density constraints
- Optional CTCF motif orientation scoring

## Project Structure

```
pMMC~/
├── src/                    C++ and CUDA source files
│   ├── main.cpp            CLI entry point
│   ├── LooperSolver.cpp    Main reconstruction engine
│   ├── Heatmap.cpp         Contact frequency matrices
│   ├── Chromosome.cpp      Bead-chain polymer model
│   ├── HierarchicalChromosome.cu   Hierarchical tree (CUDA)
│   ├── ParallelMonteCarlo*.cu      MC simulation kernels
│   ├── SyntheticGenerator.cpp      Synthetic data generation
│   ├── BenchmarkRunner.cpp         Benchmark pipeline
│   ├── MetricsFramework.cpp        Structural comparison metrics
│   ├── CifWriter.cpp / PdbWriter.cpp   Output format writers
│   └── ...
├── include/                Header files with Doxygen documentation
├── thirdparty/             Bundled dependencies (INI parser, RMSD, matrix lib)
├── test/                   Test and CI scripts (.bat)
├── MSVC++/                 Visual Studio project files
├── config.ini              Default configuration with inline docs
├── CMakeLists.txt          CMake build configuration
└── Dockerfile              NVIDIA CUDA container build
```

## Testing

Test scripts are provided in the `test/` directory:

```bash
test/build_and_test.bat        # Compile and verify basic functionality
test/ci_local.bat              # Full CI pipeline (build + validate + benchmark)
test/test_determinism.bat      # Reproducibility with fixed seeds
test/test_interchrom.bat       # Inter-chromosomal reconstruction
test/run_benchmarks.bat        # Performance benchmarks
```

## License

See the repository for license information.
