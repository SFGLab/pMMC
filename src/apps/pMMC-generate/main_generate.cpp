#include <AppCommon.hpp>
#include <Static.hpp>
#include <Settings.hpp>
#include <SyntheticGenerator.hpp>

using namespace std;

void usage() {
  printf("pMMC-generate — synthetic data generation\n\n");
  printf("Usage: pMMC-generate -s <settings.ini>\n\n");
  printf("All parameters are controlled by the settings file.\n\n");
  printf("Key settings.ini fields in [main]:\n");
  printf("  output_dir      Output directory (default: ./synthetic/)\n");
  printf("  chromosomes     Chromosome name (default: chr22)\n");
  printf("  label           Label for output files (default: synthetic)\n");
  printf("  num_loops       Number of loops (default: 100)\n");
  printf("  ensemble_size   Ensemble size (default: 20)\n");
  printf("  resolution      Resolution in bp (default: 25000)\n");
  printf("  seed            Random seed (0=time-based)\n");
}

void prepareGenerate() {
  std::string outdir = Settings::outputDir.empty() ? "./synthetic/" : Settings::outputDir;
  std::string chr = Settings::chromosomes.empty() ? "chr22" : Settings::chromosomes;
  std::string name = Settings::label.empty() ? "synthetic" : Settings::label;
  int num_loops = Settings::numLoops;
  int ens = Settings::ensembleSize;
  int res = Settings::resolution;

  // chr22 length for hg19
  int chr_length = 51304566;

  SyntheticGenerator gen;
  gen.setChromosome(chr, chr_length);
  gen.setResolution(res);
  gen.setNumLoops(num_loops);
  gen.setEnsembleSize(ens);
  gen.setOutputDir(outdir);
  gen.setLabel(name);
  gen.generate();
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

  prepareGenerate();

  printf("end\n");
  return 0;
}
