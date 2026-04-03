#include <AppCommon.hpp>
#include <Static.hpp>
#include <FlowchartGenerator.hpp>

using namespace std;

void usage() {
  printf("Usage: pMMC-flowchart [options]\n\n");
  printf("Generate algorithm flowchart (SVG).\n\n");
  printf("Options:\n");
  printf("  -o  output file or directory (default: pMMC_algorithm.svg)\n");
  args.clear();
  exit(0);
}

void prepareFlowchart() {
  string output = flag_set('o') ? get_arg('o') : "pMMC_algorithm.svg";
  // If output is a directory, append default filename
  if (!output.empty() && (output.back() == '/' || output.back() == '\\'))
    output += "pMMC_algorithm.svg";
  printf("Generating algorithm flowchart: %s\n", output.c_str());
  if (FlowchartGenerator::generate(output)) {
    printf("Flowchart written to %s\n", output.c_str());
  } else {
    error_msg("Failed to generate flowchart\n");
  }
}

int main(int argc, char **argv) {
  setbuf(stdout, NULL);
  parseArgs(argc, argv);

  if (flag_set('a') && get_arg('a') != "flowchart") {
    printf("This executable only supports the 'flowchart' action.\n");
    usage();
  }

  prepareFlowchart();

  args.clear();
  printf("end\n");
  return 0;
}
