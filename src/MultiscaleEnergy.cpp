#include <MultiscaleEnergy.h>

#include <cmath>
#include <cstdio>

// --- EnergyTerm ---

EnergyTerm::EnergyTerm() : name(""), value(0.0), previous_value(0.0), weight(1.0) {}

EnergyTerm::EnergyTerm(const std::string &name, double weight)
    : name(name), value(0.0), previous_value(0.0), weight(weight) {}

void EnergyTerm::update(double new_value) {
  previous_value = value;
  value = new_value;
}

void EnergyTerm::revert() { value = previous_value; }

double EnergyTerm::weighted() const { return value * weight; }

double EnergyTerm::delta() const {
  return (value - previous_value) * weight;
}

// --- ScaleEnergy ---

ScaleEnergy::ScaleEnergy()
    : scale_name("unnamed"), moves_proposed(0), moves_accepted(0),
      moves_rejected(0) {}

ScaleEnergy::ScaleEnergy(const std::string &name)
    : scale_name(name), moves_proposed(0), moves_accepted(0),
      moves_rejected(0) {}

void ScaleEnergy::addTerm(const std::string &name, double weight) {
  terms.push_back(EnergyTerm(name, weight));
}

EnergyTerm *ScaleEnergy::getTerm(const std::string &name) {
  for (auto &t : terms) {
    if (t.name == name)
      return &t;
  }
  return nullptr;
}

double ScaleEnergy::totalEnergy() const {
  double sum = 0.0;
  for (const auto &t : terms)
    sum += t.weighted();
  return sum;
}

double ScaleEnergy::totalDelta() const {
  double sum = 0.0;
  for (const auto &t : terms)
    sum += t.delta();
  return sum;
}

void ScaleEnergy::recordProposal() { moves_proposed++; }

void ScaleEnergy::recordAcceptance() { moves_accepted++; }

void ScaleEnergy::recordRejection() { moves_rejected++; }

void ScaleEnergy::revertAll() {
  for (auto &t : terms)
    t.revert();
}

double ScaleEnergy::acceptanceRatio() const {
  if (moves_proposed == 0)
    return 0.0;
  return static_cast<double>(moves_accepted) / moves_proposed;
}

void ScaleEnergy::resetCounters() {
  moves_proposed = 0;
  moves_accepted = 0;
  moves_rejected = 0;
}

// --- MultiscaleEnergy ---

MultiscaleEnergy::MultiscaleEnergy()
    : domain_scale("domain"), loop_scale("loop") {}

void MultiscaleEnergy::initStandardTerms() {
  // Domain scale: heatmap-driven energy terms
  domain_scale.addTerm("heatmap_fit", 1.0);
  domain_scale.addTerm("inter_chr_heatmap", 1.0);
  domain_scale.addTerm("density_fit", 1.0);

  // Loop scale: arc and smoothing energy terms
  loop_scale.addTerm("arc_distance", 1.0);
  loop_scale.addTerm("linker_stretch", 1.0);
  loop_scale.addTerm("linker_squeeze", 1.0);
  loop_scale.addTerm("angular_bending", 1.0);
  loop_scale.addTerm("ctcf_orientation", 1.0);
  loop_scale.addTerm("subanchor_heatmap", 1.0);
}

double MultiscaleEnergy::totalEnergy() const {
  return domain_scale.totalEnergy() + loop_scale.totalEnergy();
}

void MultiscaleEnergy::printDecomposition() const {
  printf("\n========================================\n");
  printf("  MULTISCALE ENERGY DECOMPOSITION\n");
  printf("========================================\n");

  printf("\n  DOMAIN SCALE (%s):\n", domain_scale.scale_name.c_str());
  printf("  %-25s %12s %12s %12s\n", "Term", "Raw", "Weight", "Weighted");
  printf("  %-25s %12s %12s %12s\n", "----", "---", "------", "--------");
  for (const auto &t : domain_scale.terms) {
    if (std::fabs(t.value) > 1e-12)
      printf("  %-25s %12.4f %12.4f %12.4f\n", t.name.c_str(), t.value,
             t.weight, t.weighted());
  }
  printf("  %-25s %12s %12s %12.4f\n", "SUBTOTAL", "", "",
         domain_scale.totalEnergy());
  printf("  Acceptance: %d/%d (%.3f)\n", domain_scale.moves_accepted,
         domain_scale.moves_proposed, domain_scale.acceptanceRatio());

  printf("\n  LOOP SCALE (%s):\n", loop_scale.scale_name.c_str());
  printf("  %-25s %12s %12s %12s\n", "Term", "Raw", "Weight", "Weighted");
  printf("  %-25s %12s %12s %12s\n", "----", "---", "------", "--------");
  for (const auto &t : loop_scale.terms) {
    if (std::fabs(t.value) > 1e-12)
      printf("  %-25s %12.4f %12.4f %12.4f\n", t.name.c_str(), t.value,
             t.weight, t.weighted());
  }
  printf("  %-25s %12s %12s %12.4f\n", "SUBTOTAL", "", "",
         loop_scale.totalEnergy());
  printf("  Acceptance: %d/%d (%.3f)\n", loop_scale.moves_accepted,
         loop_scale.moves_proposed, loop_scale.acceptanceRatio());

  printf("\n  TOTAL ENERGY: %.4f\n", totalEnergy());
  printf("========================================\n\n");
}

void MultiscaleEnergy::writeDecomposition(const std::string &path,
                                          bool append) const {
  FILE *f = fopen(path.c_str(), append ? "a" : "w");
  if (!f)
    return;

  if (!append) {
    fprintf(f, "# Multiscale Energy Decomposition\n");
    fprintf(f, "scale,term,raw_value,weight,weighted_value\n");
  }

  for (const auto &t : domain_scale.terms) {
    fprintf(f, "%s,%s,%.6f,%.6f,%.6f\n", domain_scale.scale_name.c_str(),
            t.name.c_str(), t.value, t.weight, t.weighted());
  }
  for (const auto &t : loop_scale.terms) {
    fprintf(f, "%s,%s,%.6f,%.6f,%.6f\n", loop_scale.scale_name.c_str(),
            t.name.c_str(), t.value, t.weight, t.weighted());
  }

  fclose(f);
}

void MultiscaleEnergy::writeStepLogHeader(FILE *f) {
  fprintf(f, "step,temperature,domain_energy,loop_energy,total_energy,"
             "domain_accept_ratio,loop_accept_ratio\n");
}

void MultiscaleEnergy::writeStepLog(FILE *f, int step,
                                    double temperature) const {
  fprintf(f, "%d,%.6f,%.6f,%.6f,%.6f,%.4f,%.4f\n", step, temperature,
          domain_scale.totalEnergy(), loop_scale.totalEnergy(), totalEnergy(),
          domain_scale.acceptanceRatio(), loop_scale.acceptanceRatio());
}
