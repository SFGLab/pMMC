#include <DistanceMapGenerator.h>

#include <algorithm>
#include <cmath>
#include <cstdio>

DistanceMapGenerator::DistanceMapGenerator()
    : contact_threshold(2.0f), contact_exponent(2.0f) {}

void DistanceMapGenerator::setContactThreshold(float threshold) {
  contact_threshold = threshold;
}

void DistanceMapGenerator::setContactExponent(float exponent) {
  contact_exponent = exponent;
}

Heatmap DistanceMapGenerator::computeDistanceMap(const Chromosome &chr) const {
  int n = chr.size;
  Heatmap h(n);

  for (int i = 0; i < n; ++i) {
    h.v[i][i] = 0.0f;
    for (int j = i + 1; j < n; ++j) {
      float dist = (chr.points[i] - chr.points[j]).length();
      h.v[i][j] = dist;
      h.v[j][i] = dist;
    }
  }

  return h;
}

Heatmap DistanceMapGenerator::computeContactMap(const Chromosome &chr) const {
  int n = chr.size;
  Heatmap h(n);

  float thresh2 = contact_threshold * contact_threshold;

  for (int i = 0; i < n; ++i) {
    h.v[i][i] = 1.0f;
    for (int j = i + 1; j < n; ++j) {
      vector3 diff = chr.points[i] - chr.points[j];
      float d2 = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
      float val = (d2 < thresh2) ? 1.0f : 0.0f;
      h.v[i][j] = val;
      h.v[j][i] = val;
    }
  }

  return h;
}

Heatmap
DistanceMapGenerator::computeContactFrequencyMap(const Chromosome &chr) const {
  int n = chr.size;
  Heatmap h(n);

  for (int i = 0; i < n; ++i) {
    h.v[i][i] = 0.0f;
    for (int j = i + 1; j < n; ++j) {
      float dist = (chr.points[i] - chr.points[j]).length();
      float freq = 0.0f;
      if (dist > 1e-6f)
        freq = 1.0f / std::pow(dist, contact_exponent);
      h.v[i][j] = freq;
      h.v[j][i] = freq;
    }
  }

  return h;
}

void DistanceMapGenerator::addEnsembleMember(const Chromosome &chr) {
  ensemble.push_back(chr);
}

void DistanceMapGenerator::clearEnsemble() { ensemble.clear(); }

int DistanceMapGenerator::ensembleSize() const {
  return static_cast<int>(ensemble.size());
}

Heatmap DistanceMapGenerator::computeEnsembleAverageDistance() const {
  if (ensemble.empty())
    return Heatmap();

  int n = ensemble[0].size;
  Heatmap avg(n);
  int count = static_cast<int>(ensemble.size());

  for (const auto &chr : ensemble) {
    int m = std::min(n, chr.size);
    for (int i = 0; i < m; ++i) {
      for (int j = i + 1; j < m; ++j) {
        float dist = (chr.points[i] - chr.points[j]).length();
        avg.v[i][j] += dist;
        avg.v[j][i] += dist;
      }
    }
  }

  float inv_count = 1.0f / count;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      avg.v[i][j] *= inv_count;

  return avg;
}

Heatmap DistanceMapGenerator::computeEnsembleAverageContact() const {
  if (ensemble.empty())
    return Heatmap();

  int n = ensemble[0].size;
  Heatmap avg(n);
  int count = static_cast<int>(ensemble.size());
  float thresh2 = contact_threshold * contact_threshold;

  for (const auto &chr : ensemble) {
    int m = std::min(n, chr.size);
    for (int i = 0; i < m; ++i) {
      for (int j = i + 1; j < m; ++j) {
        vector3 diff = chr.points[i] - chr.points[j];
        float d2 = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
        if (d2 < thresh2) {
          avg.v[i][j] += 1.0f;
          avg.v[j][i] += 1.0f;
        }
      }
    }
  }

  float inv_count = 1.0f / count;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      avg.v[i][j] *= inv_count;

  return avg;
}

Heatmap DistanceMapGenerator::computeEnsembleAverageFrequency() const {
  if (ensemble.empty())
    return Heatmap();

  int n = ensemble[0].size;
  Heatmap avg(n);
  int count = static_cast<int>(ensemble.size());

  for (const auto &chr : ensemble) {
    int m = std::min(n, chr.size);
    for (int i = 0; i < m; ++i) {
      for (int j = i + 1; j < m; ++j) {
        float dist = (chr.points[i] - chr.points[j]).length();
        if (dist > 1e-6f) {
          float freq = 1.0f / std::pow(dist, contact_exponent);
          avg.v[i][j] += freq;
          avg.v[j][i] += freq;
        }
      }
    }
  }

  float inv_count = 1.0f / count;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      avg.v[i][j] *= inv_count;

  return avg;
}

void DistanceMapGenerator::writeDistanceMap(
    const Chromosome &chr, const std::string &output_path) const {
  Heatmap h = computeDistanceMap(chr);
  h.toFile(output_path, false);
  printf("Distance map (%dx%d) written to %s\n", (int)h.size, (int)h.size,
         output_path.c_str());
}

void DistanceMapGenerator::writeContactMap(
    const Chromosome &chr, const std::string &output_path) const {
  Heatmap h = computeContactMap(chr);
  h.toFile(output_path, false);
  printf("Contact map (%dx%d) written to %s\n", (int)h.size, (int)h.size,
         output_path.c_str());
}
