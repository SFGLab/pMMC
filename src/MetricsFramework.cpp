#include <MetricsFramework.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>

MetricsFramework::MetricsFramework() {}

float MetricsFramework::pearsonCorrelation(const std::vector<float> &a,
                                           const std::vector<float> &b) const {
  if (a.size() != b.size() || a.empty())
    return 0.0f;

  int n = static_cast<int>(a.size());
  double sum_a = 0.0, sum_b = 0.0;
  for (int i = 0; i < n; ++i) {
    sum_a += a[i];
    sum_b += b[i];
  }
  double mean_a = sum_a / n;
  double mean_b = sum_b / n;

  double cov = 0.0, var_a = 0.0, var_b = 0.0;
  for (int i = 0; i < n; ++i) {
    double da = a[i] - mean_a;
    double db = b[i] - mean_b;
    cov += da * db;
    var_a += da * da;
    var_b += db * db;
  }

  double denom = std::sqrt(var_a * var_b);
  if (denom < 1e-12)
    return 0.0f;

  return static_cast<float>(cov / denom);
}

float MetricsFramework::distanceCorrelation(const Chromosome &a,
                                            const Chromosome &b,
                                            int subsample_step) const {
  int n = std::min(a.size, b.size);
  if (n < 2)
    return 0.0f;

  int step = std::max(1, subsample_step);

  std::vector<float> dist_a, dist_b;
  for (int i = 0; i < n; i += step) {
    for (int j = i + step; j < n; j += step) {
      dist_a.push_back((a.points[i] - a.points[j]).length());
      dist_b.push_back((b.points[i] - b.points[j]).length());
    }
  }

  return pearsonCorrelation(dist_a, dist_b);
}

float MetricsFramework::alignedRMSD(Chromosome &a,
                                     const Chromosome &b) const {
  a.align(b);
  return a.calcRMSD(b);
}

std::vector<ContactDecayPoint>
MetricsFramework::computeContactDecay(const Heatmap &heatmap,
                                      int max_separation,
                                      int bin_size) const {
  int n = static_cast<int>(heatmap.size);
  if (max_separation <= 0 || max_separation > n)
    max_separation = n;
  if (bin_size < 1)
    bin_size = 1;

  int num_bins = (max_separation + bin_size - 1) / bin_size;
  std::vector<double> sum(num_bins, 0.0);
  std::vector<int> count(num_bins, 0);

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      int sep = j - i;
      if (sep > max_separation)
        break;
      int bin = (sep - 1) / bin_size;
      if (bin >= 0 && bin < num_bins) {
        sum[bin] += heatmap.v[i][j];
        count[bin]++;
      }
    }
  }

  std::vector<ContactDecayPoint> result;
  for (int b = 0; b < num_bins; ++b) {
    ContactDecayPoint pt;
    pt.genomic_separation = (b * bin_size) + bin_size / 2 + 1;
    pt.mean_contact = (count[b] > 0) ? static_cast<float>(sum[b] / count[b]) : 0.0f;
    pt.count = count[b];
    result.push_back(pt);
  }

  return result;
}

void MetricsFramework::writeContactDecay(
    const std::vector<ContactDecayPoint> &decay,
    const std::string &output_path) const {
  FILE *f = fopen(output_path.c_str(), "w");
  if (!f) {
    printf("Error: cannot write %s\n", output_path.c_str());
    return;
  }

  fprintf(f, "genomic_separation,mean_contact,pair_count\n");
  for (const auto &pt : decay) {
    fprintf(f, "%d,%.8f,%d\n", pt.genomic_separation, pt.mean_contact,
            pt.count);
  }
  fclose(f);
  printf("Contact decay curve (%d points) written to %s\n",
         static_cast<int>(decay.size()), output_path.c_str());
}

std::vector<LoopEnrichmentResult> MetricsFramework::computeLoopEnrichment(
    const Heatmap &heatmap,
    const std::vector<std::pair<int, int>> &loop_pairs) const {

  int n = static_cast<int>(heatmap.size);

  // Precompute average contact per genomic separation
  std::vector<double> sep_sum(n, 0.0);
  std::vector<int> sep_count(n, 0);

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      int sep = j - i;
      sep_sum[sep] += heatmap.v[i][j];
      sep_count[sep]++;
    }
  }

  std::vector<float> bg_avg(n, 0.0f);
  for (int s = 1; s < n; ++s) {
    if (sep_count[s] > 0)
      bg_avg[s] = static_cast<float>(sep_sum[s] / sep_count[s]);
  }

  std::vector<LoopEnrichmentResult> results;
  for (const auto &lp : loop_pairs) {
    int a = lp.first;
    int b = lp.second;
    if (a < 0 || b < 0 || a >= n || b >= n)
      continue;
    if (a > b)
      std::swap(a, b);

    int sep = b - a;

    LoopEnrichmentResult res;
    res.anchor_a = a;
    res.anchor_b = b;
    res.loop_contact = heatmap.v[a][b];
    res.bg_contact = (sep < n) ? bg_avg[sep] : 0.0f;
    res.enrichment =
        (res.bg_contact > 1e-10f) ? res.loop_contact / res.bg_contact : 0.0f;
    results.push_back(res);
  }

  return results;
}

void MetricsFramework::writeLoopEnrichment(
    const std::vector<LoopEnrichmentResult> &results,
    const std::string &output_path) const {
  FILE *f = fopen(output_path.c_str(), "w");
  if (!f) {
    printf("Error: cannot write %s\n", output_path.c_str());
    return;
  }

  fprintf(f, "anchor_a,anchor_b,loop_contact,bg_contact,enrichment\n");
  for (const auto &r : results) {
    fprintf(f, "%d,%d,%.8f,%.8f,%.4f\n", r.anchor_a, r.anchor_b,
            r.loop_contact, r.bg_contact, r.enrichment);
  }
  fclose(f);
  printf("Loop enrichment (%d loops) written to %s\n",
         static_cast<int>(results.size()), output_path.c_str());
}

float MetricsFramework::structuralSimilarity(const Heatmap &a,
                                             const Heatmap &b,
                                             int h_param) const {
  int n = static_cast<int>(std::min(a.size, b.size));
  if (n < 2)
    return 0.0f;

  auto diag_corrs = perDiagonalCorrelation(a, b);
  if (diag_corrs.empty())
    return 0.0f;

  // Weight by number of pairs at each diagonal and distance kernel
  double weighted_sum = 0.0;
  double weight_total = 0.0;

  for (const auto &dc : diag_corrs) {
    int k = dc.first;
    float corr = dc.second;

    if (std::isnan(corr))
      continue;

    int num_pairs = n - k;
    if (num_pairs < 2)
      continue;

    // Gaussian distance weighting
    double w = num_pairs * std::exp(-static_cast<double>(k * k) /
                                    (2.0 * h_param * h_param));
    weighted_sum += w * corr;
    weight_total += w;
  }

  if (weight_total < 1e-12)
    return 0.0f;

  return static_cast<float>(weighted_sum / weight_total);
}

std::vector<std::pair<int, float>>
MetricsFramework::perDiagonalCorrelation(const Heatmap &a,
                                         const Heatmap &b) const {
  int n = static_cast<int>(std::min(a.size, b.size));
  std::vector<std::pair<int, float>> result;

  for (int k = 1; k < n; ++k) {
    std::vector<float> vals_a, vals_b;
    for (int i = 0; i + k < n; ++i) {
      vals_a.push_back(a.v[i][i + k]);
      vals_b.push_back(b.v[i][i + k]);
    }

    if (vals_a.size() < 2)
      continue;

    float corr = pearsonCorrelation(vals_a, vals_b);
    result.push_back({k, corr});
  }

  return result;
}

void MetricsFramework::writeMetricsReport(const Chromosome &chr_a,
                                          const Chromosome &chr_b,
                                          const Heatmap &heatmap_a,
                                          const Heatmap &heatmap_b,
                                          const std::string &output_path) const {
  FILE *f = fopen(output_path.c_str(), "w");
  if (!f) {
    printf("Error: cannot write %s\n", output_path.c_str());
    return;
  }

  // Distance correlation
  float dist_corr = distanceCorrelation(chr_a, chr_b,
                                        std::max(1, chr_a.size / 500));

  // RMSD (make copies for alignment)
  Chromosome a_copy = chr_a;
  Chromosome b_copy = chr_b;
  a_copy.align(b_copy);
  float rmsd = a_copy.calcRMSD(b_copy);

  // Structural similarity
  float scc = structuralSimilarity(heatmap_a, heatmap_b);

  fprintf(f, "# Structural Comparison Metrics\n");
  fprintf(f, "distance_correlation = %.6f\n", dist_corr);
  fprintf(f, "rmsd = %.6f\n", rmsd);
  fprintf(f, "structural_similarity_scc = %.6f\n", scc);
  fprintf(f, "structure_a_beads = %d\n", chr_a.size);
  fprintf(f, "structure_b_beads = %d\n", chr_b.size);
  fprintf(f, "heatmap_a_size = %d\n", static_cast<int>(heatmap_a.size));
  fprintf(f, "heatmap_b_size = %d\n", static_cast<int>(heatmap_b.size));

  fclose(f);
  printf("Metrics report written to %s\n", output_path.c_str());
}
