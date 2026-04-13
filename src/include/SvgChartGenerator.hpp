/**
 * @file SvgChartGenerator.hpp
 * @brief Built-in SVG chart generation for pMMC benchmark and diagnostic plots.
 *
 * Generates standalone browser-viewable SVG charts from CSV data, eliminating
 * the need for external plotting tools (gnuplot, matplotlib, etc.).
 *
 * Supported chart types:
 *   - Line chart (single and multi-series)
 *   - Bar chart (simple, grouped, stacked)
 *   - Histogram
 *   - Heatmap (correlation matrix)
 *   - Scatter plot
 */

#ifndef SVG_CHART_GENERATOR_H_
#define SVG_CHART_GENERATOR_H_

#include <string>
#include <vector>

struct ChartSeries {
  std::string label;
  std::vector<float> values;
  std::string color; // SVG color (e.g., "#3b82f6")
};

class SvgChartGenerator {
public:
  // Chart dimensions (defaults)
  static const int DEFAULT_WIDTH = 800;
  static const int DEFAULT_HEIGHT = 500;
  static const int MARGIN_LEFT = 80;
  static const int MARGIN_RIGHT = 30;
  static const int MARGIN_TOP = 50;
  static const int MARGIN_BOTTOM = 70;

  /**
   * Line chart: one or more series plotted against shared x-axis values.
   */
  static bool lineChart(const std::vector<float> &x_values,
                        const std::vector<ChartSeries> &series,
                        const std::string &title,
                        const std::string &x_label,
                        const std::string &y_label,
                        const std::string &output_path,
                        bool log_x = false, bool log_y = false);

  /**
   * Bar chart: one value per category.
   */
  static bool barChart(const std::vector<std::string> &categories,
                       const std::vector<float> &values,
                       const std::string &title,
                       const std::string &x_label,
                       const std::string &y_label,
                       const std::string &output_path,
                       const std::string &color = "#3b82f6");

  /**
   * Grouped bar chart: multiple series per category.
   */
  static bool groupedBarChart(const std::vector<std::string> &categories,
                              const std::vector<ChartSeries> &series,
                              const std::string &title,
                              const std::string &x_label,
                              const std::string &y_label,
                              const std::string &output_path);

  /**
   * Stacked bar chart: components stacked within each category.
   */
  static bool stackedBarChart(const std::vector<std::string> &categories,
                              const std::vector<ChartSeries> &series,
                              const std::string &title,
                              const std::string &x_label,
                              const std::string &y_label,
                              const std::string &output_path);

  /**
   * Histogram: distribution of values in bins.
   */
  static bool histogram(const std::vector<std::string> &bin_labels,
                        const std::vector<float> &counts,
                        const std::string &title,
                        const std::string &x_label,
                        const std::string &y_label,
                        const std::string &output_path,
                        const std::string &color = "#8b5cf6");

  /**
   * Heatmap: 2D matrix with color-coded cells.
   */
  static bool heatmap(const std::vector<std::vector<float>> &matrix,
                      const std::vector<std::string> &row_labels,
                      const std::vector<std::string> &col_labels,
                      const std::string &title,
                      const std::string &output_path,
                      float min_val = 0.0f, float max_val = 1.0f);

  /**
   * Scatter plot: x vs y with optional labels.
   */
  static bool scatterPlot(const std::vector<float> &x_values,
                          const std::vector<float> &y_values,
                          const std::string &title,
                          const std::string &x_label,
                          const std::string &y_label,
                          const std::string &output_path,
                          const std::string &color = "#ef4444");

  // ─── Convenience: read CSV and generate chart ─────────────

  /**
   * Read a CSV and produce an energy convergence line chart.
   * Input: *_energy.csv (step, temperature, domain_energy, loop_energy, total_energy, ...)
   */
  static bool energyConvergenceChart(const std::string &csv_path,
                                     const std::string &svg_path);

  /**
   * Read benchmark_results.csv and produce a multi-panel chart.
   * Input: benchmark_results.csv (loops, dist_correlation, rmsd, elapsed_seconds)
   */
  static bool benchmarkResultsChart(const std::string &csv_path,
                                    const std::string &svg_path);

  /**
   * Read loop_length_distribution CSV and produce histogram.
   */
  static bool loopLengthChart(const std::string &csv_path,
                              const std::string &svg_path);

  /**
   * Read energy_decomposition CSV and produce stacked bar chart.
   */
  static bool energyDecompositionChart(const std::string &csv_path,
                                       const std::string &svg_path);

  /**
   * Read similarity_matrix CSV and produce heatmap.
   */
  static bool similarityHeatmapChart(const std::string &csv_path,
                                      const std::string &svg_path);

  /**
   * Read rmsd_convergence CSV and produce line chart.
   * Input: rmsd_convergence_*.csv (phase,rmsd)
   */
  static bool rmsdConvergenceChart(const std::string &csv_path,
                                    const std::string &svg_path);

  /**
   * Generate a radius of gyration bar chart from territory report data.
   * @param chr_names  Chromosome labels.
   * @param rg_values  Rg values per chromosome.
   * @param svg_path   Output SVG file path.
   */
  static bool rgBarChart(const std::vector<std::string> &chr_names,
                         const std::vector<float> &rg_values,
                         const std::string &svg_path);

private:
  // SVG helpers
  static std::string svgHeader(int width, int height);
  static std::string svgFooter();
  static std::string svgRect(float x, float y, float w, float h,
                             const std::string &fill,
                             const std::string &stroke = "none",
                             float stroke_width = 0, float rx = 0);
  static std::string svgLine(float x1, float y1, float x2, float y2,
                             const std::string &stroke,
                             float stroke_width = 1.0f,
                             const std::string &dash = "");
  static std::string svgText(float x, float y, const std::string &text,
                             int font_size = 12,
                             const std::string &anchor = "middle",
                             const std::string &fill = "#333",
                             float rotate = 0);
  static std::string svgCircle(float cx, float cy, float r,
                               const std::string &fill);
  static std::string svgPath(const std::string &d, const std::string &stroke,
                             float stroke_width = 2.0f,
                             const std::string &fill = "none");

  // Axis/grid helpers
  static std::vector<float> niceTickValues(float min_val, float max_val,
                                           int max_ticks = 8);
  static std::string formatTickLabel(float val);
  static std::string colorScale(float val, float min_val, float max_val);

  // Default palette for multi-series
  static const char *seriesColor(int index);
};


// ============================================================================
// Implementation
// ============================================================================

#include "HeatmapImageWriter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <sstream>

// ─── Color palette ──────────────────────────────────────────────────────────

static const char *PALETTE[] = {
    "#3b82f6", // blue
    "#ef4444", // red
    "#22c55e", // green
    "#f59e0b", // amber
    "#8b5cf6", // purple
    "#06b6d4", // cyan
    "#ec4899", // pink
    "#14b8a6", // teal
    "#f97316", // orange
    "#6366f1", // indigo
};
static const int PALETTE_SIZE = 10;

inline const char *SvgChartGenerator::seriesColor(int index) {
  return PALETTE[index % PALETTE_SIZE];
}

// ─── PNG companion helper ──────────────────────────────────────────────────
// Replace .svg extension with .png
static std::string pngPath(const std::string &svg_path) {
  std::string p = svg_path;
  size_t pos = p.rfind(".svg");
  if (pos != std::string::npos)
    p.replace(pos, 4, ".png");
  else
    p += ".png";
  return p;
}

// Parse hex color "#rrggbb" to Color
static HeatmapImageWriter::Color hexColor(const std::string &hex) {
  HeatmapImageWriter::Color c = {0, 0, 0};
  if (hex.size() >= 7 && hex[0] == '#') {
    unsigned int r, g, b;
    sscanf(hex.c_str() + 1, "%02x%02x%02x", &r, &g, &b);
    c.r = (uint8_t)r;
    c.g = (uint8_t)g;
    c.b = (uint8_t)b;
  }
  return c;
}

// Render a line chart to PNG
static void writeLineChartPNG(const std::string &path,
                              const std::vector<float> &x_values,
                              const std::vector<ChartSeries> &series,
                              const std::string &title,
                              const std::string &x_label,
                              const std::string &y_label) {
  const int W = 800, H = 500;
  const int ml = 80, mr = 150, mt = 50, mb = 70;
  int pw = W - ml - mr, ph = H - mt - mb;

  std::vector<uint8_t> pixels(W * H * 3, 255);
  HeatmapImageWriter::Color black = {0, 0, 0};
  HeatmapImageWriter::Color gray = {200, 200, 200};

  // Title
  HeatmapImageWriter::drawString(pixels, W, H,
      W / 2 - (int)title.size() * 4, 10, title, black, 2);

  // Find range
  float x_min = *std::min_element(x_values.begin(), x_values.end());
  float x_max = *std::max_element(x_values.begin(), x_values.end());
  float y_min = 1e30f, y_max = -1e30f;
  for (const auto &s : series)
    for (float v : s.values) { y_min = std::min(y_min, v); y_max = std::max(y_max, v); }
  float yp = (y_max - y_min) * 0.05f;
  if (yp < 1e-6f) yp = 1.0f;
  y_min -= yp; y_max += yp;
  float xr = x_max - x_min;
  if (xr < 1e-6f) xr = 1.0f;
  float yr = y_max - y_min;

  // Axes
  HeatmapImageWriter::drawLine(pixels, W, H, ml, H - mb, W - mr, H - mb, black);
  HeatmapImageWriter::drawLine(pixels, W, H, ml, mt, ml, H - mb, black);

  // Axis labels
  HeatmapImageWriter::drawString(pixels, W, H,
      ml + pw / 2 - (int)x_label.size() * 3, H - 15, x_label, black, 1);
  HeatmapImageWriter::drawString(pixels, W, H,
      5, mt + ph / 2, y_label, black, 1);

  // Plot series
  for (size_t si = 0; si < series.size(); ++si) {
    const auto &s = series[si];
    HeatmapImageWriter::Color col = hexColor(s.color.empty() ?
        PALETTE[si % PALETTE_SIZE] : s.color.c_str());
    int n = std::min((int)x_values.size(), (int)s.values.size());
    for (int i = 1; i < n; ++i) {
      int x0 = ml + (int)(pw * (x_values[i-1] - x_min) / xr);
      int y0 = mt + (int)(ph * (1.0f - (s.values[i-1] - y_min) / yr));
      int x1 = ml + (int)(pw * (x_values[i] - x_min) / xr);
      int y1 = mt + (int)(ph * (1.0f - (s.values[i] - y_min) / yr));
      HeatmapImageWriter::drawLine(pixels, W, H, x0, y0, x1, y1, col);
      // Thicken line
      HeatmapImageWriter::drawLine(pixels, W, H, x0, y0+1, x1, y1+1, col);
    }
    // Legend
    int ly = mt + 10 + (int)si * 14;
    HeatmapImageWriter::fillRect(pixels, W, H, W - mr + 10, ly, 12, 8, col);
    HeatmapImageWriter::drawString(pixels, W, H, W - mr + 25, ly,
        s.label, black, 1);
  }

  HeatmapImageWriter::writeRawPNG(path, pixels, W, H);
  printf("[PNG] Chart written: %s\n", path.c_str());
}

// Render a bar chart to PNG
static void writeBarChartPNG(const std::string &path,
                             const std::vector<std::string> &categories,
                             const std::vector<float> &values,
                             const std::string &title,
                             const std::string &x_label,
                             const std::string &y_label,
                             const std::string &color) {
  const int W = 800, H = 500;
  const int ml = 80, mr = 30, mt = 50, mb = 70;
  int pw = W - ml - mr, ph = H - mt - mb;
  int n = (int)std::min(categories.size(), values.size());

  std::vector<uint8_t> pixels(W * H * 3, 255);
  HeatmapImageWriter::Color black = {0, 0, 0};
  HeatmapImageWriter::Color barCol = hexColor(color);

  // Title
  HeatmapImageWriter::drawString(pixels, W, H,
      W / 2 - (int)title.size() * 4, 10, title, black, 2);

  float y_max = *std::max_element(values.begin(), values.begin() + n) * 1.1f;
  if (y_max < 1e-6f) y_max = 1.0f;

  // Axes
  HeatmapImageWriter::drawLine(pixels, W, H, ml, H - mb, W - mr, H - mb, black);
  HeatmapImageWriter::drawLine(pixels, W, H, ml, mt, ml, H - mb, black);

  // Bars
  int bar_total_w = pw / n;
  int bar_w = bar_total_w * 7 / 10;
  int bar_gap = bar_total_w * 15 / 100;

  for (int i = 0; i < n; ++i) {
    int bx = ml + i * bar_total_w + bar_gap;
    int bh = (int)(ph * values[i] / y_max);
    int by = mt + ph - bh;
    HeatmapImageWriter::fillRect(pixels, W, H, bx, by, bar_w, bh, barCol);

    // Category label
    std::string lbl = categories[i];
    if (lbl.size() > 10) lbl = lbl.substr(0, 10);
    HeatmapImageWriter::drawString(pixels, W, H,
        bx + bar_w / 2 - (int)lbl.size() * 3, H - mb + 5, lbl, black, 1);
  }

  // Labels
  HeatmapImageWriter::drawString(pixels, W, H,
      ml + pw / 2 - (int)x_label.size() * 3, H - 10, x_label, black, 1);
  HeatmapImageWriter::drawString(pixels, W, H, 5, mt + ph / 2, y_label, black, 1);

  HeatmapImageWriter::writeRawPNG(path, pixels, W, H);
  printf("[PNG] Chart written: %s\n", path.c_str());
}

// Render a grouped bar chart to PNG
static void writeGroupedBarChartPNG(const std::string &path,
                                    const std::vector<std::string> &categories,
                                    const std::vector<ChartSeries> &series,
                                    const std::string &title,
                                    const std::string &x_label,
                                    const std::string &y_label) {
  const int W = 800, H = 500;
  const int ml = 80, mr = 130, mt = 50, mb = 70;
  int pw = W - ml - mr, ph = H - mt - mb;
  int nc = (int)categories.size(), ns = (int)series.size();

  std::vector<uint8_t> pixels(W * H * 3, 255);
  HeatmapImageWriter::Color black = {0, 0, 0};

  HeatmapImageWriter::drawString(pixels, W, H,
      W / 2 - (int)title.size() * 4, 10, title, black, 2);

  float y_max = 0;
  for (const auto &s : series) for (float v : s.values) y_max = std::max(y_max, v);
  y_max *= 1.1f;
  if (y_max < 1e-6f) y_max = 1.0f;

  HeatmapImageWriter::drawLine(pixels, W, H, ml, H - mb, W - mr, H - mb, black);
  HeatmapImageWriter::drawLine(pixels, W, H, ml, mt, ml, H - mb, black);

  int grp_w = pw / nc;
  int bar_w = grp_w * 8 / (10 * ns);

  for (int ci = 0; ci < nc; ++ci) {
    int gx = ml + ci * grp_w + grp_w / 10;
    for (int si = 0; si < ns; ++si) {
      float val = (ci < (int)series[si].values.size()) ? series[si].values[ci] : 0;
      int bh = (int)(ph * val / y_max);
      int bx = gx + si * bar_w;
      HeatmapImageWriter::Color col = hexColor(series[si].color.empty() ?
          PALETTE[si % PALETTE_SIZE] : series[si].color.c_str());
      HeatmapImageWriter::fillRect(pixels, W, H, bx, mt + ph - bh, bar_w - 1, bh, col);
    }
    std::string lbl = categories[ci];
    if (lbl.size() > 8) lbl = lbl.substr(0, 8);
    HeatmapImageWriter::drawString(pixels, W, H,
        ml + ci * grp_w + grp_w / 2 - (int)lbl.size() * 3, H - mb + 5, lbl, black, 1);
  }

  // Legend
  for (int si = 0; si < ns; ++si) {
    HeatmapImageWriter::Color col = hexColor(series[si].color.empty() ?
        PALETTE[si % PALETTE_SIZE] : series[si].color.c_str());
    int ly = mt + 10 + si * 14;
    HeatmapImageWriter::fillRect(pixels, W, H, W - mr + 10, ly, 12, 8, col);
    HeatmapImageWriter::drawString(pixels, W, H, W - mr + 25, ly,
        series[si].label, black, 1);
  }

  HeatmapImageWriter::drawString(pixels, W, H,
      ml + pw / 2 - (int)x_label.size() * 3, H - 10, x_label, black, 1);

  HeatmapImageWriter::writeRawPNG(path, pixels, W, H);
  printf("[PNG] Chart written: %s\n", path.c_str());
}

// Render a heatmap to PNG
static void writeHeatmapPNG(const std::string &path,
                            const std::vector<std::vector<float>> &matrix,
                            const std::vector<std::string> &row_labels,
                            const std::string &title,
                            float min_val, float max_val) {
  int nr = (int)matrix.size();
  if (nr == 0) return;
  int nc = (int)matrix[0].size();

  int cell_size = std::max(4, std::min(40, 500 / std::max(nr, nc)));
  int label_margin = 80;
  int W = label_margin + nc * cell_size + 80;
  int H = 50 + nr * cell_size + 40;

  std::vector<uint8_t> pixels(W * H * 3, 255);
  HeatmapImageWriter::Color black = {0, 0, 0};

  HeatmapImageWriter::drawString(pixels, W, H,
      W / 2 - (int)title.size() * 4, 10, title, black, 2);

  float range = max_val - min_val;
  if (range < 1e-12f) range = 1.0f;

  for (int r = 0; r < nr; ++r) {
    for (int c = 0; c < nc; ++c) {
      float t = (matrix[r][c] - min_val) / range;
      t = std::max(0.0f, std::min(1.0f, t));
      HeatmapImageWriter::Color col = HeatmapImageWriter::valueToColor(t);
      HeatmapImageWriter::fillRect(pixels, W, H,
          label_margin + c * cell_size, 50 + r * cell_size,
          cell_size, cell_size, col);
    }
    // Row label
    if (r < (int)row_labels.size()) {
      std::string rl = row_labels[r];
      if (rl.size() > 10) rl = rl.substr(0, 10);
      HeatmapImageWriter::drawString(pixels, W, H,
          label_margin - (int)rl.size() * 6 - 4,
          50 + r * cell_size + cell_size / 2 - 3, rl, black, 1);
    }
  }

  // Colorbar
  int cbx = label_margin + nc * cell_size + 15;
  int cbh = nr * cell_size;
  for (int y = 0; y < cbh; ++y) {
    float t = 1.0f - (float)y / (float)(cbh - 1);
    HeatmapImageWriter::Color col = HeatmapImageWriter::valueToColor(t);
    HeatmapImageWriter::fillRect(pixels, W, H, cbx, 50 + y, 15, 1, col);
  }

  HeatmapImageWriter::writeRawPNG(path, pixels, W, H);
  printf("[PNG] Chart written: %s\n", path.c_str());
}

// Render a scatter plot to PNG
static void writeScatterPlotPNG(const std::string &path,
                                const std::vector<float> &x_values,
                                const std::vector<float> &y_values,
                                const std::string &title,
                                const std::string &x_label,
                                const std::string &y_label,
                                const std::string &color) {
  const int W = 800, H = 500;
  const int ml = 80, mr = 30, mt = 50, mb = 70;
  int pw = W - ml - mr, ph = H - mt - mb;
  int n = (int)std::min(x_values.size(), y_values.size());

  std::vector<uint8_t> pixels(W * H * 3, 255);
  HeatmapImageWriter::Color black = {0, 0, 0};
  HeatmapImageWriter::Color dotCol = hexColor(color);

  HeatmapImageWriter::drawString(pixels, W, H,
      W / 2 - (int)title.size() * 4, 10, title, black, 2);

  float x_min = *std::min_element(x_values.begin(), x_values.begin() + n);
  float x_max = *std::max_element(x_values.begin(), x_values.begin() + n);
  float y_min = *std::min_element(y_values.begin(), y_values.begin() + n);
  float y_max = *std::max_element(y_values.begin(), y_values.begin() + n);
  float xp = (x_max - x_min) * 0.05f, yp = (y_max - y_min) * 0.05f;
  if (xp < 1e-6f) xp = 1.0f;
  if (yp < 1e-6f) yp = 1.0f;
  x_min -= xp; x_max += xp; y_min -= yp; y_max += yp;
  float xr = x_max - x_min, yr = y_max - y_min;

  HeatmapImageWriter::drawLine(pixels, W, H, ml, H - mb, W - mr, H - mb, black);
  HeatmapImageWriter::drawLine(pixels, W, H, ml, mt, ml, H - mb, black);

  for (int i = 0; i < n; ++i) {
    int px = ml + (int)(pw * (x_values[i] - x_min) / xr);
    int py = mt + (int)(ph * (1.0f - (y_values[i] - y_min) / yr));
    HeatmapImageWriter::fillRect(pixels, W, H, px - 2, py - 2, 5, 5, dotCol);
  }

  HeatmapImageWriter::drawString(pixels, W, H,
      ml + pw / 2 - (int)x_label.size() * 3, H - 10, x_label, black, 1);
  HeatmapImageWriter::drawString(pixels, W, H, 5, mt + ph / 2, y_label, black, 1);

  HeatmapImageWriter::writeRawPNG(path, pixels, W, H);
  printf("[PNG] Chart written: %s\n", path.c_str());
}

// ─── SVG primitives ─────────────────────────────────────────────────────────

inline std::string SvgChartGenerator::svgHeader(int width, int height) {
  char buf[512];
  snprintf(buf, sizeof(buf),
           "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
           "<svg xmlns=\"http://www.w3.org/2000/svg\" "
           "width=\"%d\" height=\"%d\" viewBox=\"0 0 %d %d\">\n"
           "<rect width=\"%d\" height=\"%d\" fill=\"#ffffff\"/>\n",
           width, height, width, height, width, height);
  return std::string(buf);
}

inline std::string SvgChartGenerator::svgFooter() { return "</svg>\n"; }

inline std::string SvgChartGenerator::svgRect(float x, float y, float w, float h,
                                       const std::string &fill,
                                       const std::string &stroke,
                                       float stroke_width, float rx) {
  char buf[512];
  snprintf(buf, sizeof(buf),
           "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" "
           "fill=\"%s\" stroke=\"%s\" stroke-width=\"%.1f\" rx=\"%.1f\"/>\n",
           x, y, w, h, fill.c_str(), stroke.c_str(), stroke_width, rx);
  return std::string(buf);
}

inline std::string SvgChartGenerator::svgLine(float x1, float y1, float x2, float y2,
                                       const std::string &stroke,
                                       float stroke_width,
                                       const std::string &dash) {
  char buf[256];
  if (dash.empty()) {
    snprintf(buf, sizeof(buf),
             "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" "
             "stroke=\"%s\" stroke-width=\"%.1f\"/>\n",
             x1, y1, x2, y2, stroke.c_str(), stroke_width);
  } else {
    snprintf(buf, sizeof(buf),
             "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" "
             "stroke=\"%s\" stroke-width=\"%.1f\" stroke-dasharray=\"%s\"/>\n",
             x1, y1, x2, y2, stroke.c_str(), stroke_width, dash.c_str());
  }
  return std::string(buf);
}

inline std::string SvgChartGenerator::svgText(float x, float y,
                                       const std::string &text, int font_size,
                                       const std::string &anchor,
                                       const std::string &fill, float rotate) {
  char buf[512];
  if (rotate != 0) {
    snprintf(buf, sizeof(buf),
             "<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"%s\" "
             "font-family=\"Arial,sans-serif\" font-size=\"%d\" fill=\"%s\" "
             "transform=\"rotate(%.1f %.1f %.1f)\">%s</text>\n",
             x, y, anchor.c_str(), font_size, fill.c_str(), rotate, x, y,
             text.c_str());
  } else {
    snprintf(buf, sizeof(buf),
             "<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"%s\" "
             "font-family=\"Arial,sans-serif\" font-size=\"%d\" "
             "fill=\"%s\">%s</text>\n",
             x, y, anchor.c_str(), font_size, fill.c_str(), text.c_str());
  }
  return std::string(buf);
}

inline std::string SvgChartGenerator::svgCircle(float cx, float cy, float r,
                                         const std::string &fill) {
  char buf[128];
  snprintf(buf, sizeof(buf),
           "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"%.1f\" fill=\"%s\"/>\n", cx, cy,
           r, fill.c_str());
  return std::string(buf);
}

inline std::string SvgChartGenerator::svgPath(const std::string &d,
                                       const std::string &stroke,
                                       float stroke_width,
                                       const std::string &fill) {
  char buf[64];
  snprintf(buf, sizeof(buf),
           "<path d=\"%s\" stroke=\"%s\" stroke-width=\"%.1f\" fill=\"%s\" "
           "stroke-linejoin=\"round\"/>\n",
           // d is potentially long, build separately
           "", stroke.c_str(), stroke_width, fill.c_str());
  std::string result = "<path d=\"" + d + "\" stroke=\"" + stroke + "\"";
  char attrs[128];
  snprintf(attrs, sizeof(attrs),
           " stroke-width=\"%.1f\" fill=\"%s\" stroke-linejoin=\"round\"/>\n",
           stroke_width, fill.c_str());
  result += attrs;
  return result;
}

// ─── Axis helpers ───────────────────────────────────────────────────────────

inline std::vector<float> SvgChartGenerator::niceTickValues(float min_val,
                                                     float max_val,
                                                     int max_ticks) {
  std::vector<float> ticks;
  if (max_val <= min_val) {
    ticks.push_back(min_val);
    return ticks;
  }

  float range = max_val - min_val;
  float rough_step = range / max_ticks;

  // Round to nice number
  float mag = std::pow(10.0f, std::floor(std::log10(rough_step)));
  float norm = rough_step / mag;
  float nice;
  if (norm <= 1.5f)
    nice = 1.0f;
  else if (norm <= 3.5f)
    nice = 2.0f;
  else if (norm <= 7.5f)
    nice = 5.0f;
  else
    nice = 10.0f;

  float step = nice * mag;
  float start = std::floor(min_val / step) * step;

  for (float v = start; v <= max_val + step * 0.01f; v += step) {
    if (v >= min_val - step * 0.01f)
      ticks.push_back(v);
  }
  return ticks;
}

inline std::string SvgChartGenerator::formatTickLabel(float val) {
  char buf[32];
  if (std::abs(val) >= 1e6f)
    snprintf(buf, sizeof(buf), "%.1fM", val / 1e6f);
  else if (std::abs(val) >= 1e4f)
    snprintf(buf, sizeof(buf), "%.1fk", val / 1e3f);
  else if (val == std::floor(val) && std::abs(val) < 1e6f)
    snprintf(buf, sizeof(buf), "%d", (int)val);
  else
    snprintf(buf, sizeof(buf), "%.2f", val);
  return std::string(buf);
}

inline std::string SvgChartGenerator::colorScale(float val, float min_val,
                                          float max_val) {
  float t = (max_val > min_val) ? (val - min_val) / (max_val - min_val) : 0.5f;
  t = std::max(0.0f, std::min(1.0f, t));

  // Blue (low) -> White (mid) -> Red (high)
  int r, g, b;
  if (t < 0.5f) {
    float s = t * 2.0f;
    r = (int)(59 + (255 - 59) * s);
    g = (int)(130 + (255 - 130) * s);
    b = 246;
  } else {
    float s = (t - 0.5f) * 2.0f;
    r = 239;
    g = (int)(255 - (255 - 68) * s);
    b = (int)(246 - (246 - 68) * s);
  }

  char buf[8];
  snprintf(buf, sizeof(buf), "#%02x%02x%02x", r, g, b);
  return std::string(buf);
}

// ─── Line chart ─────────────────────────────────────────────────────────────

inline bool SvgChartGenerator::lineChart(const std::vector<float> &x_values,
                                  const std::vector<ChartSeries> &series,
                                  const std::string &title,
                                  const std::string &x_label,
                                  const std::string &y_label,
                                  const std::string &output_path, bool log_x,
                                  bool log_y) {
  if (x_values.empty() || series.empty())
    return false;

  int W = DEFAULT_WIDTH, H = DEFAULT_HEIGHT;
  int ml = MARGIN_LEFT, mr = MARGIN_RIGHT + 120; // extra for legend
  int mt = MARGIN_TOP, mb = MARGIN_BOTTOM;
  float pw = (float)(W - ml - mr), ph = (float)(H - mt - mb);

  // Find data range
  float x_min = *std::min_element(x_values.begin(), x_values.end());
  float x_max = *std::max_element(x_values.begin(), x_values.end());
  float y_min = 1e30f, y_max = -1e30f;
  for (const auto &s : series) {
    for (float v : s.values) {
      if (v < y_min) y_min = v;
      if (v > y_max) y_max = v;
    }
  }
  float y_pad = (y_max - y_min) * 0.05f;
  if (y_pad < 1e-6f) y_pad = 1.0f;
  y_min -= y_pad;
  y_max += y_pad;

  std::string svg = svgHeader(W, H);

  // Title
  svg += svgText((float)W / 2, 25, title, 16, "middle", "#1e293b");

  // Grid and axes
  auto y_ticks = niceTickValues(y_min, y_max, 6);
  for (float yt : y_ticks) {
    float py = mt + ph * (1.0f - (yt - y_min) / (y_max - y_min));
    svg += svgLine((float)ml, py, (float)(W - mr), py, "#e5e7eb", 1, "4,4");
    svg += svgText((float)(ml - 8), py + 4, formatTickLabel(yt), 11, "end",
                   "#666");
  }
  auto x_ticks = niceTickValues(x_min, x_max, 8);
  for (float xt : x_ticks) {
    float px = ml + pw * (xt - x_min) / (x_max - x_min);
    svg += svgLine(px, (float)mt, px, (float)(H - mb), "#e5e7eb", 1, "4,4");
    svg += svgText(px, (float)(H - mb + 18), formatTickLabel(xt), 11, "middle",
                   "#666");
  }

  // Axis lines
  svg += svgLine((float)ml, (float)(H - mb), (float)(W - mr), (float)(H - mb),
                 "#475569", 1.5f);
  svg += svgLine((float)ml, (float)mt, (float)ml, (float)(H - mb), "#475569",
                 1.5f);

  // Labels
  svg += svgText((float)(ml + pw / 2), (float)(H - 15), x_label, 13, "middle",
                 "#333");
  svg += svgText(18, (float)(mt + ph / 2), y_label, 13, "middle", "#333",
                 -90);

  // Plot lines
  for (size_t si = 0; si < series.size(); ++si) {
    const auto &s = series[si];
    std::string color = s.color.empty() ? seriesColor((int)si) : s.color;
    int n = std::min((int)x_values.size(), (int)s.values.size());
    if (n < 2) continue;

    std::string path_d;
    for (int i = 0; i < n; ++i) {
      float px = ml + pw * (x_values[i] - x_min) / (x_max - x_min);
      float py =
          mt + ph * (1.0f - (s.values[i] - y_min) / (y_max - y_min));
      char pt[32];
      snprintf(pt, sizeof(pt), "%s%.1f %.1f", i == 0 ? "M" : " L", px, py);
      path_d += pt;
    }
    svg += svgPath(path_d, color, 2.0f);

    // Legend entry
    float ly = (float)(mt + 15 + si * 20);
    float lx = (float)(W - mr + 10);
    svg += svgLine(lx, ly, lx + 20, ly, color, 2.5f);
    svg += svgText(lx + 25, ly + 4, s.label, 11, "start", "#333");
  }

  svg += svgFooter();

  FILE *f = fopen(output_path.c_str(), "w");
  if (!f) return false;
  fprintf(f, "%s", svg.c_str());
  fclose(f);
  printf("[SVG] Chart written: %s\n", output_path.c_str());

  // PNG companion
  writeLineChartPNG(pngPath(output_path), x_values, series, title,
                    x_label, y_label);
  return true;
}

// ─── Bar chart ──────────────────────────────────────────────────────────────

inline bool SvgChartGenerator::barChart(const std::vector<std::string> &categories,
                                 const std::vector<float> &values,
                                 const std::string &title,
                                 const std::string &x_label,
                                 const std::string &y_label,
                                 const std::string &output_path,
                                 const std::string &color) {
  if (categories.empty() || values.empty()) return false;
  int n = (int)std::min(categories.size(), values.size());

  int W = DEFAULT_WIDTH, H = DEFAULT_HEIGHT;
  int ml = MARGIN_LEFT, mr = MARGIN_RIGHT, mt = MARGIN_TOP, mb = MARGIN_BOTTOM;
  float pw = (float)(W - ml - mr), ph = (float)(H - mt - mb);

  float y_max = *std::max_element(values.begin(), values.begin() + n);
  float y_min = 0.0f;
  y_max *= 1.1f;
  if (y_max < 1e-6f) y_max = 1.0f;

  std::string svg = svgHeader(W, H);
  svg += svgText((float)W / 2, 25, title, 16, "middle", "#1e293b");

  // Y-axis grid
  auto y_ticks = niceTickValues(y_min, y_max, 6);
  for (float yt : y_ticks) {
    float py = mt + ph * (1.0f - (yt - y_min) / (y_max - y_min));
    svg += svgLine((float)ml, py, (float)(W - mr), py, "#e5e7eb", 1, "4,4");
    svg += svgText((float)(ml - 8), py + 4, formatTickLabel(yt), 11, "end",
                   "#666");
  }

  // Bars
  float bar_total_w = pw / n;
  float bar_w = bar_total_w * 0.7f;
  float bar_gap = bar_total_w * 0.15f;

  for (int i = 0; i < n; ++i) {
    float bx = ml + i * bar_total_w + bar_gap;
    float bh = ph * values[i] / y_max;
    float by = mt + ph - bh;
    svg += svgRect(bx, by, bar_w, bh, color, "none", 0, 3);

    // Category label
    float tx = ml + i * bar_total_w + bar_total_w / 2;
    std::string label = categories[i];
    if (label.length() > 10)
      svg += svgText(tx, (float)(H - mb + 18), label, 10, "end", "#666", -35);
    else
      svg += svgText(tx, (float)(H - mb + 18), label, 11, "middle", "#666");

    // Value label on bar
    char vbuf[32];
    snprintf(vbuf, sizeof(vbuf), "%.1f", values[i]);
    svg += svgText(bx + bar_w / 2, by - 5, vbuf, 10, "middle", "#333");
  }

  // Axes
  svg += svgLine((float)ml, (float)(H - mb), (float)(W - mr), (float)(H - mb),
                 "#475569", 1.5f);
  svg += svgLine((float)ml, (float)mt, (float)ml, (float)(H - mb), "#475569",
                 1.5f);
  svg += svgText((float)(ml + pw / 2), (float)(H - 10), x_label, 13, "middle",
                 "#333");
  svg += svgText(18, (float)(mt + ph / 2), y_label, 13, "middle", "#333",
                 -90);

  svg += svgFooter();

  FILE *f = fopen(output_path.c_str(), "w");
  if (!f) return false;
  fprintf(f, "%s", svg.c_str());
  fclose(f);
  printf("[SVG] Chart written: %s\n", output_path.c_str());

  // PNG companion
  writeBarChartPNG(pngPath(output_path), categories, values, title,
                   x_label, y_label, color);
  return true;
}

// ─── Grouped bar chart ──────────────────────────────────────────────────────

inline bool SvgChartGenerator::groupedBarChart(
    const std::vector<std::string> &categories,
    const std::vector<ChartSeries> &series, const std::string &title,
    const std::string &x_label, const std::string &y_label,
    const std::string &output_path) {
  if (categories.empty() || series.empty()) return false;

  int nc = (int)categories.size(), ns = (int)series.size();
  int W = DEFAULT_WIDTH, H = DEFAULT_HEIGHT;
  int ml = MARGIN_LEFT, mr = MARGIN_RIGHT + 100, mt = MARGIN_TOP,
      mb = MARGIN_BOTTOM;
  float pw = (float)(W - ml - mr), ph = (float)(H - mt - mb);

  float y_max = 0;
  for (const auto &s : series)
    for (float v : s.values)
      y_max = std::max(y_max, v);
  y_max *= 1.1f;
  if (y_max < 1e-6f) y_max = 1.0f;

  std::string svg = svgHeader(W, H);
  svg += svgText((float)W / 2, 25, title, 16, "middle", "#1e293b");

  auto y_ticks = niceTickValues(0, y_max, 6);
  for (float yt : y_ticks) {
    float py = mt + ph * (1.0f - yt / y_max);
    svg += svgLine((float)ml, py, (float)(W - mr), py, "#e5e7eb", 1, "4,4");
    svg += svgText((float)(ml - 8), py + 4, formatTickLabel(yt), 11, "end",
                   "#666");
  }

  float grp_w = pw / nc;
  float bar_w = grp_w * 0.8f / ns;

  for (int ci = 0; ci < nc; ++ci) {
    float gx = ml + ci * grp_w + grp_w * 0.1f;
    for (int si = 0; si < ns; ++si) {
      float val = (ci < (int)series[si].values.size()) ? series[si].values[ci] : 0;
      float bh = ph * val / y_max;
      float bx = gx + si * bar_w;
      std::string c = series[si].color.empty() ? seriesColor(si) : series[si].color;
      svg += svgRect(bx, mt + ph - bh, bar_w - 1, bh, c, "none", 0, 2);
    }
    svg += svgText(ml + ci * grp_w + grp_w / 2, (float)(H - mb + 18),
                   categories[ci], 11, "middle", "#666");
  }

  // Legend
  for (int si = 0; si < ns; ++si) {
    float ly = (float)(mt + 15 + si * 20);
    float lx = (float)(W - mr + 10);
    std::string c = series[si].color.empty() ? seriesColor(si) : series[si].color;
    svg += svgRect(lx, ly - 6, 14, 14, c, "none", 0, 2);
    svg += svgText(lx + 18, ly + 4, series[si].label, 11, "start", "#333");
  }

  svg += svgLine((float)ml, (float)(H - mb), (float)(W - mr), (float)(H - mb),
                 "#475569", 1.5f);
  svg += svgLine((float)ml, (float)mt, (float)ml, (float)(H - mb), "#475569",
                 1.5f);
  svg += svgText((float)(ml + pw / 2), (float)(H - 10), x_label, 13, "middle",
                 "#333");
  svg += svgText(18, (float)(mt + ph / 2), y_label, 13, "middle", "#333",
                 -90);

  svg += svgFooter();
  FILE *f = fopen(output_path.c_str(), "w");
  if (!f) return false;
  fprintf(f, "%s", svg.c_str());
  fclose(f);
  printf("[SVG] Chart written: %s\n", output_path.c_str());

  // PNG companion
  writeGroupedBarChartPNG(pngPath(output_path), categories, series, title,
                          x_label, y_label);
  return true;
}

// ─── Stacked bar chart ──────────────────────────────────────────────────────

inline bool SvgChartGenerator::stackedBarChart(
    const std::vector<std::string> &categories,
    const std::vector<ChartSeries> &series, const std::string &title,
    const std::string &x_label, const std::string &y_label,
    const std::string &output_path) {
  if (categories.empty() || series.empty()) return false;

  int nc = (int)categories.size(), ns = (int)series.size();
  int W = DEFAULT_WIDTH, H = DEFAULT_HEIGHT;
  int ml = MARGIN_LEFT, mr = MARGIN_RIGHT + 120, mt = MARGIN_TOP,
      mb = MARGIN_BOTTOM;
  float pw = (float)(W - ml - mr), ph = (float)(H - mt - mb);

  // Compute stacked totals
  float y_max = 0;
  for (int ci = 0; ci < nc; ++ci) {
    float total = 0;
    for (const auto &s : series)
      total += (ci < (int)s.values.size()) ? s.values[ci] : 0;
    y_max = std::max(y_max, total);
  }
  y_max *= 1.1f;
  if (y_max < 1e-6f) y_max = 1.0f;

  std::string svg = svgHeader(W, H);
  svg += svgText((float)W / 2, 25, title, 16, "middle", "#1e293b");

  auto y_ticks = niceTickValues(0, y_max, 6);
  for (float yt : y_ticks) {
    float py = mt + ph * (1.0f - yt / y_max);
    svg += svgLine((float)ml, py, (float)(W - mr), py, "#e5e7eb", 1, "4,4");
    svg += svgText((float)(ml - 8), py + 4, formatTickLabel(yt), 11, "end",
                   "#666");
  }

  float bar_total_w = pw / nc;
  float bar_w = bar_total_w * 0.7f;

  for (int ci = 0; ci < nc; ++ci) {
    float bx = ml + ci * bar_total_w + bar_total_w * 0.15f;
    float cum = 0;
    for (int si = 0; si < ns; ++si) {
      float val = (ci < (int)series[si].values.size()) ? series[si].values[ci] : 0;
      float bh = ph * val / y_max;
      float by = mt + ph - cum * ph / y_max - bh;
      std::string c = series[si].color.empty() ? seriesColor(si) : series[si].color;
      svg += svgRect(bx, by, bar_w, bh, c, "#fff", 0.5f, 2);
      cum += val;
    }
    svg += svgText(ml + ci * bar_total_w + bar_total_w / 2,
                   (float)(H - mb + 18), categories[ci], 11, "middle", "#666");
  }

  for (int si = 0; si < ns; ++si) {
    float ly = (float)(mt + 15 + si * 20);
    float lx = (float)(W - mr + 10);
    std::string c = series[si].color.empty() ? seriesColor(si) : series[si].color;
    svg += svgRect(lx, ly - 6, 14, 14, c, "none", 0, 2);
    svg += svgText(lx + 18, ly + 4, series[si].label, 11, "start", "#333");
  }

  svg += svgLine((float)ml, (float)(H - mb), (float)(W - mr), (float)(H - mb),
                 "#475569", 1.5f);
  svg += svgLine((float)ml, (float)mt, (float)ml, (float)(H - mb), "#475569",
                 1.5f);
  svg += svgText((float)(ml + pw / 2), (float)(H - 10), x_label, 13, "middle",
                 "#333");
  svg += svgText(18, (float)(mt + ph / 2), y_label, 13, "middle", "#333",
                 -90);

  svg += svgFooter();
  FILE *f = fopen(output_path.c_str(), "w");
  if (!f) return false;
  fprintf(f, "%s", svg.c_str());
  fclose(f);
  printf("[SVG] Chart written: %s\n", output_path.c_str());

  // PNG companion (render as grouped bar for stacked)
  writeGroupedBarChartPNG(pngPath(output_path), categories, series, title,
                          x_label, y_label);
  return true;
}

// ─── Histogram ──────────────────────────────────────────────────────────────

inline bool SvgChartGenerator::histogram(const std::vector<std::string> &bin_labels,
                                  const std::vector<float> &counts,
                                  const std::string &title,
                                  const std::string &x_label,
                                  const std::string &y_label,
                                  const std::string &output_path,
                                  const std::string &color) {
  return barChart(bin_labels, counts, title, x_label, y_label, output_path,
                  color);
}

// ─── Heatmap ────────────────────────────────────────────────────────────────

inline bool SvgChartGenerator::heatmap(const std::vector<std::vector<float>> &matrix,
                                const std::vector<std::string> &row_labels,
                                const std::vector<std::string> &col_labels,
                                const std::string &title,
                                const std::string &output_path, float min_val,
                                float max_val) {
  int nr = (int)matrix.size();
  if (nr == 0) return false;
  int nc = (int)matrix[0].size();

  int label_margin = 100;
  int cell_size = std::max(15, std::min(40, 500 / std::max(nr, nc)));
  int W = label_margin + nc * cell_size + 80; // +80 for color bar
  int H = MARGIN_TOP + nr * cell_size + 40;

  std::string svg = svgHeader(W, H);
  svg += svgText((float)W / 2, 25, title, 16, "middle", "#1e293b");

  float ox = (float)label_margin, oy = (float)MARGIN_TOP;

  for (int r = 0; r < nr; ++r) {
    for (int c = 0; c < nc; ++c) {
      float val = matrix[r][c];
      std::string clr = colorScale(val, min_val, max_val);
      svg += svgRect(ox + c * cell_size, oy + r * cell_size, (float)cell_size,
                     (float)cell_size, clr, "#fff", 0.5f);

      // Value text for small matrices
      if (nr <= 15 && nc <= 15 && cell_size >= 25) {
        char vb[16];
        snprintf(vb, sizeof(vb), "%.2f", val);
        svg += svgText(ox + c * cell_size + cell_size / 2.0f,
                       oy + r * cell_size + cell_size / 2.0f + 4, vb, 9,
                       "middle", "#333");
      }
    }
    // Row label
    if (r < (int)row_labels.size()) {
      std::string rl = row_labels[r];
      if (rl.length() > 12) rl = rl.substr(0, 12);
      svg += svgText(ox - 5, oy + r * cell_size + cell_size / 2.0f + 4, rl, 10,
                     "end", "#333");
    }
  }

  // Column labels
  for (int c = 0; c < nc && c < (int)col_labels.size(); ++c) {
    std::string cl = col_labels[c];
    if (cl.length() > 12) cl = cl.substr(0, 12);
    svg += svgText(ox + c * cell_size + cell_size / 2.0f,
                   oy + nr * cell_size + 15, cl, 10, "end", "#333", -45);
  }

  // Color bar
  float cbx = ox + nc * cell_size + 20;
  float cby = oy;
  float cbh = (float)(nr * cell_size);
  int cb_steps = 20;
  for (int i = 0; i < cb_steps; ++i) {
    float t = 1.0f - (float)i / cb_steps;
    float val = min_val + t * (max_val - min_val);
    std::string clr = colorScale(val, min_val, max_val);
    svg += svgRect(cbx, cby + i * cbh / cb_steps, 15,
                   cbh / cb_steps + 1, clr);
  }
  svg += svgText(cbx + 20, cby + 4, formatTickLabel(max_val), 10, "start",
                 "#666");
  svg += svgText(cbx + 20, cby + cbh + 4, formatTickLabel(min_val), 10,
                 "start", "#666");

  svg += svgFooter();
  FILE *f = fopen(output_path.c_str(), "w");
  if (!f) return false;
  fprintf(f, "%s", svg.c_str());
  fclose(f);
  printf("[SVG] Chart written: %s\n", output_path.c_str());

  // PNG companion
  writeHeatmapPNG(pngPath(output_path), matrix, row_labels, title,
                  min_val, max_val);
  return true;
}

// ─── Scatter plot ───────────────────────────────────────────────────────────

inline bool SvgChartGenerator::scatterPlot(const std::vector<float> &x_values,
                                    const std::vector<float> &y_values,
                                    const std::string &title,
                                    const std::string &x_label,
                                    const std::string &y_label,
                                    const std::string &output_path,
                                    const std::string &color) {
  int n = (int)std::min(x_values.size(), y_values.size());
  if (n == 0) return false;

  int W = DEFAULT_WIDTH, H = DEFAULT_HEIGHT;
  int ml = MARGIN_LEFT, mr = MARGIN_RIGHT, mt = MARGIN_TOP, mb = MARGIN_BOTTOM;
  float pw = (float)(W - ml - mr), ph = (float)(H - mt - mb);

  float x_min = *std::min_element(x_values.begin(), x_values.begin() + n);
  float x_max = *std::max_element(x_values.begin(), x_values.begin() + n);
  float y_min = *std::min_element(y_values.begin(), y_values.begin() + n);
  float y_max = *std::max_element(y_values.begin(), y_values.begin() + n);
  float xp = (x_max - x_min) * 0.05f, yp = (y_max - y_min) * 0.05f;
  if (xp < 1e-6f) xp = 1.0f;
  if (yp < 1e-6f) yp = 1.0f;
  x_min -= xp; x_max += xp; y_min -= yp; y_max += yp;

  std::string svg = svgHeader(W, H);
  svg += svgText((float)W / 2, 25, title, 16, "middle", "#1e293b");

  auto y_ticks = niceTickValues(y_min, y_max, 6);
  for (float yt : y_ticks) {
    float py = mt + ph * (1.0f - (yt - y_min) / (y_max - y_min));
    svg += svgLine((float)ml, py, (float)(W - mr), py, "#e5e7eb", 1, "4,4");
    svg += svgText((float)(ml - 8), py + 4, formatTickLabel(yt), 11, "end",
                   "#666");
  }
  auto x_ticks = niceTickValues(x_min, x_max, 8);
  for (float xt : x_ticks) {
    float px = ml + pw * (xt - x_min) / (x_max - x_min);
    svg += svgText(px, (float)(H - mb + 18), formatTickLabel(xt), 11, "middle",
                   "#666");
  }

  svg += svgLine((float)ml, (float)(H - mb), (float)(W - mr), (float)(H - mb),
                 "#475569", 1.5f);
  svg += svgLine((float)ml, (float)mt, (float)ml, (float)(H - mb), "#475569",
                 1.5f);
  svg += svgText((float)(ml + pw / 2), (float)(H - 15), x_label, 13, "middle",
                 "#333");
  svg += svgText(18, (float)(mt + ph / 2), y_label, 13, "middle", "#333",
                 -90);

  for (int i = 0; i < n; ++i) {
    float px = ml + pw * (x_values[i] - x_min) / (x_max - x_min);
    float py = mt + ph * (1.0f - (y_values[i] - y_min) / (y_max - y_min));
    svg += svgCircle(px, py, 4, color);
  }

  svg += svgFooter();
  FILE *f = fopen(output_path.c_str(), "w");
  if (!f) return false;
  fprintf(f, "%s", svg.c_str());
  fclose(f);
  printf("[SVG] Chart written: %s\n", output_path.c_str());

  // PNG companion
  writeScatterPlotPNG(pngPath(output_path), x_values, y_values, title,
                      x_label, y_label, color);
  return true;
}

// ─── Convenience: CSV -> Chart ──────────────────────────────────────────────

inline bool SvgChartGenerator::energyConvergenceChart(const std::string &csv_path,
                                               const std::string &svg_path) {
  FILE *f = fopen(csv_path.c_str(), "r");
  if (!f) return false;

  std::vector<float> steps, domain_e, loop_e, total_e;
  char line[512];
  fgets(line, sizeof(line), f); // skip header

  int step;
  float temp, de, le, te, da, la;
  while (fscanf(f, "%d,%f,%f,%f,%f,%f,%f", &step, &temp, &de, &le, &te, &da,
                &la) == 7) {
    steps.push_back((float)step);
    domain_e.push_back(de);
    loop_e.push_back(le);
    total_e.push_back(te);
  }
  fclose(f);

  if (steps.empty()) return false;

  std::vector<ChartSeries> series;
  series.push_back({"Total Energy", total_e, "#1e293b"});
  series.push_back({"Domain Energy", domain_e, "#3b82f6"});
  series.push_back({"Loop Energy", loop_e, "#ef4444"});

  return lineChart(steps, series, "Energy Convergence", "MC Step", "Energy",
                   svg_path);
}

inline bool SvgChartGenerator::benchmarkResultsChart(const std::string &csv_path,
                                              const std::string &svg_path) {
  FILE *f = fopen(csv_path.c_str(), "r");
  if (!f) return false;

  std::vector<std::string> cats;
  std::vector<float> corrs, rmsds, times;
  char line[512];
  fgets(line, sizeof(line), f); // skip header

  int loops;
  float dc, rm, el;
  while (fscanf(f, "%d,%f,%f,%f", &loops, &dc, &rm, &el) == 4) {
    char buf[32];
    snprintf(buf, sizeof(buf), "%d", loops);
    cats.push_back(buf);
    corrs.push_back(dc);
    rmsds.push_back(rm);
    times.push_back(el);
  }
  fclose(f);

  if (cats.empty()) return false;

  std::vector<ChartSeries> series;
  series.push_back({"Dist. Corr.", corrs, "#3b82f6"});
  series.push_back({"RMSD", rmsds, "#ef4444"});
  series.push_back({"Time (s)", times, "#22c55e"});

  return groupedBarChart(cats, series, "Benchmark Results",
                         "Number of Loops", "Value", svg_path);
}

inline bool SvgChartGenerator::loopLengthChart(const std::string &csv_path,
                                        const std::string &svg_path) {
  FILE *f = fopen(csv_path.c_str(), "r");
  if (!f) return false;

  std::vector<std::string> labels;
  std::vector<float> counts;
  char line[512];
  fgets(line, sizeof(line), f); // skip header

  long long bs, be;
  int cnt;
  float frac;
  while (fscanf(f, "%lld,%lld,%d,%f", &bs, &be, &cnt, &frac) == 4) {
    if (cnt == 0) continue;
    char buf[64];
    if (be < 0)
      snprintf(buf, sizeof(buf), ">%lldMb", bs / 1000000);
    else if (be >= 1000000)
      snprintf(buf, sizeof(buf), "%lld-%lldMb", bs / 1000000, be / 1000000);
    else
      snprintf(buf, sizeof(buf), "%lld-%lldk", bs / 1000, be / 1000);
    labels.push_back(buf);
    counts.push_back((float)cnt);
  }
  fclose(f);

  if (labels.empty()) return false;

  return histogram(labels, counts, "Loop Length Distribution",
                   "Loop Length", "Count", svg_path, "#8b5cf6");
}

inline bool SvgChartGenerator::energyDecompositionChart(const std::string &csv_path,
                                                  const std::string &svg_path) {
  FILE *f = fopen(csv_path.c_str(), "r");
  if (!f) return false;

  std::vector<std::string> terms;
  std::vector<float> raw_vals, weighted_vals;
  char line[512];

  // Skip comment and header lines (lines starting with # or containing "scale,term")
  while (fgets(line, sizeof(line), f)) {
    if (line[0] != '#' && strstr(line, "scale") == nullptr)
      break;
  }

  // Now parse data lines - first line was already read
  char scale[64], term[64];
  float raw, weight, weighted;
  do {
    if (sscanf(line, "%[^,],%[^,],%f,%f,%f", scale, term, &raw, &weight,
               &weighted) == 5) {
      if (raw > 0.001f || weighted > 0.001f) { // skip zero-valued terms
        std::string lbl = std::string(scale) + ":" + term;
        terms.push_back(lbl);
        raw_vals.push_back(raw);
        weighted_vals.push_back(weighted);
      }
    }
  } while (fgets(line, sizeof(line), f));
  fclose(f);

  if (terms.empty()) return false;

  std::vector<ChartSeries> series;
  series.push_back({"Raw", raw_vals, "#3b82f6"});
  series.push_back({"Weighted", weighted_vals, "#ef4444"});

  return groupedBarChart(terms, series, "Energy Decomposition", "Energy Term",
                         "Value", svg_path);
}

inline bool SvgChartGenerator::similarityHeatmapChart(const std::string &csv_path,
                                                const std::string &svg_path) {
  FILE *f = fopen(csv_path.c_str(), "r");
  if (!f) return false;

  std::vector<std::string> labels;
  std::vector<std::vector<float>> matrix;
  char line[8192];

  // First line: header with cell line names
  if (!fgets(line, sizeof(line), f)) {
    fclose(f);
    return false;
  }
  // Parse header: "cell_line,name1,name2,..."
  std::istringstream hdr(line);
  std::string token;
  std::getline(hdr, token, ','); // skip "cell_line"
  while (std::getline(hdr, token, ',')) {
    while (!token.empty() && (token.back() == '\n' || token.back() == '\r'))
      token.pop_back();
    if (!token.empty())
      labels.push_back(token);
  }

  // Data rows
  while (fgets(line, sizeof(line), f)) {
    std::istringstream row(line);
    std::string name;
    std::getline(row, name, ','); // row label
    std::vector<float> vals;
    std::string val_str;
    while (std::getline(row, val_str, ',')) {
      vals.push_back(std::stof(val_str));
    }
    matrix.push_back(vals);
  }
  fclose(f);

  if (matrix.empty() || labels.empty()) return false;

  return heatmap(matrix, labels, labels, "Structural Similarity (SCC)",
                 svg_path, 0.0f, 1.0f);
}

// ---------------------------------------------------------------------------
// RMSD convergence chart from CSV: phase,rmsd
// ---------------------------------------------------------------------------
inline bool SvgChartGenerator::rmsdConvergenceChart(const std::string &csv_path,
                                              const std::string &svg_path) {
  FILE *f = fopen(csv_path.c_str(), "r");
  if (!f) return false;

  char line[1024];
  std::vector<std::string> phases;
  std::vector<float> values;

  // Skip header
  if (fgets(line, sizeof(line), f)) { /* header */ }

  while (fgets(line, sizeof(line), f)) {
    char phase[256];
    float val;
    if (sscanf(line, "%[^,],%f", phase, &val) == 2) {
      phases.push_back(phase);
      values.push_back(val);
    }
  }
  fclose(f);

  if (phases.empty()) return false;

  return barChart(phases, values,
                  "RMSD Convergence by Phase",
                  "Phase", "RMSD",
                  svg_path, "#ef4444");
}

// ---------------------------------------------------------------------------
// Rg bar chart
// ---------------------------------------------------------------------------
inline bool SvgChartGenerator::rgBarChart(const std::vector<std::string> &chr_names,
                                    const std::vector<float> &rg_values,
                                    const std::string &svg_path) {
  if (chr_names.empty() || rg_values.empty()) return false;
  return barChart(chr_names, rg_values,
                  "Radius of Gyration (Rg) by Chromosome",
                  "Chromosome", "Rg",
                  svg_path, "#10b981");
}

#endif /* SVG_CHART_GENERATOR_H_ */
