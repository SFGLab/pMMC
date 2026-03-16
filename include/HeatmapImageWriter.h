/**
 * @file HeatmapImageWriter.h
 * @brief PNG and SVG image generation for Heatmap contact maps.
 *
 * Produces publication-quality heatmap images with colorbar, axis labels,
 * and title. PNG uses an embedded minimal encoder (no external dependencies).
 * SVG generates standards-compliant XML directly.
 */

#ifndef HEATMAPIMAGEWRITER_H_
#define HEATMAPIMAGEWRITER_H_

#include <cstdint>
#include <string>
#include <vector>

class Heatmap;

/**
 * @class HeatmapImageWriter
 * @brief Writes Heatmap objects as PNG or SVG image files.
 *
 * Both formats include a blue-white-red color gradient, a colorbar on
 * the right side, axis labels (bin indices), and an optional title.
 */
class HeatmapImageWriter {
public:
  /**
   * @brief Write the heatmap as a PNG image.
   * @param heat   The heatmap to render.
   * @param filename Output file path (should end in .png).
   * @param title  Optional title displayed above the heatmap.
   * @return True on success, false on I/O or encoding error.
   */
  static bool writePNG(const Heatmap &heat, const std::string &filename,
                       const std::string &title = "");

  /**
   * @brief Write the heatmap as an SVG image.
   * @param heat   The heatmap to render.
   * @param filename Output file path (should end in .svg).
   * @param title  Optional title displayed above the heatmap.
   * @return True on success, false on I/O error.
   */
  static bool writeSVG(const Heatmap &heat, const std::string &filename,
                       const std::string &title = "");

  /** @brief RGB color triplet (public for use by helpers). */
  struct Color {
    uint8_t r, g, b;
  };

  /**
   * @brief Map a normalized value [0,1] to a blue-white-red color.
   * @param t Value in [0,1].
   * @return Corresponding RGB color.
   */
  static Color valueToColor(float t);

private:
  /**
   * @brief Render a character glyph into an RGB pixel buffer.
   * @param pixels  RGB buffer (3 bytes per pixel).
   * @param imgW    Image width in pixels.
   * @param imgH    Image height in pixels.
   * @param x0      Top-left x of the glyph bounding box.
   * @param y0      Top-left y of the glyph bounding box.
   * @param ch      Character to render.
   * @param color   Foreground color.
   * @param scale   Font scale multiplier (1 = 5x7 bitmap).
   */
  static void drawChar(std::vector<uint8_t> &pixels, int imgW, int imgH,
                       int x0, int y0, char ch, Color color, int scale = 1);

  /**
   * @brief Render a string into an RGB pixel buffer.
   * @param pixels  RGB buffer.
   * @param imgW    Image width.
   * @param imgH    Image height.
   * @param x0      Starting x position.
   * @param y0      Starting y position.
   * @param text    String to render.
   * @param color   Foreground color.
   * @param scale   Font scale multiplier.
   */
  static void drawString(std::vector<uint8_t> &pixels, int imgW, int imgH,
                          int x0, int y0, const std::string &text,
                          Color color, int scale = 1);

  /**
   * @brief Write raw RGB pixel data as a valid PNG file.
   *
   * Uses uncompressed DEFLATE (store blocks) so no zlib dependency
   * is required. Produces a fully spec-compliant PNG.
   *
   * @param filename Output path.
   * @param pixels   RGB pixel data (width*height*3 bytes).
   * @param width    Image width.
   * @param height   Image height.
   * @return True on success.
   */
  static bool writeRawPNG(const std::string &filename,
                           const std::vector<uint8_t> &pixels,
                           int width, int height);
};

#endif /* HEATMAPIMAGEWRITER_H_ */
