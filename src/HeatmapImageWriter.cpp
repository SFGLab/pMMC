/**
 * @file HeatmapImageWriter.cpp
 * @brief Implementation of PNG and SVG heatmap image generation.
 */

#include <HeatmapImageWriter.h>
#include <Heatmap.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// ============================================================
// Minimal 5x7 bitmap font for labels
// ============================================================
// Each character is 5 columns x 7 rows, stored as 7 bytes
// where each byte represents one row (bit 4=leftmost, bit 0=rightmost).

static const uint8_t FONT_5x7[][7] = {
    // space (32)
    {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00},
    // ! (33)
    {0x04, 0x04, 0x04, 0x04, 0x04, 0x00, 0x04},
    // " (34)
    {0x0A, 0x0A, 0x00, 0x00, 0x00, 0x00, 0x00},
    // # (35)
    {0x0A, 0x1F, 0x0A, 0x0A, 0x1F, 0x0A, 0x00},
    // $ (36)
    {0x04, 0x0F, 0x14, 0x0E, 0x05, 0x1E, 0x04},
    // % (37)
    {0x18, 0x19, 0x02, 0x04, 0x08, 0x13, 0x03},
    // & (38)
    {0x08, 0x14, 0x14, 0x08, 0x15, 0x12, 0x0D},
    // ' (39)
    {0x04, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00},
    // ( (40)
    {0x02, 0x04, 0x08, 0x08, 0x08, 0x04, 0x02},
    // ) (41)
    {0x08, 0x04, 0x02, 0x02, 0x02, 0x04, 0x08},
    // * (42)
    {0x00, 0x04, 0x15, 0x0E, 0x15, 0x04, 0x00},
    // + (43)
    {0x00, 0x04, 0x04, 0x1F, 0x04, 0x04, 0x00},
    // , (44)
    {0x00, 0x00, 0x00, 0x00, 0x00, 0x04, 0x08},
    // - (45)
    {0x00, 0x00, 0x00, 0x1F, 0x00, 0x00, 0x00},
    // . (46)
    {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x04},
    // / (47)
    {0x01, 0x01, 0x02, 0x04, 0x08, 0x10, 0x10},
    // 0 (48)
    {0x0E, 0x11, 0x13, 0x15, 0x19, 0x11, 0x0E},
    // 1 (49)
    {0x04, 0x0C, 0x04, 0x04, 0x04, 0x04, 0x0E},
    // 2 (50)
    {0x0E, 0x11, 0x01, 0x06, 0x08, 0x10, 0x1F},
    // 3 (51)
    {0x0E, 0x11, 0x01, 0x0E, 0x01, 0x11, 0x0E},
    // 4 (52)
    {0x02, 0x06, 0x0A, 0x12, 0x1F, 0x02, 0x02},
    // 5 (53)
    {0x1F, 0x10, 0x1E, 0x01, 0x01, 0x11, 0x0E},
    // 6 (54)
    {0x06, 0x08, 0x10, 0x1E, 0x11, 0x11, 0x0E},
    // 7 (55)
    {0x1F, 0x01, 0x02, 0x04, 0x08, 0x08, 0x08},
    // 8 (56)
    {0x0E, 0x11, 0x11, 0x0E, 0x11, 0x11, 0x0E},
    // 9 (57)
    {0x0E, 0x11, 0x11, 0x0F, 0x01, 0x02, 0x0C},
    // : (58)
    {0x00, 0x00, 0x04, 0x00, 0x00, 0x04, 0x00},
    // ; (59)
    {0x00, 0x00, 0x04, 0x00, 0x00, 0x04, 0x08},
    // < (60)
    {0x02, 0x04, 0x08, 0x10, 0x08, 0x04, 0x02},
    // = (61)
    {0x00, 0x00, 0x1F, 0x00, 0x1F, 0x00, 0x00},
    // > (62)
    {0x08, 0x04, 0x02, 0x01, 0x02, 0x04, 0x08},
    // ? (63)
    {0x0E, 0x11, 0x01, 0x02, 0x04, 0x00, 0x04},
    // @ (64)
    {0x0E, 0x11, 0x17, 0x15, 0x17, 0x10, 0x0E},
    // A (65)
    {0x0E, 0x11, 0x11, 0x1F, 0x11, 0x11, 0x11},
    // B (66)
    {0x1E, 0x11, 0x11, 0x1E, 0x11, 0x11, 0x1E},
    // C (67)
    {0x0E, 0x11, 0x10, 0x10, 0x10, 0x11, 0x0E},
    // D (68)
    {0x1E, 0x11, 0x11, 0x11, 0x11, 0x11, 0x1E},
    // E (69)
    {0x1F, 0x10, 0x10, 0x1E, 0x10, 0x10, 0x1F},
    // F (70)
    {0x1F, 0x10, 0x10, 0x1E, 0x10, 0x10, 0x10},
    // G (71)
    {0x0E, 0x11, 0x10, 0x17, 0x11, 0x11, 0x0E},
    // H (72)
    {0x11, 0x11, 0x11, 0x1F, 0x11, 0x11, 0x11},
    // I (73)
    {0x0E, 0x04, 0x04, 0x04, 0x04, 0x04, 0x0E},
    // J (74)
    {0x07, 0x02, 0x02, 0x02, 0x02, 0x12, 0x0C},
    // K (75)
    {0x11, 0x12, 0x14, 0x18, 0x14, 0x12, 0x11},
    // L (76)
    {0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x1F},
    // M (77)
    {0x11, 0x1B, 0x15, 0x15, 0x11, 0x11, 0x11},
    // N (78)
    {0x11, 0x19, 0x15, 0x13, 0x11, 0x11, 0x11},
    // O (79)
    {0x0E, 0x11, 0x11, 0x11, 0x11, 0x11, 0x0E},
    // P (80)
    {0x1E, 0x11, 0x11, 0x1E, 0x10, 0x10, 0x10},
    // Q (81)
    {0x0E, 0x11, 0x11, 0x11, 0x15, 0x12, 0x0D},
    // R (82)
    {0x1E, 0x11, 0x11, 0x1E, 0x14, 0x12, 0x11},
    // S (83)
    {0x0E, 0x11, 0x10, 0x0E, 0x01, 0x11, 0x0E},
    // T (84)
    {0x1F, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04},
    // U (85)
    {0x11, 0x11, 0x11, 0x11, 0x11, 0x11, 0x0E},
    // V (86)
    {0x11, 0x11, 0x11, 0x11, 0x0A, 0x0A, 0x04},
    // W (87)
    {0x11, 0x11, 0x11, 0x15, 0x15, 0x1B, 0x11},
    // X (88)
    {0x11, 0x11, 0x0A, 0x04, 0x0A, 0x11, 0x11},
    // Y (89)
    {0x11, 0x11, 0x0A, 0x04, 0x04, 0x04, 0x04},
    // Z (90)
    {0x1F, 0x01, 0x02, 0x04, 0x08, 0x10, 0x1F},
    // [ (91)
    {0x0E, 0x08, 0x08, 0x08, 0x08, 0x08, 0x0E},
    // backslash (92)
    {0x10, 0x10, 0x08, 0x04, 0x02, 0x01, 0x01},
    // ] (93)
    {0x0E, 0x02, 0x02, 0x02, 0x02, 0x02, 0x0E},
    // ^ (94)
    {0x04, 0x0A, 0x11, 0x00, 0x00, 0x00, 0x00},
    // _ (95)
    {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x1F},
    // ` (96)
    {0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00},
    // a (97)
    {0x00, 0x00, 0x0E, 0x01, 0x0F, 0x11, 0x0F},
    // b (98)
    {0x10, 0x10, 0x1E, 0x11, 0x11, 0x11, 0x1E},
    // c (99)
    {0x00, 0x00, 0x0E, 0x11, 0x10, 0x11, 0x0E},
    // d (100)
    {0x01, 0x01, 0x0F, 0x11, 0x11, 0x11, 0x0F},
    // e (101)
    {0x00, 0x00, 0x0E, 0x11, 0x1F, 0x10, 0x0E},
    // f (102)
    {0x06, 0x08, 0x1E, 0x08, 0x08, 0x08, 0x08},
    // g (103)
    {0x00, 0x00, 0x0F, 0x11, 0x0F, 0x01, 0x0E},
    // h (104)
    {0x10, 0x10, 0x1E, 0x11, 0x11, 0x11, 0x11},
    // i (105)
    {0x04, 0x00, 0x0C, 0x04, 0x04, 0x04, 0x0E},
    // j (106)
    {0x02, 0x00, 0x06, 0x02, 0x02, 0x12, 0x0C},
    // k (107)
    {0x10, 0x10, 0x12, 0x14, 0x18, 0x14, 0x12},
    // l (108)
    {0x0C, 0x04, 0x04, 0x04, 0x04, 0x04, 0x0E},
    // m (109)
    {0x00, 0x00, 0x1A, 0x15, 0x15, 0x15, 0x15},
    // n (110)
    {0x00, 0x00, 0x1E, 0x11, 0x11, 0x11, 0x11},
    // o (111)
    {0x00, 0x00, 0x0E, 0x11, 0x11, 0x11, 0x0E},
    // p (112)
    {0x00, 0x00, 0x1E, 0x11, 0x1E, 0x10, 0x10},
    // q (113)
    {0x00, 0x00, 0x0F, 0x11, 0x0F, 0x01, 0x01},
    // r (114)
    {0x00, 0x00, 0x16, 0x19, 0x10, 0x10, 0x10},
    // s (115)
    {0x00, 0x00, 0x0F, 0x10, 0x0E, 0x01, 0x1E},
    // t (116)
    {0x08, 0x08, 0x1E, 0x08, 0x08, 0x09, 0x06},
    // u (117)
    {0x00, 0x00, 0x11, 0x11, 0x11, 0x11, 0x0F},
    // v (118)
    {0x00, 0x00, 0x11, 0x11, 0x0A, 0x0A, 0x04},
    // w (119)
    {0x00, 0x00, 0x11, 0x11, 0x15, 0x15, 0x0A},
    // x (120)
    {0x00, 0x00, 0x11, 0x0A, 0x04, 0x0A, 0x11},
    // y (121)
    {0x00, 0x00, 0x11, 0x11, 0x0F, 0x01, 0x0E},
    // z (122)
    {0x00, 0x00, 0x1F, 0x02, 0x04, 0x08, 0x1F},
};

static const int FONT_FIRST_CHAR = 32;
static const int FONT_LAST_CHAR = 122;
static const int FONT_W = 5;
static const int FONT_H = 7;

// ============================================================
// Color mapping: blue-white-red
// ============================================================

HeatmapImageWriter::Color HeatmapImageWriter::valueToColor(float t) {
    // Clamp to [0,1]
    if (t < 0.0f) t = 0.0f;
    if (t > 1.0f) t = 1.0f;

    Color c;
    if (t < 0.5f) {
        // Blue (0) to White (0.5)
        float s = t / 0.5f; // 0..1
        c.r = static_cast<uint8_t>(s * 255.0f);
        c.g = static_cast<uint8_t>(s * 255.0f);
        c.b = 255;
    } else {
        // White (0.5) to Red (1)
        float s = (t - 0.5f) / 0.5f; // 0..1
        c.r = 255;
        c.g = static_cast<uint8_t>((1.0f - s) * 255.0f);
        c.b = static_cast<uint8_t>((1.0f - s) * 255.0f);
    }
    return c;
}

// ============================================================
// Bitmap font rendering
// ============================================================

void HeatmapImageWriter::drawChar(std::vector<uint8_t> &pixels, int imgW,
                                   int imgH, int x0, int y0, char ch,
                                   Color color, int scale) {
    int idx = static_cast<int>(ch) - FONT_FIRST_CHAR;
    if (idx < 0 || idx > (FONT_LAST_CHAR - FONT_FIRST_CHAR))
        return;

    const uint8_t *glyph = FONT_5x7[idx];
    for (int row = 0; row < FONT_H; ++row) {
        uint8_t bits = glyph[row];
        for (int col = 0; col < FONT_W; ++col) {
            if (bits & (1 << (FONT_W - 1 - col))) {
                // Fill scale x scale block
                for (int sy = 0; sy < scale; ++sy) {
                    for (int sx = 0; sx < scale; ++sx) {
                        int px = x0 + col * scale + sx;
                        int py = y0 + row * scale + sy;
                        if (px >= 0 && px < imgW && py >= 0 && py < imgH) {
                            int off = (py * imgW + px) * 3;
                            pixels[off + 0] = color.r;
                            pixels[off + 1] = color.g;
                            pixels[off + 2] = color.b;
                        }
                    }
                }
            }
        }
    }
}

void HeatmapImageWriter::drawString(std::vector<uint8_t> &pixels, int imgW,
                                     int imgH, int x0, int y0,
                                     const std::string &text, Color color,
                                     int scale) {
    int charW = (FONT_W + 1) * scale; // 1 pixel gap between chars
    for (size_t i = 0; i < text.size(); ++i) {
        drawChar(pixels, imgW, imgH, x0 + static_cast<int>(i) * charW, y0,
                 text[i], color, scale);
    }
}

// ============================================================
// CRC-32 for PNG (IEEE polynomial)
// ============================================================

static uint32_t crc32_table[256];
static bool crc32_initialized = false;

static void init_crc32() {
    if (crc32_initialized) return;
    for (uint32_t i = 0; i < 256; ++i) {
        uint32_t c = i;
        for (int k = 0; k < 8; ++k) {
            if (c & 1)
                c = 0xEDB88320u ^ (c >> 1);
            else
                c >>= 1;
        }
        crc32_table[i] = c;
    }
    crc32_initialized = true;
}

static uint32_t update_crc32(uint32_t crc, const uint8_t *data, size_t len) {
    for (size_t i = 0; i < len; ++i)
        crc = crc32_table[(crc ^ data[i]) & 0xFF] ^ (crc >> 8);
    return crc;
}

// ============================================================
// Adler-32 checksum for zlib wrapper
// ============================================================

static uint32_t adler32(const uint8_t *data, size_t len) {
    uint32_t a = 1, b = 0;
    for (size_t i = 0; i < len; ++i) {
        a = (a + data[i]) % 65521;
        b = (b + a) % 65521;
    }
    return (b << 16) | a;
}

// ============================================================
// Minimal PNG writer using uncompressed DEFLATE store blocks
// ============================================================

// Helper: write a 4-byte big-endian uint32
static void write_be32(std::vector<uint8_t> &out, uint32_t val) {
    out.push_back(static_cast<uint8_t>((val >> 24) & 0xFF));
    out.push_back(static_cast<uint8_t>((val >> 16) & 0xFF));
    out.push_back(static_cast<uint8_t>((val >> 8) & 0xFF));
    out.push_back(static_cast<uint8_t>(val & 0xFF));
}

// Write a PNG chunk (type + data, auto-computes length and CRC)
static void write_chunk(std::vector<uint8_t> &out, const char type[4],
                         const uint8_t *data, uint32_t len) {
    init_crc32();
    write_be32(out, len);
    // Type bytes
    for (int i = 0; i < 4; ++i)
        out.push_back(static_cast<uint8_t>(type[i]));

    // CRC covers type + data
    uint32_t crc = 0xFFFFFFFFu;
    crc = update_crc32(crc, reinterpret_cast<const uint8_t *>(type), 4);
    if (len > 0) {
        for (uint32_t i = 0; i < len; ++i)
            out.push_back(data[i]);
        crc = update_crc32(crc, data, len);
    }
    crc ^= 0xFFFFFFFFu;
    write_be32(out, crc);
}

bool HeatmapImageWriter::writeRawPNG(const std::string &filename,
                                      const std::vector<uint8_t> &pixels,
                                      int width, int height) {
    std::vector<uint8_t> png;

    // PNG signature
    static const uint8_t sig[] = {137, 80, 78, 71, 13, 10, 26, 10};
    png.insert(png.end(), sig, sig + 8);

    // IHDR chunk: 13 bytes
    {
        uint8_t ihdr[13];
        // Width (4 bytes BE)
        ihdr[0] = (width >> 24) & 0xFF;
        ihdr[1] = (width >> 16) & 0xFF;
        ihdr[2] = (width >> 8) & 0xFF;
        ihdr[3] = width & 0xFF;
        // Height (4 bytes BE)
        ihdr[4] = (height >> 24) & 0xFF;
        ihdr[5] = (height >> 16) & 0xFF;
        ihdr[6] = (height >> 8) & 0xFF;
        ihdr[7] = height & 0xFF;
        ihdr[8] = 8;  // bit depth
        ihdr[9] = 2;  // color type = RGB
        ihdr[10] = 0; // compression
        ihdr[11] = 0; // filter
        ihdr[12] = 0; // interlace
        write_chunk(png, "IHDR", ihdr, 13);
    }

    // IDAT chunk: Build raw image data with filter byte 0 per row,
    // then wrap in uncompressed DEFLATE blocks inside a zlib stream.

    // Raw data: for each row, 1 filter byte (0) + width*3 RGB bytes
    size_t raw_row_size = 1 + static_cast<size_t>(width) * 3;
    size_t raw_size = raw_row_size * static_cast<size_t>(height);

    std::vector<uint8_t> raw_data;
    raw_data.reserve(raw_size);
    for (int y = 0; y < height; ++y) {
        raw_data.push_back(0); // filter: None
        const uint8_t *row = &pixels[y * width * 3];
        raw_data.insert(raw_data.end(), row, row + width * 3);
    }

    // Build zlib stream with uncompressed DEFLATE store blocks.
    // zlib header: CMF=0x78 FLG=0x01 (no dict, lowest compression)
    std::vector<uint8_t> zlib_stream;
    zlib_stream.push_back(0x78);
    zlib_stream.push_back(0x01);

    // Split raw_data into DEFLATE store blocks of up to 65535 bytes
    const size_t max_block = 65535;
    size_t offset = 0;
    while (offset < raw_data.size()) {
        size_t remaining = raw_data.size() - offset;
        size_t block_len = (remaining > max_block) ? max_block : remaining;
        bool last = (offset + block_len >= raw_data.size());

        // BFINAL | BTYPE=00 (store)
        zlib_stream.push_back(last ? 0x01 : 0x00);
        // LEN (little-endian 16-bit)
        uint16_t len16 = static_cast<uint16_t>(block_len);
        zlib_stream.push_back(len16 & 0xFF);
        zlib_stream.push_back((len16 >> 8) & 0xFF);
        // NLEN (one's complement)
        uint16_t nlen16 = static_cast<uint16_t>(~len16);
        zlib_stream.push_back(nlen16 & 0xFF);
        zlib_stream.push_back((nlen16 >> 8) & 0xFF);
        // Literal data
        zlib_stream.insert(zlib_stream.end(),
                           raw_data.begin() + offset,
                           raw_data.begin() + offset + block_len);
        offset += block_len;
    }

    // Adler-32 checksum of uncompressed data (big-endian)
    uint32_t adler = adler32(raw_data.data(), raw_data.size());
    zlib_stream.push_back((adler >> 24) & 0xFF);
    zlib_stream.push_back((adler >> 16) & 0xFF);
    zlib_stream.push_back((adler >> 8) & 0xFF);
    zlib_stream.push_back(adler & 0xFF);

    write_chunk(png, "IDAT", zlib_stream.data(),
                static_cast<uint32_t>(zlib_stream.size()));

    // IEND chunk
    write_chunk(png, "IEND", nullptr, 0);

    // Write to file
    FILE *f = fopen(filename.c_str(), "wb");
    if (!f) return false;
    size_t written = fwrite(png.data(), 1, png.size(), f);
    fclose(f);
    return written == png.size();
}

// ============================================================
// PNG heatmap writer
// ============================================================

bool HeatmapImageWriter::writePNG(const Heatmap &heat,
                                   const std::string &filename,
                                   const std::string &title) {
    if (heat.size == 0) return false;

    const int n = static_cast<int>(heat.size);

    // Layout constants
    const int cellSize = std::max(1, 800 / n); // at least 1 pixel per cell
    const int matW = cellSize * n;
    const int matH = cellSize * n;

    const int fontScale = 2;
    const int charW = (FONT_W + 1) * fontScale;
    const int charH = FONT_H * fontScale;

    const int marginLeft = 60;      // space for Y axis labels
    const int marginTop = title.empty() ? 30 : (charH * 3 + 20);
    const int marginRight = 100;    // space for colorbar
    const int marginBottom = 60;    // space for X axis labels

    const int imgW = marginLeft + matW + marginRight;
    const int imgH = marginTop + matH + marginBottom;

    // Allocate white background
    std::vector<uint8_t> pixels(imgW * imgH * 3, 255);
    Color black = {0, 0, 0};

    // Compute value range
    float vmin = 1e30f, vmax = -1e30f;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            float val = heat.v[i][j];
            if (val < vmin) vmin = val;
            if (val > vmax) vmax = val;
        }
    }
    float range = vmax - vmin;
    if (range < 1e-12f) range = 1.0f;

    // Draw heatmap cells
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            float t = (heat.v[i][j] - vmin) / range;
            Color c = valueToColor(t);
            int px0 = marginLeft + j * cellSize;
            int py0 = marginTop + i * cellSize;
            for (int dy = 0; dy < cellSize; ++dy) {
                for (int dx = 0; dx < cellSize; ++dx) {
                    int px = px0 + dx;
                    int py = py0 + dy;
                    if (px < imgW && py < imgH) {
                        int off = (py * imgW + px) * 3;
                        pixels[off + 0] = c.r;
                        pixels[off + 1] = c.g;
                        pixels[off + 2] = c.b;
                    }
                }
            }
        }
    }

    // Draw title
    if (!title.empty()) {
        int titleX = marginLeft + (matW - static_cast<int>(title.size()) * charW * 2) / 2;
        if (titleX < 2) titleX = 2;
        drawString(pixels, imgW, imgH, titleX, 8, title, black, fontScale * 2);
    }

    // Axis labels: show a subset of bin indices to avoid clutter
    int labelStep = 1;
    if (n > 40) labelStep = 5;
    if (n > 100) labelStep = 10;
    if (n > 500) labelStep = 50;
    if (n > 1000) labelStep = 100;

    // X axis labels (bottom)
    for (int j = 0; j < n; j += labelStep) {
        std::string lbl = std::to_string(j);
        int lx = marginLeft + j * cellSize + cellSize / 2 -
                 static_cast<int>(lbl.size()) * charW / 2;
        int ly = marginTop + matH + 8;
        drawString(pixels, imgW, imgH, lx, ly, lbl, black, fontScale);
    }

    // Y axis labels (left)
    for (int i = 0; i < n; i += labelStep) {
        std::string lbl = std::to_string(i);
        int lx = marginLeft - static_cast<int>(lbl.size()) * charW - 4;
        if (lx < 0) lx = 0;
        int ly = marginTop + i * cellSize + cellSize / 2 - charH / 2;
        drawString(pixels, imgW, imgH, lx, ly, lbl, black, fontScale);
    }

    // Colorbar
    int cbX = marginLeft + matW + 20;
    int cbW = 20;
    int cbH = matH;
    int cbY = marginTop;

    for (int y = 0; y < cbH; ++y) {
        float t = 1.0f - static_cast<float>(y) / static_cast<float>(cbH - 1);
        Color c = valueToColor(t);
        for (int x = 0; x < cbW; ++x) {
            int px = cbX + x;
            int py = cbY + y;
            if (px < imgW && py < imgH) {
                int off = (py * imgW + px) * 3;
                pixels[off + 0] = c.r;
                pixels[off + 1] = c.g;
                pixels[off + 2] = c.b;
            }
        }
    }

    // Colorbar border
    for (int y = 0; y < cbH; ++y) {
        int offL = ((cbY + y) * imgW + cbX) * 3;
        int offR = ((cbY + y) * imgW + cbX + cbW - 1) * 3;
        if (cbY + y < imgH) {
            pixels[offL] = pixels[offL + 1] = pixels[offL + 2] = 0;
            pixels[offR] = pixels[offR + 1] = pixels[offR + 2] = 0;
        }
    }
    for (int x = 0; x < cbW; ++x) {
        int offT = (cbY * imgW + cbX + x) * 3;
        int offB = ((cbY + cbH - 1) * imgW + cbX + x) * 3;
        if (cbX + x < imgW) {
            pixels[offT] = pixels[offT + 1] = pixels[offT + 2] = 0;
            pixels[offB] = pixels[offB + 1] = pixels[offB + 2] = 0;
        }
    }

    // Colorbar labels (min at bottom, max at top)
    {
        char buf[32];
        snprintf(buf, sizeof(buf), "%.2g", vmax);
        drawString(pixels, imgW, imgH, cbX + cbW + 4, cbY, std::string(buf),
                   black, fontScale);
        snprintf(buf, sizeof(buf), "%.2g", vmin);
        drawString(pixels, imgW, imgH, cbX + cbW + 4, cbY + cbH - charH,
                   std::string(buf), black, fontScale);
        // Mid value
        snprintf(buf, sizeof(buf), "%.2g", (vmin + vmax) / 2.0f);
        drawString(pixels, imgW, imgH, cbX + cbW + 4, cbY + cbH / 2 - charH / 2,
                   std::string(buf), black, fontScale);
    }

    return writeRawPNG(filename, pixels, imgW, imgH);
}

// ============================================================
// SVG heatmap writer
// ============================================================

static std::string colorToSVG(HeatmapImageWriter::Color c) {
    char buf[24];
    snprintf(buf, sizeof(buf), "rgb(%d,%d,%d)", c.r, c.g, c.b);
    return std::string(buf);
}

bool HeatmapImageWriter::writeSVG(const Heatmap &heat,
                                   const std::string &filename,
                                   const std::string &title) {
    if (heat.size == 0) return false;

    const int n = static_cast<int>(heat.size);

    // Layout
    const int cellSize = std::max(2, 800 / n);
    const int matW = cellSize * n;
    const int matH = cellSize * n;

    const int marginLeft = 60;
    const int marginTop = title.empty() ? 30 : 60;
    const int marginRight = 120;
    const int marginBottom = 60;

    const int svgW = marginLeft + matW + marginRight;
    const int svgH = marginTop + matH + marginBottom;

    // Compute range
    float vmin = 1e30f, vmax = -1e30f;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            float val = heat.v[i][j];
            if (val < vmin) vmin = val;
            if (val > vmax) vmax = val;
        }
    }
    float range = vmax - vmin;
    if (range < 1e-12f) range = 1.0f;

    std::ostringstream svg;
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
        << "width=\"" << svgW << "\" height=\"" << svgH << "\" "
        << "viewBox=\"0 0 " << svgW << " " << svgH << "\">\n";

    // White background
    svg << "<rect width=\"" << svgW << "\" height=\"" << svgH
        << "\" fill=\"white\"/>\n";

    // Title
    if (!title.empty()) {
        int titleX = marginLeft + matW / 2;
        svg << "<text x=\"" << titleX << "\" y=\"" << 30
            << "\" font-family=\"monospace\" font-size=\"18\" "
            << "text-anchor=\"middle\" font-weight=\"bold\">"
            << title << "</text>\n";
    }

    // Heatmap cells
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            float t = (heat.v[i][j] - vmin) / range;
            Color c = valueToColor(t);
            int x = marginLeft + j * cellSize;
            int y = marginTop + i * cellSize;
            svg << "<rect x=\"" << x << "\" y=\"" << y
                << "\" width=\"" << cellSize << "\" height=\"" << cellSize
                << "\" fill=\"" << colorToSVG(c) << "\"/>\n";
        }
    }

    // Axis labels
    int labelStep = 1;
    if (n > 40) labelStep = 5;
    if (n > 100) labelStep = 10;
    if (n > 500) labelStep = 50;
    if (n > 1000) labelStep = 100;

    int fontSize = 10;
    if (cellSize * labelStep < 14) fontSize = 8;

    // X axis labels
    for (int j = 0; j < n; j += labelStep) {
        int x = marginLeft + j * cellSize + cellSize / 2;
        int y = marginTop + matH + 20;
        svg << "<text x=\"" << x << "\" y=\"" << y
            << "\" font-family=\"monospace\" font-size=\"" << fontSize
            << "\" text-anchor=\"middle\">" << j << "</text>\n";
    }

    // Y axis labels
    for (int i = 0; i < n; i += labelStep) {
        int x = marginLeft - 6;
        int y = marginTop + i * cellSize + cellSize / 2 + fontSize / 3;
        svg << "<text x=\"" << x << "\" y=\"" << y
            << "\" font-family=\"monospace\" font-size=\"" << fontSize
            << "\" text-anchor=\"end\">" << i << "</text>\n";
    }

    // Colorbar
    int cbX = marginLeft + matW + 20;
    int cbW = 20;
    int cbH = matH;
    int cbY = marginTop;
    int cbSteps = std::min(cbH, 256);

    // Colorbar gradient as stacked rects
    for (int s = 0; s < cbSteps; ++s) {
        float t = 1.0f - static_cast<float>(s) / static_cast<float>(cbSteps - 1);
        Color c = valueToColor(t);
        int ry = cbY + s * cbH / cbSteps;
        int rh = (s + 1) * cbH / cbSteps - s * cbH / cbSteps;
        svg << "<rect x=\"" << cbX << "\" y=\"" << ry
            << "\" width=\"" << cbW << "\" height=\"" << rh
            << "\" fill=\"" << colorToSVG(c) << "\"/>\n";
    }

    // Colorbar border
    svg << "<rect x=\"" << cbX << "\" y=\"" << cbY
        << "\" width=\"" << cbW << "\" height=\"" << cbH
        << "\" fill=\"none\" stroke=\"black\" stroke-width=\"1\"/>\n";

    // Colorbar labels
    {
        char buf[32];
        snprintf(buf, sizeof(buf), "%.4g", vmax);
        svg << "<text x=\"" << (cbX + cbW + 4) << "\" y=\"" << (cbY + fontSize)
            << "\" font-family=\"monospace\" font-size=\"" << fontSize
            << "\">" << buf << "</text>\n";

        snprintf(buf, sizeof(buf), "%.4g", (vmin + vmax) / 2.0f);
        svg << "<text x=\"" << (cbX + cbW + 4) << "\" y=\""
            << (cbY + cbH / 2 + fontSize / 3)
            << "\" font-family=\"monospace\" font-size=\"" << fontSize
            << "\">" << buf << "</text>\n";

        snprintf(buf, sizeof(buf), "%.4g", vmin);
        svg << "<text x=\"" << (cbX + cbW + 4) << "\" y=\"" << (cbY + cbH)
            << "\" font-family=\"monospace\" font-size=\"" << fontSize
            << "\">" << buf << "</text>\n";
    }

    svg << "</svg>\n";

    // Write file
    std::ofstream ofs(filename);
    if (!ofs.is_open()) return false;
    ofs << svg.str();
    ofs.close();
    return true;
}
