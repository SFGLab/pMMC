#ifndef FLOWCHART_GENERATOR_H
#define FLOWCHART_GENERATOR_H

#include <string>

class FlowchartGenerator {
public:
    static bool generate(const std::string& output_path);
};


// ============================================================================
// Implementation
// ============================================================================

#include <fstream>
#include <sstream>
#include <string>

// ---------------------------------------------------------------------------
// Helper: build an SVG arrowhead marker definition
// ---------------------------------------------------------------------------
static std::string arrowMarker(const std::string& id,
                               const std::string& color) {
    std::ostringstream m;
    m << "    <marker id=\"" << id << "\" markerWidth=\"10\" markerHeight=\"7\" "
      << "refX=\"10\" refY=\"3.5\" orient=\"auto\" markerUnits=\"strokeWidth\">\n"
      << "      <polygon points=\"0 0, 10 3.5, 0 7\" fill=\"" << color << "\"/>\n"
      << "    </marker>\n";
    return m.str();
}

// ---------------------------------------------------------------------------
// Helper: a rounded rectangle process box with centred multi-line text
// ---------------------------------------------------------------------------
static std::string processBox(double cx, double cy, double w, double h,
                              const std::string& fill,
                              const std::string& stroke,
                              const std::string lines[], int nlines,
                              const int fontSizes[] = nullptr) {
    std::ostringstream s;
    double x = cx - w / 2.0;
    double y = cy - h / 2.0;
    s << "  <rect x=\"" << x << "\" y=\"" << y
      << "\" width=\"" << w << "\" height=\"" << h
      << "\" rx=\"12\" ry=\"12\" fill=\"" << fill
      << "\" stroke=\"" << stroke << "\" stroke-width=\"2\"/>\n";

    // Lay out text lines centred vertically
    double lineH = 20.0;
    double startY = cy - ((nlines - 1) * lineH) / 2.0;
    for (int i = 0; i < nlines; ++i) {
        int fs = (fontSizes ? fontSizes[i] : 14);
        std::string weight = (i == 0) ? "bold" : "normal";
        s << "  <text x=\"" << cx << "\" y=\"" << (startY + i * lineH)
          << "\" text-anchor=\"middle\" dominant-baseline=\"central\" "
          << "font-family=\"Arial, Helvetica, sans-serif\" font-size=\""
          << fs << "\" font-weight=\"" << weight
          << "\" fill=\"#222\">" << lines[i] << "</text>\n";
    }
    return s.str();
}

// ---------------------------------------------------------------------------
// Helper: a diamond decision shape with centred text
// ---------------------------------------------------------------------------
static std::string diamond(double cx, double cy, double w, double h,
                           const std::string& fill,
                           const std::string& stroke,
                           const std::string lines[], int nlines) {
    std::ostringstream s;
    // Diamond via polygon
    s << "  <polygon points=\""
      << cx << "," << (cy - h / 2.0) << " "
      << (cx + w / 2.0) << "," << cy << " "
      << cx << "," << (cy + h / 2.0) << " "
      << (cx - w / 2.0) << "," << cy
      << "\" fill=\"" << fill << "\" stroke=\"" << stroke
      << "\" stroke-width=\"2\"/>\n";

    double lineH = 17.0;
    double startY = cy - ((nlines - 1) * lineH) / 2.0;
    for (int i = 0; i < nlines; ++i) {
        s << "  <text x=\"" << cx << "\" y=\"" << (startY + i * lineH)
          << "\" text-anchor=\"middle\" dominant-baseline=\"central\" "
          << "font-family=\"Arial, Helvetica, sans-serif\" font-size=\"13\" "
          << "font-weight=\"bold\" fill=\"#222\">"
          << lines[i] << "</text>\n";
    }
    return s.str();
}

// ---------------------------------------------------------------------------
// Helper: a vertical arrow between two y-positions at the given x
// ---------------------------------------------------------------------------
static std::string arrow(double x, double y1, double y2,
                         const std::string& markerId,
                         const std::string& color) {
    std::ostringstream s;
    s << "  <line x1=\"" << x << "\" y1=\"" << y1
      << "\" x2=\"" << x << "\" y2=\"" << y2
      << "\" stroke=\"" << color << "\" stroke-width=\"2\" "
      << "marker-end=\"url(#" << markerId << ")\"/>\n";
    return s.str();
}

// ---------------------------------------------------------------------------
// Helper: a right-angle loop-back arrow (goes right, up, then left back in)
// ---------------------------------------------------------------------------
static std::string loopbackArrow(double xCenter, double yStart,
                                 double yEnd, double xOffset,
                                 const std::string& markerId,
                                 const std::string& color) {
    std::ostringstream s;
    double xRight = xCenter + xOffset;
    s << "  <path d=\"M " << xCenter << " " << yStart
      << " L " << xRight << " " << yStart
      << " L " << xRight << " " << yEnd
      << " L " << xCenter << " " << yEnd
      << "\" fill=\"none\" stroke=\"" << color
      << "\" stroke-width=\"2\" marker-end=\"url(#" << markerId << ")\"/>\n";
    return s.str();
}

// ---------------------------------------------------------------------------
// FlowchartGenerator::generate
// ---------------------------------------------------------------------------
inline bool FlowchartGenerator::generate(const std::string& output_path) {
    // Layout constants
    const double W       = 800.0;   // canvas width
    const double boxW    = 360.0;   // process box width
    const double boxH    = 68.0;    // process box height
    const double diaW    = 260.0;   // diamond width
    const double diaH    = 80.0;    // diamond height
    const double gap     = 32.0;    // vertical gap between shapes
    const double cx      = W / 2.0; // centre x

    // Colours per hierarchy level
    const std::string colInput   = "#dbeafe"; // light blue
    const std::string colTree    = "#e0e7ff"; // light indigo
    const std::string colL01     = "#fef3c7"; // amber
    const std::string colSphere  = "#fce7f3"; // pink
    const std::string colL2      = "#d1fae5"; // green
    const std::string colL3      = "#ede9fe"; // violet
    const std::string colOutput  = "#ffedd5"; // orange
    const std::string colDecis   = "#f3f4f6"; // gray
    const std::string colDone    = "#d1d5db"; // dark gray
    const std::string strokeClr  = "#475569"; // slate

    // Y positions for each element (top to bottom)
    double y = 50.0;  // title
    double yTitle = y;

    y += 55.0;  // input box top-centre
    double yInput = y;

    y += boxH / 2.0 + gap;
    y += boxH / 2.0;
    double yTree = y;

    y += boxH / 2.0 + gap;
    y += diaH / 2.0;
    double yEnsStart = y;

    y += diaH / 2.0 + gap;
    y += boxH / 2.0;
    double yL01 = y;

    y += boxH / 2.0 + gap;
    y += diaH / 2.0;
    double ySphere = y;

    y += diaH / 2.0 + gap;
    y += boxH / 2.0;
    double yL2 = y;

    y += boxH / 2.0 + gap;
    y += boxH / 2.0;
    double yL3 = y;

    y += boxH / 2.0 + gap;
    y += boxH / 2.0;
    double yOut = y;

    y += boxH / 2.0 + gap;
    y += diaH / 2.0;
    double yEnsEnd = y;

    y += diaH / 2.0 + gap;
    y += 20.0;       // done pill
    double yDone = y;

    double H = yDone + 50.0; // canvas height

    // Start building SVG
    std::ostringstream svg;
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
        << "width=\"" << W << "\" height=\"" << H << "\" "
        << "viewBox=\"0 0 " << W << " " << H << "\">\n";

    // Background
    svg << "  <rect width=\"" << W << "\" height=\"" << H
        << "\" fill=\"#ffffff\"/>\n";

    // Defs: arrowheads
    svg << "  <defs>\n";
    svg << arrowMarker("arrowGray", strokeClr);
    svg << arrowMarker("arrowGreen", "#16a34a");
    svg << "  </defs>\n";

    // Title
    svg << "  <text x=\"" << cx << "\" y=\"" << yTitle
        << "\" text-anchor=\"middle\" dominant-baseline=\"central\" "
        << "font-family=\"Arial, Helvetica, sans-serif\" font-size=\"24\" "
        << "font-weight=\"bold\" fill=\"#1e293b\">pMMC Algorithm Flowchart</text>\n";

    // ---- Boxes ----

    // 1. Input Data Loading
    {
        std::string lines[] = {
            "Input Data Loading",
            "BED anchors, BEDPE clusters,",
            "singletons, heatmaps"
        };
        svg << processBox(cx, yInput, boxW, boxH, colInput, strokeClr, lines, 3);
    }

    // Arrow input -> tree
    svg << arrow(cx, yInput + boxH / 2.0, yTree - boxH / 2.0, "arrowGray", strokeClr);

    // 2. Build Hierarchical Tree
    {
        std::string lines[] = {
            "Build Hierarchical Tree",
            "4 levels: Chr &gt; Segment &gt; Anchor &gt; Subanchor"
        };
        svg << processBox(cx, yTree, boxW, boxH, colTree, strokeClr, lines, 2);
    }

    // Arrow tree -> ensemble start
    svg << arrow(cx, yTree + boxH / 2.0, yEnsStart - diaH / 2.0, "arrowGray", strokeClr);

    // 3. Ensemble Loop Start (diamond)
    {
        std::string lines[] = { "Ensemble Loop", "m &gt; 1?" };
        svg << diamond(cx, yEnsStart, diaW, diaH, colDecis, strokeClr, lines, 2);
    }

    // Arrow ensemble start -> L0-1
    svg << arrow(cx, yEnsStart + diaH / 2.0, yL01 - boxH / 2.0, "arrowGray", strokeClr);

    // 4. Level 0-1: Heatmap MC
    {
        std::string lines[] = {
            "Level 0-1: Heatmap Monte Carlo",
            "Chromosome &amp; Segment scale",
            "Simulated annealing"
        };
        svg << processBox(cx, yL01, boxW, boxH, colL01, strokeClr, lines, 3);
    }

    // Arrow L0-1 -> sphere
    svg << arrow(cx, yL01 + boxH / 2.0, ySphere - diaH / 2.0, "arrowGray", strokeClr);

    // 5. Spherical Container (diamond)
    {
        std::string lines[] = { "Spherical", "Container?" };
        svg << diamond(cx, ySphere, diaW, diaH, colSphere, "#db2777", lines, 2);
    }
    // Side label "If enabled: enforce boundary penalty"
    svg << "  <text x=\"" << (cx + diaW / 2.0 + 12) << "\" y=\"" << ySphere
        << "\" text-anchor=\"start\" dominant-baseline=\"central\" "
        << "font-family=\"Arial, Helvetica, sans-serif\" font-size=\"12\" "
        << "font-style=\"italic\" fill=\"#9f1239\">"
        << "If enabled: enforce boundary penalty</text>\n";

    // Arrow sphere -> L2
    svg << arrow(cx, ySphere + diaH / 2.0, yL2 - boxH / 2.0, "arrowGray", strokeClr);

    // 6. Level 2: Arc-Distance MC
    {
        std::string lines[] = {
            "Level 2: Arc-Distance Monte Carlo",
            "Anchor scale, PET cluster arcs",
            "Simulated annealing"
        };
        svg << processBox(cx, yL2, boxW, boxH, colL2, strokeClr, lines, 3);
    }

    // Arrow L2 -> L3
    svg << arrow(cx, yL2 + boxH / 2.0, yL3 - boxH / 2.0, "arrowGray", strokeClr);

    // 7. Level 3: Smoothing MC
    {
        std::string lines[] = {
            "Level 3: Smoothing Monte Carlo",
            "Subanchor scale, polymer springs",
            "Angular + distance constraints"
        };
        svg << processBox(cx, yL3, boxW, boxH, colL3, strokeClr, lines, 3);
    }

    // Arrow L3 -> output
    svg << arrow(cx, yL3 + boxH / 2.0, yOut - boxH / 2.0, "arrowGray", strokeClr);

    // 8. Output Generation
    {
        std::string lines[] = {
            "Output Generation",
            "HCM, PDB, CIF formats",
            "Visualization scripts (.pml, .cxc)"
        };
        svg << processBox(cx, yOut, boxW, boxH, colOutput, strokeClr, lines, 3);
    }

    // Arrow output -> ensemble end
    svg << arrow(cx, yOut + boxH / 2.0, yEnsEnd - diaH / 2.0, "arrowGray", strokeClr);

    // 9. Ensemble Loop End (diamond)
    {
        std::string lines[] = { "More models?" };
        svg << diamond(cx, yEnsEnd, diaW, diaH, colDecis, strokeClr, lines, 1);
    }

    // "Yes" loop-back arrow from ensemble-end right side back up to L0-1 top
    svg << loopbackArrow(cx + diaW / 2.0, yEnsEnd,
                         yL01 - boxH / 2.0 + 5, 80,
                         "arrowGreen", "#16a34a");
    // "Yes" label
    svg << "  <text x=\"" << (cx + diaW / 2.0 + 8) << "\" y=\"" << (yEnsEnd - 8)
        << "\" text-anchor=\"start\" dominant-baseline=\"auto\" "
        << "font-family=\"Arial, Helvetica, sans-serif\" font-size=\"13\" "
        << "font-weight=\"bold\" fill=\"#16a34a\">Yes</text>\n";

    // "No" arrow down to Done
    svg << arrow(cx, yEnsEnd + diaH / 2.0, yDone - 18, "arrowGray", strokeClr);
    // "No" label
    svg << "  <text x=\"" << (cx + 10) << "\" y=\"" << (yEnsEnd + diaH / 2.0 + 16)
        << "\" text-anchor=\"start\" dominant-baseline=\"auto\" "
        << "font-family=\"Arial, Helvetica, sans-serif\" font-size=\"13\" "
        << "font-weight=\"bold\" fill=\"" << strokeClr << "\">No</text>\n";

    // 10. Done (rounded pill)
    svg << "  <rect x=\"" << (cx - 50) << "\" y=\"" << (yDone - 18)
        << "\" width=\"100\" height=\"36\" rx=\"18\" ry=\"18\" fill=\""
        << colDone << "\" stroke=\"" << strokeClr << "\" stroke-width=\"2\"/>\n";
    svg << "  <text x=\"" << cx << "\" y=\"" << yDone
        << "\" text-anchor=\"middle\" dominant-baseline=\"central\" "
        << "font-family=\"Arial, Helvetica, sans-serif\" font-size=\"15\" "
        << "font-weight=\"bold\" fill=\"#1e293b\">Done</text>\n";

    // Legend
    double legendX = 20;
    double legendY = H - 120;
    double swatchSize = 14;
    double legendLineH = 22;
    struct LegendItem { std::string color; std::string label; };
    LegendItem legend[] = {
        { colInput,  "Input / Output" },
        { colTree,   "Tree Construction" },
        { colL01,    "Level 0-1 (Heatmap MC)" },
        { colL2,     "Level 2 (Arc-Distance MC)" },
        { colL3,     "Level 3 (Smoothing MC)" },
    };
    int nLegend = 5;

    svg << "  <text x=\"" << legendX << "\" y=\"" << (legendY - 10)
        << "\" font-family=\"Arial, Helvetica, sans-serif\" font-size=\"13\" "
        << "font-weight=\"bold\" fill=\"#475569\">Legend</text>\n";

    for (int i = 0; i < nLegend; ++i) {
        double ly = legendY + i * legendLineH;
        svg << "  <rect x=\"" << legendX << "\" y=\"" << (ly - swatchSize / 2.0)
            << "\" width=\"" << swatchSize << "\" height=\"" << swatchSize
            << "\" rx=\"3\" ry=\"3\" fill=\"" << legend[i].color
            << "\" stroke=\"" << strokeClr << "\" stroke-width=\"1\"/>\n";
        svg << "  <text x=\"" << (legendX + swatchSize + 8) << "\" y=\"" << ly
            << "\" dominant-baseline=\"central\" "
            << "font-family=\"Arial, Helvetica, sans-serif\" font-size=\"12\" "
            << "fill=\"#334155\">" << legend[i].label << "</text>\n";
    }

    svg << "</svg>\n";

    // Write to file
    std::ofstream out(output_path);
    if (!out.is_open()) {
        return false;
    }
    out << svg.str();
    out.close();
    return out.good() || !out.fail();
}

#endif // FLOWCHART_GENERATOR_H
