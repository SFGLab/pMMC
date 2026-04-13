#!/usr/bin/env python3
"""
Generate a manuscript-quality SVG flowchart for the pMMC algorithm.

Features beyond the C++ FlowchartGenerator:
  - Dual-column layout (main flow + GPU/CPU execution detail)
  - CPU/GPU decision diamonds at each MC phase
  - Persistent GPU Cache sidebar
  - GPU OOM fallback arrows
  - MDS initialization optional path
  - Energy function annotations
  - Color-coded hierarchy levels

Usage:
    python generate_flowchart.py [-o OUTPUT_PATH]

Default output: doc/pMMC_algorithm_flowchart.svg
No external dependencies — pure SVG string generation.
"""

import argparse
import os

# ============================================================================
# SVG Helper Functions
# ============================================================================

def arrow_marker(marker_id, color):
    return (
        f'    <marker id="{marker_id}" markerWidth="10" markerHeight="7" '
        f'refX="10" refY="3.5" orient="auto" markerUnits="strokeWidth">\n'
        f'      <polygon points="0 0, 10 3.5, 0 7" fill="{color}"/>\n'
        f'    </marker>\n'
    )

def rounded_rect(x, y, w, h, rx, fill, stroke, sw=2):
    return (
        f'  <rect x="{x}" y="{y}" width="{w}" height="{h}" '
        f'rx="{rx}" ry="{rx}" fill="{fill}" stroke="{stroke}" '
        f'stroke-width="{sw}"/>\n'
    )

def text_elem(x, y, content, size=14, weight="normal", fill="#222",
              anchor="middle", baseline="central", font_style="normal"):
    return (
        f'  <text x="{x}" y="{y}" text-anchor="{anchor}" '
        f'dominant-baseline="{baseline}" '
        f'font-family="Arial, Helvetica, sans-serif" font-size="{size}" '
        f'font-weight="{weight}" font-style="{font_style}" '
        f'fill="{fill}">{content}</text>\n'
    )

def process_box(cx, cy, w, h, fill, stroke, lines, font_sizes=None):
    """Rounded rectangle with centered multi-line text."""
    s = rounded_rect(cx - w/2, cy - h/2, w, h, 12, fill, stroke)
    line_h = 20.0
    start_y = cy - ((len(lines) - 1) * line_h) / 2.0
    for i, line in enumerate(lines):
        fs = font_sizes[i] if font_sizes else (14 if i > 0 else 14)
        fw = "bold" if i == 0 else "normal"
        s += text_elem(cx, start_y + i * line_h, line, size=fs, weight=fw)
    return s

def diamond(cx, cy, w, h, fill, stroke, lines):
    """Diamond decision shape with centered text."""
    pts = (f"{cx},{cy - h/2} {cx + w/2},{cy} "
           f"{cx},{cy + h/2} {cx - w/2},{cy}")
    s = (f'  <polygon points="{pts}" fill="{fill}" '
         f'stroke="{stroke}" stroke-width="2"/>\n')
    line_h = 17.0
    start_y = cy - ((len(lines) - 1) * line_h) / 2.0
    for i, line in enumerate(lines):
        s += text_elem(cx, start_y + i * line_h, line, size=13, weight="bold")
    return s

def v_arrow(x, y1, y2, marker_id, color, sw=2):
    """Vertical arrow."""
    return (
        f'  <line x1="{x}" y1="{y1}" x2="{x}" y2="{y2}" '
        f'stroke="{color}" stroke-width="{sw}" '
        f'marker-end="url(#{marker_id})"/>\n'
    )

def h_arrow(x1, y, x2, marker_id, color, sw=2):
    """Horizontal arrow."""
    return (
        f'  <line x1="{x1}" y1="{y}" x2="{x2}" y2="{y}" '
        f'stroke="{color}" stroke-width="{sw}" '
        f'marker-end="url(#{marker_id})"/>\n'
    )

def path_arrow(points, marker_id, color, sw=2, dashed=False):
    """Polyline arrow through a list of (x, y) points."""
    d = "M " + " L ".join(f"{x} {y}" for x, y in points)
    dash = ' stroke-dasharray="6,4"' if dashed else ""
    return (
        f'  <path d="{d}" fill="none" stroke="{color}" '
        f'stroke-width="{sw}" marker-end="url(#{marker_id})"{dash}/>\n'
    )

def loopback_arrow(x_center, y_start, y_end, x_offset, marker_id, color):
    """Right-angle loopback arrow."""
    xr = x_center + x_offset
    return (
        f'  <path d="M {x_center} {y_start} L {xr} {y_start} '
        f'L {xr} {y_end} L {x_center} {y_end}" fill="none" '
        f'stroke="{color}" stroke-width="2" '
        f'marker-end="url(#{marker_id})"/>\n'
    )

def side_box(cx, cy, w, h, fill, stroke, lines, font_sizes=None,
             rx=8, stroke_width=1.5):
    """Smaller annotation box for GPU/CPU detail column."""
    s = rounded_rect(cx - w/2, cy - h/2, w, h, rx, fill, stroke,
                     sw=stroke_width)
    line_h = 16.0
    start_y = cy - ((len(lines) - 1) * line_h) / 2.0
    for i, line in enumerate(lines):
        fs = font_sizes[i] if font_sizes else (12 if i > 0 else 12)
        fw = "bold" if i == 0 else "normal"
        s += text_elem(cx, start_y + i * line_h, line, size=fs, weight=fw,
                       fill="#333")
    return s


# ============================================================================
# Colors
# ============================================================================
COL_INPUT  = "#dbeafe"   # light blue
COL_TREE   = "#e0e7ff"   # light indigo
COL_L01    = "#fef3c7"   # amber
COL_L2     = "#d1fae5"   # green
COL_L3     = "#ede9fe"   # violet
COL_OUTPUT = "#ffedd5"   # orange
COL_DECIS  = "#f3f4f6"   # gray
COL_DONE   = "#d1d5db"   # dark gray
COL_GPU    = "#bbf7d0"   # light green (GPU)
COL_CPU    = "#bfdbfe"   # light blue  (CPU)
COL_CACHE  = "#fce7f3"   # pink (GPU cache)
COL_SPHERE = "#fce7f3"   # pink
COL_MDS    = "#e0f2fe"   # sky blue
STROKE     = "#475569"   # slate
STROKE_GPU = "#16a34a"   # green
STROKE_CPU = "#2563eb"   # blue


# ============================================================================
# Main Flowchart Generation
# ============================================================================

def generate_flowchart():
    # Layout constants
    W = 1100.0          # total canvas width
    cx_main = 320.0     # center of main (left) column
    cx_detail = 830.0   # center of detail (right) column
    box_w = 340.0       # main box width
    box_h = 68.0        # main box height
    dia_w = 240.0       # diamond width
    dia_h = 76.0        # diamond height
    side_w = 280.0      # side box width
    side_h = 56.0       # side box height
    gap = 30.0          # vertical gap

    # Y positions (computed top-down)
    y = 50.0
    y_title = y

    y += 55.0
    y_input = y         # Input Data box

    y += box_h/2 + gap + box_h/2
    y_tree = y          # Build Tree

    y += box_h/2 + gap + dia_h/2
    y_ens_start = y     # Ensemble diamond

    y += dia_h/2 + gap + dia_h/2
    y_mds = y           # MDS init diamond

    y += dia_h/2 + gap + box_h/2
    y_l01 = y           # Level 0-1 Heatmap MC

    y += box_h/2 + gap + dia_h/2
    y_sphere = y        # Spherical Container diamond

    y += dia_h/2 + gap + box_h/2
    y_l2 = y            # Level 2 Arc-Distance MC

    # Extra gap between L2 and L3 to fit GPU cache in the right column
    y += box_h/2 + 120 + box_h/2
    y_l3 = y            # Level 3 Smoothing MC

    y += box_h/2 + gap + box_h/2
    y_out = y           # Output Generation

    y += box_h/2 + gap + dia_h/2
    y_ens_end = y       # More models? diamond

    y += dia_h/2 + gap + 20
    y_done = y          # Done pill

    H = y_done + 55.0   # canvas height

    # GPU/CPU detail y-positions (aligned with main boxes)
    y_gpu_l01 = y_l01
    y_gpu_l2  = y_l2
    y_gpu_l3  = y_l3
    # GPU cache placed between L2 and L3 detail sections
    y_cache   = y_l2 + box_h/2 + 65

    # ---- Build SVG ----
    svg = []
    svg.append('<?xml version="1.0" encoding="UTF-8"?>\n')
    svg.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{W}" height="{H}" viewBox="0 0 {W} {H}">\n'
    )
    # Background
    svg.append(f'  <rect width="{W}" height="{H}" fill="#ffffff"/>\n')

    # Defs: arrowheads
    svg.append('  <defs>\n')
    svg.append(arrow_marker("arrowGray", STROKE))
    svg.append(arrow_marker("arrowGreen", STROKE_GPU))
    svg.append(arrow_marker("arrowBlue", STROKE_CPU))
    svg.append(arrow_marker("arrowOrange", "#ea580c"))
    svg.append('  </defs>\n')

    # ---- Title ----
    svg.append(text_elem(W/2, y_title, "pMMC Algorithm Flowchart",
                         size=24, weight="bold", fill="#1e293b"))
    svg.append(text_elem(W/2, y_title + 26,
                         "Parallel Multiscale Monte Carlo for 3D Genome Reconstruction",
                         size=13, fill="#64748b"))

    # ==================================================================
    # LEFT COLUMN: Main Algorithm Flow
    # ==================================================================

    # 1. Input Data Loading
    svg.append(process_box(cx_main, y_input, box_w, box_h, COL_INPUT, STROKE, [
        "Input Data Loading",
        "BED anchors, BEDPE clusters,",
        "singletons (or Hi-C contact matrix)"
    ]))

    # Arrow: input -> tree
    svg.append(v_arrow(cx_main, y_input + box_h/2, y_tree - box_h/2,
                       "arrowGray", STROKE))

    # 2. Build Hierarchical Tree
    svg.append(process_box(cx_main, y_tree, box_w, box_h, COL_TREE, STROKE, [
        "Build Hierarchical Tree",
        "4 levels: Chr &gt; Segment &gt; Anchor &gt; Subanchor"
    ]))

    # Arrow: tree -> ensemble
    svg.append(v_arrow(cx_main, y_tree + box_h/2, y_ens_start - dia_h/2,
                       "arrowGray", STROKE))

    # 3. Ensemble Loop Start
    svg.append(diamond(cx_main, y_ens_start, dia_w, dia_h, COL_DECIS, STROKE, [
        "Ensemble Loop", "m = 1..M"
    ]))

    # Arrow: ensemble -> MDS
    svg.append(v_arrow(cx_main, y_ens_start + dia_h/2, y_mds - dia_h/2,
                       "arrowGray", STROKE))

    # 4. MDS Initialization Diamond
    svg.append(diamond(cx_main, y_mds, dia_w, dia_h, COL_MDS, "#0284c7", [
        "MDS Init?"
    ]))
    # Side label
    svg.append(text_elem(cx_main + dia_w/2 + 10, y_mds - 10,
                         "Yes: Classical MDS from",
                         size=11, weight="normal", fill="#0369a1",
                         anchor="start", font_style="italic"))
    svg.append(text_elem(cx_main + dia_w/2 + 10, y_mds + 8,
                         "distance heatmap (Torgerson)",
                         size=11, weight="normal", fill="#0369a1",
                         anchor="start", font_style="italic"))

    # Arrow: MDS -> L0-1
    svg.append(v_arrow(cx_main, y_mds + dia_h/2, y_l01 - box_h/2,
                       "arrowGray", STROKE))

    # 5. Level 0-1: Heatmap Monte Carlo
    svg.append(process_box(cx_main, y_l01, box_w, box_h, COL_L01, STROKE, [
        "Level 0\u20131: Heatmap Monte Carlo",
        "Chromosome &amp; Segment scale",
        "E = \u03a3((d\u1d62\u2c7c \u2212 D\u1d62\u2c7c)/D\u1d62\u2c7c)\u00b2 + E_boundary"
    ]))

    # Arrow: L0-1 -> sphere
    svg.append(v_arrow(cx_main, y_l01 + box_h/2, y_sphere - dia_h/2,
                       "arrowGray", STROKE))

    # 6. Spherical Container Diamond
    svg.append(diamond(cx_main, y_sphere, dia_w, dia_h, COL_SPHERE,
                       "#db2777", ["Spherical", "Container?"]))
    svg.append(text_elem(cx_main + dia_w/2 + 10, y_sphere,
                         "If enabled: quadratic boundary penalty",
                         size=11, fill="#9f1239", anchor="start",
                         font_style="italic"))

    # Arrow: sphere -> L2
    svg.append(v_arrow(cx_main, y_sphere + dia_h/2, y_l2 - box_h/2,
                       "arrowGray", STROKE))

    # 7. Level 2: Arc-Distance MC
    svg.append(process_box(cx_main, y_l2, box_w, box_h, COL_L2, STROKE, [
        "Level 2: Arc-Distance Monte Carlo",
        "Anchor scale, PET cluster arcs",
        "E = E_arc + E_loop + E_CTCF"
    ]))

    # Arrow: L2 -> L3
    svg.append(v_arrow(cx_main, y_l2 + box_h/2, y_l3 - box_h/2,
                       "arrowGray", STROKE))

    # 8. Level 3: Smoothing MC
    svg.append(process_box(cx_main, y_l3, box_w, box_h, COL_L3, STROKE, [
        "Level 3: Smoothing Monte Carlo",
        "Subanchor scale, polymer springs",
        "E = E_linker + E_angle"
    ]))

    # Arrow: L3 -> output
    svg.append(v_arrow(cx_main, y_l3 + box_h/2, y_out - box_h/2,
                       "arrowGray", STROKE))

    # 9. Output Generation
    svg.append(process_box(cx_main, y_out, box_w, box_h, COL_OUTPUT, STROKE, [
        "Output Generation",
        "HCM, PDB, mmCIF formats",
        "Visualization scripts (.pml, .cxc)"
    ]))

    # Arrow: output -> ensemble end
    svg.append(v_arrow(cx_main, y_out + box_h/2, y_ens_end - dia_h/2,
                       "arrowGray", STROKE))

    # 10. Ensemble Loop End
    svg.append(diamond(cx_main, y_ens_end, dia_w, dia_h, COL_DECIS, STROKE, [
        "More models?"
    ]))

    # "Yes" loopback arrow
    svg.append(loopback_arrow(cx_main + dia_w/2, y_ens_end,
                              y_mds - dia_h/2 + 5, 75, "arrowGreen",
                              STROKE_GPU))
    svg.append(text_elem(cx_main + dia_w/2 + 8, y_ens_end - 8,
                         "Yes", size=13, weight="bold", fill=STROKE_GPU,
                         anchor="start"))

    # "No" arrow to Done
    svg.append(v_arrow(cx_main, y_ens_end + dia_h/2, y_done - 18,
                       "arrowGray", STROKE))
    svg.append(text_elem(cx_main + 10, y_ens_end + dia_h/2 + 16,
                         "No", size=13, weight="bold", fill=STROKE,
                         anchor="start"))

    # Done pill
    svg.append(rounded_rect(cx_main - 50, y_done - 18, 100, 36, 18,
                            COL_DONE, STROKE))
    svg.append(text_elem(cx_main, y_done, "Done", size=15, weight="bold",
                         fill="#1e293b"))

    # ==================================================================
    # RIGHT COLUMN: GPU/CPU Execution Detail
    # ==================================================================

    # Column header
    svg.append(text_elem(cx_detail, y_title + 4,
                         "GPU / CPU Execution Detail",
                         size=16, weight="bold", fill="#334155"))

    # ---- L0-1 GPU/CPU detail ----
    # Horizontal arrow from L0-1 box to detail
    svg.append(h_arrow(cx_main + box_w/2, y_gpu_l01,
                       cx_detail - side_w/2 - 50, "arrowGray", STROKE))

    # GPU box
    y_gpu_top = y_gpu_l01 - side_h - 8
    svg.append(side_box(cx_detail, y_gpu_top, side_w, side_h,
                        COL_GPU, STROKE_GPU, [
        "GPU: Multi-Warp Parallel SA",
        "N warps \u00d7 independent SA trajectories",
        "half3 positions (FP16 memory savings)",
    ]))

    # CPU box
    y_cpu_bot = y_gpu_l01 + side_h + 8
    svg.append(side_box(cx_detail, y_cpu_bot, side_w, side_h,
                        COL_CPU, STROKE_CPU, [
        "CPU: OpenMP Metropolis SA",
        "Single-thread or parallel regions",
        "Full float32 precision",
    ]))

    # Decision mini-diamond between GPU/CPU
    mini_dia_w = 120
    mini_dia_h = 44
    svg.append(diamond(cx_detail - side_w/2 - 25, y_gpu_l01,
                       mini_dia_w, mini_dia_h, COL_DECIS, STROKE, [
        "CUDA?"
    ]))

    # Arrows from mini-diamond to GPU and CPU
    svg.append(path_arrow([
        (cx_detail - side_w/2 - 25, y_gpu_l01 - mini_dia_h/2),
        (cx_detail - side_w/2 - 25, y_gpu_top),
        (cx_detail - side_w/2, y_gpu_top)
    ], "arrowGreen", STROKE_GPU))
    svg.append(text_elem(cx_detail - side_w/2 - 45, y_gpu_l01 - mini_dia_h/2 - 5,
                         "Yes", size=11, weight="bold", fill=STROKE_GPU))

    svg.append(path_arrow([
        (cx_detail - side_w/2 - 25, y_gpu_l01 + mini_dia_h/2),
        (cx_detail - side_w/2 - 25, y_cpu_bot),
        (cx_detail - side_w/2, y_cpu_bot)
    ], "arrowBlue", STROKE_CPU))
    svg.append(text_elem(cx_detail - side_w/2 - 45, y_gpu_l01 + mini_dia_h/2 + 14,
                         "No", size=11, weight="bold", fill=STROKE_CPU))

    # ---- L2 GPU/CPU detail ----
    svg.append(h_arrow(cx_main + box_w/2, y_gpu_l2,
                       cx_detail - side_w/2 - 10, "arrowGray", STROKE))

    svg.append(side_box(cx_detail, y_gpu_l2 - 34, side_w, side_h,
                        COL_GPU, STROKE_GPU, [
        "GPU: Warp-Level Arc Scoring",
        "N-body energy via warp reduction",
        "Best-of-warps selection on host",
    ]))

    svg.append(side_box(cx_detail, y_gpu_l2 + 34, side_w, 40,
                        COL_CPU, STROKE_CPU, [
        "CPU Fallback: Local O(N) scoring",
        "Incremental energy updates",
    ]))

    # OOM fallback dashed arrow
    svg.append(path_arrow([
        (cx_detail + side_w/2, y_gpu_l2 - 34 + side_h/2),
        (cx_detail + side_w/2 + 20, y_gpu_l2 - 34 + side_h/2),
        (cx_detail + side_w/2 + 20, y_gpu_l2 + 34),
        (cx_detail + side_w/2, y_gpu_l2 + 34)
    ], "arrowOrange", "#ea580c", dashed=True))
    svg.append(text_elem(cx_detail + side_w/2 + 25, y_gpu_l2,
                         "OOM", size=10, weight="bold", fill="#ea580c",
                         anchor="start"))
    svg.append(text_elem(cx_detail + side_w/2 + 25, y_gpu_l2 + 13,
                         "fallback", size=10, fill="#ea580c",
                         anchor="start"))

    # ---- Persistent GPU Cache (between L2 and L3) ----
    cache_h = 55
    svg.append(rounded_rect(cx_detail - side_w/2, y_cache - cache_h/2,
                            side_w, cache_h, 8, COL_CACHE, "#db2777", 2))
    svg.append(text_elem(cx_detail, y_cache - 14,
                         "Persistent GPU Resource Cache",
                         size=12, weight="bold", fill="#831843"))
    svg.append(text_elem(cx_detail, y_cache + 4,
                         "curandState, device buffers, isDone flags",
                         size=11, fill="#9f1239"))
    svg.append(text_elem(cx_detail, y_cache + 20,
                         "Eliminates per-block cudaMalloc overhead",
                         size=10, fill="#9f1239", font_style="italic"))

    # ---- L3 GPU/CPU detail ----
    svg.append(h_arrow(cx_main + box_w/2, y_gpu_l3,
                       cx_detail - side_w/2 - 10, "arrowGray", STROKE))

    svg.append(side_box(cx_detail, y_gpu_l3 - 28, side_w, side_h,
                        COL_GPU, STROKE_GPU, [
        "GPU: Multi-Warp Smooth Kernel",
        "Per-warp independent SA copies",
        "Spring length + angular bending",
    ]))

    svg.append(side_box(cx_detail, y_gpu_l3 + 34, side_w, 36,
                        COL_CPU, STROKE_CPU, [
        "CPU Fallback: Sequential smooth",
    ]))

    # ==================================================================
    # LEGEND
    # ==================================================================
    legend_x = 20
    legend_y = H - 140
    swatch = 14
    line_h = 22

    items = [
        (COL_INPUT,  STROKE,     "Input / Output"),
        (COL_TREE,   STROKE,     "Tree Construction"),
        (COL_L01,    STROKE,     "Level 0\u20131 (Heatmap MC)"),
        (COL_L2,     STROKE,     "Level 2 (Arc-Distance MC)"),
        (COL_L3,     STROKE,     "Level 3 (Smoothing MC)"),
        (COL_GPU,    STROKE_GPU, "GPU Execution Path"),
        (COL_CPU,    STROKE_CPU, "CPU Execution Path"),
        (COL_CACHE,  "#db2777",  "Persistent GPU Cache"),
    ]

    svg.append(text_elem(legend_x + 5, legend_y - 12, "Legend",
                         size=13, weight="bold", fill=STROKE, anchor="start"))

    for i, (color, stroke_c, label) in enumerate(items):
        ly = legend_y + i * line_h
        svg.append(rounded_rect(legend_x, ly - swatch/2, swatch, swatch,
                                3, color, stroke_c, 1))
        svg.append(text_elem(legend_x + swatch + 8, ly, label,
                             size=12, fill="#334155", anchor="start"))

    # Close SVG
    svg.append('</svg>\n')

    return ''.join(svg)


# ============================================================================
# Entry Point
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate pMMC algorithm flowchart (SVG)")
    parser.add_argument("-o", "--output",
                        default=None,
                        help="Output SVG path (default: doc/pMMC_algorithm_flowchart.svg)")
    args = parser.parse_args()

    # Default output path: doc/ relative to script location
    if args.output is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(script_dir)
        output_path = os.path.join(project_root, "doc",
                                   "pMMC_algorithm_flowchart.svg")
    else:
        output_path = args.output

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    svg_content = generate_flowchart()

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(svg_content)

    print(f"Flowchart written to {output_path}")


if __name__ == "__main__":
    main()
