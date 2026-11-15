"""Tests for automated planar diagram layout, SVG/TikZ rendering, and the `/pics` gallery."""

from __future__ import annotations

import json
import math
from pathlib import Path
import shutil
import subprocess
import xml.etree.ElementTree as ET

from surfaceknots.ch11_census import disjoint_union
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.drawing import (
    compute_diagram_layout,
    label_to_stem,
    render_svg,
    render_tikz,
)
from surfaceknots.yoshikawa_data import YOSHIKAWA_BY_LABEL


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "code" / "data"
PICS_DIR = REPO_ROOT / "pics"


def test_layout_geometry_and_split_free_loops() -> None:
    d0 = MarkedDiagram.unknot_sphere()
    lay0 = compute_diagram_layout(d0)
    assert len(lay0.vertex_positions) == 0
    assert len(lay0.circles) == 1
    assert lay0.num_surface_components == 1

    d_hopf = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["6_1^0,1"].canonical_epd)
    d_split = disjoint_union(d_hopf, d0)
    lay_sp = compute_diagram_layout(d_split)
    assert len(lay_sp.vertex_positions) == 6
    assert len(lay_sp.markers) == 2
    assert len(lay_sp.circles) == 1
    assert lay_sp.num_surface_components == 3
    for x, y in lay_sp.vertex_positions:
        assert math.isfinite(x) and math.isfinite(y)

    # Distinct vertex coordinates
    for i in range(len(lay_sp.vertex_positions)):
        for j in range(i + 1, len(lay_sp.vertex_positions)):
            dist = math.hypot(
                lay_sp.vertex_positions[i][0] - lay_sp.vertex_positions[j][0],
                lay_sp.vertex_positions[i][1] - lay_sp.vertex_positions[j][1],
            )
            assert dist > 1e-3

    svg_txt = render_svg(d_split, title="6_1^0,1 u 0_1", subtitle="split test")
    root = ET.fromstring(svg_txt)
    assert root.tag.endswith("svg")


def test_canonical_gallery_files_exist_and_parse(tmp_path: Path) -> None:
    yosh = json.loads((DATA_DIR / "yoshikawa_ch10.json").read_text())
    ch11 = json.loads((DATA_DIR / "expanded_ch11.json").read_text())["entries"]

    expected_labels = [r["label"] for r in yosh] + [r["label"] for r in ch11]
    for r in ch11:
        if r["num_move_components"] > 1:
            expected_labels.extend(r["move_component_ids"])

    stems = [label_to_stem(lbl) for lbl in expected_labels]
    assert len(set(stems)) == len(stems), "Collision in label_to_stem!"
    assert len(stems) == 36

    for lbl, stem in zip(expected_labels, stems):
        svg_file = PICS_DIR / "svg" / f"{stem}.svg"
        tex_file = PICS_DIR / "tikz" / f"{stem}.tex"
        assert svg_file.exists(), f"Missing SVG for {lbl}: {svg_file}"
        assert tex_file.exists(), f"Missing TikZ for {lbl}: {tex_file}"

        root = ET.fromstring(svg_file.read_text())
        assert root.tag.endswith("svg")
        tex_txt = tex_file.read_text()
        assert r"\begin{tikzpicture}" in tex_txt
        assert r"\end{tikzpicture}" in tex_txt

    # Live pdflatex compilation check on a multi-component marked diagram (`10_1^0,0,1`)
    if shutil.which("pdflatex") is not None:
        sample_tex = render_tikz(
            MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_1^0,0,1"].canonical_epd),
            label="10_1^0,0,1",
            standalone=True,
        )
        tex_p = tmp_path / "sample.tex"
        tex_p.write_text(sample_tex)
        res = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", f"-output-directory={tmp_path}", str(tex_p)],
            capture_output=True,
            text=True,
            check=False,
        )
        assert res.returncode == 0, f"pdflatex failed:\n{res.stdout}\n{res.stderr}"
        assert (tmp_path / "sample.pdf").exists()


def test_splitting_movie_and_extended_gallery() -> None:
    from surfaceknots.drawing import compute_splitting_layouts, render_splitting_movie_svg

    d = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_1^0,1"].canonical_epd)
    lay_m, lay_0, lay_p = compute_splitting_layouts(d)

    # Pixel-registered: identical vertex coordinates and bounding box across all 3 slices
    assert lay_m.vertex_positions == lay_0.vertex_positions == lay_p.vertex_positions
    assert lay_m.bbox == lay_0.bbox == lay_p.bbox
    # Each marked vertex adds 2 hyperbolic bridge Beziers in L_- and L_+
    assert len(lay_m.beziers) == len(lay_0.beziers) + 2 * d.num_marked
    assert len(lay_p.beziers) == len(lay_0.beziers) + 2 * d.num_marked

    movie_svg = render_splitting_movie_svg(d, title="10_1^0,1", subtitle="movie test")
    root = ET.fromstring(movie_svg)
    assert root.tag.endswith("svg")
    assert "Hyperbolic splitting cobordism" in movie_svg

    # Verify exported extended gallery files in pics/families and pics/movies
    assert (PICS_DIR / "families" / "svg" / "Tube_4_1.svg").exists()
    assert (PICS_DIR / "families" / "tikz" / "Tube_4_1.tex").exists()
    assert (PICS_DIR / "movies" / "11_1_sup_0_0_m2_movie.svg").exists()
    assert (PICS_DIR / "movies" / "Mixed_L6a4_01_movie.svg").exists()


def test_one_vertex_self_loop_and_cut_vertex_layouts() -> None:
    from surfaceknots.drawing import compute_splitting_layouts, render_splitting_movie_svg

    # 1. Single-vertex planar figure-eight / double-loop diagrams (m == 1):
    for epd in ("X[1,1,2,2]", "Y[1,1,2,2]", "Z[1,2,2,1]"):
        d1 = MarkedDiagram.from_epd(epd)
        lay1 = compute_diagram_layout(d1)
        assert len(lay1.vertex_positions) == 1
        assert math.isfinite(lay1.vertex_positions[0][0]) and math.isfinite(lay1.vertex_positions[0][1])
        assert len(lay1.beziers) == 4
        svg1 = render_svg(d1, title=epd)
        assert ET.fromstring(svg1).tag.endswith("svg")
        tex1 = render_tikz(d1, label=epd)
        assert r"\begin{tikzpicture}" in tex1 and r"\end{tikzpicture}" in tex1
        m_lay, z_lay, p_lay = compute_splitting_layouts(d1)
        assert m_lay.vertex_positions == z_lay.vertex_positions == p_lay.vertex_positions
        assert ET.fromstring(render_splitting_movie_svg(d1, title=epd)).tag.endswith("svg")

    # 2. Planar cut-vertex blocks & R1 self-loop dumbbells:
    cut_cases = [
        "X[1,2,3,4], Y[2,5,5,1], Y[4,6,6,3]",
        "X[1,1,2,3], X[3,2,4,4]",
        "X[1,2,3,4], X[2,5,6,1], X[6,7,7,5], X[4,8,9,3], X[9,10,10,8]",
    ]
    for epd in cut_cases:
        dc = MarkedDiagram.from_epd(epd)
        lay_c = compute_diagram_layout(dc)
        assert len(lay_c.vertex_positions) == dc.num_vertices
        for i in range(len(lay_c.vertex_positions)):
            x1, y1 = lay_c.vertex_positions[i]
            assert math.isfinite(x1) and math.isfinite(y1)
            for j in range(i + 1, len(lay_c.vertex_positions)):
                x2, y2 = lay_c.vertex_positions[j]
                assert math.hypot(x1 - x2, y1 - y2) > 1e-3, (epd, i, j, lay_c.vertex_positions)
        assert ET.fromstring(render_svg(dc, title="cut_test")).tag.endswith("svg")
        assert r"\begin{tikzpicture}" in render_tikz(dc, label="cut_test")

    # 3. Iterated 4-R1-stabilized 10_1 (nested self-loops on 14 vertices):
    d10 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_1"].canonical_epd)
    for _ in range(4):
        opp = list(d10.opp) + [0, 0, 0, 0]
        vtype = list(d10.vtype) + ["X0"]
        nv = d10.num_vertices
        d0 = 0
        d1 = opp[d0]
        opp[4 * nv + 0] = 4 * nv + 1
        opp[4 * nv + 1] = 4 * nv + 0
        opp[4 * nv + 2] = d1
        opp[d1] = 4 * nv + 2
        opp[4 * nv + 3] = d0
        opp[d0] = 4 * nv + 3
        d10 = MarkedDiagram(num_vertices=nv + 1, opp=tuple(opp), vtype=tuple(vtype), free_loops=d10.free_loops)
    lay_r1 = compute_diagram_layout(d10)
    assert len(lay_r1.vertex_positions) == 14
    assert ET.fromstring(render_svg(d10, title="10_1_r1x4")).tag.endswith("svg")

