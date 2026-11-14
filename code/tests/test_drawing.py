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
