"""Tests for the self-contained offline interactive SurfaceKnots Web Atlas (`atlas.py`)."""

from __future__ import annotations

import json
from pathlib import Path

from surfaceknots.atlas import _normalize_entry, build_atlas_payload, render_interactive_atlas_html
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.yoshikawa_data import YOSHIKAWA_BY_LABEL, compute_full_invariant_record


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "code" / "data"
EXPO_ATLAS = REPO_ROOT / "expo" / "surfaceknots_atlas.html"


def _add_r1_kink(d: MarkedDiagram) -> MarkedDiagram:
    opp = list(d.opp) + [0, 0, 0, 0]
    vtype = list(d.vtype) + ["X0"]
    nv = d.num_vertices
    d0 = 0
    d1 = opp[d0]
    opp[4 * nv + 0] = 4 * nv + 1
    opp[4 * nv + 1] = 4 * nv + 0
    opp[4 * nv + 2] = d1
    opp[d1] = 4 * nv + 2
    opp[4 * nv + 3] = d0
    opp[d0] = 4 * nv + 3
    return MarkedDiagram(num_vertices=nv + 1, opp=tuple(opp), vtype=tuple(vtype), free_loops=d.free_loops)


def test_atlas_payload_and_pairwise_separation() -> None:
    payload = build_atlas_payload(DATA_DIR)
    counts = payload["counts"]
    assert counts["yoshikawa_ch10"] == 23
    assert counts["expanded_ch11"] == 10
    ch12_path = DATA_DIR / "expanded_ch12.json"
    if ch12_path.exists():
        ch12_expected = len(json.loads(ch12_path.read_text())["entries"])
        assert counts["expanded_ch12"] == ch12_expected
        assert ch12_expected >= 25
        assert counts["spoke_total"] >= 37
        assert counts["total_distinct_classes"] >= 95
    else:
        assert counts["spoke_total"] >= 44
        assert counts["total_distinct_classes"] >= 77
    total_n = counts["total_distinct_classes"]

    entries = payload["entries"]
    assert len(entries) == total_n

    # All unified Hub + Spoke entries must have pairwise-distinct 10-layer SHA-256 hash vectors,
    # valid nearest-neighbor witnesses, clean pi_2 display, and explicit epistemic rigor badges!
    hash_tuples = [tuple(e["lh"]) for e in entries]
    assert len(set(hash_tuples)) == total_n
    for e in entries:
        assert e["nw"].startswith("Nearest: "), e["label"]
        assert "rk\u2265?" not in e["ld"][8] and "(None)" not in e["ld"][8], e["label"]
        assert e["rigor_badge"].startswith(("PROVABLE (", "PACKET (", "HEURISTIC (", "BOUNDED ([")), e["label"]

    html_str = render_interactive_atlas_html(DATA_DIR)
    assert "tailwindcss.min.js" not in html_str
    assert "<script src=" not in html_str
    assert "<style>" in html_str
    assert 'id="atlas-data"' in html_str
    assert EXPO_ATLAS.exists()
    assert len(EXPO_ATLAS.read_text()) == len(html_str)


def test_atlas_layer_hashes_o4_and_stabilization_invariance() -> None:
    # 10_1 vs 4-R1-stabilized 10_1 (where bounded-BFS constructive_so4_amphicheiral drops out)
    d10 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_1"].canonical_epd)
    d10_stab = d10
    for _ in range(4):
        d10_stab = _add_r1_kink(d10_stab)
    e10 = _normalize_entry(
        compute_full_invariant_record(d10, label="10_1"),
        tier_code="yoshikawa_ch10",
        tier_label="Yoshikawa (ch <= 10)",
        hub_reals=[],
    )
    e10_stab = _normalize_entry(
        compute_full_invariant_record(d10_stab, label="10_1_r1x4"),
        tier_code="yoshikawa_ch10",
        tier_label="Yoshikawa (ch <= 10)",
        hub_reals=[],
    )
    assert e10["lh"] == e10_stab["lh"]

    # O(4) mirror / orientation reversal invariance on chiral RP^2 (e=+2 vs e=-2) and 10_2
    for lbl in ("2_1^-1", "10_2", "8_1^-1,-1"):
        d = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL[lbl].canonical_epd)
        eb = _normalize_entry(
            compute_full_invariant_record(d, label=lbl),
            tier_code="yoshikawa_ch10",
            tier_label="Yoshikawa (ch <= 10)",
            hub_reals=[],
        )
        em = _normalize_entry(
            compute_full_invariant_record(d.mirror_3d(), label=f"{lbl}_mir"),
            tier_code="yoshikawa_ch10",
            tier_label="Yoshikawa (ch <= 10)",
            hub_reals=[],
        )
        ez = _normalize_entry(
            compute_full_invariant_record(d.reverse_orientation_zt(), label=f"{lbl}_rzt"),
            tier_code="yoshikawa_ch10",
            tier_label="Yoshikawa (ch <= 10)",
            hub_reals=[],
        )
        assert eb["lh"] == em["lh"] == ez["lh"], lbl
