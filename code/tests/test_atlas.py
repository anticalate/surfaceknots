"""Tests for the self-contained 74-entry interactive SurfaceKnots Web Atlas (`atlas.py`)."""

from __future__ import annotations

from pathlib import Path

from surfaceknots.atlas import build_atlas_payload, render_interactive_atlas_html


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "code" / "data"
EXPO_ATLAS = REPO_ROOT / "expo" / "surfaceknots_atlas.html"


def test_atlas_payload_and_pairwise_separation() -> None:
    payload = build_atlas_payload(DATA_DIR)
    counts = payload["counts"]
    assert counts["yoshikawa_ch10"] == 23
    assert counts["expanded_ch11"] == 10
    assert counts["exact_ch12"] == 5
    assert counts["spoke_total"] == 41
    assert counts["total_distinct_classes"] == 74
    assert counts["hub_realized_count"] == 15

    entries = payload["entries"]
    assert len(entries) == 74

    # All 74 unified Hub + Spoke entries must have pairwise-distinct 10-layer SHA-256 hash vectors!
    hash_tuples = [tuple(e["lh"]) for e in entries]
    assert len(set(hash_tuples)) == 74

    html_str = render_interactive_atlas_html(DATA_DIR)
    assert "tailwindcss.min.js" in html_str
    assert 'id="atlas-data"' in html_str
    assert EXPO_ATLAS.exists()
    assert len(EXPO_ATLAS.read_text()) == len(html_str)
