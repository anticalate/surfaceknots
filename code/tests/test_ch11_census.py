"""Tests for the expanded ch = 11 surface-link census, dataset integrity, and composite/split operations."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from surfaceknots.ch11_census import disjoint_union, planar_connected_sum
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.moves import (
    find_planar_connected_sum_decomposition,
    has_certified_cyclic_unknotted_tube,
)
from surfaceknots.unlink import certify_classical_unlink
from surfaceknots.yoshikawa_data import YOSHIKAWA_BY_LABEL, invariant_discrimination_key


DATA_DIR = Path(__file__).resolve().parents[1] / "data"


def test_planar_connected_sum_and_disjoint_union() -> None:
    d_8_1 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["8_1"].canonical_epd)
    d_t2 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^1"].canonical_epd)

    csum = planar_connected_sum(d_8_1, 0, d_t2, 0, flip_splice=False)
    if not (
        certify_classical_unlink(csum.resolve(1)).is_unlink
        and certify_classical_unlink(csum.resolve(-1)).is_unlink
    ):
        csum = planar_connected_sum(d_8_1, 0, d_t2, 0, flip_splice=True)

    assert csum.ch_index == 10
    assert certify_classical_unlink(csum.resolve(1)).is_unlink
    assert certify_classical_unlink(csum.resolve(-1)).is_unlink
    assert find_planar_connected_sum_decomposition(csum) is not None

    d_hopf = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["6_1^0,1"].canonical_epd)
    sp = disjoint_union(d_hopf, d_t2)
    assert sp.ch_index == 8
    assert sp.surface_topology().num_components == 3


def test_expanded_ch11_dataset_integrity_and_distinctness() -> None:
    yosh_path = DATA_DIR / "yoshikawa_ch10.json"
    ch11_json = DATA_DIR / "expanded_ch11.json"
    ch11_csv = DATA_DIR / "expanded_ch11.csv"
    assert ch11_json.exists()
    assert ch11_csv.exists()

    yosh_entries = json.loads(yosh_path.read_text())
    yosh_keys = {invariant_discrimination_key(r) for r in yosh_entries}
    assert len(yosh_keys) == 23

    payload = json.loads(ch11_json.read_text())
    entries = payload["entries"]
    assert payload["ch_index"] == 11
    assert payload["certified_new_prime_invariant_classes"] == len(entries)
    assert len(entries) > 0

    with ch11_csv.open() as f:
        csv_rows = list(csv.DictReader(f))
    assert len(csv_rows) == len(entries)

    seen_ch11_keys = set()
    for rec in entries:
        assert rec["label"].startswith("11_")
        m_comps = rec["num_move_components"]
        assert m_comps == rec["num_diagram_representatives"]
        if m_comps == 1:
            assert rec["cell_status"] == "EXACT_SINGLETON"
            assert rec["move_component_ids"] == [rec["label"]]
        else:
            assert rec["cell_status"] == "NEW_PACKET"
            assert len(rec["move_component_ids"]) == m_comps
            assert rec["move_component_ids"][0] == f"{rec['label']}.a"

        wit = rec["nearest_separating_witness"]
        assert wit["nearest_label"] != rec["label"]
        assert 1 <= wit["separating_layer_index"] <= 10

        d = MarkedDiagram.from_epd(rec["epd"])
        assert d.ch_index == 11
        assert certify_classical_unlink(d.resolve(1)).is_unlink
        assert certify_classical_unlink(d.resolve(-1)).is_unlink
        assert find_planar_connected_sum_decomposition(d) is None
        assert not has_certified_cyclic_unknotted_tube(d)

        dkey = invariant_discrimination_key(rec)
        assert dkey not in yosh_keys, f"{rec['label']} collided with a ch <= 10 Yoshikawa entry!"
        assert dkey not in seen_ch11_keys, f"Duplicate invariant key among certified ch=11 classes on {rec['label']}!"
        seen_ch11_keys.add(dkey)

    assert len(entries) == 10
    assert payload["exact_singletons_count"] == 9
    assert payload["new_packets_count"] == 1
    assert payload["exact_singletons_count"] + payload["new_packets_count"] == len(entries)

    by_lbl = {r["label"]: r for r in entries}
    assert "11_1^0,0,-2" in by_lbl
    assert by_lbl["11_1^0,0,-2"]["sym_quandle_q6_col"] == 96
    assert by_lbl["11_2^0,-2"]["sigma2_b1"] == 1
    assert by_lbl["11_3^0,-2"]["sigma2_b1"] == 0
    assert by_lbl["11_3^0,-2"]["sigma2_h1_torsion"] == [4]

    for rc in payload["residual_collisions"]:
        assert rc["cell_status"] == "AMBIGUOUS_SHADOW"
        assert len(rc["collision_matches"]) >= 1
        d_rc = MarkedDiagram.from_epd(rc["epd"])
        assert d_rc.ch_index == 11
        assert certify_classical_unlink(d_rc.resolve(1)).is_unlink
        assert certify_classical_unlink(d_rc.resolve(-1)).is_unlink

