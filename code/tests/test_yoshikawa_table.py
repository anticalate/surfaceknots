"""Tests for the machine-verified Yoshikawa (ch <= 10) table, Table II errata, and distinguishing invariants."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.invariants import (
    compute_alexander_invariants,
    double_linking_mod2,
    mochizuki_3cocycle_invariant,
    peripheral_sl2f3_profile,
    sublink_summary,
    wirtinger_presentation,
)
from surfaceknots.unlink import certify_classical_unlink
from surfaceknots.yoshikawa_data import (
    YOSHIKAWA_BY_LABEL,
    YOSHIKAWA_TABLE,
    check_e1_ideal_against_reference,
    invariant_discrimination_key,
)


DATA_DIR = Path(__file__).resolve().parents[1] / "data"


def test_yoshikawa_datasets_complete_and_collision_free() -> None:
    json_path = DATA_DIR / "yoshikawa_ch10.json"
    csv_path = DATA_DIR / "yoshikawa_ch10.csv"
    assert json_path.exists()
    assert csv_path.exists()

    payload = json.loads(json_path.read_text())
    entries = payload["entries"] if isinstance(payload, dict) else payload
    assert len(entries) == 23
    assert len(YOSHIKAWA_TABLE) == 23

    with csv_path.open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 23

    seen_keys = set()
    errata_labels = []
    chiral_labels = []

    by_label = {}
    for rec in entries:
        by_label[rec["label"]] = rec
        assert rec["cell_status"] == "EXACT_SINGLETON"
        assert rec["num_move_components"] == 1
        assert rec["move_component_ids"] == [rec["label"]]
        wit = rec["nearest_separating_witness"]
        assert wit["nearest_label"] != rec["label"]
        assert 1 <= wit["separating_layer_index"] <= 10
        if rec["has_yoshikawa_table_ii_erratum"]:
            errata_labels.append(rec["label"])
        if rec["certified_chiral_by_c3"]:
            chiral_labels.append(rec["label"])
        key = invariant_discrimination_key(rec)
        assert key not in seen_keys, f"Invariant collision on {rec['label']}"
        seen_keys.add(key)

    assert len(seen_keys) == 23
    assert errata_labels == ["9_1^1,-2", "10_1^0,1", "10_2^0,-2"]
    assert chiral_labels == ["10_2"]

    # 8_1^{1,1} and 10_1^{1,1} agree on Layers 1..8 and are separated at Layer 9 (Sanderson Dlk)!
    w_811 = by_label["8_1^1,1"]["nearest_separating_witness"]
    assert w_811["nearest_label"] == "10_1^1,1"
    assert w_811["separating_layer"] == "double_triple_linking"
    assert w_811["separating_layer_index"] == 9

    # Kamada-Oshiro Q_6 symmetric 3-cocycle obstruction on 8_1^{-1,-1} vs 10_1^{-2,-2}
    assert by_label["8_1^-1,-1"]["certified_positive_triple_points_by_q6"] is True
    assert by_label["8_1^-1,-1"]["sym_quandle_q6_nontrivial_weights"] == 24
    assert by_label["10_1^-2,-2"]["certified_positive_triple_points_by_q6"] is False
    assert by_label["10_1^-2,-2"]["sym_quandle_q6_trivial_weights"] == 24


def test_yoshikawa_three_table_ii_errata_live() -> None:
    for label in ("9_1^1,-2", "10_1^0,1", "10_2^0,-2"):
        entry = YOSHIKAWA_BY_LABEL[label]
        d = MarkedDiagram.from_epd(entry.canonical_epd)
        assert certify_classical_unlink(d.resolve(1)).is_unlink
        assert certify_classical_unlink(d.resolve(-1)).is_unlink
        printed_ok, _ = check_e1_ideal_against_reference(d, entry.yoshikawa_printed_e1)
        corrected_ok, _ = check_e1_ideal_against_reference(d, entry.corrected_e1)
        assert printed_ok is False
        assert corrected_ok is True

    # Specifically verify H_1(S^4 \ F) = Z + Z_2 on 9_1^{1,-2}
    d_9 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["9_1^1,-2"].canonical_epd)
    wp_9 = wirtinger_presentation(d_9)
    assert wp_9.abelianization_invariant_factors() == (1, (2,))
    al_9 = compute_alexander_invariants(d_9)
    assert al_9.orientable_flags == (True, False)
    assert any(g in ("2", "2*y") for g in al_9.e1_generators)


def test_distinguishing_4d_invariants_on_classical_collision_pairs() -> None:
    # 1. Chirality / non-invertibility of 2-twist spun trefoil 10_2 via Mochizuki 3-cocycle over R_3
    d_10_2 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_2"].canonical_epd)
    c3_oriented = mochizuki_3cocycle_invariant(d_10_2)
    assert c3_oriented is not None
    assert set(c3_oriented) == {"3 + 6*u", "3 + 6*u^2"}

    # 2. Sanderson mod-2 double-linking separates 8_1^{1,1} (Dlk = 0) from 10_1^{1,1} (Dlk = 1)
    d_8_11 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["8_1^1,1"].canonical_epd)
    d_10_11 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_1^1,1"].canonical_epd)
    assert [row[2] for row in double_linking_mod2(d_8_11)] == [0, 0]
    assert [row[2] for row in double_linking_mod2(d_10_11)] == [1, 1]

    # 3. Peripheral SL(2, F_3) longitudinal holonomy separates 10_1^1 = Tube(3_1) from 8_1 = tau^0(3_1)
    d_8_1 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["8_1"].canonical_epd)
    d_10_1_1 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_1^1"].canonical_epd)
    hom8, per8 = peripheral_sl2f3_profile(d_8_1)
    hom10, per10 = peripheral_sl2f3_profile(d_10_1_1)
    assert hom8 == 72 and len(per8) == 0
    assert hom10 == 72 and len(per10) == 24

    # 4. Exact sublink structure of 3-component link 10_1^{0,0,1} (NOT Brunnian!)
    d_3comp = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_1^0,0,1"].canonical_epd)
    subs = sublink_summary(d_3comp)
    assert subs["all_singletons_unknotted"] is True
    assert subs["is_brunnian"] is False
    proper = subs["proper_sublinks"]
    pairs = [p for p in proper if len(p["components"]) == 2]
    s2_s2 = [p for p in pairs if list(p["types"]) == ["S^2", "S^2"]]
    s2_t2 = [p for p in pairs if sorted(p["types"]) == ["S^2", "T^2"]]
    assert len(s2_s2) == 1 and s2_s2[0]["certified_split"] is True
    assert len(s2_t2) == 2 and all(
        p["min_ch_found"] == 6 and p["certified_split"] is False for p in s2_t2
    )
