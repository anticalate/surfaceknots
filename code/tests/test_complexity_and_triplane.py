"""Tests for topological complexity floors, excess ch-index, targeted C99 enumeration,
and bridge-trisection tri-plane diagrams / censuses.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from surfaceknots.complexity import (
    MINIMAL_KLEIN_BOTTLE_E0_EPD,
    analyze_catalog_complexity,
    analyze_connected_sum_collapses,
    build_complexity_comparison,
    certify_klein_bottle_subadditivity,
    component_unknot_floor,
    construct_minimal_unknot_diagram,
    diagram_c_m_split,
    excess_ch_index,
    unknot_ch_floor,
)
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.fast_enum import (
    batch_check_tangle_pairs_unlink,
    enumerate_targeted_prime_ch_census,
    generate_prime_shadows_c_list,
)
from surfaceknots.triplane import (
    TangleDiagram,
    TriPlaneDiagram,
    compute_triplane_invariant_record,
    concentrated_triplane_to_marked_diagram,
    generate_catalan_matchings,
    standard_rp2_triplane,
    standard_torus_triplane,
    standard_unknot_s2_triplane,
    triplane_o4_key,
    union_two_tangles_to_classical_link,
)
from surfaceknots.triplane_census import (
    enumerate_minimal_braided_tangles,
    verify_general_planar_tangles_b4,
)
from surfaceknots.unlink import certify_classical_unlink
from surfaceknots.yoshikawa_data import YOSHIKAWA_BY_LABEL, compute_full_invariant_record


DATA_DIR = Path(__file__).resolve().parents[1] / "data"


def test_exact_unknot_ch_floors_and_realizations() -> None:
    # Orientable surfaces Sigma_g (g = 0..4, e = 0)
    for g in range(5):
        ch0, c0, m0 = component_unknot_floor(g, 0)
        assert (ch0, c0, m0) == (2 * g, 0, 2 * g)
        d = construct_minimal_unknot_diagram(g, 0)
        assert d.ch_index == ch0
        assert d.num_classical == c0
        assert d.num_marked == m0
        st = d.surface_topology()
        assert st.num_components == 1
        assert st.components[0].orientable is True
        assert st.components[0].genus_symbol == g
        assert st.components[0].normal_euler == 0

    # Non-orientable unknots #_k RP^2 across all Massey/Whitney-valid normal Euler numbers
    for k in range(1, 6):
        for e in range(-2 * k, 2 * k + 1, 4):
            if (e - 2 * k) % 4 != 0:
                continue
            ch0, c0, m0 = component_unknot_floor(-k, e)
            assert c0 == max(1, abs(e) // 2)
            assert m0 == k
            assert ch0 == k + max(1, abs(e) // 2)
            d = construct_minimal_unknot_diagram(-k, e)
            assert d.ch_index == ch0
            assert d.num_classical == c0
            assert d.num_marked == m0
            st = d.surface_topology()
            assert st.num_components == 1
            assert st.components[0].orientable is False
            assert st.components[0].genus_symbol == -k
            assert st.components[0].normal_euler == e

    # Whitney congruence violation raises ValueError
    with pytest.raises(ValueError):
        component_unknot_floor(-1, 0)
    with pytest.raises(ValueError):
        component_unknot_floor(1, 2)


def test_klein_bottle_subadditivity_and_hopf_excess_ch() -> None:
    cert = certify_klein_bottle_subadditivity()
    assert cert["sum_of_ch"] == 4
    assert cert["true_ch_of_connected_sum"] == 3
    assert cert["defect"] == 1
    assert cert["minimal_epd"] == MINIMAL_KLEIN_BOTTLE_E0_EPD
    assert cert["total_yoshikawa_steps"] == 4
    for epd in cert["move_path_epds"]:  # type: ignore[union-attr]
        d = MarkedDiagram.from_epd(str(epd))
        assert certify_classical_unlink(d.resolve(1)).is_unlink
        assert certify_classical_unlink(d.resolve(-1)).is_unlink

    # Verify all 6 Hopf-type 2-links in Yoshikawa's table have excess_ch == 4
    hopf_labels = ("6_1^0,1", "7_1^0,-2", "8_1^1,1", "8_1^-1,-1", "9_1^1,-2", "10_1^-2,-2")
    for lbl in hopf_labels:
        d = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL[lbl].canonical_epd)
        cm = diagram_c_m_split(d)
        assert cm["excess_ch"] == 4, (lbl, cm)
        st = d.surface_topology()
        assert excess_ch_index(
            d.ch_index,
            [c.genus_symbol for c in st.components],
            [c.normal_euler for c in st.components],
        ) == 4


def test_viro_collapse_and_rp2_stabilization() -> None:
    rows = {r["expression"]: r for r in analyze_connected_sum_collapses()}
    assert rows["10_3 # P+"]["viro_barden_collapse_to_unknot"] is True
    assert rows["10_3 # P+"]["certified_true_ch"] == 2
    assert rows["10_3 # P+"]["ch_collapse_amount"] == 10
    # 8_1 # P+, 9_1 # P+, and 10_2 # P+ all retain S_3 quotient (hom_S3 == 10)
    for expr in ("8_1 # P+", "9_1 # P+", "10_2 # P+"):
        assert rows[expr]["viro_barden_collapse_to_unknot"] is False
        assert rows[expr]["hom_S3"] == 10


def test_c99_shadow_generator_and_targeted_census() -> None:
    expected_a113201 = {4: 2, 5: 3, 6: 9, 7: 18, 8: 62, 9: 198, 10: 803}
    for n, count in expected_a113201.items():
        assert len(generate_prime_shadows_c_list(n)) == count

    # Targeted ch = 8, m = 2 connected S^2 scan recovers 8_1 cleanly
    cands8, _ = enumerate_targeted_prime_ch_census(
        ch=8,
        m_min=2,
        m_max=2,
        target_components=1,
        target_orientability=1,
        target_chi_min=2,
        target_chi_max=2,
        num_workers=1,
    )
    assert len(cands8) >= 1
    d8 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["8_1"].canonical_epd)
    rec_target = compute_full_invariant_record(d8, label="8_1", include_sublinks=False)
    rec_found = compute_full_invariant_record(cands8[0], label="cand8", include_sublinks=False)
    assert triplane_o4_key(rec_target) == triplane_o4_key(rec_found)


def test_triplane_invariants_compiler_and_meander_obstruction() -> None:
    for lbl, y_lbl, tp in (
        ("0_1", "0_1", standard_unknot_s2_triplane()),
        ("2_1^-1+", "2_1^-1", standard_rp2_triplane(1)),
        ("2_1^-1-", "2_1^-1", standard_rp2_triplane(-1)),
        ("2_1^1", "2_1^1", standard_torus_triplane()),
    ):
        assert tp.is_certified_valid()
        trec = compute_triplane_invariant_record(tp, lbl)
        md = concentrated_triplane_to_marked_diagram(tp)
        assert md is not None
        mrec = compute_full_invariant_record(md, label=lbl, include_sublinks=False)
        yrec = compute_full_invariant_record(
            MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL[y_lbl].canonical_epd),
            label=y_lbl,
            include_sublinks=False,
        )
        assert triplane_o4_key(trec) == triplane_o4_key(mrec) == triplane_o4_key(yrec)

    # Explicit b = 4 Meander Interleaving Obstruction witness:
    # Cycle A = {0, 3, 4, 7} and Cycle B = {1, 2, 5, 6} alternate A-B-A-B around S^1!
    m_inter_1 = ((0, 7), (1, 6), (2, 5), (3, 4))
    m_inter_3 = ((0, 3), (1, 2), (4, 7), (5, 6))
    tp_inter = TriPlaneDiagram(
        TangleDiagram(4, m_inter_1, ()),
        TangleDiagram(4, m_inter_1, ()),
        TangleDiagram(4, m_inter_3, ()),
    )
    assert tp_inter.is_certified_valid()
    assert concentrated_triplane_to_marked_diagram(tp_inter) is None


def test_braided_tangles_and_exported_complexity_datasets() -> None:
    by_c, _ = enumerate_minimal_braided_tangles(4, max_crossings=2)
    assert [len(by_c[0]), len(by_c[1]), len(by_c[2])] == [14, 96, 582]

    sub = by_c[0] + by_c[1][:20]
    raw = batch_check_tangle_pairs_unlink(sub, sub)
    n = len(sub)
    for i in range(n):
        for j in range(n):
            got = bool(raw[i * n + j])
            exp = certify_classical_unlink(union_two_tangles_to_classical_link(sub[i], sub[j])).is_unlink
            assert got == exp

    t12_path = DATA_DIR / "complexity" / "targeted_ch12_13_census.json"
    tp_path = DATA_DIR / "complexity" / "triplane_b4_census.json"
    comp_path = DATA_DIR / "complexity" / "complexity_comparison.json"
    assert t12_path.exists()
    assert tp_path.exists()
    assert comp_path.exists()

    t12 = json.loads(t12_path.read_text())
    assert len(t12["ch12_m2_s2_knots"]["new_prime_classes"]) >= 7

    tp_data = json.loads(tp_path.read_text())
    glmo = tp_data["gong_lewis_monkman_osnes_2026_verification"]
    assert glmo["conjecture_3_6_spun_trefoil_c_tri_eq_6_at_b4"] is True
    assert glmo["corollary_3_5_two_twist_spun_trefoil_c_tri_eq_6_at_b4"] is True
    assert glmo["conjecture_3_7_only_0_1_8_1_10_2_for_s2_up_to_c6_at_b4"] is True
    assert glmo["rank_inversion_9_1_absent_at_c6"] is True

    comp_data = json.loads(comp_path.read_text())
    assert comp_data["total_distinct_classes_compared"] >= 70


def test_general_planar_tangles_b4_and_dlk2_o4_key() -> None:
    res1 = verify_general_planar_tangles_b4(max_crossings=1)
    assert res1["all_general_planar_trivial_tangles_realized_by_braided_up_to_max_c"] is True
    assert res1["per_crossing_counts"]["1"]["non_braided_trivial_diagrams"] == 0  # type: ignore[index]

    res2 = verify_general_planar_tangles_b4(max_crossings=2)
    assert res2["per_crossing_counts"]["2"]["free_pi1_rank4_diagrams"] == 224  # type: ignore[index]
    assert res2["per_crossing_counts"]["2"]["certified_isotopic_to_braided_le_c"] == 208  # type: ignore[index]
    assert res2["per_crossing_counts"]["2"]["non_braided_trivial_diagrams"] == 16  # type: ignore[index]

    res3 = verify_general_planar_tangles_b4(max_crossings=3)
    assert res3["per_crossing_counts"]["3"]["free_pi1_rank4_diagrams"] == 1232  # type: ignore[index]
    assert res3["per_crossing_counts"]["3"]["certified_isotopic_to_braided_le_c"] == 1084  # type: ignore[index]
    assert res3["per_crossing_counts"]["3"]["non_braided_trivial_diagrams"] == 148  # type: ignore[index]

    # (8_1^1,1, 10_1^1,1) collide on the base 6-tuple triplane_o4_key (both b >= 6),
    # and separate cleanly when include_double_linking=True:
    r8 = compute_full_invariant_record(
        MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["8_1^1,1"].canonical_epd),
        label="8_1^1,1",
        include_sublinks=False,
    )
    r10 = compute_full_invariant_record(
        MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_1^1,1"].canonical_epd),
        label="10_1^1,1",
        include_sublinks=False,
    )
    assert triplane_o4_key(r8) == triplane_o4_key(r10)
    assert triplane_o4_key(r8, include_double_linking=True) != triplane_o4_key(r10, include_double_linking=True)


def test_live_build_complexity_comparison() -> None:
    live = build_complexity_comparison(DATA_DIR)
    rows_by_lbl = {r["label"]: r for r in live["rows"]}  # type: ignore[union-attr]

    # Exact bridge & c_tri cross-join checks via witness_triplane recomputation:
    assert rows_by_lbl["0_1"]["b4_census_bridge"] == 1
    assert rows_by_lbl["0_1"]["c_tri_b4_braided"] == 0
    assert rows_by_lbl["2_1^1"]["bridge_bounds"] == [3, 3]
    assert rows_by_lbl["2_1^1"]["b4_census_bridge"] == 3
    assert rows_by_lbl["2_1^1"]["c_tri_b4_braided"] == 0
    assert rows_by_lbl["2_1^-1"]["b4_census_bridge"] == 2
    assert rows_by_lbl["2_1^-1"]["c_tri_b4_braided"] == 1
    assert int(rows_by_lbl["8_1"]["c_tri_b4_braided"]) <= 6
    assert rows_by_lbl["10_2"]["c_tri_b4_braided"] == 6

    # All 9 b <= 4 census bridge matches must resolve cleanly:
    b4_matched = [r["label"] for r in live["rows"] if r["b4_census_bridge"] is not None]  # type: ignore[union-attr]
    assert len(b4_matched) >= 9

    # Explicit epistemic certification tiers and zero duplicate promoted Spoke rows:
    valid_tiers = {"PROVEN_EXACT", "HEURISTIC_INVARIANT_MATCH", "UPPER_BOUND_ONLY"}
    for r in live["rows"]:  # type: ignore[union-attr]
        assert r["certification_tier"] in valid_tiers
    ch12_path = DATA_DIR / "expanded_ch12.json"
    if ch12_path.exists():
        ch12_data = json.loads(ch12_path.read_text())
        for pr in ch12_data.get("promoted_exact_ch12_spokes", []):
            assert pr["family_id"] not in rows_by_lbl
