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
    assert errata_labels == ["9_1^1,-2", "10_2^0,-2"]
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


def test_yoshikawa_two_table_ii_errata_live() -> None:
    for label in ("9_1^1,-2", "10_2^0,-2"):
        entry = YOSHIKAWA_BY_LABEL[label]
        d = MarkedDiagram.from_epd(entry.canonical_epd)
        assert certify_classical_unlink(d.resolve(1)).is_unlink
        assert certify_classical_unlink(d.resolve(-1)).is_unlink
        printed_ok, _ = check_e1_ideal_against_reference(d, entry.yoshikawa_printed_e1)
        corrected_ok, _ = check_e1_ideal_against_reference(d, entry.corrected_e1)
        assert printed_ok is False
        assert corrected_ok is True

    # Verify 10_1^{0,1}'s printed E_1 = (x - 1, y - 1)(xy + 1) matches (x - 1, y - 1)(x + y) under y -> y^{-1}
    e_101 = YOSHIKAWA_BY_LABEL["10_1^0,1"]
    d_101 = MarkedDiagram.from_epd(e_101.canonical_epd)
    p101_ok, _ = check_e1_ideal_against_reference(d_101, e_101.yoshikawa_printed_e1)
    assert p101_ok is True

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


def test_o4_mirror_invariance_linking_and_euler_pairing() -> None:
    from surfaceknots.ch11_census import disjoint_union, planar_connected_sum
    from surfaceknots.invariants import exact_sublink_invariant_profile, mod2_triple_linking_profile
    from surfaceknots.yoshikawa_data import _unoriented_euler_key, compute_full_invariant_record

    p_plus = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^-1"].canonical_epd)
    d_811 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["8_1^1,1"].canonical_epd)

    # 1. (8_1^{1,1}) # 2_1^{-1} has mixed T^2(e=0) u N_3(e=+-2) with nontrivial Sanderson Dlk:
    d_mix = planar_connected_sum(d_811, 0, p_plus, 0, flip_splice=False)
    d_mix_mir = d_mix.mirror_3d()
    assert double_linking_mod2(d_mix) == double_linking_mod2(d_mix_mir)
    assert mod2_triple_linking_profile(d_mix) == mod2_triple_linking_profile(d_mix_mir)
    rec_pos = compute_full_invariant_record(d_mix, label="mix_pos", include_sublinks=False)
    rec_neg = compute_full_invariant_record(d_mix_mir, label="mix_neg", include_sublinks=False)
    assert invariant_discrimination_key(rec_pos) == invariant_discrimination_key(rec_neg)

    # 2. Check ch=12 entry 12_1^1,-2 (T^2 u K^2(|e|=4)) if present:
    ch12_path = DATA_DIR / "expanded_ch12.json"
    if ch12_path.exists():
        ch12_map = {e["label"]: e for e in json.loads(ch12_path.read_text())["entries"]}
        if "12_1^1,-2" in ch12_map:
            d_12 = MarkedDiagram.from_epd(ch12_map["12_1^1,-2"]["epd"])
            d_12_mir = d_12.mirror_3d()
            assert double_linking_mod2(d_12) == double_linking_mod2(d_12_mir)
            r12_pos = compute_full_invariant_record(d_12, label="12_pos", include_sublinks=False)
            r12_neg = compute_full_invariant_record(d_12_mir, label="12_neg", include_sublinks=False)
            assert invariant_discrimination_key(r12_pos) == invariant_discrimination_key(r12_neg)

    # 3. Paired (genus, normal_euler) with sum(e_list) == 0:
    #    RP^2(e=+2) u RP^2(e=+2) u K^2(e=-4) vs its mirror AND vs RP^2(+2) u RP^2(-2) u K^2(0)
    k_minus4 = planar_connected_sum(
        p_plus.mirror_3d(), 0, p_plus.mirror_3d(), 0, flip_splice=False
    )
    if k_minus4.surface_topology().components[0].normal_euler != -4:
        k_minus4 = planar_connected_sum(
            p_plus.mirror_3d(), 0, p_plus.mirror_3d(), 0, flip_splice=True
        )
    k_zero = planar_connected_sum(p_plus, 0, p_plus.mirror_3d(), 0, flip_splice=False)
    if k_zero.surface_topology().components[0].normal_euler != 0:
        k_zero = planar_connected_sum(p_plus, 0, p_plus.mirror_3d(), 0, flip_splice=True)

    d_a = disjoint_union(disjoint_union(p_plus, p_plus), k_minus4)
    d_a_mir = d_a.mirror_3d()
    d_b = disjoint_union(disjoint_union(p_plus, p_plus.mirror_3d()), k_zero)

    assert exact_sublink_invariant_profile(d_a) == exact_sublink_invariant_profile(d_a_mir)
    assert exact_sublink_invariant_profile(d_a) != exact_sublink_invariant_profile(d_b)

    rec_legacy_pos = {"yoshikawa_genera": ["-1", "-1", "-2"], "normal_euler_numbers": [2, 2, -4]}
    rec_legacy_neg = {"yoshikawa_genera": ["-1", "-1", "-2"], "normal_euler_numbers": [-2, -2, 4]}
    rec_legacy_alt = {"yoshikawa_genera": ["-1", "-1", "-2"], "normal_euler_numbers": [2, -2, 0]}
    assert _unoriented_euler_key(rec_legacy_pos) == _unoriented_euler_key(rec_legacy_neg)
    assert _unoriented_euler_key(rec_legacy_pos) != _unoriented_euler_key(rec_legacy_alt)

    # Also verify paired (genus, euler) prevents cross-component swap (-1, +2), (-2, -2) vs (-1, -2), (-2, +2):
    rec_swap_1 = {"yoshikawa_genera": ["-1", "-2"], "normal_euler_numbers": [2, -2]}
    rec_swap_2 = {"yoshikawa_genera": ["-1", "-2"], "normal_euler_numbers": [-2, 2]}
    assert _unoriented_euler_key(rec_swap_1) == _unoriented_euler_key(rec_swap_2)
    rec_diff_1 = {"yoshikawa_genera": ["-1", "-2", "0"], "normal_euler_numbers": [2, 2, -4]}
    rec_diff_2 = {"yoshikawa_genera": ["-1", "-2", "0"], "normal_euler_numbers": [2, -4, 2]}
    assert _unoriented_euler_key(rec_diff_1) != _unoriented_euler_key(rec_diff_2)


def test_sublink_summary_and_triple_linking_canonicalization() -> None:
    from surfaceknots.ch11_census import disjoint_union
    from surfaceknots.invariants import mod2_triple_linking_profile

    # 1. [Sub6 H4] 4-component link 6_1^{0,1} u 0_1 u 0_1:
    #    3-component proper sublink (6_1^{0,1} u 0_1) is split (2 pieces) yet NOT a trivial unlink!
    d_hopf = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["6_1^0,1"].canonical_epd)
    d_unk = MarkedDiagram.from_epd("O[1]")
    d_4 = disjoint_union(disjoint_union(d_hopf, d_unk), d_unk)
    summary_4 = sublink_summary(d_4, max_states=64)
    triples_t2_s2_s2 = {
        p["components"]: p
        for p in summary_4["proper_sublinks"]
        if len(p["components"]) == 3 and sorted(p["types"]) == ["S^2", "S^2", "T^2"]
    }
    assert set(triples_t2_s2_s2) == {(0, 1, 2), (0, 1, 3), (0, 2, 3)}
    for hopf_sub in ((0, 1, 2), (0, 1, 3)):
        assert triples_t2_s2_s2[hopf_sub]["certified_split"] is True
        assert triples_t2_s2_s2[hopf_sub]["certified_trivial_unlink"] is False
    assert triples_t2_s2_s2[(0, 2, 3)]["certified_split"] is True
    assert triples_t2_s2_s2[(0, 2, 3)]["certified_trivial_unlink"] is True

    # 2. [Sub6 M1] mod2_triple_linking_profile canonicalizes over component permutations,
    #    free_loops, and O(4) mirror:
    d_1001 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_1^0,0,1"].canonical_epd)
    prof_base = mod2_triple_linking_profile(d_1001)
    assert prof_base == mod2_triple_linking_profile(d_1001.mirror_3d())
    d_with_free = MarkedDiagram(
        d_1001.num_vertices, d_1001.opp, d_1001.vtype, d_1001.free_loops + 2
    )
    assert mod2_triple_linking_profile(d_with_free) == prof_base


def test_algebraic_and_biquandle_soundness_regressions() -> None:
    import pytest

    import surfaceknots.invariants as inv_pkg
    from surfaceknots.ch11_census import disjoint_union
    from surfaceknots.invariants import (
        FiniteGroup,
        GroupPresentation,
        alexander_biquandle_invariants,
        count_homomorphisms,
    )

    # 1. Public exports in surfaceknots.invariants:
    for sym_name in (
        "FiniteGroup",
        "STANDARD_TEST_GROUPS",
        "canonical_zh1_groebner_key",
        "count_homomorphisms",
        "exact_sublink_invariant_profile",
        "mochizuki_3cocycle_orientation_pairs",
        "strong_groebner_zz",
        "tietze_simplify",
        "zh1_groebner_basis",
    ):
        assert hasattr(inv_pkg, sym_name), f"Missing export {sym_name}"

    # 2. GroupPresentation.to_human_string validation & reduction, and _CONJ_CLASSES_CACHE key:
    assert GroupPresentation(1, ((1, -1),), ()).to_human_string() == "<x1 | 1 = 1>"
    with pytest.raises(ValueError):
        GroupPresentation(1, ((2,),), ()).to_human_string()

    g_z2 = FiniteGroup("SameName", 2, ((0, 1), (1, 0)), (0, 1))
    g_z3 = FiniteGroup("SameName", 3, ((0, 1, 2), (1, 2, 0), (2, 0, 1)), (0, 2, 1))
    p_z6 = GroupPresentation(1, ((1, 1, 1, 1, 1, 1),), ())
    assert count_homomorphisms(p_z6, g_z2) == 2
    assert count_homomorphisms(p_z6, g_z3) == 3

    # 3. compute_alexander_invariants: g == 0 2-tuple, rational meridian ValueError, odd torsion ValueError:
    al_empty = compute_alexander_invariants(MarkedDiagram(0, (), (), 0), pres=GroupPresentation(0, (), ()))
    assert isinstance(al_empty.evaluation_fingerprint, tuple) and len(al_empty.evaluation_fingerprint) == 2

    d_s2 = MarkedDiagram.from_epd("O[1]")
    p_rat = GroupPresentation(2, ((1, 1, -2),), ((2,),))  # x^2 = mu -> x = mu^{1/2}
    with pytest.raises(ValueError):
        compute_alexander_invariants(d_s2, pres=p_rat)

    d_s2_rp2 = disjoint_union(d_s2, MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^-1"].canonical_epd))
    p_z_z3 = GroupPresentation(2, ((2, 2, 2),), ((1,), (2,)))
    with pytest.raises(ValueError):
        compute_alexander_invariants(d_s2_rp2, pres=p_z_z3)

    # 4. Canonical single_var_alex under t <-> t^{-1} on 9_1:
    d_9_1 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["9_1"].canonical_epd)
    sv_pos = compute_alexander_invariants(d_9_1).single_var_alex
    assert compute_alexander_invariants(d_9_1.mirror_3d()).single_var_alex == sv_pos
    assert compute_alexander_invariants(d_9_1.reverse_orientation_zt()).single_var_alex == sv_pos

    # 5. Genuine Groebner keys enabled for r + s == 3 ((3, 0) and (2, 1)):
    d_30 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_1^0,0,1"].canonical_epd)
    al_30 = compute_alexander_invariants(d_30)
    assert len(al_30.evaluation_fingerprint[0][0]) > 0
    d_21 = disjoint_union(
        MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["6_1^0,1"].canonical_epd),
        MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^-1"].canonical_epd),
    )
    al_21 = compute_alexander_invariants(d_21)
    assert len(al_21.evaluation_fingerprint[0][0]) > 0

    # 6. Tietze consequence invariance past total_basis > 160:
    wp_91 = wirtinger_presentation(d_9_1, simplify=False)
    r0 = wp_91.relators[0]
    high_deg_conseq = (1,) * 16 + r0 + (-1,) * 16
    wp_aug = GroupPresentation(
        wp_91.num_generators, wp_91.relators + (high_deg_conseq,), wp_91.component_meridians
    )
    assert (
        compute_alexander_invariants(d_9_1, pres=wp_aug).evaluation_fingerprint
        == compute_alexander_invariants(d_9_1, pres=wp_91).evaluation_fingerprint
    )

    # 7. alexander_biquandle_invariants row dedup + rank check (rem_g = 9 returns exact ('0',)):
    from surfaceknots.moves import list_omega2_pokes

    d_811 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["8_1^1,1"].canonical_epd)
    d_dup = disjoint_union(d_811, d_811)
    d_dup_5 = MarkedDiagram(d_dup.num_vertices, d_dup.opp, d_dup.vtype, d_dup.free_loops + 5)
    d_one_7 = MarkedDiagram(d_811.num_vertices, d_811.opp, d_811.vtype, d_811.free_loops + 7)
    assert alexander_biquandle_invariants(d_dup_5)[1] == ("0",)
    assert alexander_biquandle_invariants(d_one_7)[1] == ("0",)
    d_poked = list_omega2_pokes(list_omega2_pokes(d_one_7)[0])[0]
    assert alexander_biquandle_invariants(d_poked) == alexander_biquandle_invariants(d_one_7)

