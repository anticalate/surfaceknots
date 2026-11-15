"""Tests for parametric surface-link family compilers and Hub-and-Spoke cross-join."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from surfaceknots.ch11_census import planar_connected_sum
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.families import (
    CLASSICAL_PLAT_BY_LABEL,
    RIBBON_1FUSION_TABLE,
    compile_1fusion_ribbon_word,
    compile_artin_spin,
    compile_mixed_plat_link,
    compile_plat_surface,
    compile_satoh_tube,
    compile_symmetric_union_2knot,
    fox_1fusion_alexander_forms,
    fox_1fusion_alexander_poly,
    parse_ribbon_word,
)
from surfaceknots.invariants.alexander import compute_alexander_invariants
from surfaceknots.invariants.linking import peripheral_sl2f3_profile
from surfaceknots.invariants.tier3 import lomonaco_pi2_invariants
from surfaceknots.moves import _is_planar_4reg_diagram
from surfaceknots.unlink import certify_classical_unlink
from surfaceknots.yoshikawa_data import YOSHIKAWA_BY_LABEL


DATA_DIR = Path(__file__).resolve().parents[1] / "data"


def test_plat_compiler_planarity_and_admissibility() -> None:
    """Twisted Ribbon 1-Handle Pinch Theorem: every compiled plat surface is planar and admissible."""
    for lbl in ("0_1", "3_1", "4_1", "5_1", "L2a1", "L4a1", "L5a1", "L6a4"):
        spec = CLASSICAL_PLAT_BY_LABEL[lbl]
        d_tube = compile_satoh_tube(spec)
        expected_ch = 2 * len(spec.braid_word) + 2 * spec.bridge_number
        assert d_tube.ch_index == expected_ch
        assert _is_planar_4reg_diagram(d_tube)
        assert certify_classical_unlink(d_tube.resolve(1)).is_unlink
        assert certify_classical_unlink(d_tube.resolve(-1)).is_unlink
        topo = d_tube.surface_topology()
        assert topo.num_components == spec.num_components
        assert all(c.homeomorphism_type == "T^2" for c in topo.components)

        if spec.num_components == 1:
            d_spin = compile_artin_spin(spec)
            assert _is_planar_4reg_diagram(d_spin)
            assert certify_classical_unlink(d_spin.resolve(1)).is_unlink
            assert certify_classical_unlink(d_spin.resolve(-1)).is_unlink
            assert tuple(c.homeomorphism_type for c in d_spin.surface_topology().components) == ("S^2",)

            for m in (-2, -1, 1, 2):
                d_sym = compile_symmetric_union_2knot(spec, twist=m)
                assert _is_planar_4reg_diagram(d_sym)
                assert certify_classical_unlink(d_sym.resolve(1)).is_unlink
                assert certify_classical_unlink(d_sym.resolve(-1)).is_unlink
                if spec.bridge_number >= 2:
                    assert tuple(c.homeomorphism_type for c in d_sym.surface_topology().components) == ("S^2",)

    # Mixed and negative-twist Klein-bottle check on Solomon's knot L4a1:
    solomon = CLASSICAL_PLAT_BY_LABEL["L4a1"]
    for twists, expected_types in (
        ({0: 0}, ("S^2", "T^2")),
        ({0: 1}, ("K^2", "S^2")),
        ({0: -1, 1: 2}, ("K^2", "T^2")),
    ):
        dm = compile_mixed_plat_link(solomon, twists)
        assert dm.ch_index == 2 * len(solomon.braid_word) + sum(2 + abs(v) for v in twists.values())
        assert _is_planar_4reg_diagram(dm)
        assert certify_classical_unlink(dm.resolve(1)).is_unlink
        assert certify_classical_unlink(dm.resolve(-1)).is_unlink
        assert tuple(sorted(c.homeomorphism_type for c in dm.surface_topology().components)) == expected_types

    # Idle-strand edge case check:
    d_idle = compile_plat_surface(4, (), {0: 1})
    assert _is_planar_4reg_diagram(d_idle)
    assert certify_classical_unlink(d_idle.resolve(1)).is_unlink
    assert certify_classical_unlink(d_idle.resolve(-1)).is_unlink


def test_1fusion_ribbon_word_compiler_and_fox_formula() -> None:
    """NORTH/SOUTH 1-fusion compiler produces planar 4|w|+2 diagrams matching Fox's formula."""
    extra_words = [
        ("R_extra_1", parse_ribbon_word("a^3 b^-2 a^-1")),
        ("R_extra_2", parse_ribbon_word("b^-2 a^2 b")),
    ]
    all_words = [(r.family_id, r.word) for r in RIBBON_1FUSION_TABLE] + extra_words
    for fid, word in all_words:
        d = compile_1fusion_ribbon_word(word)
        assert d.ch_index == 4 * len(word) + 2, f"Unexpected raw ch_index on {fid}"
        assert _is_planar_4reg_diagram(d), f"Non-planar 1-fusion diagram on {fid}"
        assert certify_classical_unlink(d.resolve(1)).is_unlink, f"L_+ non-unlink on {fid}"
        assert certify_classical_unlink(d.resolve(-1)).is_unlink, f"L_- non-unlink on {fid}"
        topo = d.surface_topology()
        ctypes = tuple(c.homeomorphism_type for c in topo.components)
        assert ctypes == ("S^2",), f"Expected S^2 on {fid}, got {ctypes}"

        valid_alex = fox_1fusion_alexander_forms(word)
        actual_alex = compute_alexander_invariants(d).single_var_alex
        assert actual_alex in valid_alex, f"{fid}: expected {valid_alex}, got {actual_alex}"


def test_tube_4_1_separated_from_composite_10_1_conn_2_1_1() -> None:
    """Tube(4_1) and 10_1 # 2_1^1 both have ch=12, T^2, and Delta=x^2-3x+1, separated by Tier 3 & Layer 5."""
    d_tube = compile_satoh_tube(CLASSICAL_PLAT_BY_LABEL["4_1"])
    d_10_1 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_1"].canonical_epd)
    d_2_1_1 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^1"].canonical_epd)
    d_comp = planar_connected_sum(d_10_1, 0, d_2_1_1, 0)

    assert d_tube.ch_index == 12
    assert d_comp.ch_index == 12
    assert tuple(c.homeomorphism_type for c in d_tube.surface_topology().components) == ("T^2",)
    assert tuple(c.homeomorphism_type for c in d_comp.surface_topology().components) == ("T^2",)
    assert compute_alexander_invariants(d_tube).single_var_alex == "t**2 - 3*t + 1"
    assert compute_alexander_invariants(d_comp).single_var_alex == "t**2 - 3*t + 1"

    pi2_tube = lomonaco_pi2_invariants(d_tube)
    pi2_comp = lomonaco_pi2_invariants(d_comp)
    assert pi2_tube.cover_betti_profile != pi2_comp.cover_betti_profile
    assert dict(pi2_tube.cover_betti_profile)["A4"] == ((4, (2, 2), 30, 3),)
    assert dict(pi2_comp.cover_betti_profile)["A4"] == ((4, (2, 2), 33, 6),)

    _, periph_tube = peripheral_sl2f3_profile(d_tube)
    _, periph_comp = peripheral_sl2f3_profile(d_comp)
    assert len(periph_tube) == 24
    assert len(periph_comp) == 0


def test_exported_family_catalog_and_hub_realizations() -> None:
    """Verify exported family_catalog.{json,csv} and hub_realizations.json integrity."""
    fam_json = DATA_DIR / "families" / "family_catalog.json"
    fam_csv = DATA_DIR / "families" / "family_catalog.csv"
    hub_json = DATA_DIR / "families" / "hub_realizations.json"
    assert fam_json.exists()
    assert fam_csv.exists()
    assert hub_json.exists()

    payload = json.loads(fam_json.read_text())
    hub_rel = json.loads(hub_json.read_text())
    with fam_csv.open() as f:
        csv_rows = list(csv.DictReader(f))

    assert len(csv_rows) == payload["total_family_instances"]
    assert len(payload["entries"]) == payload["total_family_instances"]

    expected_hub_labels = {
        "0_1",
        "2_1^1",
        "6_1^0,1",
        "7_1^0,-2",
        "8_1",
        "8_1^1,1",
        "9_1",
        "9_1^1,-2",
        "10_1",
        "10_1^1",
        "10_1^0,1",
        "10_1^-2,-2",
        "11_1",
        "11_2",
        "11_3^0,-2",
    }
    assert expected_hub_labels.issubset(set(hub_rel.keys()))

    exact_12 = set(payload["exact_ch_12_certified_instances"])
    for expected_exact in ("Spin(5_1)", "Spin(5_2)", "Tube(4_1)", "Tube(L4a1)", "Mixed(L5a1)"):
        assert expected_exact in exact_12

    by_id = {e["family_id"]: e for e in payload["entries"]}
    assert by_id["SymUnion(4_1,1)"]["hub_match"] == "11_1"
    assert by_id["R(b^-1*a*b^-1*a^-1)"]["hub_match"] == "11_1"
    assert by_id["R(b^-1*a*b*a)"]["hub_match"] == "11_2"
    assert by_id["KleinMixed(L4a1,1)"]["hub_match"] == "11_3^0,-2"
    expected_unmatched_lower = int(payload.get("hub_max_ch", 11)) + 1
    assert by_id["Tube(L6a4)"]["ch_lower"] == expected_unmatched_lower
    assert by_id["Mixed(L6a4,01)"]["ch_lower"] == expected_unmatched_lower
