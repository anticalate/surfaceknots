"""Tests for marked ch-diagrams, Reidemeister unlink movies, Yoshikawa moves, 4D flypes, and tri-plane diagrams."""

from __future__ import annotations

from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.invariants import (
    compute_alexander_invariants,
    double_linking_mod2,
    finite_quotient_profile,
    peripheral_sl2f3_profile,
    wirtinger_presentation,
)
from surfaceknots.moves import (
    greedy_yoshikawa_reduce,
    list_4d_tait_flypes,
    list_level_yoshikawa_neighbors,
    list_omega1_6_stabilizations,
    list_omega2_pokes,
)
from surfaceknots.triplane import (
    estimate_bridge_bounds_from_marked,
    generate_catalan_matchings,
    standard_rp2_triplane,
    standard_torus_triplane,
    standard_unknot_s2_triplane,
)
from surfaceknots.unlink import certify_classical_unlink
from surfaceknots.yoshikawa_data import YOSHIKAWA_BY_LABEL


def test_epd_roundtrip_and_topology() -> None:
    d_rp2 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^-1"].canonical_epd)
    topo_rp2 = d_rp2.surface_topology()
    assert topo_rp2.num_components == 1
    assert topo_rp2.euler_characteristic == 1
    assert topo_rp2.orientability_label == "Non"
    assert abs(topo_rp2.components[0].normal_euler) == 2

    # 3D mirror negates normal Euler number on non-orientable surfaces
    d_rp2_mir = d_rp2.mirror_3d()
    assert d_rp2_mir.surface_topology().components[0].normal_euler == -topo_rp2.components[0].normal_euler

    # Round-trip EPD parsing
    d_hopf = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["6_1^0,1"].canonical_epd)
    d_rt = MarkedDiagram.from_epd(d_hopf.to_epd())
    assert d_rt.num_vertices == 6
    assert d_rt.surface_topology().yoshikawa_genera == d_hopf.surface_topology().yoshikawa_genera


def test_invariance_under_yoshikawa_moves_stabilizations_and_flypes() -> None:
    d0 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_1^1,1"].canonical_epd)
    base_dlk = [r[2] for r in double_linking_mod2(d0)]
    base_hom = finite_quotient_profile(wirtinger_presentation(d0), ("S3", "Q8"))
    base_e1 = compute_alexander_invariants(d0).evaluation_fingerprint

    # 1. Level Yoshikawa neighbors
    lvl = list_level_yoshikawa_neighbors(d0)
    for d_nbr in lvl[:4]:
        assert certify_classical_unlink(d_nbr.resolve(1)).is_unlink
        assert certify_classical_unlink(d_nbr.resolve(-1)).is_unlink
        assert [r[2] for r in double_linking_mod2(d_nbr)] == base_dlk
        assert finite_quotient_profile(wirtinger_presentation(d_nbr), ("S3", "Q8")) == base_hom
        assert compute_alexander_invariants(d_nbr).evaluation_fingerprint == base_e1

    # 2. Omega_1+ / Omega_6+ stabilizations reduce back via greedy_yoshikawa_reduce
    stabs = list_omega1_6_stabilizations(d0)
    assert len(stabs) > 0
    d_stab = stabs[0]
    assert d_stab.num_vertices == d0.num_vertices + 1
    assert [r[2] for r in double_linking_mod2(d_stab)] == base_dlk
    assert greedy_yoshikawa_reduce(d_stab).num_vertices <= d0.num_vertices

    # 3. Omega_2+ poke invariance on Tube(3_1) peripheral SL(2, F_3) profile
    d_tube = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_1^1"].canonical_epd)
    base_per = peripheral_sl2f3_profile(d_tube)
    pokes = list_omega2_pokes(d_tube)
    assert len(pokes) > 0
    assert peripheral_sl2f3_profile(pokes[0]) == base_per
    assert greedy_yoshikawa_reduce(pokes[0]).num_vertices <= d_tube.num_vertices

    # 4. 4D Tait flypes preserve admissibility and invariants when present
    for entry_lbl in ("9_1", "10_1", "10_2"):
        df = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL[entry_lbl].canonical_epd)
        for fl in list_4d_tait_flypes(df)[:2]:
            assert certify_classical_unlink(fl.resolve(1)).is_unlink
            assert certify_classical_unlink(fl.resolve(-1)).is_unlink
            assert compute_alexander_invariants(fl).evaluation_fingerprint == (
                compute_alexander_invariants(df).evaluation_fingerprint
            )


def test_triplane_diagrams_and_bridge_bounds() -> None:
    assert [len(generate_catalan_matchings(b)) for b in (1, 2, 3, 4)] == [1, 2, 5, 14]

    tp_s2 = standard_unknot_s2_triplane()
    assert tp_s2.is_certified_valid()
    assert tp_s2.bridge_number == 1
    assert tp_s2.crossing_number == 0
    assert tp_s2.euler_characteristic() == 2

    tp_rp2 = standard_rp2_triplane(1)
    assert tp_rp2.is_certified_valid()
    assert tp_rp2.bridge_number == 2
    assert tp_rp2.crossing_number == 1
    assert tp_rp2.euler_characteristic() == 1

    tp_t2 = standard_torus_triplane()
    assert tp_t2.is_certified_valid()
    assert tp_t2.bridge_number == 3
    assert tp_t2.crossing_number == 0
    assert tp_t2.euler_characteristic() == 0

    d_8_1 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["8_1"].canonical_epd)
    b_low, b_up = estimate_bridge_bounds_from_marked(d_8_1)
    assert b_low == 4 and b_up >= b_low


def test_free_loop_stabilization_invariance() -> None:
    from surfaceknots.ch11_census import disjoint_union, invariant_discrimination_key
    from surfaceknots.yoshikawa_data import compute_full_invariant_record

    def key(d: MarkedDiagram) -> tuple:
        return invariant_discrimination_key(compute_full_invariant_record(d))

    # A 1-crossing kinked unknotted sphere `0_1` (`n=1, free_loops=0`) vs `MarkedDiagram.unknot_sphere()` (`n=0, free_loops=1`)
    kink_s2 = MarkedDiagram(num_vertices=1, vtype=("X0",), opp=(1, 0, 3, 2), free_loops=0)
    free_s2 = MarkedDiagram.unknot_sphere()
    t2 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^1"].canonical_epd)

    assert key(free_s2) == key(kink_s2)
    assert key(disjoint_union(free_s2, free_s2)) == key(disjoint_union(kink_s2, free_s2))
    assert key(disjoint_union(free_s2, t2)) == key(disjoint_union(kink_s2, t2))

