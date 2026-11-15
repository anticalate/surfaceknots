"""Tests for marked ch-diagrams, Reidemeister unlink movies, Yoshikawa moves, 4D flypes, and tri-plane diagrams."""

from __future__ import annotations

import pytest

import surfaceknots.unlink as unlink_mod
from surfaceknots.ch11_census import disjoint_union, planar_connected_sum
from surfaceknots.diagrams import ClassicalLinkDiagram, MarkedDiagram
from surfaceknots.families import (
    ClassicalPlatSpec,
    compile_1fusion_ribbon_word,
    compile_artin_spin,
    conway_4plat_word,
    parse_ribbon_word,
)
from surfaceknots.invariants import (
    compute_alexander_invariants,
    double_linking_mod2,
    finite_quotient_profile,
    peripheral_sl2f3_profile,
    wirtinger_presentation,
)
from surfaceknots.moves import (
    _is_planar_4reg_diagram,
    _list_omega8_neighbors,
    canonical_marked_key,
    canonical_unoriented_so4_key,
    find_planar_connected_sum_decomposition,
    greedy_yoshikawa_reduce,
    has_certified_cyclic_unknotted_tube,
    list_4d_tait_flypes,
    list_level_yoshikawa_neighbors,
    list_omega1_6_stabilizations,
    list_omega2_pokes,
)
from surfaceknots.triplane import (
    TangleDiagram,
    estimate_bridge_bounds_from_marked,
    generate_catalan_matchings,
    standard_rp2_triplane,
    standard_torus_triplane,
    standard_unknot_s2_triplane,
    union_two_tangles_to_classical_link,
)
from surfaceknots.unlink import (
    _canonical_classical_key,
    _list_r1_plus_steps,
    certify_classical_unlink,
    find_unlink_movie,
)
from surfaceknots.yoshikawa_data import YOSHIKAWA_BY_LABEL


def test_epd_roundtrip_and_topology() -> None:
    d_rp2 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^-1"].canonical_epd)
    topo_rp2 = d_rp2.surface_topology()
    assert topo_rp2.num_components == 1
    assert topo_rp2.euler_characteristic == 1
    assert topo_rp2.orientability_label == "Non"
    # Under the time-first orientation R_t x R^3 (e = w(L_-) - w(L_+)), standard 2_1^-1 has normal_euler == +2:
    assert topo_rp2.components[0].normal_euler == 2

    # 3D mirror negates normal Euler number on non-orientable surfaces
    d_rp2_mir = d_rp2.mirror_3d()
    assert d_rp2_mir.surface_topology().components[0].normal_euler == -topo_rp2.components[0].normal_euler

    # Round-trip EPD parsing
    d_hopf = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["6_1^0,1"].canonical_epd)
    d_rt = MarkedDiagram.from_epd(d_hopf.to_epd())
    assert d_rt.num_vertices == 6
    assert d_rt.surface_topology().yoshikawa_genera == d_hopf.surface_topology().yoshikawa_genera


def test_epd_empty_and_free_loop_numbering() -> None:
    # Empty 0-component diagram round-trips losslessly via "EMPTY":
    empty = MarkedDiagram(num_vertices=0, opp=(), vtype=(), free_loops=0)
    assert empty.to_epd() == "EMPTY"
    empty_rt = MarkedDiagram.from_epd("EMPTY")
    assert empty_rt.num_vertices == 0 and empty_rt.free_loops == 0

    # Legacy "" and "0_1" continue to parse to the 1-component unknotted sphere:
    for s in ("", "0_1"):
        u = MarkedDiagram.from_epd(s)
        assert u.num_vertices == 0 and u.free_loops == 1

    # Multiple free loops with 0 vertices get distinct strand IDs O[1], O[2], O[3]:
    three_s2 = MarkedDiagram(num_vertices=0, opp=(), vtype=(), free_loops=3)
    assert three_s2.to_epd() == "O[1], O[2], O[3]"
    assert MarkedDiagram.from_epd(three_s2.to_epd()).free_loops == 3

    # Mixed crossing/marked vertices + multiple free loops allocate disjoint O[k] strand IDs:
    d_base = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^-1"].canonical_epd)
    d_mixed = MarkedDiagram(d_base.num_vertices, d_base.opp, d_base.vtype, free_loops=2)
    epd_mixed = d_mixed.to_epd()
    assert "O[5], O[6]" in epd_mixed
    d_mixed_rt = MarkedDiagram.from_epd(epd_mixed)
    assert d_mixed_rt.num_vertices == d_base.num_vertices
    assert d_mixed_rt.free_loops == 2


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
    from surfaceknots.ch11_census import invariant_discrimination_key
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


def test_canonical_classical_key_slot_shift_and_unlink_movie_r1_plus(monkeypatch: pytest.MonkeyPatch) -> None:
    # 1. Connected & disconnected slot-rotated classical crossings deduplicate identically:
    kink_a = ClassicalLinkDiagram(1, (1, 0, 3, 2), ("X0",))
    kink_b = ClassicalLinkDiagram(1, (3, 2, 1, 0), ("X1",))  # same +1 kink rotated by +1 slot
    assert _canonical_classical_key(kink_a) == _canonical_classical_key(kink_b)

    disc_aa = ClassicalLinkDiagram(2, (1, 0, 3, 2, 5, 4, 7, 6), ("X0", "X0"))
    disc_ab = ClassicalLinkDiagram(2, (1, 0, 3, 2, 7, 6, 5, 4), ("X0", "X1"))
    assert _canonical_classical_key(disc_aa) == _canonical_classical_key(disc_ab)

    # 2. Negative classical link checks (Hopf link, trefoil 3_1, Whitehead link 5_1^2, Borromean rings 6_2^3):
    hopf = MarkedDiagram.from_epd("X[1,4,2,3], X[3,2,4,1]").resolve(1)
    trefoil = MarkedDiagram.from_epd("X[1,4,2,5], X[3,6,4,1], X[5,2,6,3]").resolve(1)
    t_id3 = TangleDiagram(3, ((0, 5), (1, 4), (2, 3)), ())
    whitehead = union_two_tangles_to_classical_link(
        TangleDiagram(3, ((0, 5), (1, 4), (2, 3)), (1, -2, 1, -2, 1)), t_id3
    )
    borromean = union_two_tangles_to_classical_link(
        TangleDiagram(3, ((0, 5), (1, 4), (2, 3)), (1, -2, 1, -2, 1, -2)), t_id3
    )
    assert all(v == 0 for v in whitehead.pairwise_linking_numbers_doubled().values())
    assert all(v == 0 for v in borromean.pairwise_linking_numbers_doubled().values())
    for non_unlink in (hopf, trefoil, whitehead, borromean):
        assert certify_classical_unlink(non_unlink).is_unlink is False
        assert find_unlink_movie(non_unlink) is None

    # 3. _list_r1_plus_steps + find_unlink_movie R1+ stabilization fallback:
    cl_r3 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["8_1^-1,-1"].canonical_epd).resolve(1)
    r1p = _list_r1_plus_steps(cl_r3)
    assert len(r1p) > 0 and all(s.kind == "R1+" and s.after.num_crossings == cl_r3.num_crossings + 1 for s in r1p)

    monkeypatch.setattr(unlink_mod, "_list_r2_plus_steps", lambda _d: [])
    movie = find_unlink_movie(cl_r3, max_r3_states=0)
    assert movie is not None and len(movie) > 0
    assert any(step.kind == "R1+" for step in movie)
    assert movie[0].before == cl_r3
    for i in range(len(movie) - 1):
        assert movie[i].after == movie[i + 1].before
    assert movie[-1].after.num_crossings == 0


def test_unoriented_so4_canonical_key() -> None:
    d_rp2 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^-1"].canonical_epd)
    d_rp2_rev = d_rp2.reverse_orientation_zt()
    d_rp2_mir = d_rp2.mirror_3d()

    # R_{zt} = R_z o R_t has det = +1 in SO(4), so unoriented SO(4) key equates d and reverse_orientation_zt(d),
    # while separating chiral 4D mirrors (2_1^-1 vs its mirror with normal_euler = -2):
    assert canonical_unoriented_so4_key(d_rp2) == canonical_unoriented_so4_key(d_rp2_rev)
    assert canonical_unoriented_so4_key(d_rp2) != canonical_unoriented_so4_key(d_rp2_mir)

    # Full O(4) key equates all 3:
    assert canonical_marked_key(d_rp2, up_to_mirror_and_reverse=True) == canonical_marked_key(
        d_rp2_mir, up_to_mirror_and_reverse=True
    )

    # Check an asymmetric diagram where oriented spatial key distinguishes +F from -F = R_{zt}(F):
    d_asym = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["9_1"].canonical_epd)
    assert canonical_unoriented_so4_key(d_asym) == canonical_unoriented_so4_key(d_asym.reverse_orientation_zt())


def test_planar_2cut_detects_trivial_alexander_composite_via_a5() -> None:
    # 1-fusion ribbon 2-knot R(b^2 a^-1 b^-1 a^2) has trivial Alexander ideal E_1 = (1)
    # and cyclic solvable quotients, but |Hom(pi_1, A_5)| = 420 != 60:
    rw = compile_1fusion_ribbon_word(parse_ribbon_word("b^2 a^-1 b^-1 a^2"))
    assert tuple(compute_alexander_invariants(rw).e1_generators) == ("1",)
    comp = planar_connected_sum(rw, 0, rw, 0)
    decomp = find_planar_connected_sum_decomposition(comp)
    assert decomp is not None
    d1, d2 = decomp
    assert d1.num_vertices == rw.num_vertices and d2.num_vertices == rw.num_vertices


def test_cyclic_unknotted_tube_soundness_and_omega8() -> None:
    # Genuine cyclic-pi_1 unknotted tubes (n >= 5) return True; knotted Tube(3_1) returns False:
    t2 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^1"].canonical_epd)
    rp2 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^-1"].canonical_epd)
    tube_31 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_1^1"].canonical_epd)
    t2_rp2_stab = list_omega1_6_stabilizations(planar_connected_sum(rp2, 0, t2, 0))[0]
    rp2_4 = planar_connected_sum(planar_connected_sum(rp2, 0, rp2, 0), 0, planar_connected_sum(rp2, 0, rp2, 0), 0)
    rp2_4_tube = planar_connected_sum(rp2_4, 0, t2, 0)  # n=10, base n=8 (cp=True and red.num_vertices=8 <= 9)
    assert has_certified_cyclic_unknotted_tube(t2_rp2_stab) is True
    assert has_certified_cyclic_unknotted_tube(rp2_4_tube) is True
    assert has_certified_cyclic_unknotted_tube(tube_31) is False

    # Soundness regression [F1]: attaching a cross-neck 1-handle across tau^0(6_2) # tau^0(6_2) (ch=28,
    # det(6_2)=11 so S3/A4/F21/A5 counts match Z, yet base is composite knotted!) must NOT be certified True:
    spec_62 = ClassicalPlatSpec("6_2", "6_2", 1, 2, 6, conway_4plat_word((3, 1, 2)), (3, 1, 2))
    k62 = compile_artin_spin(spec_62)
    d0 = planar_connected_sum(k62, 0, k62, 0)
    al0 = compute_alexander_invariants(d0)
    n = d0.num_vertices
    m1 = k62.num_vertices
    vis_d = [False] * (4 * n)
    cross_neck_torus: MarkedDiagram | None = None
    for sd in range(4 * n):
        if vis_d[sd]:
            continue
        fd: list[int] = []
        cur = sd
        while not vis_d[cur]:
            vis_d[cur] = True
            fd.append(cur)
            od = d0.opp[cur]
            cur = ((od >> 2) << 2) | (((od & 3) - 1) & 3)
        da_cands = [x for x in fd if (x >> 2) < m1 and (d0.opp[x] >> 2) < m1]
        db_cands = [x for x in fd if (x >> 2) >= m1 and (d0.opp[x] >> 2) >= m1]
        for da in da_cands:
            for db in db_cands:
                oda, odb = d0.opp[da], d0.opp[db]
                u, w = n, n + 1
                nopp = list(d0.opp) + [0] * 8
                nopp[4 * u + 0], nopp[4 * w + 3] = 4 * w + 3, 4 * u + 0
                nopp[4 * u + 1], nopp[4 * w + 2] = 4 * w + 2, 4 * u + 1
                nopp[4 * u + 2], nopp[da] = da, 4 * u + 2
                nopp[4 * u + 3], nopp[oda] = oda, 4 * u + 3
                nopp[4 * w + 0], nopp[db] = db, 4 * w + 0
                nopp[4 * w + 1], nopp[odb] = odb, 4 * w + 1
                for ut, wt in (("Y", "Z"), ("Z", "Y")):
                    cand_d = MarkedDiagram(n + 2, tuple(nopp), d0.vtype + (ut, wt))
                    if not _is_planar_4reg_diagram(cand_d):
                        continue
                    st = cand_d.surface_topology()
                    if st.num_components != 1 or not st.components[0].orientable:
                        continue
                    if not (
                        certify_classical_unlink(cand_d.resolve(1)).is_unlink
                        and certify_classical_unlink(cand_d.resolve(-1)).is_unlink
                    ):
                        continue
                    if find_planar_connected_sum_decomposition(cand_d) is not None:
                        continue
                    if compute_alexander_invariants(cand_d).e1_generators != al0.e1_generators:
                        cross_neck_torus = cand_d
                        break
                if cross_neck_torus is not None:
                    break
            if cross_neck_torus is not None:
                break
        if cross_neck_torus is not None:
            break

    assert cross_neck_torus is not None
    assert has_certified_cyclic_unknotted_tube(cross_neck_torus) is False

    # Omega_8 local move admissibility and invariant preservation:
    d_o8 = MarkedDiagram.from_epd("Y[1,2,3,3], Y[4,5,6,6], X[7,8,9,1], X[7,10,11,8], X[10,5,12,11], X[12,4,2,9]")
    o8_nbrs = _list_omega8_neighbors(d_o8)
    assert len(o8_nbrs) >= 1
    for nb in o8_nbrs:
        assert _is_planar_4reg_diagram(nb)
        assert certify_classical_unlink(nb.resolve(1)).is_unlink
        assert certify_classical_unlink(nb.resolve(-1)).is_unlink
        assert nb.surface_topology() == d_o8.surface_topology()
        assert compute_alexander_invariants(nb).e1_generators == compute_alexander_invariants(d_o8).e1_generators
