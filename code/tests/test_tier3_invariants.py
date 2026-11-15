"""Unit tests for Tier-3 purely 4D invariants (`invariants/tier3.py`):
- Total mod-2 double branched cover `Sigma_2(S^4, F)` homology & Hirzebruch-Atiyah-Singer signature
- Lomonaco `pi_2(E)` & Poincare-Lefschetz regular finite cover Betti/torsion profiles
- Kamada-Oshiro octahedral symmetric quandle `(Q_6, *, rho)` colorings & 3-cocycle state sums
- 4D symmetry profiles (`F`, `-F`, `F*`) & `free_loop` / Yoshikawa move invariance
"""

from __future__ import annotations

import json
from pathlib import Path

from surfaceknots.ch11_census import disjoint_union, planar_connected_sum
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.invariants.tier3 import (
    double_branched_cover_invariants,
    lomonaco_pi2_invariants,
    octahedral_symmetric_quandle_invariants,
    symmetry_profile,
)
from surfaceknots.moves import (
    list_level_yoshikawa_neighbors,
    list_omega1_6_stabilizations,
    list_omega2_pokes,
)
from surfaceknots.yoshikawa_data import YOSHIKAWA_BY_LABEL


def _diag(label: str) -> MarkedDiagram:
    return MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL[label].canonical_epd)


def test_double_branched_cover_canonical_values() -> None:
    # 0_1: Sigma_2(S^4, S^2_std) ~= S^4 (chi=2? Wait: 2*2 - 2 = 2? No: chi(S^4)=2, 2*2 - 2 = 2! Wait: 4 - chi(S^2) = 4 - 2 = 2 == chi(S^4)!)
    s0 = double_branched_cover_invariants(_diag("0_1"))
    assert (s0.euler_characteristic, s0.b1, s0.h1_torsion, s0.b2, s0.signature) == (2, 0, (), 0, 0)
    assert s0.is_rational_homology_4sphere

    # 8_1 (spun trefoil): Q-homology 4-sphere with H_1 = Z_3
    s8 = double_branched_cover_invariants(_diag("8_1"))
    assert (s8.b1, s8.h1_torsion, s8.b2, s8.signature) == (0, (3,), 0, 0)
    assert s8.is_rational_homology_4sphere

    # 2_1^-1 = P_+ (e = +2, chi = 1): Sigma_2(S^4, P_+) ~= bar(CP^2), sigma = -1, b2+ = 0, b2- = 1
    p_plus = _diag("2_1^-1")
    sp_plus = double_branched_cover_invariants(p_plus)
    assert (sp_plus.euler_characteristic, sp_plus.b1, sp_plus.h1_torsion, sp_plus.b2) == (3, 0, (), 1)
    assert (sp_plus.signature, sp_plus.b2_plus, sp_plus.b2_minus) == (-1, 0, 1)

    # Spatial mirror P_- (e = -2): Sigma_2(S^4, P_-) ~= CP^2, sigma = +1, b2+ = 1, b2- = 0, same o4_key!
    sp_minus = double_branched_cover_invariants(p_plus.mirror_3d())
    assert (sp_minus.signature, sp_minus.b2_plus, sp_minus.b2_minus) == (1, 1, 0)
    assert sp_plus.o4_key() == sp_minus.o4_key()

    # 2_1^1 = T^2_std (chi = 0, kei_rk = 1): Sigma_2 ~= S^2 x S^2 (# copies -> b1=0, b2=2, sigma=0, b2+=1, b2-=1)
    st = double_branched_cover_invariants(_diag("2_1^1"))
    assert (st.euler_characteristic, st.b1, st.b2, st.b2_plus, st.b2_minus) == (4, 0, 2, 1, 1)

    # Check ch=11 separation between 11_2^0,-2 and 11_3^0,-2:
    ch11_path = Path(__file__).resolve().parents[1] / "data" / "expanded_ch11.json"
    ch11_entries = {e["label"]: e for e in json.loads(ch11_path.read_text())["entries"]}
    d_11_2 = MarkedDiagram.from_epd(ch11_entries["11_2^0,-2"]["epd"])
    d_11_3 = MarkedDiagram.from_epd(ch11_entries["11_3^0,-2"]["epd"])
    s_11_2 = double_branched_cover_invariants(d_11_2)
    s_11_3 = double_branched_cover_invariants(d_11_3)
    assert (s_11_2.b1, s_11_2.h1_torsion, s_11_2.b2) == (1, (), 2)
    assert (s_11_3.b1, s_11_3.h1_torsion, s_11_3.b2) == (0, (4,), 0)
    assert s_11_3.is_rational_homology_4sphere

    # Check 11_1^-2 (knotted K^2 with e = 0, chi = 0): H_1(Sigma_2) = Z_3, b2 = 2, sigma = 0
    d_11_k = MarkedDiagram.from_epd(ch11_entries["11_1^-2"]["epd"])
    s_11_k = double_branched_cover_invariants(d_11_k)
    assert (s_11_k.euler_characteristic, s_11_k.b1, s_11_k.h1_torsion, s_11_k.b2, s_11_k.signature) == (
        4,
        0,
        (3,),
        2,
        0,
    )

    # Check extremal Klein bottle P_+ # P_+ (e = +4, chi = 0): Sigma_2 ~= 2 bar(CP^2), sigma = -2!
    d_k_ext = planar_connected_sum(p_plus, 0, p_plus, 0, flip_splice=False)
    if d_k_ext.surface_topology().components[0].normal_euler != 4:
        d_k_ext = planar_connected_sum(p_plus, 0, p_plus, 0, flip_splice=True)
    s_k_ext = double_branched_cover_invariants(d_k_ext)
    assert (s_k_ext.euler_characteristic, s_k_ext.b2, s_k_ext.signature, s_k_ext.b2_plus, s_k_ext.b2_minus) == (
        4,
        2,
        -2,
        0,
        2,
    )


def test_lomonaco_pi2_and_cover_betti_canonical_values() -> None:
    # 0_1: pi_1 = Z, pi_2 = 0
    r_0 = lomonaco_pi2_invariants(_diag("0_1"))
    assert r_0.certified_pi2_z_rank == 0

    # 2_1^-1: pi_1 = Z_2, universal cover CP^2_bar \ RP^2 has b_2 = 1 -> rk_Z pi_2 = 1
    r_rp2 = lomonaco_pi2_invariants(_diag("2_1^-1"))
    assert r_rp2.certified_pi2_z_rank == 1
    assert r_rp2.to_dict()["Z2"] == ((0, (), 1, 0),)

    # 10_2 = tau^2(3_1) and 10_3 = tau^3(3_1) have infinite pi_1(E), so universal-cover
    # certified_pi2_z_rank is None while regular finite cover Betti profiles remain exact:
    r_10_2 = lomonaco_pi2_invariants(_diag("10_2"))
    assert r_10_2.certified_pi2_z_rank is None
    assert r_10_2.to_dict()["S3"] == ((1, (), 2, 2),)

    r_10_3 = lomonaco_pi2_invariants(_diag("10_3"))
    assert r_10_3.certified_pi2_z_rank is None
    assert r_10_3.to_dict()["SL2F3"] == ((1, (), 3, 3), (1, (), 7, 7))

    # Non-orientable links with finite non-abelian pi_1 (exact universal covers!):
    r_8_nn = lomonaco_pi2_invariants(_diag("8_1^-1,-1"))
    assert r_8_nn.certified_pi2_z_rank == 0
    assert r_8_nn.to_dict()["Q8"] == ((0, (), 0, 1),)

    r_10_kk = lomonaco_pi2_invariants(_diag("10_1^-2,-2"))
    assert r_10_kk.certified_pi2_z_rank == 16
    assert r_10_kk.to_dict()["Q8"] == ((0, (), 16, 1),)

    r_10_nn = lomonaco_pi2_invariants(_diag("10_1^-1,-1"))
    assert r_10_nn.certified_pi2_z_rank == 2
    assert r_10_nn.to_dict()["Dic16"] == ((0, (), 2, 3),)

    # Purely 4D cover-Betti separation of Tube(3_1) = 10_1^1 from 8_1 # 2_1^1!
    d_tube = _diag("10_1^1")
    d_stab = planar_connected_sum(_diag("8_1"), 0, _diag("2_1^1"), 0, flip_splice=False)
    if d_stab.surface_topology().yoshikawa_genera != (1,):
        d_stab = planar_connected_sum(_diag("8_1"), 0, _diag("2_1^1"), 0, flip_splice=True)
    assert lomonaco_pi2_invariants(d_tube).to_dict()["SL2F3"] == ((4, (), 54, 3),)
    assert lomonaco_pi2_invariants(d_stab).to_dict()["SL2F3"] == ((4, (), 54, 3), (4, (), 58, 7))


def test_octahedral_symmetric_quandle_and_kamada_oshiro_obstruction() -> None:
    # 2_1^-1 has pi_1 = Z_2 (order 2 meridians cannot map to order-4 elements of Q_6 subset Q_8)
    q_rp2 = octahedral_symmetric_quandle_invariants(_diag("2_1^-1"))
    assert q_rp2.col_q6 == 0
    assert not q_rp2.certifies_positive_triple_points

    # Compare 8_1^-1,-1 (t = 2, non-ribbon) vs 10_1^-2,-2 (t = 0, ribbon), which BOTH have pi_1 ~= Q_8!
    q_8 = octahedral_symmetric_quandle_invariants(_diag("8_1^-1,-1"))
    q_10 = octahedral_symmetric_quandle_invariants(_diag("10_1^-2,-2"))
    assert q_8.col_q6 == 24
    assert q_10.col_q6 == 24
    assert q_8.trivial_weight_count == 0 and q_8.nontrivial_weight_count == 24
    assert q_8.certifies_positive_triple_points
    assert q_10.trivial_weight_count == 24 and q_10.nontrivial_weight_count == 0
    assert not q_10.certifies_positive_triple_points


def test_symmetry_profile_and_tier3_yoshikawa_move_invariance() -> None:
    # 8_1 = tau^0(3_1) is constructively amphicheiral in R^4 even though classical 3_1 is chiral!
    sym_8_1 = symmetry_profile(_diag("8_1"))
    assert sym_8_1.constructive_so4_amphicheiral
    assert not sym_8_1.chiral_by_normal_euler
    assert not sym_8_1.oriented_chiral_by_cocycle

    # 10_2 = tau^2(3_1) is minus-amphicheiral (+F* ~= -F != +F): unoriented SO(4)-amphicheiral,
    # while non-invertible and oriented-chiral by Mochizuki's 3-cocycle!
    sym_10_2 = symmetry_profile(_diag("10_2"))
    assert sym_10_2.oriented_chiral_by_cocycle
    assert sym_10_2.non_invertible_by_cocycle
    assert sym_10_2.constructive_so4_amphicheiral
    assert not sym_10_2.constructive_plus_amphicheiral
    assert sym_10_2.constructive_minus_amphicheiral
    assert (
        sym_10_2.symmetry_summary
        == "Minus-amphicheiral & non-invertible (+F* ~= -F != +F; unoriented SO(4)-amphicheiral)"
    )

    # Disjoint union (+10_2) u (+10_2)* is globally invertible (reversing both components swaps
    # their cocycle states up to permutation) and unoriented SO(4)-amphicheiral:
    d_pair_10_2 = disjoint_union(_diag("10_2"), _diag("10_2").mirror_3d())
    sym_pair_10_2 = symmetry_profile(d_pair_10_2)
    assert not sym_pair_10_2.non_invertible_by_cocycle
    assert not sym_pair_10_2.oriented_chiral_by_cocycle
    assert sym_pair_10_2.constructive_so4_amphicheiral

    # Paired (genus, normal_euler) check for c >= 2 so e_list sum == 0 never decouples components:
    p_plus = _diag("2_1^-1")
    k_minus4 = planar_connected_sum(
        p_plus.mirror_3d(), 0, p_plus.mirror_3d(), 0, flip_splice=False
    )
    if k_minus4.surface_topology().components[0].normal_euler != -4:
        k_minus4 = planar_connected_sum(
            p_plus.mirror_3d(), 0, p_plus.mirror_3d(), 0, flip_splice=True
        )
    d_3link = disjoint_union(disjoint_union(p_plus, p_plus), k_minus4)
    sym_3link = symmetry_profile(d_3link, max_bfs_states=0)
    assert sym_3link.chiral_by_normal_euler
    assert sym_3link.invariant_flags() == symmetry_profile(
        d_3link.mirror_3d(), max_bfs_states=0
    ).invariant_flags()

    # 10_3 = tau^3(3_1) is constructively unoriented-amphicheiral in R^4
    sym_10_3 = symmetry_profile(_diag("10_3"))
    assert sym_10_3.constructive_so4_amphicheiral

    # 2_1^-1 is chiral by normal Euler (+2 != -2)
    sym_rp2 = symmetry_profile(_diag("2_1^-1"))
    assert sym_rp2.chiral_by_normal_euler

    # Verify strict Yoshikawa move + stabilization + free_loop invariance of ALL Tier-3 invariants!
    for lbl in ("8_1", "10_2", "8_1^-1,-1", "7_1^0,-2"):
        base = _diag(lbl)
        s_base = double_branched_cover_invariants(base).o4_key()
        p_base = lomonaco_pi2_invariants(base)
        q_base = octahedral_symmetric_quandle_invariants(base)
        f_base = symmetry_profile(base, max_bfs_states=0).invariant_flags()

        neighbors = (
            list_level_yoshikawa_neighbors(base)[:4]
            + list_omega1_6_stabilizations(base)[:3]
            + list_omega2_pokes(base)[:3]
            + [base.mirror_3d(), base.mirror_time(), base.reverse_orientation_zt()]
        )
        for nb in neighbors:
            assert double_branched_cover_invariants(nb).o4_key() == s_base
            assert lomonaco_pi2_invariants(nb) == p_base
            assert octahedral_symmetric_quandle_invariants(nb) == q_base
            assert symmetry_profile(nb, max_bfs_states=0).invariant_flags() == f_base

    # Also test free_loop stabilization invariance (`O[1]` vs `Omega_1+(O[1])` disjointly unioned with `8_1`)
    u_free = MarkedDiagram.from_epd("O[1]")
    u_kink = MarkedDiagram(num_vertices=1, opp=(1, 0, 3, 2), vtype=("X0",), free_loops=0)
    d_a = disjoint_union(_diag("8_1"), u_free)
    d_b = disjoint_union(_diag("8_1"), u_kink)
    assert double_branched_cover_invariants(d_a).o4_key() == double_branched_cover_invariants(d_b).o4_key()
    assert lomonaco_pi2_invariants(d_a) == lomonaco_pi2_invariants(d_b)
    assert octahedral_symmetric_quandle_invariants(d_a) == octahedral_symmetric_quandle_invariants(d_b)


def test_lomonaco_pi2_multistep_and_commutator_ribbon_soundness() -> None:
    from surfaceknots.families import compile_1fusion_ribbon_word, parse_ribbon_word

    def _permute_vertices(d: MarkedDiagram, shift: int) -> MarkedDiagram:
        n = d.num_vertices
        perm = [(i + shift) % n for i in range(n)]
        inv_p = [0] * n
        for ov, nv in enumerate(perm):
            inv_p[nv] = ov
        new_opp = [0] * (4 * n)
        new_vt = [d.vtype[inv_p[nv]] for nv in range(n)]
        for nv in range(n):
            ov = inv_p[nv]
            for s in range(4):
                od = d.opp[4 * ov + s]
                new_opp[4 * nv + s] = 4 * perm[od >> 2] + (od & 3)
        return MarkedDiagram(n, tuple(new_opp), tuple(new_vt), d.free_loops)

    # 1. 10_2 # R(a^4 b a^-4 b^-1) has identical S_3 regular cover Betti to 10_2,
    #    yet pi_2(E) has infinite Z-rank: certified_pi2_z_rank MUST remain None across relabelings!
    d_10_2 = _diag("10_2")
    d_rib = compile_1fusion_ribbon_word(parse_ribbon_word("a^4 b a^-4 b^-1"))
    d_sum = planar_connected_sum(d_10_2, 0, d_rib, 0, flip_splice=False)
    lom_base = lomonaco_pi2_invariants(d_sum, include_sl2f3=False)
    assert lom_base.certified_pi2_z_rank is None
    for sh in (5, 17, 35, 42):
        lom_sh = lomonaco_pi2_invariants(_permute_vertices(d_sum, sh), include_sl2f3=False)
        assert lom_sh == lom_base
        assert lom_sh.certified_pi2_z_rank is None

    # 2. Multi-step Yoshikawa walk on 10_2 (where meridian elimination can introduce higher
    #    powers of mu) preserves lomonaco_pi2_invariants identically:
    p_10_2 = lomonaco_pi2_invariants(d_10_2)
    cur = d_10_2
    for step_fn in (list_omega1_6_stabilizations, list_omega2_pokes, list_level_yoshikawa_neighbors):
        nbs = step_fn(cur)
        if nbs:
            cur = nbs[0]
            assert lomonaco_pi2_invariants(cur) == p_10_2


def test_r1_plus_movie_step_invariance() -> None:
    from surfaceknots.invariants.linking import _transport_component_tags_forward
    from surfaceknots.invariants.quandle import _evaluate_classical_movie_cocycle
    from surfaceknots.invariants.tier3 import _enumerate_q6_graph_colorings, _eval_q6_movie_bitmask
    from surfaceknots.unlink import UnlinkMovieStep, _list_r1_plus_steps, _list_r3_steps, find_unlink_movie

    d8 = _diag("8_1")
    lp8 = d8.resolve(1)
    base_mov = find_unlink_movie(lp8)
    assert base_mov is not None

    r1_steps = _list_r1_plus_steps(lp8)
    assert r1_steps
    st0 = r1_steps[0]
    step_r1_minus = UnlinkMovieStep(
        kind="R1-",
        before=st0.after,
        after=lp8,
        meta={"old_to_new": {v: v for v in range(lp8.num_crossings)}},
    )
    aug_cancel = [st0, step_r1_minus] + base_mov

    aug_r3 = None
    for st_cand in r1_steps:
        r3_cands = _list_r3_steps(st_cand.after)
        if r3_cands:
            tail = find_unlink_movie(r3_cands[0].after)
            if tail is not None:
                aug_r3 = [st_cand, r3_cands[0]] + tail
                break
    assert aug_r3 is not None

    n_c = lp8.num_crossings
    cycs, _ = lp8.circle_components()
    ddir = [0] * (4 * n_c)
    for cyc in cycs:
        for hd in cyc:
            ddir[hd] = 1
            ddir[hd ^ 2] = -1

    q6_cols = _enumerate_q6_graph_colorings(d8)
    cl_vs = [v for v in range(d8.num_vertices) if d8.vtype[v] in ("X0", "X1")]

    for mov in (aug_cancel, aug_r3):
        tags = _transport_component_tags_forward(lp8, mov)
        assert tags[-1] == ()
        for c_const in range(3):
            assert _evaluate_classical_movie_cocycle(
                mov, [c_const] * (4 * n_c), ddir, 3
            ) == _evaluate_classical_movie_cocycle(base_mov, [c_const] * (4 * n_c), ddir, 3)
        for qc in q6_cols[:6]:
            ic = [qc[4 * vo + s] for vo in cl_vs for s in range(4)]
            assert _eval_q6_movie_bitmask(mov, ic) == _eval_q6_movie_bitmask(base_mov, ic)

