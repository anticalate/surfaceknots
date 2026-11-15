"""Audit-remediation soundness tests, blind Yoshikawa census recovery, Tietze/surgery/Dlk
regressions, and pure-Python verification of the ch = 11 / ch = 12 headline theorems
(including the two Cappell-Shaneson 2-knots CS_A0 and CS_A1 at ch = 12, m = 2).
"""

from __future__ import annotations

import json
from pathlib import Path
import time

import pytest
from sympy.combinatorics.free_groups import free_group
from sympy.combinatorics.fp_groups import FpGroup

from surfaceknots.ch11_census import build_reference_catalog, planar_connected_sum
from surfaceknots.complexity import (
    certify_klein_bottle_subadditivity,
    construct_minimal_unknot_diagram,
)
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.families import CLASSICAL_PLAT_BY_LABEL, compile_1fusion_ribbon_word, compile_satoh_tube, parse_ribbon_word
from surfaceknots.fast_enum import (
    _filter_cyclic_tubes_parallel,
    cascade_tunnel_diagrams_c,
    enumerate_targeted_census_level_candidates,
)
from surfaceknots.invariants.alexander import compute_alexander_invariants
from surfaceknots.invariants.fundamental_group import GroupPresentation, tietze_simplify, wirtinger_presentation
from surfaceknots.invariants.linking import double_linking_mod2
from surfaceknots.invariants.tier3 import lomonaco_pi2_invariants
from surfaceknots.move_search import _indep_is_single_local_move, verify_move_path_independent
from surfaceknots.moves import (
    _list_marked_tube_flip_neighbors,
    _list_omega3_and_4_neighbors,
    _list_omega5_neighbors,
    _list_omega7_neighbors,
    _list_omega8_neighbors,
    list_4d_tait_flypes,
    list_omega1_6_stabilizations,
    list_omega2_pokes,
)
from surfaceknots.unlink import certify_classical_unlink
from surfaceknots.yoshikawa_data import (
    YOSHIKAWA_BY_LABEL,
    YOSHIKAWA_TABLE,
    check_e1_ideal_against_reference,
    compute_full_invariant_record,
    invariant_discrimination_key,
)


DATA_DIR = Path(__file__).resolve().parents[1] / "data"

CS_A0_EPD = (
    "X[1,2,3,4], X[5,6,4,7], X[6,8,9,1], Y[9,10,11,2], "
    "X[12,13,14,3], Y[14,15,16,7], X[5,17,18,8], X[17,19,10,18], "
    "X[11,20,21,12], X[13,21,22,15], X[23,24,19,16], X[24,23,22,20]"
)
CS_A1_EPD = (
    "X[1,2,3,4], X[4,5,6,1], X[7,8,9,6], Y[9,10,11,2], "
    "X[3,12,13,14], Y[14,15,16,17], X[18,7,5,17], X[8,18,19,10], "
    "X[20,21,12,11], X[21,22,15,13], X[23,24,19,16], X[24,23,22,20]"
)


def _blind_prime_keys_at_ch(ch: int, num_workers: int = 4) -> set[tuple]:
    diags, _ = enumerate_targeted_census_level_candidates(
        ch=ch,
        m_min=1,
        m_max=ch - 1,
        max_bfs_states=1024,
        num_workers=num_workers,
    )
    tun = cascade_tunnel_diagrams_c(
        ch=ch,
        diagrams=diags,
        depth_plus1=192,
        depth_plus2=128,
        max_bfs_states=2048,
        num_workers=num_workers,
    )
    surv = _filter_cyclic_tubes_parallel(tun, num_workers=num_workers)
    found_keys = set()
    for d in surv:
        rec = compute_full_invariant_record(d, label=f"blind_{ch}")
        found_keys.add(invariant_discrimination_key(rec))
    return found_keys


def test_blind_ch_le_9_recovers_yoshikawa_prime_table() -> None:
    """Blind enumeration across ch = 6..9 (m = 1 .. ch - 1) recovers 100% of Yoshikawa's prime
    classes at each crossing-and-saddle index with zero spurious non-reference prime classes.
    """
    ref_cat = build_reference_catalog(max_composite_ch=9)
    for ch in (6, 7, 8, 9):
        expected_entries = [e for e in YOSHIKAWA_TABLE if e.ch_index == ch]
        expected_keys = {
            invariant_discrimination_key(
                compute_full_invariant_record(MarkedDiagram.from_epd(e.canonical_epd), label=e.label)
            )
            for e in expected_entries
        }
        blind_keys = _blind_prime_keys_at_ch(ch, num_workers=4)
        assert expected_keys.issubset(blind_keys), f"Missing Yoshikawa prime class at ch={ch}"

        # Any additional blind survivor must match a strictly smaller prime, unknot, split, or composite key
        lower_or_comp_keys = {
            k for k, tags in ref_cat.items() if any(t not in {e.label for e in expected_entries} for t in tags)
        }
        unexpected = blind_keys - expected_keys - lower_or_comp_keys
        assert not unexpected, f"Spurious unknown invariant class at ch={ch}: {len(unexpected)}"


@pytest.mark.slow
def test_blind_ch_10_recovers_all_twelve_yoshikawa_classes() -> None:
    expected_entries = [e for e in YOSHIKAWA_TABLE if e.ch_index == 10]
    assert len(expected_entries) == 12
    expected_keys = {
        invariant_discrimination_key(
            compute_full_invariant_record(MarkedDiagram.from_epd(e.canonical_epd), label=e.label)
        )
        for e in expected_entries
    }
    blind_keys = _blind_prime_keys_at_ch(10, num_workers=8)
    assert expected_keys.issubset(blind_keys)


def test_groebner_e1_meridian_pinning_and_tietze_invariance() -> None:
    """Verify two-way strong Groebner Z-ideal check on all 23 Yoshikawa entries, fast simplify=False
    runtime, and strict peripheral-preserving Tietze invariance on mixed-orientability links.
    """
    strict_matches = 0
    errata_labels = []
    for entry in YOSHIKAWA_TABLE:
        d = MarkedDiagram.from_epd(entry.canonical_epd)
        ok_print, _ = check_e1_ideal_against_reference(d, entry.yoshikawa_printed_e1)
        ok_corr, _ = check_e1_ideal_against_reference(d, entry.corrected_e1)
        assert ok_corr is True, f"Corrected E_1 failed two-way Groebner equality on {entry.label}"
        if ok_print:
            strict_matches += 1
        else:
            errata_labels.append(entry.label)

    assert strict_matches == 21
    assert errata_labels == ["9_1^1,-2", "10_2^0,-2"]

    ch11 = {e["label"]: e for e in json.loads((DATA_DIR / "expanded_ch11.json").read_text())["entries"]}
    d_11_2 = MarkedDiagram.from_epd(ch11["11_2^0,-2"]["epd"])
    d_11_3 = MarkedDiagram.from_epd(ch11["11_3^0,-2"]["epd"])

    t0 = time.perf_counter()
    al2_raw = compute_alexander_invariants(d_11_2)
    al3_raw = compute_alexander_invariants(d_11_3)
    elapsed = time.perf_counter() - t0
    assert elapsed < 1.5, f"Strong Groebner Alexander evaluation took too long: {elapsed:.2f}s"

    assert check_e1_ideal_against_reference(d_11_2, ("x**2 + y", "(x + 1)*(y - 1)"))[0] is True
    assert check_e1_ideal_against_reference(d_11_3, ("x**2 + y", "(x - 1)*(y - 1)"))[0] is True
    assert check_e1_ideal_against_reference(d_11_2, ("x**2 + y", "(x - 1)*(y - 1)"))[0] is False
    assert al2_raw.evaluation_fingerprint != al3_raw.evaluation_fingerprint

    # Adversarial peripheral-preserving Nielsen substitution g_2 -> g_1 * g_2
    # that previously swapped E_1(11_2^{0,-2}) and E_1(11_3^{0,-2}) before meridian pinning:
    for d_link, expected_fp in ((d_11_2, al2_raw.evaluation_fingerprint), (d_11_3, al3_raw.evaluation_fingerprint)):
        pres = tietze_simplify(wirtinger_presentation(d_link))
        assert pres.num_generators == 2

        def _sub_word(w: tuple[int, ...]) -> tuple[int, ...]:
            out: list[int] = []
            for lit in w:
                seq = (1, 2) if lit == 2 else ((-2, -1) if lit == -2 else (lit,))
                for x in seq:
                    if out and out[-1] == -x:
                        out.pop()
                    else:
                        out.append(x)
            return tuple(out)

        from dataclasses import replace

        adv_pres = replace(
            pres,
            relators=tuple(_sub_word(r) for r in pres.relators),
            component_meridians=tuple(_sub_word(m) for m in pres.component_meridians),
        )
        al_adv = compute_alexander_invariants(d_link, pres=adv_pres)
        assert al_adv.evaluation_fingerprint == expected_fp


def test_independent_move_verifier_rejects_surgery_and_accepts_genuine_moves() -> None:
    """Verify surgical 2-vertex splice Tube(3_1) -> 8_1 # T^2 is rejected while genuine local moves pass."""
    d_tube = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_1^1"].canonical_epd)
    d_8_1 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["8_1"].canonical_epd)
    d_t2 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^1"].canonical_epd)
    d_sum = planar_connected_sum(d_8_1, 0, d_t2, 0, flip_splice=False)
    if d_sum.surface_topology().yoshikawa_genera != (1,):
        d_sum = planar_connected_sum(d_8_1, 0, d_t2, 0, flip_splice=True)

    # Also test the exact 2-marked-vertex surgical flip turning 10_1^1 into a composite-group diagram:
    epd_tube31 = YOSHIKAWA_BY_LABEL["10_1^1"].canonical_epd
    assert _indep_is_single_local_move(d_tube, d_sum, allow_o4=True) is None
    res_false = verify_move_path_independent([epd_tube31, d_sum.to_epd()], allow_o4=True)
    assert res_false.is_valid is False

    # Construct surgical 2-crossing flip on Tube(3_1) preserving opp, ch=10, T^2 topology, and admissibility:
    vtype_flipped = list(d_tube.vtype)
    vtype_flipped[2] = "X1"
    vtype_flipped[5] = "X1"
    d_surg = MarkedDiagram(num_vertices=d_tube.num_vertices, opp=d_tube.opp, vtype=tuple(vtype_flipped))
    assert certify_classical_unlink(d_surg.resolve(1)).is_unlink
    assert certify_classical_unlink(d_surg.resolve(-1)).is_unlink
    assert d_surg.surface_topology().yoshikawa_genera == d_tube.surface_topology().yoshikawa_genera
    assert invariant_discrimination_key(compute_full_invariant_record(d_surg, label="surg")) != (
        invariant_discrimination_key(compute_full_invariant_record(d_tube, label="tube"))
    )
    assert _indep_is_single_local_move(d_tube, d_surg, allow_o4=True) is None
    assert verify_move_path_independent([d_tube.to_epd(), d_surg.to_epd()], allow_o4=True).is_valid is False

    # Exact 3-vertex adversarial surgical pair (Tube(3_1) -> 8_1 # T^2) from third-party audit:
    rest_adv = (
        "X[6,7,8,9], X[9,10,11,2], X[12,13,14,3], X[15,16,5,14], "
        "Z[16,17,18,7], X[18,19,10,8], Y[11,19,20,12], X[13,20,17,15]"
    )
    epd_adv_a = f"X[2,21,21,22], Z[1,22,3,4], Y[4,5,6,1], {rest_adv}"
    epd_adv_b = f"Y[2,21,22,23], Z[23,22,21,3], X[5,24,24,6], {rest_adv}"
    d_adv_a = MarkedDiagram.from_epd(epd_adv_a)
    d_adv_b = MarkedDiagram.from_epd(epd_adv_b)
    assert _indep_is_single_local_move(d_adv_a, d_adv_b, allow_o4=True) is None
    assert verify_move_path_independent([epd_adv_a, epd_adv_b], allow_o4=True).is_valid is False

    # Positive validation on Klein-bottle subadditivity certificate:
    kb_cert = certify_klein_bottle_subadditivity()
    res_kb = verify_move_path_independent([str(x) for x in kb_cert["move_path_epds"]], allow_o4=True)  # type: ignore[union-attr]
    assert res_kb.is_valid is True

    # Positive validation across genuine move families (Omega1/6, Omega2, Omega3/4, Omega5, Omega7, Omega8, flype, tube-flip):
    for nb in (
        list_omega1_6_stabilizations(d_8_1)[:2]
        + list_omega2_pokes(d_8_1)[:2]
        + _list_omega3_and_4_neighbors(d_8_1)[:2]
        + _list_omega5_neighbors(d_8_1)[:2]
        + list_4d_tait_flypes(d_8_1)[:2]
    ):
        assert _indep_is_single_local_move(d_8_1, nb, allow_o4=True) is not None
    for nb in _list_marked_tube_flip_neighbors(d_tube)[:2]:
        assert _indep_is_single_local_move(d_tube, nb, allow_o4=True) is not None

    # Deterministic short stabilization walk to exercise genuine Omega_7 and Omega_8 moves:
    import random

    rng = random.Random(0)
    pool = [MarkedDiagram.from_epd(e.canonical_epd) for e in YOSHIKAWA_TABLE if e.canonical_epd and e.ch_index >= 6]
    walkers = [
        list_omega1_6_stabilizations,
        list_omega2_pokes,
        _list_omega3_and_4_neighbors,
        _list_omega5_neighbors,
        _list_omega7_neighbors,
    ]
    seen_o7 = 0
    seen_o8 = 0
    for _ in range(120):
        if seen_o7 >= 3 and seen_o8 >= 3:
            break
        cur = rng.choice(pool)
        n0 = cur.num_vertices
        for _step in range(10):
            if seen_o7 < 3:
                o7_nbs = _list_omega7_neighbors(cur)
                if o7_nbs:
                    nb7 = rng.choice(o7_nbs)
                    assert _indep_is_single_local_move(cur, nb7, allow_o4=True) is not None
                    seen_o7 += 1
            if seen_o8 < 3:
                o8_nbs = _list_omega8_neighbors(cur)
                if o8_nbs:
                    nb8 = rng.choice(o8_nbs)
                    assert _indep_is_single_local_move(cur, nb8, allow_o4=True) is not None
                    seen_o8 += 1
            w_fn = rng.choice(walkers)
            if cur.num_vertices >= n0 + 3 and w_fn in walkers[:2]:
                w_fn = walkers[2]
            w_nbs = w_fn(cur)
            if w_nbs:
                cur = rng.choice(w_nbs)
    assert seen_o7 >= 1 and seen_o8 >= 1, f"Expected genuine Omega_7 and Omega_8 instances (got {seen_o7}, {seen_o8})"


def test_sanderson_dlk_regression_and_high_marker_invariance() -> None:
    """Verify Dlk counterexample regression (0_1 u K^2(e=0), 9_1^{1,-2}) and m > 12 stabilization invariance."""
    from surfaceknots.ch11_census import disjoint_union

    d_s2 = construct_minimal_unknot_diagram(0, 0)
    d_k2 = construct_minimal_unknot_diagram(-2, 0)
    d_split = disjoint_union(d_s2, d_k2)
    assert len(double_linking_mod2(d_split)) == 0

    # lambda_{T^2, S^2} == 0 on 9_1^{0,1} -> returns (); lambda_{T^2, K^2} != 0 on 9_1^{1,-2} -> invariant 0!
    d_9_01 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["9_1^0,1"].canonical_epd)
    assert len(double_linking_mod2(d_9_01)) == 0
    d_9_err = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["9_1^1,-2"].canonical_epd)
    assert [row[2] for row in double_linking_mod2(d_9_err)] == [0]

    for lbl, expected_vals in (("8_1^1,1", [0, 0]), ("10_1^1,1", [1, 1])):
        cur = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL[lbl].canonical_epd)
        assert [row[2] for row in double_linking_mod2(cur)] == expected_vals
        # Stabilize with marked saddles (Omega_6^+) until num_marked >= 14 (> 12 threshold):
        while cur.num_marked < 14:
            stabs = [nb for nb in list_omega1_6_stabilizations(cur) if nb.num_marked > cur.num_marked]
            assert stabs, f"Expected Omega_6^+ stabilization on {lbl}"
            cur = stabs[0]
        assert cur.num_marked >= 14
        assert [row[2] for row in double_linking_mod2(cur)] == expected_vals


def test_ch11_packet_group_isomorphism_and_2O_klein_bottle() -> None:
    """Verify the explicit isomorphism between R(b a b^-1 a^-1) and R(b a b^-1 a) in 11_1.{a,b,c},
    and pi_1(S^4 \\ 11_1^-2) ~= 2O (binary octahedral group of order 48, meridian order 8).
    """
    # 1. Explicit mutually inverse group isomorphisms between
    #    G_b = <a, b | w_b^-1 a w_b b^-1 = 1>, w_b = b a b^-1 a^-1  (11_1.b = R(b a b^-1 a^-1))
    #    G_c = <a, b | w_c^-1 a w_c b^-1 = 1>, w_c = b a b^-1 a      (11_1.c = R(b a b^-1 a))
    #    given by:
    #      phi(a) = a^-1, phi(b) = (b^-1 a) b^-1 (a^-1 b)
    #      psi(a) = a^-1, psi(b) = (a^-1 b a) b^-1 (a^-1 b^-1 a)
    F_ab, a, b = free_group("a, b")
    wb = b * a * b**-1 * a**-1
    wc = b * a * b**-1 * a
    rb = wb**-1 * a * wb * b**-1
    rc = wc**-1 * a * wc * b**-1

    def _sub(word, im_a, im_b):
        res = F_ab.identity
        a_sym = a.array_form[0][0]
        for sym, exp in word.array_form:
            res = res * ((im_a if sym == a_sym else im_b) ** exp)
        return res

    phi_a = a**-1
    phi_b = (b**-1 * a) * b**-1 * (a**-1 * b)
    psi_a = a**-1
    psi_b = (a**-1 * b * a) * b**-1 * (a**-1 * b**-1 * a)

    # Pure free-group identities: psi(phi(x)) == x and phi(psi(x)) == x ALREADY in F(a, b) or modulo relators!
    def _in_normal_closure(word, rel) -> bool:
        """BFS up to 3 conjugate relator insertions in F(a, b) to prove word == 1 rigorously."""
        if word == F_ab.identity:
            return True
        gens_pm = [a, a**-1, b, b**-1]
        conjs = [rel, rel**-1]
        for g1 in gens_pm:
            conjs.extend([g1 * rel * g1**-1, g1 * rel**-1 * g1**-1])
            for g2 in gens_pm:
                if g1 * g2 != F_ab.identity:
                    c = g1 * g2
                    conjs.extend([c * rel * c**-1, c * rel**-1 * c**-1])
        # Also allow inserting rel^{+-1} at any position of the current word (1 step = 1 relator rewrite):
        frontier = {word}
        visited = {word}
        for _ in range(3):
            nxt = set()
            for w in frontier:
                letters = []
                for sym, exp in w.array_form:
                    g = (a if sym == a.array_form[0][0] else b)
                    letters.extend([g if exp > 0 else g**-1] * abs(exp))
                for pos in range(len(letters) + 1):
                    pre = F_ab.identity
                    for x in letters[:pos]:
                        pre = pre * x
                    suf = F_ab.identity
                    for x in letters[pos:]:
                        suf = suf * x
                    for r_pm in (rel, rel**-1):
                        cand = pre * r_pm * suf
                        if cand == F_ab.identity:
                            return True
                        if len(cand) <= len(word) + 4 and cand not in visited:
                            visited.add(cand)
                            nxt.add(cand)
            frontier = nxt
        return False

    assert _in_normal_closure(_sub(rb, phi_a, phi_b), rc)
    assert _in_normal_closure(_sub(rc, psi_a, psi_b), rb)
    assert _in_normal_closure(_sub(phi_a, psi_a, psi_b) * a**-1, rb)
    assert _in_normal_closure(_sub(phi_b, psi_a, psi_b) * b**-1, rb)
    assert _in_normal_closure(_sub(psi_a, phi_a, phi_b) * a**-1, rc)
    assert _in_normal_closure(_sub(psi_b, phi_a, phi_b) * b**-1, rc)

    # 2. Verify pi_1(S^4 \ 11_1^{-2}) ~= 2O (order 48, unique involution, meridian order 8):
    ch11 = {e["label"]: e for e in json.loads((DATA_DIR / "expanded_ch11.json").read_text())["entries"]}
    d_k = MarkedDiagram.from_epd(ch11["11_1^-2"]["epd"])
    pres_k = tietze_simplify(wirtinger_presentation(d_k))
    syms = ",".join(f"g{i}" for i in range(1, pres_k.num_generators + 1))
    fg_tuple = free_group(syms)
    F_k, gens_k = fg_tuple[0], fg_tuple[1:]

    def _to_fw(rel: tuple[int, ...]):
        w = F_k.identity
        for lit in rel:
            g = gens_k[abs(lit) - 1]
            w = w * (g if lit > 0 else g**-1)
        return w

    fp_k = FpGroup(F_k, [_to_fw(r) for r in pres_k.relators])
    assert fp_k.order() == 48
    perm_res = fp_k._to_perm_group()
    perm_g = perm_res[0] if isinstance(perm_res, tuple) else perm_res
    assert perm_g.order() == 48
    # Unique involution distinguishes 2O from GL(2, 3):
    elems = list(perm_g.generate())
    involutions = [x for x in elems if x.order() == 2]
    assert len(involutions) == 1
    # Meridian has order 8 (> 2 => no standard RP^2 or unknotted K^2 summand):
    m_lit = pres_k.component_meridians[0][0]
    m_perm = perm_g.generators[abs(m_lit) - 1]
    assert m_perm.order() == 8


def test_ch12_two_cappell_shaneson_2knots_and_cover_separation() -> None:
    """Verify admissibility, cubic unit-norm Alexander polynomials, and regular F21 / F57 cover-homology
    separation of both ch = 12, m = 2 Cappell-Shaneson 2-knots (CS_A0 and CS_A1) from their 1-fusion
    ribbon Alexander-twins (11_2 = R(b^-1 a b a) and R(b^-1 a b a^-1 b a)).
    """
    d_cs0 = MarkedDiagram.from_epd(CS_A0_EPD)
    d_cs1 = MarkedDiagram.from_epd(CS_A1_EPD)

    for d_cs in (d_cs0, d_cs1):
        assert d_cs.ch_index == 12
        assert d_cs.num_marked == 2
        assert certify_classical_unlink(d_cs.resolve(1)).is_unlink
        assert certify_classical_unlink(d_cs.resolve(-1)).is_unlink
        topo = d_cs.surface_topology()
        assert topo.num_components == 1
        assert topo.components[0].homeomorphism_type == "S^2"

    # CS_A0 vs ribbon 11_2 = R(b^-1 a b a): separated by 2-primary second-commutator torsion in F21 cover!
    ch11 = {e["label"]: e for e in json.loads((DATA_DIR / "expanded_ch11.json").read_text())["entries"]}
    d_rib0 = MarkedDiagram.from_epd(ch11["11_2"]["epd"])
    al_cs0 = compute_alexander_invariants(d_cs0).single_var_alex
    al_rib0 = compute_alexander_invariants(d_rib0).single_var_alex
    assert al_cs0 in ("t**3 - 2*t**2 + t - 1", "t**3 - t**2 + 2*t - 1")
    assert al_rib0 in ("t**3 - 2*t**2 + t - 1", "t**3 - t**2 + 2*t - 1")

    cov_cs0 = lomonaco_pi2_invariants(d_cs0).to_dict()
    cov_rib0 = lomonaco_pi2_invariants(d_rib0).to_dict()
    assert cov_cs0["F21"] == ((1, (7,), 6, 6),)
    assert cov_rib0["F21"] == ((1, (2, 2, 2, 2, 2, 14), 6, 6),)

    rec_cs0 = compute_full_invariant_record(d_cs0, label="CS_A0")
    rec_rib0 = compute_full_invariant_record(d_rib0, label="11_2")
    assert invariant_discrimination_key(rec_cs0) != invariant_discrimination_key(rec_rib0)

    # CS_A1 vs ribbon R(b^-1 a b a^-1 b a): separated by 11-primary second-commutator torsion in F57 cover!
    d_rib1 = compile_1fusion_ribbon_word(parse_ribbon_word("b^-1 a b a^-1 b a"))
    al_cs1 = compute_alexander_invariants(d_cs1).single_var_alex
    al_rib1 = compute_alexander_invariants(d_rib1).single_var_alex
    assert al_cs1 in ("t**3 - 2*t**2 + 3*t - 1", "t**3 - 3*t**2 + 2*t - 1")
    assert al_rib1 in ("t**3 - 2*t**2 + 3*t - 1", "t**3 - 3*t**2 + 2*t - 1")

    cov_cs1 = lomonaco_pi2_invariants(d_cs1).to_dict()
    cov_rib1 = lomonaco_pi2_invariants(d_rib1).to_dict()
    assert cov_cs1["F57"] == ((1, (19,), 18, 18),)
    assert cov_rib1["F57"] == ((1, (11, 11, 209), 18, 18),)

    rec_cs1 = compute_full_invariant_record(d_cs1, label="CS_A1")
    rec_rib1 = compute_full_invariant_record(d_rib1, label="R_A1")
    assert invariant_discrimination_key(rec_cs1) != invariant_discrimination_key(rec_rib1)

    # If expanded_ch12.json is present, verify both CS_A0 and CS_A1 appear as certified EXACT_SINGLETON classes!
    ch12_path = DATA_DIR / "expanded_ch12.json"
    if ch12_path.exists():
        ch12_payload = json.loads(ch12_path.read_text())
        ch12_by_key = {invariant_discrimination_key(e): e for e in ch12_payload["entries"]}
        k0 = invariant_discrimination_key(rec_cs0)
        k1 = invariant_discrimination_key(rec_cs1)
        assert k0 in ch12_by_key, "CS_A0 must appear as a new certified ch=12 prime class!"
        assert k1 in ch12_by_key, "CS_A1 must appear as a new certified ch=12 prime class!"
        assert ch12_by_key[k0]["cell_status"] == "EXACT_SINGLETON"
        assert ch12_by_key[k1]["cell_status"] == "EXACT_SINGLETON"


def test_rp2_cert_worker_and_tunnel_discard_regressions() -> None:
    """Regression test for [H-1] `_rp2_cert_worker` signature/status and [H-2] `tunnel_diagrams_c` discard check."""
    from surfaceknots.complexity import MINIMAL_KLEIN_BOTTLE_E0_EPD
    from surfaceknots.families import compile_mixed_plat_link, simplify_compiled_diagram
    from surfaceknots.fast_enum import tunnel_diagrams_c
    from surfaceknots.move_search import are_prime_diagrams_tunnel_merged_c
    from surfaceknots.targeted_census import _rp2_cert_worker

    # [H-1] _rp2_cert_worker must run cleanly without TypeError and distinguish reducible vs irreducible statuses:
    p2_epd = YOSHIKAWA_BY_LABEL["2_1^-1"].canonical_epd
    res_p2 = _rp2_cert_worker(p2_epd)
    cert_p2 = res_p2["reducible_certificate"]
    assert isinstance(cert_p2, dict)
    assert cert_p2["status"] == "STANDARD_P2"
    assert cert_p2["path_verified"] is True

    d_8_1 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["8_1"].canonical_epd)
    d_p2 = MarkedDiagram.from_epd(p2_epd)
    d_8_1_p2 = planar_connected_sum(d_8_1, 0, d_p2, 0, flip_splice=False)
    res_split = _rp2_cert_worker(d_8_1_p2.to_epd())
    cert_split = res_split["reducible_certificate"]
    assert isinstance(cert_split, dict)
    assert cert_split["status"] == "SPLIT_CONNECTED_SUM"

    res_irred = _rp2_cert_worker(d_8_1.to_epd())
    cert_irred = res_irred["reducible_certificate"]
    assert isinstance(cert_irred, dict)
    assert res_irred["meridian_order_lower_bound"] > 2  # type: ignore[operator]
    assert cert_irred["status"] == "CERTIFIED_IRREDUCIBLE"

    # [H-2] tunnel_diagrams_c discards 2-edge-cut composites; are_prime_diagrams_tunnel_merged_c must NOT conflate:
    klein_tube_3_1 = simplify_compiled_diagram(
        compile_mixed_plat_link(CLASSICAL_PLAT_BY_LABEL["3_1"], {0: 0, 1: 1}),
        max_bfs_states=32,
    )
    k2_0 = MarkedDiagram.from_epd(MINIMAL_KLEIN_BOTTLE_E0_EPD)
    g0 = planar_connected_sum(d_8_1, 0, k2_0, 0, flip_splice=False)
    assert klein_tube_3_1.num_vertices == 11 and g0.num_vertices == 11
    assert len(tunnel_diagrams_c(11, [g0], depth_plus1=16, depth_plus2=8, max_bfs_states=128)) == 0
    assert (
        are_prime_diagrams_tunnel_merged_c(
            11, klein_tube_3_1, g0, depth_plus1=16, depth_plus2=8, max_bfs_states=128
        )
        is False
    )
    assert (
        are_prime_diagrams_tunnel_merged_c(
            8, d_8_1, d_8_1.mirror_3d(), depth_plus1=8, depth_plus2=4, max_bfs_states=64
        )
        is True
    )


def test_chiral_noninvertible_connected_sums_and_nonsplit_certificate() -> None:
    """Regression test for [H-3] 4-variant + flip_splice connected sums and [H-4] non-split certificate."""
    from surfaceknots.ch11_census import _oriented_diagram_variants, disjoint_union
    from surfaceknots.ch12_census import build_def5_1_high_composite_matches
    from surfaceknots.families import (
        _is_certified_nonsplit_multi_component,
        compile_mixed_plat_link,
        simplify_compiled_diagram,
    )

    # [H-3a] 8_1^-1,-1 has mirror_3d isomorphic to itself, requiring mirror_time (*t) to reach its second oriented class:
    d_8_mm = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["8_1^-1,-1"].canonical_epd)
    vars_8_mm = _oriented_diagram_variants("8_1^-1,-1", d_8_mm)
    assert [lbl for lbl, _ in vars_8_mm] == ["8_1^-1,-1", "8_1^-1,-1*t"]

    # [H-3b] Non-invertible 11_2 connected sums fl=False vs fl=True have distinct Alexander polynomials:
    ch11 = {e["label"]: e for e in json.loads((DATA_DIR / "expanded_ch11.json").read_text())["entries"]}
    d_11_2_s2 = MarkedDiagram.from_epd(ch11["11_2"]["epd"])
    cs_11_f = planar_connected_sum(d_11_2_s2, 0, d_11_2_s2, 0, flip_splice=False)
    cs_11_t = planar_connected_sum(d_11_2_s2, 0, d_11_2_s2, 0, flip_splice=True)
    assert compute_alexander_invariants(cs_11_f).single_var_alex != (
        compute_alexander_invariants(cs_11_t).single_var_alex
    )

    # [H-3c] Chiral 10_2 # 10_2 (fl=True) vs 10_2 # (-10_2^*) (fl=False) separated by Mochizuki c_3 & matched in Def 5.1:
    d_10_2 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_2"].canonical_epd)
    cs_10_opp = planar_connected_sum(d_10_2, 0, d_10_2, 0, flip_splice=False)
    cs_10_same = planar_connected_sum(d_10_2, 0, d_10_2, 0, flip_splice=True)
    rec_opp = compute_full_invariant_record(cs_10_opp, label="10_2_opp")
    rec_same = compute_full_invariant_record(cs_10_same, label="10_2_same")
    dk_opp = invariant_discrimination_key(rec_opp)
    dk_same = invariant_discrimination_key(rec_same)
    assert dk_opp != dk_same
    assert rec_opp["mochizuki_c3_mod3"] != rec_same["mochizuki_c3_mod3"]
    def51_matches = build_def5_1_high_composite_matches(12, [rec_opp, rec_same], DATA_DIR, num_workers=2)
    assert dk_opp in def51_matches
    assert dk_same in def51_matches

    # [H-4] _is_certified_nonsplit_multi_component must reject ALL 8 Non/Mixed split unions (0 false positives):
    s2 = MarkedDiagram.unknot_sphere()
    t2 = construct_minimal_unknot_diagram(1, 0)
    p_raw = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^-1"].canonical_epd)
    p_pos = p_raw if p_raw.surface_topology().components[0].normal_euler == 2 else p_raw.mirror_3d()
    k2_0 = construct_minimal_unknot_diagram(-2, 0)
    k2_4 = construct_minimal_unknot_diagram(-2, 4)
    d_8_1 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["8_1"].canonical_epd)

    split_diags = [
        disjoint_union(s2, p_pos),
        disjoint_union(p_pos, p_pos),
        disjoint_union(d_8_1, p_pos),
        disjoint_union(t2, p_pos),
        disjoint_union(k2_0, p_pos),
        disjoint_union(k2_4, k2_4),
        disjoint_union(d_8_1, k2_0),
        disjoint_union(d_8_1, k2_4),
        disjoint_union(s2, s2),
        disjoint_union(d_8_1, t2),
    ]
    for sp_d in split_diags:
        rec_sp = compute_full_invariant_record(sp_d, label="split_test", include_sublinks=False)
        assert _is_certified_nonsplit_multi_component(rec_sp) is False

    # Genuine non-split Mixed/Ori links must still certify True:
    km_l5a1 = simplify_compiled_diagram(
        compile_mixed_plat_link(CLASSICAL_PLAT_BY_LABEL["L5a1"], {1: 1}), max_bfs_states=32
    )
    rec_km = compute_full_invariant_record(km_l5a1, label="KleinMixed(L5a1,1)", include_sublinks=False)
    assert _is_certified_nonsplit_multi_component(rec_km) is True

    tube_borr = simplify_compiled_diagram(
        compile_satoh_tube(CLASSICAL_PLAT_BY_LABEL["L6a4"]), max_bfs_states=32
    )
    rec_borr = compute_full_invariant_record(tube_borr, label="Tube(L6a4)", include_sublinks=True)
    assert _is_certified_nonsplit_multi_component(rec_borr) is True


def test_bigon_10_1_relabeling_invariance_and_so4_verifier() -> None:
    """Regression test for [M-2] 10_1 bigon relabeling + B3 cord key and [M-3/L-2] SO(4) & flype checks."""
    from surfaceknots.families import simplify_compiled_diagram
    from surfaceknots.move_search import analyze_marked_bigon_1handles

    def _permute_vertices(d: MarkedDiagram, perm: list[int]) -> MarkedDiagram:
        n = d.num_vertices
        inv_p = [0] * n
        for ov, nv in enumerate(perm):
            inv_p[nv] = ov
        new_vt = tuple(d.vtype[inv_p[nv]] for nv in range(n))
        new_opp = [0] * (4 * n)
        for nv in range(n):
            ov = inv_p[nv]
            for s in range(4):
                od = d.opp[(ov << 2) | s]
                new_opp[(nv << 2) | s] = (perm[od >> 2] << 2) | (od & 3)
        return MarkedDiagram(num_vertices=n, opp=tuple(new_opp), vtype=new_vt, free_loops=d.free_loops)

    # [M-2a] All 24 cyclic + reversed-cyclic relabelings of Tube(4_1) must classify as base "10_1":
    tube_4_1 = simplify_compiled_diagram(
        compile_satoh_tube(CLASSICAL_PLAT_BY_LABEL["4_1"]), max_bfs_states=32
    )
    n12 = tube_4_1.num_vertices
    perms24 = [[(v + k) % n12 for v in range(n12)] for k in range(n12)] + [
        [(k - v) % n12 for v in range(n12)] for k in range(n12)
    ]
    for p in perms24:
        recs = analyze_marked_bigon_1handles(_permute_vertices(tube_4_1, p))
        assert any(
            r.base_identification == "10_1"
            and r.cord_double_coset_key == "4_1:norm_r=0:triv=False"
            and r.is_trivial_handle is False
            for r in recs
        )

    # [M-2b] ch=11 residual non-orientable longitude handle over 8_1 has exact cord key B3:('r=0', -1, -6):
    c11_residuals = json.loads((DATA_DIR / "expanded_ch11.json").read_text())["residual_collisions"]
    matched_b3_longitude = [
        hb.cord_double_coset_key
        for rc in c11_residuals
        for hb in analyze_marked_bigon_1handles(MarkedDiagram.from_epd(str(rc["epd"])))
        if hb.base_identification == "8_1" and not hb.handle_orientable
    ]
    assert "B3:('r=0', -1, -6)" in matched_b3_longitude

    # [M-3] reverse_orientation_zt() has det = +1 in SO(4): accepted under allow_o4=False, mirror_3d rejected!
    d_10_2 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_2"].canonical_epd)
    ver_so4 = verify_move_path_independent(
        [d_10_2.to_epd(), d_10_2.reverse_orientation_zt().to_epd()], allow_o4=False
    )
    assert ver_so4.is_valid is True
    assert ver_so4.step_kinds == ("ISO_SO4",)

    ver_mirror_so4 = verify_move_path_independent(
        [d_10_2.to_epd(), d_10_2.mirror_3d().to_epd()], allow_o4=False
    )
    assert ver_mirror_so4.is_valid is False

    ver_mirror_o4 = verify_move_path_independent(
        [d_10_2.to_epd(), d_10_2.mirror_3d().to_epd()], allow_o4=True
    )
    assert ver_mirror_o4.is_valid is True
    assert ver_mirror_o4.step_kinds == ("MIRROR_O4",)

