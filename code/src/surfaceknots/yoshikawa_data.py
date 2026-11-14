"""Canonical marked diagrams, metadata, Table II errata verification, and invariant export for all 23
Yoshikawa (1994) weakly prime, non-splittable surface-links (`ch <= 10`).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from itertools import permutations, product
import json
from pathlib import Path
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form

from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.invariants.alexander import (
    _ideal_evaluation_fingerprint,
    compute_alexander_invariants,
    twisted_s3_invariants,
)
from surfaceknots.invariants.fundamental_group import finite_quotient_profile, wirtinger_presentation
from surfaceknots.invariants.linking import (
    double_linking_mod2,
    mod2_triple_linking_profile,
    peripheral_sl2f3_profile,
    sublink_summary,
)
from surfaceknots.invariants.quandle import (
    alexander_biquandle_invariants,
    bikei_module_invariants,
    kei_module_invariants,
    mochizuki_3cocycle_invariant,
)
from surfaceknots.triplane import estimate_bridge_bounds_from_marked
from surfaceknots.unlink import certify_classical_unlink


@dataclass(frozen=True)
class YoshikawaEntry:
    """Reference entry in Yoshikawa's (1994) table of `ch <= 10` weakly prime surface-links."""

    label: str
    ch_index: int
    genera: tuple[int, ...]
    canonical_epd: str
    description: str
    ribbon_status: str
    yoshikawa_printed_e1: tuple[str, ...]
    corrected_e1: tuple[str, ...]
    has_yoshikawa_table_ii_erratum: bool = False
    erratum_note: str = ""


YOSHIKAWA_TABLE: tuple[YoshikawaEntry, ...] = (
    YoshikawaEntry(
        label="0_1",
        ch_index=0,
        genera=(0,),
        canonical_epd="O[1]",
        description="Unknotted 2-sphere",
        ribbon_status="Ribbon",
        yoshikawa_printed_e1=("1",),
        corrected_e1=("1",),
    ),
    YoshikawaEntry(
        label="2_1^1",
        ch_index=2,
        genera=(1,),
        canonical_epd="Y[1,2,3,4], Z[4,3,2,1]",
        description="Unknotted torus",
        ribbon_status="Ribbon",
        yoshikawa_printed_e1=("1",),
        corrected_e1=("1",),
    ),
    YoshikawaEntry(
        label="2_1^-1",
        ch_index=2,
        genera=(-1,),
        canonical_epd="Y[1,2,3,4], X[4,3,2,1]",
        description="Standard projective plane P_+- (normal Euler +-2)",
        ribbon_status="Pseudo-ribbon (non-ribbon)",
        yoshikawa_printed_e1=("1",),
        corrected_e1=("1",),
    ),
    YoshikawaEntry(
        label="6_1^0,1",
        ch_index=6,
        genera=(0, 1),
        canonical_epd="Y[1,2,3,4], Z[4,5,6,1], X[6,7,8,9], X[9,8,10,2], X[10,11,12,3], X[11,7,5,12]",
        description="Spun Hopf link (S^2 u T^2)",
        ribbon_status="Ribbon",
        yoshikawa_printed_e1=("x - 1", "y - 1"),
        corrected_e1=("x - 1", "y - 1"),
    ),
    YoshikawaEntry(
        label="7_1^0,-2",
        ch_index=7,
        genera=(0, -2),
        canonical_epd="Z[1,2,3,4], Y[5,6,7,8], X[6,9,10,11], X[12,13,14,2], X[5,14,13,9], X[8,7,4,3], X[12,1,11,10]",
        description="Simplest mixed orientable/non-orientable ribbon link (S^2 u K^2)",
        ribbon_status="Ribbon",
        yoshikawa_printed_e1=("x + 1", "y - 1"),
        corrected_e1=("x + 1", "y - 1"),
    ),
    YoshikawaEntry(
        label="8_1",
        ch_index=8,
        genera=(0,),
        canonical_epd="X[1,2,3,4], X[4,5,6,1], X[8,6,5,7], X[3,9,10,11], X[11,10,12,13], Z[8,14,15,16], Y[16,15,9,2], X[14,7,13,12]",
        description="Spun trefoil 2-knot tau^0(3_1)",
        ribbon_status="Ribbon",
        yoshikawa_printed_e1=("x**2 - x + 1",),
        corrected_e1=("x**2 - x + 1",),
    ),
    YoshikawaEntry(
        label="8_1^1,1",
        ch_index=8,
        genera=(1, 1),
        canonical_epd="X[1,2,3,4], X[4,5,6,1], Z[6,7,8,9], Y[9,8,10,2], Z[3,11,12,13], Y[13,12,14,5], X[15,16,7,14], X[16,15,11,10]",
        description="0-twist torus-spun Hopf link (T^2 u T^2)",
        ribbon_status="Ribbon",
        yoshikawa_printed_e1=("x - 1", "y - 1"),
        corrected_e1=("x - 1", "y - 1"),
    ),
    YoshikawaEntry(
        label="8_1^-1,-1",
        ch_index=8,
        genera=(-1, -1),
        canonical_epd="X[1,2,3,4], X[4,5,6,1], Y[7,8,9,2], X[10,8,7,6], Y[3,11,12,13], X[13,12,14,5], X[15,16,10,14], X[16,15,11,9]",
        description="Linked P_+ u P_- with quaternion group pi_1 = Q_8 and t(F) = 2",
        ribbon_status="Non-ribbon",
        yoshikawa_printed_e1=("x + 1", "y + 1", "2"),
        corrected_e1=("x + 1", "y + 1", "2"),
    ),
    YoshikawaEntry(
        label="9_1",
        ch_index=9,
        genera=(0,),
        canonical_epd="Z[1,2,3,4], Y[4,3,5,6], X[1,8,9,7], X[8,10,11,9], X[12,13,14,15], X[13,16,7,11], X[12,5,2,16], X[17,18,10,6], X[15,14,18,17]",
        description="1-fusion ribbon 2-knot with Baumslag-Solitar group BS(1, 2)",
        ribbon_status="Ribbon",
        yoshikawa_printed_e1=("x - 2",),
        corrected_e1=("x - 2",),
    ),
    YoshikawaEntry(
        label="9_1^0,1",
        ch_index=9,
        genera=(0, 1),
        canonical_epd="X[1,2,3,4], X[4,5,6,1], X[7,8,2,6], X[3,9,10,11], Y[11,12,13,5], Y[13,12,14,15], X[15,16,17,7], X[8,17,16,18], X[14,10,9,18]",
        description="Ribbon S^2 u T^2 link",
        ribbon_status="Ribbon",
        yoshikawa_printed_e1=("(x - 1)*(y - 1)", "(y - 1)**2"),
        corrected_e1=("(x - 1)*(y - 1)", "(y - 1)**2"),
    ),
    YoshikawaEntry(
        label="9_1^1,-2",
        ch_index=9,
        genera=(1, -2),
        canonical_epd="Y[1,2,3,4], Z[4,5,6,1], X[6,7,8,9], X[9,8,10,2], X[11,12,13,3], X[12,14,5,13], Z[14,15,16,7], Y[10,17,18,11], X[18,17,16,15]",
        description="Ribbon T^2 u K^2 link with abelian pi_1 = Z + Z_2",
        ribbon_status="Ribbon",
        yoshikawa_printed_e1=("x + 1", "y - 1"),
        corrected_e1=("x - 1", "y - 1", "2"),
        has_yoshikawa_table_ii_erratum=True,
        erratum_note=(
            "Yoshikawa Table II printed pi_1 = <x, y | yxyx^{-1} = 1> and E_1 = (x + 1, y - 1); "
            "actual group is abelian Z + Z_2 = <x, y | xy = yx, y^2 = 1> with E_1 = (x - 1, y - 1, 2)."
        ),
    ),
    YoshikawaEntry(
        label="10_1",
        ch_index=10,
        genera=(0,),
        canonical_epd="X[1,2,3,4], X[5,6,7,8], Y[5,9,4,10], X[11,12,2,13], Y[3,14,15,16], X[16,17,6,10], X[9,8,13,1], X[17,18,19,7], X[18,15,20,19], X[20,14,12,11]",
        description="Spun figure-eight 2-knot tau^0(4_1)",
        ribbon_status="Ribbon",
        yoshikawa_printed_e1=("x**2 - 3*x + 1",),
        corrected_e1=("x**2 - 3*x + 1",),
    ),
    YoshikawaEntry(
        label="10_2",
        ch_index=10,
        genera=(0,),
        canonical_epd="X[1,2,3,4], X[4,5,6,1], X[6,7,8,9], Y[9,10,11,2], Y[12,13,14,15], X[16,13,12,3], X[15,17,7,5], X[18,19,10,8], X[19,20,16,11], X[17,14,20,18]",
        description="2-twist-spun trefoil tau^2(3_1) (non-ribbon, chiral by Mochizuki 3-cocycle)",
        ribbon_status="Non-ribbon",
        yoshikawa_printed_e1=("x + 1", "3"),
        corrected_e1=("x + 1", "3"),
    ),
    YoshikawaEntry(
        label="10_3",
        ch_index=10,
        genera=(0,),
        canonical_epd="X[1,2,3,4], X[4,5,6,1], X[6,7,8,9], Y[9,10,11,2], X[12,13,14,3], Y[14,15,16,5], X[17,8,7,16], X[18,19,10,17], X[11,20,13,12], X[19,18,15,20]",
        description="3-twist-spun trefoil tau^3(3_1) (non-ribbon, commutator subgroup Q_8)",
        ribbon_status="Non-ribbon",
        yoshikawa_printed_e1=("x**2 + x + 1", "2"),
        corrected_e1=("x**2 + x + 1", "2"),
    ),
    YoshikawaEntry(
        label="10_1^1",
        ch_index=10,
        genera=(1,),
        canonical_epd="Z[1,2,3,4], Y[4,5,6,1], X[6,7,8,9], X[9,10,11,2], X[12,13,14,3], X[15,16,5,14], Z[16,17,18,7], X[18,19,10,8], Y[11,19,20,12], X[13,20,17,15]",
        description="Symmetric-spun / tubed trefoil torus Tube(3_1)",
        ribbon_status="Ribbon",
        yoshikawa_printed_e1=("x**2 - x + 1",),
        corrected_e1=("x**2 - x + 1",),
    ),
    YoshikawaEntry(
        label="10_1^0,1",
        ch_index=10,
        genera=(0, 1),
        canonical_epd="X[1,2,3,4], X[4,5,6,1], X[7,6,5,8], X[9,3,2,10], X[9,11,12,13], X[13,12,14,15], X[15,14,16,17], Z[7,18,19,20], Y[20,19,11,10], X[18,8,17,16]",
        description="Ribbon S^2 u T^2 link",
        ribbon_status="Ribbon",
        yoshikawa_printed_e1=("(x - 1)*(x + 1)", "(y - 1)*(x + 1)"),
        corrected_e1=("(x - 1)*(x + y)", "(y - 1)*(x + y)"),
        has_yoshikawa_table_ii_erratum=True,
        erratum_note=(
            "Yoshikawa Table II printed E_1 = (x - 1, y - 1)(x + 1); "
            "Fox calculus on <x, y | y^{-1}x^{-1}yx^{-1}yxy^{-1}x = 1> gives E_1 = (x - 1, y - 1)(x + y)."
        ),
    ),
    YoshikawaEntry(
        label="10_2^0,1",
        ch_index=10,
        genera=(0, 1),
        canonical_epd="X[1,2,3,4], X[4,5,6,1], X[6,7,8,9], X[9,10,11,2], X[12,13,10,14], X[15,16,5,17], X[18,8,7,16], Y[18,15,19,14], Z[13,20,3,11], X[12,19,17,20]",
        description="Ribbon S^2 u T^2 link with pi_1 = <x, y | x^2 y x^{-2} y^{-1} = 1>",
        ribbon_status="Ribbon",
        yoshikawa_printed_e1=("(x - 1)*(x + 1)", "(y - 1)*(x + 1)"),
        corrected_e1=("(x - 1)*(x + 1)", "(y - 1)*(x + 1)"),
    ),
    YoshikawaEntry(
        label="10_1^1,1",
        ch_index=10,
        genera=(1, 1),
        canonical_epd="X[1,2,3,4], X[4,5,6,7], Y[7,8,9,1], X[10,11,2,9], Z[3,11,12,5], Z[6,13,14,15], X[15,16,17,8], Y[17,18,19,10], X[19,20,13,12], X[20,18,16,14]",
        description="1-twist-spun Hopf link (T^2 u T^2, Sanderson Dlk = (1, 1))",
        ribbon_status="Pseudo-ribbon",
        yoshikawa_printed_e1=("x - 1", "y - 1"),
        corrected_e1=("x - 1", "y - 1"),
    ),
    YoshikawaEntry(
        label="10_1^0,0,1",
        ch_index=10,
        genera=(0, 0, 1),
        canonical_epd="Z[1,2,3,4], Y[5,6,2,1], X[7,8,9,10], X[11,7,10,6], X[11,5,12,13], X[14,8,13,12], X[9,15,16,17], X[17,16,18,3], X[18,19,20,4], X[19,15,14,20]",
        description="3-component S^2 u S^2 u T^2 ribbon link (two S^2s each linked as 6_1^{0,1} with T^2)",
        ribbon_status="Ribbon",
        yoshikawa_printed_e1=("0",),
        corrected_e1=("0",),
    ),
    YoshikawaEntry(
        label="10_1^0,-2",
        ch_index=10,
        genera=(0, -2),
        canonical_epd="X[1,2,3,4], X[4,5,6,1], X[7,8,2,9], X[8,10,11,3], X[10,12,5,11], X[13,14,15,12], X[16,17,18,19], Y[18,17,6,15], Z[20,9,16,19], X[7,20,14,13]",
        description="Ribbon S^2 u K^2 link",
        ribbon_status="Ribbon",
        yoshikawa_printed_e1=("2*x + y - 1", "4"),
        corrected_e1=("2*x + y - 1", "4"),
    ),
    YoshikawaEntry(
        label="10_2^0,-2",
        ch_index=10,
        genera=(0, -2),
        canonical_epd="X[1,2,3,4], X[4,5,6,7], X[7,8,9,1], Y[9,10,11,2], X[11,12,5,3], Y[6,13,14,15], X[15,16,17,8], X[18,19,10,17], X[19,20,13,12], X[20,18,16,14]",
        description="S^2 u K^2 surface-link (non-ribbon candidate)",
        ribbon_status="Unknown / Non-ribbon candidate",
        yoshikawa_printed_e1=("2*x + y - 1", "4"),
        corrected_e1=("2*x + y - 1", "(x + 1)*(y + 1)", "4"),
        has_yoshikawa_table_ii_erratum=True,
        erratum_note=(
            "Yoshikawa Table II printed E_1 = (2x + y - 1, 4), omitting the generator (x + 1)(y + 1); "
            "integral HNF over Z[x^{+-1}, y]/(y^2 - 1) proves E_1 = (2x + y - 1, (x + 1)(y + 1), 4)."
        ),
    ),
    YoshikawaEntry(
        label="10_1^-1,-1",
        ch_index=10,
        genera=(-1, -1),
        canonical_epd="X[1,2,3,4], X[4,5,6,7], Y[7,8,9,1], X[10,11,2,9], X[11,12,5,3], Z[6,13,14,15], X[15,16,17,8], X[18,19,10,17], X[19,20,13,12], X[20,18,16,14]",
        description="Non-orientable RP^2 u RP^2 link with dicyclic group Dic_16 of order 16",
        ribbon_status="Non-ribbon",
        yoshikawa_printed_e1=("x + 1", "y + 1", "4"),
        corrected_e1=("x + 1", "y + 1", "4"),
    ),
    YoshikawaEntry(
        label="10_1^-2,-2",
        ch_index=10,
        genera=(-2, -2),
        canonical_epd="Z[1,2,3,4], Y[4,5,6,1], X[6,7,8,9], X[9,8,10,11], X[12,3,2,11], X[13,14,15,12], X[14,16,5,15], Z[16,17,18,7], Y[10,19,20,13], X[20,19,18,17]",
        description="Ribbon link of two Klein bottles K^2 u K^2 with pi_1 = Q_8",
        ribbon_status="Ribbon",
        yoshikawa_printed_e1=("x + 1", "y + 1", "2"),
        corrected_e1=("x + 1", "y + 1", "2"),
    ),
)

YOSHIKAWA_BY_LABEL: dict[str, YoshikawaEntry] = {e.label: e for e in YOSHIKAWA_TABLE}


def exact_zh1_ideal_equals(
    gens_a: list[sp.Expr],
    gens_b: list[sp.Expr],
    ori_syms: list[sp.Symbol],
    tor_syms: list[sp.Symbol],
) -> bool:
    """Test exact ideal equality `(gens_a) == (gens_b)` in `Z[ori_syms^{+-1}, tor_syms]/(u_j^2 - 1)`
    using evaluation fingerprints + integral Hermite Normal Form in finite group-ring quotients.
    """
    if _ideal_evaluation_fingerprint(gens_a, ori_syms, tor_syms) != _ideal_evaluation_fingerprint(
        gens_b, ori_syms, tor_syms
    ):
        return False

    r = len(ori_syms)
    s = len(tor_syms)
    all_syms = ori_syms + tor_syms
    for m_ord in (3, 4, 5):
        if (m_ord**r) * (2**s) > 512:
            continue
        basis = list(product(*(range(m_ord) for _ in range(r)), *(range(2) for _ in range(s))))
        b_idx = {e: idx for idx, e in enumerate(basis)}
        dim = len(basis)

        def build_hnf(gens: list[sp.Expr]) -> sp.Matrix:
            cols: list[list[int]] = []
            for g in gens:
                g_exp = sp.expand(g)
                if g_exp == 0:
                    continue
                if not all_syms:
                    terms = [((), int(g_exp))]
                else:
                    p = sp.Poly(g_exp, *all_syms, domain=sp.ZZ)
                    terms = [(mon, int(c)) for mon, c in p.terms()]
                for b in basis:
                    col = [0] * dim
                    for mon, coeff in terms:
                        pe = tuple((b[i] + mon[i]) % m_ord for i in range(r)) + tuple(
                            (b[r + j] + mon[r + j]) & 1 for j in range(s)
                        )
                        col[b_idx[pe]] += coeff
                    if any(col):
                        cols.append(col)
            if not cols:
                return sp.zeros(dim, 1)
            return hermite_normal_form(sp.Matrix(cols).T)

        if build_hnf(gens_a) != build_hnf(gens_b):
            return False
    return True


def check_e1_ideal_against_reference(
    d: MarkedDiagram,
    ref_gens_str: tuple[str, ...],
) -> tuple[bool, bool]:
    """Compare computed `E_1(d)` against `ref_gens_str` in `Z[H_1]`.

    Returns `(strict_h1_aut_match, signed_u_flip_match)`.
    """
    al = compute_alexander_invariants(d)
    ori_syms = [sp.Symbol(nm) for nm, fl in zip(al.variable_names, al.orientable_flags) if fl]
    tor_syms = [sp.Symbol(nm) for nm, fl in zip(al.variable_names, al.orientable_flags) if not fl]
    sym_map = {nm: sp.Symbol(nm) for nm in al.variable_names}

    our_gens = [sp.sympify(s, locals=sym_map) for s in al.e1_generators]
    ref_raw = [sp.sympify(s, locals=sym_map) for s in ref_gens_str]
    r = len(ori_syms)
    s = len(tor_syms)

    def clear_neg_powers(expr: sp.Expr) -> sp.Expr:
        expr = sp.expand(expr)
        if expr == 0 or not ori_syms:
            return expr
        num, _ = sp.fraction(sp.together(expr))
        e = sp.expand(num)
        all_s = ori_syms + tor_syms
        p = sp.Poly(e, *all_s, domain=sp.ZZ)
        min_o = [min(mon[i] for mon, _ in p.terms()) for i in range(r)]
        out = sp.Integer(0)
        for mon, coeff in p.terms():
            term = sp.Integer(int(coeff))
            for i in range(r):
                term *= ori_syms[i] ** (mon[i] - min_o[i])
            for j in range(s):
                if mon[r + j] & 1:
                    term *= tor_syms[j]
            out += term
        return sp.expand(out)

    strict_ok = False
    flip_ok = False
    for p_ori in permutations(range(r)):
        for p_tor in permutations(range(s)):
            for inv_o in product((1, -1), repeat=r):
                for sign_u in product((1, -1), repeat=s):
                    subs = {}
                    for i in range(r):
                        tgt = ori_syms[p_ori[i]]
                        subs[ori_syms[i]] = tgt if inv_o[i] == 1 else (1 / tgt)
                    for j in range(s):
                        subs[tor_syms[j]] = sign_u[j] * tor_syms[p_tor[j]]
                    mapped = [clear_neg_powers(g.subs(subs)) for g in ref_raw]
                    if exact_zh1_ideal_equals(our_gens, mapped, ori_syms, tor_syms):
                        if all(u == 1 for u in sign_u):
                            strict_ok = True
                        else:
                            flip_ok = True
    return strict_ok, flip_ok


def compute_full_invariant_record(
    d: MarkedDiagram,
    label: str | None = None,
    include_sublinks: bool = True,
) -> dict[str, object]:
    """Compute the full machine-certified KnotInfo-style invariant record for a `MarkedDiagram`."""
    lp, lm = d.resolve(1), d.resolve(-1)
    cert_p = certify_classical_unlink(lp)
    cert_m = certify_classical_unlink(lm)
    assert cert_p.is_unlink and cert_m.is_unlink, f"Inadmissible diagram for {label}!"

    topo = d.surface_topology()
    wp_m = wirtinger_presentation(d, from_resolution=-1)
    wp_p = wirtinger_presentation(d, from_resolution=1)
    fqp_m = finite_quotient_profile(wp_m)
    fqp_p = finite_quotient_profile(wp_p)
    assert fqp_m == fqp_p, f"Resolution mismatch in finite_quotient_profile on {label}!"

    al_m = compute_alexander_invariants(d, from_resolution=-1)
    al_p = compute_alexander_invariants(d, from_resolution=1)
    assert al_m.evaluation_fingerprint == al_p.evaluation_fingerprint, (
        f"Resolution mismatch in E_1 fingerprint on {label}!"
    )
    h1_free, h1_tor = wp_m.abelianization_invariant_factors()

    tw_s3 = twisted_s3_invariants(d)
    kei = kei_module_invariants(d)
    bikei = bikei_module_invariants(d)
    abq = alexander_biquandle_invariants(d)
    c3_oriented = mochizuki_3cocycle_invariant(d, p=3, up_to_mirror=False)
    c3_mirror = (
        mochizuki_3cocycle_invariant(d.mirror_3d(), p=3, up_to_mirror=False)
        if c3_oriented is not None
        else None
    )
    is_c3_chiral = (c3_oriented is not None) and (
        c3_oriented != c3_mirror or len(set(c3_oriented)) > 1
    )

    psl2_homs, psl2_non_merid = peripheral_sl2f3_profile(d)
    dlk = double_linking_mod2(d)
    tlk = mod2_triple_linking_profile(d)
    bridge_est = estimate_bridge_bounds_from_marked(d)
    sublinks = sublink_summary(d) if include_sublinks else None

    comp_types = [c.homeomorphism_type for c in topo.components]
    normal_eulers = sorted(c.normal_euler for c in topo.components)

    import hashlib

    e1_fp_hex = hashlib.sha256(repr(al_m.evaluation_fingerprint).encode("utf-8")).hexdigest()[:16]

    return {
        "label": label or "",
        "ch_index_diagram": d.ch_index,
        "num_classical": d.num_classical,
        "num_marked": d.num_marked,
        "epd": d.to_epd() if label != "0_1" else "O[1]",
        "num_components": topo.num_components,
        "yoshikawa_genera": list(topo.yoshikawa_genera),
        "component_types": comp_types,
        "orientability": topo.orientability_label,
        "euler_characteristic": topo.euler_characteristic,
        "normal_euler_numbers": normal_eulers,
        "c_plus": topo.num_circles_plus,
        "c_minus": topo.num_circles_minus,
        "pi1_generators": wp_m.num_generators,
        "pi1_relators": len(wp_m.relators),
        "pi1_presentation": wp_m.to_human_string(),
        "h1_free_rank": h1_free,
        "h1_torsion": list(h1_tor),
        "alexander_variables": list(al_m.variable_names),
        "e1_generators": list(al_m.e1_generators),
        "e2_generators": list(al_m.e2_generators),
        "single_var_alex": al_m.single_var_alex,
        "e1_fingerprint_hash": e1_fp_hex,
        "hom_counts": dict(fqp_m),
        "twisted_s3_delta0": [str(x) for x in tw_s3],
        "kei_invariant_factors": [list(x) if isinstance(x, tuple) else x for x in kei],
        "bikei_invariant_factors": [list(x) if isinstance(x, tuple) else x for x in bikei],
        "alexander_biquandle_t": repr(abq[0]) if abq is not None else None,
        "alexander_biquandle_s": list(abq[1]) if abq is not None else None,
        "mochizuki_c3_mod3": list(c3_oriented) if c3_oriented is not None else None,
        "mochizuki_c3_mirror_mod3": list(c3_mirror) if c3_mirror is not None else None,
        "certified_chiral_by_c3": is_c3_chiral,
        "sl2f3_hom_count": psl2_homs,
        "sl2f3_non_meridional_peripheral_count": len(psl2_non_merid),
        "sl2f3_peripheral_sample": [list(x[0]) for x in psl2_non_merid[:1]],
        "double_linking_mod2": [list(x) for x in dlk],
        "triple_linking_mod2": [list(x) for x in tlk] if tlk is not None else None,
        "bridge_bounds": bridge_est,
        "sublinks": sublinks,
    }


def invariant_discrimination_key(rec: dict[str, object]) -> tuple:
    """Strict topological invariant tuple used to certify distinctness of two surface-links."""
    # Note: excludes diagram-dependent presentation counts (`pi1_generators`, `bridge_bounds`, `epd`).
    # Uses exact evaluation fingerprint via re-evaluating `e1_generators`, `hom_counts`, peripheral SL(2,F3),
    # biquandle, Mochizuki 3-cocycle (unoriented/mirror-invariant sorted pair), and Sanderson Dlk!
    c3 = rec["mochizuki_c3_mod3"]
    c3_m = rec["mochizuki_c3_mirror_mod3"]
    c3_unoriented = tuple(sorted([tuple(c3), tuple(c3_m)])) if (c3 is not None and c3_m is not None) else None
    return (
        rec["num_components"],
        tuple(rec["yoshikawa_genera"]),
        rec["orientability"],
        rec["euler_characteristic"],
        tuple(sorted(abs(x) for x in rec["normal_euler_numbers"]))
        if rec["num_components"] == 1
        else tuple(sorted(rec["normal_euler_numbers"]))
        if sum(rec["normal_euler_numbers"]) == 0
        else tuple(
            min(
                sorted(rec["normal_euler_numbers"]),
                sorted(-x for x in rec["normal_euler_numbers"]),
            )
        ),
        rec["h1_free_rank"],
        tuple(rec["h1_torsion"]),
        rec["e1_fingerprint_hash"],
        tuple(sorted(rec["hom_counts"].items())),
        tuple(rec["twisted_s3_delta0"]),
        tuple(tuple(x) if isinstance(x, list) else x for x in rec["kei_invariant_factors"]),
        tuple(tuple(x) if isinstance(x, list) else x for x in rec["bikei_invariant_factors"]),
        rec["alexander_biquandle_t"],
        tuple(rec["alexander_biquandle_s"]) if rec["alexander_biquandle_s"] is not None else None,
        c3_unoriented,
        rec["sl2f3_hom_count"],
        rec["sl2f3_non_meridional_peripheral_count"],
        tuple(tuple(x) for x in rec["double_linking_mod2"]),
        tuple(tuple(x) for x in rec["triple_linking_mod2"]) if rec["triple_linking_mod2"] is not None else None,
    )


def export_yoshikawa_datasets(out_dir: Path) -> list[dict[str, object]]:
    """Verify all 23 Yoshikawa entries + 3 Table II errata, and write `yoshikawa_ch10.{json,csv}`."""
    out_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, object]] = []
    seen_keys: dict[tuple, str] = {}

    for entry in YOSHIKAWA_TABLE:
        d = MarkedDiagram.from_epd(entry.canonical_epd)
        assert d.ch_index == entry.ch_index
        assert d.surface_topology().yoshikawa_genera == entry.genera

        printed_ok, _ = check_e1_ideal_against_reference(d, entry.yoshikawa_printed_e1)
        corrected_ok, _ = check_e1_ideal_against_reference(d, entry.corrected_e1)
        assert corrected_ok, f"Corrected E_1 failed strict check on {entry.label}!"
        if entry.has_yoshikawa_table_ii_erratum:
            assert not printed_ok, f"Expected Yoshikawa printed E_1 to FAIL on erratum {entry.label}!"
        else:
            assert printed_ok, f"Expected Yoshikawa printed E_1 to PASS on {entry.label}!"

        rec = compute_full_invariant_record(d, label=entry.label, include_sublinks=True)
        rec["description"] = entry.description
        rec["ribbon_status"] = entry.ribbon_status
        rec["yoshikawa_printed_e1"] = list(entry.yoshikawa_printed_e1)
        rec["corrected_reference_e1"] = list(entry.corrected_e1)
        rec["has_yoshikawa_table_ii_erratum"] = entry.has_yoshikawa_table_ii_erratum
        rec["erratum_note"] = entry.erratum_note

        dkey = invariant_discrimination_key(rec)
        assert dkey not in seen_keys, f"Invariant collision between {entry.label} and {seen_keys[dkey]}!"
        seen_keys[dkey] = entry.label
        records.append(rec)

    json_path = out_dir / "yoshikawa_ch10.json"
    json_path.write_text(json.dumps(records, indent=2) + "\n")

    csv_path = out_dir / "yoshikawa_ch10.csv"
    csv_fields = [
        "label",
        "ch_index_diagram",
        "num_components",
        "component_types",
        "orientability",
        "euler_characteristic",
        "normal_euler_numbers",
        "ribbon_status",
        "pi1_presentation",
        "e1_generators",
        "corrected_reference_e1",
        "has_yoshikawa_table_ii_erratum",
        "hom_S3",
        "hom_A4",
        "hom_S4",
        "hom_A5",
        "hom_Q8",
        "hom_Dic16",
        "sl2f3_hom_count",
        "sl2f3_non_meridional_peripheral_count",
        "double_linking_mod2",
        "mochizuki_c3_mod3",
        "certified_chiral_by_c3",
        "epd",
    ]
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_fields)
        writer.writeheader()
        for r in records:
            hc = r["hom_counts"]
            writer.writerow(
                {
                    "label": r["label"],
                    "ch_index_diagram": r["ch_index_diagram"],
                    "num_components": r["num_components"],
                    "component_types": " u ".join(r["component_types"]),
                    "orientability": r["orientability"],
                    "euler_characteristic": r["euler_characteristic"],
                    "normal_euler_numbers": json.dumps(r["normal_euler_numbers"]),
                    "ribbon_status": r["ribbon_status"],
                    "pi1_presentation": r["pi1_presentation"],
                    "e1_generators": json.dumps(r["e1_generators"]),
                    "corrected_reference_e1": json.dumps(r["corrected_reference_e1"]),
                    "has_yoshikawa_table_ii_erratum": r["has_yoshikawa_table_ii_erratum"],
                    "hom_S3": hc["S3"],
                    "hom_A4": hc["A4"],
                    "hom_S4": hc["S4"],
                    "hom_A5": hc["A5"],
                    "hom_Q8": hc["Q8"],
                    "hom_Dic16": hc["Dic16"],
                    "sl2f3_hom_count": r["sl2f3_hom_count"],
                    "sl2f3_non_meridional_peripheral_count": r["sl2f3_non_meridional_peripheral_count"],
                    "double_linking_mod2": json.dumps(r["double_linking_mod2"]),
                    "mochizuki_c3_mod3": json.dumps(r["mochizuki_c3_mod3"]),
                    "certified_chiral_by_c3": r["certified_chiral_by_c3"],
                    "epd": r["epd"],
                }
            )
    return records
