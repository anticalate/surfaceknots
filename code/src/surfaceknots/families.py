"""Parametric family compilers and Hub-and-Spoke cross-join for knotted surfaces in R^4.

Implements constructive, planar, admissibility-certified `MarkedDiagram` compilers for major
infinite families of 2-knots and surface-links, together with automated invariant cross-joining
against the exhaustive `ch <= 11` Hub catalog:

- **Spoke 1 (Spun & Symmetric-Union 2-Knots)**:
  - `compile_artin_spin`: Artin's spun 2-knot `tau^0(K)` of any `b`-bridge classical knot `K`
    with constructive upper bound `ch(tau^0(K)) <= 2 * c_plat(K) + 2 * (b - 1)`.
  - `compile_symmetric_union_2knot`: Eisermann-Lamm `m`-twisted symmetric-union ribbon 2-knot
    `U_m(K)` with constructive upper bound `ch(U_m(K)) <= 2 * c_plat(K) + 2 * (b - 1) + |m|`.
- **Spoke 2 (1-Fusion Ribbon 2-Knots)**:
  - `compile_1fusion_ribbon_word`: General NORTH/SOUTH alternating planar compiler for the
    1-fusion ribbon 2-knot `R(w)` of any reduced word `w in F(a, b)`, realizing group
    `pi_1 = <a, b | w^{-1} a w = b>` and exact Fox Alexander polynomial
    `Delta_{R(w)}(t) = 1 + (t - 1) * (dw / da)|_{a=b=t}` in `4 * |w| + 2` vertices.
- **Spoke 3 (Satoh Tubes, Mixed Sphere-Torus Links, & Twisted Klein-Bottle Tubes)**:
  - `compile_satoh_tube`: Satoh's ribbon torus-link `Tube(L)` of any `b`-bridge classical link `L`
    with constructive upper bound `ch(Tube(L)) <= 2 * c_plat(L) + 2 * b`.
  - `compile_mixed_plat_link`: Mixed sphere / torus / Klein-bottle surface-links obtained by
    attaching `m_k`-twisted `Z-Y` 1-handles across a subset of middle bridge necks of `L`.
"""

from __future__ import annotations

from collections import defaultdict
import csv
from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

from surfaceknots.ch11_census import build_reference_catalog
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.moves import (
    _is_planar_4reg_diagram,
    canonical_marked_key,
    greedy_yoshikawa_reduce,
    list_4d_tait_flypes,
    list_level_yoshikawa_neighbors,
    simplify_and_explore_yoshikawa,
)
from surfaceknots.unlink import certify_classical_unlink
from surfaceknots.yoshikawa_data import (
    assign_nearest_separating_witnesses,
    compute_full_invariant_record,
    invariant_discrimination_key,
)


@dataclass(frozen=True)
class ClassicalPlatSpec:
    """Plat specification of a classical `b`-bridge knot or link on `2*b` strands."""

    label: str
    display_name: str
    num_components: int
    bridge_number: int
    crossing_number: int
    braid_word: tuple[int, ...]
    conway_cf: tuple[int, ...] | None = None

    @property
    def num_strands(self) -> int:
        return 2 * self.bridge_number

    @property
    def rational_fraction(self) -> str | None:
        if not self.conway_cf:
            return "0/1" if self.crossing_number == 0 else None
        val = Fraction(self.conway_cf[-1], 1)
        for a in reversed(self.conway_cf[:-1]):
            val = Fraction(a, 1) + Fraction(1, val)
        return f"{val.numerator}/{val.denominator}"


@dataclass(frozen=True)
class RibbonWordSpec:
    """Specification of a 1-fusion ribbon 2-knot `R(w)` by a word `w` in `F(a, b)`."""

    family_id: str
    word_str: str
    word: tuple[tuple[str, int], ...]
    notes: str = ""


def conway_4plat_word(continued_fraction: tuple[int, ...]) -> tuple[int, ...]:
    """Convert an odd-length additive Conway continued fraction `[a_1, ..., a_{2k+1}]` into a
    4-strand plat braid word alternating generator `2` (strands `1, 2`) and `-1` (strands `0, 1`).

    Starting and ending on generator `2` guarantees that end twists never untwist across the
    top/bottom bridge caps `(0, 1)` and `(2, 3)`.
    """
    if not continued_fraction:
        return ()
    if len(continued_fraction) % 2 != 1:
        raise ValueError(f"Conway 4-plat continued fraction must have odd length: {continued_fraction}")
    out: list[int] = []
    for idx, c in enumerate(continued_fraction):
        if c == 0:
            raise ValueError(f"Zero entry in Conway continued fraction: {continued_fraction}")
        gen = 2 if (idx % 2 == 0) else -1
        s = gen if c > 0 else -gen
        out.extend([s] * abs(c))
    return tuple(out)


def parse_ribbon_word(word_str: str) -> tuple[tuple[str, int], ...]:
    """Parse a human-readable word like `'b^-1 a b^2 a^-1'` or `'1'` into `(('b', -1), ...)`."""
    s = word_str.strip()
    if not s or s == "1":
        return ()
    out: list[tuple[str, int]] = []
    for tok in s.split():
        if "^" in tok:
            base, exp_s = tok.split("^", 1)
            exp = int(exp_s)
        else:
            base, exp = tok, 1
        if base not in ("a", "b") or exp == 0:
            raise ValueError(f"Invalid ribbon word token {tok!r} in {word_str!r}")
        sgn = 1 if exp > 0 else -1
        out.extend([(base, sgn)] * abs(exp))
    return tuple(out)


def fox_1fusion_alexander_forms(word: tuple[tuple[str, int], ...]) -> tuple[str, str]:
    """Return both `+- t^k`-normalized representatives `(p(t), t^deg * p(t^{-1}))` of the
    single-variable Alexander polynomial `Delta_{R(w)}(t) = 1 + (t - 1) * (dw / da)|_{a=b=t}`
    under `Aut(H_1) = {t -> t, t -> t^{-1}}`, ordered canonically `(min, max)`.
    """
    t = sp.Symbol("t")
    deriv = sp.Integer(0)
    prefix_pow = 0
    for gen, sgn in word:
        if gen == "a":
            if sgn == 1:
                deriv += t**prefix_pow
            else:
                deriv -= t ** (prefix_pow - 1)
        prefix_pow += sgn

    raw = sp.expand(1 + (t - 1) * deriv)
    num, _ = sp.fraction(sp.together(raw))
    poly = sp.Poly(sp.expand(num), t, domain=sp.ZZ)
    terms = [(m[0], int(c)) for m, c in poly.terms() if c != 0]
    min_e = min(e for e, _ in terms)
    max_e = max(e for e, _ in terms)
    deg = max_e - min_e
    fwd = sorted(((e - min_e, c) for e, c in terms), reverse=True)
    if fwd[0][1] < 0:
        fwd = [(e, -c) for e, c in fwd]
    rev = sorted(((deg - (e - min_e), c) for e, c in terms), reverse=True)
    if rev[0][1] < 0:
        rev = [(e, -c) for e, c in rev]
    low, high = (fwd, rev) if fwd <= rev else (rev, fwd)
    s_low = str(sp.expand(sum(sp.Integer(c) * (t**e) for e, c in low)))
    s_high = str(sp.expand(sum(sp.Integer(c) * (t**e) for e, c in high)))
    return (s_low, s_high)


def fox_1fusion_alexander_poly(word: tuple[tuple[str, int], ...]) -> str:
    """Compute the canonical single-variable Alexander polynomial of `R(w)` up to `+- t^{+-k}`."""
    return fox_1fusion_alexander_forms(word)[0]


def compile_plat_surface(
    num_strands: int,
    braid_word: tuple[int, ...],
    neck_twists: dict[int, int],
) -> MarkedDiagram:
    """Compile the planar doubled `b`-bridge plat diagram `Plat(beta * beta^{-1})` with
    `m_k`-twisted horizontal `Z-Y` 1-handle gadgets inserted at `k in neck_twists`.

    Theorem (Twisted Ribbon 1-Handle Pinch):
    For every `b >= 1`, every braid word `beta in B_{2b}`, and every assignment `neck_twists[k] = m_k`:
    1. The resulting 4-regular graph `D` is strictly planar (`_is_planar_4reg_diagram(D) == True`)
       with `ch(D) = 2 * len(beta) + sum_k (2 + |m_k|)`.
    2. In both resolutions `L_+(D)` and `L_-(D)`, each horizontal `Z-Y` neck resolves to one
       vertical pass-through and one chain of `|m_k|` Reidemeister I kinks, which R1-reduces to
       parallel vertical strands so `beta` cancels `beta^{-1}` by Reidemeister II moves to the
       `b`-component classical unlink (`D` is 100% admissible without search).
    """
    if num_strands <= 0 or num_strands % 2 != 0:
        raise ValueError(f"num_strands must be positive and even, got {num_strands}")
    b = num_strands // 2
    for k in neck_twists:
        if not (0 <= k < b):
            raise ValueError(f"Neck cap index {k} out of range 0..{b - 1}")

    vtypes: list[str] = []
    adj: list[list[tuple[int, int] | None]] = []

    def new_v(vt: str) -> int:
        vid = len(vtypes)
        vtypes.append(vt)
        adj.append([None, None, None, None])
        return vid

    def link(va: int, sa: int, vb: int, sb: int) -> None:
        adj[va][sa] = (vb, sb)
        adj[vb][sb] = (va, sa)

    curr: list[tuple[int | None, int]] = [(None, i) for i in range(num_strands)]
    top_stubs: list[tuple[int, int] | None] = [None] * num_strands

    def feed_strand(s_idx: int, target: tuple[int, int]) -> None:
        if curr[s_idx][0] is None:
            top_stubs[s_idx] = target
        else:
            va, sa = curr[s_idx]
            assert va is not None
            link(va, sa, target[0], target[1])

    def append_gen(g: int) -> None:
        i = abs(g) - 1
        if not (0 <= i < num_strands - 1):
            raise ValueError(f"Braid generator {g} out of range for {num_strands} strands")
        vt = "X0" if g > 0 else "X1"
        v = new_v(vt)
        # CCW cyclic order around v: 0=NE (i+1 in), 1=NW (i in), 2=SW (i out), 3=SE (i+1 out)
        feed_strand(i + 1, (v, 0))
        feed_strand(i, (v, 1))
        curr[i] = (v, 2)
        curr[i + 1] = (v, 3)

    for g in braid_word:
        append_gen(g)

    for k in range(b):
        if k in neck_twists:
            m = neck_twists[k]
            sL, sR = 2 * k, 2 * k + 1
            vz = new_v("Z")
            vy = new_v("Y")
            feed_strand(sL, (vz, 1))
            feed_strand(sR, (vy, 0))
            prev_top = (vz, 0)
            prev_bot = (vz, 3)
            vt_tw = "X0" if m > 0 else "X1"
            for _ in range(abs(m)):
                vx = new_v(vt_tw)
                link(prev_top[0], prev_top[1], vx, 1)
                link(prev_bot[0], prev_bot[1], vx, 2)
                prev_top = (vx, 0)
                prev_bot = (vx, 3)
            link(prev_top[0], prev_top[1], vy, 1)
            link(prev_bot[0], prev_bot[1], vy, 2)
            curr[sL] = (vz, 2)
            curr[sR] = (vy, 3)

    for g in reversed(braid_word):
        append_gen(-g)

    idle_loops = 0
    for k in range(b):
        sL, sR = 2 * k, 2 * k + 1
        tL, tR = top_stubs[sL], top_stubs[sR]
        if tL is None and tR is None:
            idle_loops += 1
        elif tL is None:
            assert tR is not None and curr[sR][0] is not None
            link(tR[0], tR[1], curr[sR][0], curr[sR][1])
        elif tR is None:
            assert tL is not None and curr[sL][0] is not None
            link(tL[0], tL[1], curr[sL][0], curr[sL][1])
        else:
            assert curr[sL][0] is not None and curr[sR][0] is not None
            link(curr[sL][0], curr[sL][1], curr[sR][0], curr[sR][1])
            link(tL[0], tL[1], tR[0], tR[1])

    n = len(vtypes)
    if n == 0:
        return MarkedDiagram(num_vertices=0, vtype=(), opp=(), free_loops=idle_loops)
    opp = [0] * (4 * n)
    for v in range(n):
        for s in range(4):
            nb = adj[v][s]
            assert nb is not None
            opp[4 * v + s] = 4 * nb[0] + nb[1]
    return MarkedDiagram(num_vertices=n, vtype=tuple(vtypes), opp=tuple(opp), free_loops=idle_loops)


def compile_artin_spin(spec: ClassicalPlatSpec) -> MarkedDiagram:
    """Compile Artin's spun 2-knot `tau^0(K)` of a `b`-bridge classical knot `spec`."""
    if spec.num_components != 1:
        raise ValueError(f"compile_artin_spin expects a 1-component knot, got {spec.label}")
    necks = {k: 0 for k in range(spec.bridge_number - 1)}
    return compile_plat_surface(spec.num_strands, spec.braid_word, necks)


def compile_symmetric_union_2knot(spec: ClassicalPlatSpec, twist: int = 1) -> MarkedDiagram:
    """Compile the `twist`-half-twisted symmetric-union ribbon 2-knot `U_m(K)` of `spec`."""
    if spec.num_components != 1:
        raise ValueError(f"compile_symmetric_union_2knot expects a 1-component knot, got {spec.label}")
    necks = {k: (twist if k == 0 else 0) for k in range(spec.bridge_number - 1)}
    return compile_plat_surface(spec.num_strands, spec.braid_word, necks)


def compile_satoh_tube(spec: ClassicalPlatSpec) -> MarkedDiagram:
    """Compile Satoh's ribbon torus-link `Tube(L)` of a `b`-bridge classical knot or link `spec`."""
    necks = {k: 0 for k in range(spec.bridge_number)}
    return compile_plat_surface(spec.num_strands, spec.braid_word, necks)


def compile_mixed_plat_link(spec: ClassicalPlatSpec, neck_twists: dict[int, int]) -> MarkedDiagram:
    """Compile a mixed sphere / torus / Klein-bottle surface-link from `spec` with `neck_twists`."""
    return compile_plat_surface(spec.num_strands, spec.braid_word, neck_twists)


def compile_1fusion_ribbon_word(word: tuple[tuple[str, int], ...]) -> MarkedDiagram:
    """Compile a reduced word `w` in `F(a, b)` into a planar admissible `MarkedDiagram` of the
    1-fusion ribbon 2-knot `R(w)` with `pi_1(S^4 \\ R(w)) = <a, b | w^{-1} a w = b>`.

    Uses NORTH/SOUTH half-plane alternation across the West tip of `O_a` and East tip of `O_b`
    so the 2-rail ribbon never encloses either component regardless of exponent lengths.
    """
    vtypes: list[str] = []
    adj: list[list[tuple[int, int] | None]] = []

    def new_v(vt: str) -> int:
        vid = len(vtypes)
        vtypes.append(vt)
        adj.append([None, None, None, None])
        return vid

    def link(va: int, sa: int, vb: int, sb: int) -> None:
        adj[va][sa] = (vb, sb)
        adj[vb][sb] = (va, sa)

    vz = new_v("Z")
    vy = new_v("Y")
    a_crossings: list[tuple[int, int, int, int]] = []
    b_crossings: list[tuple[int, int, int, int]] = []

    cur_west: tuple[int, int] = (vz, 0)
    cur_east: tuple[int, int] = (vz, 3)
    side = "NORTH"

    for gen, sgn in word:
        if gen not in ("a", "b") or sgn not in (1, -1):
            raise ValueError(f"Invalid ribbon word generator {(gen, sgn)!r}")
        entry_vt = "X1" if sgn > 0 else "X0"
        exit_vt = "X0" if sgn > 0 else "X1"

        if side == "NORTH":
            vt_W, vt_E = new_v(entry_vt), new_v(entry_vt)
            vb_W, vb_E = new_v(exit_vt), new_v(exit_vt)
            if gen == "b":
                link(cur_east[0], cur_east[1], vt_W, 0)
                link(cur_west[0], cur_west[1], vt_E, 0)
            else:
                link(cur_west[0], cur_west[1], vt_E, 0)
                link(cur_east[0], cur_east[1], vt_W, 0)
            link(vt_W, 2, vb_W, 0)
            link(vt_E, 2, vb_E, 0)
            cur_west = (vb_W, 2)
            cur_east = (vb_E, 2)
            side = "SOUTH"
        else:
            vt_W, vt_E = new_v(exit_vt), new_v(exit_vt)
            vb_W, vb_E = new_v(entry_vt), new_v(entry_vt)
            if gen == "b":
                link(cur_east[0], cur_east[1], vb_W, 2)
                link(cur_west[0], cur_west[1], vb_E, 2)
            else:
                link(cur_west[0], cur_west[1], vb_E, 2)
                link(cur_east[0], cur_east[1], vb_W, 2)
            link(vb_W, 0, vt_W, 2)
            link(vb_E, 0, vt_E, 2)
            cur_west = (vt_W, 0)
            cur_east = (vt_E, 0)
            side = "NORTH"

        link(vt_W, 3, vt_E, 1)
        link(vb_W, 3, vb_E, 1)

        if gen == "a":
            a_crossings.insert(0, (vt_W, vt_E, vb_W, vb_E))
        else:
            b_crossings.append((vt_W, vt_E, vb_W, vb_E))

    if not a_crossings:
        link(vz, 1, vz, 2)
    else:
        link(a_crossings[0][0], 1, a_crossings[0][2], 1)
        for idx in range(len(a_crossings) - 1):
            link(a_crossings[idx][1], 3, a_crossings[idx + 1][0], 1)
            link(a_crossings[idx][3], 3, a_crossings[idx + 1][2], 1)
        link(a_crossings[-1][1], 3, vz, 1)
        link(a_crossings[-1][3], 3, vz, 2)

    if side == "NORTH":
        link(cur_east[0], cur_east[1], vy, 0)
        link(cur_west[0], cur_west[1], vy, 3)
    else:
        link(cur_east[0], cur_east[1], vy, 3)
        link(cur_west[0], cur_west[1], vy, 0)

    if not b_crossings:
        link(vy, 1, vy, 2)
    else:
        link(b_crossings[0][0], 1, b_crossings[0][2], 1)
        for idx in range(len(b_crossings) - 1):
            link(b_crossings[idx][1], 3, b_crossings[idx + 1][0], 1)
            link(b_crossings[idx][3], 3, b_crossings[idx + 1][2], 1)
        link(b_crossings[-1][1], 3, vy, 1)
        link(b_crossings[-1][3], 3, vy, 2)

    n = len(vtypes)
    opp = [0] * (4 * n)
    for v in range(n):
        for s in range(4):
            nb = adj[v][s]
            assert nb is not None
            opp[4 * v + s] = 4 * nb[0] + nb[1]
    return MarkedDiagram(num_vertices=n, vtype=tuple(vtypes), opp=tuple(opp))


def _plat_2b(label: str, name: str, comps: int, cf: tuple[int, ...]) -> ClassicalPlatSpec:
    bw = conway_4plat_word(cf)
    return ClassicalPlatSpec(
        label=label,
        display_name=name,
        num_components=comps,
        bridge_number=1 if not cf and comps == 1 else 2,
        crossing_number=len(bw),
        braid_word=bw,
        conway_cf=cf,
    )


# Exhaustive table of prime classical knots up to 7 crossings (all 2-bridge) + classical links up to 6 crossings:
CLASSICAL_PLAT_TABLE: tuple[ClassicalPlatSpec, ...] = (
    ClassicalPlatSpec("0_1", "Unknot", 1, 1, 0, (), ()),
    _plat_2b("3_1", "Trefoil knot", 1, (1, 1, 1)),
    _plat_2b("4_1", "Figure-eight knot", 1, (1, 1, 2)),
    _plat_2b("5_1", "(2,5)-torus cinquefoil knot", 1, (1, 3, 1)),
    _plat_2b("5_2", "Three-twist knot 5_2", 1, (1, 1, 3)),
    _plat_2b("6_1", "Stevedore knot 6_1", 1, (1, 1, 4)),
    _plat_2b("6_2", "Knot 6_2", 1, (1, 1, 1, 2, 1)),
    _plat_2b("6_3", "Knot 6_3", 1, (1, 1, 1, 1, 2)),
    _plat_2b("7_1", "(2,7)-torus septafoil knot", 1, (1, 5, 1)),
    _plat_2b("7_2", "Five-twist knot 7_2", 1, (1, 1, 5)),
    _plat_2b("7_3", "Knot 7_3", 1, (1, 2, 4)),
    _plat_2b("7_4", "Knot 7_4", 1, (1, 2, 1, 2, 1)),
    _plat_2b("7_5", "Knot 7_5", 1, (1, 1, 2, 2, 1)),
    _plat_2b("7_6", "Knot 7_6", 1, (1, 1, 1, 2, 2)),
    _plat_2b("7_7", "Knot 7_7", 1, (1, 1, 1, 1, 1, 1, 1)),
    _plat_2b("L2a1", "Hopf link", 2, (2,)),
    _plat_2b("L4a1", "Solomon's knot (2,4)-torus link", 2, (1, 2, 1)),
    _plat_2b("L5a1", "Whitehead link", 2, (1, 1, 1, 1, 1)),
    _plat_2b("L6a1", "Link L6a1 (6_3^2)", 2, (1, 1, 2, 1, 1)),
    _plat_2b("L6a2", "Link L6a2 (6_2^2)", 2, (1, 2, 3)),
    _plat_2b("L6a3", "(2,6)-torus link L6a3 (6_1^2)", 2, (1, 4, 1)),
    ClassicalPlatSpec(
        label="L6a4",
        display_name="Borromean rings (Brunnian 3-bridge link)",
        num_components=3,
        bridge_number=3,
        crossing_number=6,
        braid_word=(4, 2, -3, 2, -3, 2, -3, -4),
        conway_cf=None,
    ),
)

CLASSICAL_PLAT_BY_LABEL: dict[str, ClassicalPlatSpec] = {s.label: s for s in CLASSICAL_PLAT_TABLE}


RIBBON_1FUSION_TABLE: tuple[RibbonWordSpec, ...] = (
    RibbonWordSpec("R(1)", "1", (), "Trivial 1-fusion unknotted sphere"),
    RibbonWordSpec("R(b*a)", "b a", parse_ribbon_word("b a"), "Spun trefoil 8_1 via 2-syllable 1-fusion"),
    RibbonWordSpec("R(b*a^-1)", "b a^-1", parse_ribbon_word("b a^-1"), "Ribbon 2-knot 9_1 (Delta = 2x - 1)"),
    RibbonWordSpec(
        "R(b^-1*a*b*a^-1)",
        "b^-1 a b a^-1",
        parse_ribbon_word("b^-1 a b a^-1"),
        "Spun figure-eight 10_1 via commutator [b, a^-1]",
    ),
    RibbonWordSpec(
        "R(b^-1*a*b^-1*a^-1)",
        "b^-1 a b^-1 a^-1",
        parse_ribbon_word("b^-1 a b^-1 a^-1"),
        "ch=11 ribbon 2-knot 11_1 (Delta = 2x^2 - 2x + 1)",
    ),
    RibbonWordSpec(
        "R(b^-1*a*b*a)",
        "b^-1 a b a",
        parse_ribbon_word("b^-1 a b a"),
        "ch=11 asymmetric ribbon 2-knot 11_2 (Delta = x^3 - x^2 + 2x - 1)",
    ),
    RibbonWordSpec(
        "R(b^2*a)",
        "b^2 a",
        parse_ribbon_word("b^2 a"),
        "Non-reciprocal degree-3 1-fusion ribbon 2-knot (Delta = x^3 - x + 1)",
    ),
    RibbonWordSpec(
        "R(b^2*a^-1)",
        "b^2 a^-1",
        parse_ribbon_word("b^2 a^-1"),
        "Golden-ratio non-reciprocal 1-fusion ribbon 2-knot (Delta = x^2 - x - 1)",
    ),
    RibbonWordSpec(
        "R(b^3*a)",
        "b^3 a",
        parse_ribbon_word("b^3 a"),
        "Non-reciprocal degree-4 1-fusion ribbon 2-knot (Delta = x^4 - x + 1)",
    ),
    RibbonWordSpec(
        "R(b^3*a^-1)",
        "b^3 a^-1",
        parse_ribbon_word("b^3 a^-1"),
        "Non-reciprocal degree-3 1-fusion ribbon 2-knot (Delta = x^3 - x^2 - 1)",
    ),
    RibbonWordSpec(
        "R(b^2*a^2)",
        "b^2 a^2",
        parse_ribbon_word("b^2 a^2"),
        "Symmetric degree-4 1-fusion ribbon 2-knot (Delta = x^4 - x^2 + 1)",
    ),
    RibbonWordSpec(
        "R(b^2*a^-2)",
        "b^2 a^-2",
        parse_ribbon_word("b^2 a^-2"),
        "Non-monic non-reciprocal 1-fusion ribbon 2-knot (Delta = x^2 - 2)",
    ),
    RibbonWordSpec(
        "R(b^-1*a^2*b*a^-1)",
        "b^-1 a^2 b a^-1",
        parse_ribbon_word("b^-1 a^2 b a^-1"),
        "4-syllable non-reciprocal 1-fusion ribbon 2-knot (Delta = x^3 - 2x^2 - x + 1)",
    ),
    RibbonWordSpec(
        "R(b^2*a*b^-1*a)",
        "b^2 a b^-1 a",
        parse_ribbon_word("b^2 a b^-1 a"),
        "Non-monic non-reciprocal 1-fusion ribbon 2-knot (Delta = x^3 - 2x^2 + 2)",
    ),
)


def simplify_compiled_diagram(d: MarkedDiagram, max_bfs_states: int = 6) -> MarkedDiagram:
    """Simplify an admissible planar `MarkedDiagram` via greedy Yoshikawa reduction, bounded
    level-move exploration, and deduplicated 4D Tait flypes.
    """
    if not _is_planar_4reg_diagram(d):
        raise ValueError("Compiled diagram failed strict 4-regular planarity check!")
    if not (
        certify_classical_unlink(d.resolve(1)).is_unlink
        and certify_classical_unlink(d.resolve(-1)).is_unlink
    ):
        raise ValueError("Compiled diagram failed hyperbolic splitting unlink admissibility!")

    best = greedy_yoshikawa_reduce(d)
    best = simplify_and_explore_yoshikawa(best, max_states=max_bfs_states)[0]
    seen_fl: set[tuple] = {canonical_marked_key(best, up_to_mirror_and_reverse=True)}
    checked_flypes = 0
    for fl in list_4d_tait_flypes(best):
        fk = canonical_marked_key(fl, up_to_mirror_and_reverse=True)
        if fk in seen_fl:
            continue
        seen_fl.add(fk)
        g_fl = greedy_yoshikawa_reduce(fl)
        if g_fl.ch_index < best.ch_index:
            best = simplify_and_explore_yoshikawa(g_fl, max_states=max_bfs_states)[0]
            continue
        if best.ch_index <= 12 and checked_flypes < 2:
            checked_flypes += 1
            cand = simplify_and_explore_yoshikawa(g_fl, max_states=4)[0]
            if cand.ch_index < best.ch_index:
                best = cand
    return best


def _enumerate_spoke_specs() -> list[dict[str, object]]:
    """Assemble the complete catalog of Spoke 1, Spoke 2, and Spoke 3 family instances."""
    specs: list[dict[str, object]] = []

    # Spoke 1a: Artin spins `tau^0(K)` for all 15 knots up to 7 crossings
    for kspec in CLASSICAL_PLAT_TABLE:
        if kspec.num_components != 1:
            continue
        specs.append(
            {
                "family_id": f"Spin({kspec.label})",
                "spoke": "Spoke 1: Spun & Symmetric-Union 2-Knots",
                "subfamily": "artin_spin",
                "source_label": kspec.label,
                "parameters": {
                    "classical_knot": kspec.label,
                    "bridge_number": kspec.bridge_number,
                    "crossing_number": kspec.crossing_number,
                    "conway_cf": list(kspec.conway_cf) if kspec.conway_cf is not None else None,
                    "rational_fraction": kspec.rational_fraction,
                    "twist": 0,
                },
                "diagram": compile_artin_spin(kspec),
                "description": f"Artin spin tau^0({kspec.label}) ({kspec.display_name})",
            }
        )

    # Spoke 1b: 1-twisted symmetric-union ribbon 2-knots `U_1(K)` for knots `3_1 .. 6_3`
    for k_lbl in ("3_1", "4_1", "5_1", "5_2", "6_1", "6_2", "6_3"):
        kspec = CLASSICAL_PLAT_BY_LABEL[k_lbl]
        specs.append(
            {
                "family_id": f"SymUnion({kspec.label},1)",
                "spoke": "Spoke 1: Spun & Symmetric-Union 2-Knots",
                "subfamily": "symmetric_union_2knot",
                "source_label": kspec.label,
                "parameters": {
                    "classical_knot": kspec.label,
                    "bridge_number": kspec.bridge_number,
                    "crossing_number": kspec.crossing_number,
                    "twist": 1,
                },
                "diagram": compile_symmetric_union_2knot(kspec, twist=1),
                "description": f"1-twisted symmetric-union ribbon 2-knot U_1({kspec.label})",
            }
        )

    # Spoke 2: 1-fusion ribbon 2-knots `R(w)`
    for rspec in RIBBON_1FUSION_TABLE:
        specs.append(
            {
                "family_id": rspec.family_id,
                "spoke": "Spoke 2: 1-Fusion Ribbon 2-Knots",
                "subfamily": "ribbon_1fusion",
                "source_label": rspec.word_str,
                "parameters": {
                    "word": rspec.word_str,
                    "word_length": len(rspec.word),
                    "fox_alexander_poly": fox_1fusion_alexander_poly(rspec.word),
                },
                "diagram": compile_1fusion_ribbon_word(rspec.word),
                "description": rspec.notes,
            }
        )

    # Spoke 3a: Satoh tubes `Tube(L)` for classical knots (`0_1 .. 6_1`) and links (`L2a1 .. L6a4`)
    for t_lbl in ("0_1", "3_1", "4_1", "5_1", "5_2", "6_1", "L2a1", "L4a1", "L5a1", "L6a1", "L6a2", "L6a3", "L6a4"):
        lspec = CLASSICAL_PLAT_BY_LABEL[t_lbl]
        specs.append(
            {
                "family_id": f"Tube({lspec.label})",
                "spoke": "Spoke 3: Satoh Tubes & Mixed Plat Surface-Links",
                "subfamily": "satoh_tube",
                "source_label": lspec.label,
                "parameters": {
                    "classical_base": lspec.label,
                    "num_base_components": lspec.num_components,
                    "bridge_number": lspec.bridge_number,
                    "neck_twists": {str(k): 0 for k in range(lspec.bridge_number)},
                },
                "diagram": compile_satoh_tube(lspec),
                "description": f"Satoh ribbon torus-link Tube({lspec.label}) ({lspec.display_name})",
            }
        )

    # Spoke 3b: Mixed Sphere-Torus & Non-Orientable Klein-Bottle Twisted Plat Surface-Links
    mixed_configs: list[tuple[str, str, dict[int, int], str]] = [
        ("Mixed(L2a1)", "L2a1", {0: 0}, "Mixed S^2 u T^2 Hopf surface-link"),
        ("Mixed(L4a1)", "L4a1", {0: 0}, "Mixed S^2 u T^2 Solomon surface-link (identifies 10_1^0,1)"),
        ("Mixed(L5a1)", "L5a1", {0: 0}, "Mixed S^2 u T^2 Whitehead surface-link (exact ch = 12)"),
        ("Mixed(L6a1)", "L6a1", {0: 0}, "Mixed S^2 u T^2 L6a1 surface-link"),
        ("Mixed(L6a2)", "L6a2", {0: 0}, "Mixed S^2 u T^2 L6a2 surface-link"),
        ("Mixed(L6a3)", "L6a3", {0: 0}, "Mixed S^2 u T^2 (2,6)-torus surface-link"),
        ("Mixed(L6a4,01)", "L6a4", {0: 0, 1: 0}, "Brunnian 3-component S^2 u T^2 u T^2 Borromean surface-link"),
        ("KleinMixed(L2a1,1)", "L2a1", {0: 1}, "Mixed S^2 u K^2 Hopf surface-link (identifies 7_1^0,-2)"),
        ("KleinMixed(L2a1,10)", "L2a1", {0: 1, 1: 0}, "Mixed T^2 u K^2 Hopf surface-link (identifies 9_1^1,-2)"),
        ("KleinTube(L2a1,11)", "L2a1", {0: 1, 1: 1}, "Double Klein-bottle K^2 u K^2 Hopf link (identifies 10_1^-2,-2)"),
        ("KleinMixed(L4a1,1)", "L4a1", {0: 1}, "Mixed S^2 u K^2 Solomon surface-link (identifies 11_3^0,-2)"),
        ("KleinMixed(L5a1,1)", "L5a1", {0: 1}, "Mixed S^2 u K^2 Whitehead surface-link (ch in [12, 13])"),
    ]
    for fid, base_lbl, ntw, desc in mixed_configs:
        lspec = CLASSICAL_PLAT_BY_LABEL[base_lbl]
        subf = "mixed_sphere_torus" if all(v == 0 for v in ntw.values()) else "klein_twisted_plat"
        specs.append(
            {
                "family_id": fid,
                "spoke": "Spoke 3: Satoh Tubes & Mixed Plat Surface-Links",
                "subfamily": subf,
                "source_label": base_lbl,
                "parameters": {
                    "classical_base": base_lbl,
                    "num_base_components": lspec.num_components,
                    "bridge_number": lspec.bridge_number,
                    "neck_twists": {str(k): v for k, v in ntw.items()},
                },
                "diagram": compile_mixed_plat_link(lspec, ntw),
                "description": desc,
            }
        )

    return specs


def _is_certified_nonsplit_multi_component(rec: dict[str, object]) -> bool:
    """Certify that a multi-component surface-link record cannot be a nontrivial split union."""
    n_comp = int(rec["num_components"])  # type: ignore[arg-type]
    if n_comp <= 1:
        return True
    # Any split union `F_1 u F_2` of orientable components with `b_1 >= 2` has `E_1(F) = (0)`.
    e1 = list(rec["e1_generators"])  # type: ignore[arg-type]
    if e1 != ["0"]:
        return True
    # For 3-component Brunnian links (where E_1 = (0)), compare |Hom(pi_1, G)| against free-product
    # factorization across proper sublinks if available (`L6a4` has non-free-product S_4 / SL(2,F_3) / A_4 counts!).
    hc = rec["hom_counts"]  # type: ignore[assignment]
    if n_comp == 3 and tuple(sorted(str(x) for x in rec["component_types"])) in (  # type: ignore[arg-type]
        ("S^2", "T^2", "T^2"),
        ("T^2", "T^2", "T^2"),
    ):
        # Any split 3-component link of `S^2` and `T^2` at ch <= 11 has at least one unknotted component
        # split off (`pi_1 = Z * G_2`), forcing `|Hom(pi_1, S_3)|` to be divisible by `|S_3| = 6`, AND if all
        # 2-component sublinks are trivial unlinks (Brunnian), a split link would be the 3-component unlink
        # (`|Hom(Z^3_free, S_3)| = 6^3 = 216`, `|Hom(Z^3_free, A_4)| = 12^3 = 1728`), whereas Borromean differs!
        sub = rec.get("sublinks")
        if isinstance(sub, dict) and sub.get("is_brunnian") and hc["A4"] != 12**3:  # type: ignore[index]
            return True
    return False


def build_family_catalog(data_dir: Path) -> dict[str, object]:
    """Compile all Spoke 1-3 family entries, compute full Tier 1-3 invariants, cross-join against
    the exhaustive `ch <= 11` Hub (`yoshikawa_ch10.json` + `expanded_ch11.json`) and composite/split
    reference catalog (`max_composite_ch = 12`), and certify exact `[ch_lower, ch_upper]` intervals.
    """
    yosh_path = data_dir / "yoshikawa_ch10.json"
    ch11_path = data_dir / "expanded_ch11.json"
    if not yosh_path.exists() or not ch11_path.exists():
        raise FileNotFoundError(
            f"Hub datasets {yosh_path} and {ch11_path} must exist before building family catalog."
        )

    yosh_entries: list[dict[str, object]] = json.loads(yosh_path.read_text())
    ch11_payload: dict[str, object] = json.loads(ch11_path.read_text())
    ch11_entries: list[dict[str, object]] = ch11_payload["entries"]  # type: ignore[assignment]
    ch11_residuals: list[dict[str, object]] = ch11_payload["residual_collisions"]  # type: ignore[assignment]

    # Index Hub entries by invariant_discrimination_key AND by canonical_marked_key of known EPDs:
    hub_by_inv: dict[tuple, dict[str, object]] = {}
    hub_epd_keys: dict[str, set[tuple]] = {}
    for he in yosh_entries + ch11_entries:
        lbl = str(he["label"])
        hub_by_inv[invariant_discrimination_key(he)] = he
        epd_set: set[tuple] = set()
        epds_list = [str(he["epd"])] + [str(x) for x in he.get("move_component_epds", [])]  # type: ignore[union-attr]
        for epd_s in epds_list:
            dh = MarkedDiagram.from_epd(epd_s)
            epd_set.add(canonical_marked_key(dh, up_to_mirror_and_reverse=True))
            for fl in list_4d_tait_flypes(dh):
                epd_set.add(canonical_marked_key(fl, up_to_mirror_and_reverse=True))
        hub_epd_keys[lbl] = epd_set

    residual_by_inv: dict[tuple, list[str]] = defaultdict(list)
    for re_rec in ch11_residuals:
        rk = invariant_discrimination_key(re_rec)
        for m_tag in re_rec.get("collision_matches", []):  # type: ignore[union-attr]
            if str(m_tag) not in residual_by_inv[rk]:
                residual_by_inv[rk].append(str(m_tag))

    ref_catalog = build_reference_catalog(max_composite_ch=12)

    spoke_specs = _enumerate_spoke_specs()
    entries_out: list[dict[str, object]] = []
    hub_realizations: dict[str, list[dict[str, object]]] = defaultdict(list)

    for sp_item in spoke_specs:
        fid = str(sp_item["family_id"])
        raw_d: MarkedDiagram = sp_item["diagram"]  # type: ignore[assignment]
        raw_ch = raw_d.ch_index
        simp_d = simplify_compiled_diagram(raw_d, max_bfs_states=6)
        compiled_ch = simp_d.ch_index

        need_sublinks = simp_d.surface_topology().num_components >= 3
        inv_rec = compute_full_invariant_record(
            simp_d,
            label=fid,
            include_sublinks=need_sublinks,
            compute_so4_symmetry_bfs=(compiled_ch <= 12),
        )

        # Verify Fox's 1-fusion Alexander polynomial formula on every Spoke 2 word!
        params = sp_item["parameters"]  # type: ignore[assignment]
        if sp_item["subfamily"] == "ribbon_1fusion":
            valid_forms = fox_1fusion_alexander_forms(parse_ribbon_word(str(params["word"])))  # type: ignore[index]
            actual_alex = inv_rec["single_var_alex"]
            assert actual_alex in valid_forms, (
                f"Fox formula mismatch on {fid}: expected one of {valid_forms}, got {actual_alex}"
            )

        dkey = invariant_discrimination_key(inv_rec)
        hub_hit = hub_by_inv.get(dkey)
        ref_hits = [t for t in ref_catalog.get(dkey, []) if not t.startswith("ch<=10:")]
        res_hits = residual_by_inv.get(dkey, [])

        hub_match: str | None = None
        hub_match_mode: str | None = None
        ch_lower: int
        ch_upper: int
        exact_ch_certified = False
        complexity_reason: str

        if hub_hit is not None:
            hub_match = str(hub_hit["label"])
            hub_ch = int(hub_hit["ch_index_diagram"])  # type: ignore[arg-type]
            hub_status = str(hub_hit["cell_status"])

            # Level moves, 4D Tait flypes, and C mountain-pass tunneling check when compiled_ch == hub_ch:
            move_matched = False
            if compiled_ch == hub_ch:
                simp_keys: set[tuple] = {canonical_marked_key(simp_d, up_to_mirror_and_reverse=True)}
                if simp_keys & hub_epd_keys[hub_match]:
                    move_matched = True
                else:
                    bfs_q: list[MarkedDiagram] = [simp_d]
                    q_idx = 0
                    while q_idx < len(bfs_q) and len(simp_keys) < 24 and not move_matched:
                        cur_d = bfs_q[q_idx]
                        q_idx += 1
                        for nb_d in list_level_yoshikawa_neighbors(cur_d) + list_4d_tait_flypes(cur_d):
                            nk = canonical_marked_key(nb_d, up_to_mirror_and_reverse=True)
                            if nk in hub_epd_keys[hub_match]:
                                move_matched = True
                                break
                            if nk not in simp_keys:
                                simp_keys.add(nk)
                                bfs_q.append(nb_d)
                                if len(simp_keys) >= 24:
                                    break
                if not move_matched and 3 <= compiled_ch <= 11 and simp_d.free_loops == 0:
                    from surfaceknots.fast_enum import tunnel_diagrams_c

                    dh = MarkedDiagram.from_epd(str(hub_hit["epd"]))
                    t_merged = tunnel_diagrams_c(
                        compiled_ch, [simp_d, dh], depth_plus1=16, depth_plus2=8, max_bfs_states=512
                    )
                    if len(t_merged) == 1:
                        move_matched = True

            other_res = [r for r in res_hits if r != f"ch<=10:{hub_match}"]
            if move_matched:
                hub_match_mode = "MOVE_CERTIFIED"
                ch_lower = hub_ch
                ch_upper = hub_ch
                exact_ch_certified = True
                complexity_reason = (
                    f"Move-equivalent to Hub entry {hub_match} (ch = {hub_ch}, {hub_status}); "
                    f"disjoint from all lower-ch/composite/split references."
                )
            else:
                hub_match_mode = "INVARIANT_MATCH"
                ch_lower = hub_ch if not ref_hits and not other_res else 0
                ch_upper = compiled_ch
                exact_ch_certified = ch_lower == ch_upper
                complexity_reason = (
                    f"Full 10-layer invariant match with Hub entry {hub_match} (ch = {hub_ch}, "
                    f"{hub_status}); compiled diagram has ch = {compiled_ch}."
                )

            hub_realizations[hub_match].append(
                {
                    "family_id": fid,
                    "spoke": sp_item["spoke"],
                    "subfamily": sp_item["subfamily"],
                    "match_mode": hub_match_mode,
                    "compiled_ch": compiled_ch,
                    "raw_ch": raw_ch,
                }
            )
        elif ref_hits or res_hits:
            # Collides with a known unknot / split / composite (`ch <= 12`) or ch=11 residual shadow:
            ch_lower = 0
            ch_upper = compiled_ch
            exact_ch_certified = compiled_ch == 0
            complexity_reason = (
                f"Invariant collides with non-prime/residual reference(s): {ref_hits + res_hits}."
            )
        else:
            # Disjoint from ALL ch <= 11 prime entries, ch=11 residual collisions, AND ch <= 12
            # unknot / split (`ch <= 11`) / composite (`ch <= 12`) references!
            nonsplit_ok = _is_certified_nonsplit_multi_component(inv_rec)
            if nonsplit_ok:
                ch_lower = 12
                ch_upper = compiled_ch
                exact_ch_certified = ch_lower == ch_upper
                if exact_ch_certified:
                    complexity_reason = (
                        "PROVEN EXACT ch = 12: constructive diagram has ch = 12 and full invariant "
                        "key is disjoint from the exhaustive ch <= 11 prime census, residual shadows, "
                        "and all ch <= 12 composite/split/unknot classes."
                    )
                else:
                    complexity_reason = (
                        f"Certified ch >= 12 (disjoint from exhaustive ch <= 11 prime, residual, and "
                        f"composite/split catalogs); constructive diagram gives ch <= {compiled_ch}."
                    )
            else:
                ch_lower = 0
                ch_upper = compiled_ch
                exact_ch_certified = False
                complexity_reason = "Multi-component link without algebraic non-split certificate."

        record = {
            "family_id": fid,
            "spoke": sp_item["spoke"],
            "subfamily": sp_item["subfamily"],
            "source_label": sp_item["source_label"],
            "description": sp_item["description"],
            "parameters": params,
            "raw_ch": raw_ch,
            "compiled_ch": compiled_ch,
            "ch_lower": ch_lower,
            "ch_upper": ch_upper,
            "exact_ch_certified": exact_ch_certified,
            "hub_match": hub_match,
            "hub_match_mode": hub_match_mode,
            "reference_collisions": ref_hits + res_hits,
            "complexity_certification_reason": complexity_reason,
            **inv_rec,
        }
        if hub_hit is not None and "nearest_separating_witness" in hub_hit:
            record["nearest_separating_witness"] = hub_hit["nearest_separating_witness"]
        entries_out.append(record)

    # Assign nearest separating witnesses across the union of the 33-entry Hub + new `ch >= 12` Spoke entries:
    new_spoke_pool = [e for e in entries_out if e["hub_match"] is None and not e["reference_collisions"]]
    # Deduplicate `new_spoke_pool` by invariant key so distinct Spoke constructions of the SAME `ch >= 12`
    # surface (e.g. if a ribbon word equals a spun/symmetric-union knot) don't self-collide in witness pool:
    unique_new_by_inv: dict[tuple, dict[str, object]] = {}
    for e in new_spoke_pool:
        k = invariant_discrimination_key(e)
        if k not in unique_new_by_inv:
            e["label"] = e["family_id"]
            unique_new_by_inv[k] = e

    combined_pool = yosh_entries + ch11_entries + list(unique_new_by_inv.values())
    assign_nearest_separating_witnesses(list(unique_new_by_inv.values()), reference_pool=combined_pool)
    for e in entries_out:
        k = invariant_discrimination_key(e)
        if k in unique_new_by_inv and "nearest_separating_witness" not in e:
            e["nearest_separating_witness"] = unique_new_by_inv[k]["nearest_separating_witness"]

    exact_ch12_entries = [
        e["family_id"]
        for e in entries_out
        if e["exact_ch_certified"] and e["ch_lower"] == 12 and e["ch_upper"] == 12
    ]

    return {
        "hub_max_ch": 11,
        "total_family_instances": len(entries_out),
        "hub_matched_instances": sum(1 for e in entries_out if e["hub_match"] is not None),
        "distinct_hub_entries_realized": len(hub_realizations),
        "certified_ch_ge_12_instances": sum(1 for e in entries_out if e["ch_lower"] >= 12),
        "distinct_ch_ge_12_invariant_classes": len(unique_new_by_inv),
        "exact_ch_12_certified_instances": exact_ch12_entries,
        "hub_realizations": {k: hub_realizations[k] for k in sorted(hub_realizations)},
        "entries": entries_out,
    }


def export_family_datasets(data_dir: Path) -> dict[str, object]:
    """Build and write `code/data/families/family_catalog.{json,csv}` and `hub_realizations.json`."""
    fam_dir = data_dir / "families"
    fam_dir.mkdir(parents=True, exist_ok=True)

    payload = build_family_catalog(data_dir)
    json_path = fam_dir / "family_catalog.json"
    json_path.write_text(json.dumps(payload, indent=2) + "\n")

    hub_rel_path = fam_dir / "hub_realizations.json"
    hub_rel_path.write_text(json.dumps(payload["hub_realizations"], indent=2) + "\n")

    csv_path = fam_dir / "family_catalog.csv"
    csv_fields = [
        "family_id",
        "spoke",
        "subfamily",
        "source_label",
        "raw_ch",
        "compiled_ch",
        "ch_lower",
        "ch_upper",
        "exact_ch_certified",
        "hub_match",
        "hub_match_mode",
        "num_components",
        "component_types",
        "orientability",
        "normal_euler_numbers",
        "e1_generators",
        "sigma2_b1",
        "sigma2_h1_torsion",
        "sigma2_b2",
        "sigma2_abs_signature",
        "certified_pi2_z_rank",
        "sl2f3_non_meridional_peripheral_count",
        "symmetry_summary",
        "nearest_separating_label",
        "nearest_separating_layer",
        "epd",
    ]
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_fields)
        writer.writeheader()
        for e in payload["entries"]:  # type: ignore[union-attr]
            wit = e.get("nearest_separating_witness") or {}
            sym_d = e["symmetry"]  # type: ignore[index]
            writer.writerow(
                {
                    "family_id": e["family_id"],
                    "spoke": e["spoke"],
                    "subfamily": e["subfamily"],
                    "source_label": e["source_label"],
                    "raw_ch": e["raw_ch"],
                    "compiled_ch": e["compiled_ch"],
                    "ch_lower": e["ch_lower"],
                    "ch_upper": e["ch_upper"],
                    "exact_ch_certified": e["exact_ch_certified"],
                    "hub_match": e["hub_match"] or "",
                    "hub_match_mode": e["hub_match_mode"] or "",
                    "num_components": e["num_components"],
                    "component_types": " u ".join(e["component_types"]),  # type: ignore[arg-type]
                    "orientability": e["orientability"],
                    "normal_euler_numbers": json.dumps(e["normal_euler_numbers"]),
                    "e1_generators": json.dumps(e["e1_generators"]),
                    "sigma2_b1": e["sigma2_b1"],
                    "sigma2_h1_torsion": json.dumps(e["sigma2_h1_torsion"]),
                    "sigma2_b2": e["sigma2_b2"],
                    "sigma2_abs_signature": e["sigma2_abs_signature"],
                    "certified_pi2_z_rank": e["certified_pi2_z_rank"],
                    "sl2f3_non_meridional_peripheral_count": e["sl2f3_non_meridional_peripheral_count"],
                    "symmetry_summary": sym_d["summary"],  # type: ignore[index]
                    "nearest_separating_label": wit.get("nearest_label", ""),
                    "nearest_separating_layer": wit.get("separating_layer", ""),
                    "epd": e["epd"],
                }
            )

    return payload
