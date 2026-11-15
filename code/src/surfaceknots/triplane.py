"""Bridge trisections, tri-plane diagrams (T_1, T_2, T_3), invariants, and ch-diagram compilation.

Implements:
- `TangleDiagram` & `generate_catalan_matchings`: `b`-bridge trivial tangles on `2b` boundary points
  `0..2b-1` generated from a crossingless planar arch matching capped by an explicit braid word in `B_{2b}`.
- `TriPlaneDiagram`: triple `(T_1, T_2, T_3)` of `b`-bridge tangles with:
  - exact pairwise classical link construction `D_{ij} = T_i union mirror(T_j)` (`ClassicalLinkDiagram`),
  - constructive Reidemeister unlink certification of all three pairwise unions via `certify_classical_unlink`,
  - exact patch counts `(p_12, p_23, p_31)`, Euler characteristic `chi = p_12 + p_23 + p_31 - b`,
  - full `surface_topology()` (connected components via `<tau_1, tau_2, tau_3>` orbits, per-component
    orientability via bridge-point bipartiteness, and per-component normal Euler number
    `e(F_c) = -(w_self(L_12|_c) + w_self(L_23|_c) + w_self(L_31|_c))` via Joseph-Meier-Miller-Zupan 2022),
  - tri-plane crossing number `c_tri = c(T_1) + c(T_2) + c(T_3)` and concentration flag.
- Direct tri-plane algebraic topology:
  - `triplane_fundamental_group`: exact Seifert-van Kampen presentation of `pi_1(S^4 \\ F)` on `b` arch
    generators and `2b` cap relations (orientable, non-orientable, and multi-component surface-links),
  - `triplane_double_branched_cover`: exact `Sigma_2(S^4, F)` Betti/torsion/signature invariants via the
    affine reflection module `{-1, +1} |>< Z^b`,
  - `compute_triplane_invariant_record`: multivariable Alexander ideals `E_1, E_2`, `Sigma_2`, finite
    quotient hom counts, and `S_3`-twisted homology directly from `(T_1, T_2, T_3)`.
- `concentrated_triplane_to_marked_diagram`: exact zero-extra-crossing compiler from concentrated
  tri-plane diagrams `(M_1, T_2, M_3)` with non-interleaving meander cycles to planar `MarkedDiagram`s
  of `ch = c(T_2) + b - p_31`.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

from surfaceknots.diagrams import (
    ClassicalLinkDiagram,
    ComponentTopology,
    MarkedDiagram,
    SurfaceTopology,
    VertexKind,
)
from surfaceknots.invariants.alexander import compute_alexander_invariants, twisted_s3_invariants
from surfaceknots.invariants.fundamental_group import (
    GroupPresentation,
    Word,
    cyclically_reduce_word,
    finite_quotient_profile,
    invert_word,
    reduce_word,
    tietze_simplify,
)
from surfaceknots.invariants.tier3 import DoubleBranchedCoverResult
from surfaceknots.unlink import certify_classical_unlink


@dataclass(frozen=True)
class TangleDiagram:
    """A `b`-bridge tangle diagram in a disk with `2b` boundary points `0..2b-1` along the top line.

    Represented combinatorially by:
    - `num_bridges`: `b >= 1`,
    - `bottom_matching`: planar non-crossing pairing `tuple[tuple[int, int], ...]` of `0..2b-1`
      at the bottom of the tangle,
    - `braid_word`: tuple of signed Artin generators `+-k` (`1 <= k <= 2b - 1`) read from bottom
      to top (`+k` puts strand `k-1` OVER strand `k` from bottom-left to top-right; `-k` puts
      strand `k-1` UNDER strand `k`).
    """

    num_bridges: int
    bottom_matching: tuple[tuple[int, int], ...]
    braid_word: tuple[int, ...] = ()

    @property
    def num_crossings(self) -> int:
        return len(self.braid_word)


def _trace_bottom_origins(t: TangleDiagram) -> list[int]:
    """Return `wire_origin[0..2b-1]` where `wire_origin[x]` is the bottom port arriving at top port `x`."""
    w_count = 2 * t.num_bridges
    wire_origin = list(range(w_count))
    for g in t.braid_word:
        k = abs(g) - 1
        wire_origin[k], wire_origin[k + 1] = wire_origin[k + 1], wire_origin[k]
    return wire_origin


def tangle_endpoint_pairing(t: TangleDiagram) -> tuple[int, ...]:
    """Return the fixed-point-free involution `partner[0..2b-1]` induced by `t`'s arcs on `0..2b-1`."""
    w_count = 2 * t.num_bridges
    wire_origin = _trace_bottom_origins(t)
    origin_to_top = [0] * w_count
    for top_x, bot_i in enumerate(wire_origin):
        origin_to_top[bot_i] = top_x
    partner = [0] * w_count
    for u, v in t.bottom_matching:
        tx, ty = origin_to_top[u], origin_to_top[v]
        partner[tx] = ty
        partner[ty] = tx
    return tuple(partner)


def union_two_tangles_to_classical_link(
    t_top: TangleDiagram,
    t_bot: TangleDiagram,
    top_port_scomp: tuple[int, ...] | None = None,
) -> ClassicalLinkDiagram:
    """Glue `t_top` and `mirror(t_bot)` along their `2b` boundary points `0..2b-1` to form a
    closed `ClassicalLinkDiagram` `D = T_top union bar(T_bot)`.
    """
    if t_top.num_bridges != t_bot.num_bridges:
        raise ValueError("Bridge counts must match")
    b = t_top.num_bridges
    w_count = 2 * b

    bot_scomp: list[int] | None = None
    if top_port_scomp is not None:
        orig_bot = _trace_bottom_origins(t_bot)
        bot_scomp = [0] * w_count
        for top_x, bot_i in enumerate(orig_bot):
            bot_scomp[bot_i] = top_port_scomp[top_x]

    combined_braid = list(t_bot.braid_word) + [-g for g in reversed(t_top.braid_word)]
    c = len(combined_braid)

    if c == 0:
        adj_bot = {}
        for u, v in t_bot.bottom_matching:
            adj_bot[u] = v
            adj_bot[v] = u
        adj_top = {}
        for u, v in t_top.bottom_matching:
            adj_top[u] = v
            adj_top[v] = u
        seen = set()
        loops = 0
        fl_scomps: list[int] = []
        for s in range(w_count):
            if s in seen:
                continue
            loops += 1
            if top_port_scomp is not None:
                fl_scomps.append(top_port_scomp[s])
            cur = s
            while cur not in seen:
                seen.add(cur)
                nxt = adj_top[cur]
                seen.add(nxt)
                cur = adj_bot[nxt]
        return ClassicalLinkDiagram(
            num_crossings=0,
            opp=(),
            vtype=(),
            free_loops=loops,
            dart_surface_comp=() if top_port_scomp is not None else None,
            free_loop_surface_comp=tuple(fl_scomps) if top_port_scomp is not None else None,
        )

    vtypes: list[str] = []
    wire_cur: list[tuple[str, int]] = [("bot", i) for i in range(w_count)]
    wire_sc: list[int] = list(bot_scomp) if bot_scomp is not None else [0] * w_count
    dart_sc: list[int] = [0] * (4 * c)
    opp_map: dict[tuple[str, int], tuple[str, int]] = {}

    def connect(p1: tuple[str, int], p2: tuple[str, int]) -> None:
        opp_map[p1] = p2
        opp_map[p2] = p1

    for v, g in enumerate(combined_braid):
        k = abs(g) - 1
        vtypes.append("X1" if g > 0 else "X0")
        p_left = wire_cur[k]
        p_right = wire_cur[k + 1]
        connect(p_left, ("dart", 4 * v + 0))
        connect(p_right, ("dart", 4 * v + 1))
        wire_cur[k] = ("dart", 4 * v + 3)
        wire_cur[k + 1] = ("dart", 4 * v + 2)
        if bot_scomp is not None:
            sc_l, sc_r = wire_sc[k], wire_sc[k + 1]
            dart_sc[4 * v + 0] = sc_l
            dart_sc[4 * v + 2] = sc_l
            dart_sc[4 * v + 1] = sc_r
            dart_sc[4 * v + 3] = sc_r
            wire_sc[k], wire_sc[k + 1] = sc_r, sc_l

    for u, w in t_top.bottom_matching:
        connect(wire_cur[u], wire_cur[w])

    bot_partner = {}
    for u, w in t_bot.bottom_matching:
        bot_partner[("bot", u)] = ("bot", w)
        bot_partner[("bot", w)] = ("bot", u)

    opp = [0] * (4 * c)
    visited_bots = set()
    for d in range(4 * c):
        cur = opp_map[("dart", d)]
        while cur[0] == "bot":
            visited_bots.add(cur)
            nxt_b = bot_partner[cur]
            visited_bots.add(nxt_b)
            cur = opp_map[nxt_b]
        opp[d] = cur[1]

    free_loops = 0
    fl_scomps_c: list[int] = []
    for u in range(w_count):
        bp = ("bot", u)
        if bp in visited_bots:
            continue
        free_loops += 1
        if bot_scomp is not None:
            fl_scomps_c.append(bot_scomp[u])
        cur = bp
        while cur not in visited_bots:
            visited_bots.add(cur)
            nxt_b = bot_partner[cur]
            visited_bots.add(nxt_b)
            cur = opp_map[nxt_b]

    return ClassicalLinkDiagram(
        num_crossings=c,
        opp=tuple(opp),
        vtype=tuple(vtypes),  # type: ignore[arg-type]
        free_loops=free_loops,
        dart_surface_comp=tuple(dart_sc) if bot_scomp is not None else None,
        free_loop_surface_comp=tuple(fl_scomps_c) if bot_scomp is not None else None,
    )


@dataclass(frozen=True)
class TriPlaneDiagram:
    """A `(b; p_12, p_23, p_31)` tri-plane diagram `(T_1, T_2, T_3)` in `R^4` (Meier-Zupan 2017)."""

    t1: TangleDiagram
    t2: TangleDiagram
    t3: TangleDiagram

    @property
    def bridge_number(self) -> int:
        assert self.t1.num_bridges == self.t2.num_bridges == self.t3.num_bridges
        return self.t1.num_bridges

    @property
    def crossing_number(self) -> int:
        return self.t1.num_crossings + self.t2.num_crossings + self.t3.num_crossings

    @property
    def crossing_triple(self) -> tuple[int, int, int]:
        return (self.t1.num_crossings, self.t2.num_crossings, self.t3.num_crossings)

    @property
    def is_concentrated(self) -> bool:
        """True iff at least two of the three tangles have 0 crossings."""
        return sum(1 for t in (self.t1, self.t2, self.t3) if t.num_crossings == 0) >= 2

    def pairwise_links(
        self,
        top_port_scomp: tuple[int, ...] | None = None,
    ) -> tuple[ClassicalLinkDiagram, ClassicalLinkDiagram, ClassicalLinkDiagram]:
        return (
            union_two_tangles_to_classical_link(self.t1, self.t2, top_port_scomp=top_port_scomp),
            union_two_tangles_to_classical_link(self.t2, self.t3, top_port_scomp=top_port_scomp),
            union_two_tangles_to_classical_link(self.t3, self.t1, top_port_scomp=top_port_scomp),
        )

    def is_certified_valid(self) -> bool:
        """Certify that all three pairwise tangle unions `T_i union bar(T_j)` are classical unlinks."""
        for l_ij in self.pairwise_links():
            if not certify_classical_unlink(l_ij).is_unlink:
                return False
        return True

    def patch_counts(self) -> tuple[int, int, int]:
        l12, l23, l31 = self.pairwise_links()
        return (l12.num_circles(), l23.num_circles(), l31.num_circles())

    def euler_characteristic(self) -> int:
        p12, p23, p31 = self.patch_counts()
        return p12 + p23 + p31 - self.bridge_number

    def _port_orbits_and_orientability(self) -> tuple[tuple[int, ...], list[list[int]], list[bool]]:
        w_count = 2 * self.bridge_number
        tau1 = tangle_endpoint_pairing(self.t1)
        tau2 = tangle_endpoint_pairing(self.t2)
        tau3 = tangle_endpoint_pairing(self.t3)

        port_scomp = [-1] * w_count
        port_sign = [0] * w_count
        orbits: list[list[int]] = []
        orientable_flags: list[bool] = []

        for start_p in range(w_count):
            if port_scomp[start_p] != -1:
                continue
            sc = len(orbits)
            orb: list[int] = []
            is_ori = True
            port_scomp[start_p] = sc
            port_sign[start_p] = 1
            stack = [start_p]
            while stack:
                p = stack.pop()
                orb.append(p)
                expected = -port_sign[p]
                for tau in (tau1, tau2, tau3):
                    q = tau[p]
                    if port_scomp[q] == -1:
                        port_scomp[q] = sc
                        port_sign[q] = expected
                        stack.append(q)
                    elif port_sign[q] != expected:
                        is_ori = False
            orbits.append(sorted(orb))
            orientable_flags.append(is_ori)

        return tuple(port_scomp), orbits, orientable_flags

    def surface_topology(self) -> SurfaceTopology:
        """Compute exact `SurfaceTopology` (components, orientability, `chi`, and normal Euler `e(F_c)`)
        directly from `(T_1, T_2, T_3)` via bridge-point bipartiteness and Joseph-Meier-Miller-Zupan's
        tri-plane writhe formula `e(F_c) = -(w_self(L_12|_c) + w_self(L_23|_c) + w_self(L_31|_c))`.
        """
        port_scomp, orbits, orientable_flags = self._port_orbits_and_orientability()
        num_sc = len(orbits)
        tau1 = tangle_endpoint_pairing(self.t1)
        tau2 = tangle_endpoint_pairing(self.t2)
        tau3 = tangle_endpoint_pairing(self.t3)

        def count_meander_cycles(ta: tuple[int, ...], tb: tuple[int, ...]) -> list[int]:
            counts = [0] * num_sc
            seen = [False] * (2 * self.bridge_number)
            for p in range(2 * self.bridge_number):
                if seen[p]:
                    continue
                counts[port_scomp[p]] += 1
                cur = p
                while not seen[cur]:
                    seen[cur] = True
                    nxt = ta[cur]
                    seen[nxt] = True
                    cur = tb[nxt]
            return counts

        p12 = count_meander_cycles(tau1, tau2)
        p23 = count_meander_cycles(tau2, tau3)
        p31 = count_meander_cycles(tau3, tau1)

        l12, l23, l31 = self.pairwise_links(top_port_scomp=port_scomp)
        w12 = l12.self_writhe_by_surface_component(num_sc)
        w23 = l23.self_writhe_by_surface_component(num_sc)
        w31 = l31.self_writhe_by_surface_component(num_sc)

        comps: list[ComponentTopology] = []
        for sc in range(num_sc):
            b_sc = len(orbits[sc]) // 2
            chi_sc = p12[sc] + p23[sc] + p31[sc] - b_sc
            ori = orientable_flags[sc]
            e_num = -(w12[sc] + w23[sc] + w31[sc])
            if ori:
                g_sym = (2 - chi_sc) // 2
            else:
                g_sym = -(2 - chi_sc)
            comps.append(
                ComponentTopology(
                    orientable=ori,
                    euler_characteristic=chi_sc,
                    genus_symbol=g_sym,
                    normal_euler=e_num,
                )
            )

        return SurfaceTopology(
            components=tuple(comps),
            num_circles_plus=sum(p12),
            num_circles_minus=sum(p23),
        )


def _propagate_braid_words(bottom_words: list[Word], braid_word: tuple[int, ...]) -> list[Word]:
    """Propagate Wirtinger meridian words `(w_0, ..., w_{2b-1})` upward through `braid_word`."""
    words = list(bottom_words)
    for g in braid_word:
        if g > 0:
            k = g - 1
            u, v = words[k], words[k + 1]
            words[k] = reduce_word(u + v + invert_word(u))
            words[k + 1] = u
        else:
            k = -g - 1
            u, v = words[k], words[k + 1]
            words[k] = v
            words[k + 1] = reduce_word(invert_word(v) + u + v)
    return words


def triplane_fundamental_group(tp: TriPlaneDiagram, simplify: bool = True) -> GroupPresentation:
    """Compute an exact presentation of `pi_1(S^4 \\ F(tp))` via Seifert-van Kampen on `H_1 u H_2 u H_3`.

    Uses `b` generators (one per bottom arch of `T_1`) and `2b` cap relations (`b` from `T_2`, `b` from `T_3`),
    working identically for orientable, non-orientable, and multi-component surface-links.
    """
    b = tp.bridge_number
    w_count = 2 * b
    port_scomp, orbits, _ = tp._port_orbits_and_orientability()
    num_sc = len(orbits)

    orig_t1 = _trace_bottom_origins(tp.t1)
    bot_t1_scomp = [0] * w_count
    for top_x, bot_i in enumerate(orig_t1):
        bot_t1_scomp[bot_i] = port_scomp[top_x]

    w_t1_bot: list[Word] = [()] * w_count
    comp_meridian: list[Word | None] = [None] * num_sc
    for r, (u, v) in enumerate(tp.t1.bottom_matching):
        u_lo, v_hi = (u, v) if u < v else (v, u)
        gen_1based = r + 1
        w_t1_bot[u_lo] = (gen_1based,)
        w_t1_bot[v_hi] = (-gen_1based,)
        sc = bot_t1_scomp[u_lo]
        if comp_meridian[sc] is None:
            comp_meridian[sc] = (gen_1based,)

    m_top = _propagate_braid_words(w_t1_bot, tp.t1.braid_word)

    relators: list[Word] = []
    for ti in (tp.t2, tp.t3):
        inv_braid = tuple(-g for g in reversed(ti.braid_word))
        w_bot = _propagate_braid_words(m_top, inv_braid)
        for u, v in ti.bottom_matching:
            u_lo, v_hi = (u, v) if u < v else (v, u)
            rel = cyclically_reduce_word(w_bot[u_lo] + w_bot[v_hi])
            if rel:
                relators.append(rel)

    pres = GroupPresentation(
        num_generators=b,
        relators=tuple(relators),
        component_meridians=tuple(m if m is not None else () for m in comp_meridian),
    )
    return tietze_simplify(pres) if simplify else pres


def triplane_double_branched_cover(
    tp: TriPlaneDiagram,
    topo: SurfaceTopology | None = None,
) -> DoubleBranchedCoverResult:
    """Compute exact 4-manifold invariants of the mod-2 double branched cover `Sigma_2(S^4, F(tp))`
    directly from `tp` via the affine reflection module `{-1, +1} |>< Z^b` and G-signature.
    """
    if topo is None:
        topo = tp.surface_topology()
    b = tp.bridge_number
    raw_pres = triplane_fundamental_group(tp, simplify=False)

    rows: list[list[int]] = []
    for rel in raw_pres.relators:
        sgn = 1
        vec = [0] * b
        for lit in rel:
            j = abs(lit) - 1
            vec[j] += sgn
            sgn = -sgn
        assert sgn == 1, f"Odd-length relator in pi_1(S^4 \\ F): {rel}"
        if any(vec):
            rows.append(vec)

    if not rows:
        r_kei = b
        t_kei: tuple[int, ...] = ()
    else:
        snf = smith_normal_form(sp.Matrix(rows), domain=sp.ZZ)
        diag = [abs(int(snf[i, i])) for i in range(min(snf.rows, snf.cols))]
        nonzero = [d for d in diag if d > 0]
        r_kei = b - len(nonzero)
        t_kei = tuple(sorted(d for d in nonzero if d > 1))

    chi_f = topo.euler_characteristic
    total_e = sum(c.normal_euler for c in topo.components)
    chi_sigma2 = 4 - chi_f
    b1 = r_kei - 1
    b2 = 2 + 2 * b1 - chi_f
    sigma = -total_e // 2
    b2_plus = (b2 + sigma) // 2
    b2_minus = (b2 - sigma) // 2
    return DoubleBranchedCoverResult(
        euler_characteristic=chi_sigma2,
        b1=b1,
        h1_torsion=t_kei,
        b2=b2,
        signature=sigma,
        abs_signature=abs(sigma),
        b2_plus=b2_plus,
        b2_minus=b2_minus,
        is_rational_homology_4sphere=(b1 == 0 and b2 == 0),
    )


def compute_triplane_invariant_record(tp: TriPlaneDiagram, label: str = "") -> dict[str, object]:
    """Compute a complete topological, group, Alexander, and branched-cover invariant record on `tp`."""
    import hashlib

    topo = tp.surface_topology()
    pres = triplane_fundamental_group(tp, simplify=True)
    h1_rank, h1_torsion = pres.abelianization_invariant_factors()
    alex = compute_alexander_invariants(pres=pres, topo=topo)
    dbc = triplane_double_branched_cover(tp, topo=topo)
    hom_counts = finite_quotient_profile(pres)
    tw_s3 = twisted_s3_invariants(None, pres=pres)

    comp_types = [c.homeomorphism_type for c in topo.components]
    normal_eulers = sorted(c.normal_euler for c in topo.components)
    e1_fp_hex = hashlib.sha256(repr(alex.evaluation_fingerprint).encode("utf-8")).hexdigest()[:16]

    compiled_dlk2: list[list[object]] | None = None
    if tp.is_concentrated:
        compiled_md = concentrated_triplane_to_marked_diagram(tp)
        if compiled_md is not None:
            from surfaceknots.invariants.linking import double_linking_mod2

            compiled_dlk2 = [list(x) for x in double_linking_mod2(compiled_md)]

    return {
        "label": label,
        "bridge_number": tp.bridge_number,
        "c_tri": tp.crossing_number,
        "crossing_triple": list(tp.crossing_triple),
        "is_concentrated": tp.is_concentrated,
        "patch_counts": list(tp.patch_counts()),
        "num_components": topo.num_components,
        "yoshikawa_genera": list(topo.yoshikawa_genera),
        "component_types": comp_types,
        "orientability": topo.orientability_label,
        "euler_characteristic": topo.euler_characteristic,
        "normal_euler_numbers": normal_eulers,
        "abs_normal_euler_multiset": sorted(abs(x) for x in normal_eulers),
        "pi1_generators": pres.num_generators,
        "pi1_relators": len(pres.relators),
        "pi1_presentation": pres.to_human_string(),
        "h1_free_rank": h1_rank,
        "h1_torsion": list(h1_torsion),
        "alexander_variables": list(alex.variable_names),
        "e1_generators": list(alex.e1_generators),
        "e2_generators": list(alex.e2_generators),
        "single_var_alex": alex.single_var_alex,
        "e1_fingerprint_hash": e1_fp_hex,
        "hom_counts": dict(hom_counts),
        "twisted_s3_delta0": [str(x) for x in tw_s3],
        "double_linking_mod2": compiled_dlk2,
        "sigma2_euler": dbc.euler_characteristic,
        "sigma2_b1": dbc.b1,
        "sigma2_h1_torsion": list(dbc.h1_torsion),
        "sigma2_b2": dbc.b2,
        "sigma2_signature": dbc.signature,
        "sigma2_abs_signature": dbc.abs_signature,
        "sigma2_b2_unoriented_pair": [min(dbc.b2_plus, dbc.b2_minus), max(dbc.b2_plus, dbc.b2_minus)],
        "sigma2_is_q_homology_4sphere": dbc.is_rational_homology_4sphere,
    }


def _canonicalize_dlk2_o4(dlk_raw: object) -> tuple:
    """Canonicalize `double_linking_mod2` up to global `O(4)` mirror `e -> -e`."""
    import re

    if not isinstance(dlk_raw, (list, tuple)):
        return ()

    def _flip_e(s: str) -> str:
        return re.sub(r"e=(-?\d+)", lambda m: f"e={-int(m.group(1))}", s)

    pos_rows = tuple(sorted((str(r[0]), str(r[1]), int(r[2])) for r in dlk_raw))
    neg_rows = tuple(sorted((_flip_e(str(r[0])), _flip_e(str(r[1])), int(r[2])) for r in dlk_raw))
    return min(pos_rows, neg_rows)


def triplane_o4_key(rec: dict[str, object], *, include_double_linking: bool = False) -> tuple:
    """`O(4)`-invariant equivalence key compatible with both `compute_triplane_invariant_record`
    and `MarkedDiagram` / JSON catalog records (`yoshikawa_canonical.json`, `expanded_ch11.json`,
    `expanded_ch12.json`, `family_catalog.json`, `targeted_ch12_13_census.json`).

    Invariant boundary note:
    - By default (`include_double_linking=False`), returns the 6-tuple key over invariants defined
      directly on every (concentrated or unconcentrated) `TriPlaneDiagram` without requiring a planar
      `MarkedDiagram` compilation. Across `b <= 5`, this 6-tuple is completely free of the `(8_1^1,1,
      10_1^1,1)` `T^2 u T^2` collision because every `T^2` component requires `>= 3` bridges (`chi(T^2)
      = 0` cannot be realized with `<= 2` bridges), forcing `b(T^2 u T^2) >= 6`.
    - When `include_double_linking=True` (e.g. comparing `b >= 6` concentrated compilable tri-plane
      records against `MarkedDiagram` records where `rec.get("double_linking_mod2")` is non-`None` on
      both sides), appends `O(4)`-canonicalized Carter-Saito-Satoh `double_linking_mod2`, separating
      `8_1^1,1` (`dlk_2 = 0`) from `10_1^1,1` (`dlk_2 = 1`).
    """
    e_list: list[int] = [int(x) for x in rec["normal_euler_numbers"]]  # type: ignore[union-attr]
    if int(rec["num_components"]) == 1:  # type: ignore[arg-type]
        euler_key = tuple(sorted(abs(x) for x in e_list))
    elif sum(e_list) == 0:
        euler_key = tuple(sorted(e_list))
    else:
        euler_key = tuple(min(sorted(e_list), sorted(-x for x in e_list)))

    base = (
        (
            int(rec["num_components"]),  # type: ignore[arg-type]
            tuple(int(g) for g in rec["yoshikawa_genera"]),  # type: ignore[union-attr]
            str(rec["orientability"]),
            int(rec["euler_characteristic"]),  # type: ignore[arg-type]
            euler_key,
        ),
        (int(rec["h1_free_rank"]), tuple(int(t) for t in rec["h1_torsion"])),  # type: ignore[arg-type, union-attr]
        str(rec["e1_fingerprint_hash"]),
        (
            int(rec["sigma2_euler"]),  # type: ignore[arg-type]
            int(rec["sigma2_b1"]),  # type: ignore[arg-type]
            tuple(int(t) for t in rec["sigma2_h1_torsion"]),  # type: ignore[union-attr]
            int(rec["sigma2_b2"]),  # type: ignore[arg-type]
            int(rec["sigma2_abs_signature"]),  # type: ignore[arg-type]
        ),
        tuple(sorted((str(k), int(v)) for k, v in rec["hom_counts"].items())),  # type: ignore[union-attr]
        tuple(str(x) for x in rec["twisted_s3_delta0"]),  # type: ignore[union-attr]
    )
    if include_double_linking and rec.get("double_linking_mod2") is not None:
        return base + (_canonicalize_dlk2_o4(rec["double_linking_mod2"]),)
    return base


def concentrated_triplane_to_marked_diagram(tp: TriPlaneDiagram) -> MarkedDiagram | None:
    """Compile a concentrated `TriPlaneDiagram` (at least two crossingless tangles) into a planar
    `MarkedDiagram` with **zero extra classical crossings** (`ch(D) = c(tp) + b - p_cap`) whenever
    the two crossingless cap matchings have non-interleaving meander cycles on `S^1`.

    Returns `None` if `tp` is not concentrated or if the two cap matchings exhibit the meander
    interleaving obstruction in `D_+` (requiring classical R2 crossings to project into a single plane).
    """
    if not tp.is_concentrated:
        return None

    # Cyclically permute `(t1, t2, t3)` (even permutation -> preserves S^4 orientation & normal Euler sign!)
    # so `t1_c` and `t3_c` have 0 crossings and `t2_c` holds all classical crossings:
    if tp.t1.num_crossings == 0 and tp.t3.num_crossings == 0:
        t1_c, t2_c, t3_c = tp.t1, tp.t2, tp.t3
    elif tp.t2.num_crossings == 0 and tp.t1.num_crossings == 0:
        t1_c, t2_c, t3_c = tp.t2, tp.t3, tp.t1
    else:
        assert tp.t3.num_crossings == 0 and tp.t2.num_crossings == 0
        t1_c, t2_c, t3_c = tp.t3, tp.t1, tp.t2

    b = tp.bridge_number
    w_count = 2 * b

    # Partner maps for `M_minus = t3_c.bottom_matching` (realized at sgn = -1) and
    # `M_plus = t1_c.bottom_matching` (realized at sgn = +1):
    m_minus: dict[int, int] = {}
    for u, v in t3_c.bottom_matching:
        m_minus[u] = v
        m_minus[v] = u
    m_plus: dict[int, int] = {}
    for u, v in t1_c.bottom_matching:
        m_plus[u] = v
        m_plus[v] = u

    vtypes: list[VertexKind] = []
    wire_cur: list[tuple[str, int]] = [("bot", i) for i in range(w_count)]
    opp_map: dict[tuple[str, int], tuple[str, int]] = {}

    def connect(p1: tuple[str, int], p2: tuple[str, int]) -> None:
        opp_map[p1] = p2
        opp_map[p2] = p1

    # 1. Build `t2_c.braid_word` from bottom `("bot", 0..2b-1)` up to equator `wire_cur[0..2b-1]`:
    for g in t2_c.braid_word:
        v = len(vtypes)
        k = abs(g) - 1
        vtypes.append("X1" if g > 0 else "X0")
        p_left = wire_cur[k]
        p_right = wire_cur[k + 1]
        connect(p_left, ("dart", 4 * v + 0))
        connect(p_right, ("dart", 4 * v + 1))
        wire_cur[k] = ("dart", 4 * v + 3)
        wire_cur[k + 1] = ("dart", 4 * v + 2)

    # 2. Reduce `(m_minus, m_plus)` in the upper disk `D_+` on circular active port list `active`:
    active = list(range(w_count))
    port_end: dict[int, tuple[str, int]] = {i: wire_cur[i] for i in range(w_count)}
    next_virt_id = w_count

    while active:
        k_act = len(active)
        # Rule 1: Common cyclic-adjacent arch `(p_a, p_b)` in BOTH `m_minus` and `m_plus`
        found_common = -1
        for idx in range(k_act):
            pa = active[idx]
            pb = active[(idx + 1) % k_act]
            if m_minus[pa] == pb and m_plus[pa] == pb:
                found_common = idx
                break
        if found_common != -1:
            pa = active[found_common]
            pb = active[(found_common + 1) % k_act]
            connect(port_end[pa], port_end[pb])
            del m_minus[pa], m_minus[pb], m_plus[pa], m_plus[pb]
            if found_common + 1 < k_act:
                del active[found_common + 1]
                del active[found_common]
            else:
                active.pop()
                active.pop(0)
            continue

        # Rule 2: Shared-endpoint cyclic-adjacent triple `(p0, p1, p2)` (`k_act >= 4`)
        if k_act < 4:
            return None
        found_triple: tuple[int, VertexKind] | None = None
        for idx in range(k_act):
            p0 = active[idx]
            p1 = active[(idx + 1) % k_act]
            p2 = active[(idx + 2) % k_act]
            if m_minus[p0] == p1 and m_plus[p1] == p2:
                found_triple = (idx, "Y")
                break
            if m_plus[p0] == p1 and m_minus[p1] == p2:
                found_triple = (idx, "Z")
                break

        if found_triple is None:
            # Meander interleaving obstruction in `D_+`!
            return None

        idx, vt = found_triple
        p0 = active[idx]
        p1 = active[(idx + 1) % k_act]
        p2 = active[(idx + 2) % k_act]
        v = len(vtypes)
        vtypes.append(vt)
        connect(port_end[p0], ("dart", 4 * v + 0))
        connect(port_end[p1], ("dart", 4 * v + 1))
        connect(port_end[p2], ("dart", 4 * v + 2))

        merged_id = next_virt_id
        next_virt_id += 1
        port_end[merged_id] = ("dart", 4 * v + 3)

        if vt == "Y":
            # sgn = -1 caps (p0, p1) and sends p2 -> merged_id;
            # sgn = +1 caps (p1, p2) and sends p0 -> merged_id.
            q_minus = m_minus[p2]
            q_plus = m_plus[p0]
            del m_minus[p0], m_minus[p1], m_minus[p2]
            del m_plus[p0], m_plus[p1], m_plus[p2]
        else:
            # vt == "Z": sgn = -1 caps (p1, p2) and sends p0 -> merged_id;
            #            sgn = +1 caps (p0, p1) and sends p2 -> merged_id.
            q_minus = m_minus[p0]
            q_plus = m_plus[p2]
            del m_minus[p0], m_minus[p1], m_minus[p2]
            del m_plus[p0], m_plus[p1], m_plus[p2]

        m_minus[merged_id] = q_minus
        m_minus[q_minus] = merged_id
        m_plus[merged_id] = q_plus
        m_plus[q_plus] = merged_id

        # Replace `(p0, p1, p2)` at cyclic positions `idx, idx+1, idx+2` with `[merged_id]`:
        rem_set = {(idx) % k_act, (idx + 1) % k_act, (idx + 2) % k_act}
        new_active: list[int] = []
        inserted = False
        for j_pos in range(k_act):
            if j_pos == idx:
                new_active.append(merged_id)
                inserted = True
            elif j_pos not in rem_set:
                new_active.append(active[j_pos])
        assert inserted and len(new_active) == k_act - 2
        active = new_active

    # 3. Splice bottom caps `t2_c.bottom_matching` to extract `opp` and `free_loops`:
    n_verts = len(vtypes)
    bot_partner = {}
    for u, w in t2_c.bottom_matching:
        bot_partner[("bot", u)] = ("bot", w)
        bot_partner[("bot", w)] = ("bot", u)

    opp = [0] * (4 * n_verts)
    visited_bots = set()
    for d in range(4 * n_verts):
        cur = opp_map[("dart", d)]
        while cur[0] == "bot":
            visited_bots.add(cur)
            nxt_b = bot_partner[cur]
            visited_bots.add(nxt_b)
            cur = opp_map[nxt_b]
        opp[d] = cur[1]

    free_loops = 0
    for u in range(w_count):
        bp = ("bot", u)
        if bp in visited_bots:
            continue
        free_loops += 1
        cur = bp
        while cur not in visited_bots:
            visited_bots.add(cur)
            nxt_b = bot_partner[cur]
            visited_bots.add(nxt_b)
            cur = opp_map[nxt_b]

    return MarkedDiagram(
        num_vertices=n_verts,
        opp=tuple(opp),
        vtype=tuple(vtypes),
        free_loops=free_loops,
    )


@lru_cache(maxsize=16)
def generate_catalan_matchings(b: int) -> tuple[tuple[tuple[int, int], ...], ...]:
    """Generate all `C_b = (2b)! / ((b+1)! b!)` planar non-crossing matchings of `0..2b-1`."""
    if b == 0:
        return ((),)

    def rec(pts: tuple[int, ...]) -> list[list[tuple[int, int]]]:
        if not pts:
            return [[]]
        first = pts[0]
        res: list[list[tuple[int, int]]] = []
        for idx in range(1, len(pts), 2):
            partner = pts[idx]
            inside = pts[1:idx]
            outside = pts[idx + 1 :]
            for m_in in rec(inside):
                for m_out in rec(outside):
                    res.append([(first, partner)] + m_in + m_out)
        return res

    raw = rec(tuple(range(2 * b)))
    return tuple(tuple(sorted((min(u, v), max(u, v)) for u, v in m)) for m in raw)


def standard_unknot_s2_triplane() -> TriPlaneDiagram:
    """Minimal `(1; 1, 1, 1)` crossingless tri-plane diagram of `0_1` (`S^2`, `b = 1, c_tri = 0`)."""
    t = TangleDiagram(1, ((0, 1),), ())
    return TriPlaneDiagram(t, t, t)


def standard_rp2_triplane(sign: int = 1) -> TriPlaneDiagram:
    """Minimal `(2; 1, 1, 1)` 1-crossing tri-plane diagram of `2_1^-1` (`RP^2`, `b = 2, c_tri = 1`)."""
    m_a = ((0, 1), (2, 3))
    m_b = ((0, 3), (1, 2))
    t1 = TangleDiagram(2, m_a, (2 if sign >= 0 else -2,))
    t2 = TangleDiagram(2, m_a, ())
    t3 = TangleDiagram(2, m_b, ())
    return TriPlaneDiagram(t1, t2, t3)


def standard_torus_triplane() -> TriPlaneDiagram:
    """Minimal `(3; 1, 1, 1)` ZERO-crossing tri-plane diagram of `2_1^1` (`T^2`, `b = 3, c_tri = 0`)."""
    t1 = TangleDiagram(3, ((0, 3), (1, 2), (4, 5)), ())
    t2 = TangleDiagram(3, ((0, 5), (1, 4), (2, 3)), ())
    t3 = TangleDiagram(3, ((0, 1), (2, 5), (3, 4)), ())
    return TriPlaneDiagram(t1, t2, t3)


def estimate_bridge_bounds_from_marked(d: MarkedDiagram) -> tuple[int, int]:
    """Return rigorous `(b_lower, b_upper)` bounds on bridge number `b(F(D))` from topology + `D`."""
    topo = d.surface_topology()
    c = topo.num_components
    chi = topo.euler_characteristic
    b_low = max(c, 3 * c - chi)
    if d.num_vertices >= 6:
        b_low = max(b_low, 4)

    if d.num_vertices == 0 and d.free_loops == 1:
        return (1, 1)
    if d.num_vertices == 2:
        return (b_low, b_low)

    c_minus = topo.num_circles_minus
    c_plus = topo.num_circles_plus
    m = d.num_marked
    c_cl = d.num_classical
    b_up = max(b_low, min(c_minus, c_plus) + m + (c_cl + 1) // 2)
    return (b_low, b_up)
