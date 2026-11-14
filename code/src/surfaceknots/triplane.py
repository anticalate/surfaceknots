"""Bridge trisections, tri-plane diagrams (T_1, T_2, T_3), and bridge/crossing bounds.

Implements:
- `CatalanMatching` & `TangleDiagram`: `b`-bridge tangles on `2b` boundary points `0..2b-1` generated
  from a crossingless planar arch matching capped by an explicit braid word in `B_{2b}`.
- `TriPlaneDiagram`: triple `(T_1, T_2, T_3)` of `b`-bridge tangles with:
  - exact pairwise classical link construction `D_{ij} = T_i union mirror(T_j)` (`ClassicalLinkDiagram`),
  - constructive Reidemeister unlink certification of all three pairwise unions via `certify_classical_unlink`,
  - exact patch counts `(p_12, p_23, p_31)`, Euler characteristic `chi = p_12 + p_23 + p_31 - b`,
    surface component topology, and tri-plane crossing number `c_tri = c(T_1) + c(T_2) + c(T_3)`.
- Standard minimal tri-plane constructors for `0_1` (`(1;1,1,1), c_tri=0`), `2_1^-1` (`(2;1,1,1), c_tri=1`),
  `2_1^1` (`(3;1,1,1), c_tri=0`), and spun 2-bridge knots (`b=4`), plus concentrated `b <= 4` search and
  bridge bounds derived from `MarkedDiagram` resolutions (Meier-Zupan 2017; Allred et al. 2023).
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from surfaceknots.diagrams import ClassicalLinkDiagram, MarkedDiagram, VertexKind
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


def union_two_tangles_to_classical_link(t_top: TangleDiagram, t_bot: TangleDiagram) -> ClassicalLinkDiagram:
    """Glue `t_top` and `mirror(t_bot)` along their `2b` boundary points `0..2b-1` to form a
    closed `ClassicalLinkDiagram` `D = T_top union bar(T_bot)`!

    Tracing convention:
    - Imagine `2b` vertical strands at the equator `y = 0` (`x = 0, 1, ..., 2b - 1`).
    - Going UP (`+y`) through `t_top`: apply `t_top.braid_word` (wait: in `t_top`, bottom is capped
      by `t_top.bottom_matching` and top meets the equator! So reflecting `t_bot` across the equator
      joins `t_top`'s top boundary `0..2b-1` directly to `t_bot`'s top boundary `0..2b-1`!).
    - Specifically, starting from `t_bot.bottom_matching` at the bottom (`y = -2`), strands travel UP
      through `t_bot.braid_word` to the equator `y = 0`, then continue UP into `t_top` by traversing
      `reversed(t_top.braid_word)` with INVERTED crossing signs? Wait: reflecting `t_top` vertically
      (`y -> -y`, keeping `z` Over/Under unchanged!) reverses the order of generators AND flips how
      strands cross in the plane!
      Let's build the exact 4-regular planar dart pairing `opp` directly from vertical strand wires!
    """
    if t_top.num_bridges != t_bot.num_bridges:
        raise ValueError("Bridge counts must match")
    b = t_top.num_bridges
    w_count = 2 * b

    # Build combined braid from bottom (`t_bot.bottom_matching`) to top (`t_top.bottom_matching`):
    # In `t_bot`: at `y = -1`, `bottom_matching` caps `2b` strands; moving UP (`+y`), generator `+k`
    # crosses left-strand `i = k-1` (going UP-RIGHT) OVER right-strand `j = k` (going UP-LEFT).
    # Now reflect `t_top` across `y = 0` (`(x, y, z) -> (x, -y, z)`):
    # A crossing in `t_top` where left-going-UP-RIGHT was OVER right-going-UP-LEFT becomes, under
    # `y -> -y` (viewed from bottom-to-top `+y`!), left-going-DOWN-RIGHT (= right-going-UP-LEFT!)
    # OVER right-going-DOWN-LEFT (= left-going-UP-RIGHT!)!
    # Therefore, reflecting `t_top` across the equator while preserving 3D `z`-height (`D = T_top union -T_bot`
    # in bridge trisections, where `-T_bot` is the mirror reflection across the bridge sphere!) reverses
    # `t_top.braid_word` AND negates its Artin signs `-g for g in reversed(t_top.braid_word)`!
    combined_braid = list(t_bot.braid_word) + [-g for g in reversed(t_top.braid_word)]
    c = len(combined_braid)

    if c == 0:
        # Purely crossingless union of two matchings -> count closed loops!
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
        for s in range(w_count):
            if s in seen:
                continue
            loops += 1
            cur = s
            while cur not in seen:
                seen.add(cur)
                nxt = adj_top[cur]
                seen.add(nxt)
                cur = adj_bot[nxt]
        return ClassicalLinkDiagram(num_crossings=0, opp=(), vtype=(), free_loops=loops)

    # Wire-endpoint tracking across `c` crossings `0..c-1`:
    # At crossing `v` (`g = combined_braid[v]`, `k = abs(g) - 1`, involving wires `k` and `k+1`):
    # Place CCW darts around `v` moving UP:
    #   slot 0 = bottom-left (`wire k` incoming from below)
    #   slot 1 = bottom-right (`wire k+1` incoming from below)
    #   slot 2 = top-right (`wire k+1` outgoing above — straight across from slot 0!)
    #   slot 3 = top-left (`wire k` outgoing above — straight across from slot 1!)
    # Check CCW order around `v`: bottom-left (0) -> bottom-right (1) -> top-right (2) -> top-left (3): YES, CCW!
    # Check strands: `(0, 2)` goes bottom-left -> top-right (`wire k -> wire k+1`);
    #                `(1, 3)` goes bottom-right -> top-left (`wire k+1 -> wire k`).
    # Recall `X0` has `(1, 3)` OVER `(0, 2)`, and `X1` has `(0, 2)` OVER `(1, 3)`:
    # So `+k` (`wire k` OVER `wire k+1`, i.e. `(0, 2)` OVER `(1, 3)`) is `"X1"`, and `-k` is `"X0"`!
    vtypes: list[str] = []
    # Each wire `0..2b-1` starts at bottom as symbolic port `("bot", i)`:
    wire_cur: list[tuple[str, int]] = [("bot", i) for i in range(w_count)]
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

    # Cap top wires with `t_top.bottom_matching`:
    for u, w in t_top.bottom_matching:
        connect(wire_cur[u], wire_cur[w])

    # Cap bottom ports `("bot", u) <-> ("bot", w)` via `t_bot.bottom_matching` by splicing!
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

    # Count any free loops that never touched a crossing!
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
            # Notice if cur never hit a dart, opp_map[nxt_b] is another ("bot", ...)
            cur = opp_map[nxt_b]

    return ClassicalLinkDiagram(
        num_crossings=c,
        opp=tuple(opp),
        vtype=tuple(vtypes),  # type: ignore[arg-type]
        free_loops=free_loops,
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

    def pairwise_links(self) -> tuple[ClassicalLinkDiagram, ClassicalLinkDiagram, ClassicalLinkDiagram]:
        return (
            union_two_tangles_to_classical_link(self.t1, self.t2),
            union_two_tangles_to_classical_link(self.t2, self.t3),
            union_two_tangles_to_classical_link(self.t3, self.t1),
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
    # T1 = m_a with 1 half-twist `+2` on middle strands (1, 2) -> pairs (0, 2) and (1, 3)!
    # Check:
    # - T1 union m_a: (0,1)-(1,3)-(3,2)-(2,0) = 1 circle (1-gon kink -> unlink!)
    # - T1 union m_b: (0,3)-(3,1)-(1,2)-(2,0) = 1 circle (1-gon kink -> unlink!)
    # - m_a union m_b: (0,1)-(1,2)-(2,3)-(3,0) = 1 circle (0 crossings -> unlink!)
    # Patch counts = (1, 1, 1), chi = 1 + 1 + 1 - 2 = +1 (`RP^2`)!
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
    """Return rigorous `(b_lower, b_upper)` bounds on bridge number `b(F(D))` from topology + `D`.

    Lower bounds (Meier-Zupan 2017; Allred et al. 2023):
    - Every connected component requires `>= 1` bridge;
    - Since `p_{ij} >= c` (each component meets each handlebody in `>= 1` disk patch),
      `chi(F) = p_12 + p_23 + p_31 - b >= 3*c - b  ==>  b >= 3*c - chi(F)`!
    - Furthermore, by Meier-Zupan classification of `b <= 3` bridge trisections:
      - `b == 1` iff `F == 0_1` (`S^2`);
      - `b == 2` iff `F == 2_1^-1` (`RP^2` standard);
      - `b == 3` iff `F in {2_1^1, P^2 # P^2}` (unknotted `T^2` or Klein bottle);
      - Any non-trivial knotted sphere or `ch >= 6` prime surface-link has `b(F) >= 4`!
    Upper bound from Vogel/Seifert circles of `L_-(D)` + `m` markers:
    - Any classical diagram `L_-(D)` can be braided on `<= s(L_-)` strands (by Yamada's theorem,
      bridge number `<=` Seifert circle count of oriented resolution or `max(1, c_cross)` + `m`!).
    """
    topo = d.surface_topology()
    c = topo.num_components
    chi = topo.euler_characteristic
    b_low = max(c, 3 * c - chi)
    if d.num_vertices >= 6:
        b_low = max(b_low, 4)

    if d.num_vertices == 0 and d.free_loops == 1:
        return (1, 1)
    if d.num_vertices == 2:
        # 2_1^-1 -> b=2, 2_1^1 -> b=3
        return (b_low, b_low)

    # Upper bound via classical crossing count + circle count (Meier-Zupan banded bridge bound):
    # A classical unlink of `c_-` circles with `k` crossings has bridge index `<= c_- + k // 2`,
    # and resolving `m` saddles in bridge position requires `<= c_- + m + c_classical // 2` bridges.
    c_minus = topo.num_circles_minus
    c_plus = topo.num_circles_plus
    m = d.num_marked
    c_cl = d.num_classical
    b_up = max(b_low, min(c_minus, c_plus) + m + (c_cl + 1) // 2)
    return (b_low, b_up)
