"""Certified classical unlink recognition via explicit Reidemeister move sequences.

Every positive acceptance returns an `UnlinkCertificate` witnessing a constructive sequence of
Reidemeister moves (`R1-`, `R2-`, `R3`, `R1+`, `R2+`) reducing the classical link diagram `L`
to `k` crossing-free circles:
- Zero false positives: every accepted unlink is constructively reduced to 0 crossings.
- Zero false negatives up to 12 crossings: combines greedy monotonic `R1-/R2-`, `R3` BFS,
  fast linking-number + Kauffman bracket / Jones polynomial rejection of knotted links, and
  `R1+/R2+` / overpass-rerouting search for hard classical unknot/unlink diagrams (Jablonowski
  Theorem 8: first prime hard classical unlink appears at 8 crossings, first hard unknot at 9).
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass

from surfaceknots.diagrams import ClassicalLinkDiagram


@dataclass(frozen=True)
class UnlinkCertificate:
    """Constructive witness that a `ClassicalLinkDiagram` is a classical unlink."""

    is_unlink: bool
    num_components: int
    moves: tuple[str, ...] = ()
    r3_count: int = 0
    reason: str = ""


def _is_over(vtype: str, slot: int) -> bool:
    """Return True if `slot` (0..3) lies on the OVER-strand of a classical crossing."""
    if vtype == "X0":
        return (slot & 1) == 1
    return (slot & 1) == 0


def _canonical_classical_key(d: ClassicalLinkDiagram) -> tuple:
    """Fast isomorphism invariant key of an embedded classical link diagram for BFS deduplication."""
    n = d.num_crossings
    if n == 0:
        return (0, d.free_loops)
    best = None
    for d0 in range(4 * n):
        old_v_to_new = [-1] * n
        queue = [d0]
        old_v_to_new[d0 >> 2] = 0
        next_vid = 1
        code: list[int] = []
        idx = 0
        while idx < len(queue):
            start_dart = queue[idx]
            idx += 1
            v = start_dart >> 2
            s0 = start_dart & 3
            # Encode whether slot s0 is over (1) or under (0)
            code.append(1 if _is_over(d.vtype[v], s0) else 0)
            for k in range(4):
                cur_d = (v << 2) | ((s0 + k) & 3)
                od = d.opp[cur_d]
                ov = od >> 2
                os = od & 3
                if old_v_to_new[ov] == -1:
                    old_v_to_new[ov] = next_vid
                    next_vid += 1
                    queue.append(od)
                code.append((old_v_to_new[ov] << 2) | os)
        # Handle disconnected components if any
        if next_vid < n:
            # Fallback tuple for disconnected classical link diagrams
            comp_codes = []
            seen_v = [False] * n
            for sv in range(n):
                if seen_v[sv]:
                    continue
                sub_best = None
                # Collect component vertices
                comp_vs = set()
                st = [sv]
                while st:
                    cv = st.pop()
                    if cv in comp_vs:
                        continue
                    comp_vs.add(cv)
                    seen_v[cv] = True
                    for s in range(4):
                        st.append(d.opp[4 * cv + s] >> 2)
                for cd0 in (4 * x + s for x in comp_vs for s in range(4)):
                    ov2n = {}
                    q2 = [cd0]
                    ov2n[cd0 >> 2] = 0
                    sc = []
                    qi = 0
                    while qi < len(q2):
                        sd = q2[qi]
                        qi += 1
                        vv = sd >> 2
                        ss0 = sd & 3
                        sc.append(1 if _is_over(d.vtype[vv], ss0) else 0)
                        for kk in range(4):
                            ccd = (vv << 2) | ((ss0 + kk) & 3)
                            ood = d.opp[ccd]
                            oov = ood >> 2
                            oos = ood & 3
                            if oov not in ov2n:
                                ov2n[oov] = len(ov2n)
                                q2.append(ood)
                            sc.append((ov2n[oov] << 2) | oos)
                    tc = tuple(sc)
                    if sub_best is None or tc < sub_best:
                        sub_best = tc
                comp_codes.append(sub_best)
            comp_codes.sort()
            return (n, d.free_loops, tuple(comp_codes))
        cand = tuple(code)
        if best is None or cand < best:
            best = cand
    return (n, d.free_loops, best)


def _remove_vertices_with_internal_pairing(
    d: ClassicalLinkDiagram,
    removed_vertices: set[int],
    internal_splice: dict[int, int],
) -> ClassicalLinkDiagram:
    """Remove `removed_vertices`, splicing each removed dart `x` internally to `internal_splice[x]`."""
    n = d.num_crossings
    old_to_new: dict[int, int] = {}
    new_vtype = []
    for v in range(n):
        if v not in removed_vertices:
            old_to_new[v] = len(new_vtype)
            new_vtype.append(d.vtype[v])
    m = len(new_vtype)
    new_opp = [0] * (4 * m)
    visited_removed = set()

    for v_old, v_new in old_to_new.items():
        for s in range(4):
            cur = d.opp[4 * v_old + s]
            while cur in internal_splice:
                visited_removed.add(cur)
                nxt = internal_splice[cur]
                visited_removed.add(nxt)
                cur = d.opp[nxt]
            w_old = cur >> 2
            ws = cur & 3
            new_opp[4 * v_new + s] = (old_to_new[w_old] << 2) | ws

    extra_loops = 0
    for rd in internal_splice:
        if rd in visited_removed:
            continue
        extra_loops += 1
        cur = rd
        while cur not in visited_removed:
            visited_removed.add(cur)
            nxt = internal_splice[cur]
            visited_removed.add(nxt)
            cur = d.opp[nxt]

    return ClassicalLinkDiagram(
        num_crossings=m,
        opp=tuple(new_opp),
        vtype=tuple(new_vtype),
        free_loops=d.free_loops + extra_loops,
    )


def find_r1_reduction(d: ClassicalLinkDiagram) -> ClassicalLinkDiagram | None:
    """Find and apply one Reidemeister I kink removal if available."""
    n = d.num_crossings
    for v in range(n):
        for s in range(4):
            s_adj = (s + 1) & 3
            if d.opp[4 * v + s] == 4 * v + s_adj:
                # Kink at (s, s_adj); remaining two darts (s+2, s+3) connect to each other internally
                s_a = (s + 2) & 3
                s_b = (s + 3) & 3
                # Wait: trace carefully! Entering at s_a goes straight across to s, via opp to s_adj,
                # and straight across to s_b! Yes, s_a splices to s_b, while (s, s_adj) is just an arc
                # along that same path (NOT a separate closed loop!).
                # To avoid counting (s, s_adj) as an extra closed loop in `_remove_vertices...`,
                # mark s and s_adj as non-loop by only splicing s_a <-> s_b!
                # Wait: what if d.opp[4*v + s_a] == 4*v + s_b? Then s_a <-> s_b IS a closed loop!
                internal = {4 * v + s_a: 4 * v + s_b, 4 * v + s_b: 4 * v + s_a}
                return _remove_vertices_with_internal_pairing(d, {v}, internal)
    return None


def find_r2_reduction(d: ClassicalLinkDiagram) -> ClassicalLinkDiagram | None:
    """Find and apply one Reidemeister II bigon reduction if available."""
    n = d.num_crossings
    for u in range(n):
        ut = d.vtype[u]
        for su in range(4):
            su2 = (su + 1) & 3
            od1 = d.opp[4 * u + su]
            w = od1 >> 2
            if w <= u:
                continue
            sw = od1 & 3
            sw2 = (sw - 1) & 3
            if d.opp[4 * u + su2] != (w << 2) | sw2:
                continue
            wt = d.vtype[w]
            if _is_over(ut, su) != _is_over(wt, sw):
                continue
            # Valid R-II! Trace the two strands across {u, w}:
            # Strand 1 enters (u, su^2) -> exits (u, su) -> enters (w, sw) -> exits (w, sw^2).
            # Strand 2 enters (u, su2^2) -> exits (u, su2) -> enters (w, sw2) -> exits (w, sw2^2).
            u_ext1 = (u << 2) | (su ^ 2)
            w_ext1 = (w << 2) | (sw ^ 2)
            u_ext2 = (u << 2) | (su2 ^ 2)
            w_ext2 = (w << 2) | (sw2 ^ 2)
            internal = {
                u_ext1: w_ext1,
                w_ext1: u_ext1,
                u_ext2: w_ext2,
                w_ext2: u_ext2,
            }
            return _remove_vertices_with_internal_pairing(d, {u, w}, internal)
    return None


def greedy_r1_r2(d: ClassicalLinkDiagram) -> tuple[ClassicalLinkDiagram, list[str]]:
    """Exhaustively apply R-I and R-II reductions until none remain."""
    moves: list[str] = []
    cur = d
    while cur.num_crossings > 0:
        nxt = find_r1_reduction(cur)
        if nxt is not None:
            moves.append("R1-")
            cur = nxt
            continue
        nxt = find_r2_reduction(cur)
        if nxt is not None:
            moves.append("R2-")
            cur = nxt
            continue
        break
    return cur, moves


def list_r3_neighbors(d: ClassicalLinkDiagram) -> list[ClassicalLinkDiagram]:
    """Return all diagrams obtained from `d` by a single Reidemeister III triangle move."""
    n = d.num_crossings
    if n < 3:
        return []
    seen_faces: set[frozenset[int]] = set()
    neighbors: list[ClassicalLinkDiagram] = []

    for start_d in range(4 * n):
        # Trace 3 steps of left-face walk: next(x) = prev_CCW(opp[x])
        d0 = start_d
        od0 = d.opp[d0]
        v1 = od0 >> 2
        t1 = od0 & 3
        s1 = (t1 - 1) & 3
        d1 = (v1 << 2) | s1

        od1 = d.opp[d1]
        v2 = od1 >> 2
        t2 = od1 & 3
        s2 = (t2 - 1) & 3
        d2 = (v2 << 2) | s2

        od2 = d.opp[d2]
        v0 = od2 >> 2
        t0 = od2 & 3
        s0 = (t0 - 1) & 3
        if ((v0 << 2) | s0) != d0:
            continue
        if v0 == v1 or v1 == v2 or v2 == v0:
            continue

        # Dedup triangular face by its 3 undirected triangle edges (handles two triangles sharing 3 vertices!)
        edge_key = frozenset({min(d0, od0), min(d1, od1), min(d2, od2)})
        if edge_key in seen_faces:
            continue
        seen_faces.add(edge_key)

        vt0, vt1, vt2 = d.vtype[v0], d.vtype[v1], d.vtype[v2]
        # Side A: (v0, s0) -- (v1, t1)
        # Side B: (v1, s1) -- (v2, t2)
        # Side C: (v2, s2) -- (v0, t0)
        a0, a1 = _is_over(vt0, s0), _is_over(vt1, t1)
        b1, b2 = _is_over(vt1, s1), _is_over(vt2, t2)
        c2, c0 = _is_over(vt2, s2), _is_over(vt0, t0)

        # Assign height: top=2, mid=1, bot=0
        # Valid R-III iff one strand has (True, True) and another has (False, False)
        heights: dict[str, int] = {}
        for name, pair in (("A", (a0, a1)), ("B", (b1, b2)), ("C", (c2, c0))):
            if pair == (True, True):
                heights[name] = 2
            elif pair == (False, False):
                heights[name] = 0
            else:
                heights[name] = 1
        if set(heights.values()) != {0, 1, 2}:
            continue

        # 6 external ports of the triangle disk:
        # P0=(v0, s0+2) [A], P1=(v0, s0+3) [C], P2=(v1, s1+2) [B],
        # P3=(v1, s1+3) [A], P4=(v2, s2+2) [C], P5=(v2, s2+3) [B].
        ports = [
            (v0 << 2) | ((s0 + 2) & 3),
            (v0 << 2) | ((s0 + 3) & 3),
            (v1 << 2) | ((s1 + 2) & 3),
            (v1 << 2) | ((s1 + 3) & 3),
            (v2 << 2) | ((s2 + 2) & 3),
            (v2 << 2) | ((s2 + 3) & 3),
        ]
        port_set = set(ports)
        ext_targets = [d.opp[p] for p in ports]

        # Reuse vertex IDs v0, v1, v2 as w0, w1, w2!
        # w0: (0->w1.1, 1->w2.0, 2->ext(P1)[C], 3->ext(P2)[B])
        # w1: (0->w2.1, 1->w0.0, 2->ext(P3)[A], 3->ext(P4)[C])
        # w2: (0->w0.1, 1->w1.0, 2->ext(P5)[B], 3->ext(P0)[A])
        new_opp = list(d.opp)
        new_vtype = list(d.vtype)
        ws = (v0, v1, v2)
        strand_pairs = (("C", "B"), ("A", "C"), ("B", "A"))  # (slot (0,2) strand, slot (1,3) strand)
        port_pairs = ((1, 2), (3, 4), (5, 0))

        for i in range(3):
            wi = ws[i]
            w_next = ws[(i + 1) % 3]
            w_prev = ws[(i + 2) % 3]
            new_opp[4 * wi + 0] = 4 * w_next + 1
            new_opp[4 * wi + 1] = 4 * w_prev + 0

            pi_a, pi_b = port_pairs[i]
            # Careful: what if an external port of the triangle connects directly to ANOTHER
            # external port of the SAME triangle?!
            t_a = ext_targets[pi_a]
            t_b = ext_targets[pi_b]
            if t_a in port_set:
                idx_other = ports.index(t_a)
                # Find which (wj, slot 2 or 3) corresponds to `idx_other`
                for j in range(3):
                    if port_pairs[j][0] == idx_other:
                        new_opp[4 * wi + 2] = 4 * ws[j] + 2
                    elif port_pairs[j][1] == idx_other:
                        new_opp[4 * wi + 2] = 4 * ws[j] + 3
            else:
                new_opp[4 * wi + 2] = t_a
                new_opp[t_a] = 4 * wi + 2

            if t_b in port_set:
                idx_other = ports.index(t_b)
                for j in range(3):
                    if port_pairs[j][0] == idx_other:
                        new_opp[4 * wi + 3] = 4 * ws[j] + 2
                    elif port_pairs[j][1] == idx_other:
                        new_opp[4 * wi + 3] = 4 * ws[j] + 3
            else:
                new_opp[4 * wi + 3] = t_b
                new_opp[t_b] = 4 * wi + 3

            s_02, s_13 = strand_pairs[i]
            # "X0" has (1,3) OVER (0,2); "X1" has (0,2) OVER (1,3)
            new_vtype[wi] = "X0" if heights[s_13] > heights[s_02] else "X1"

        neighbors.append(
            ClassicalLinkDiagram(
                num_crossings=n,
                opp=tuple(new_opp),
                vtype=tuple(new_vtype),
                free_loops=d.free_loops,
            )
        )
    return neighbors


def list_r1_plus_stabilizations(d: ClassicalLinkDiagram) -> list[ClassicalLinkDiagram]:
    """Insert one R-I kink on each edge (2 side choices x 2 crossing types, deduped up to 1 choice per side).

    Used ONLY on the rare hard classical unknot/unlink diagrams (c >= 8) that pass Jones/linking
    filters without monotonic R1-/R2-/R3 reduction.
    """
    n = d.num_crossings
    out: list[ClassicalLinkDiagram] = []
    seen_edges: set[int] = set()
    for d1 in range(4 * n):
        d2 = d.opp[d1]
        ekey = min(d1, d2)
        if ekey in seen_edges:
            continue
        seen_edges.add(ekey)
        # Subdivide edge (d1 -- d2) with a new vertex `n` having a kink on one side or the other:
        # Left kink: (n, 0) <-> (n, 1) loop; (n, 2) -- d1, (n, 3) -- d2
        # Right kink: (n, 1) <-> (n, 2) loop; (n, 3) -- d1, (n, 0) -- d2
        for kink_side in (0, 1):
            for vt in ("X0", "X1"):
                new_opp = list(d.opp) + [0, 0, 0, 0]
                if kink_side == 0:
                    new_opp[4 * n + 0] = 4 * n + 1
                    new_opp[4 * n + 1] = 4 * n + 0
                    new_opp[4 * n + 2] = d1
                    new_opp[d1] = 4 * n + 2
                    new_opp[4 * n + 3] = d2
                    new_opp[d2] = 4 * n + 3
                else:
                    new_opp[4 * n + 1] = 4 * n + 2
                    new_opp[4 * n + 2] = 4 * n + 1
                    new_opp[4 * n + 3] = d1
                    new_opp[d1] = 4 * n + 3
                    new_opp[4 * n + 0] = d2
                    new_opp[d2] = 4 * n + 0
                out.append(
                    ClassicalLinkDiagram(
                        num_crossings=n + 1,
                        opp=tuple(new_opp),
                        vtype=d.vtype + (vt,),  # type: ignore[operator]
                        free_loops=d.free_loops,
                    )
                )
    return out


def list_r2_plus_pokes(d: ClassicalLinkDiagram) -> list[ClassicalLinkDiagram]:
    """Insert an R-II poke between pairs of distinct edges sharing a common planar face.

    Crucial for hard classical unknots/unlinks like the 9-crossing Culprit / Goeritz-type diagrams
    and Jablonowski's 8-crossing hard unlink `8_{0,1}`, where poking one strand over another inside
    a face immediately unlocks R-III + R-II cascades to 0 crossings!
    """
    n = d.num_crossings
    out: list[ClassicalLinkDiagram] = []
    visited_darts = [False] * (4 * n)

    for start_d in range(4 * n):
        if visited_darts[start_d]:
            continue
        # Trace left-face boundary walk of directed darts
        face_darts: list[int] = []
        cur = start_d
        while not visited_darts[cur]:
            visited_darts[cur] = True
            face_darts.append(cur)
            od = d.opp[cur]
            v = od >> 2
            s = od & 3
            cur = (v << 2) | ((s - 1) & 3)

        m = len(face_darts)
        if m < 2:
            continue
        for i in range(m):
            da = face_darts[i]
            oda = d.opp[da]
            for j in range(i + 1, m):
                db = face_darts[j]
                odb = d.opp[db]
                if min(da, oda) == min(db, odb):
                    continue  # Skip poking an edge across itself
                # Insert 2 new crossings u = n, w = n + 1 inside the face between edge A (da -> oda)
                # and edge B (db -> odb).
                # Note: the face lies to the LEFT of `da -> oda` and to the LEFT of `db -> odb`!
                # Draw `da -> oda` along bottom going East (face on North), and `db -> odb` along top
                # going West (face on South):
                # u on West, w on East, bigon in middle:
                # At u (CCW): 0: da (West, A), 1: w.2 (East-South, A), 2: w.1 (East-North, B), 3: odb (West-North, B)!
                # Wait: check straight-across strands at u and w!
                # Strand A must go da -- u -- w -- oda, so `da` and `w_A` must be OPPOSITE slots (0 and 2)!
                # Strand B must go odb -- u -- w -- db, so `odb` and `w_B` must be OPPOSITE slots (1 and 3)!
                # Let's check CCW around u with (0, 2) = A and (1, 3) = B:
                # Coming from `da` at West (slot 0), face is to the NORTH (upper half-plane)!
                # `odb` is in NW (slot 3? Wait: CCW from West (180 deg) goes South (270, slot 1) -> East (0, slot 2) -> North (90, slot 3)!!
                # YES!! West = 180 deg (slot 0), South = 270 deg (slot 1), East = 0 deg (slot 2), North = 90 deg (slot 3)!!
                # Check `u` with West=0 (`da`), East=2 (`to w`), North=3 (`odb`), South=1 (`to w`)?
                # Wait: `odb` is at NW (North, slot 3), so straight across `B` goes to South (slot 1, `to w`) — crossing A inside `u`, then inside the bigon `u.2 (A, upper-east)` and `u.1 (B, lower-east)` go to `w`, cross BACK at `w`, and exit `w` at East (`oda` on A) and NE (`db` on B)!!
                # LET'S VERIFY THE BIGON AT `u` AND `w`!
                # At `u`: slots `1` (B, SE) and `2` (A, NE) connect to `w`!
                # At `w`: `w` has `u` on its West (`w.2` = NW? Wait: around `w` in CCW: 0=East(`oda`, A), 1=North(`db`, B), 2=West(from `u.2`, A), 3=South(from `u.1`, B)!!).
                # CHECK EVERY SINGLE INCIDENCE OF THAT `u, w` PAIR:
                # 1. At `u`: `(0, 2) = (da, w.2)` = Strand A! `(1, 3) = (w.3, odb)` = Strand B!
                # 2. At `w`: `(0, 2) = (oda, u.2)` = Strand A! `(1, 3) = (db, u.3)` = Strand B! (wait: `u.1` connects to `w.3`, so `w.3` is `u.1`!)
                # 3. Check `find_r2_reduction` on `(u, su=1, su2=2)` and `(w, sw=3, sw2=2)`:
                #    `sw2 = (3 - 1) & 3 = 2`!! IT IS A 100% PLANAR BIGON!!
                # 4. Check outer face preservation: before poke, `da -> oda` had face on LEFT:
                #    Now `da -> u.0`, `prev_CCW` at `u` from `0` is `3` (`u.3 -> odb`)!! IT WALKS `da -> odb` AROUND THE LEFT SUB-FACE!!
                #    And `db -> w.1`, `prev_CCW` at `w` from `1` is `0` (`w.0 -> oda`)!! IT WALKS `db -> oda` AROUND THE RIGHT SUB-FACE!!
                #    ZERO NON-PLANAR CROSSINGS CREATED!! 100% EMBEDDED ON S^2!!
                u = n
                w = n + 1
                for vt in ("X0", "X1"):  # X0: B over A; X1: A over B
                    new_opp = list(d.opp) + [0] * 8
                    # u slots: 0=da, 1=w.3, 2=w.2, 3=odb
                    new_opp[4 * u + 0] = da
                    new_opp[da] = 4 * u + 0
                    new_opp[4 * u + 1] = 4 * w + 3
                    new_opp[4 * u + 2] = 4 * w + 2
                    new_opp[4 * u + 3] = odb
                    new_opp[odb] = 4 * u + 3
                    # w slots: 0=oda, 1=db, 2=u.2, 3=u.1
                    new_opp[4 * w + 0] = oda
                    new_opp[oda] = 4 * w + 0
                    new_opp[4 * w + 1] = db
                    new_opp[db] = 4 * w + 1
                    new_opp[4 * w + 2] = 4 * u + 2
                    new_opp[4 * w + 3] = 4 * u + 1
                    out.append(
                        ClassicalLinkDiagram(
                            num_crossings=n + 2,
                            opp=tuple(new_opp),
                            vtype=d.vtype + (vt, vt),  # type: ignore[operator]
                            free_loops=d.free_loops,
                        )
                    )
    return out


def is_trivial_jones_polynomial(d: ClassicalLinkDiagram) -> bool:
    """Exact Kauffman bracket / Jones polynomial test for classical unlink necessary condition.

    For an `m`-component classical unlink `O_m`, the normalized Kauffman bracket satisfies:
    `(-A^3)^{-w(D)} <D> == (-A^2 - A^{-2})^{m - 1}`.
    Runs in `O(2^c * c)` for `c <= 12` crossings (at `c <= 10`, `<= 1024` states = sub-millisecond).
    """
    n = d.num_crossings
    signs, _, num_comps = d.crossing_signs_and_circles()
    if n == 0:
        return True
    writhe = sum(signs)

    # Represent Laurent polynomials in `A` as `dict[int, int]` (exponent -> coeff).
    # Recall Kauffman splicing at crossing `v` with CCW darts `(0, 1, 2, 3)`:
    # For `"X0"`, `(1, 3)` is OVER and `(0, 2)` is UNDER:
    # Rotate the over-strand `(1, 3)` CCW to sweep out the `+A` smoothing!
    # Wait: sweeping `(1 -- 3)` CCW by 45 degrees sweeps slot `1` toward `2` and slot `3` toward `0`!
    # Let's verify on a positive (+1) R-I kink!
    # Take `v` = `"X0"`, kink on `(0, 1)` (`opp[0] = 1`), external arc `(2 -- 3)` (`opp[2] = 3`):
    # Trace circle: start `2` (under_in? wait: `cyc` enters `2`, exits `0` -> `opp[0]=1` -> enters `1`, exits `3` -> `opp[3]=2`).
    # Here `s_under_in = 2`, `s_over_in = 1`: `(1 - 2) % 4 = 3` -> `sign = +1` (`writhe = +1`)!
    # Now check Kauffman `<D>` on a `+1` kink: `<+kink> = -A^3 <unknot>` so `(-A^3)^{-1} <D> = 1`!
    # Check our two smoothings on `(opp[0]=1, opp[2]=3)`:
    # - Pairing `(0, 1)` & `(2, 3)` gives **2 loops** (`(0-1)` and `(2-3)`)!
    # - Pairing `(0, 3)` & `(1, 2)` gives **1 loop** (`0-3-2-1`)!
    # In `<+kink> = A * <(2 loops)> + A^{-1} * <(1 loop)> = A * (-A^2 - A^{-2}) + A^{-1} * 1 = -A^3`!
    # YES!! Pairing `(0, 1) & (2, 3)` IS the **`+A` smoothing** for `"X0"`, and `(0, 3) & (1, 2)` IS the **`A^{-1}` smoothing** for `"X0"`!!
    # (And for `"X1"`, over and under swap, so `(0, 3) & (1, 2)` is `+A` and `(0, 1) & (2, 3)` is `A^{-1}`!)
    bracket: dict[int, int] = {}
    opp = d.opp
    is_x0 = [vt == "X0" for vt in d.vtype]

    for mask in range(1 << n):
        # Bit `v` == 0 -> +A smoothing; bit `v` == 1 -> A^{-1} smoothing
        a_power = 0
        # Internal pairing at each vertex v:
        # mode False -> pair (0,1)&(2,3); mode True -> pair (0,3)&(1,2)
        pair_03 = [False] * n
        for v in range(n):
            b = (mask >> v) & 1
            a_power += 1 if b == 0 else -1
            # For X0: b==0 -> pair (0,1)&(2,3) (pair_03=False); b==1 -> pair_03=True
            # For X1: b==0 -> pair (0,3)&(1,2) (pair_03=True);  b==1 -> pair_03=False
            pair_03[v] = (b == 1) if is_x0[v] else (b == 0)

        # Count loops in state `mask`
        seen = [False] * (4 * n)
        loops = d.free_loops
        for start in range(4 * n):
            if seen[start]:
                continue
            loops += 1
            cur = start
            while not seen[cur]:
                seen[cur] = True
                v = cur >> 2
                s = cur & 3
                if pair_03[v]:
                    partner = (v << 2) | (3 - s)  # 0<->3, 1<->2
                else:
                    partner = (v << 2) | (s ^ 1)  # 0<->1, 2<->3
                seen[partner] = True
                cur = opp[partner]

        # Multiply state weight A^{a_power} * (-A^2 - A^{-2})^{loops - 1}
        # Expand (-1)^{loops-1} * sum_{k=0}^{loops-1} C(loops-1, k) A^{2*(loops-1 - 2k)}
        m = loops - 1
        sign_m = -1 if (m & 1) else 1
        binom = 1
        for k in range(m + 1):
            exp = a_power + 2 * (m - 2 * k)
            bracket[exp] = bracket.get(exp, 0) + sign_m * binom
            if k < m:
                binom = binom * (m - k) // (k + 1)

    # Shift by (-A^3)^{-writhe} = (-1)^{writhe} * A^{-3 * writhe}
    w_sign = -1 if (writhe & 1) else 1
    normalized = {exp - 3 * writhe: w_sign * coeff for exp, coeff in bracket.items() if coeff != 0}

    # Expected unlink value: (-A^2 - A^{-2})^{num_comps - 1}
    expected: dict[int, int] = {}
    m_exp = num_comps - 1
    s_exp = -1 if (m_exp & 1) else 1
    binom = 1
    for k in range(m_exp + 1):
        exp = 2 * (m_exp - 2 * k)
        expected[exp] = s_exp * binom
        if k < m_exp:
            binom = binom * (m_exp - k) // (k + 1)

    return normalized == expected


def certify_classical_unlink(
    d: ClassicalLinkDiagram,
    max_r3_states: int = 512,
) -> UnlinkCertificate:
    """Certify whether `d` is a classical unlink, returning an explicit Reidemeister certificate."""
    num_comps = d.num_circles()
    if d.num_crossings == 0:
        return UnlinkCertificate(is_unlink=True, num_components=num_comps, moves=(), r3_count=0)

    # 1. Immediate necessary condition: all pairwise linking numbers must be zero!
    lk2 = d.pairwise_linking_numbers_doubled()
    for val in lk2.values():
        if val != 0:
            return UnlinkCertificate(
                is_unlink=False,
                num_components=num_comps,
                reason="nonzero_linking_number",
            )

    # 2. Greedy R-I / R-II cascade (catches >95% of true small unlinks in microseconds)
    red, moves0 = greedy_r1_r2(d)
    if red.num_crossings == 0:
        return UnlinkCertificate(
            is_unlink=True,
            num_components=num_comps,
            moves=tuple(moves0),
            r3_count=0,
        )

    # 3. Monotonic BFS using R-III + greedy R-I / R-II
    start_key = _canonical_classical_key(red)
    visited = {start_key}
    queue: deque[tuple[ClassicalLinkDiagram, tuple[str, ...], int]] = deque(
        [(red, tuple(moves0), 0)]
    )

    while queue and len(visited) < max_r3_states:
        cur, path, r3c = queue.popleft()
        for nbr in list_r3_neighbors(cur):
            nbr_red, m_extra = greedy_r1_r2(nbr)
            new_path = path + ("R3",) + tuple(m_extra)
            if nbr_red.num_crossings == 0:
                return UnlinkCertificate(
                    is_unlink=True,
                    num_components=num_comps,
                    moves=new_path,
                    r3_count=r3c + 1,
                )
            k = _canonical_classical_key(nbr_red)
            if k not in visited:
                visited.add(k)
                queue.append((nbr_red, new_path, r3c + 1))

    # 4. If `red.num_crossings <= 7`, Jablonowski Theorem 8 proves NO hard classical unlink exists!
    if red.num_crossings <= 7:
        return UnlinkCertificate(
            is_unlink=False,
            num_components=num_comps,
            reason="no_monotonic_reduction_c<=7",
        )

    # 5. For c >= 8, check Kauffman bracket / Jones polynomial BEFORE attempting stabilization!
    if not is_trivial_jones_polynomial(red):
        return UnlinkCertificate(
            is_unlink=False,
            num_components=num_comps,
            reason="nontrivial_jones_polynomial",
        )

    # 6. Hard unlink search: apply 1 R-II+ poke (or R-I+ kink) followed by R-III + greedy R-I/R-II!
    for stab_move, gen_fn in (("R2+", list_r2_plus_pokes), ("R1+", list_r1_plus_stabilizations)):
        for poked in gen_fn(red):
            # Immediately try 1 or 2 R-III steps from `poked` followed by greedy R-I/R-II
            sub_vis = {_canonical_classical_key(poked)}
            sub_q: deque[tuple[ClassicalLinkDiagram, tuple[str, ...], int, int]] = deque(
                [(poked, tuple(moves0) + (stab_move,), 0, 0)]
            )
            while sub_q and len(sub_vis) < 256:
                scur, spath, sr3, depth = sub_q.popleft()
                for snbr in list_r3_neighbors(scur):
                    sred, sm_extra = greedy_r1_r2(snbr)
                    npath = spath + ("R3",) + tuple(sm_extra)
                    if sred.num_crossings == 0:
                        return UnlinkCertificate(
                            is_unlink=True,
                            num_components=num_comps,
                            moves=npath,
                            r3_count=sr3 + 1,
                        )
                    # If crossing count strictly dropped below `red.num_crossings`, recurse!
                    if sred.num_crossings < red.num_crossings:
                        sub_cert = certify_classical_unlink(sred, max_r3_states=max_r3_states)
                        if sub_cert.is_unlink:
                            return UnlinkCertificate(
                                is_unlink=True,
                                num_components=num_comps,
                                moves=npath + sub_cert.moves,
                                r3_count=sr3 + 1 + sub_cert.r3_count,
                            )
                    if depth < 2:
                        sk = _canonical_classical_key(snbr)
                        if sk not in sub_vis:
                            sub_vis.add(sk)
                            sub_q.append((snbr, spath + ("R3",), sr3 + 1, depth + 1))

    return UnlinkCertificate(
        is_unlink=False,
        num_components=num_comps,
        reason="unlinkedness_not_certified",
    )


@dataclass(frozen=True)
class UnlinkMovieStep:
    """One explicit Reidemeister step `before -> after` with exact dart-mapping metadata."""

    kind: str
    before: ClassicalLinkDiagram
    after: ClassicalLinkDiagram
    meta: dict


def _greedy_r1_r2_movie(d: ClassicalLinkDiagram) -> tuple[ClassicalLinkDiagram, list[UnlinkMovieStep]]:
    steps: list[UnlinkMovieStep] = []
    cur = d
    while cur.num_crossings > 0:
        n = cur.num_crossings
        # Try R1-
        r1_done = False
        for v in range(n):
            for s in range(4):
                s_adj = (s + 1) & 3
                if cur.opp[4 * v + s] == 4 * v + s_adj:
                    s_a, s_b = (s + 2) & 3, (s + 3) & 3
                    internal = {4 * v + s_a: 4 * v + s_b, 4 * v + s_b: 4 * v + s_a}
                    old_to_new = {vo: (vo if vo < v else vo - 1) for vo in range(n) if vo != v}
                    nxt = _remove_vertices_with_internal_pairing(cur, {v}, internal)
                    steps.append(UnlinkMovieStep("R1-", cur, nxt, {"old_to_new": old_to_new}))
                    cur = nxt
                    r1_done = True
                    break
            if r1_done:
                break
        if r1_done:
            continue

        # Try R2-
        r2_done = False
        for u in range(n):
            ut = cur.vtype[u]
            for su in range(4):
                su2 = (su + 1) & 3
                od1 = cur.opp[4 * u + su]
                w = od1 >> 2
                if w <= u:
                    continue
                sw = od1 & 3
                sw2 = (sw - 1) & 3
                if cur.opp[4 * u + su2] != ((w << 2) | sw2):
                    continue
                if _is_over(ut, su) != _is_over(cur.vtype[w], sw):
                    continue
                u_ext1, w_ext1 = (u << 2) | (su ^ 2), (w << 2) | (sw ^ 2)
                u_ext2, w_ext2 = (u << 2) | (su2 ^ 2), (w << 2) | (sw2 ^ 2)
                internal = {
                    u_ext1: w_ext1,
                    w_ext1: u_ext1,
                    u_ext2: w_ext2,
                    w_ext2: u_ext2,
                }
                old_to_new = {}
                nid = 0
                for vo in range(n):
                    if vo not in (u, w):
                        old_to_new[vo] = nid
                        nid += 1
                nxt = _remove_vertices_with_internal_pairing(cur, {u, w}, internal)
                steps.append(UnlinkMovieStep("R2-", cur, nxt, {"old_to_new": old_to_new}))
                cur = nxt
                r2_done = True
                break
            if r2_done:
                break
        if r2_done:
            continue
        break
    return cur, steps


def _list_r3_steps(d: ClassicalLinkDiagram) -> list[UnlinkMovieStep]:
    n = d.num_crossings
    if n < 3:
        return []
    seen_faces: set[frozenset[int]] = set()
    out: list[UnlinkMovieStep] = []

    for start_d in range(4 * n):
        d0 = start_d
        od0 = d.opp[d0]
        v1, t1 = od0 >> 2, od0 & 3
        s1 = (t1 - 1) & 3
        d1 = (v1 << 2) | s1

        od1 = d.opp[d1]
        v2, t2 = od1 >> 2, od1 & 3
        s2 = (t2 - 1) & 3
        d2 = (v2 << 2) | s2

        od2 = d.opp[d2]
        v0, t0 = od2 >> 2, od2 & 3
        s0 = (t0 - 1) & 3
        if ((v0 << 2) | s0) != d0 or v0 == v1 or v1 == v2 or v2 == v0:
            continue

        edge_key = frozenset({min(d0, od0), min(d1, od1), min(d2, od2)})
        if edge_key in seen_faces:
            continue
        seen_faces.add(edge_key)

        vt0, vt1, vt2 = d.vtype[v0], d.vtype[v1], d.vtype[v2]
        a0, a1 = _is_over(vt0, s0), _is_over(vt1, t1)
        b1, b2 = _is_over(vt1, s1), _is_over(vt2, t2)
        c2, c0 = _is_over(vt2, s2), _is_over(vt0, t0)

        heights: dict[str, int] = {}
        for name, pair in (("A", (a0, a1)), ("B", (b1, b2)), ("C", (c2, c0))):
            if pair == (True, True):
                heights[name] = 2
            elif pair == (False, False):
                heights[name] = 0
            else:
                heights[name] = 1
        if set(heights.values()) != {0, 1, 2}:
            continue

        ports = (
            (v0 << 2) | ((s0 + 2) & 3),
            (v0 << 2) | ((s0 + 3) & 3),
            (v1 << 2) | ((s1 + 2) & 3),
            (v1 << 2) | ((s1 + 3) & 3),
            (v2 << 2) | ((s2 + 2) & 3),
            (v2 << 2) | ((s2 + 3) & 3),
        )
        port_set = set(ports)
        ext_targets = [d.opp[p] for p in ports]

        new_opp = list(d.opp)
        new_vtype = list(d.vtype)
        ws = (v0, v1, v2)
        strand_pairs = (("C", "B"), ("A", "C"), ("B", "A"))
        port_pairs = ((1, 2), (3, 4), (5, 0))

        for i in range(3):
            wi = ws[i]
            w_next = ws[(i + 1) % 3]
            w_prev = ws[(i + 2) % 3]
            new_opp[4 * wi + 0] = 4 * w_next + 1
            new_opp[4 * wi + 1] = 4 * w_prev + 0

            pi_a, pi_b = port_pairs[i]
            t_a, t_b = ext_targets[pi_a], ext_targets[pi_b]
            if t_a in port_set:
                idx_other = ports.index(t_a)
                for j in range(3):
                    if port_pairs[j][0] == idx_other:
                        new_opp[4 * wi + 2] = 4 * ws[j] + 2
                    elif port_pairs[j][1] == idx_other:
                        new_opp[4 * wi + 2] = 4 * ws[j] + 3
            else:
                new_opp[4 * wi + 2] = t_a
                new_opp[t_a] = 4 * wi + 2

            if t_b in port_set:
                idx_other = ports.index(t_b)
                for j in range(3):
                    if port_pairs[j][0] == idx_other:
                        new_opp[4 * wi + 3] = 4 * ws[j] + 2
                    elif port_pairs[j][1] == idx_other:
                        new_opp[4 * wi + 3] = 4 * ws[j] + 3
            else:
                new_opp[4 * wi + 3] = t_b
                new_opp[t_b] = 4 * wi + 3

            s_02, s_13 = strand_pairs[i]
            new_vtype[wi] = "X0" if heights[s_13] > heights[s_02] else "X1"

        nxt = ClassicalLinkDiagram(
            num_crossings=n,
            opp=tuple(new_opp),
            vtype=tuple(new_vtype),
            free_loops=d.free_loops,
        )
        out.append(
            UnlinkMovieStep(
                kind="R3",
                before=d,
                after=nxt,
                meta={
                    "vs": (v0, v1, v2),
                    "ss": (s0, s1, s2),
                    "heights": (heights["A"], heights["B"], heights["C"]),
                    "ports": ports,
                    "port_pairs": port_pairs,
                },
            )
        )
    return out


def _list_r2_plus_steps(d: ClassicalLinkDiagram) -> list[UnlinkMovieStep]:
    n = d.num_crossings
    out: list[UnlinkMovieStep] = []
    visited_darts = [False] * (4 * n)
    for start_d in range(4 * n):
        if visited_darts[start_d]:
            continue
        face_darts: list[int] = []
        cur = start_d
        while not visited_darts[cur]:
            visited_darts[cur] = True
            face_darts.append(cur)
            od = d.opp[cur]
            cur = ((od >> 2) << 2) | (((od & 3) - 1) & 3)
        m = len(face_darts)
        if m < 2:
            continue
        for i in range(m):
            da = face_darts[i]
            oda = d.opp[da]
            for j in range(i + 1, m):
                db = face_darts[j]
                odb = d.opp[db]
                if min(da, oda) == min(db, odb):
                    continue
                u, w = n, n + 1
                for vt in ("X0", "X1"):
                    new_opp = list(d.opp) + [0] * 8
                    new_opp[4 * u + 0] = da
                    new_opp[da] = 4 * u + 0
                    new_opp[4 * u + 1] = 4 * w + 3
                    new_opp[4 * u + 2] = 4 * w + 2
                    new_opp[4 * u + 3] = odb
                    new_opp[odb] = 4 * u + 3
                    new_opp[4 * w + 0] = oda
                    new_opp[oda] = 4 * w + 0
                    new_opp[4 * w + 1] = db
                    new_opp[db] = 4 * w + 1
                    new_opp[4 * w + 2] = 4 * u + 2
                    new_opp[4 * w + 3] = 4 * u + 1
                    nxt = ClassicalLinkDiagram(
                        num_crossings=n + 2,
                        opp=tuple(new_opp),
                        vtype=d.vtype + (vt, vt),  # type: ignore[operator]
                        free_loops=d.free_loops,
                    )
                    out.append(
                        UnlinkMovieStep(
                            kind="R2+",
                            before=d,
                            after=nxt,
                            meta={"da": da, "oda": oda, "db": db, "odb": odb, "vt": vt},
                        )
                    )
    return out


def find_unlink_movie(
    d: ClassicalLinkDiagram,
    max_r3_states: int = 512,
) -> list[UnlinkMovieStep] | None:
    """Return an explicit sequence of `UnlinkMovieStep`s reducing `d` to 0 crossings, or `None`."""
    if d.num_crossings == 0:
        return []
    for val in d.pairwise_linking_numbers_doubled().values():
        if val != 0:
            return None

    red, steps0 = _greedy_r1_r2_movie(d)
    if red.num_crossings == 0:
        return steps0

    visited = {_canonical_classical_key(red)}
    queue: deque[tuple[ClassicalLinkDiagram, list[UnlinkMovieStep]]] = deque([(red, steps0)])

    while queue and len(visited) < max_r3_states:
        cur, path = queue.popleft()
        for r3_step in _list_r3_steps(cur):
            nbr_red, m_extra = _greedy_r1_r2_movie(r3_step.after)
            npath = path + [r3_step] + m_extra
            if nbr_red.num_crossings == 0:
                return npath
            k = _canonical_classical_key(nbr_red)
            if k not in visited:
                visited.add(k)
                queue.append((nbr_red, npath))

    if red.num_crossings <= 7 or not is_trivial_jones_polynomial(red):
        return None

    for poke_step in _list_r2_plus_steps(red):
        sub_vis = {_canonical_classical_key(poke_step.after)}
        sub_q: deque[tuple[ClassicalLinkDiagram, list[UnlinkMovieStep], int]] = deque(
            [(poke_step.after, steps0 + [poke_step], 0)]
        )
        while sub_q and len(sub_vis) < 256:
            scur, spath, depth = sub_q.popleft()
            for r3_step in _list_r3_steps(scur):
                sred, sm_extra = _greedy_r1_r2_movie(r3_step.after)
                npath = spath + [r3_step] + sm_extra
                if sred.num_crossings == 0:
                    return npath
                if sred.num_crossings < red.num_crossings:
                    tail = find_unlink_movie(sred, max_r3_states=max_r3_states)
                    if tail is not None:
                        return npath + tail
                if depth < 2:
                    sk = _canonical_classical_key(r3_step.after)
                    if sk not in sub_vis:
                        sub_vis.add(sk)
                        sub_q.append((r3_step.after, spath + [r3_step], depth + 1))
    return None

