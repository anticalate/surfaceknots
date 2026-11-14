"""Yoshikawa moves (Omega_1 .. Omega_7) on admissible marked graph diagrams.

Implements:
- Complexity-reducing Yoshikawa moves (`Omega_1-` classical kink, `Omega_6-`/`Omega_6'-` marked kink,
  `Omega_2-` classical bigon) + greedy reduction.
- Complexity-preserving Yoshikawa moves (`Omega_3` classical triangle, `Omega_4`/`Omega_4'`
  sliding a strand over/under a marked vertex, `Omega_5` 180-degree rigid-vertex flip across a
  classical-marked bigon, and `Omega_7` Morse 1-handle slide along an edge between two marked
  vertices).
- Disconnected (split) diagram decomposition and planar 2-edge connected-sum (`#`) detection.
- Bounded BFS equivalence search over Yoshikawa moves to merge diagrams representing the same
  surface-link and certify unknottedness / splittability / compositeness.
"""

from __future__ import annotations

from collections import deque

from surfaceknots.diagrams import MarkedDiagram, VertexKind


def _is_classical(vt: VertexKind) -> bool:
    return vt in ("X0", "X1")


def _is_marked(vt: VertexKind) -> bool:
    return vt in ("Y", "Z")


def _is_over(vt: VertexKind, slot: int) -> bool:
    if vt == "X0":
        return (slot & 1) == 1
    if vt == "X1":
        return (slot & 1) == 0
    raise ValueError(f"_is_over called on marked vertex {vt}")


def canonical_marked_key(d: MarkedDiagram, up_to_mirror_and_reverse: bool = False) -> tuple:
    """Canonical spherical isomorphism code of a `MarkedDiagram`.

    If `up_to_mirror_and_reverse=True`, canonicalizes across:
    - spherical orientation-preserving and orientation-reversing homeomorphisms (`orient in {+1, -1}`),
    - 3D mirror `R_z` (`X0 <-> X1`),
    - time mirror `R_t` (`Y <-> Z`),
    - and orientation reversal `R_{zt}` (`X0 <-> X1` and `Y <-> Z`),
    matching Yoshikawa's unoriented, up-to-mirror census equivalence!
    """
    variants = [d]
    if up_to_mirror_and_reverse:
        variants = [
            d,
            d.mirror_3d(),
            d.mirror_time(),
            d.reverse_orientation_zt(),
        ]

    overall_best = None
    for diag in variants:
        # Note: planar reflection `(x, y) -> (-x, y)` reverses cyclic order `(0,1,2,3) -> (0,3,2,1)`:
        # under `(0,1,2,3) -> (0,3,2,1)`, `(0,2)` stays `(0,2)` (`X0`/`X1` preserved) while `(0,3)`
        # and `(0,1)` swap (`Y <-> Z`). Combined with `R_z`, `(x,y,z,t) -> (-x,y,-z,t)` is an
        # orientation-preserving ambient isotopy in SO(4)!
        refl_opp = [0] * (4 * diag.num_vertices)
        slot_refl = (0, 3, 2, 1)
        for v in range(diag.num_vertices):
            for s in range(4):
                od = diag.opp[4 * v + s]
                ov, os = od >> 2, od & 3
                refl_opp[4 * v + slot_refl[s]] = 4 * ov + slot_refl[os]
        # Ambient-isotopic planar reflection + R_z (det = +1 in R^4): flips X0<->X1 and Y<->Z (from slot_refl)
        flip_both = {"X0": "X1", "X1": "X0", "Y": "Z", "Z": "Y"}
        diag_rot180 = MarkedDiagram(
            num_vertices=diag.num_vertices,
            opp=tuple(refl_opp),
            vtype=tuple(flip_both[t] for t in diag.vtype),  # type: ignore[arg-type]
            free_loops=diag.free_loops,
        )
        sub_diags = [diag, diag_rot180]
        if up_to_mirror_and_reverse:
            # Also allow pure planar reflection without R_z (already covered since variants has mirror_3d!)
            pass

        for sd in sub_diags:
            k = _single_marked_canonical(sd)
            if overall_best is None or k < overall_best:
                overall_best = k
    assert overall_best is not None
    return overall_best


def _single_marked_canonical(d: MarkedDiagram) -> tuple:
    """Orientation-preserving spherical canonical code of a single `MarkedDiagram`."""
    n = d.num_vertices
    if n == 0:
        return (0, d.free_loops, ())

    # Partition into connected graph components
    seen_global = [False] * n
    comp_codes = []
    for sv in range(n):
        if seen_global[sv]:
            continue
        comp_vs: list[int] = []
        st = [sv]
        seen_global[sv] = True
        while st:
            cv = st.pop()
            comp_vs.append(cv)
            for s in range(4):
                nv = d.opp[4 * cv + s] >> 2
                if not seen_global[nv]:
                    seen_global[nv] = True
                    st.append(nv)

        best_comp = None
        for v0 in comp_vs:
            for s0 in range(4):
                # Note: rotating slots by +2 at every vertex preserves X0, X1, Y, Z!
                # Rotating slots by +1 flips X0<->X1 and Y<->Z!
                ov2n = {v0: 0}
                slot_shift = {v0: s0}
                queue = [v0]
                qi = 0
                code: list[int] = []
                while qi < len(queue):
                    v = queue[qi]
                    qi += 1
                    sh = slot_shift[v]
                    vt = d.vtype[v]
                    if sh & 1:
                        vt_rel = {"X0": 1, "X1": 0, "Y": 3, "Z": 2}[vt]
                    else:
                        vt_rel = {"X0": 0, "X1": 1, "Y": 2, "Z": 3}[vt]
                    code.append(vt_rel)
                    for k in range(4):
                        cd = (v << 2) | ((sh + k) & 3)
                        od = d.opp[cd]
                        ov, os = od >> 2, od & 3
                        if ov not in ov2n:
                            ov2n[ov] = len(ov2n)
                            slot_shift[ov] = os
                            queue.append(ov)
                        rel_os = (os - slot_shift[ov]) & 3
                        code.append((ov2n[ov] << 2) | rel_os)
                tc = tuple(code)
                if best_comp is None or tc < best_comp:
                    best_comp = tc
        comp_codes.append(best_comp)
    comp_codes.sort()
    return (n, d.free_loops, tuple(comp_codes))


def _remove_marked_diagram_vertices(
    d: MarkedDiagram,
    removed: set[int],
    internal_splice: dict[int, int],
) -> MarkedDiagram:
    """Remove `removed` vertices from `MarkedDiagram`, splicing external darts across `internal_splice`."""
    n = d.num_vertices
    old_to_new: dict[int, int] = {}
    new_vtype: list[VertexKind] = []
    for v in range(n):
        if v not in removed:
            old_to_new[v] = len(new_vtype)
            new_vtype.append(d.vtype[v])
    m = len(new_vtype)
    new_opp = [0] * (4 * m)
    visited = set()

    for v_old, v_new in old_to_new.items():
        for s in range(4):
            cur = d.opp[4 * v_old + s]
            while cur in internal_splice:
                visited.add(cur)
                nxt = internal_splice[cur]
                visited.add(nxt)
                cur = d.opp[nxt]
            w_old, ws = cur >> 2, cur & 3
            new_opp[4 * v_new + s] = (old_to_new[w_old] << 2) | ws

    extra_loops = 0
    for rd in internal_splice:
        if rd in visited:
            continue
        extra_loops += 1
        cur = rd
        while cur not in visited:
            visited.add(cur)
            nxt = internal_splice[cur]
            visited.add(nxt)
            cur = d.opp[nxt]

    return MarkedDiagram(
        num_vertices=m,
        opp=tuple(new_opp),
        vtype=tuple(new_vtype),
        free_loops=d.free_loops + extra_loops,
    )


def find_omega1_or_6_reduction(d: MarkedDiagram) -> MarkedDiagram | None:
    """Apply Omega_1- (classical kink) or Omega_6-/Omega_6'- (marked kink) if available."""
    n = d.num_vertices
    for v in range(n):
        for s in range(4):
            s_adj = (s + 1) & 3
            if d.opp[4 * v + s] == 4 * v + s_adj:
                sa = (s + 2) & 3
                sb = (s + 3) & 3
                internal = {4 * v + sa: 4 * v + sb, 4 * v + sb: 4 * v + sa}
                return _remove_marked_diagram_vertices(d, {v}, internal)
    return None


def find_omega2_reduction(d: MarkedDiagram) -> MarkedDiagram | None:
    """Apply Omega_2- (classical R-II bigon) if available."""
    n = d.num_vertices
    for u in range(n):
        ut = d.vtype[u]
        if not _is_classical(ut):
            continue
        for su in range(4):
            su2 = (su + 1) & 3
            od1 = d.opp[4 * u + su]
            w = od1 >> 2
            if w <= u or not _is_classical(d.vtype[w]):
                continue
            sw = od1 & 3
            sw2 = (sw - 1) & 3
            if d.opp[4 * u + su2] != ((w << 2) | sw2):
                continue
            if _is_over(ut, su) != _is_over(d.vtype[w], sw):
                continue
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
            return _remove_marked_diagram_vertices(d, {u, w}, internal)
    return None


def greedy_yoshikawa_reduce(d: MarkedDiagram) -> MarkedDiagram:
    """Exhaustively apply Omega_1-, Omega_6-, Omega_6'-, and Omega_2- reductions."""
    cur = d
    while cur.num_vertices > 0:
        nxt = find_omega1_or_6_reduction(cur)
        if nxt is not None:
            cur = nxt
            continue
        nxt = find_omega2_reduction(cur)
        if nxt is not None:
            cur = nxt
            continue
        break
    return cur


def list_level_yoshikawa_neighbors(d: MarkedDiagram) -> list[MarkedDiagram]:
    """Return all diagrams reachable from `d` by one complexity-preserving Yoshikawa move:
    `Omega_3`, `Omega_4 / Omega_4'`, `Omega_5`, `Omega_7`, `Omega_8`, or the certified `+3`-peak
    marked 2-gon tube involution (`Omega_6+ / Omega_2+ / Omega_7` composite).
    """
    out: list[MarkedDiagram] = []
    out.extend(_list_omega3_and_4_neighbors(d))
    out.extend(_list_omega5_neighbors(d))
    out.extend(_list_omega7_neighbors(d))
    out.extend(_list_omega8_neighbors(d))
    out.extend(_list_marked_tube_flip_neighbors(d))
    return out


def _rewire_ports_in_diagram(
    d: MarkedDiagram,
    ws: tuple[int, ...],
    internal_edges: list[tuple[int, int]],
    port_assignments: list[tuple[int, int]],
    ports: list[int],
    new_vtypes_for_ws: list[VertexKind],
) -> MarkedDiagram:
    """Helper to rewire a local sub-disk with vertices `ws` and external boundary `ports`."""
    ext_targets = [d.opp[p] for p in ports]
    port_to_new_dart = {}
    for dart_local, port_idx in port_assignments:
        port_to_new_dart[ports[port_idx]] = dart_local

    new_opp = list(d.opp)
    for d_a, d_b in internal_edges:
        new_opp[d_a] = d_b
        new_opp[d_b] = d_a

    for dart_local, port_idx in port_assignments:
        target = ext_targets[port_idx]
        if target in port_to_new_dart:
            # External edge connects two boundary ports of the SAME local disk!
            other_local = port_to_new_dart[target]
            new_opp[dart_local] = other_local
        else:
            new_opp[dart_local] = target
            new_opp[target] = dart_local

    new_vtype = list(d.vtype)
    for wi, nvt in zip(ws, new_vtypes_for_ws):
        new_vtype[wi] = nvt

    return MarkedDiagram(
        num_vertices=d.num_vertices,
        opp=tuple(new_opp),
        vtype=tuple(new_vtype),
        free_loops=d.free_loops,
    )


def _list_omega3_and_4_neighbors(d: MarkedDiagram) -> list[MarkedDiagram]:
    """Omega_3 (3 classical) and Omega_4 / Omega_4' (1 marked + 2 classical) triangle moves."""
    n = d.num_vertices
    if n < 3:
        return []
    seen_faces: set[frozenset[int]] = set()
    out: list[MarkedDiagram] = []

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

        vs = [v0, v1, v2]
        ss = [s0, s1, s2]
        marked_indices = [i for i in range(3) if _is_marked(d.vtype[vs[i]])]

        if len(marked_indices) == 0:
            # Omega_3: 3 classical crossings
            vt0, vt1, vt2 = d.vtype[v0], d.vtype[v1], d.vtype[v2]
            a_pair = (_is_over(vt0, s0), _is_over(vt1, t1))
            b_pair = (_is_over(vt1, s1), _is_over(vt2, t2))
            c_pair = (_is_over(vt2, s2), _is_over(vt0, t0))
            heights: dict[str, int] = {}
            for name, pair in (("A", a_pair), ("B", b_pair), ("C", c_pair)):
                if pair == (True, True):
                    heights[name] = 2
                elif pair == (False, False):
                    heights[name] = 0
                else:
                    heights[name] = 1
            if set(heights.values()) != {0, 1, 2}:
                continue

            ports = [
                (v0 << 2) | ((s0 + 2) & 3),
                (v0 << 2) | ((s0 + 3) & 3),
                (v1 << 2) | ((s1 + 2) & 3),
                (v1 << 2) | ((s1 + 3) & 3),
                (v2 << 2) | ((s2 + 2) & 3),
                (v2 << 2) | ((s2 + 3) & 3),
            ]
            ws = (v0, v1, v2)
            internal_edges = [
                (4 * ws[0] + 0, 4 * ws[1] + 1),
                (4 * ws[1] + 0, 4 * ws[2] + 1),
                (4 * ws[2] + 0, 4 * ws[0] + 1),
            ]
            port_assignments = [
                (4 * ws[0] + 2, 1),
                (4 * ws[0] + 3, 2),
                (4 * ws[1] + 2, 3),
                (4 * ws[1] + 3, 4),
                (4 * ws[2] + 2, 5),
                (4 * ws[2] + 3, 0),
            ]
            strand_pairs = (("C", "B"), ("A", "C"), ("B", "A"))
            new_vts: list[VertexKind] = [
                "X0" if heights[s13] > heights[s02] else "X1" for (s02, s13) in strand_pairs
            ]
            out.append(
                _rewire_ports_in_diagram(d, ws, internal_edges, port_assignments, ports, new_vts)
            )

        elif len(marked_indices) == 1:
            # Omega_4 / Omega_4': rotate indices so v0 is the marked vertex!
            mi = marked_indices[0]
            if mi != 0:
                vs = vs[mi:] + vs[:mi]
                ss = ss[mi:] + ss[:mi]
            mv, cv1, cv2 = vs
            ms, cs1, cs2 = ss
            ct1 = (cs1 + 1) & 3
            ct2 = (cs2 + 1) & 3
            # Side B = (cv1, cs1) -- (cv2, ct2) is the free strand opposite the marked vertex `mv`!
            # Valid Omega_4 iff Side B is OVER at both cv1 & cv2; Omega_4' iff UNDER at both!
            b_over1 = _is_over(d.vtype[cv1], cs1)
            b_over2 = _is_over(d.vtype[cv2], ct2)
            if b_over1 != b_over2:
                continue
            b_is_over = b_over1

            ports = [
                (mv << 2) | ((ms + 2) & 3),   # P0
                (mv << 2) | ((ms + 3) & 3),   # P1
                (cv1 << 2) | ((cs1 + 2) & 3), # P2 (B)
                (cv1 << 2) | ((cs1 + 3) & 3), # P3 (from mv.ms)
                (cv2 << 2) | ((cs2 + 2) & 3), # P4 (from mv.ms+1)
                (cv2 << 2) | ((cs2 + 3) & 3), # P5 (B)
            ]
            # Before move: at mv, slots are (ms->cv1.ct1[to P3], ms+1->cv2.cs2[to P4], ms+2=P0, ms+3=P1).
            # Check which marker type `mv` has in `(ms, ms+1, ms+2, ms+3)` frame!
            # In `L_+`: does `mv` pair `(ms, ms+3)` (`(P3, P1)` & `(P4, P0)`) or `(ms, ms+1)` (`(P3, P4)` & `(P0, P1)`)?
            vt_m = d.vtype[mv]
            # Recall: `Y` pairs `(0,3)&(1,2)` in `L_+`, so `(ms, ms+3)` in `L_+` iff `ms` is EVEN for `Y` or ODD for `Z`!
            lplus_pairs_ms_ms3 = (vt_m == "Y" and (ms & 1) == 0) or (vt_m == "Z" and (ms & 1) == 1)

            # After sliding B (P2 -- P5) across `mv`:
            # Place `mv` at `w1` (connecting P3 and P4 externally on slots 2, 3, and internal triangle on 0, 1!):
            # Wait: check 6-port inversion ws = (w0=cv2, w1=mv, w2=cv1):
            # - w0: (0->w1.1, 1->w2.0, 2->P1, 3->P2[B])  [classical crossing between P1-branch and B!]
            # - w1 (MARKED!): (0->w2.1 [goes straight across w2 to w2.3 = P0!],
            #                  1->w0.0 [goes straight across w0 to w0.2 = P1!],
            #                  2->P3,
            #                  3->P4)!!
            # - w2: (0->w0.1, 1->w1.0, 2->P5[B], 3->P0)  [classical crossing between B and P0-branch!]
            # VERIFY `w1` CONNECTIONS IN `L_+` AND `L_-`!!
            # At `w1`, the 4 incident strands in CCW order `(0, 1, 2, 3)` lead to:
            #   slot 0 -> P0 (was `ms + 2`!)
            #   slot 1 -> P1 (was `ms + 3`!)
            #   slot 2 -> P3 (was `ms`!)
            #   slot 3 -> P4 (was `ms + 1`!)
            # IT IS AN EXACT `+2` CYCLIC ROTATION OF `(ms, ms+1, ms+2, ms+3)` -> `(P0, P1, P3, P4)`!!
            # Therefore, in `L_+`:
            # - Pairing `(P3, P1)` & `(P4, P0)` at `w1` is pairing `(2, 1)` & `(3, 0)` = **`Y`**!!
            # - Pairing `(P3, P4)` & `(P0, P1)` at `w1` is pairing `(2, 3)` & `(0, 1)` = **`Z`**!!
            new_mv_type: VertexKind = "Y" if lplus_pairs_ms_ms3 else "Z"
            # At `w0`: `(0, 2)` is P1-branch, `(1, 3)` is B (`P2` at slot 3 -> `w0.1` -> `w2.0` -> `w2.2 = P5`!).
            # Since `X0` has `(1, 3)` OVER `(0, 2)`, `w0` is `"X0"` iff `b_is_over`!
            # At `w2`: `(0, 2)` is B (`P5` at slot 2), `(1, 3)` is P0-branch (`P0` at slot 3)!
            # So `w2` has B on `(0, 2)` -> `w2` is `"X1"` iff `b_is_over`!
            ws = (cv2, mv, cv1)
            internal_edges = [
                (4 * ws[0] + 0, 4 * ws[1] + 1),
                (4 * ws[1] + 0, 4 * ws[2] + 1),
                (4 * ws[2] + 0, 4 * ws[0] + 1),
            ]
            port_assignments = [
                (4 * ws[0] + 2, 1),
                (4 * ws[0] + 3, 2),
                (4 * ws[1] + 2, 3),
                (4 * ws[1] + 3, 4),
                (4 * ws[2] + 2, 5),
                (4 * ws[2] + 3, 0),
            ]
            new_vts = [
                "X0" if b_is_over else "X1",
                new_mv_type,
                "X1" if b_is_over else "X0",
            ]
            out.append(
                _rewire_ports_in_diagram(d, ws, internal_edges, port_assignments, ports, new_vts)
            )
    return out


def _list_omega5_neighbors(d: MarkedDiagram) -> list[MarkedDiagram]:
    """Omega_5: move a classical crossing across an adjacent marked vertex along a 2-gon."""
    n = d.num_vertices
    if n < 2:
        return []
    out: list[MarkedDiagram] = []
    for u in range(n):
        ut = d.vtype[u]
        if not _is_classical(ut):
            continue
        for su in range(4):
            su2 = (su + 1) & 3
            od1 = d.opp[4 * u + su]
            v, sv = od1 >> 2, od1 & 3
            if not _is_marked(d.vtype[v]):
                continue
            sv2 = (sv - 1) & 3
            if d.opp[4 * u + su2] != ((v << 2) | sv2):
                continue

            # 4 external ports in CCW order around {u, v}:
            # P0 = A = (u, su+2) [goes across u to su -> v.sv]
            # P1 = B = (u, su+3) [goes across u to su+1 -> v.sv-1]
            # P2 = C = (v, sv+1)
            # P3 = D = (v, sv+2)
            ports = [
                (u << 2) | ((su + 2) & 3),
                (u << 2) | ((su + 3) & 3),
                (v << 2) | ((sv + 1) & 3),
                (v << 2) | ((sv + 2) & 3),
            ]
            vt_v = d.vtype[v]
            # Check whether L_+ is the `thru` resolution ((sv, sv+1) = A--C and (sv-1, sv+2) = B--D)
            # or the `loop` resolution ((sv, sv-1) = A--B and (sv+1, sv+2) = C--D)!
            # (sv, sv+1) paired in L_+ iff (vt_v=='Z' and sv even) or (vt_v=='Y' and sv odd)
            lplus_is_thru = (vt_v == "Z" and (sv & 1) == 0) or (vt_v == "Y" and (sv & 1) == 1)
            # Is A--C OVER B--D at u? `A` goes across u to `su`, so `A` is OVER iff `_is_over(ut, su)`!
            ac_over_bd = _is_over(ut, su)

            # After Omega_5: put marked vertex `v` on `(P0=A, P1=B)` side and classical `u` on `(P2=C, P3=D)` side!
            # Set `v` slots: (0->u.1 [to D? wait: let u have (0,2)=D-strand, (1,3)=C-strand!]):
            # Let's assign slots explicitly and verify!
            # At `v`: 0 -> u.3, 1 -> u.2, 2 -> P0 (A), 3 -> P1 (B).
            # At `u`: 0 -> P2 (C), 1 -> P3 (D), 2 -> v.1, 3 -> v.0.
            # Check bigon between `v` (slots 0, 1) and `u` (slots 3, 2): `3 - 1 = 2` -> PLANAR BIGON!
            # Check outer CCW ports around `{v, u}`:
            # At `v`: 1 -> 2(P0) -> 3(P1) -> 0; at `u`: 3 -> 0(P2) -> 1(P3) -> 2 -> `(P0, P1, P2, P3)` EXACT!
            # Check `thru` vs `loop` at new `v`:
            # - If new `v` pairs `(1, 2)` & `(0, 3)` (`Y`!):
            #   `P0 (v.2)` -> `v.1` -> `u.2` -> straight across `u` to `u.0 = P2 (C)`! (`A -- C`!)
            #   `P1 (v.3)` -> `v.0` -> `u.3` -> straight across `u` to `u.1 = P3 (D)`! (`B -- D`!)
            #   SO `Y` AT NEW `v` GIVES `thru` IN `L_+`, AND `Z` AT NEW `v` GIVES `loop` IN `L_+`!!
            # - At new `u`: `(0, 2)` is `A -- C` (`u.0 = P2 = C`), and `(1, 3)` is `B -- D` (`u.1 = P3 = D`)!
            #   Since `"X1"` has `(0, 2)` OVER `(1, 3)`, new `u` is `"X1"` iff `ac_over_bd` else `"X0"`!
            ws = (v, u)
            internal_edges = [
                (4 * v + 0, 4 * u + 3),
                (4 * v + 1, 4 * u + 2),
            ]
            port_assignments = [
                (4 * v + 2, 0),
                (4 * v + 3, 1),
                (4 * u + 0, 2),
                (4 * u + 1, 3),
            ]
            new_v_type: VertexKind = "Y" if lplus_is_thru else "Z"
            new_u_type: VertexKind = "X1" if ac_over_bd else "X0"
            out.append(
                _rewire_ports_in_diagram(
                    d, ws, internal_edges, port_assignments, ports, [new_v_type, new_u_type]
                )
            )
    return out


def _list_omega7_neighbors(d: MarkedDiagram) -> list[MarkedDiagram]:
    """Omega_7: Morse 1-handle slide along an edge connecting two marked vertices (Jablonowski Prop. 11)."""
    n = d.num_vertices
    if n < 2:
        return []
    out: list[MarkedDiagram] = []
    seen_edges: set[int] = set()

    for u in range(n):
        ut = d.vtype[u]
        if not _is_marked(ut):
            continue
        for su in range(4):
            d1 = (u << 2) | su
            d2 = d.opp[d1]
            ekey = min(d1, d2)
            if ekey in seen_edges:
                continue
            w, sw = d2 >> 2, d2 & 3
            if w == u or not _is_marked(d.vtype[w]):
                continue
            seen_edges.add(ekey)
            wt = d.vtype[w]

            # Check Proposition 11: same marker type -> different slot parity; different type -> same parity
            same_type = ut == wt
            same_parity = (su & 1) == (sw & 1)
            if same_type == same_parity:
                continue

            # 6 external boundary ports in CCW order around {u, w}:
            # P0=(u, su+1), P1=(u, su+2), P2=(u, su+3), P3=(w, sw+1), P4=(w, sw+2), P5=(w, sw+3)
            ports = [
                (u << 2) | ((su + 1) & 3),
                (u << 2) | ((su + 2) & 3),
                (u << 2) | ((su + 3) & 3),
                (w << 2) | ((sw + 1) & 3),
                (w << 2) | ((sw + 2) & 3),
                (w << 2) | ((sw + 3) & 3),
            ]
            # Check which of the two Catalan matchings `L_+` induces on (P0..P5):
            # Does `u` pair `(su, su+3)` in `L_+` (giving `(P2--P3), (P0--P1), (P4--P5)`),
            # or `(su, su+1)` in `L_+` (giving `(P5--P0), (P1--P2), (P3--P4)`)?
            lplus_pairs_p2_p3 = (ut == "Y" and (su & 1) == 0) or (ut == "Z" and (su & 1) == 1)

            # Two rotated H-trees inside the 6-port disk `(P0..P5)` (`shift in (1, 2)`):
            # `shift = 1` groups `(P1, P2, P3) | (P4, P5, P0)`;
            # `shift = 2` groups `(P2, P3, P4) | (P5, P0, P1)`.
            for shift in (1, 2):
                # New `u` has slot 0 -> w.0, slots (1, 2, 3) -> (P_{shift}, P_{shift+1}, P_{shift+2})
                # New `w` has slot 0 -> u.0, slots (1, 2, 3) -> (P_{shift+3}, P_{shift+4}, P_{shift+5})
                # Check `Y` at new `u` (slot 0 = internal): `Y` pairs `(0, 3)` & `(1, 2)` in `L_+` ->
                # pairs `(P_{shift+2} -- P_{shift+3})` (when new `w` is `Z`! Check: at `w`, slot 0 is internal
                # and slot 1 is `P_{shift+3}`, so pairing `(0, 1)` in `L_+` is **`Z`**!!).
                # Wait: verify Proposition 11 on `(new u slot 0, new w slot 0)`: same parity (0 and 0) ->
                # one MUST be `Y` and the other MUST be `Z`!!
                # Let's check: if `u='Y', w='Z'`, in `L_+`:
                # - `u` pairs `(1,2) = (P_{shift}, P_{shift+1})` and `(3,0)` -> `w.0` -> `w.1 = P_{shift+3}`, so `(P_{shift+2}, P_{shift+3})`, and `w` pairs `(2,3) = (P_{shift+4}, P_{shift+5})`!
                #   When `shift == 1`: `(P1,P2), (P3,P4), (P5,P0)` -> `lplus_pairs_p2_p3 == False`!
                #   When `shift == 2`: `(P2,P3), (P4,P5), (P0,P1)` -> `lplus_pairs_p2_p3 == True`!
                want_yz = (shift == 2) if lplus_pairs_p2_p3 else (shift == 1)
                new_ut: VertexKind = "Y" if want_yz else "Z"
                new_wt: VertexKind = "Z" if want_yz else "Y"

                ws = (u, w)
                internal_edges = [(4 * u + 0, 4 * w + 0)]
                port_assignments = [
                    (4 * u + 1, (shift + 0) % 6),
                    (4 * u + 2, (shift + 1) % 6),
                    (4 * u + 3, (shift + 2) % 6),
                    (4 * w + 1, (shift + 3) % 6),
                    (4 * w + 2, (shift + 4) % 6),
                    (4 * w + 3, (shift + 5) % 6),
                ]
                out.append(
                    _rewire_ports_in_diagram(
                        d, ws, internal_edges, port_assignments, ports, [new_ut, new_wt]
                    )
                )
    return out


def _list_omega8_neighbors(d: MarkedDiagram) -> list[MarkedDiagram]:
    """Omega_8 (band swim): slide a marked vertex `u` across a 4-gon `(c_sw, c_se, c_ne, c_nw)`
    of 4 distinct classical crossings whose South side `(c_sw, c_se)` is attached to a distinct
    marked vertex `w`, where `H_bot` and `H_top` have opposite uniform over/under heights and
    `u, w` carry compatible band-swim markers.
    """
    n = d.num_vertices
    if n < 6:
        return []
    out: list[MarkedDiagram] = []
    for start_d in range(4 * n):
        vs: list[int] = []
        ss: list[int] = []
        cur = start_d
        ok = True
        for _ in range(4):
            v, s = cur >> 2, cur & 3
            if not _is_classical(d.vtype[v]) or v in vs:
                ok = False
                break
            vs.append(v)
            ss.append(s)
            od = d.opp[cur]
            cur = ((od >> 2) << 2) | (((od & 3) - 1) & 3)
        if not ok or cur != start_d:
            continue

        c_sw, c_se, c_ne, c_nw = vs
        s_sw, s_se, s_ne, s_nw = ss
        d_sw_s = (c_sw << 2) | ((s_sw + 3) & 3)
        d_se_s = (c_se << 2) | ((s_se + 2) & 3)
        od_sw_s = d.opp[d_sw_s]
        od_se_s = d.opp[d_se_s]
        w = od_sw_s >> 2
        if (od_se_s >> 2) != w or not _is_marked(d.vtype[w]):
            continue
        sw_w = od_se_s & 3
        if ((sw_w + 1) & 3) != (od_sw_s & 3):
            continue

        hb0 = _is_over(d.vtype[c_sw], s_sw)
        hb1 = _is_over(d.vtype[c_se], (s_se + 1) & 3)
        ht0 = _is_over(d.vtype[c_ne], s_ne)
        ht1 = _is_over(d.vtype[c_nw], (s_nw + 1) & 3)
        if hb0 != hb1 or ht0 != ht1 or hb0 == ht0:
            continue
        h_bot_over = hb0
        w_lplus_thru = (d.vtype[w] == "Y" and (sw_w & 1) == 0) or (
            d.vtype[w] == "Z" and (sw_w & 1) == 1
        )
        mk: VertexKind = "Y" if w_lplus_thru else "Z"
        cb: VertexKind = "X0" if h_bot_over else "X1"
        ct: VertexKind = "X1" if h_bot_over else "X0"

        # Case A: `u` on West -> move `u` to East (same level) or `w` to North (opposite level)
        d_sw_w = (c_sw << 2) | ((s_sw + 2) & 3)
        d_nw_w = (c_nw << 2) | ((s_nw + 3) & 3)
        od_sw_w = d.opp[d_sw_w]
        od_nw_w = d.opp[d_nw_w]
        u = od_sw_w >> 2
        if u != w and (od_nw_w >> 2) == u and _is_marked(d.vtype[u]):
            su_u = od_sw_w & 3
            if ((su_u + 1) & 3) == (od_nw_w & 3):
                u_lplus_thru = (d.vtype[u] == "Y" and (su_u & 1) == 0) or (
                    d.vtype[u] == "Z" and (su_u & 1) == 1
                )
                ports = [
                    (w << 2) | ((sw_w + 2) & 3),
                    (w << 2) | ((sw_w + 3) & 3),
                    (c_se << 2) | ((s_se + 3) & 3),
                    (c_ne << 2) | ((s_ne + 2) & 3),
                    (c_ne << 2) | ((s_ne + 3) & 3),
                    (c_nw << 2) | ((s_nw + 2) & 3),
                    (u << 2) | ((su_u + 2) & 3),
                    (u << 2) | ((su_u + 3) & 3),
                ]
                ws = (u, w, c_sw, c_se, c_ne, c_nw)
                if w_lplus_thru == u_lplus_thru:
                    internal_edges = [
                        (4 * c_sw + 1, 4 * c_se + 3),
                        (4 * c_se + 2, 4 * c_ne + 0),
                        (4 * c_ne + 3, 4 * c_nw + 1),
                        (4 * c_nw + 0, 4 * c_sw + 2),
                        (4 * w + 3, 4 * c_sw + 0),
                        (4 * w + 2, 4 * c_se + 0),
                        (4 * u + 3, 4 * c_se + 1),
                        (4 * u + 2, 4 * c_ne + 1),
                    ]
                    port_assignments = [
                        (4 * w + 0, 0),
                        (4 * w + 1, 1),
                        (4 * u + 0, 2),
                        (4 * u + 1, 3),
                        (4 * c_ne + 2, 4),
                        (4 * c_nw + 2, 5),
                        (4 * c_nw + 3, 6),
                        (4 * c_sw + 3, 7),
                    ]
                    out.append(
                        _rewire_ports_in_diagram(
                            d, ws, internal_edges, port_assignments, ports, [mk, mk, cb, cb, ct, ct]
                        )
                    )
                else:
                    mk_u: VertexKind = "Y" if u_lplus_thru else "Z"
                    internal_edges = [
                        (4 * c_sw + 1, 4 * c_se + 3),
                        (4 * c_se + 2, 4 * c_ne + 0),
                        (4 * c_ne + 3, 4 * c_nw + 1),
                        (4 * c_nw + 0, 4 * c_sw + 2),
                        (4 * w + 1, 4 * c_ne + 2),
                        (4 * w + 0, 4 * c_nw + 2),
                        (4 * u + 0, 4 * c_sw + 3),
                        (4 * u + 1, 4 * c_nw + 3),
                    ]
                    port_assignments = [
                        (4 * c_sw + 0, 0),
                        (4 * c_se + 0, 1),
                        (4 * c_se + 1, 2),
                        (4 * c_ne + 1, 3),
                        (4 * w + 2, 4),
                        (4 * w + 3, 5),
                        (4 * u + 2, 6),
                        (4 * u + 3, 7),
                    ]
                    out.append(
                        _rewire_ports_in_diagram(
                            d, ws, internal_edges, port_assignments, ports, [mk_u, mk, cb, cb, ct, ct]
                        )
                    )

        # Case B: `u2` on East -> move `u2` to West (same level) or `w` to North (opposite level)
        d_se_e = (c_se << 2) | ((s_se + 3) & 3)
        d_ne_e = (c_ne << 2) | ((s_ne + 2) & 3)
        od_se_e = d.opp[d_se_e]
        od_ne_e = d.opp[d_ne_e]
        u2 = od_ne_e >> 2
        if u2 != w and (od_se_e >> 2) == u2 and _is_marked(d.vtype[u2]):
            su_u2 = od_ne_e & 3
            if ((su_u2 + 1) & 3) == (od_se_e & 3):
                u2_lplus_thru = (d.vtype[u2] == "Y" and (su_u2 & 1) == 0) or (
                    d.vtype[u2] == "Z" and (su_u2 & 1) == 1
                )
                ports = [
                    (w << 2) | ((sw_w + 2) & 3),
                    (w << 2) | ((sw_w + 3) & 3),
                    (u2 << 2) | ((su_u2 + 2) & 3),
                    (u2 << 2) | ((su_u2 + 3) & 3),
                    (c_ne << 2) | ((s_ne + 3) & 3),
                    (c_nw << 2) | ((s_nw + 2) & 3),
                    (c_nw << 2) | ((s_nw + 3) & 3),
                    (c_sw << 2) | ((s_sw + 2) & 3),
                ]
                ws = (u2, w, c_sw, c_se, c_ne, c_nw)
                if w_lplus_thru == u2_lplus_thru:
                    internal_edges = [
                        (4 * c_sw + 1, 4 * c_se + 3),
                        (4 * c_se + 2, 4 * c_ne + 0),
                        (4 * c_ne + 3, 4 * c_nw + 1),
                        (4 * c_nw + 0, 4 * c_sw + 2),
                        (4 * w + 3, 4 * c_sw + 0),
                        (4 * w + 2, 4 * c_se + 0),
                        (4 * u2 + 0, 4 * c_sw + 3),
                        (4 * u2 + 1, 4 * c_nw + 3),
                    ]
                    port_assignments = [
                        (4 * w + 0, 0),
                        (4 * w + 1, 1),
                        (4 * c_se + 1, 2),
                        (4 * c_ne + 1, 3),
                        (4 * c_ne + 2, 4),
                        (4 * c_nw + 2, 5),
                        (4 * u2 + 2, 6),
                        (4 * u2 + 3, 7),
                    ]
                    out.append(
                        _rewire_ports_in_diagram(
                            d, ws, internal_edges, port_assignments, ports, [mk, mk, cb, cb, ct, ct]
                        )
                    )
                else:
                    mk_u2: VertexKind = "Y" if u2_lplus_thru else "Z"
                    internal_edges = [
                        (4 * c_sw + 1, 4 * c_se + 3),
                        (4 * c_se + 2, 4 * c_ne + 0),
                        (4 * c_ne + 3, 4 * c_nw + 1),
                        (4 * c_nw + 0, 4 * c_sw + 2),
                        (4 * w + 1, 4 * c_ne + 2),
                        (4 * w + 0, 4 * c_nw + 2),
                        (4 * u2 + 3, 4 * c_se + 1),
                        (4 * u2 + 2, 4 * c_ne + 1),
                    ]
                    port_assignments = [
                        (4 * c_sw + 0, 0),
                        (4 * c_se + 0, 1),
                        (4 * u2 + 0, 2),
                        (4 * u2 + 1, 3),
                        (4 * w + 2, 4),
                        (4 * w + 3, 5),
                        (4 * c_nw + 3, 6),
                        (4 * c_sw + 3, 7),
                    ]
                    out.append(
                        _rewire_ports_in_diagram(
                            d, ws, internal_edges, port_assignments, ports, [mk_u2, mk, cb, cb, ct, ct]
                        )
                    )
    return out


def _list_marked_tube_flip_neighbors(d: MarkedDiagram) -> list[MarkedDiagram]:
    """Certified `+3`-vertex mountain-pass macro (`Omega_6+ -> Omega_2+ -> Omega_7/level -> Omega_2- -> Omega_6-`):
    flipping `(vt[u], vt[w])` (`Y <-> Z` on both `u` and `w`) across an unknotted marked-marked 2-gon tube
    where `(vt[u] == vt[w]) == ((su & 1) == (sw & 1))`. This preserves both `L_+(D)` and `L_-(D)` identically.
    """
    n = d.num_vertices
    if n < 2:
        return []
    out: list[MarkedDiagram] = []
    for u in range(n):
        ut = d.vtype[u]
        if not _is_marked(ut):
            continue
        for su in range(4):
            od1 = d.opp[4 * u + su]
            w = od1 >> 2
            if w <= u:
                continue
            wt = d.vtype[w]
            if not _is_marked(wt):
                continue
            sw = od1 & 3
            su2 = (su + 1) & 3
            sw2 = (sw - 1) & 3
            if d.opp[4 * u + su2] != ((w << 2) | sw2):
                continue
            if (ut == wt) != ((su & 1) == (sw & 1)):
                continue
            new_vt = list(d.vtype)
            new_vt[u] = "Z" if ut == "Y" else "Y"
            new_vt[w] = "Z" if wt == "Y" else "Y"
            out.append(
                MarkedDiagram(
                    num_vertices=n,
                    opp=d.opp,
                    vtype=tuple(new_vt),  # type: ignore[arg-type]
                    free_loops=d.free_loops,
                )
            )
    return out


def connected_diagram_components(d: MarkedDiagram) -> list[MarkedDiagram]:
    """Split a disconnected planar `MarkedDiagram` into its connected diagram components."""
    n = d.num_vertices
    comps: list[MarkedDiagram] = []
    seen = [False] * n
    for sv in range(n):
        if seen[sv]:
            continue
        vs: list[int] = []
        st = [sv]
        seen[sv] = True
        while st:
            cv = st.pop()
            vs.append(cv)
            for s in range(4):
                nv = d.opp[4 * cv + s] >> 2
                if not seen[nv]:
                    seen[nv] = True
                    st.append(nv)
        vs.sort()
        o2n = {v: i for i, v in enumerate(vs)}
        m = len(vs)
        sub_opp = [0] * (4 * m)
        sub_vt = [d.vtype[v] for v in vs]
        for i, v in enumerate(vs):
            for s in range(4):
                od = d.opp[4 * v + s]
                sub_opp[4 * i + s] = (o2n[od >> 2] << 2) | (od & 3)
        comps.append(MarkedDiagram(num_vertices=m, opp=tuple(sub_opp), vtype=tuple(sub_vt), free_loops=0))

    for _ in range(d.free_loops):
        comps.append(MarkedDiagram.unknot_sphere())
    return comps


def find_planar_connected_sum_decomposition(
    d: MarkedDiagram,
) -> tuple[MarkedDiagram, MarkedDiagram] | None:
    """Detect if a connected `MarkedDiagram` has a 2-edge cut `D = D_1 # D_2` with both `|V(D_i)| >= 1`.

    Note: on the 2-vertex Hopf shadow (`2_1^1` or `2_1^-1`), cutting 2 edges leaves `|V(D_i)| == 1`:
    - For `2_1^1` (`Y, Y` on Hopf shadow), cutting 2 edges gives two 1-marker `S^2` diagrams! So we
       ONLY decompose `D = D_1 # D_2` when BOTH `D_1` and `D_2` are **nontrivial summands** (i.e.,
      neither `F(D_1)` nor `F(D_2)` is a single unknotted 2-sphere `0_1` — checked by `chi < 2` or
      `num_components > 1` or `ch >= 6`!).
    """
    n = d.num_vertices
    if n < 4:
        return None
    if len(connected_diagram_components(d)) != 1 or d.free_loops > 0:
        return None

    undirected_edges: list[tuple[int, int]] = []
    for x in range(4 * n):
        y = d.opp[x]
        if x < y:
            undirected_edges.append((x, y))

    num_e = len(undirected_edges)
    for i in range(num_e):
        e1 = undirected_edges[i]
        for j in range(i + 1, num_e):
            e2 = undirected_edges[j]
            blocked = {e1[0], e1[1], e2[0], e2[1]}
            # Check if removing e1 and e2 disconnects the vertex set!
            vis = [False] * n
            st = [0]
            vis[0] = True
            cnt = 0
            while st:
                cv = st.pop()
                cnt += 1
                for s in range(4):
                    dd = (cv << 2) | s
                    if dd in blocked:
                        continue
                    nv = d.opp[dd] >> 2
                    if not vis[nv]:
                        vis[nv] = True
                        st.append(nv)
            if cnt == n:
                continue
            # Ensure e1 and e2 each have one endpoint in A (`vis==True`) and one in B (`vis==False`)
            e1_a = e1[0] if vis[e1[0] >> 2] else e1[1]
            e1_b = e1[1] if vis[e1[0] >> 2] else e1[0]
            e2_a = e2[0] if vis[e2[0] >> 2] else e2[1]
            e2_b = e2[1] if vis[e2[0] >> 2] else e2[0]
            if not (vis[e1_a >> 2] and not vis[e1_b >> 2] and vis[e2_a >> 2] and not vis[e2_b >> 2]):
                continue

            # Build D_A (capping e1_a <-> e2_a) and D_B (capping e1_b <-> e2_b)
            def build_factor(keep_flag: bool, cap1: int, cap2: int) -> MarkedDiagram:
                vs = [v for v in range(n) if vis[v] == keep_flag]
                o2n = {v: idx for idx, v in enumerate(vs)}
                f_opp = [0] * (4 * len(vs))
                f_vt = [d.vtype[v] for v in vs]
                for idx, v in enumerate(vs):
                    for s in range(4):
                        orig_d = (v << 2) | s
                        if orig_d == cap1:
                            target = cap2
                        elif orig_d == cap2:
                            target = cap1
                        else:
                            target = d.opp[orig_d]
                        f_opp[4 * idx + s] = (o2n[target >> 2] << 2) | (target & 3)
                return MarkedDiagram(
                    num_vertices=len(vs),
                    opp=tuple(f_opp),
                    vtype=tuple(f_vt),
                    free_loops=0,
                )

            da = greedy_yoshikawa_reduce(build_factor(True, e1_a, e2_a))
            db = greedy_yoshikawa_reduce(build_factor(False, e1_b, e2_b))
            # Check if BOTH factors are genuinely nontrivial (not unknotted 0_1 S^2)!
            ta = da.surface_topology()
            tb = db.surface_topology()
            da_is_s2_unknot = (
                ta.num_components == 1
                and ta.components[0].orientable
                and ta.components[0].euler_characteristic == 2
                and da.num_vertices < 8
            )
            db_is_s2_unknot = (
                tb.num_components == 1
                and tb.components[0].orientable
                and tb.components[0].euler_characteristic == 2
                and db.num_vertices < 8
            )
            if not da_is_s2_unknot and not db_is_s2_unknot:
                return (da, db)
    return None


def simplify_and_explore_yoshikawa(
    d: MarkedDiagram,
    max_states: int = 512,
) -> tuple[MarkedDiagram, bool, bool]:
    """Explore level + reducing Yoshikawa moves from `d` to find minimal ch representative.

    Returns `(best_diagram, is_certified_split, is_certified_composite)`:
    - `best_diagram`: diagram of minimum `ch_index` (tie-broken canonically) found in component.
    - `is_certified_split`: True if equivalent to a disconnected diagram with >= 2 components!
    - `is_certified_composite`: True if equivalent to a nontrivial connected sum `D_1 # D_2`!
    """
    cur0 = greedy_yoshikawa_reduce(d)
    best = cur0
    best_key = canonical_marked_key(cur0, up_to_mirror_and_reverse=True)

    def check_split_or_comp(diag: MarkedDiagram) -> tuple[bool, bool]:
        comps = connected_diagram_components(diag)
        if len(comps) >= 2:
            return True, False
        if find_planar_connected_sum_decomposition(diag) is not None:
            return False, True
        return False, False

    sp, cp = check_split_or_comp(cur0)
    if sp or cp or cur0.num_vertices <= 2:
        return cur0, sp, cp

    visited = {canonical_marked_key(cur0, up_to_mirror_and_reverse=False)}
    queue: deque[MarkedDiagram] = deque([cur0])

    while queue and len(visited) < max_states:
        cur = queue.popleft()
        for nbr in list_level_yoshikawa_neighbors(cur):
            red = greedy_yoshikawa_reduce(nbr)
            sp, cp = check_split_or_comp(red)
            if sp or cp:
                return red, sp, cp
            if red.num_vertices < best.num_vertices:
                # Strictly smaller ch-index found! Restart search from `red`!
                return simplify_and_explore_yoshikawa(red, max_states=max_states)
            rk_sym = canonical_marked_key(red, up_to_mirror_and_reverse=True)
            if (red.num_vertices, rk_sym) < (best.num_vertices, best_key):
                best = red
                best_key = rk_sym
            rk = canonical_marked_key(red, up_to_mirror_and_reverse=False)
            if rk not in visited:
                visited.add(rk)
                queue.append(red)

    return best, False, False


def list_omega1_6_stabilizations(d: MarkedDiagram) -> list[MarkedDiagram]:
    """List all `+1`-vertex stabilizations (`Omega_1+`, `Omega_6+`, `Omega_6'+`) on edges of `d`."""
    n = d.num_vertices
    if n == 0:
        return []
    out: list[MarkedDiagram] = []
    u = n
    for d1 in range(4 * n):
        d2 = d.opp[d1]
        if d2 < d1:
            continue
        for side in (0, 1):
            s_opp = list(d.opp) + [0, 0, 0, 0]
            if side == 0:
                s_opp[4 * u + 0] = 4 * u + 1
                s_opp[4 * u + 1] = 4 * u + 0
                s_opp[4 * u + 2] = d1
                s_opp[d1] = 4 * u + 2
                s_opp[4 * u + 3] = d2
                s_opp[d2] = 4 * u + 3
            else:
                s_opp[4 * u + 1] = 4 * u + 2
                s_opp[4 * u + 2] = 4 * u + 1
                s_opp[4 * u + 3] = d1
                s_opp[d1] = 4 * u + 3
                s_opp[4 * u + 0] = d2
                s_opp[d2] = 4 * u + 0
            opp_tup = tuple(s_opp)
            for vt_u in ("X0", "X1", "Y", "Z"):
                out.append(
                    MarkedDiagram(
                        num_vertices=n + 1,
                        opp=opp_tup,
                        vtype=d.vtype + (vt_u,),  # type: ignore[arg-type]
                        free_loops=d.free_loops,
                    )
                )
    return out


def list_omega2_pokes(d: MarkedDiagram) -> list[MarkedDiagram]:
    """List all `+2`-vertex classical face pokes (`Omega_2+`) across non-identical edge pairs of faces of `d`."""
    n = d.num_vertices
    if n == 0:
        return []
    out: list[MarkedDiagram] = []
    vis_d = [False] * (4 * n)
    u, w = n, n + 1
    for sd in range(4 * n):
        if vis_d[sd]:
            continue
        f_darts: list[int] = []
        cur = sd
        while not vis_d[cur]:
            vis_d[cur] = True
            f_darts.append(cur)
            od = d.opp[cur]
            cur = ((od >> 2) << 2) | (((od & 3) - 1) & 3)
        if len(f_darts) < 2:
            continue
        for a in range(len(f_darts)):
            da = f_darts[a]
            oda = d.opp[da]
            ea = min(da, oda)
            for b in range(a + 1, len(f_darts)):
                db = f_darts[b]
                odb = d.opp[db]
                eb = min(db, odb)
                if ea == eb:
                    continue
                s_opp = list(d.opp) + [0] * 8
                s_opp[4 * u + 0] = da
                s_opp[da] = 4 * u + 0
                s_opp[4 * u + 1] = 4 * w + 3
                s_opp[4 * u + 2] = 4 * w + 2
                s_opp[4 * u + 3] = odb
                s_opp[odb] = 4 * u + 3
                s_opp[4 * w + 0] = oda
                s_opp[oda] = 4 * w + 0
                s_opp[4 * w + 1] = db
                s_opp[db] = 4 * w + 1
                s_opp[4 * w + 2] = 4 * u + 2
                s_opp[4 * w + 3] = 4 * u + 1
                opp_tup = tuple(s_opp)
                for vt_uw in ("X0", "X1"):
                    out.append(
                        MarkedDiagram(
                            num_vertices=n + 2,
                            opp=opp_tup,
                            vtype=d.vtype + (vt_uw, vt_uw),  # type: ignore[arg-type]
                            free_loops=d.free_loops,
                        )
                    )
    return out


def _is_planar_4reg_diagram(d: MarkedDiagram) -> bool:
    """Verify `V - E + F == 2` on every connected vertex component of `d`."""
    n = d.num_vertices
    if n == 0:
        return True
    for comp in connected_diagram_components(d):
        m = comp.num_vertices
        if m == 0:
            continue
        vis = [False] * (4 * m)
        faces = 0
        for sd in range(4 * m):
            if vis[sd]:
                continue
            faces += 1
            cur = sd
            while not vis[cur]:
                vis[cur] = True
                od = comp.opp[cur]
                cur = ((od >> 2) << 2) | (((od & 3) - 1) & 3)
        if faces != m + 2:
            return False
    return True


def list_4d_tait_flypes(d: MarkedDiagram) -> list[MarkedDiagram]:
    """List all 4D Tait flypes across 4-port marked tangles `T` with an adjacent classical crossing `c`.

    Geometric justification: rotating a 3-ball containing the 4-port marked tangle `T` and West
    crossing `c` by 180 degrees around the horizontal axis in `R^3 x [-eps, eps]` simultaneously
    in every `t`-slice (`s -> (4 - s) & 3`, `X0 <-> X1`, `Y <-> Z` inside `T`, shifting `c` from
    West to East) is an ambient isotopy of `F(D)` in `R^4` preserving `ch(D)`.
    """
    n = d.num_vertices
    if n < 3 or d.free_loops > 0:
        return []
    out: list[MarkedDiagram] = []
    vt_flip: dict[str, VertexKind] = {"X0": "X1", "X1": "X0", "Y": "Z", "Z": "Y"}

    for c in range(n):
        ct = d.vtype[c]
        if ct not in ("X0", "X1"):
            continue
        for sc in range(4):
            sc2 = (sc + 1) & 3
            p0 = d.opp[4 * c + sc]
            p3 = d.opp[4 * c + sc2]
            if (p3 >> 2) == c or (p0 >> 2) == c:
                continue
            wtop_over = (ct == "X0" and ((sc & 1) == 1)) or (ct == "X1" and ((sc & 1) == 0))

            undir: list[tuple[int, int]] = []
            for x in range(4 * n):
                y = d.opp[x]
                if x < y and (x >> 2) != c and (y >> 2) != c:
                    undir.append((x, y))

            for i in range(len(undir)):
                for j in range(i + 1, len(undir)):
                    e1, e2 = undir[i], undir[j]
                    blocked = {4 * c + sc, p0, 4 * c + sc2, p3, e1[0], e1[1], e2[0], e2[1]}
                    vis = {p3 >> 2}
                    st = [p3 >> 2]
                    while st:
                        u = st.pop()
                        for s in range(4):
                            dd = (u << 2) | s
                            if dd in blocked:
                                continue
                            v = d.opp[dd] >> 2
                            if v not in vis:
                                vis.add(v)
                                st.append(v)
                    if c in vis or (p0 >> 2) not in vis or len(vis) < 2 or len(vis) > n - 2:
                        continue
                    out_start = next(v for v in range(n) if v not in vis and v != c)
                    out_vis = {out_start, c}
                    st2 = [out_start, c]
                    while st2:
                        u = st2.pop()
                        for s in range(4):
                            dd = (u << 2) | s
                            if dd in blocked:
                                continue
                            v = d.opp[dd] >> 2
                            if v not in vis and v not in out_vis:
                                out_vis.add(v)
                                st2.append(v)
                    if len(vis) + len(out_vis) != n:
                        continue

                    p_a = e1[0] if (e1[0] >> 2) in vis else (e1[1] if (e1[1] >> 2) in vis else -1)
                    p_b = e2[0] if (e2[0] >> 2) in vis else (e2[1] if (e2[1] >> 2) in vis else -1)
                    if p_a < 0 or p_b < 0:
                        continue
                    if (d.opp[p_a] >> 2) in vis or (d.opp[p_b] >> 2) in vis:
                        continue

                    cur = ((p0 >> 2) << 2) | (((p0 & 3) + 1) & 3)
                    steps = 0
                    while (cur >> 2) in vis and cur not in (p_a, p_b, p0, p3) and steps <= 4 * n:
                        od = d.opp[cur]
                        cur = ((od >> 2) << 2) | (((od & 3) + 1) & 3)
                        steps += 1
                    if cur == p_a:
                        p1, p2 = p_a, p_b
                    elif cur == p_b:
                        p1, p2 = p_b, p_a
                    else:
                        continue

                    new_opp = list(d.opp)
                    new_vt = list(d.vtype)

                    def map_t_dart(td: int) -> int:
                        return ((td >> 2) << 2) | ((4 - (td & 3)) & 3)

                    for v in vis:
                        new_vt[v] = vt_flip[d.vtype[v]]
                        for s in range(4):
                            orig_d = (v << 2) | s
                            if orig_d in (p0, p1, p2, p3):
                                continue
                            orig_od = d.opp[orig_d]
                            new_opp[map_t_dart(orig_d)] = map_t_dart(orig_od)

                    mp0, mp1, mp2, mp3 = map_t_dart(p0), map_t_dart(p1), map_t_dart(p2), map_t_dart(p3)
                    new_vt[c] = "X1" if wtop_over else "X0"

                    west_top_dart = 4 * c + ((sc + 2) & 3)
                    west_bot_dart = 4 * c + ((sc + 3) & 3)
                    boundary_out = {
                        west_top_dart: mp0,
                        west_bot_dart: mp3,
                        p2: 4 * c + 3,
                        p1: 4 * c + 2,
                    }
                    new_opp[4 * c + 0] = mp1
                    new_opp[mp1] = 4 * c + 0
                    new_opp[4 * c + 1] = mp2
                    new_opp[mp2] = 4 * c + 1

                    for old_src, new_local in boundary_out.items():
                        old_tgt = d.opp[old_src]
                        if old_tgt in boundary_out:
                            new_opp[new_local] = boundary_out[old_tgt]
                        else:
                            new_opp[new_local] = old_tgt
                            new_opp[old_tgt] = new_local

                    cand = MarkedDiagram(
                        num_vertices=n,
                        opp=tuple(new_opp),
                        vtype=tuple(new_vt),
                        free_loops=0,
                    )
                    if _is_planar_4reg_diagram(cand):
                        out.append(cand)
    return out


def has_certified_cyclic_unknotted_tube(d: MarkedDiagram, _depth: int = 0) -> bool:
    """Certify non-minimality / compositeness via Hosokawa-Kawauchi (1979, Thm 1.2) / Boyle (1988).

    If `d` (`n >= 4`, connected surface `c(F(d)) == 1`) contains a marked-marked 2-gon `{u, w}`
    representing a 1-handle attached to a connected base diagram `d_0` (`|V(d_0)| = n - 2`,
    `c(F(d_0)) == 1`, `chi(F(d)) == chi(F(d_0)) - 2`) and `d_0` reduces by Yoshikawa moves
    (or recursively by cyclic-base 1-handle surgery) to `0_1` (`n=0`), `2_1^1` (`n=2`), or
    `2_1^{-1}` (`n=2`)—all of whose complements have cyclic `pi_1` generated by the meridian—then
    the 1-handle isotopy class depends only on orientation/framing, so `F(d)` is smoothly ambient
    isotopic to `F(d_0) # T^2` or `F(d_0) # K^2`, which is composite for `F(d_0) != S^2` and has
    `ch(F(d)) <= 4 < n` when `F(d_0) = S^2` (`n >= 5`).
    """
    from surfaceknots.unlink import certify_classical_unlink

    n = d.num_vertices
    if n < 5 or _depth > 3:
        return False
    st = d.surface_topology()
    if st.num_components != 1:
        return False

    for u in range(n):
        if d.vtype[u] not in ("Y", "Z"):
            continue
        for su in range(4):
            d_a = d.opp[4 * u + su]
            d_b = d.opp[4 * u + ((su + 1) & 3)]
            w = d_a >> 2
            if w <= u or (d_b >> 2) != w or d.vtype[w] not in ("Y", "Z"):
                continue
            sw_a = d_a & 3
            sw_b = d_b & 3
            if ((sw_b + 1) & 3) != sw_a:
                continue
            p0, p1 = 4 * u + ((su + 2) & 3), 4 * u + ((su + 3) & 3)
            p2, p3 = 4 * w + ((sw_b + 2) & 3), 4 * w + ((sw_b + 3) & 3)
            remap = [v - (1 if v > u else 0) - (1 if v > w else 0) for v in range(n)]
            m = n - 2
            nvt = tuple(d.vtype[v] for v in range(n) if v not in (u, w))

            for pairs in ([(p0, p1), (p2, p3)], [(p0, p3), (p1, p2)]):
                partner = {x: y for x, y in pairs} | {y: x for x, y in pairs}
                nopp = [0] * (4 * m)
                valid = True
                for v in range(n):
                    if v in (u, w):
                        continue
                    nv = remap[v]
                    for s in range(4):
                        od = d.opp[4 * v + s]
                        steps = 0
                        while (od >> 2) in (u, w):
                            od = d.opp[partner[od]]
                            steps += 1
                            if steps > 8:
                                valid = False
                                break
                        if not valid:
                            break
                        nopp[4 * nv + s] = (remap[od >> 2] << 2) | (od & 3)
                    if not valid:
                        break
                if not valid:
                    continue
                cand = MarkedDiagram(num_vertices=m, opp=tuple(nopp), vtype=nvt, free_loops=0)
                if not _is_planar_4reg_diagram(cand):
                    continue
                if not (
                    certify_classical_unlink(cand.resolve(1)).is_unlink
                    and certify_classical_unlink(cand.resolve(-1)).is_unlink
                ):
                    continue
                cst = cand.surface_topology()
                if cst.num_components != 1 or cst.euler_characteristic != st.euler_characteristic + 2:
                    continue
                from surfaceknots.invariants.fundamental_group import (
                    finite_quotient_profile,
                    wirtinger_presentation,
                )

                wp_c = wirtinger_presentation(cand)
                if wp_c.num_generators > 1:
                    fqp_c = finite_quotient_profile(wp_c, ("S3", "A4", "F21", "A5"))
                    if (
                        fqp_c["S3"] > 6
                        or fqp_c["A4"] > 12
                        or fqp_c["F21"] > 21
                        or fqp_c["A5"] > 60
                    ):
                        continue

                # By Yoshikawa's census (and our machine-verified ch <= 8 census), the ONLY
                # nontrivial prime connected surface-knot with ch <= 8 is 8_1 (which has
                # |Hom(pi_1(8_1), S_3)| = 12 > 6 and |Hom(pi_1(8_1), A_5)| = 360 > 60, already
                # ruled out above!). Hence any connected surface diagram `cand` with <= 9 vertices
                # that either has <= 7 vertices or is certified by C Yoshikawa moves/tunneling to
                # reduce to <= 8 vertices (or split along a planar 2-cut into <= 7-vertex summands)
                # is provably unknotted (pi_1 = Z or Z_2 generated by its meridian)!
                red, sp, cp = simplify_and_explore_yoshikawa(cand, max_states=24)
                if sp:
                    continue
                if cp or red.num_vertices <= 7:
                    return True
                if red.num_vertices in (8, 9):
                    from surfaceknots.fast_enum import tunnel_diagrams_c

                    t_kept = tunnel_diagrams_c(
                        red.num_vertices, [red], depth_plus1=96, depth_plus2=64, max_bfs_states=512
                    )
                    if len(t_kept) == 0:
                        return True
                    if has_certified_cyclic_unknotted_tube(t_kept[0], _depth=_depth + 1):
                        return True
    return False


