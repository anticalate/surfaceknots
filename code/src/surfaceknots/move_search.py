"""Constructive Yoshikawa move-path search, independent step-by-step path verification,
RP^2 reducibility / meridian-order certificates, and faithful 1-handle cord calculus
(Boyle 1988, Kamada 2014) for bases `8_1`, `9_1`, and `10_1`.
"""

from __future__ import annotations

from collections import deque
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from fractions import Fraction
import re
from typing import Callable

from surfaceknots.diagrams import MarkedDiagram, VertexKind
from surfaceknots.moves import (
    _list_marked_tube_flip_neighbors,
    _list_omega3_and_4_neighbors,
    _list_omega5_neighbors,
    _list_omega7_neighbors,
    _list_omega8_neighbors,
    canonical_marked_key,
    connected_diagram_components,
    find_omega1_or_6_reduction,
    find_omega2_reduction,
    find_planar_connected_sum_decomposition,
    greedy_yoshikawa_reduce,
    list_4d_tait_flypes,
    list_omega1_6_stabilizations,
    list_omega2_pokes,
)
from surfaceknots.unlink import certify_classical_unlink


# ============================================================================
# 1. Step-by-step Yoshikawa move-path validator (standalone parser, ribbon
#    canonicalizer, kink/R-II/port-tangle checks; candidate level moves from moves.py)
# ============================================================================


@dataclass(frozen=True)
class MovePathVerification:
    """Result of independent step-by-step verification of a Yoshikawa move sequence."""

    is_valid: bool
    num_steps: int
    step_kinds: tuple[str, ...]
    error: str = ""


def _indep_parse_epd(epd: str) -> tuple[int, tuple[int, ...], tuple[str, ...], int]:
    """Standalone EPD parser + planar Euler characteristic check (shares zero code with moves.py)."""
    s = epd.strip()
    if not s or s == "0_1":
        return 0, (), (), 1
    tok_re = re.compile(
        r"([XYZ])\s*\[\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\]|O\s*\[\s*(\d+)\s*\]"
    )
    items: list[tuple[str, int, int, int, int]] = []
    free_loops = 0
    pos = 0
    while pos < len(s):
        while pos < len(s) and s[pos] in " \t\r\n,;":
            pos += 1
        if pos >= len(s):
            break
        m = tok_re.match(s, pos)
        if m is None:
            raise ValueError(f"Invalid EPD token at {pos}: {epd!r}")
        if m.group(1) is not None:
            items.append((m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4)), int(m.group(5))))
        else:
            free_loops += 1
        pos = m.end()
    n = len(items)
    if n == 0:
        if free_loops <= 0:
            raise ValueError("Empty EPD")
        return 0, (), (), free_loops
    occ: dict[int, list[int]] = {}
    vtype: list[str] = []
    for v, (k, a, b, c, d) in enumerate(items):
        vtype.append("X0" if k == "X" else k)
        for sl, lbl in enumerate((a, b, c, d)):
            occ.setdefault(lbl, []).append(4 * v + sl)
    opp = [0] * (4 * n)
    for lbl, ds in occ.items():
        if len(ds) != 2:
            raise ValueError(f"Edge {lbl} occurs {len(ds)} times")
        opp[ds[0]] = ds[1]
        opp[ds[1]] = ds[0]
    for hd in range(4 * n):
        if opp[opp[hd]] != hd or opp[hd] == hd:
            raise ValueError("Non-involution opp")
    # Planar Euler check per graph component:
    vis_v = [False] * n
    comp_id = [-1] * n
    sizes: list[int] = []
    for sv in range(n):
        if not vis_v[sv]:
            cid = len(sizes)
            vis_v[sv] = True
            comp_id[sv] = cid
            stk = [sv]
            cnt = 0
            while stk:
                u = stk.pop()
                cnt += 1
                for sl in range(4):
                    w = opp[4 * u + sl] >> 2
                    if not vis_v[w]:
                        vis_v[w] = True
                        comp_id[w] = cid
                        stk.append(w)
            sizes.append(cnt)
    fvis = [False] * (4 * n)
    faces = [0] * len(sizes)
    for hd in range(4 * n):
        if not fvis[hd]:
            faces[comp_id[hd >> 2]] += 1
            cur = hd
            while not fvis[cur]:
                fvis[cur] = True
                od = opp[cur]
                cur = ((od >> 2) << 2) | (((od & 3) - 1) & 3)
    for cid, sz in enumerate(sizes):
        if faces[cid] != sz + 2:
            raise ValueError(f"Non-planar component V={sz} F={faces[cid]}")
    return n, tuple(opp), tuple(vtype), free_loops


def _indep_wirtinger_relations(
    n: int,
    opp: tuple[int, ...],
    vtype: tuple[str, ...],
    free_loops: int,
) -> tuple[int, dict[int, tuple[int, int]], list[tuple]]:
    """Independent Wirtinger relation builder on sheet-generators of `(n, opp, vtype)`."""
    if n == 0:
        return free_loops, {}, []
    parent = list(range(4 * n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: int, y: int) -> None:
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[ry] = rx

    # Orient edges arbitrarily, then build sheet generators on over-strands & marked vertices:
    # Actually even simpler: orientable/non-orientable Wirtinger from L_- resolution!
    # Every circle in L_- gives a generator? No: every ARC of L_- gives a Wirtinger generator,
    # and every saddle (marked vertex) equates merging/splitting arcs!
    for hd in range(4 * n):
        union(hd, opp[hd])
    for v in range(n):
        vt = vtype[v]
        if vt == "X0":
            union(4 * v + 1, 4 * v + 3)
        elif vt == "X1":
            union(4 * v + 0, 4 * v + 2)
        else:
            # In L_-, Y pairs (0,1)&(2,3), Z pairs (0,3)&(1,2); saddle adds 1 relation, while
            # the 4 darts around a marked vertex lie on the same sheet up to orientation!
            pass

    # Use L_- resolution directly: trace oriented circles of L_- to orient darts consistently!
    seen_d = [False] * (4 * n)
    dart_sign = [0] * (4 * n)  # +1 incoming to vertex, -1 outgoing from vertex
    for hd in range(4 * n):
        if seen_d[hd]:
            continue
        cur = hd
        while not seen_d[cur]:
            seen_d[cur] = True
            dart_sign[cur] = 1
            v, s = cur >> 2, cur & 3
            vt = vtype[v]
            if vt in ("X0", "X1"):
                ex = (v << 2) | (s ^ 2)
            elif vt == "Y":
                ex = (v << 2) | (s ^ 1)
            else:
                ex = (v << 2) | (3 - s)
            seen_d[ex] = True
            dart_sign[ex] = -1
            cur = opp[ex]

    # Over-strands in L_- connect across over-crossings; group darts into L_- over-arcs:
    arc_p = list(range(4 * n))

    def afind(x: int) -> int:
        while arc_p[x] != x:
            arc_p[x] = arc_p[arc_p[x]]
            x = arc_p[x]
        return x

    def aunion(x: int, y: int) -> None:
        rx, ry = afind(x), afind(y)
        if rx != ry:
            arc_p[ry] = rx

    for hd in range(4 * n):
        aunion(hd, opp[hd])
    for v in range(n):
        vt = vtype[v]
        if vt == "X0":
            aunion(4 * v + 1, 4 * v + 3)
        elif vt == "X1":
            aunion(4 * v + 0, 4 * v + 2)
        elif vt == "Y":
            aunion(4 * v + 0, 4 * v + 1)
            aunion(4 * v + 2, 4 * v + 3)
        else:
            aunion(4 * v + 0, 4 * v + 3)
            aunion(4 * v + 1, 4 * v + 2)

    roots = sorted({afind(hd) for hd in range(4 * n)})
    r2id = {r: i for i, r in enumerate(roots)}
    eid: dict[int, tuple[int, int]] = {}
    for hd in range(4 * n):
        # Meridian orientation via right-hand rule on L_- tangent `dart_sign`:
        eid[hd] = (r2id[afind(hd)], dart_sign[hd])

    rels: list[tuple] = []
    for v in range(n):
        vt = vtype[v]
        if vt in ("X0", "X1"):
            us = (0, 2) if vt == "X0" else (1, 3)
            os_ = (1, 3) if vt == "X0" else (0, 2)
            u_in = 4 * v + (us[0] if dart_sign[4 * v + us[0]] == 1 else us[1])
            u_out = 4 * v + (us[1] if dart_sign[4 * v + us[0]] == 1 else us[0])
            o_in = 4 * v + (os_[0] if dart_sign[4 * v + os_[0]] == 1 else os_[1])
            eps = 1 if (((o_in & 3) - (u_in & 3)) & 3) == 3 else -1
            # Relation: eid(u_out) * O * eid(u_in) * O^-1 = 1 (note eid(u_out) has dart_sign = -1!)
            rels.append(("conj", u_out, u_in, o_in, eps))
        else:
            # Marked saddle equates the two bands of L_+:
            if vt == "Y":
                # L_- pairs (0,1) and (2,3); saddle connects slot 0 to slot 3:
                rels.append(("eq", 4 * v + 0, 4 * v + 3))
            else:
                # L_- pairs (0,3) and (1,2); saddle connects slot 0 to slot 1:
                rels.append(("eq", 4 * v + 0, 4 * v + 1))

    return len(roots) + free_loops, eid, rels


def _indep_s3_hom_count(n: int, opp: tuple[int, ...], vtype: tuple[str, ...], free_loops: int) -> int:
    """Independent |Hom(pi_1, S_3)| counter directly from `_indep_wirtinger_relations`."""
    from itertools import product

    s3 = [(0, 1, 2), (1, 0, 2), (0, 2, 1), (2, 1, 0), (1, 2, 0), (2, 0, 1)]
    idx_of = {p: i for i, p in enumerate(s3)}
    mul = [[idx_of[tuple(s3[a][s3[b][k]] for k in range(3))] for b in range(6)] for a in range(6)]
    inv = [next(b for b in range(6) if mul[a][b] == 0) for a in range(6)]

    k, eid, rels = _indep_wirtinger_relations(n, opp, vtype, free_loops)
    if k > 10:
        # Use propagation from independent generators:
        pass
    # Propagate using eq relations first to reduce search variables:
    adj: list[list[tuple[int, int]]] = [[] for _ in range(k)]
    conj_rels = []
    for r in rels:
        if r[0] == "eq":
            g1, s1 = eid[r[1]]
            g2, s2 = eid[r[2]]
            # val(r[1]) * val(r[2]) == 1 => g1^{s1} = g2^{-s2}
            rel_exp = -s1 * s2  # g2 = g1^{rel_exp}
            adj[g1].append((g2, rel_exp))
            adj[g2].append((g1, rel_exp))
        else:
            conj_rels.append(r)

    comp_root = [-1] * k
    comp_exp = [1] * k
    roots: list[int] = []
    self_inv_roots: set[int] = set()
    for i in range(k):
        if comp_root[i] != -1:
            continue
        rid = len(roots)
        roots.append(i)
        comp_root[i] = rid
        comp_exp[i] = 1
        stk = [i]
        while stk:
            u = stk.pop()
            for w, rexp in adj[u]:
                want = comp_exp[u] * rexp
                if comp_root[w] == -1:
                    comp_root[w] = rid
                    comp_exp[w] = want
                    stk.append(w)
                elif comp_root[w] == rid and comp_exp[w] != want:
                    self_inv_roots.add(rid)

    if len(roots) > 9:
        return -1

    domain_for_r = [
        [g for g in range(6) if inv[g] == g] if rid in self_inv_roots else list(range(6))
        for rid in range(len(roots))
    ]
    cnt = 0
    for tup in product(*domain_for_r):
        gval = [tup[comp_root[i]] if comp_exp[i] == 1 else inv[tup[comp_root[i]]] for i in range(k)]
        ok = True
        for _, u_out, u_in, o_in, eps in conj_rels:
            go, so = eid[o_in]
            O = gval[go] if (so * eps > 0) else inv[gval[go]]
            gi, si = eid[u_in]
            X = gval[gi] if si > 0 else inv[gval[gi]]
            gu, su = eid[u_out]
            Y = gval[gu] if su > 0 else inv[gval[gu]]
            # Relation: Y * O * X * O^-1 == 1 => O * X * O^-1 == inv[Y]
            if mul[mul[O][X]][inv[O]] != inv[Y]:
                ok = False
                break
        if ok:
            cnt += 1
    return cnt


def _indep_canon_codes(
    n: int,
    opp: tuple[int, ...],
    vtype: tuple[str, ...],
    free_loops: int,
    allow_o4: bool,
) -> set[tuple]:
    """Independent canonical ribbon-graph codes across all connected graph components.

    When `allow_o4=False`, includes all 4 orientation-preserving `SO(4)` diagram symmetries
    (`id`, `(y,z)` half-turn, `(z,t)` half-turn `reverse_orientation_zt`, and `(y,t)` half-turn,
    all with `det = +1` in `SO(4)`). When `allow_o4=True`, includes all 8 `O(4)` mirror/reverse
    symmetries (`det = +-1`).
    """
    if n == 0:
        return {(0, free_loops, ())}

    seen_global = [False] * n
    comps: list[list[int]] = []
    for sv in range(n):
        if seen_global[sv]:
            continue
        comp_vs: list[int] = []
        stk = [sv]
        seen_global[sv] = True
        while stk:
            cv = stk.pop()
            comp_vs.append(cv)
            for s in range(4):
                nv = opp[(cv << 2) | s] >> 2
                if not seen_global[nv]:
                    seen_global[nv] = True
                    stk.append(nv)
        comps.append(comp_vs)

    transforms: list[tuple[bool, bool, bool]] = [
        (False, False, False),
        (True, True, False),
        (False, True, True),
        (True, False, True),
    ]
    if allow_o4:
        transforms = [
            (refl, flip_x, flip_yz)
            for refl in (False, True)
            for flip_x in (False, True)
            for flip_yz in (False, True)
        ]

    codes: set[tuple] = set()
    for refl, flip_x, flip_yz in transforms:
        s_step = -1 if refl else 1
        comp_bests: list[tuple] = []
        for comp_vs in comps:
            best_comp: tuple | None = None
            for v0 in comp_vs:
                for s0_init in range(4):
                    v2n = [-1] * n
                    ref_slot = [0] * n
                    v2n[v0] = 0
                    ref_slot[v0] = s0_init
                    q = [(v0 << 2) | s0_init]
                    nxt = 1
                    seq: list[object] = []
                    idx = 0
                    while idx < len(q):
                        sd = q[idx]
                        idx += 1
                        u, s0 = sd >> 2, sd & 3
                        vt_raw = vtype[u]
                        if vt_raw in ("X0", "X1"):
                            ut = ("X1" if vt_raw == "X0" else "X0") if flip_x else vt_raw
                            over0 = ((s0 & 1) == 0) if ut == "X1" else ((s0 & 1) == 1)
                            seq.append("Xo" if over0 else "Xu")
                        else:
                            ut = ("Z" if vt_raw == "Y" else "Y") if flip_yz else vt_raw
                            s1_raw = (s0 + s_step) & 3
                            p01 = ((s0 >> 1) == (s1_raw >> 1)) if ut == "Y" else ((s0 >> 1) != (s1_raw >> 1))
                            seq.append("M01" if p01 else "M03")
                        for k in range(4):
                            cd = (u << 2) | ((s0 + k * s_step) & 3)
                            od = opp[cd]
                            ov, os = od >> 2, od & 3
                            if v2n[ov] == -1:
                                v2n[ov] = nxt
                                ref_slot[ov] = os
                                nxt += 1
                                q.append(od)
                            rel_os = ((os - ref_slot[ov]) * s_step) & 3
                            seq.append(v2n[ov])
                            seq.append(rel_os)
                    tup = tuple(seq)
                    if best_comp is None or tup < best_comp:
                        best_comp = tup
            assert best_comp is not None
            comp_bests.append(best_comp)
        codes.add((n, free_loops, tuple(sorted(comp_bests))))
    return codes


def _indep_check_kink_splice(
    d_big: MarkedDiagram,
    d_small: MarkedDiagram,
    allow_o4: bool,
) -> str | None:
    """Independently check whether removing one adjacent self-loop vertex from `d_big` yields `d_small`."""
    n_big = d_big.num_vertices
    if n_big != d_small.num_vertices + 1:
        return None
    target_codes = _indep_canon_codes(
        d_small.num_vertices, d_small.opp, d_small.vtype, d_small.free_loops, allow_o4=allow_o4
    )
    for v in range(n_big):
        for s in range(4):
            da = (v << 2) | s
            db = (v << 2) | ((s + 1) & 3)
            if d_big.opp[da] != db:
                continue
            dc = (v << 2) | ((s + 2) & 3)
            dd = (v << 2) | ((s + 3) & 3)
            keep = [u for u in range(n_big) if u != v]
            o2n = {u: idx for idx, u in enumerate(keep)}
            new_vt = tuple(d_big.vtype[u] for u in keep)
            new_opp = [0] * (4 * len(keep))
            if d_big.opp[dc] == dd:
                new_fl = d_big.free_loops + 1
                for u in keep:
                    nu = o2n[u]
                    for su in range(4):
                        od = d_big.opp[(u << 2) | su]
                        new_opp[(nu << 2) | su] = (o2n[od >> 2] << 2) | (od & 3)
            else:
                new_fl = d_big.free_loops
                p = d_big.opp[dc]
                q = d_big.opp[dd]
                rp = (o2n[p >> 2] << 2) | (p & 3)
                rq = (o2n[q >> 2] << 2) | (q & 3)
                new_opp[rp] = rq
                new_opp[rq] = rp
                for u in keep:
                    nu = o2n[u]
                    for su in range(4):
                        orig = (u << 2) | su
                        if orig in (p, q):
                            continue
                        od = d_big.opp[orig]
                        new_opp[(nu << 2) | su] = (o2n[od >> 2] << 2) | (od & 3)
            cand_codes = _indep_canon_codes(
                len(keep), tuple(new_opp), new_vt, new_fl, allow_o4=allow_o4
            )
            if cand_codes & target_codes:
                return "OMEGA_1" if d_big.vtype[v] in ("X0", "X1") else "OMEGA_6"
    return None


def _indep_check_r2_splice(
    d_big: MarkedDiagram,
    d_small: MarkedDiagram,
    allow_o4: bool,
) -> bool:
    """Independently check whether removing one classical R-II bigon from `d_big` yields `d_small`."""
    n_big = d_big.num_vertices
    if n_big != d_small.num_vertices + 2:
        return False
    target_codes = _indep_canon_codes(
        d_small.num_vertices, d_small.opp, d_small.vtype, d_small.free_loops, allow_o4=allow_o4
    )
    for u in range(n_big):
        ut = d_big.vtype[u]
        if ut not in ("X0", "X1"):
            continue
        for su in range(4):
            od0 = d_big.opp[(u << 2) | su]
            w, sw1 = od0 >> 2, od0 & 3
            if w <= u:
                continue
            wt = d_big.vtype[w]
            if wt not in ("X0", "X1"):
                continue
            sw = (sw1 - 1) & 3
            su1 = (su + 1) & 3
            if d_big.opp[(u << 2) | su1] != ((w << 2) | sw):
                continue
            over_u = ((su1 & 1) == 1) if ut == "X0" else ((su1 & 1) == 0)
            over_w = ((sw & 1) == 1) if wt == "X0" else ((sw & 1) == 0)
            if over_u != over_w:
                continue

            a1 = (u << 2) | ((su + 3) & 3)
            b1 = (w << 2) | ((sw + 2) & 3)
            a2 = (u << 2) | ((su + 2) & 3)
            b2 = (w << 2) | ((sw + 3) & 3)
            thru = {a1: b1, b1: a1, a2: b2, b2: a2}

            keep = [v for v in range(n_big) if v not in (u, w)]
            o2n = {v: idx for idx, v in enumerate(keep)}
            new_vt = tuple(d_big.vtype[v] for v in keep)
            new_opp = [0] * (4 * len(keep))
            for v in keep:
                nv = o2n[v]
                for sv in range(4):
                    cur = d_big.opp[(v << 2) | sv]
                    while (cur >> 2) in (u, w):
                        cur = d_big.opp[thru[cur]]
                    new_opp[(nv << 2) | sv] = (o2n[cur >> 2] << 2) | (cur & 3)

            extra_loops = 0
            vis_int: set[int] = set()
            for start_x in (a1, b1, a2, b2):
                if start_x in vis_int:
                    continue
                comp_x: list[int] = []
                stk = [start_x]
                vis_int.add(start_x)
                closed = True
                while stk:
                    cx = stk.pop()
                    comp_x.append(cx)
                    tx = thru[cx]
                    if tx not in vis_int:
                        vis_int.add(tx)
                        stk.append(tx)
                    ox = d_big.opp[cx]
                    if (ox >> 2) in (u, w):
                        if ox not in vis_int:
                            vis_int.add(ox)
                            stk.append(ox)
                    else:
                        closed = False
                if closed:
                    extra_loops += 1

            cand_codes = _indep_canon_codes(
                len(keep), tuple(new_opp), new_vt, d_big.free_loops + extra_loops, allow_o4=allow_o4
            )
            if cand_codes & target_codes:
                return True
    return False


def _indep_verify_port_preserving_local_rewire(
    d1: MarkedDiagram,
    nb: MarkedDiagram,
    max_mod_vertices: int,
) -> bool:
    """Independently verify that `nb` modifies a vertex subset `mod_vs` of size `<= max_mod_vertices`
    in `d1`, leaves the untouched complement and its labeled port darts bijectively intact, and
    induces the exact same `L_+` and `L_-` port-to-port tangle matchings and internal loop counts.
    """
    from itertools import combinations

    if d1.num_vertices != nb.num_vertices or d1.free_loops != nb.free_loops:
        return False
    n = d1.num_vertices
    must_in = {v for v in range(n) if d1.vtype[v] != nb.vtype[v]}
    if len(must_in) > max_mod_vertices:
        return False
    changed_vs = [
        v
        for v in range(n)
        if v not in must_in
        and any(d1.opp[(v << 2) | s] != nb.opp[(v << 2) | s] for s in range(4))
    ]
    if not must_in and not changed_vs:
        return False

    def check_subset(mod_vs: set[int]) -> bool:
        return _indep_verify_port_tangle_subset(d1, nb, mod_vs)

    max_rem = min(len(changed_vs), max_mod_vertices - len(must_in))
    min_rem = 1 if not must_in else 0
    for r_sz in range(max_rem, min_rem - 1, -1):
        for extra in combinations(changed_vs, r_sz):
            if check_subset(must_in | set(extra)):
                return True
    return False


def _indep_verify_port_tangle_subset(
    d1: MarkedDiagram,
    nb: MarkedDiagram,
    mod_vs: set[int],
    max_ext_ports: int | None = None,
) -> bool:
    """Verify on a fixed vertex subset `mod_vs` that the complement and its labeled port darts
    are untouched and that `d1` and `nb` induce identical `L_+` and `L_-` port-to-port tangle
    matchings and internal loop counts.
    """
    n = d1.num_vertices
    ext_ports: list[int] = []
    for v in range(n):
        if v in mod_vs:
            continue
        for s in range(4):
            hd = (v << 2) | s
            o1 = d1.opp[hd]
            o2 = nb.opp[hd]
            in1 = (o1 >> 2) in mod_vs
            in2 = (o2 >> 2) in mod_vs
            if in1 != in2:
                return False
            if in1:
                ext_ports.append(hd)
            elif o1 != o2:
                return False

    if max_ext_ports is not None and len(ext_ports) > max_ext_ports:
        return False
    if not ext_ports:
        return True
    ext_set = set(ext_ports)

    def trace_ports(diag: MarkedDiagram, sign: int) -> tuple[dict[int, int], int] | None:
        mapping: dict[int, int] = {}
        used_internal: set[int] = set()

        def step_res(cur_d: int) -> int:
            u, s = cur_d >> 2, cur_d & 3
            vt = diag.vtype[u]
            if vt in ("X0", "X1"):
                return (u << 2) | (s ^ 2)
            if (vt == "Y" and sign == -1) or (vt == "Z" and sign == 1):
                return (u << 2) | (s ^ 1)
            return (u << 2) | (3 - s)

        for ep in ext_ports:
            cur = diag.opp[ep]
            steps = 0
            while (cur >> 2) in mod_vs:
                used_internal.add(cur)
                ex = step_res(cur)
                used_internal.add(ex)
                cur = diag.opp[ex]
                steps += 1
                if steps > 4 * len(mod_vs) + 2:
                    return None
            if cur not in ext_set:
                return None
            mapping[ep] = cur

        int_cycles = 0
        for v in mod_vs:
            for s in range(4):
                sd = (v << 2) | s
                if sd in used_internal:
                    continue
                int_cycles += 1
                cur = sd
                steps = 0
                while cur not in used_internal:
                    if (cur >> 2) not in mod_vs:
                        return None
                    used_internal.add(cur)
                    ex = step_res(cur)
                    used_internal.add(ex)
                    cur = diag.opp[ex]
                    steps += 1
                    if steps > 4 * len(mod_vs) + 2:
                        return None
        return mapping, int_cycles

    for sgn in (1, -1):
        m1 = trace_ports(d1, sgn)
        m2 = trace_ports(nb, sgn)
        if m1 is None or m2 is None or m1 != m2:
            return False
    return True


def _indep_verify_tait_flype_rewire(d1: MarkedDiagram, nb: MarkedDiagram) -> bool:
    """Independently verify that a 4D Tait flype candidate `d1 -> nb` decomposes into a 4-port
    sub-tangle `T` (`|T| >= 2`) undergoing exact half-turn (`s -> (4 - s) & 3`, `X0 <-> X1`,
    `Y <-> Z`) plus one adjacent classical crossing `c`, leaving the external complement untouched
    and preserving both `L_+` and `L_-` port-to-port tangle matchings and internal loop counts.
    """
    if d1.num_vertices != nb.num_vertices or d1.free_loops != 0 or nb.free_loops != 0:
        return False
    n = d1.num_vertices
    if n < 3:
        return False
    vt_flip = {"X0": "X1", "X1": "X0", "Y": "Z", "Z": "Y"}
    must_in = {v for v in range(n) if d1.vtype[v] != nb.vtype[v]}
    if len(must_in) < 2:
        return False
    changed_other = [
        v
        for v in range(n)
        if v not in must_in
        and any(d1.opp[(v << 2) | s] != nb.opp[(v << 2) | s] for s in range(4))
    ]
    cand_tc: list[tuple[set[int], int]] = []
    for c in changed_other:
        if d1.vtype[c] in ("X0", "X1") and nb.vtype[c] in ("X0", "X1"):
            cand_tc.append((must_in, c))
    if not changed_other:
        for c in range(n):
            if c not in must_in and d1.vtype[c] in ("X0", "X1") and nb.vtype[c] in ("X0", "X1"):
                cand_tc.append((must_in, c))
    for c in must_in:
        if d1.vtype[c] in ("X0", "X1") and nb.vtype[c] in ("X0", "X1") and len(must_in) >= 3:
            cand_tc.append((must_in - {c}, c))

    for t_set, c in cand_tc:
        mod_vs = t_set | {c}
        if len(t_set) < 2 or len(t_set) > n - 2:
            continue
        if any(nb.vtype[v] != vt_flip[d1.vtype[v]] for v in t_set):
            continue
        t_boundary = 0
        ok_t = True
        for v in t_set:
            for s in range(4):
                hd = (v << 2) | s
                od = d1.opp[hd]
                if (od >> 2) in t_set:
                    m_hd = (v << 2) | ((4 - s) & 3)
                    m_od = ((od >> 2) << 2) | ((4 - (od & 3)) & 3)
                    if nb.opp[m_hd] != m_od:
                        ok_t = False
                        break
                else:
                    t_boundary += 1
            if not ok_t:
                break
        if not ok_t or t_boundary != 4:
            continue
        if _indep_verify_port_tangle_subset(d1, nb, mod_vs, max_ext_ports=4):
            return True
    return False


def _indep_is_single_local_move(
    d1: MarkedDiagram,
    d2: MarkedDiagram,
    allow_o4: bool = True,
) -> str | None:
    """Classify a single step `d1 -> d2`.

    Standalone checks (sharing zero code with `moves.py`):
    - Ribbon-graph symmetry canonicalization (`_indep_canon_codes`) for `SO(4)` (`allow_o4=False`,
      including `reverse_orientation_zt`) and `O(4)` (`allow_o4=True`, up to ambient mirror).
    - Exact 1-vertex self-loop kink splice (`Omega_1`, `Omega_6`) and 2-vertex R-II bigon splice
      (`Omega_2`).
    - Port-labeled `L_+`/`L_-` tangle-matching & internal-loop verification on every `Delta V = 0`
      move (`_indep_verify_port_preserving_local_rewire` and `_indep_verify_tait_flype_rewire`).

    Candidate enumeration for `Delta V = 0` planar templates (`Omega_3`, `Omega_4`, `Omega_5`,
    `Omega_7`, `Omega_8`, `MARKED_TUBE_FLIP`, `TAIT_FLYPE_4D`) uses the generator functions in
    `moves.py`, followed by standalone port-tangle and canonical-code verification.
    """
    n1, opp1, vt1, fl1 = d1.num_vertices, d1.opp, d1.vtype, d1.free_loops
    n2, opp2, vt2, fl2 = d2.num_vertices, d2.opp, d2.vtype, d2.free_loops
    diff = n2 - n1
    if abs(diff) > 2:
        return None

    if diff == 0 and fl1 == fl2:
        c1 = _indep_canon_codes(n1, opp1, vt1, fl1, allow_o4=False)
        c2 = _indep_canon_codes(n2, opp2, vt2, fl2, allow_o4=False)
        if c1 & c2:
            return "ISO_SO4"
        if allow_o4:
            co1 = _indep_canon_codes(n1, opp1, vt1, fl1, allow_o4=True)
            co2 = _indep_canon_codes(n2, opp2, vt2, fl2, allow_o4=True)
            if co1 & co2:
                return "MIRROR_O4"

    if diff == -1:
        kink = _indep_check_kink_splice(d1, d2, allow_o4=allow_o4)
        return f"{kink}_MINUS" if kink is not None else None
    if diff == 1:
        kink = _indep_check_kink_splice(d2, d1, allow_o4=allow_o4)
        return f"{kink}_PLUS" if kink is not None else None
    if diff == -2:
        return "OMEGA_2_MINUS" if _indep_check_r2_splice(d1, d2, allow_o4=allow_o4) else None
    if diff == 2:
        return "OMEGA_2_PLUS" if _indep_check_r2_splice(d2, d1, allow_o4=allow_o4) else None

    # diff == 0: check exact port-labeled planar templates + independent port-tangle verification
    if fl1 != fl2:
        return None
    target_codes = _indep_canon_codes(n2, opp2, vt2, fl2, allow_o4=allow_o4)

    template_families: tuple[tuple[str, Callable[[MarkedDiagram], list[MarkedDiagram]], int], ...] = (
        ("OMEGA_3_4", _list_omega3_and_4_neighbors, 3),
        ("OMEGA_5", _list_omega5_neighbors, 2),
        ("OMEGA_7", _list_omega7_neighbors, 2),
        ("OMEGA_8", _list_omega8_neighbors, 6),
        ("MARKED_TUBE_FLIP", _list_marked_tube_flip_neighbors, 2),
    )
    for tag, gen_fn, max_k in template_families:
        for nb in gen_fn(d1):
            nb_codes = _indep_canon_codes(
                nb.num_vertices, nb.opp, nb.vtype, nb.free_loops, allow_o4=allow_o4
            )
            if (nb_codes & target_codes) and _indep_verify_port_preserving_local_rewire(d1, nb, max_k):
                return tag

    for nb in list_4d_tait_flypes(d1):
        nb_codes = _indep_canon_codes(
            nb.num_vertices, nb.opp, nb.vtype, nb.free_loops, allow_o4=allow_o4
        )
        if (nb_codes & target_codes) and _indep_verify_tait_flype_rewire(d1, nb):
            return "TAIT_FLYPE_4D"

    return None


def verify_move_path_independent(
    epd_sequence: list[str] | tuple[str, ...],
    allow_o4: bool = True,
) -> MovePathVerification:
    """Validate an EPD move path step-by-step with port-labeled local move verification and full invariant checks.

    Parameters
    ----------
    allow_o4 : bool, default True
        When `True`, validates topological equivalence in `R^4` up to `O(4)` ambient reflection
        (mirror images `MIRROR_O4` allowed; normal Euler numbers compared up to global sign).
        When `False`, enforces orientation-preserving `SO(4)` ambient isotopy (`ISO_SO4`, including
        `reverse_orientation_zt`, with signed normal Euler numbers preserved).

    Verifies:
    1. Every EPD parses cleanly via standalone `_indep_parse_epd` with valid involution `opp` and genus-0 planar faces.
    2. Every intermediate diagram is admissible (`L_+` and `L_-` certified classical unlinks).
    3. Topological, Alexander, and peripheral invariants (`num_components`, `yoshikawa_genera`,
       `euler_characteristic`, `normal_euler_numbers`, `H_1`, finite quotient counts,
       canonical Alexander `evaluation_fingerprint`, peripheral `SL(2, F_3)` count, and standalone
       `|Hom(pi_1, S_3)|`) remain strictly invariant at every step.
    4. Every transition `epd[i] -> epd[i+1]` passes `_indep_is_single_local_move` (standalone kink/R-II
       splices or `moves.py` `Delta V = 0` candidate templates verified by standalone port-tangle &
       ribbon-code checks).
    """
    if not epd_sequence:
        return MovePathVerification(False, 0, (), "Empty path")

    diags: list[MarkedDiagram] = []
    indep_s3: list[int] = []
    for idx, epd in enumerate(epd_sequence):
        try:
            n, opp, vt, fl = _indep_parse_epd(epd)
            d = MarkedDiagram.from_epd(epd)
        except Exception as exc:
            return MovePathVerification(False, idx, (), f"Step {idx} parse/planarity error: {exc}")
        cp = certify_classical_unlink(d.resolve(1))
        cm = certify_classical_unlink(d.resolve(-1))
        if not (cp.is_unlink and cm.is_unlink):
            return MovePathVerification(
                False,
                idx,
                (),
                f"Step {idx} inadmissible: L+={cp.reason}, L-={cm.reason}",
            )
        diags.append(d)
        indep_s3.append(_indep_s3_hom_count(n, opp, vt, fl))

    from surfaceknots.invariants.alexander import compute_alexander_invariants
    from surfaceknots.invariants.fundamental_group import (
        finite_quotient_profile,
        tietze_simplify,
        wirtinger_presentation,
    )
    from surfaceknots.invariants.linking import peripheral_sl2f3_profile

    def quick_inv(d: MarkedDiagram) -> tuple:
        topo = d.surface_topology()
        es = [c.normal_euler for c in topo.components]
        if allow_o4:
            e_key = (
                tuple(sorted(abs(x) for x in es))
                if topo.num_components == 1
                else tuple(min(sorted(es), sorted(-x for x in es)))
            )
        else:
            e_key = tuple(sorted(es))
        wp = tietze_simplify(wirtinger_presentation(d))
        h1 = wp.abelianization_invariant_factors()
        fqp = finite_quotient_profile(wp)
        alex_fp = compute_alexander_invariants(d, pres=wp).evaluation_fingerprint
        sl2_tot, sl2_periph = peripheral_sl2f3_profile(d)
        return (
            topo.num_components,
            topo.yoshikawa_genera,
            topo.euler_characteristic,
            e_key,
            h1,
            tuple(sorted(fqp.items())),
            alex_fp,
            sl2_tot,
            len(sl2_periph),
        )

    base_inv = quick_inv(diags[0])
    base_s3 = indep_s3[0]
    step_kinds: list[str] = []
    for idx in range(len(diags) - 1):
        cur_inv = quick_inv(diags[idx + 1])
        if cur_inv != base_inv:
            return MovePathVerification(
                False,
                idx + 1,
                tuple(step_kinds),
                f"Invariant changed at step {idx}->{idx+1}: {base_inv} != {cur_inv}",
            )
        if base_s3 != -1 and indep_s3[idx + 1] != -1 and indep_s3[idx + 1] != base_s3:
            return MovePathVerification(
                False,
                idx + 1,
                tuple(step_kinds),
                f"Independent |Hom(pi_1, S_3)| changed at step {idx}->{idx+1}: {base_s3} != {indep_s3[idx+1]}",
            )
        kind = _indep_is_single_local_move(diags[idx], diags[idx + 1], allow_o4=allow_o4)
        if kind is None:
            return MovePathVerification(
                False,
                idx + 1,
                tuple(step_kinds),
                f"Unrecognized non-local transition at step {idx}->{idx+1}",
            )
        step_kinds.append(kind)

    return MovePathVerification(True, len(step_kinds), tuple(step_kinds))


# ============================================================================
# 2. Constructive single-step move expansion & parallel bidirectional search
# ============================================================================


def list_labeled_single_moves(
    d: MarkedDiagram,
    max_ch: int,
    include_stabilizations: bool = True,
    include_flypes: bool = True,
) -> list[tuple[str, MarkedDiagram]]:
    """Return admissible single-step move neighbors `(move_tag, nbr)` of `d`."""
    out: list[tuple[str, MarkedDiagram]] = []

    def _adm(cand: MarkedDiagram) -> bool:
        return (
            certify_classical_unlink(cand.resolve(1)).is_unlink
            and certify_classical_unlink(cand.resolve(-1)).is_unlink
        )

    r16 = find_omega1_or_6_reduction(d)
    if r16 is not None and _adm(r16):
        out.append(("Omega_1_6-", r16))
    r2 = find_omega2_reduction(d)
    if r2 is not None and _adm(r2):
        out.append(("Omega_2-", r2))

    for nb in _list_omega3_and_4_neighbors(d):
        if _adm(nb):
            out.append(("Omega_3_4", nb))
    for nb in _list_omega5_neighbors(d):
        if _adm(nb):
            out.append(("Omega_5", nb))
    for nb in _list_omega7_neighbors(d):
        if _adm(nb):
            out.append(("Omega_7", nb))
    for nb in _list_omega8_neighbors(d):
        if _adm(nb):
            out.append(("Omega_8", nb))
    for nb in _list_marked_tube_flip_neighbors(d):
        if _adm(nb):
            out.append(("TubeFlip", nb))
    if include_flypes:
        for nb in list_4d_tait_flypes(d):
            if _adm(nb):
                out.append(("TaitFlype4D", nb))

    if include_stabilizations:
        if d.num_vertices + 1 <= max_ch:
            for nb in list_omega1_6_stabilizations(d):
                if _adm(nb):
                    out.append(("Omega_1_6+", nb))
        if d.num_vertices + 2 <= max_ch:
            for nb in list_omega2_pokes(d):
                if _adm(nb):
                    out.append(("Omega_2+", nb))

    return out


def _greedy_path_from(d: MarkedDiagram) -> tuple[MarkedDiagram, list[tuple[str, str]]]:
    """Greedily apply `Omega_1/6-` and `Omega_2-` while recording every intermediate `(move, epd)`."""
    cur = d
    steps: list[tuple[str, str]] = []
    while cur.num_vertices > 0:
        n1 = find_omega1_or_6_reduction(cur)
        if n1 is not None:
            steps.append(("Omega_1_6-", n1.to_epd()))
            cur = n1
            continue
        n2 = find_omega2_reduction(cur)
        if n2 is not None:
            steps.append(("Omega_2-", n2.to_epd()))
            cur = n2
            continue
        break
    return cur, steps


def _expand_wave_chunk(
    args: tuple[list[tuple[tuple, str]], int, bool, bool],
) -> list[tuple[tuple, tuple, list[tuple[str, str]], int]]:
    """Worker helper for parallel bidirectional BFS wave expansion."""
    chunk, max_ch, allow_o4, include_stabs = args
    results: list[tuple[tuple, tuple, list[tuple[str, str]], int]] = []
    seen_local: set[tuple] = set()
    for src_key, src_epd in chunk:
        d = MarkedDiagram.from_epd(src_epd)
        for mv_tag, nb in list_labeled_single_moves(
            d,
            max_ch=max_ch,
            include_stabilizations=include_stabs,
            include_flypes=True,
        ):
            nb_epd = nb.to_epd()
            red, red_steps = _greedy_path_from(nb)
            # Record both the immediate neighbor (if stabilized) and its greedy reduction:
            for target_d, path_seg in (
                [(nb, [(mv_tag, nb_epd)]), (red, [(mv_tag, nb_epd)] + red_steps)]
                if red_steps
                else [(nb, [(mv_tag, nb_epd)])]
            ):
                tk = canonical_marked_key(target_d, up_to_mirror_and_reverse=allow_o4)
                if tk not in seen_local:
                    seen_local.add(tk)
                    results.append((src_key, tk, path_seg, target_d.num_vertices))
    return results


@dataclass(frozen=True)
class BidirectionalSearchResult:
    """Outcome of a 64-core parallel bidirectional mountain-pass Yoshikawa move search."""

    connected: bool
    max_ch: int
    waves_run: int
    forward_states: int
    backward_states: int
    epd_path: tuple[str, ...] = ()
    move_tags: tuple[str, ...] = ()
    independent_verification: MovePathVerification | None = None


def parallel_bidirectional_move_search(
    start_epds: list[str],
    goal_epds: list[str],
    *,
    max_ch: int = 15,
    max_waves: int = 8,
    beam_width: int = 20000,
    workers: int = 16,
    allow_o4: bool = True,
) -> BidirectionalSearchResult:
    """Search for an explicit Yoshikawa move path connecting `start_epds` to `goal_epds` up to `max_ch`."""
    f_vis: dict[tuple, tuple[tuple | None, list[tuple[str, str]], str]] = {}
    b_vis: dict[tuple, tuple[tuple | None, list[tuple[str, str]], str]] = {}
    f_front: list[tuple[tuple, str]] = []
    b_front: list[tuple[tuple, str]] = []

    for epd in start_epds:
        d = MarkedDiagram.from_epd(epd)
        k = canonical_marked_key(d, up_to_mirror_and_reverse=allow_o4)
        if k not in f_vis:
            f_vis[k] = (None, [], d.to_epd())
            f_front.append((k, d.to_epd()))

    for epd in goal_epds:
        d = MarkedDiagram.from_epd(epd)
        k = canonical_marked_key(d, up_to_mirror_and_reverse=allow_o4)
        if k not in b_vis:
            b_vis[k] = (None, [], d.to_epd())
            b_front.append((k, d.to_epd()))

    def reconstruct(meet_key: tuple) -> BidirectionalSearchResult:
        # Forward chain: root -> ... -> meet_key
        f_segs: list[list[tuple[str, str]]] = []
        cur: tuple | None = meet_key
        root_epd = ""
        while cur is not None:
            prev_k, seg, rep_epd = f_vis[cur]
            if prev_k is None:
                root_epd = rep_epd
            else:
                f_segs.append(seg)
            cur = prev_k
        f_segs.reverse()
        epds = [root_epd]
        tags: list[str] = []
        for seg in f_segs:
            for tg, e_str in seg:
                tags.append(tg)
                epds.append(e_str)

        # Backward chain: meet_key -> ... -> goal_root
        b_chains: list[tuple[str, list[tuple[str, str]]]] = []
        cur_b: tuple | None = meet_key
        while cur_b is not None:
            prev_b, seg_b, rep_b = b_vis[cur_b]
            b_chains.append((rep_b, seg_b))
            cur_b = prev_b
        # b_chains[0] has rep_b isomorphic (mod O(4)) to epds[-1]!
        # For each step i -> i+1 in b_chains, in the backward tree `prev_b` expanded to `cur_b` via `seg_b`:
        # so reversing `seg_b` walks from `cur_b`'s `rep_b` back to `prev_b`'s `rep_b`!
        if b_chains[0][0] != epds[-1]:
            tags.append("ISO_JOIN")
            epds.append(b_chains[0][0])
        for idx in range(len(b_chains) - 1):
            _, seg_from_parent = b_chains[idx]
            parent_rep, _ = b_chains[idx + 1]
            # Intermediate EPDs in seg_from_parent were parent_rep -> s_0 -> s_1 -> ... -> child_rep:
            inter_epds = [parent_rep] + [e_str for _, e_str in seg_from_parent[:-1]]
            inter_tags = [f"rev({tg})" for tg, _ in reversed(seg_from_parent)]
            for tg, target_e in zip(inter_tags, reversed(inter_epds)):
                tags.append(tg)
                epds.append(target_e)

        ver = verify_move_path_independent(epds, allow_o4=allow_o4)
        return BidirectionalSearchResult(
            connected=True,
            max_ch=max_ch,
            waves_run=waves_done,
            forward_states=len(f_vis),
            backward_states=len(b_vis),
            epd_path=tuple(epds),
            move_tags=tuple(tags),
            independent_verification=ver,
        )

    waves_done = 0
    inter0 = set(f_vis.keys()) & set(b_vis.keys())
    if inter0:
        return reconstruct(next(iter(inter0)))

    import multiprocessing as mp

    with ProcessPoolExecutor(max_workers=max(1, workers), mp_context=mp.get_context("forkserver")) as ex:
        for wave in range(1, max_waves + 1):
            if not f_front or not b_front:
                break
            waves_done = wave
            expand_forward = len(f_front) <= len(b_front)
            active_front = f_front if expand_forward else b_front
            active_vis = f_vis if expand_forward else b_vis
            other_vis = b_vis if expand_forward else f_vis

            if len(active_front) > beam_width:
                active_front.sort(key=lambda item: (MarkedDiagram.from_epd(item[1]).num_vertices, item[1]))
                active_front = active_front[:beam_width]

            chunk_sz = max(1, (len(active_front) + workers * 4 - 1) // (workers * 4))
            tasks = [
                (active_front[i : i + chunk_sz], max_ch, allow_o4, True)
                for i in range(0, len(active_front), chunk_sz)
            ]
            next_front_candidates: list[tuple[int, tuple, str]] = []
            meet_found: tuple | None = None
            for batch in ex.map(_expand_wave_chunk, tasks):
                for src_k, tk, path_seg, nv in batch:
                    if tk not in active_vis:
                        rep_epd = path_seg[-1][1]
                        active_vis[tk] = (src_k, path_seg, rep_epd)
                        next_front_candidates.append((nv, tk, rep_epd))
                        if tk in other_vis and meet_found is None:
                            meet_found = tk
            if meet_found is not None:
                return reconstruct(meet_found)

            next_front_candidates.sort(key=lambda x: (x[0], x[2]))
            next_front = [(tk, rep) for _, tk, rep in next_front_candidates[:beam_width]]
            if expand_forward:
                f_front = next_front
            else:
                b_front = next_front

    return BidirectionalSearchResult(
        connected=False,
        max_ch=max_ch,
        waves_run=waves_done,
        forward_states=len(f_vis),
        backward_states=len(b_vis),
    )


# ============================================================================
# 3. Constructive RP^2 reducibility & meridian-order certificates (AUDIT §9.4)
# ============================================================================


@dataclass(frozen=True)
class RP2ReducibilityCertificate:
    """Constructive reducibility / irreducibility analysis for a connected RP^2 diagram."""

    status: str  # "STANDARD_P2", "REDUCED_LOWER_CH", "SPLIT_CONNECTED_SUM", "CERTIFIED_IRREDUCIBLE", "UNRESOLVED"
    start_ch: int
    terminal_ch: int
    epd_path: tuple[str, ...]
    factor_epds: tuple[str, ...] = ()
    meridian_order_lower_bound: int = 2
    meridian_witness_group: str = ""
    path_verified: bool = False


def meridian_order_lower_bound(d: MarkedDiagram) -> tuple[int, str]:
    """Compute a rigorous lower bound on the order of meridians in `pi_1(S^4 \\ F)` via finite quotients.

    By Yoshikawa (1994) Prop. 5.3 and Hughes-Kim-Miller-Nahm (2026), any connected surface with a
    standard `P_+` or `P_-` connected-sum factor (or unknotted `K^2` factor) has a meridian of order 2.
    Thus `meridian_order > 2` on a connected `RP^2` (or single-meridian-conjugacy-class surface)
    immediately certifies irreducibility!
    """
    from sympy.combinatorics.fp_groups import FpGroup, coset_enumeration_r
    from sympy.combinatorics.free_groups import free_group
    from surfaceknots.invariants.fundamental_group import wirtinger_presentation

    wp = wirtinger_presentation(d)
    if wp.num_generators <= 1:
        return (2 if d.surface_topology().orientability_label != "Ori" else 0), "abelian"

    # Test if adding `x_1^2 = 1` gives order > 2 AND check order of `x_1` in `G / <<x_1^k>>`!
    g_names = ", ".join(f"x{i+1}" for i in range(wp.num_generators))
    fg_tuple = free_group(g_names)
    F = fg_tuple[0]
    fgens = fg_tuple[1:]
    gmap = {}
    for i, fg in enumerate(fgens, start=1):
        gmap[i] = fg
        gmap[-i] = fg**-1

    rels = []
    for r in wp.relators:
        w = F.identity
        for lit in r:
            w = w * gmap[lit]
        rels.append(w)

    # First check if `G` itself is small finite (like `2O` of order 48 on `11_1^{-2}`!):
    try:
        G0 = FpGroup(F, rels)
        ct0 = coset_enumeration_r(G0, [], max_cosets=2000)
        ct0.compress()
        # Trace permutation of generator x1 on cosets to get its exact order in G!
        perm = [ct0.table[c][0] for c in range(len(ct0.table))]
        cur = 0
        ord_x1 = 1
        cur = perm[0]
        while cur != 0 and ord_x1 <= len(perm) + 1:
            cur = perm[cur]
            ord_x1 += 1
        return ord_x1, f"exact_pi1(order={len(ct0.table)})"
    except Exception:
        pass

    best_ord = 2
    best_src = "H_1"
    for k in (4, 6, 8):
        try:
            Gk = FpGroup(F, rels + [fgens[0] ** k])
            ct = coset_enumeration_r(Gk, [], max_cosets=1500)
            ct.compress()
            perm = [ct.table[c][0] for c in range(len(ct.table))]
            cur = perm[0]
            ord_k = 1
            while cur != 0 and ord_k <= k + 1:
                cur = perm[cur]
                ord_k += 1
            if ord_k > best_ord:
                best_ord = ord_k
                best_src = f"G/<<x1^{k}>>(order={len(ct.table)})"
        except Exception:
            continue
    return best_ord, best_src


def _find_any_2edge_cut_factors(d: MarkedDiagram) -> tuple[MarkedDiagram, MarkedDiagram] | None:
    """Find ANY planar 2-edge cut splitting `d` into two positive-vertex factors `(d_a, d_b)`
    (specifically catching `K # P_{+-}` where `P_{+-}` has only 2 vertices!).
    """
    n = d.num_vertices
    if n < 4 or len(connected_diagram_components(d)) != 1:
        return None
    undir = [(x, d.opp[x]) for x in range(4 * n) if x < d.opp[x]]
    for i in range(len(undir)):
        e1 = undir[i]
        for j in range(i + 1, len(undir)):
            e2 = undir[j]
            blocked = {e1[0], e1[1], e2[0], e2[1]}
            vis = [False] * n
            vis[0] = True
            st = [0]
            cnt = 0
            while st:
                u = st.pop()
                cnt += 1
                for s in range(4):
                    dd = (u << 2) | s
                    if dd in blocked:
                        continue
                    w = d.opp[dd] >> 2
                    if not vis[w]:
                        vis[w] = True
                        st.append(w)
            if cnt < 2 or cnt > n - 2:
                continue
            e1_a = e1[0] if vis[e1[0] >> 2] else e1[1]
            e1_b = e1[1] if vis[e1[0] >> 2] else e1[0]
            e2_a = e2[0] if vis[e2[0] >> 2] else e2[1]
            e2_b = e2[1] if vis[e2[0] >> 2] else e2[0]
            if not (vis[e1_a >> 2] and not vis[e1_b >> 2] and vis[e2_a >> 2] and not vis[e2_b >> 2]):
                continue

            def mk_fac(flag: bool, c1: int, c2: int) -> MarkedDiagram:
                vs = [v for v in range(n) if vis[v] == flag]
                o2n = {v: idx for idx, v in enumerate(vs)}
                f_opp = [0] * (4 * len(vs))
                f_vt = [d.vtype[v] for v in vs]
                for idx, v in enumerate(vs):
                    for s in range(4):
                        od = (v << 2) | s
                        tgt = c2 if od == c1 else (c1 if od == c2 else d.opp[od])
                        f_opp[4 * idx + s] = (o2n[tgt >> 2] << 2) | (tgt & 3)
                return MarkedDiagram(len(vs), tuple(f_opp), tuple(f_vt), 0)

            da = mk_fac(True, e1_a, e2_a)
            db = mk_fac(False, e1_b, e2_b)
            # Reject trivial S^2 bubble cuts where one side is an unknotted S^2:
            ta, tb = da.surface_topology(), db.surface_topology()
            da_red = greedy_yoshikawa_reduce(da)
            db_red = greedy_yoshikawa_reduce(db)
            if da_red.num_vertices == 0 or db_red.num_vertices == 0:
                continue
            if (ta.euler_characteristic == 1 and tb.euler_characteristic == 2 and db_red.num_vertices >= 8) or (
                tb.euler_characteristic == 1 and ta.euler_characteristic == 2 and da_red.num_vertices >= 8
            ):
                return (da_red, db_red)
    return None


def find_rp2_reducible_certificate(
    d: MarkedDiagram,
    max_ch: int | None = None,
    max_states: int = 3000,
) -> RP2ReducibilityCertificate:
    """Constructively certify an `RP^2` diagram as `2_1^{+-1}`, reducible to `< ch`, split `K # P_{+-}`,
    or irreducible via meridian order > 2.
    """
    start_ch = d.num_vertices
    start_epd = d.to_epd()
    m_ord, m_src = meridian_order_lower_bound(d)
    if m_ord > 2:
        return RP2ReducibilityCertificate(
            status="CERTIFIED_IRREDUCIBLE",
            start_ch=start_ch,
            terminal_ch=start_ch,
            epd_path=(start_epd,),
            meridian_order_lower_bound=m_ord,
            meridian_witness_group=m_src,
            path_verified=True,
        )

    if start_ch <= 2:
        return RP2ReducibilityCertificate(
            status="STANDARD_P2",
            start_ch=start_ch,
            terminal_ch=start_ch,
            epd_path=(start_epd,),
            meridian_order_lower_bound=m_ord,
            meridian_witness_group=m_src,
            path_verified=True,
        )

    red0, steps0 = _greedy_path_from(d)
    if red0.num_vertices < start_ch:
        epds = [start_epd] + [e for _, e in steps0]
        ver = verify_move_path_independent(epds)
        return RP2ReducibilityCertificate(
            status="STANDARD_P2" if red0.num_vertices <= 2 else "REDUCED_LOWER_CH",
            start_ch=start_ch,
            terminal_ch=red0.num_vertices,
            epd_path=tuple(epds),
            meridian_order_lower_bound=m_ord,
            meridian_witness_group=m_src,
            path_verified=ver.is_valid,
        )

    cut0 = _find_any_2edge_cut_factors(d)
    if cut0 is not None:
        return RP2ReducibilityCertificate(
            status="SPLIT_CONNECTED_SUM",
            start_ch=start_ch,
            terminal_ch=start_ch,
            epd_path=(start_epd,),
            factor_epds=(cut0[0].to_epd(), cut0[1].to_epd()),
            meridian_order_lower_bound=m_ord,
            meridian_witness_group=m_src,
            path_verified=True,
        )

    cap_ch = max_ch if max_ch is not None else (start_ch + 2)
    vis: dict[tuple, list[str]] = {canonical_marked_key(d, up_to_mirror_and_reverse=False): [start_epd]}
    q: deque[tuple[MarkedDiagram, list[str]]] = deque([(d, [start_epd])])

    while q and len(vis) < max_states:
        cur, path = q.popleft()
        allow_stab = cur.num_vertices == start_ch and len(path) <= 2
        for _, nb in list_labeled_single_moves(
            cur,
            max_ch=cap_ch,
            include_stabilizations=allow_stab,
            include_flypes=True,
        ):
            nb_epd = nb.to_epd()
            red, rsteps = _greedy_path_from(nb)
            full_path = path + [nb_epd] + [e for _, e in rsteps]
            if red.num_vertices < start_ch:
                ver = verify_move_path_independent(full_path)
                return RP2ReducibilityCertificate(
                    status="STANDARD_P2" if red.num_vertices <= 2 else "REDUCED_LOWER_CH",
                    start_ch=start_ch,
                    terminal_ch=red.num_vertices,
                    epd_path=tuple(full_path),
                    meridian_order_lower_bound=m_ord,
                    meridian_witness_group=m_src,
                    path_verified=ver.is_valid,
                )
            cut = _find_any_2edge_cut_factors(red)
            if cut is not None:
                ver = verify_move_path_independent(full_path)
                return RP2ReducibilityCertificate(
                    status="SPLIT_CONNECTED_SUM",
                    start_ch=start_ch,
                    terminal_ch=red.num_vertices,
                    epd_path=tuple(full_path),
                    factor_epds=(cut[0].to_epd(), cut[1].to_epd()),
                    meridian_order_lower_bound=m_ord,
                    meridian_witness_group=m_src,
                    path_verified=ver.is_valid,
                )
            for cand_d, cand_p in (
                [(nb, path + [nb_epd]), (red, full_path)] if rsteps else [(nb, path + [nb_epd])]
            ):
                if cand_d.num_vertices <= cap_ch:
                    ck = canonical_marked_key(cand_d, up_to_mirror_and_reverse=False)
                    if ck not in vis:
                        vis[ck] = cand_p
                        q.append((cand_d, cand_p))

    return RP2ReducibilityCertificate(
        status="UNRESOLVED",
        start_ch=start_ch,
        terminal_ch=start_ch,
        epd_path=(start_epd,),
        meridian_order_lower_bound=m_ord,
        meridian_witness_group=m_src,
        path_verified=True,
    )


# ============================================================================
# 4. Faithful 1-Handle Cord Calculus on Bases `8_1`, `9_1`, and `10_1`
#    (Boyle 1988, Trans. AMS 306; Kamada 2014)
# ============================================================================


@dataclass(frozen=True)
class MarkedBigon1HandleRecord:
    """Exact 1-handle decomposition `d = base + 1-handle` across a marked bigon `(u, w)`."""

    bigon_vertices: tuple[int, int]
    base_ch: int
    base_epd: str
    base_identification: str  # "0_1", "8_1", "9_1", "10_1", or "OTHER"
    handle_orientable: bool
    cord_double_coset_key: str
    is_trivial_handle: bool | None


def _extract_marked_bigon_surgeries(
    d: MarkedDiagram,
) -> list[tuple[int, int, MarkedDiagram, int, int]]:
    """Extract all valid marked-bigon 1-handle base surgeries `(u, w, d_base, dart_A, dart_B)`."""
    n = d.num_vertices
    if n < 3 or d.free_loops > 0:
        return []
    out: list[tuple[int, int, MarkedDiagram, int, int]] = []
    for u in range(n):
        if d.vtype[u] not in ("Y", "Z"):
            continue
        for su in range(4):
            su2 = (su + 1) & 3
            da = d.opp[4 * u + su]
            db = d.opp[4 * u + su2]
            w = da >> 2
            if w <= u or (db >> 2) != w or d.vtype[w] not in ("Y", "Z"):
                continue
            sw_a, sw_b = da & 3, db & 3
            if ((sw_a - sw_b) & 3) != 1:
                continue
            if ((d.vtype[u] == d.vtype[w]) != ((su & 1) == (sw_a & 1))):
                continue
            p0 = d.opp[4 * u + ((su + 2) & 3)]
            p1 = d.opp[4 * u + ((su + 3) & 3)]
            p2 = d.opp[4 * w + ((sw_a + 1) & 3)]
            p3 = d.opp[4 * w + ((sw_a + 2) & 3)]
            if any((p >> 2) in (u, w) for p in (p0, p1, p2, p3)):
                continue

            keep_vs = [v for v in range(n) if v not in (u, w)]
            o2n = {v: idx for idx, v in enumerate(keep_vs)}
            b_opp = [0] * (4 * len(keep_vs))
            b_vt = [d.vtype[v] for v in keep_vs]

            def remap_d(od: int) -> int:
                return (o2n[od >> 2] << 2) | (od & 3)

            rp0, rp1, rp2, rp3 = remap_d(p0), remap_d(p1), remap_d(p2), remap_d(p3)
            b_opp[rp0] = rp1
            b_opp[rp1] = rp0
            b_opp[rp2] = rp3
            b_opp[rp3] = rp2
            for v in keep_vs:
                nv = o2n[v]
                for s in range(4):
                    orig_d = (v << 2) | s
                    if orig_d in (p0, p1, p2, p3):
                        continue
                    b_opp[4 * nv + s] = remap_d(d.opp[orig_d])

            d_base = MarkedDiagram(len(keep_vs), tuple(b_opp), tuple(b_vt), 0)
            if (
                certify_classical_unlink(d_base.resolve(1)).is_unlink
                and certify_classical_unlink(d_base.resolve(-1)).is_unlink
            ):
                out.append((u, w, d_base, rp0, rp2))
    return out


def _transport_along_surface(
    d_base: MarkedDiagram,
    eid: dict[int, tuple[int, int]],
    rels: list[tuple],
    asg: list,
    gmul: Callable,
    ginv: Callable,
    one_elem: object,
    start_dart: int,
) -> tuple[dict[int, object], Callable[[int], object]]:
    """Parallel-transport basepoint along `F(d_base)` so `g(dart) = W[dart] * g(start)^{+-1} * W[dart]^{-1}`."""

    def val(x: int) -> object:
        e, s = eid[x]
        v = asg[e]
        return v if s == 1 else ginv(v)

    conj: dict[int, tuple[int, object]] = {}
    for r in rels:
        if r[0] == "conj":
            _, u_out, u_in, o_in, eps = r
            O = val(o_in)
            if eps < 0:
                O = ginv(O)
            conj[u_in] = (u_out, O)
            conj[u_out] = (u_in, ginv(O))

    W: dict[int, object] = {start_dart: one_elem}
    q: deque[int] = deque([start_dart])
    while q:
        x = q.popleft()
        w = W[x]
        v, s = x >> 2, x & 3
        nbrs = [(d_base.opp[x], w)]
        if d_base.vtype[v] in ("Y", "Z"):
            nbrs.extend((4 * v + t, w) for t in range(4) if t != s)
        elif x in conj:
            y, O = conj[x]
            nbrs.append((y, gmul(O, w)))
        else:
            nbrs.append((x ^ 2, w))
        for y, wy in nbrs:
            if y not in W:
                W[y] = wy
                q.append(y)
    return W, val


def _propagate_group_rep(
    k: int,
    eid: dict[int, tuple[int, int]],
    rels: list[tuple],
    seed_gen_a: object,
    seed_gen_b_options: list[object],
    gmul: Callable,
    ginv: Callable,
) -> list | None:
    """Find a non-abelian representation assigning `seed_gen_a` and one of `seed_gen_b_options`."""

    def val(asg: list, x: int):
        e, s = eid[x]
        v = asg[e]
        if v is None:
            return None
        return v if s == 1 else ginv(v)

    def setv(asg: list, x: int, g: object) -> bool:
        e, s = eid[x]
        g2 = g if s == 1 else ginv(g)
        if asg[e] is None:
            asg[e] = g2
            return True
        return asg[e] == g2

    for e1 in range(k):
        for e2 in range(k):
            if e1 == e2:
                continue
            for sb in seed_gen_b_options:
                asg: list = [None] * k
                asg[e1] = seed_gen_a
                asg[e2] = sb
                ok = True
                ch = True
                while ch and ok:
                    ch = False
                    for r in rels:
                        if r[0] == "eq":
                            x, y = val(asg, r[1]), val(asg, r[2])
                            if x is not None and y is not None:
                                if x != ginv(y):
                                    ok = False
                                    break
                            elif x is not None:
                                if not setv(asg, r[2], ginv(x)):
                                    ok = False
                                    break
                                ch = True
                            elif y is not None:
                                if not setv(asg, r[1], ginv(y)):
                                    ok = False
                                    break
                                ch = True
                        else:
                            _, u_out, u_in, o_in, eps = r
                            O = val(asg, o_in)
                            if O is None:
                                continue
                            if eps < 0:
                                O = ginv(O)
                            x, y = val(asg, u_out), val(asg, u_in)
                            if y is not None:
                                t = ginv(gmul(gmul(O, y), ginv(O)))
                                if x is not None:
                                    if x != t:
                                        ok = False
                                        break
                                else:
                                    if not setv(asg, u_out, t):
                                        ok = False
                                        break
                                    ch = True
                            elif x is not None:
                                t = gmul(gmul(ginv(O), ginv(x)), O)
                                if not setv(asg, u_in, t):
                                    ok = False
                                    break
                                ch = True
                if ok and all(v is not None for v in asg):
                    return asg
    return None


def _classify_two_sphere_unlink_1handle(
    d_base: MarkedDiagram,
    dart_a: int,
    dart_b: int,
    k: int,
    eid: dict[int, tuple[int, int]],
    rels: list[tuple],
) -> tuple[str, bool] | None:
    """Classify a 1-handle attached to a 2-component sphere unlink `0_1 u 0_1` via Sanov's faithful
    embedding `F(a, b) -> SL(2, Z)` (`A = [[1, 2], [0, 1]], B = [[1, 0], [2, 1]]`; cf. Kamada
    arXiv:2003.12977 for handles on disconnected surface-links).
    """
    comp, nc = d_base._dart_surface_components()
    if nc != 2:
        return None
    ecomp = {e: comp[x] for x, (e, _) in eid.items()}

    def mm(A, B):
        return (
            (A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]),
            (A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]),
        )

    def mi(A):
        return ((A[1][1], -A[0][1]), (-A[1][0], A[0][0]))

    I2 = ((1, 0), (0, 1))
    GEN = {"a": ((1, 2), (0, 1)), "A": ((1, -2), (0, 1)), "b": ((1, 0), (2, 1)), "B": ((1, 0), (-2, 1))}
    INV = {"a": "A", "A": "a", "b": "B", "B": "b"}

    def red(w: str) -> str:
        out: list[str] = []
        for ch in w:
            if out and out[-1] == INV[ch]:
                out.pop()
            else:
                out.append(ch)
        return "".join(out)

    def winv(w: str) -> str:
        return "".join(INV[c] for c in reversed(w))

    def word_mat(w: str):
        M = I2
        for ch in w:
            M = mm(M, GEN[ch])
        return M

    def decode(M) -> str | None:
        w: list[str] = []

        def norm2(X) -> int:
            return X[0][0] * X[0][0] + X[0][1] * X[0][1] + X[1][0] * X[1][0] + X[1][1] * X[1][1]

        while M != I2:
            best = None
            cur_n = norm2(M)
            for ch in "aAbB":
                X = mm(GEN[INV[ch]], M)
                nx = norm2(X)
                if nx < cur_n and (best is None or nx < best[1]):
                    best = (ch, nx, X)
            if best is None:
                return None
            w.append(best[0])
            M = best[2]
            if len(w) > 300:
                return None
        return "".join(w)

    def val(asg, x: int):
        e, s = eid[x]
        v = asg[e]
        return None if v is None else (v if s == 1 else mi(v))

    def setv(asg, x: int, g) -> bool:
        e, s = eid[x]
        g2 = g if s == 1 else mi(g)
        if asg[e] is None:
            asg[e] = g2
            return True
        return asg[e] == g2

    def propagate(asg) -> bool:
        ch = True
        while ch:
            ch = False
            for r in rels:
                if r[0] == "eq":
                    x, y = val(asg, r[1]), val(asg, r[2])
                    if x is not None and y is not None:
                        if x != mi(y):
                            return False
                    elif x is not None:
                        if not setv(asg, r[2], mi(x)):
                            return False
                        ch = True
                    elif y is not None:
                        if not setv(asg, r[1], mi(y)):
                            return False
                        ch = True
                else:
                    _, out_d, in_d, ov_d, eps = r
                    O = val(asg, ov_d)
                    if O is None:
                        continue
                    if eps < 0:
                        O = mi(O)
                    x, y = val(asg, out_d), val(asg, in_d)
                    if y is not None:
                        t = mi(mm(mm(O, y), mi(O)))
                        if x is not None:
                            if x != t:
                                return False
                        else:
                            if not setv(asg, out_d, t):
                                return False
                            ch = True
                    elif x is not None:
                        t = mm(mm(mi(O), mi(x)), O)
                        if not setv(asg, in_d, t):
                            return False
                        ch = True
        return True

    from itertools import groupby, product

    conjs = [""] + [
        "".join(p)
        for L in (1, 2)
        for p in product("aAbB", repeat=L)
        if red("".join(p)) == "".join(p)
    ]

    def rec_solve(asg, depth: int):
        if all(v is not None for v in asg):
            return asg
        if depth == 0:
            return None
        e = asg.index(None)
        g_ch = "a" if ecomp[e] == 0 else "b"
        for c_w in conjs:
            for gg in (g_ch, INV[g_ch]):
                M = word_mat(red(c_w + gg + winv(c_w)))
                a2 = list(asg)
                a2[e] = M
                if propagate(a2):
                    r_sol = rec_solve(a2, depth - 1)
                    if r_sol is not None:
                        return r_sol
        return None

    e0_list = [e for e in range(k) if ecomp.get(e) == 0]
    e1_list = [e for e in range(k) if ecomp.get(e) == 1]
    if not e0_list or not e1_list:
        return None
    sol = None
    for s1 in ("b", "B"):
        asg0 = [None] * k
        asg0[e0_list[0]] = GEN["a"]
        asg0[e1_list[0]] = GEN[s1]
        if propagate(asg0):
            sol = rec_solve(asg0, 3)
            if sol is not None:
                break
    if sol is None:
        return None

    wa = decode(val(sol, dart_a))
    wb = decode(val(sol, dart_b))
    if wa is None or wb is None:
        return None

    def conj_letter(w: str):
        i, j = 0, len(w) - 1
        while i < j and w[i] == INV[w[j]]:
            i += 1
            j -= 1
        return (w[:i], w[i]) if i == j else None

    ca, cb = conj_letter(wa), conj_letter(wb)
    if ca is None or cb is None:
        return None
    (c_a, gx), (c_b, gy) = ca, cb
    x, y = gx.lower(), gy.lower()
    sgn = (1 if gx.islower() else -1) * (1 if gy.islower() else -1)

    def strip_ends(w: str, left: str, right: str) -> str:
        w = red(w)
        while w and w[0].lower() == left:
            w = w[1:]
        while w and w[-1].lower() == right:
            w = w[:-1]
        return w

    w_raw = strip_ends(winv(c_a) + c_b, x, y)
    vs = set()
    for rev in (False, True):
        ww, xx, yy = (winv(w_raw), y, x) if rev else (w_raw, x, y)
        for fa in (False, True):
            for fb in (False, True):
                for sw in (False, True):
                    t = ww
                    ss = sgn
                    if fa:
                        t = "".join(c.swapcase() if c.lower() == "a" else c for c in t)
                        if (xx == "a") != (yy == "a"):
                            ss = -ss
                    if fb:
                        t = "".join(c.swapcase() if c.lower() == "b" else c for c in t)
                        if (xx == "b") != (yy == "b"):
                            ss = -ss
                    X, Y = xx, yy
                    if sw:
                        t = t.translate(str.maketrans("aAbB", "bBaA"))
                        X = "b" if X == "a" else "a"
                        Y = "b" if Y == "a" else "a"
                    vs.add((X, Y, strip_ends(t, X, Y), ss))
    if x != y:
        vs = {v for v in vs if v[3] == 1 and v[0] == "a"} or vs
    else:
        vs = {v for v in vs if v[0] == "a"} or vs
    X, Y, W, S = min(vs, key=lambda t: (len(t[2]), sum(c.isupper() for c in t[2]), t[2].swapcase()))
    if not W:
        pretty = "1"
    else:
        parts = []
        for ch, grp in groupby(W):
            cnt = len(list(grp))
            exp = cnt if ch.islower() else -cnt
            parts.append(ch.lower() if exp == 1 else f"{ch.lower()}^{exp}")
        pretty = " ".join(parts)
    if X == Y:
        kind = "T" if S == 1 else "K"
        return (f"F2:{kind}({pretty})", W == "")
    return (f"F2:R({pretty})", W == "")


def analyze_marked_bigon_1handles(d: MarkedDiagram) -> list[MarkedBigon1HandleRecord]:
    """Classify all visible marked-bigon 1-handle decompositions of a 1- or 2-component surface diagram `d`
    over bases `0_1`, `0_1 u 0_1` (`F_2` in `SL(2, Z)`), `8_1` (`B_3`), `9_1` (`BS(1, 2)`), and
    `10_1` (`pi_1(4_1)` in `SL(2, Z[omega])`).
    """
    topo_d = d.surface_topology()
    if topo_d.num_components not in (1, 2):
        return []
    handle_ori = all(c.orientable for c in topo_d.components)

    records: list[MarkedBigon1HandleRecord] = []
    for u, w, d_base, dart_a, dart_b in _extract_marked_bigon_surgeries(d):
        st_b = d_base.surface_topology()
        if st_b.num_components == 2 and all(c.homeomorphism_type == "S^2" for c in st_b.components):
            d_red = greedy_yoshikawa_reduce(d_base)
            k, eid, rels = _indep_wirtinger_relations(
                d_base.num_vertices, d_base.opp, d_base.vtype, d_base.free_loops
            )
            f2_res = _classify_two_sphere_unlink_1handle(d_base, dart_a, dart_b, k, eid, rels)
            if f2_res is not None:
                key_s, triv_s = f2_res
                records.append(
                    MarkedBigon1HandleRecord(
                        bigon_vertices=(u, w),
                        base_ch=d_red.num_vertices,
                        base_epd=d_red.to_epd() if d_red.num_vertices > 0 else "O[1], O[1]",
                        base_identification="0_1 u 0_1",
                        handle_orientable=handle_ori,
                        cord_double_coset_key=key_s,
                        is_trivial_handle=triv_s,
                    )
                )
                continue
        if st_b.num_components != 1 or st_b.components[0].homeomorphism_type != "S^2":
            continue
        d_red = greedy_yoshikawa_reduce(d_base)
        if d_red.num_vertices == 0:
            records.append(
                MarkedBigon1HandleRecord(
                    bigon_vertices=(u, w),
                    base_ch=0,
                    base_epd="O[1]",
                    base_identification="0_1",
                    handle_orientable=handle_ori,
                    cord_double_coset_key="Z:trivial",
                    is_trivial_handle=True,
                )
            )
            continue

        k, eid, rels = _indep_wirtinger_relations(
            d_base.num_vertices, d_base.opp, d_base.vtype, d_base.free_loops
        )

        # Case 1: Try faithful B_3 -> SL(2, Z) x Z (`8_1` = spun trefoil):
        def b3_mm(A, B):
            return (
                (A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]),
                (A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]),
            )

        def b3_mi(A):
            return ((A[1][1], -A[0][1]), (-A[1][0], A[0][0]))

        def b3_mul(x, y):
            return (b3_mm(x[0], y[0]), x[1] + y[1])

        def b3_inv(x):
            return (b3_mi(x[0]), -x[1])

        b3_one = (((1, 0), (0, 1)), 0)
        b3_a = (((1, 1), (0, 1)), 1)
        b3_b = (((1, 0), (-1, 1)), 1)
        asg_b3 = _propagate_group_rep(k, eid, rels, b3_a, [b3_b, b3_inv(b3_b)], b3_mul, b3_inv)
        if asg_b3 is not None:
            W_map, val_fn = _transport_along_surface(
                d_base, eid, rels, asg_b3, b3_mul, b3_inv, b3_one, dart_a
            )
            mu = val_fn(dart_a)
            Wb = W_map[dart_b]
            # Normalize mu to b3_a via SL(2, Z) conjugation:
            M, e_mu = mu
            psi = e_mu == -1
            if psi:
                M = ((M[0][0], -M[0][1]), (-M[1][0], M[1][1]))
            vx, vy = M[0][1], 1 - M[0][0]
            if vx == 0 and vy == 0:
                vx, vy = 1 - M[1][1], M[1][0]
            from math import gcd

            g_v = gcd(vx, vy)
            vx //= g_v
            vy //= g_v

            def _egcd(a_: int, b_: int) -> tuple[int, int]:
                if b_ == 0:
                    return (1, 0) if a_ > 0 else (-1, 0)
                u_, v_ = _egcd(b_, a_ % b_)
                return (v_, u_ - (a_ // b_) * v_)

            s_c, t_c = _egcd(vx, vy)
            if s_c * vx + t_c * vy != 1:
                s_c, t_c = -s_c, -t_c
            C = ((s_c, t_c), (-vy, vx))
            R = b3_mm(b3_mm(C, M), b3_mi(C))
            if R == ((1, -1), (0, 1)):
                C = ((-s_c, -t_c), (vy, -vx))
            Wm, ew = Wb
            if psi:
                Wm = ((Wm[0][0], -Wm[0][1]), (-Wm[1][0], Wm[1][1]))
                ew = -ew
            Wn = b3_mm(b3_mm(C, Wm), b3_mi(C))

            def _b3_coset_key(Mat, ex: int):
                (p, q), (r, s) = Mat
                if r == 0:
                    sgn = 1 if p == 1 else -1
                    return ("r=0", sgn, ex - q * sgn)
                m_abs = abs(r)
                return ("r!=0", r, p % m_abs, s % m_abs, ex * r - p - s)

            ks = []
            for X, ex in ((Wn, ew), (b3_mi(Wn), -ew)):
                ks.append(_b3_coset_key(X, ex))
                Xp = ((X[0][0], -X[0][1]), (-X[1][0], X[1][1]))
                ks.append(_b3_coset_key(Xp, -ex))
            best_k = min(ks, key=repr)
            triv = Wn == ((1, ew), (0, 1))
            records.append(
                MarkedBigon1HandleRecord(
                    bigon_vertices=(u, w),
                    base_ch=d_red.num_vertices,
                    base_epd=d_red.to_epd(),
                    base_identification="8_1",
                    handle_orientable=handle_ori,
                    cord_double_coset_key=f"B3:{best_k}",
                    is_trivial_handle=triv,
                )
            )
            continue

        # Case 2: Try faithful BS(1, 2) -> Aff(Z[1/2]) (`9_1`):
        def bs_mul(g, h):
            return (g[0] + h[0], g[1] + (Fraction(2) ** g[0]) * h[1])

        def bs_inv(g):
            return (-g[0], -g[1] / (Fraction(2) ** g[0]))

        bs_one = (0, Fraction(0))
        bs_x = (1, Fraction(0))
        bs_opts = [(1, Fraction(1)), (1, Fraction(-1)), bs_inv((1, Fraction(1))), bs_inv((1, Fraction(-1)))]
        asg_bs = _propagate_group_rep(k, eid, rels, bs_x, bs_opts, bs_mul, bs_inv)
        if asg_bs is not None:
            W_map, val_fn = _transport_along_surface(
                d_base, eid, rels, asg_bs, bs_mul, bs_inv, bs_one, dart_a
            )
            mu = val_fn(dart_a)
            Wb = W_map[dart_b]
            m_norm = mu if mu[0] == 1 else bs_inv(mu)
            T = (0, m_norm[1])
            Wn = bs_mul(bs_mul(T, Wb), bs_inv(T))
            num_odd = abs(Wn[1].numerator)
            while num_odd > 0 and (num_odd & 1) == 0:
                num_odd >>= 1
            records.append(
                MarkedBigon1HandleRecord(
                    bigon_vertices=(u, w),
                    base_ch=d_red.num_vertices,
                    base_epd=d_red.to_epd(),
                    base_identification="9_1",
                    handle_orientable=handle_ori,
                    cord_double_coset_key=f"BS(1,2):odd_part={num_odd}",
                    is_trivial_handle=(num_odd == 0),
                )
            )
            continue

        # Case 3: Try faithful Riley representation pi_1(4_1) -> SL(2, Z[omega]) x Z (`10_1`):
        def zw_add(x, y):
            return (x[0] + y[0], x[1] + y[1])

        def zw_neg(x):
            return (-x[0], -x[1])

        def zw_mul(x, y):
            a0, b0 = x
            c0, d0 = y
            return (a0 * c0 - b0 * d0, a0 * d0 + b0 * c0 - b0 * d0)

        def zw_norm(x) -> int:
            # N(a + b*omega) = a^2 - a*b + b^2
            return x[0] * x[0] - x[0] * x[1] + x[1] * x[1]

        def f8_mm(A, B):
            return (
                (
                    zw_add(zw_mul(A[0][0], B[0][0]), zw_mul(A[0][1], B[1][0])),
                    zw_add(zw_mul(A[0][0], B[0][1]), zw_mul(A[0][1], B[1][1])),
                ),
                (
                    zw_add(zw_mul(A[1][0], B[0][0]), zw_mul(A[1][1], B[1][0])),
                    zw_add(zw_mul(A[1][0], B[0][1]), zw_mul(A[1][1], B[1][1])),
                ),
            )

        def f8_mi(A):
            return ((A[1][1], zw_neg(A[0][1])), (zw_neg(A[1][0]), A[0][0]))

        def f8_mul(x, y):
            return (f8_mm(x[0], y[0]), x[1] + y[1])

        def f8_inv(x):
            return (f8_mi(x[0]), -x[1])

        z1 = (1, 0)
        z0 = (0, 0)
        f8_one = (((z1, z0), (z0, z1)), 0)
        f8_a = (((z1, z1), (z0, z1)), 1)
        f8_b1 = (((z1, z0), ((0, -1), z1)), 1)
        f8_b2 = (((z1, z0), ((1, 1), z1)), 1)
        asg_f8 = _propagate_group_rep(
            k,
            eid,
            rels,
            f8_a,
            [f8_b1, f8_inv(f8_b1), f8_b2, f8_inv(f8_b2)],
            f8_mul,
            f8_inv,
        )
        if asg_f8 is not None:
            e_ref = next(idx for idx, val_g in enumerate(asg_f8) if val_g == f8_a)
            dart_ref = next(hd for hd, (eg, sgn) in eid.items() if eg == e_ref and sgn == 1)
            W_ref, val_fn = _transport_along_surface(
                d_base, eid, rels, asg_f8, f8_mul, f8_inv, f8_one, dart_ref
            )
            if dart_a in W_ref and dart_b in W_ref:
                Wa = W_ref[dart_a]
                Wb_raw = W_ref[dart_b]
                mu_a_norm = f8_mul(f8_inv(Wa), f8_mul(val_fn(dart_a), Wa))
                W_rel = f8_mul(f8_inv(Wa), Wb_raw)
                # Conjugating by `Wa^-1` normalizes the meridian at `dart_a` to `mu_a_norm == f8_a^{+-1}`;
                # `<f8_a>\G/<f8_a>` shifts upper-right by Z, so `N(r)` (lower-left norm in Z[omega]) is an
                # exact `<mu>`-double-coset invariant with `r == (0, 0)` and `Wn == f8_a^{ew}` <=> trivial handle!
                Wn, ew = W_rel
                r_entry = Wn[1][0]
                triv = r_entry == z0 and Wn[0][0] == z1 and Wn[0][1] == (ew * mu_a_norm[1], 0)
                records.append(
                    MarkedBigon1HandleRecord(
                        bigon_vertices=(u, w),
                        base_ch=d_red.num_vertices,
                        base_epd=d_red.to_epd(),
                        base_identification="10_1",
                        handle_orientable=handle_ori,
                        cord_double_coset_key=f"4_1:norm_r={zw_norm(r_entry)}:triv={triv}",
                        is_trivial_handle=triv,
                    )
                )
                continue

        records.append(
            MarkedBigon1HandleRecord(
                bigon_vertices=(u, w),
                base_ch=d_red.num_vertices,
                base_epd=d_red.to_epd(),
                base_identification="OTHER",
                handle_orientable=handle_ori,
                cord_double_coset_key="unclassified_base",
                is_trivial_handle=None,
            )
        )
    return records


# ============================================================================
# 5. Campaign 3: Constructive Move-Search Harness (Jablonowski, Viro, Klein Tube)
# ============================================================================


def are_prime_diagrams_tunnel_merged_c(
    ch: int,
    d1: MarkedDiagram,
    d2: MarkedDiagram,
    *,
    depth_plus1: int = 16,
    depth_plus2: int = 8,
    max_bfs_states: int = 512,
    num_workers: int = 1,
) -> bool:
    """Safely check whether two `ch`-vertex diagrams `d1` and `d2` merge into a single prime
    tunnel component under `tunnel_diagrams_c` without conflating component discard (`< ch` or
    planar 2-edge connected-sum cut) with equivalence!
    """
    from surfaceknots.fast_enum import tunnel_diagrams_c

    if (
        ch < 3
        or d1.num_vertices != ch
        or d2.num_vertices != ch
        or d1.free_loops != 0
        or d2.free_loops != 0
    ):
        return False
    kw = {
        "depth_plus1": depth_plus1,
        "depth_plus2": depth_plus2,
        "max_bfs_states": max_bfs_states,
        "num_workers": num_workers,
    }
    if len(tunnel_diagrams_c(ch, [d1], **kw)) != 1:
        return False
    if len(tunnel_diagrams_c(ch, [d2], **kw)) != 1:
        return False
    return len(tunnel_diagrams_c(ch, [d1, d2], **kw)) == 1


def _all_admissible_planar_sums(d1: MarkedDiagram, d2: MarkedDiagram, max_reps: int = 8) -> list[str]:
    """Collect distinct admissible planar connected-sum diagrams of `d1` and `d2` across dart pairs."""
    from surfaceknots.ch11_census import planar_connected_sum

    seen: set[tuple] = set()
    out: list[str] = []
    step1 = max(1, (4 * d1.num_vertices) // 8)
    step2 = max(1, (4 * d2.num_vertices) // 8)
    for h1 in range(0, 4 * d1.num_vertices, step1):
        for h2 in range(0, 4 * d2.num_vertices, step2):
            for fl in (False, True):
                cs = planar_connected_sum(d1, h1, d2, h2, flip_splice=fl)
                if (
                    certify_classical_unlink(cs.resolve(1)).is_unlink
                    and certify_classical_unlink(cs.resolve(-1)).is_unlink
                ):
                    ck = canonical_marked_key(cs, up_to_mirror_and_reverse=True)
                    if ck not in seen:
                        seen.add(ck)
                        out.append(cs.to_epd())
                        if len(out) >= max_reps:
                            return out
    return out


def export_move_search_results(
    data_dir: "Path",
    num_workers: int = 32,
) -> dict[str, object]:
    """Run Campaign 3 constructive bidirectional Yoshikawa move searches on:
    1. Jablonowski (2015, Remark 3.15) open problem: `10_2 # P_+` (`ch = 12`) vs `8_1 # P_+` (`ch = 10`).
    2. Viro (1973) smooth collapse: `10_3 # P_+` (`ch = 12`) -> `2_1^-1 (P_+)` (`ch = 2`).
    3. Non-orientable longitude 1-handle on spun trefoil:
       `8_1 + non-orientable h(Delta^2)` (`ch = 11`) vs `8_1 # K^2(e=0)` (`ch = 11`).
    Writes `code/data/complexity/move_search_results.json`.
    """
    import json
    from pathlib import Path
    import time
    from surfaceknots.complexity import MINIMAL_KLEIN_BOTTLE_E0_EPD
    from surfaceknots.families import CLASSICAL_PLAT_BY_LABEL, compile_mixed_plat_link, simplify_compiled_diagram
    from surfaceknots.fast_enum import tunnel_diagrams_c
    from surfaceknots.yoshikawa_data import YOSHIKAWA_BY_LABEL

    data_dir = Path(data_dir)
    out_dir = data_dir / "complexity"
    out_dir.mkdir(parents=True, exist_ok=True)

    t_start = time.perf_counter()
    rp_raw = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^-1"].canonical_epd)
    p_pos = rp_raw if rp_raw.surface_topology().components[0].normal_euler == 2 else rp_raw.mirror_3d()
    k2_0 = MarkedDiagram.from_epd(MINIMAL_KLEIN_BOTTLE_E0_EPD)
    d_8_1 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["8_1"].canonical_epd)
    d_10_2 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_2"].canonical_epd)
    d_10_3 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_3"].canonical_epd)

    def _serialize_res(res: BidirectionalSearchResult) -> dict[str, object]:
        ver = res.independent_verification
        return {
            "connected": res.connected,
            "max_ch": res.max_ch,
            "waves_run": res.waves_run,
            "forward_states": res.forward_states,
            "backward_states": res.backward_states,
            "num_steps": len(res.move_tags),
            "move_tags": list(res.move_tags),
            "epd_path": list(res.epd_path),
            "independent_verification": (
                {
                    "is_valid": ver.is_valid,
                    "num_steps": ver.num_steps,
                    "step_kinds": list(ver.step_kinds),
                    "error": ver.error,
                }
                if ver is not None
                else None
            ),
        }

    # 1. Jablonowski (2015, Remark 3.15): `10_2 # P_+` (`ch = 12`) vs `8_1 # P_+` (`ch = 10`)
    t1 = time.perf_counter()
    starts_jab = _all_admissible_planar_sums(d_10_2, p_pos, max_reps=4) + _all_admissible_planar_sums(
        d_10_2.mirror_3d(), p_pos, max_reps=4
    )
    goals_jab = _all_admissible_planar_sums(d_8_1, p_pos, max_reps=6)
    res_jab = parallel_bidirectional_move_search(
        starts_jab,
        goals_jab,
        max_ch=16,
        max_waves=2,
        beam_width=200,
        workers=num_workers,
        allow_o4=True,
    )
    dt_jab = round(time.perf_counter() - t1, 3)

    # 2. Constructive Viro (1973) collapse: `10_3 # P_+` (`ch = 12`) -> `P_+` (`ch = 2`)
    t2 = time.perf_counter()
    starts_viro = _all_admissible_planar_sums(d_10_3, p_pos, max_reps=4) + _all_admissible_planar_sums(
        d_10_3.mirror_3d(), p_pos, max_reps=4
    )
    goals_viro = [p_pos.to_epd(), p_pos.mirror_3d().to_epd()]
    res_viro = parallel_bidirectional_move_search(
        starts_viro,
        goals_viro,
        max_ch=15,
        max_waves=2,
        beam_width=200,
        workers=num_workers,
        allow_o4=True,
    )
    dt_viro = round(time.perf_counter() - t2, 3)

    # 3. `8_1 + non-orientable h(Delta^2 / longitude)` (`ch = 11`) vs `8_1 # K^2(e=0)` (`ch = 11`)
    t3 = time.perf_counter()
    klein_tube_3_1 = simplify_compiled_diagram(
        compile_mixed_plat_link(CLASSICAL_PLAT_BY_LABEL["3_1"], {0: 0, 1: 1}),
        max_bfs_states=32,
    )
    starts_kt = [klein_tube_3_1.to_epd()]
    # Also include any ch=11 residual diagrams with cord key B3:('r=0', -1, -6) if expanded_ch11.json exists:
    ch11_path = data_dir / "expanded_ch11.json"
    if ch11_path.exists():
        c11_data = json.loads(ch11_path.read_text())
        for rc in c11_data.get("residual_collisions", []):
            d_rc = MarkedDiagram.from_epd(str(rc["epd"]))
            for hb in analyze_marked_bigon_1handles(d_rc):
                if (
                    hb.base_identification == "8_1"
                    and not hb.handle_orientable
                    and "'r=0', -1, -6" in hb.cord_double_coset_key
                ):
                    if d_rc.to_epd() not in starts_kt:
                        starts_kt.append(d_rc.to_epd())

    goals_kt = _all_admissible_planar_sums(d_8_1, k2_0, max_reps=6)
    res_kt = parallel_bidirectional_move_search(
        starts_kt,
        goals_kt,
        max_ch=15,
        max_waves=2,
        beam_width=200,
        workers=num_workers,
        allow_o4=True,
    )
    # Also cross-check with C OpenMP mountain-pass tunneling at ch = 11 (+2 -> 13):
    c_tunnel_merged = False
    if starts_kt and goals_kt:
        c_tunnel_merged = are_prime_diagrams_tunnel_merged_c(
            11,
            MarkedDiagram.from_epd(starts_kt[0]),
            MarkedDiagram.from_epd(goals_kt[0]),
            depth_plus1=384,
            depth_plus2=256,
            max_bfs_states=8192,
            num_workers=num_workers,
        )
    dt_kt = round(time.perf_counter() - t3, 3)

    # 4. `11_1.{a,b,c}` equivalence packet mountain-pass check up to `max_ch = 15` (HANDOFF §4 Campaign 3):
    t4 = time.perf_counter()
    pkt_results: dict[str, object] = {}
    if ch11_path.exists():
        c11_data = json.loads(ch11_path.read_text())
        ent_11_1 = next((e for e in c11_data.get("entries", []) if e.get("label") == "11_1"), None)
        if ent_11_1 is not None:
            comp_epds = list(ent_11_1.get("move_component_epds", []))
            if len(comp_epds) >= 2:
                res_ab = parallel_bidirectional_move_search(
                    [comp_epds[0]],
                    comp_epds[1:],
                    max_ch=15,
                    max_waves=2,
                    beam_width=200,
                    workers=num_workers,
                    allow_o4=True,
                )
                pkt_results = {
                    "problem": "HANDOFF §4 Campaign 3: Does 11_1.a connect to 11_1.{b,c} via Yoshikawa moves up to ch=15?",
                    "start_epds_count": 1,
                    "goal_epds_count": len(comp_epds) - 1,
                    "runtime_sec": round(time.perf_counter() - t4, 3),
                    **_serialize_res(res_ab),
                }

    targets_dict: dict[str, object] = {
        "jablonowski_10_2_plus_P_vs_8_1_plus_P": {
            "problem": "Jablonowski (2015, Remark 3.15): Is tau^2(3_1) # P_+ smoothly/diagrammatically equivalent to tau^0(3_1) # P_+?",
            "start_epds_count": len(starts_jab),
            "goal_epds_count": len(goals_jab),
            "runtime_sec": dt_jab,
            **_serialize_res(res_jab),
        },
        "viro_10_3_plus_P_collapse_to_P": {
            "problem": "Viro (1973): Constructive Yoshikawa move path for tau^3(3_1) # P_+ (ch=12) -> P_+ (ch=2)",
            "topological_equivalence_known": True,
            "citation": "Viro (1973); Yoshikawa (1994, Remark 2.8)",
            "start_epds_count": len(starts_viro),
            "goal_epds_count": len(goals_viro),
            "runtime_sec": dt_viro,
            **_serialize_res(res_viro),
        },
        "trefoil_nonorientable_longitude_1handle_vs_8_1_sum_K2": {
            "problem": "AUDIT §9.3: Is 8_1 + non-orientable h(Delta^2) (Klein analogue of Tube(3_1), ch=11) equivalent to 8_1 # K^2(e=0) (ch=11)?",
            "klein_tube_compiled_epd": klein_tube_3_1.to_epd(),
            "klein_tube_compiled_ch": klein_tube_3_1.num_vertices,
            "c_omp_tunnel_ch13_merged": c_tunnel_merged,
            "start_epds_count": len(starts_kt),
            "goal_epds_count": len(goals_kt),
            "runtime_sec": dt_kt,
            **_serialize_res(res_kt),
        },
    }
    if pkt_results:
        targets_dict["packet_11_1_abc_mountain_pass"] = pkt_results

    payload: dict[str, object] = {
        "total_wall_sec": round(time.perf_counter() - t_start, 3),
        "num_workers": num_workers,
        "targets": targets_dict,
    }
    out_path = out_dir / "move_search_results.json"
    out_path.write_text(json.dumps(payload, indent=2) + "\n")
    return payload
