"""Surface-link sublink extraction, peripheral SL(2,F_3) longitudes, and mod-2 linking invariants.

Implements:
1. `extract_subsurface_diagram` and `sublink_summary`: exact diagram projection onto any nonempty
   subset of connected components `F_{i_1} u ... u F_{i_k}`, certifying proper sublink triviality,
   splittability, and Brunnian surface-link structure (e.g. `10_1^{0,0,1}`).
2. `peripheral_sl2f3_profile`: evaluates the image of the peripheral subgroup
   `P_i = <iota_*(pi_1(F_i)), mu_i> <= pi_1(S^4 \\ F)` under all homomorphisms
   `rho: pi_1(S^4 \\ F) -> SL(2, F_3)`, separating ribbon-torus knots such as `10_1^1 = Tube(3_1)`
   (`24` representations with `|rho(P_0)| = 6 > 3 = |<rho(mu_0)>|`) from `8_1` and `8_1 # 2_1^1` (`0`).
3. `mod2_triple_linking_profile`: Carter-Saito-Satoh mod-2 triple point linking invariants
   `t_{ijk}(F) in Z_2` (`|{i, j, k}| = 3`) extracted from explicit Reidemeister unlink movies `(M_-, M_+)`.
4. `double_linking_mod2`: Rohlin-Matsumoto / Carter-Saito-Sato mod-2 torus-linking invariant
   `q_{F_i}(PD(lambda_{ij})) in {0, 1}` evaluating the pin^- membrane quadratic form on the unique
   nonzero class `PD(lambda_{ij}) in ker(lambda_{ij}) subset H_1(F_i; Z_2)` across mixed saddle slices,
   separating `8_1^{1,1}` (`(0, 0)`) from the 1-twist-spun Hopf link `10_1^{1,1}` (`(1, 1)`).
5. `exact_sublink_invariant_profile`: presentation- and BFS-independent `O(4)`-invariant profile of
   all 1-component (and, for `c >= 3`, 2-component) proper sublinks, separating links with unknotted
   components (such as `11_2^{0,1}` and `11_1^{0,-2}`) from connected sums sharing the same ambient
   link group and `E_1` ideal (such as `6_1^{0,1} # 8_1` and `7_1^{0,-2} # 8_1`).
"""

from __future__ import annotations

from itertools import combinations, permutations, product
import math
import sympy as sp

from surfaceknots.diagrams import ClassicalLinkDiagram, MarkedDiagram
from surfaceknots.invariants.fundamental_group import cyclically_reduce_word
from surfaceknots.moves import connected_diagram_components, simplify_and_explore_yoshikawa
from surfaceknots.unlink import UnlinkMovieStep, _is_over, find_unlink_movie


def extract_subsurface_diagram(d: MarkedDiagram, keep_components: tuple[int, ...]) -> MarkedDiagram:
    """Return the `MarkedDiagram` of the sub-link consisting of components in `keep_components`.

    Vertices of omitted components are removed, and mixed classical crossings between a kept
    component and an omitted component are spliced straight along the kept strand.
    """
    keep_set = set(keep_components)
    n = d.num_vertices
    dart_comp, total_comps = d._dart_surface_components()
    num_vert_comps = (max(dart_comp) + 1) if n > 0 else 0

    kept_free = sum(1 for c in range(num_vert_comps, total_comps) if c in keep_set)
    if n == 0:
        return MarkedDiagram(num_vertices=0, opp=(), vtype=(), free_loops=kept_free)

    surviving_vs: list[int] = []
    old_to_new: dict[int, int] = {}
    internal_splice: dict[int, int] = {}

    for v in range(n):
        vt = d.vtype[v]
        if vt in ("Y", "Z"):
            if dart_comp[4 * v] in keep_set:
                old_to_new[v] = len(surviving_vs)
                surviving_vs.append(v)
        else:
            c02 = dart_comp[4 * v] in keep_set
            c13 = dart_comp[4 * v + 1] in keep_set
            if c02 and c13:
                old_to_new[v] = len(surviving_vs)
                surviving_vs.append(v)
            elif c02:
                internal_splice[4 * v + 0] = 4 * v + 2
                internal_splice[4 * v + 2] = 4 * v + 0
            elif c13:
                internal_splice[4 * v + 1] = 4 * v + 3
                internal_splice[4 * v + 3] = 4 * v + 1

    m = len(surviving_vs)
    new_opp = [0] * (4 * m)
    new_vtype = tuple(d.vtype[v] for v in surviving_vs)
    visited_spliced: set[int] = set()

    for v_old in surviving_vs:
        v_new = old_to_new[v_old]
        for s in range(4):
            cur = d.opp[4 * v_old + s]
            while cur in internal_splice:
                visited_spliced.add(cur)
                nxt = internal_splice[cur]
                visited_spliced.add(nxt)
                cur = d.opp[nxt]
            new_opp[4 * v_new + s] = (old_to_new[cur >> 2] << 2) | (cur & 3)

    extra_loops = 0
    for sd in internal_splice:
        if sd in visited_spliced:
            continue
        extra_loops += 1
        cur = sd
        while cur not in visited_spliced:
            visited_spliced.add(cur)
            nxt = internal_splice[cur]
            visited_spliced.add(nxt)
            cur = d.opp[nxt]

    return MarkedDiagram(
        num_vertices=m,
        opp=tuple(new_opp),
        vtype=new_vtype,
        free_loops=kept_free + extra_loops,
    )


def _is_certified_standard_unknot_single(d: MarkedDiagram) -> bool:
    """Check whether a simplified 1-component diagram `d` is a standard minimal unknotted surface."""
    topo = d.surface_topology()
    if topo.num_components != 1:
        return False
    c = topo.components[0]
    if c.orientable:
        # Standard unknotted Sigma_g has ch = 2*g and 0 classical crossings!
        return d.num_classical == 0 and d.num_vertices == 2 * c.genus_symbol
    # Standard unknotted #_k RP^2 with normal Euler e has minimal ch = k + max(1, |e|/2)
    # (e.g. ch = 2 for RP^2, ch = 3 for K^2(e=0), ch = 4 for K^2(e=+-4)) and cyclic group!
    k_non = abs(c.genus_symbol)
    if d.num_vertices <= max(2 * k_non, k_non + max(1, abs(c.normal_euler) // 2)):
        from surfaceknots.invariants.fundamental_group import wirtinger_presentation

        wp = wirtinger_presentation(d)
        return wp.num_generators <= 1
    return False


def sublink_summary(d: MarkedDiagram, max_states: int = 256) -> dict[str, object]:
    """Summarize all proper nonempty sub-surface-links of `d`."""
    topo = d.surface_topology()
    c = topo.num_components
    if c <= 1:
        return {
            "num_components": c,
            "all_singletons_unknotted": None,
            "is_brunnian": False,
            "proper_sublinks": [],
        }

    singletons_trivial = True
    all_proper_trivial_or_split = True
    records: list[dict[str, object]] = []

    for r_sz in range(1, c):
        for comb in combinations(range(c), r_sz):
            sub_d = extract_subsurface_diagram(d, comb)
            best, is_split, is_comp = simplify_and_explore_yoshikawa(sub_d, max_states=max_states)
            comp_types = tuple(topo.components[i].homeomorphism_type for i in comb)
            if r_sz == 1:
                triv = _is_certified_standard_unknot_single(best)
                if not triv:
                    singletons_trivial = False
                    all_proper_trivial_or_split = False
                records.append(
                    {
                        "components": comb,
                        "types": comp_types,
                        "min_ch_found": best.num_vertices,
                        "certified_unknot": triv,
                    }
                )
            else:
                # For r_sz >= 2, check if `best` decomposes into `r_sz` certified unknotted singletons!
                pieces = connected_diagram_components(best)
                completely_trivial = bool(pieces) and all(
                    p.surface_topology().num_components == 1
                    and _is_certified_standard_unknot_single(p)
                    for p in pieces
                )
                if not completely_trivial:
                    all_proper_trivial_or_split = False
                records.append(
                    {
                        "components": comb,
                        "types": comp_types,
                        "min_ch_found": best.num_vertices,
                        "certified_split": is_split or (best.num_vertices == 0),
                        "certified_composite": is_comp,
                        "certified_trivial_unlink": completely_trivial,
                    }
                )

    return {
        "num_components": c,
        "all_singletons_unknotted": singletons_trivial,
        "is_brunnian": (c >= 3 and all_proper_trivial_or_split),
        "proper_sublinks": records,
    }


# ---------------------------------------------------------------------------
# Precomputed multiplication table for SL(2, F_3) (order 24)
# ---------------------------------------------------------------------------


def _build_sl2_f3() -> tuple[int, list[list[int]], list[int]]:
    elems: list[tuple[int, int, int, int]] = []
    for a in range(3):
        for b in range(3):
            for c in range(3):
                for d_el in range(3):
                    if (a * d_el - b * c) % 3 == 1:
                        elems.append((a, b, c, d_el))
    idx = {m: i for i, m in enumerate(elems)}
    n = len(elems)
    mul = [[0] * n for _ in range(n)]
    inv = [0] * n
    for i, (a, b, c, d_el) in enumerate(elems):
        inv[i] = idx[(d_el % 3, (-b) % 3, (-c) % 3, a % 3)]
        for j, (e, f, g, h) in enumerate(elems):
            mul[i][j] = idx[
                (
                    (a * e + b * g) % 3,
                    (a * f + b * h) % 3,
                    (c * e + d_el * g) % 3,
                    (c * f + d_el * h) % 3,
                )
            ]
    return idx[(1, 0, 0, 1)], mul, inv


_SL2_ID, _SL2_MUL, _SL2_INV = _build_sl2_f3()


def _sl2_subgroup_size(gens: list[int]) -> int:
    seen = {_SL2_ID}
    st = [_SL2_ID]
    g_all = list(set(gens) | {_SL2_INV[g] for g in gens})
    while st:
        x = st.pop()
        for g in g_all:
            y = _SL2_MUL[x][g]
            if y not in seen:
                seen.add(y)
                st.append(y)
    return len(seen)


def peripheral_sl2f3_profile(
    d: MarkedDiagram,
) -> tuple[int, tuple[tuple[tuple[str, int, int], ...], ...]]:
    """Compute `(num_homs, non_meridional_peripheral_profiles)` into `SL(2, F_3)`.

    For each homomorphism `rho: pi_1(S^4 \\ F) -> SL(2, F_3)` and each orientable component `F_i`,
    compares `|<rho(mu_i)>|` against `|<rho(iota_*(pi_1(F_i))), rho(mu_i)>|`.
    Records `(homeomorphism_type, ord_meridian, ord_peripheral)` whenever `ord_peripheral > ord_meridian`.
    """
    n = d.num_vertices
    topo = d.surface_topology()
    dc, num_scomps = d._dart_surface_components()
    if n == 0:
        return (24 ** d.free_loops, ())

    is_cl = [vt in ("X0", "X1") for vt in d.vtype]

    def is_und(v: int, s: int) -> bool:
        return ((s & 1) == 0) if d.vtype[v] == "X0" else ((s & 1) == 1)

    def mpart(v: int, s: int) -> int:
        # L_- resolution (`from_resolution = -1`)
        return (3 - s) if d.vtype[v] == "Z" else (s ^ 1)

    dart_dir = [0] * (4 * n)
    cycles: list[list[int]] = []
    for sd in range(4 * n):
        if dart_dir[sd] != 0:
            continue
        cyc: list[int] = []
        cur = sd
        while dart_dir[cur] == 0:
            dart_dir[cur] = 1
            cyc.append(cur)
            v, s = cur >> 2, cur & 3
            exs = (s ^ 2) if is_cl[v] else mpart(v, s)
            out_d = (v << 2) | exs
            dart_dir[out_d] = -1
            cur = d.opp[out_d]
        cycles.append(cyc)

    dart_arc = [-1] * (4 * n)
    g = 0
    for cyc in cycles:
        upos = [k for k, hd in enumerate(cyc) if is_cl[hd >> 2] and is_und(hd >> 2, hd & 3)]
        if not upos:
            for hd in cyc:
                v, s = hd >> 2, hd & 3
                exs = (s ^ 2) if is_cl[v] else mpart(v, s)
                dart_arc[hd] = g
                dart_arc[(v << 2) | exs] = g
            g += 1
        else:
            m = len(upos)
            for idx_u in range(m):
                p0 = upos[idx_u]
                p1 = upos[(idx_u + 1) % m]
                v0, s0 = cyc[p0] >> 2, cyc[p0] & 3
                dart_arc[(v0 << 2) | (s0 ^ 2)] = g
                k = (p0 + 1) % len(cyc)
                while k != p1:
                    hd = cyc[k]
                    v, s = hd >> 2, hd & 3
                    exs = (s ^ 2) if is_cl[v] else mpart(v, s)
                    dart_arc[hd] = g
                    dart_arc[(v << 2) | exs] = g
                    k = (k + 1) % len(cyc)
                dart_arc[cyc[p1]] = g
                g += 1

    # Factor `d.free_loops` analytically (`24 ** free_loops`) instead of branching 24x per free loop!
    free_mult = 24 ** d.free_loops
    if g == 0:
        return (free_mult, ())

    rels: list[tuple[int, ...]] = []
    under_step_lit: dict[tuple[int, int], int] = {}

    for v in range(n):
        vt = d.vtype[v]
        if vt in ("X0", "X1"):
            us = (0, 2) if vt == "X0" else (1, 3)
            os = (1, 3) if vt == "X0" else (0, 2)
            u_in = us[0] if dart_dir[(v << 2) | us[0]] == 1 else us[1]
            u_out = u_in ^ 2
            o_in = os[0] if dart_dir[(v << 2) | os[0]] == 1 else os[1]
            x = dart_arc[(v << 2) | u_in] + 1
            y = dart_arc[(v << 2) | u_out] + 1
            z = dart_arc[(v << 2) | o_in] + 1
            diff = (o_in - u_in) & 3
            if diff == 3:
                w_rel = cyclically_reduce_word((z, x, -z, -y))
                w_lit = -z
            else:
                w_rel = cyclically_reduce_word((-z, x, z, -y))
                w_lit = z
            if w_rel:
                rels.append(w_rel)
            under_step_lit[((v << 2) | u_in, (v << 2) | u_out)] = w_lit
            under_step_lit[((v << 2) | u_out, (v << 2) | u_in)] = -w_lit
        else:
            s0 = 1 if vt == "Z" else 0
            d0 = (v << 2) | s0
            d2 = (v << 2) | ((s0 + 2) & 3)
            a1 = dart_arc[d0] + 1
            a2 = dart_arc[d2] + 1
            lit1 = a1 if dart_dir[d0] == 1 else -a1
            lit2 = -a2 if dart_dir[d2] == 1 else a2
            w_rel = cyclically_reduce_word((lit1, lit2))
            if w_rel:
                rels.append(w_rel)

    # Extract peripheral loop generators on N_0(F_i) for each orientable component:
    comp_periph: list[tuple[str, bool, int, list[tuple[int, ...]]]] = []
    for ci in range(num_scomps - d.free_loops):
        ct = topo.components[ci]
        if not ct.orientable:
            continue
        c_darts = [hd for hd in range(4 * n) if dc[hd] == ci]
        if not c_darts:
            continue
        root = c_darts[0]
        mu_gen = dart_arc[root] + 1
        adj: dict[int, list[tuple[int, tuple[int, ...], tuple]]] = {hd: [] for hd in c_darts}
        for hd in c_darts:
            od = d.opp[hd]
            if hd < od:
                adj[hd].append((od, (), ("e", hd, od)))
                adj[od].append((hd, (), ("e", hd, od)))
        for v in range(n):
            if is_cl[v]:
                for s0, s2 in ((0, 2), (1, 3)):
                    h0, h2 = 4 * v + s0, 4 * v + s2
                    if dc[h0] == ci:
                        w02 = (under_step_lit[(h0, h2)],) if (h0, h2) in under_step_lit else ()
                        w20 = (-w02[0],) if w02 else ()
                        adj[h0].append((h2, w02, ("cl", v, s0)))
                        adj[h2].append((h0, w20, ("cl", v, s0)))
            else:
                if dc[4 * v] == ci:
                    for s in range(4):
                        h0, h1 = 4 * v + s, 4 * v + ((s + 1) & 3)
                        adj[h0].append((h1, (), ("sd", v, s)))
                        adj[h1].append((h0, (), ("sd", v, s)))

        word_to: dict[int, tuple[int, ...]] = {root: ()}
        used_ekeys = set()
        bfs_q = [root]
        for u in bfs_q:
            for nb, w_step, ekey in adj[u]:
                if nb not in word_to:
                    word_to[nb] = word_to[u] + w_step
                    used_ekeys.add(ekey)
                    bfs_q.append(nb)

        loop_words: list[tuple[int, ...]] = []
        seen_ekeys = set()
        for u in c_darts:
            for nb, w_step, ekey in adj[u]:
                if ekey not in used_ekeys and ekey not in seen_ekeys:
                    seen_ekeys.add(ekey)
                    raw_lw = word_to[u] + w_step + tuple(-x for x in reversed(word_to[nb]))
                    st_red: list[int] = []
                    for lit in raw_lw:
                        if st_red and st_red[-1] == -lit:
                            st_red.pop()
                        else:
                            st_red.append(lit)
                    if st_red:
                        loop_words.append(tuple(st_red))
        comp_periph.append((ct.homeomorphism_type, ct.normal_euler == 0, mu_gen, loop_words))

    # Precompute the 7 conjugacy classes `(rep, size)` of SL(2, F_3) (inner automorphisms preserve
    # all subgroup orders `mu_ord` and `p_ord` identically!):
    sl2_seen = [False] * 24
    sl2_conj_classes: list[tuple[int, int]] = []
    for x0 in range(24):
        if sl2_seen[x0]:
            continue
        cls_s = set()
        for h_el in range(24):
            cx = _SL2_MUL[_SL2_MUL[h_el][x0]][_SL2_INV[h_el]]
            cls_s.add(cx)
            sl2_seen[cx] = True
        sl2_conj_classes.append((x0, len(cls_s)))

    # Solve Hom(pi_1_arcs, SL(2, F_3)) with unit/domain propagation + root conjugacy symmetry on `0..g-1`:
    rel_gens = [tuple(abs(x) - 1 for x in r) for r in rels]
    gen_occ = [0] * g
    for rg in rel_gens:
        for gi in rg:
            gen_occ[gi] += 1
    first_gen = max(range(g), key=lambda i: (gen_occ[i], -i))

    assign = [-1] * g
    weighted_homs: list[tuple[tuple[int, ...], int]] = []

    def eval_lit(lit: int) -> int:
        v = assign[abs(lit) - 1]
        return v if lit > 0 else _SL2_INV[v]

    def backtrack(mult: int, is_root: bool = False) -> None:
        changed: list[int] = []
        best_domain: tuple[int, int, list[int]] | None = None
        while True:
            progressed = False
            best_domain = None
            for r, rg in zip(rels, rel_gens):
                unassigned_pos = [k for k, gi in enumerate(rg) if assign[gi] == -1]
                if not unassigned_pos:
                    cur = _SL2_ID
                    for lit in r:
                        cur = _SL2_MUL[cur][eval_lit(lit)]
                    if cur != _SL2_ID:
                        for gi in changed:
                            assign[gi] = -1
                        return
                elif len(unassigned_pos) == 1:
                    pos = unassigned_pos[0]
                    gi = rg[pos]
                    pref = _SL2_ID
                    for lit in r[:pos]:
                        pref = _SL2_MUL[pref][eval_lit(lit)]
                    suff = _SL2_ID
                    for lit in r[pos + 1 :]:
                        suff = _SL2_MUL[suff][eval_lit(lit)]
                    val_lit = _SL2_MUL[_SL2_INV[pref]][_SL2_INV[suff]]
                    assign[gi] = val_lit if r[pos] > 0 else _SL2_INV[val_lit]
                    changed.append(gi)
                    progressed = True
                else:
                    u_gens = {rg[k] for k in unassigned_pos}
                    if len(u_gens) == 1:
                        gi = next(iter(u_gens))
                        cands: list[int] = []
                        for cand in range(24):
                            assign[gi] = cand
                            cur = _SL2_ID
                            for lit in r:
                                cur = _SL2_MUL[cur][eval_lit(lit)]
                            if cur == _SL2_ID:
                                cands.append(cand)
                        assign[gi] = -1
                        if not cands:
                            for g_chg in changed:
                                assign[g_chg] = -1
                            return
                        if len(cands) == 1:
                            assign[gi] = cands[0]
                            changed.append(gi)
                            progressed = True
                        elif best_domain is None or len(cands) < best_domain[0]:
                            best_domain = (len(cands), gi, cands)
            if not progressed:
                break

        rem_gens = [i for i in range(g) if assign[i] == -1]
        if not rem_gens:
            weighted_homs.append((tuple(assign), mult))
        elif is_root:
            rg0 = max(rem_gens, key=lambda i: (gen_occ[i], -i))
            for rep_c, cls_sz in sl2_conj_classes:
                assign[rg0] = rep_c
                backtrack(mult * cls_sz, False)
            assign[rg0] = -1
        elif best_domain is not None:
            _, nxt_g, cands = best_domain
            for cand in cands:
                assign[nxt_g] = cand
                backtrack(mult, False)
            assign[nxt_g] = -1
        else:
            # Pick unassigned generator participating in the most partially-assigned / shortest relators:
            scores = {i: [0, 0, gen_occ[i], -i] for i in rem_gens}
            for rg in rel_gens:
                u_set = [gi for gi in set(rg) if assign[gi] == -1]
                if len(u_set) == 2:
                    for gi in u_set:
                        scores[gi][0] += 1
                elif 0 < len(u_set) < len(set(rg)):
                    for gi in u_set:
                        scores[gi][1] += 1
            nxt_g = max(rem_gens, key=lambda i: scores[i])
            for cand in range(24):
                assign[nxt_g] = cand
                backtrack(mult, False)
            assign[nxt_g] = -1

        for gi in changed:
            assign[gi] = -1

    # Also factor across disconnected vertex components if `d` is a split union (`num_v_comps > 1`)!
    from surfaceknots.moves import connected_diagram_components

    v_comps = [c for c in connected_diagram_components(d) if c.num_vertices > 0]
    if len(v_comps) > 1:
        # Combine exact `(count, Counter-like multiset of peripheral records)` across split components:
        from collections import Counter

        tot_homs = free_mult
        dist: Counter[tuple[tuple[str, int, int], ...]] = Counter({(): free_mult})
        for vc in v_comps:
            vc_cnt, vc_per = peripheral_sl2f3_profile(vc)
            tot_homs *= vc_cnt
            vc_dist: Counter[tuple[tuple[str, int, int], ...]] = Counter(vc_per)
            vc_triv = vc_cnt - len(vc_per)
            if vc_triv > 0:
                vc_dist[()] += vc_triv
            new_dist: Counter[tuple[tuple[str, int, int], ...]] = Counter()
            for r1, m1 in dist.items():
                for r2, m2 in vc_dist.items():
                    comb_r = tuple(sorted(r1 + r2))
                    new_dist[comb_r] += m1 * m2
            dist = new_dist
        out_recs: list[tuple[tuple[str, int, int], ...]] = []
        for r_tup in sorted(k for k in dist if k):
            out_recs.extend([r_tup] * dist[r_tup])
        return (tot_homs, tuple(out_recs))

    backtrack(free_mult, True)

    periph_records: list[tuple[tuple[str, int, int], ...]] = []
    total_hom_count = 0
    for h, w_mult in weighted_homs:
        total_hom_count += w_mult

        def ev_w(w: tuple[int, ...]) -> int:
            c = _SL2_ID
            for lit in w:
                v = h[abs(lit) - 1]
                c = _SL2_MUL[c][v if lit > 0 else _SL2_INV[v]]
            return c

        rec: list[tuple[str, int, int]] = []
        for htype, is_e0, mu_g, lws in comp_periph:
            mu_el = h[mu_g - 1]
            mu_ord = _sl2_subgroup_size([mu_el])
            lw_els = [ev_w(lw) for lw in lws]
            if is_e0:
                for le in lw_els:
                    assert _SL2_MUL[mu_el][le] == _SL2_MUL[le][mu_el], (
                        "Meridian failed to commute with peripheral longitude!"
                    )
            p_ord = _sl2_subgroup_size([mu_el] + lw_els)
            if p_ord > mu_ord:
                rec.append((htype, mu_ord, p_ord))
        if rec:
            rec.sort()
            periph_records.extend([tuple(rec)] * w_mult)

    periph_records.sort()
    return (total_hom_count, tuple(periph_records))


# ---------------------------------------------------------------------------
# Movie tag propagation, mod-2 triple linking, and Sanderson mod-2 double linking
# ---------------------------------------------------------------------------


def _transport_component_tags_forward(
    L0: ClassicalLinkDiagram,
    movie: list[UnlinkMovieStep],
) -> list[tuple[int, ...]]:
    """Propagate `dart_surface_comp` forward along every step of `movie` (`0..len(movie)`)."""
    assert L0.dart_surface_comp is not None
    tags: list[tuple[int, ...]] = [L0.dart_surface_comp]
    cur_sc = L0.dart_surface_comp
    for step in movie:
        bef, aft = step.before, step.after
        m = aft.num_crossings
        nxt_sc = [0] * (4 * m)
        if step.kind in ("R1-", "R2-"):
            o2n = step.meta["old_to_new"]
            for vo, vn in o2n.items():
                for s in range(4):
                    nxt_sc[4 * vn + s] = cur_sc[4 * vo + s]
        elif step.kind == "R1+":
            n = bef.num_crossings
            for hd in range(4 * n):
                nxt_sc[hd] = cur_sc[hd]
            d_ext = next(aft.opp[4 * n + s] for s in range(4) if aft.opp[4 * n + s] < 4 * n)
            c_kink = cur_sc[d_ext]
            for s in range(4):
                nxt_sc[4 * n + s] = c_kink
        elif step.kind == "R2+":
            n = bef.num_crossings
            for hd in range(4 * n):
                nxt_sc[hd] = cur_sc[hd]
            u, w = n, n + 1
            ca = cur_sc[step.meta["da"]]
            cb = cur_sc[step.meta["db"]]
            nxt_sc[4 * u + 0] = nxt_sc[4 * u + 2] = ca
            nxt_sc[4 * w + 0] = nxt_sc[4 * w + 2] = ca
            nxt_sc[4 * u + 1] = nxt_sc[4 * u + 3] = cb
            nxt_sc[4 * w + 1] = nxt_sc[4 * w + 3] = cb
        elif step.kind == "R3":
            vs = step.meta["vs"]
            ss = step.meta["ss"]
            nxt_sc = list(cur_sc)
            cA = cur_sc[(vs[0] << 2) | ss[0]]
            cB = cur_sc[(vs[1] << 2) | ss[1]]
            cC = cur_sc[(vs[2] << 2) | ss[2]]
            strand_comps = ((cC, cB), (cA, cC), (cB, cA))
            for idx_i, wi in enumerate(vs):
                s02, s13 = strand_comps[idx_i]
                nxt_sc[4 * wi + 0] = nxt_sc[4 * wi + 2] = s02
                nxt_sc[4 * wi + 1] = nxt_sc[4 * wi + 3] = s13
        else:
            raise ValueError(f"Unknown movie step kind {step.kind!r}")
        cur_sc = tuple(nxt_sc)
        tags.append(cur_sc)
    return tags


def mod2_triple_linking_profile(d: MarkedDiagram) -> tuple[tuple[str, str, str], ...] | None:
    """Compute Carter-Saito-Satoh mod-2 triple point linking invariants `t_{ijk}(F) in Z_2`.

    Counts `R3` triple points `(Top, Mid, Bot) = (F_i, F_j, F_k)` mod 2 across `(M_-, M_+)`
    for all ordered triples of 3 distinct components (`|{i, j, k}| = 3`), canonicalized up to
    component permutation preserving `(homeomorphism_type, normal_euler)`.
    """
    topo = d.surface_topology()
    c = topo.num_components
    if c < 3:
        return ()

    lp, lm = d.resolve(1), d.resolve(-1)
    mp = find_unlink_movie(lp)
    mm = find_unlink_movie(lm)
    if mp is None or mm is None:
        return None

    counts: dict[tuple[int, int, int], int] = {}
    for L0, movie in ((lp, mp), (lm, mm)):
        sc_tags = _transport_component_tags_forward(L0, movie)
        for idx, step in enumerate(movie):
            if step.kind != "R3":
                continue
            sc = sc_tags[idx]
            vs = step.meta["vs"]
            ss = step.meta["ss"]
            hs = step.meta["heights"]
            strands = [
                (hs[0], sc[(vs[0] << 2) | ss[0]]),
                (hs[1], sc[(vs[1] << 2) | ss[1]]),
                (hs[2], sc[(vs[2] << 2) | ss[2]]),
            ]
            strands.sort(reverse=True)  # Top (2), Mid (1), Bot (0)
            ci, cj, ck = strands[0][1], strands[1][1], strands[2][1]
            if ci != cj and cj != ck and ci != ck:
                counts[(ci, cj, ck)] = counts.get((ci, cj, ck), 0) ^ 1

    nonzero = [k for k, v in counts.items() if v == 1]
    if not nonzero:
        return ()

    # Restrict to active components participating in nonzero triples (invariant under free_loops!)
    # and canonicalize over label-preserving permutations AND O(4) mirror flip `e -> -e`:
    active_comps = sorted({u for trip in nonzero for u in trip})
    out_deg = {u: sum(1 for (a, _, _) in nonzero if a == u) for u in active_comps}
    mid_deg = {u: sum(1 for (_, b, _) in nonzero if b == u) for u in active_comps}
    in_deg = {u: sum(1 for (_, _, cc) in nonzero if cc == u) for u in active_comps}

    best: tuple[tuple[str, str, str], ...] | None = None
    for sgn in (1, -1):
        lbl_map = {
            u: f"{topo.components[u].homeomorphism_type}(e={sgn * topo.components[u].normal_euler})"
            for u in active_comps
        }
        unique_lbls = sorted(set(lbl_map.values()))
        buckets = [[u for u in active_comps if lbl_map[u] == lb] for lb in unique_lbls]
        total_perms = math.prod(math.factorial(len(b)) for b in buckets)
        if total_perms > 5040:
            # Refine buckets canonically by directed incidence degree triple `(out, mid, in)`
            refined: list[list[int]] = []
            for b in buckets:
                by_deg: dict[tuple[int, int, int], list[int]] = {}
                for u in b:
                    by_deg.setdefault((out_deg[u], mid_deg[u], in_deg[u]), []).append(u)
                for dk in sorted(by_deg):
                    refined.append(by_deg[dk])
            buckets = refined
            total_perms = math.prod(math.factorial(len(b)) for b in buckets)

        if total_perms <= 5040:
            bucket_iters = [permutations(b) for b in buckets]
            for p_buckets in product(*bucket_iters):
                order = [u for grp in p_buckets for u in grp]
                u_to_r = {u: r for r, u in enumerate(order)}
                cand = tuple(
                    sorted(
                        (
                            f"{u_to_r[i]}:{lbl_map[i]}",
                            f"{u_to_r[j]}:{lbl_map[j]}",
                            f"{u_to_r[k]}:{lbl_map[k]}",
                        )
                        for (i, j, k) in nonzero
                    )
                )
                if best is None or cand < best:
                    best = cand
        else:
            order = [u for b in buckets for u in b]
            u_to_r = {u: r for r, u in enumerate(order)}
            cand = tuple(
                sorted(
                    (
                        f"{u_to_r[i]}:{lbl_map[i]}",
                        f"{u_to_r[j]}:{lbl_map[j]}",
                        f"{u_to_r[k]}:{lbl_map[k]}",
                    )
                    for (i, j, k) in nonzero
                )
            )
            if best is None or cand < best:
                best = cand

    return best if best is not None else ()


def classical_knot_arf(c_diag: ClassicalLinkDiagram) -> int:
    """Compute Arf invariant in `{0, 1}` of a 1-component classical knot diagram via `|det(K)| mod 8`."""
    cycs, total = c_diag.circle_components()
    assert total == 1, f"Expected 1 circle component, got {total}"
    n = c_diag.num_crossings
    if n <= 2:
        return 0
    cyc = cycs[0]
    under_pos = [k for k, hd in enumerate(cyc) if not _is_over(c_diag.vtype[hd >> 2], hd & 3)]
    g = len(under_pos)
    if g <= 1:
        return 0
    dart_arc = [0] * (4 * n)
    for idx, p0 in enumerate(under_pos):
        p1 = under_pos[(idx + 1) % g]
        dart_arc[cyc[p0] ^ 2] = idx
        k = (p0 + 1) % len(cyc)
        while k != p1:
            hd = cyc[k]
            dart_arc[hd] = idx
            dart_arc[hd ^ 2] = idx
            k = (k + 1) % len(cyc)
        dart_arc[cyc[p1]] = idx

    rows = []
    for v in range(n):
        vt = c_diag.vtype[v]
        us = (0, 2) if vt == "X0" else (1, 3)
        os = (1, 3) if vt == "X0" else (0, 2)
        a_ov = dart_arc[4 * v + os[0]]
        a_u1 = dart_arc[4 * v + us[0]]
        a_u2 = dart_arc[4 * v + us[1]]
        row = [0] * g
        row[a_ov] += 2
        row[a_u1] -= 1
        row[a_u2] -= 1
        rows.append(row)
    mat = sp.Matrix([r[:-1] for r in rows[:-1]])
    det_val = abs(int(mat.det()))
    assert det_val % 2 == 1, f"Expected odd determinant for knot, got {det_val}"
    return 1 if (det_val % 8) in (3, 5) else 0


def _sky_potential_mod2_on_marked(
    d: MarkedDiagram,
    sign: int,
    comp_i: int,
    comp_j: int,
) -> list[int] | None:
    """Compute backward sky-potential `p[hd] in {0, 1}` on darts of `d` on `comp_i` wrt `comp_j`."""
    assert comp_i != comp_j
    L0 = d.resolve(1 if sign > 0 else -1)
    movie = find_unlink_movie(L0)
    if movie is None:
        return None
    sc_tags = _transport_component_tags_forward(L0, movie)

    pot_cur: list[int] = []
    for idx in range(len(movie) - 1, -1, -1):
        step = movie[idx]
        bef = step.before
        sc_bef = sc_tags[idx]
        n = bef.num_crossings
        pot_bef = [0] * (4 * n)

        def is_i_under_j(v: int, s: int) -> int:
            if sc_bef[4 * v + s] != comp_i:
                return 0
            if _is_over(bef.vtype[v], s):
                return 0
            other_c = sc_bef[4 * v + ((s + 1) & 3)]
            return 1 if other_c == comp_j else 0

        if step.kind == "R1-":
            o2n = step.meta["old_to_new"]
            removed_v = next(vo for vo in range(n) if vo not in o2n)
            for vo, vn in o2n.items():
                for s in range(4):
                    pot_bef[4 * vo + s] = pot_cur[4 * vn + s]
            ext_val = 0
            for s in range(4):
                od = bef.opp[4 * removed_v + s]
                if (od >> 2) != removed_v:
                    ext_val = pot_bef[od]
                    break
            for s in range(4):
                pot_bef[4 * removed_v + s] = ext_val
        elif step.kind == "R2-":
            o2n = step.meta["old_to_new"]
            rem = [vo for vo in range(n) if vo not in o2n]
            for vo, vn in o2n.items():
                for s in range(4):
                    pot_bef[4 * vo + s] = pot_cur[4 * vn + s]
            u, w = rem[0], rem[1]
            for rv in (u, w):
                for s in range(4):
                    od = bef.opp[4 * rv + s]
                    if (od >> 2) not in (u, w):
                        pot_bef[4 * rv + s] = pot_bef[od]
            for _ in range(4):
                for rv in (u, w):
                    for s in range(4):
                        hd = 4 * rv + s
                        hd_opp_slot = 4 * rv + (s ^ 2)
                        od = bef.opp[hd]
                        if (od >> 2) not in (u, w):
                            pot_bef[hd] = pot_bef[od]
                            pot_bef[hd_opp_slot] = pot_bef[hd] ^ is_i_under_j(rv, s)
                            pot_bef[bef.opp[hd_opp_slot]] = pot_bef[hd_opp_slot]
            for su in range(4):
                su2 = (su + 1) & 3
                od1 = bef.opp[4 * u + su]
                if (od1 >> 2) == w and bef.opp[4 * u + su2] == ((w << 2) | (((od1 & 3) - 1) & 3)):
                    sw = od1 & 3
                    sw2 = (sw - 1) & 3
                    for s_in, sw_in in ((su, sw), (su2, sw2)):
                        u_ext = (u << 2) | (s_in ^ 2)
                        w_ext = (w << 2) | (sw_in ^ 2)
                        if bef.opp[u_ext] == w_ext:
                            pot_bef[u_ext] = 0
                            pot_bef[w_ext] = 0
                            pot_bef[4 * u + s_in] = is_i_under_j(u, s_in)
                            pot_bef[4 * w + sw_in] = is_i_under_j(u, s_in)
                    break
        elif step.kind == "R2+":
            for hd in range(4 * n):
                if hd not in (step.meta["da"], step.meta["oda"], step.meta["db"], step.meta["odb"]):
                    pot_bef[hd] = pot_cur[hd]
            u, w = n, n + 1
            pot_bef[step.meta["da"]] = pot_cur[4 * u + 0]
            pot_bef[step.meta["oda"]] = pot_cur[4 * w + 0]
            pot_bef[step.meta["db"]] = pot_cur[4 * w + 1]
            pot_bef[step.meta["odb"]] = pot_cur[4 * u + 3]
        elif step.kind == "R3":
            vs = set(step.meta["vs"])
            ws = step.meta["vs"]
            ports = step.meta["ports"]
            port_pairs = step.meta["port_pairs"]
            for v in range(n):
                if v not in vs:
                    for s in range(4):
                        pot_bef[4 * v + s] = pot_cur[4 * v + s]
            for idx_i in range(3):
                wi = ws[idx_i]
                pa, pb = port_pairs[idx_i]
                pot_bef[ports[pa]] = pot_cur[4 * wi + 2]
                pot_bef[ports[pb]] = pot_cur[4 * wi + 3]
            for rv in ws:
                for s in range(4):
                    hd = 4 * rv + s
                    if hd in ports:
                        hd_in = 4 * rv + (s ^ 2)
                        pot_bef[hd_in] = pot_bef[hd] ^ is_i_under_j(rv, s)
        pot_cur = pot_bef

    dc, _ = d._dart_surface_components()
    nd = d.num_vertices
    pot_d = [0] * (4 * nd)
    classical_to_l0: dict[int, int] = {}
    cid = 0
    for v in range(nd):
        if d.vtype[v] in ("X0", "X1"):
            classical_to_l0[v] = cid
            cid += 1

    seeded = [False] * (4 * nd)
    stack: list[int] = []
    for v, lv in classical_to_l0.items():
        for s in range(4):
            hd = 4 * v + s
            if dc[hd] == comp_i:
                pot_d[hd] = pot_cur[4 * lv + s]
                seeded[hd] = True
                stack.append(hd)

    def marked_partner(v: int, s: int) -> int:
        vt = d.vtype[v]
        p03 = (vt == "Y" and sign > 0) or (vt == "Z" and sign < 0)
        return (3 - s) if p03 else (s ^ 1)

    while stack:
        hd = stack.pop()
        nhd = d.opp[hd]
        if not seeded[nhd]:
            pot_d[nhd] = pot_d[hd]
            seeded[nhd] = True
            stack.append(nhd)
        v, s = hd >> 2, hd & 3
        if d.vtype[v] in ("Y", "Z"):
            shd = 4 * v + marked_partner(v, s)
            if not seeded[shd]:
                pot_d[shd] = pot_d[hd]
                seeded[shd] = True
                stack.append(shd)
    return pot_d


def compute_torus_dlk_pair(d: MarkedDiagram, comp_i: int, comp_j: int) -> int | None:
    """Compute the Rohlin-Matsumoto / Carter-Saito-Sato mod-2 torus-linking invariant `q_{F_i}(PD(lambda_{ij}))`
    in `{0, 1}` when `F_i` is an `e = 0` torus and the mod-2 linking class `lambda_{ij} in H^1(F_i; Z_2)`
    with `F_j` is non-zero.

    Mathematical proof of ambient isotopy invariance:
    1. By Alexander duality in `S^4`, mod-2 linking with `F_j` defines a canonical ambient-isotopy-invariant
       linear functional `lambda_{ij} : H_1(F_i; Z_2) -> Z_2`. Because `L_-` and `L_+` are classical unlinks
       (all pairwise linking numbers zero), the indicator cochain of under-crossings of `F_i` beneath `F_j`
       vanishes on every `L_-` circle, `L_+` circle, and marked-vertex 4-cycle, so it is an exact cellular
       1-cocycle representing `lambda_{ij}`.
    2. When `F_i ~ T^2` and `lambda_{ij} != 0`, `ker(lambda_{ij}) subset H_1(F_i; Z_2) ~ Z_2^2` is 1-dimensional
       and contains a unique non-zero class `alpha_0 = PD(lambda_{ij})` (the homology class of the double-point
       curve `F_i cap W_j` for any mod-2 3-chain `W_j` with `partial W_j = F_j`, since the mod-2 intersection
       form on `T^2` is alternating).
    3. On an orientable surface `F_i subset R^4`, the Rohlin-Matsumoto pin^- membrane quadratic form
       `q_{F_i} : H_1(F_i; Z_2) -> Z_2` satisfies `q(x + y) = q(x) + q(y) + (x . y) mod 2`. For any circle `C`
       in a regular horizontal cross-section `F_i cap (R^3 x {t_0}) = C sqcup C'_1 sqcup ... sqcup C'_r`, an
       embedded Seifert surface `M subset R^3 x {t_0}` for `C` has zero interior self-intersections and zero
       normal framing obstruction relative to `partial/partial t in T F_i`, so
       `q_{F_i}([C]) = sum_{s=1}^r lk_{R^3}(C, C'_s) mod 2` (using signed 3D linking numbers, with no Arf term).
    4. When `lambda_{ij} == 0` (an ambient-isotopy-invariant condition, e.g. on `9_1^{0,1}`, `10_1^{0,1}`,
       `10_2^{0,1}`, and `11_1^{1,-2}`), `ker(lambda_{ij}) = H_1(F_i; Z_2)` has no distinguished non-zero class,
       so this function returns `None` by design and the pair is omitted from `double_linking_mod2`.
    """
    topo = d.surface_topology()
    ci = topo.components[comp_i]
    if not (ci.orientable and ci.genus_symbol == 1 and ci.normal_euler == 0):
        return None

    dc, _ = d._dart_surface_components()
    nd = d.num_vertices
    darts_i = [hd for hd in range(4 * nd) if dc[hd] == comp_i]
    if not darts_i:
        return None

    eidx: dict[tuple[int, int, str], int] = {}

    def get_edge(a: int, b: int, tag: str) -> int:
        key = (a, b, tag) if a <= b else (b, a, tag)
        res = eidx.get(key)
        if res is None:
            res = len(eidx)
            eidx[key] = res
        return res

    lam: dict[int, int] = {}
    for hd in darts_i:
        get_edge(hd, d.opp[hd], "e")
    for v in range(nd):
        vt = d.vtype[v]
        if vt in ("X0", "X1"):
            for s in (0, 1):
                if dc[4 * v + s] == comp_i:
                    ei = get_edge(4 * v + s, 4 * v + s + 2, "c")
                    is_und_j = (not _is_over(vt, s)) and (dc[4 * v + ((s + 1) & 3)] == comp_j)
                    lam[ei] = 1 if is_und_j else 0
        elif dc[4 * v] == comp_i:
            for s in range(4):
                get_edge(4 * v + s, 4 * v + ((s + 1) & 3), "m")

    marked_vs = [v for v in range(nd) if d.vtype[v] in ("Y", "Z") and dc[4 * v] == comp_i]

    def slice_cycles(s_plus: set[int]) -> list[list[tuple[int, int]]]:
        seen: set[int] = set()
        cycs: list[list[tuple[int, int]]] = []
        for hd in darts_i:
            if hd in seen:
                continue
            cyc: list[tuple[int, int]] = []
            cur = hd
            while cur not in seen:
                seen.add(cur)
                v, s = cur >> 2, cur & 3
                vt = d.vtype[v]
                if vt in ("X0", "X1"):
                    ex = (v << 2) | (s ^ 2)
                else:
                    sgn = 1 if v in s_plus else -1
                    p03 = (vt == "Y" and sgn > 0) or (vt == "Z" and sgn < 0)
                    ex = (v << 2) | ((3 - s) if p03 else (s ^ 1))
                seen.add(ex)
                cyc.append((cur, ex))
                cur = d.opp[ex]
            cycs.append(cyc)
        return cycs

    def cyc_mask(cyc: list[tuple[int, int]]) -> int:
        m = 0
        for a, b in cyc:
            v = a >> 2
            tag = "c" if d.vtype[v] in ("X0", "X1") else "m"
            m ^= 1 << get_edge(a, b, tag)
            m ^= 1 << get_edge(b, d.opp[b], "e")
        return m

    boundaries: list[int] = []
    for s_plus in (set(), set(marked_vs)):
        for cyc in slice_cycles(s_plus):
            boundaries.append(cyc_mask(cyc))
    for v in marked_vs:
        m_sq = 0
        for s in range(4):
            m_sq ^= 1 << get_edge(4 * v + s, 4 * v + ((s + 1) & 3), "m")
        boundaries.append(m_sq)

    basis: dict[int, int] = {}
    for b_vec in boundaries:
        x = b_vec
        while x:
            hb = x.bit_length() - 1
            if hb in basis:
                x ^= basis[hb]
            else:
                basis[hb] = x
                break

    pivots_desc = sorted(basis.keys(), reverse=True)

    def canon_class(x: int) -> int:
        for hb in pivots_desc:
            if (x >> hb) & 1:
                x ^= basis[hb]
        return x

    def eval_lam(m: int) -> int:
        tot = 0
        for ei, w in lam.items():
            if w and ((m >> ei) & 1):
                tot ^= 1
        return tot

    for b_vec in boundaries:
        if eval_lam(b_vec) != 0:
            return None

    subsets: list[set[int]] = []
    if len(marked_vs) <= 12:
        for r_sz in range(len(marked_vs) + 1):
            for comb in combinations(marked_vs, r_sz):
                subsets.append(set(comb))
    else:
        # Exact O(m) spanning-tree + co-tree fundamental cycle construction over G_- and G_+:
        # Merging all L_- disks along a spanning tree `tree_minus` leaves a single disk D_0 subset T^2,
        # and each co-tree saddle `e not in tree_minus` attaches a 1-handle whose boundary curves
        # represent its fundamental homology class in H_1(T^2; F_2) (surjective onto F_2^2!).
        cycs_0 = slice_cycles(set())
        which_0: dict[int, int] = {}
        for idx_c, cyc in enumerate(cycs_0):
            for a, b in cyc:
                which_0[a] = idx_c
                which_0[b] = idx_c
        uf_p = list(range(len(cycs_0)))

        def _find_u(x: int) -> int:
            while uf_p[x] != x:
                uf_p[x] = uf_p[uf_p[x]]
                x = uf_p[x]
            return x

        tree_minus: set[int] = set()
        for v in marked_vs:
            cs_v = list({which_0[4 * v + s] for s in range(4)})
            if len(cs_v) == 2:
                ra, rb = _find_u(cs_v[0]), _find_u(cs_v[1])
                if ra != rb:
                    uf_p[ra] = rb
                    tree_minus.add(v)
        subsets.append(set(tree_minus))
        cotree_minus = [v for v in marked_vs if v not in tree_minus]
        for v in cotree_minus:
            subsets.append(tree_minus | {v})

        # Dual spanning tree `tree_plus` on L_+ disks as redundant cross-check:
        all_mk = set(marked_vs)
        cycs_p = slice_cycles(all_mk)
        which_p: dict[int, int] = {}
        for idx_c, cyc in enumerate(cycs_p):
            for a, b in cyc:
                which_p[a] = idx_c
                which_p[b] = idx_c
        uf_pp = list(range(len(cycs_p)))

        def _find_up(x: int) -> int:
            while uf_pp[x] != x:
                uf_pp[x] = uf_pp[uf_pp[x]]
                x = uf_pp[x]
            return x

        tree_plus: set[int] = set()
        for v in marked_vs:
            cs_v = list({which_p[4 * v + s] for s in range(4)})
            if len(cs_v) == 2:
                ra, rb = _find_up(cs_v[0]), _find_up(cs_v[1])
                if ra != rb:
                    uf_pp[ra] = rb
                    tree_plus.add(v)
        for v in marked_vs:
            if v not in tree_plus:
                subsets.append(all_mk - (tree_plus | {v}))

    self_class_vs = [
        v
        for v in range(nd)
        if d.vtype[v] in ("X0", "X1") and dc[4 * v] == comp_i and dc[4 * v + 1] == comp_i
    ]

    qvals: dict[int, set[int]] = {}
    lamvals: dict[int, set[int]] = {}

    for s_plus in subsets:
        cycs = slice_cycles(s_plus)
        which: dict[int, int] = {}
        indart: set[int] = set()
        for k_idx, cyc in enumerate(cycs):
            for a, b in cyc:
                which[a] = k_idx
                which[b] = k_idx
                indart.add(a)
        lk2: dict[tuple[int, int], int] = {}
        for v in self_class_vs:
            ka = which[4 * v]
            kb = which[4 * v + 1]
            if ka == kb:
                continue
            vt = d.vtype[v]
            us = (0, 2) if vt == "X0" else (1, 3)
            os_ = (1, 3) if vt == "X0" else (0, 2)
            u_in = us[0] if (4 * v + us[0]) in indart else us[1]
            o_in = os_[0] if (4 * v + os_[0]) in indart else os_[1]
            sgn = 1 if ((o_in - u_in) & 3) == 3 else -1
            key = (ka, kb) if ka < kb else (kb, ka)
            lk2[key] = lk2.get(key, 0) + sgn

        for k_idx, cyc in enumerate(cycs):
            raw_m = cyc_mask(cyc)
            cl = canon_class(raw_m)
            if cl == 0:
                continue
            q_val = 0
            for (ca, cb), val2 in lk2.items():
                if k_idx == ca or k_idx == cb:
                    q_val ^= (val2 // 2) & 1
            qvals.setdefault(cl, set()).add(q_val)
            lamvals.setdefault(cl, set()).add(eval_lam(raw_m))

    if not qvals:
        return None
    if not (all(len(s) == 1 for s in qvals.values()) and all(len(s) == 1 for s in lamvals.values())):
        return None

    ker_classes = [cl for cl in qvals if lamvals[cl] == {0}]
    nonker_classes = [cl for cl in qvals if lamvals[cl] == {1}]
    if not nonker_classes:
        return None
    if len(ker_classes) == 1:
        return next(iter(qvals[ker_classes[0]]))
    if len(ker_classes) == 0 and len(nonker_classes) == 2:
        a_cl, b_cl = nonker_classes
        return next(iter(qvals[a_cl])) ^ next(iter(qvals[b_cl])) ^ 1
    return None


def double_linking_mod2(d: MarkedDiagram) -> tuple[tuple[str, str, int], ...]:
    """Return sorted `(type_i, type_j, Dlk_{ij})` across all component pairs `(i, j)` with `F_i ~ T^2`
    whose mod-2 linking functional `lambda_{ij}: H_1(F_i; F_2) -> F_2` with `F_j` is non-zero.

    By design, pairs `(i, j)` where `lambda_{ij} == 0` on `H_1(T^2; F_2)` (such as `9_1^{0,1}`,
    `10_1^{0,1}`, `10_2^{0,1}`, `11_1^{1,-2}`, or split unions) have no distinguished non-zero kernel
    line `ker(lambda_{ij}) \\ {0}` and are omitted (`compute_torus_dlk_pair` returns `None`), which is
    itself a topological invariant of the surface-link.
    """
    topo = d.surface_topology()
    c = topo.num_components
    if c < 2:
        return ()
    raw: list[tuple[str, str, int, int]] = []
    for i in range(c):
        ci = topo.components[i]
        if not (ci.orientable and ci.genus_symbol == 1 and ci.normal_euler == 0):
            continue
        for j in range(c):
            if i == j:
                continue
            val = compute_torus_dlk_pair(d, i, j)
            if val is not None:
                cj = topo.components[j]
                raw.append((ci.homeomorphism_type, cj.homeomorphism_type, cj.normal_euler, val))
    if not raw:
        return ()
    prof_pos = tuple(sorted((ti, f"{tj}(e={ej})", v) for ti, tj, ej, v in raw))
    prof_neg = tuple(sorted((ti, f"{tj}(e={-ej})", v) for ti, tj, ej, v in raw))
    return min(prof_pos, prof_neg)


def exact_sublink_invariant_profile(d: MarkedDiagram) -> tuple:
    """Return an exact, BFS-free `O(4)`-invariant tuple of proper sublink invariants `(1-comp, 2-comp)`.

    For `c <= 1` components returns `()`. For `c >= 2`, computes the sorted multiset of 1-component
    profiles `(homeomorphism_type, abs(normal_euler), hom_counts, e1_hash)` over `0 <= i < c`; for
    `c >= 3`, also includes the sorted multiset of 2-component profiles over `0 <= i < j < c`.
    """
    import hashlib
    from surfaceknots.invariants.alexander import compute_alexander_invariants
    from surfaceknots.invariants.fundamental_group import finite_quotient_profile, wirtinger_presentation

    topo = d.surface_topology()
    c = topo.num_components
    if c <= 1:
        return ()

    singles: list[tuple] = []
    for i in range(c):
        di = extract_subsurface_diagram(d, (i,))
        cti = topo.components[i]
        wpi = wirtinger_presentation(di)
        fqpi = tuple(sorted(finite_quotient_profile(wpi).items()))
        ali = compute_alexander_invariants(di)
        e1i = hashlib.sha256(repr(ali.evaluation_fingerprint).encode("utf-8")).hexdigest()[:16]
        singles.append((cti.homeomorphism_type, abs(cti.normal_euler), fqpi, e1i))
    singles.sort()

    if c == 2:
        return (tuple(singles),)

    pairs: list[tuple] = []
    for i, j in combinations(range(c), 2):
        dij = extract_subsurface_diagram(d, (i, j))
        st_ij = dij.surface_topology()
        ge_pos = tuple(sorted((ct.genus_symbol, ct.normal_euler) for ct in st_ij.components))
        ge_neg = tuple(sorted((ct.genus_symbol, -ct.normal_euler) for ct in st_ij.components))
        t_key = min(ge_pos, ge_neg)
        wpij = wirtinger_presentation(dij)
        fqpij = tuple(sorted(finite_quotient_profile(wpij).items()))
        alij = compute_alexander_invariants(dij)
        e1ij = hashlib.sha256(repr(alij.evaluation_fingerprint).encode("utf-8")).hexdigest()[:16]
        pairs.append((t_key, fqpij, e1ij))
    pairs.sort()
    return (tuple(singles), tuple(pairs))

