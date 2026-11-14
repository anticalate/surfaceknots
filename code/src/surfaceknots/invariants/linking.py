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
4. `double_linking_mod2`: Sanderson / Carter-Saito-Satoh / Takase-Tanaka mod-2 double linking invariant
   `Dlk_{ij}(F) in Z_2` via backward sky-potential transport of `lambda_{ij}: H_1(F_i; Z_2) -> Z_2`
   and Rokhlin's spin quadratic form `q_i: H_1(F_i; Z_2) -> Z_2`, separating `8_1^{1,1}` (`(0, 0)`)
   from the 1-twist-spun Hopf link `10_1^{1,1}` (`(1, 1)`).
"""

from __future__ import annotations

from itertools import combinations, permutations
import sympy as sp

from surfaceknots.diagrams import ClassicalLinkDiagram, MarkedDiagram
from surfaceknots.invariants.fundamental_group import cyclically_reduce_word
from surfaceknots.moves import simplify_and_explore_yoshikawa
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
    # Standard unknotted #_k RP^2 with normal Euler e has ch = 2*k (minimal ch = 2 for RP^2, 4 for K^2)
    # and cyclic Z_2 fundamental group!
    k_non = abs(c.genus_symbol)
    if d.num_vertices == 2 * k_non:
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
                # For r_sz >= 2, check if completely split into unknotted singletons!
                completely_trivial = False
                if is_split or best.num_vertices == 0:
                    # Check each singleton of `comb`
                    completely_trivial = all(
                        _is_certified_standard_unknot_single(
                            simplify_and_explore_yoshikawa(extract_subsurface_diagram(d, (idx_i,)), max_states=max_states)[0]
                        )
                        for idx_i in comb
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

    # Solve Hom(pi_1_arcs, SL(2, F_3)) with unit-propagation + root conjugacy symmetry on `0..g-1`:
    rel_gens = [tuple(abs(x) - 1 for x in r) for r in rels]
    assign = [-1] * g
    weighted_homs: list[tuple[tuple[int, ...], int]] = []

    def eval_lit(lit: int) -> int:
        v = assign[abs(lit) - 1]
        return v if lit > 0 else _SL2_INV[v]

    def backtrack(mult: int) -> None:
        changed: list[int] = []
        while True:
            progressed = False
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
            if not progressed:
                break

        nxt_g = next((i for i in range(g) if assign[i] == -1), -1)
        if nxt_g == -1:
            weighted_homs.append((tuple(assign), mult))
        elif nxt_g == 0:
            for rep_c, cls_sz in sl2_conj_classes:
                assign[0] = rep_c
                backtrack(mult * cls_sz)
            assign[0] = -1
        else:
            for cand in range(24):
                assign[nxt_g] = cand
                backtrack(mult)
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

    backtrack(free_mult)

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

    # Canonicalize under permutations preserving `(homeomorphism_type, normal_euler)` if c <= 6:
    comp_labels = [f"{ct.homeomorphism_type}(e={ct.normal_euler})" for ct in topo.components]
    if c <= 6:
        best = None
        for perm in permutations(range(c)):
            if any(comp_labels[i] != comp_labels[perm[i]] for i in range(c)):
                continue
            cand = tuple(
                sorted(
                    (comp_labels[perm[i]], comp_labels[perm[j]], comp_labels[perm[k]])
                    for (i, j, k) in nonzero
                )
            )
            if best is None or cand < best:
                best = cand
        return best if best is not None else ()
    return tuple(sorted((comp_labels[i], comp_labels[j], comp_labels[k]) for (i, j, k) in nonzero))


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
    """Compute Sanderson's mod-2 double linking `Dlk_{ij}(F) in {0, 1}` when `F_i` is an `e = 0` torus."""
    topo = d.surface_topology()
    ci = topo.components[comp_i]
    if not (ci.orientable and ci.genus_symbol == 1 and ci.normal_euler == 0):
        return None

    pot_minus = _sky_potential_mod2_on_marked(d, -1, comp_i, comp_j)
    pot_plus = _sky_potential_mod2_on_marked(d, 1, comp_i, comp_j)
    if pot_minus is None or pot_plus is None:
        return None

    dc, _ = d._dart_surface_components()
    nd = d.num_vertices
    edge_id_of_dart: dict[int, int] = {}
    num_e = 0
    for hd in range(4 * nd):
        if dc[hd] == comp_i and hd not in edge_id_of_dart:
            od = d.opp[hd]
            edge_id_of_dart[hd] = num_e
            edge_id_of_dart[od] = num_e
            num_e += 1
    if num_e == 0:
        return None

    def trace_res_circles(sgn: int) -> list[list[tuple[int, int]]]:
        seen = set()
        cycs = []
        for hd in range(4 * nd):
            if dc[hd] != comp_i or hd in seen:
                continue
            c = []
            cur = hd
            while cur not in seen:
                seen.add(cur)
                v, s = cur >> 2, cur & 3
                vt = d.vtype[v]
                if vt in ("X0", "X1"):
                    ex = (v << 2) | (s ^ 2)
                else:
                    p03 = (vt == "Y" and sgn > 0) or (vt == "Z" and sgn < 0)
                    ex = (v << 2) | ((3 - s) if p03 else (s ^ 1))
                seen.add(ex)
                c.append((cur, ex))
                cur = d.opp[ex]
            cycs.append(c)
        return cycs

    cycs_minus = trace_res_circles(-1)
    cycs_plus = trace_res_circles(1)
    v_minus_0 = num_e
    v_plus_0 = num_e + len(cycs_minus)
    num_kv = v_plus_0 + len(cycs_plus)

    k_edges: list[tuple[int, int, int]] = []
    transition_cell: dict[tuple[int, int], list[int]] = {}

    def add_k_edge(u: int, w: int, wt: int) -> int:
        idx = len(k_edges)
        k_edges.append((u, w, wt & 1))
        return idx

    for v in range(nd):
        vt = d.vtype[v]
        if vt in ("X0", "X1"):
            for s0, s2 in ((0, 2), (1, 3)):
                if dc[4 * v + s0] == comp_i:
                    u = edge_id_of_dart[4 * v + s0]
                    w = edge_id_of_dart[4 * v + s2]
                    is_und_j = (not _is_over(vt, s0)) and (dc[4 * v + ((s0 + 1) & 3)] == comp_j)
                    ce = add_k_edge(u, w, 1 if is_und_j else 0)
                    transition_cell[(4 * v + s0, 4 * v + s2)] = [ce]
                    transition_cell[(4 * v + s2, 4 * v + s0)] = [ce]
        else:
            if dc[4 * v] == comp_i:
                for s in range(4):
                    u = edge_id_of_dart[4 * v + s]
                    w = edge_id_of_dart[4 * v + ((s + 1) & 3)]
                    ce = add_k_edge(u, w, 0)
                    transition_cell[(4 * v + s, 4 * v + ((s + 1) & 3))] = [ce]
                    transition_cell[(4 * v + ((s + 1) & 3), 4 * v + s)] = [ce]

    spokes_minus: list[list[int]] = []
    for r, cyc in enumerate(cycs_minus):
        spokes_minus.append(
            [add_k_edge(edge_id_of_dart[d_in], v_minus_0 + r, pot_minus[d_in]) for d_in, _ in cyc]
        )

    spokes_plus: list[list[int]] = []
    for s_idx, cyc in enumerate(cycs_plus):
        spokes_plus.append(
            [add_k_edge(edge_id_of_dart[d_in], v_plus_0 + s_idx, pot_plus[d_in]) for d_in, _ in cyc]
        )

    num_ke = len(k_edges)
    two_cells: list[list[int]] = []
    for v in range(nd):
        if d.vtype[v] in ("Y", "Z") and dc[4 * v] == comp_i:
            two_cells.append([transition_cell[(4 * v + s, 4 * v + ((s + 1) & 3))][0] for s in range(4)])
    for r, cyc in enumerate(cycs_minus):
        sp_r = spokes_minus[r]
        m_len = len(cyc)
        for k in range(m_len):
            d_in, d_out = cyc[k]
            two_cells.append([sp_r[k]] + transition_cell[(d_in, d_out)] + [sp_r[(k + 1) % m_len]])
    for s_idx, cyc in enumerate(cycs_plus):
        sp_s = spokes_plus[s_idx]
        m_len = len(cyc)
        for k in range(m_len):
            d_in, d_out = cyc[k]
            two_cells.append([sp_s[k]] + transition_cell[(d_in, d_out)] + [sp_s[(k + 1) % m_len]])

    for fc in two_cells:
        w_sum = 0
        for ce in fc:
            w_sum ^= k_edges[ce][2]
        assert w_sum == 0, "1-cocycle boundary check failed on K_i!"

    parent_v = [-1] * num_kv
    parent_e = [-1] * num_kv
    parent_v[0] = 0
    q_bfs = [0]
    tree_edges = set()
    for u in q_bfs:
        for ei, (a, b, _) in enumerate(k_edges):
            nb = b if a == u else (a if b == u else -1)
            if nb != -1 and parent_v[nb] == -1:
                parent_v[nb] = u
                parent_e[nb] = ei
                tree_edges.add(ei)
                q_bfs.append(nb)

    def path_to_root_mask(v: int) -> int:
        mask = 0
        cur = v
        while cur != 0:
            mask ^= 1 << parent_e[cur]
            cur = parent_v[cur]
        return mask

    z1_cycles = [
        (1 << ei) ^ path_to_root_mask(a) ^ path_to_root_mask(b)
        for ei, (a, b, _) in enumerate(k_edges)
        if ei not in tree_edges
    ]

    def gf2_rref(vecs: list[int]) -> dict[int, int]:
        basis: dict[int, int] = {}
        for v in vecs:
            x = v
            while x:
                hb = x.bit_length() - 1
                if hb in basis:
                    x ^= basis[hb]
                else:
                    basis[hb] = x
                    break
        return basis

    def gf2_reduce(x: int, basis: dict[int, int]) -> int:
        while x:
            hb = x.bit_length() - 1
            if hb in basis:
                x ^= basis[hb]
            else:
                break
        return x

    b1_vecs = []
    for fc in two_cells:
        m_fc = 0
        for ce in fc:
            m_fc ^= 1 << ce
        b1_vecs.append(m_fc)
    b1_basis = gf2_rref(b1_vecs)

    h1_reps: list[int] = []
    aug_basis = dict(b1_basis)
    for cyc in z1_cycles:
        rem = gf2_reduce(cyc, aug_basis)
        if rem:
            aug_basis[rem.bit_length() - 1] = rem
            h1_reps.append(cyc)
    if len(h1_reps) != 2:
        return None

    def h1_coords(cycle_mask: int) -> tuple[int, int]:
        for c0 in (0, 1):
            for c1 in (0, 1):
                cand = cycle_mask ^ (h1_reps[0] if c0 else 0) ^ (h1_reps[1] if c1 else 0)
                if gf2_reduce(cand, b1_basis) == 0:
                    return (c0, c1)
        raise AssertionError("Cycle not in Z_1!")

    def eval_lambda(cycle_mask: int) -> int:
        tot = 0
        for ei in range(num_ke):
            if (cycle_mask >> ei) & 1:
                tot ^= k_edges[ei][2]
        return tot

    l0 = eval_lambda(h1_reps[0])
    l1 = eval_lambda(h1_reps[1])
    if (l0, l1) == (0, 0):
        return None

    ker_vec = next(
        (c0, c1)
        for c0 in (0, 1)
        for c1 in (0, 1)
        if (c0, c1) != (0, 0) and ((c0 * l0) ^ (c1 * l1)) == 0
    )

    marked_vs = [v for v in range(nd) if d.vtype[v] in ("Y", "Z") and dc[4 * v] == comp_i]
    self_class_vs = [
        v
        for v in range(nd)
        if d.vtype[v] in ("X0", "X1") and dc[4 * v] == comp_i and dc[4 * v + 1] == comp_i
    ]

    q_vals: dict[tuple[int, int], int] = {(0, 0): 0}
    subsets: list[set[int]] = []
    if len(marked_vs) <= 8:
        for r_sz in range(len(marked_vs) + 1):
            for comb in combinations(marked_vs, r_sz):
                subsets.append(set(comb))
    else:
        for r_sz in (0, 1, 2):
            for comb in combinations(marked_vs, r_sz):
                subsets.append(set(comb))
                subsets.append(set(marked_vs) - set(comb))

    for S_plus in subsets:
        seen_d = set()
        slice_cycs: list[list[tuple[int, int]]] = []
        dart_to_scyc: dict[int, int] = {}
        for hd in range(4 * nd):
            if dc[hd] != comp_i or hd in seen_d:
                continue
            sc_idx = len(slice_cycs)
            c = []
            cur = hd
            while cur not in seen_d:
                seen_d.add(cur)
                v, s = cur >> 2, cur & 3
                vt = d.vtype[v]
                if vt in ("X0", "X1"):
                    ex = (v << 2) | (s ^ 2)
                else:
                    sgn = 1 if v in S_plus else -1
                    p03 = (vt == "Y" and sgn > 0) or (vt == "Z" and sgn < 0)
                    ex = (v << 2) | ((3 - s) if p03 else (s ^ 1))
                seen_d.add(ex)
                dart_to_scyc[cur] = sc_idx
                dart_to_scyc[ex] = sc_idx
                c.append((cur, ex))
                cur = d.opp[ex]
            slice_cycs.append(c)

        cross_ab: dict[tuple[int, int], int] = {}
        for v in self_class_vs:
            ca = dart_to_scyc[4 * v]
            cb = dart_to_scyc[4 * v + 1]
            if ca != cb:
                key = (ca, cb) if ca < cb else (cb, ca)
                cross_ab[key] = cross_ab.get(key, 0) + 1

        for sc_idx, cyc in enumerate(slice_cycs):
            cmask = 0
            for d_in, d_out in cyc:
                for ce in transition_cell[(d_in, d_out)]:
                    cmask ^= 1 << ce
            coords = h1_coords(cmask)
            if coords == (0, 0) or coords in q_vals:
                continue

            c_self_vs = [
                v for v in self_class_vs if dart_to_scyc[4 * v] == sc_idx and dart_to_scyc[4 * v + 1] == sc_idx
            ]
            if len(c_self_vs) <= 2:
                arf_c = 0
            else:
                v2n = {v: idx for idx, v in enumerate(c_self_vs)}
                c_opp = [0] * (4 * len(c_self_vs))
                c_vt = tuple(d.vtype[v] for v in c_self_vs)
                for v in c_self_vs:
                    nv = v2n[v]
                    for s in range(4):
                        cur = d.opp[4 * v + s]
                        while (cur >> 2) not in v2n:
                            wv, ws = cur >> 2, cur & 3
                            wvt = d.vtype[wv]
                            if wvt in ("X0", "X1"):
                                ex = (wv << 2) | (ws ^ 2)
                            else:
                                sgn = 1 if wv in S_plus else -1
                                p03 = (wvt == "Y" and sgn > 0) or (wvt == "Z" and sgn < 0)
                                ex = (wv << 2) | ((3 - ws) if p03 else (ws ^ 1))
                            cur = d.opp[ex]
                        c_opp[4 * nv + s] = (v2n[cur >> 2] << 2) | (cur & 3)
                arf_c = classical_knot_arf(ClassicalLinkDiagram(len(c_self_vs), tuple(c_opp), c_vt, 0))

            lk_sum = 0
            for o_idx in range(len(slice_cycs)):
                if o_idx != sc_idx:
                    key = (sc_idx, o_idx) if sc_idx < o_idx else (o_idx, sc_idx)
                    lk_sum ^= (cross_ab.get(key, 0) // 2) & 1

            q_vals[coords] = (arf_c ^ lk_sum) & 1
            if len(q_vals) == 3:
                break
        if len(q_vals) >= 3:
            break

    nonzero = [(1, 0), (0, 1), (1, 1)]
    known_nz = [x for x in nonzero if x in q_vals]
    if len(known_nz) < 2:
        return None
    if len(known_nz) == 2:
        a, b = known_nz[0], known_nz[1]
        c_vec = (a[0] ^ b[0], a[1] ^ b[1])
        q_vals[c_vec] = q_vals[a] ^ q_vals[b] ^ 1
    return q_vals[ker_vec]


def double_linking_mod2(d: MarkedDiagram) -> tuple[tuple[str, str, int], ...]:
    """Return sorted `(type_i, type_j, Dlk_{ij})` across all component pairs `(i, j)` with `F_i ~ T^2`."""
    topo = d.surface_topology()
    c = topo.num_components
    if c < 2:
        return ()
    out: list[tuple[str, str, int]] = []
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
                out.append((ci.homeomorphism_type, f"{cj.homeomorphism_type}(e={cj.normal_euler})", val))
    out.sort()
    return tuple(out)
