"""Tier-3 purely 4-dimensional homotopy, branched-cover, symmetric quandle, and symmetry invariants.

Implements:
1. `double_branched_cover_invariants(d)`: Exact closed oriented 4-manifold invariants of the total
   mod-2 double branched cover `Sigma_2(S^4, F)` (`chi`, `b_1`, `H_1` torsion, `b_2`, `sigma`,
   `b_2^+`, `b_2^-`) via Fox-Plans homology and the Hirzebruch-Atiyah-Singer G-Signature Theorem.
2. `lomonaco_pi2_invariants(d)`: Exact regular finite-cover homology `(b_1, Tors H_1, b_2, b_3)` of
   `widetilde{E}_rho -> E = S^4 \\ nu(F)` via Schreier-Tietze lifting and Poincare-Lefschetz duality
   (Theorem 5.5 of `math/foundations_and_invariants.md`), together with certified `rk_Z pi_2(E)` on
   finite-`pi_1` and cyclic-kernel finite-commutator surface-knots (`0` on `0_1`, `1` on `2_1^-1`,
   `2` on `10_2`, `7` on `10_3`, `0` on `8_1^-1,-1`, `16` on `10_1^-2,-2`, `2` on `10_1^-1,-1`).
3. `octahedral_symmetric_quandle_invariants(d)`: Kamada-Oshiro octahedral symmetric quandle
   `(Q_6, *, rho)` coloring counts and 12-dimensional `H^3_{Q, rho}(Q_6; Z_2)` symmetric 3-cocycle
   state-sum invariant directly from constructive Reidemeister unlink movies `(M_-, M_+)`, applicable
   to ALL surface-links (orientable, non-orientable, and mixed) without choosing orientations.
4. `symmetry_profile(d)`: 4D chirality, orientation-sensitivity (`+F` vs `-F`), normal-Euler mirror
   obstruction, and constructive `SO(4)` Yoshikawa amphicheirality certificates.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.invariants.fundamental_group import (
    FiniteGroup,
    GroupPresentation,
    STANDARD_TEST_GROUPS,
    Word,
    _group_from_permutations,
    cyclically_reduce_word,
    finite_quotient_profile,
    reduce_word,
    tietze_simplify,
    wirtinger_presentation,
)
from surfaceknots.invariants.linking import _SL2_ID, _SL2_INV, _SL2_MUL
from surfaceknots.invariants.quandle import kei_module_invariants, mochizuki_3cocycle_invariant
from surfaceknots.moves import (
    canonical_marked_key,
    greedy_yoshikawa_reduce,
    list_level_yoshikawa_neighbors,
)
from surfaceknots.unlink import _is_over as cl_is_over, find_unlink_movie


# ============================================================================
# 1. Total Mod-2 Double Branched Cover Sigma_2(S^4, F)
# ============================================================================


@dataclass(frozen=True)
class DoubleBranchedCoverResult:
    """Exact 4-manifold invariants of the total mod-2 double branched cover `Sigma_2(S^4, F)`.

    Under an ambient orientation-reversing mirror `F -> F*` in `O(4) \\ SO(4)`, `signature` negates
    (`b_2^+ <-> b_2^-`); `unoriented_betti_pair = (min(b2_plus, b2_minus), max(b2_plus, b2_minus))`
    and `abs_signature = abs(signature)` are invariant under all of `O(4)`.
    """

    euler_characteristic: int
    b1: int
    h1_torsion: tuple[int, ...]
    b2: int
    signature: int
    abs_signature: int
    b2_plus: int
    b2_minus: int
    is_rational_homology_4sphere: bool

    def o4_key(self) -> tuple:
        """`O(4)`-invariant tuple suitable for unmirrored census stratification."""
        return (
            self.euler_characteristic,
            self.b1,
            self.h1_torsion,
            self.b2,
            self.abs_signature,
        )


def double_branched_cover_invariants(d: MarkedDiagram) -> DoubleBranchedCoverResult:
    """Compute exact homology and signature invariants of `Sigma_2(S^4, F(D))`."""
    topo = d.surface_topology()
    if topo.num_components == 0:
        # Branched cover over empty set: S^4 u S^4 (here treat base connected convention if nonempty)
        return DoubleBranchedCoverResult(
            euler_characteristic=4,
            b1=0,
            h1_torsion=(),
            b2=0,
            signature=0,
            abs_signature=0,
            b2_plus=0,
            b2_minus=0,
            is_rational_homology_4sphere=False,
        )

    r_kei, t_kei = kei_module_invariants(d)
    chi_f = topo.euler_characteristic
    total_normal_euler = sum(c.normal_euler for c in topo.components)
    assert total_normal_euler % 2 == 0, f"Odd total normal Euler number {total_normal_euler}!"

    chi_sigma2 = 4 - chi_f
    b1 = r_kei - 1
    b2 = 2 + 2 * b1 - chi_f
    sigma = -total_normal_euler // 2
    assert (b2 + sigma) % 2 == 0 and b2 >= abs(sigma), (
        f"Invalid 4-manifold Betti/signature: b2={b2}, sigma={sigma}"
    )
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


# ============================================================================
# 2. Lomonaco pi_2(E) & Poincare-Lefschetz Regular Finite Cover Homology
# ============================================================================


def _build_cover_target_groups() -> tuple[FiniteGroup, ...]:
    z2 = _group_from_permutations("Z2", [(1, 0)])
    z3 = _group_from_permutations("Z3", [(1, 2, 0)])
    sl2_order = 24
    perm_to_0 = list(range(sl2_order))
    if _SL2_ID != 0:
        perm_to_0[0], perm_to_0[_SL2_ID] = _SL2_ID, 0
    inv_p0 = {old: new for new, old in enumerate(perm_to_0)}
    sl2_mul = tuple(
        tuple(inv_p0[_SL2_MUL[perm_to_0[i]][perm_to_0[j]]] for j in range(sl2_order))
        for i in range(sl2_order)
    )
    sl2_inv = tuple(inv_p0[_SL2_INV[perm_to_0[i]]] for i in range(sl2_order))
    sl2_grp = FiniteGroup("SL2F3", 24, sl2_mul, sl2_inv)

    by_name = {g.name: g for g in STANDARD_TEST_GROUPS}
    return (
        z2,
        z3,
        by_name["S3"],
        by_name["Q8"],
        by_name["A4"],
        by_name["Dic16"],
        sl2_grp,
    )


_COVER_TARGET_GROUPS: tuple[FiniteGroup, ...] = _build_cover_target_groups()


@dataclass(frozen=True)
class CoverHomologyResult:
    """Regular finite cover Betti/torsion profiles `(b_1, Tors H_1, b_2, b_3)` and `pi_2(E)` data."""

    cover_betti_profile: tuple[tuple[str, tuple[tuple[int, tuple[int, ...], int, int], ...]], ...]
    certified_pi2_z_rank: int | None
    pi2_certification_reason: str | None

    def to_dict(self) -> dict[str, tuple]:
        return dict(self.cover_betti_profile)


def _extract_tietze_peripheral_system(
    d: MarkedDiagram,
) -> tuple[int, tuple[Word, ...], list[tuple[Word, list[Word]]]]:
    """Extract Tietze-simplified `pi_1(E)` presentation together with basepoint-consistent
    peripheral generating words `(mu_i, [sky_longitudes_i])` for EVERY component `F_i`
    (orientable or non-orientable).
    """
    n = d.num_vertices
    dc, num_scomps = d._dart_surface_components()
    if n == 0:
        periph0 = [((i + 1,), []) for i in range(d.free_loops)]
        return (d.free_loops, (), periph0)

    is_cl = [vt in ("X0", "X1") for vt in d.vtype]

    def is_und(v: int, s: int) -> bool:
        return ((s & 1) == 0) if d.vtype[v] == "X0" else ((s & 1) == 1)

    def mpart(v: int, s: int) -> int:
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

    rels: list[Word] = []
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

    raw_periph: list[tuple[int, list[Word]]] = []
    graph_scomps = num_scomps - d.free_loops
    for ci in range(graph_scomps):
        c_darts = [hd for hd in range(4 * n) if dc[hd] == ci]
        root = c_darts[0]
        mu_gen = dart_arc[root] + 1
        adj: dict[int, list[tuple[int, Word, tuple]]] = {hd: [] for hd in c_darts}
        for hd in c_darts:
            od = d.opp[hd]
            if hd < od:
                adj[hd].append((od, (), ("e", hd, od)))
                adj[od].append((hd, (), ("e", hd, od)))
            elif hd == od:
                adj[hd].append((od, (), ("e_loop", hd)))
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

        tree_word: dict[int, Word] = {root: ()}
        tree_edges: set[tuple] = set()
        q = deque([root])
        while q:
            u = q.popleft()
            for nb, w_edge, eid in adj[u]:
                if nb not in tree_word:
                    tree_word[nb] = tree_word[u] + w_edge
                    tree_edges.add(eid)
                    q.append(nb)

        cyc_words: list[Word] = []
        seen_non_tree: set[tuple] = set()
        for u in c_darts:
            for nb, w_edge, eid in adj[u]:
                if eid in tree_edges or eid in seen_non_tree:
                    continue
                seen_non_tree.add(eid)
                # CRITICAL: use `reduce_word` (preserving basepoint `root`!), NOT `cyclically_reduce_word`!
                cw = reduce_word(
                    tree_word[u] + w_edge + tuple(-x for x in reversed(tree_word[nb]))
                )
                if cw:
                    cyc_words.append(cw)
        raw_periph.append((mu_gen, cyc_words))

    for _ in range(d.free_loops):
        g += 1
        raw_periph.append((g, []))

    flat_words: list[Word] = []
    slices: list[tuple[int, int]] = []
    for mu_gen, cyc_words in raw_periph:
        start_idx = len(flat_words)
        flat_words.append((mu_gen,))
        flat_words.extend(cyc_words)
        slices.append((start_idx, len(flat_words)))

    simp = tietze_simplify(
        GroupPresentation(
            num_generators=g,
            relators=tuple(rels),
            component_meridians=tuple(flat_words),
        )
    )
    simp_periph: list[tuple[Word, list[Word]]] = []
    for a_idx, b_idx in slices:
        mu_w = simp.component_meridians[a_idx]
        cycs_w = [w for w in simp.component_meridians[a_idx + 1 : b_idx] if w]
        simp_periph.append((mu_w, cycs_w))

    return (simp.num_generators, simp.relators, simp_periph)


def _exact_rational_rank(rows: list[list[int]], ncols: int) -> int:
    """Compute exact rank over Q of an integer matrix via Sylvester-Bareiss fraction-free elimination."""
    mat = [list(r) for r in rows if any(r)]
    m = len(mat)
    if m == 0 or ncols == 0:
        return 0
    rank = 0
    prev_piv = 1
    for col in range(ncols):
        piv_row = -1
        for i in range(rank, m):
            if mat[i][col] != 0:
                piv_row = i
                break
        if piv_row == -1:
            continue
        if piv_row != rank:
            mat[rank], mat[piv_row] = mat[piv_row], mat[rank]
        row_r = mat[rank]
        piv = row_r[col]
        for i in range(rank + 1, m):
            row_i = mat[i]
            factor = row_i[col]
            if factor == 0:
                if piv != prev_piv:
                    for j in range(col + 1, ncols):
                        if row_i[j]:
                            row_i[j] = (piv * row_i[j]) // prev_piv
            else:
                row_i[col] = 0
                for j in range(col + 1, ncols):
                    row_i[j] = (piv * row_i[j] - factor * row_r[j]) // prev_piv
        prev_piv = piv
        rank += 1
        if rank == m:
            break
    return rank


def _cover_homology_for_surjection(
    num_gens: int,
    relators: tuple[Word, ...],
    periph: list[tuple[Word, list[Word]]],
    chi_f: int,
    grp: FiniteGroup,
    gen_images: tuple[int, ...],
) -> tuple[int, tuple[int, ...], int, int]:
    """Compute exact `(b_1, Tors H_1, b_2, b_3)` of the connected `d`-sheeted regular cover
    `widetilde{E}_rho -> E` via Schreier lifting and Poincare-Lefschetz duality (Theorem 5.5).
    """
    d = grp.order
    mul = grp.mul
    inv = grp.inv

    parent_edge: list[int] = [-1] * d
    visited = [False] * d
    visited[0] = True
    q = deque([0])
    while q:
        u = q.popleft()
        for j in range(num_gens):
            v = mul[u][gen_images[j]]
            if not visited[v]:
                visited[v] = True
                parent_edge[v] = u * num_gens + j
                q.append(v)

    tree_edge_set = {pe for pe in parent_edge if pe >= 0}
    edge_to_var = [-1] * (d * num_gens)
    n_vars = 0
    for e_idx in range(d * num_gens):
        if e_idx not in tree_edge_set:
            edge_to_var[e_idx] = n_vars
            n_vars += 1

    def trace_word_from_sheet(start_s: int, w: Word) -> tuple[int, list[int]]:
        cur = start_s
        vec = [0] * n_vars
        for lit in w:
            if lit > 0:
                j = lit - 1
                e_idx = cur * num_gens + j
                v_id = edge_to_var[e_idx]
                if v_id >= 0:
                    vec[v_id] += 1
                cur = mul[cur][gen_images[j]]
            else:
                j = (-lit) - 1
                nxt = mul[cur][inv[gen_images[j]]]
                e_idx = nxt * num_gens + j
                v_id = edge_to_var[e_idx]
                if v_id >= 0:
                    vec[v_id] -= 1
                cur = nxt
        return cur, vec

    rel_rows: list[list[int]] = []
    for r in relators:
        for s in range(d):
            end_s, vec = trace_word_from_sheet(s, r)
            assert end_s == s
            if any(vec):
                rel_rows.append(vec)

    if n_vars == 0:
        b1 = 0
        torsion: tuple[int, ...] = ()
        rk_r = 0
    elif not rel_rows:
        b1 = n_vars
        torsion = ()
        rk_r = 0
    else:
        mat_r = sp.Matrix(rel_rows)
        snf = smith_normal_form(mat_r, domain=sp.ZZ)
        diag = [abs(int(snf[i, i])) for i in range(min(snf.rows, snf.cols))]
        nonzero = [x for x in diag if x > 0]
        rk_r = len(nonzero)
        b1 = n_vars - rk_r
        torsion = tuple(sorted(x for x in nonzero if x > 1))

    b0_bdry = 0
    periph_cycle_rows: list[list[int]] = []
    for mu_w, cyc_words in periph:
        w_list: list[Word] = ([mu_w] if mu_w else []) + list(cyc_words)
        w_imgs = [trace_word_from_sheet(0, w)[0] for w in w_list]
        comp_vis = [False] * d
        for root_s in range(d):
            if comp_vis[root_s]:
                continue
            b0_bdry += 1
            comp_vis[root_s] = True
            path_word: dict[int, Word] = {root_s: ()}
            q_p = deque([root_s])
            comp_nodes = [root_s]
            tree_p_edges: set[tuple[int, int]] = set()
            while q_p:
                u = q_p.popleft()
                for w_idx, g_im in enumerate(w_imgs):
                    v = mul[u][g_im]
                    if not comp_vis[v]:
                        comp_vis[v] = True
                        comp_nodes.append(v)
                        path_word[v] = path_word[u] + w_list[w_idx]
                        tree_p_edges.add((u, w_idx))
                        q_p.append(v)
            for u in comp_nodes:
                for w_idx, g_im in enumerate(w_imgs):
                    if (u, w_idx) in tree_p_edges:
                        continue
                    v = mul[u][g_im]
                    loop_w = path_word[u] + w_list[w_idx] + tuple(-x for x in reversed(path_word[v]))
                    end_s, p_vec = trace_word_from_sheet(root_s, loop_w)
                    assert end_s == root_s
                    if any(p_vec):
                        periph_cycle_rows.append(p_vec)

    if b1 == 0 or not periph_cycle_rows:
        rk_im = 0
    else:
        rk_comb = _exact_rational_rank(rel_rows + periph_cycle_rows, n_vars)
        rk_im = rk_comb - rk_r

    b3 = (b0_bdry - 1) + b1 - rk_im
    b2 = d * (2 - chi_f) - 1 + b1 + b3
    assert b3 >= 0 and b2 >= 0, f"Invalid cover Betti numbers b2={b2}, b3={b3}!"
    return (b1, torsion, b2, b3)


# Known finite fundamental groups occurring in low-ch censuses (`pi_1(E) ~= Gamma` literally)
_CERTIFIED_FINITE_PI1_QUOTIENT_PROFILES: dict[dict[str, int] | tuple[tuple[str, int], ...], str] = {}


def lomonaco_pi2_invariants(
    d: MarkedDiagram,
    include_sl2f3: bool = True,
) -> CoverHomologyResult:
    """Compute exact regular finite cover homology `(b_1, Tors H_1, b_2, b_3)` across
    `(Z2, Z3, S3, Q8, A4, Dic16, SL2F3)` and infer `rk_Z pi_2(S^4 \\ F)` where certified.
    """
    topo = d.surface_topology()
    if topo.num_components == 0:
        return CoverHomologyResult(
            cover_betti_profile=(),
            certified_pi2_z_rank=0,
            pi2_certification_reason="empty_link",
        )

    num_gens, rels, periph = _extract_tietze_peripheral_system(d)
    chi_f = topo.euler_characteristic
    orientable_count = sum(1 for c in topo.components if c.orientable)

    profile_items: list[tuple[str, tuple[tuple[int, tuple[int, ...], int, int], ...]]] = []

    for grp in _COVER_TARGET_GROUPS:
        if grp.name == "SL2F3":
            # Condition strictly on topological invariant `topo.num_components` (`|pi_0(F)|`), never `free_loops`!
            if not include_sl2f3 or topo.num_components >= 3:
                continue

        order = grp.order
        mul = grp.mul
        inv = grp.inv

        adj_g: list[set[int]] = [set() for _ in range(num_gens)]
        for r in rels:
            vs_r = {abs(lit) - 1 for lit in r}
            for a in vs_r:
                adj_g[a] |= vs_r - {a}
        var_order: list[int] = []
        seen_v = [False] * num_gens
        for s_v in range(num_gens):
            if not seen_v[s_v]:
                seen_v[s_v] = True
                qv = deque([s_v])
                while qv:
                    u = qv.popleft()
                    var_order.append(u)
                    for wv in sorted(adj_g[u]):
                        if not seen_v[wv]:
                            seen_v[wv] = True
                            qv.append(wv)
        pos_of = {v_orig: step_i for step_i, v_orig in enumerate(var_order)}
        step_checks: list[list[Word]] = [[] for _ in range(num_gens)]
        forced_at: list[tuple[int, int] | None] = [None] * num_gens
        for r in rels:
            if len(r) == 2 and abs(r[0]) != abs(r[1]):
                v0, v1 = abs(r[0]) - 1, abs(r[1]) - 1
                s0, s1 = (1 if r[0] > 0 else -1), (1 if r[1] > 0 else -1)
                if pos_of[v1] > pos_of[v0] and forced_at[pos_of[v1]] is None:
                    forced_at[pos_of[v1]] = (v0, -s0 * s1)
                elif pos_of[v0] > pos_of[v1] and forced_at[pos_of[v0]] is None:
                    forced_at[pos_of[v0]] = (v1, -s0 * s1)
            mx = max(pos_of[abs(lit) - 1] for lit in r)
            step_checks[mx].append(r)

        conj_reps: list[int] = []
        seen_conj: set[int] = set()
        for x in range(order):
            if x in seen_conj:
                continue
            conj_reps.append(x)
            for h in range(order):
                seen_conj.add(mul[mul[h][x]][inv[h]])

        assign = [0] * num_gens
        surj_homs: list[tuple[int, ...]] = []

        def is_surjective(imgs: tuple[int, ...]) -> bool:
            seen_el = {0}
            st = [0]
            gens_el = set(imgs)
            while st:
                x = st.pop()
                for ge in gens_el:
                    y = mul[x][ge]
                    if y not in seen_el:
                        seen_el.add(y)
                        st.append(y)
            return len(seen_el) == order

        def dfs(step_i: int) -> None:
            if step_i == num_gens:
                tup = tuple(assign)
                if is_surjective(tup):
                    surj_homs.append(tup)
                return
            v_orig = var_order[step_i]
            fc = forced_at[step_i]
            if fc is not None:
                src_v, sgn = fc
                cand = (assign[src_v] if sgn == 1 else inv[assign[src_v]],)
            elif step_i == 0:
                cand = tuple(conj_reps)
            else:
                cand = tuple(range(order))

            for val in cand:
                assign[v_orig] = val
                ok = True
                for r in step_checks[step_i]:
                    acc = 0
                    for lit in r:
                        gv = assign[lit - 1] if lit > 0 else inv[assign[-lit - 1]]
                        acc = mul[acc][gv]
                    if acc != 0:
                        ok = False
                        break
                if ok:
                    dfs(step_i + 1)

        if num_gens > 0:
            dfs(0)

        if surj_homs:
            # Deduplicate by exact kernel (`Aut(grp)`-orbit) via canonical BFS Cayley signature
            uniq_homs: list[tuple[int, ...]] = []
            seen_sigs: set[tuple[int, ...]] = set()
            for hom in surj_homs:
                remap = [-1] * order
                remap[0] = 0
                bfs_order = [0]
                nxt_id = 1
                sig_list: list[int] = []
                for u in bfs_order:
                    for g_im in hom:
                        v = mul[u][g_im]
                        if remap[v] == -1:
                            remap[v] = nxt_id
                            nxt_id += 1
                            bfs_order.append(v)
                        sig_list.append(remap[v])
                sig = tuple(sig_list)
                if sig not in seen_sigs:
                    seen_sigs.add(sig)
                    uniq_homs.append(hom)

            bettis = tuple(
                sorted(
                    {
                        _cover_homology_for_surjection(num_gens, rels, periph, chi_f, grp, hom)
                        for hom in uniq_homs
                    }
                )
            )
            profile_items.append((grp.name, bettis))

    prof_dict = dict(profile_items)

    # Certified `rk_Z pi_2(E)` detection where presentation / cover is exact universal or cyclic-universal:
    pi2_rk: int | None = None
    reason: str | None = None

    if num_gens == 1 and not rels:
        # Exact infinite cyclic fundamental group `pi_1(E) == Z`!
        # Wait: when `pi_1(E) == Z` and `D` is an unknotted orientable surface `Sigma_g` (`c=1`):
        # For `S^2` (`chi = 2`), `E ~= S^1` -> `pi_2(E) = 0` (`rk_Z = 0`).
        # For `g >= 1` (`chi = 2 - 2g < 2`), `pi_2(E) ~= Z[t^{+-1}]^{2g}` has infinite Z-rank!
        if chi_f == 2 and topo.num_components == 1:
            pi2_rk = 0
            reason = "exact_pi1_Z_sphere"
    elif num_gens == 1 and len(rels) == 1 and abs(rels[0][0]) == 1 and len(set(rels[0])) == 1:
        # Cyclic finite fundamental group `Z_m` (`m = len(rels[0])`, e.g. `Z_2` on `2_1^-1`):
        m_ord = len(rels[0])
        if m_ord == 2 and "Z2" in prof_dict and len(prof_dict["Z2"]) == 1:
            b1_u, tors_u, b2_u, _ = prof_dict["Z2"][0]
            if b1_u == 0 and not tors_u:
                pi2_rk = b2_u
                reason = "universal_cover_Z2"
    elif orientable_count == 0:
        # Check if one of the non-abelian finite groups `Q8` or `Dic16` is the universal cover
        # (`b1 == 0, tors == ()` on a 2-generator finite presentation known to have order 8 or 16):
        # Verify finite order directly by Todd-Coxeter coset enumeration (with small cap `<= 64`)!
        tc_order = _todd_coxeter_group_order(num_gens, rels, max_cosets=64)
        if tc_order is not None:
            for gname, g_ord in (("Z2", 2), ("S3", 6), ("Q8", 8), ("A4", 12), ("Dic16", 16), ("SL2F3", 24)):
                if tc_order == g_ord and gname in prof_dict and len(prof_dict[gname]) == 1:
                    b1_u, tors_u, b2_u, _ = prof_dict[gname][0]
                    if b1_u == 0 and not tors_u:
                        pi2_rk = b2_u
                        reason = f"universal_cover_{gname}"
                        break
    elif topo.num_components == 1 and chi_f == 2 and orientable_count == 1:
        # Connected 2-sphere knot `S^2`: test if commutator subgroup `[pi_1, pi_1]` is finite by
        # killing `k * mu_0` (`mu_0^k = 1` for `k in (2, 3)`) via Todd-Coxeter!
        mu_w = periph[0][0]
        for k_mer, expected_g in ((1, None), (2, "S3"), (3, "SL2F3")):
            extra_rel = mu_w * k_mer
            q_ord = _todd_coxeter_group_order(num_gens, rels + (extra_rel,), max_cosets=64)
            if k_mer == 1 and q_ord == 1:
                # `pi_1(E) / <mu_0> = 1` always; check if `pi_1` itself is `Z`:
                pass
            elif expected_g is not None and q_ord is not None and expected_g in prof_dict:
                grp_expected_ord = 6 if expected_g == "S3" else 24
                if q_ord == grp_expected_ord:
                    # `pi_1 / <mu_0^k> ~= Gamma` of exact order `q_ord`! Wait: check if `ker(pi_1 -> Gamma)` is `Z`!
                    # Verify by checking `|pi_1 / <mu_0^{2*k_mer}>| == 2 * q_ord` AND `b1 == 1, tors == ()`!
                    q_ord2 = _todd_coxeter_group_order(num_gens, rels + (mu_w * (2 * k_mer),), max_cosets=128)
                    if q_ord2 == 2 * grp_expected_ord:
                        cands = [b2_c for (b1_c, t_c, b2_c, _) in prof_dict[expected_g] if b1_c == 1 and not t_c]
                        if cands:
                            pi2_rk = max(cands)
                            reason = f"finite_commutator_cover_{expected_g}"
                            break

    return CoverHomologyResult(
        cover_betti_profile=tuple(profile_items),
        certified_pi2_z_rank=pi2_rk,
        pi2_certification_reason=reason,
    )


def _todd_coxeter_group_order(
    num_gens: int,
    relators: tuple[Word, ...],
    max_cosets: int = 128,
) -> int | None:
    """Exact Todd-Coxeter coset enumeration over the trivial subgroup (returns `|G|` if `<= max_cosets`,
    else `None`).
    """
    if num_gens == 0:
        return 1
    if not relators:
        return None

    # Columns `0..2*num_gens-1`: `2*j` is `+g_{j+1}`, `2*j + 1` is `-g_{j+1}`.
    ncols = 2 * num_gens
    inv_col = [c ^ 1 for c in range(ncols)]
    table: list[list[int]] = [[-1] * ncols]
    parent = [0]
    live_count = 1
    coinc_q: deque[tuple[int, int]] = deque()

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def merge(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            coinc_q.append((ra, rb))

    def process_coincidences() -> None:
        nonlocal live_count
        while coinc_q:
            a, b = coinc_q.popleft()
            ra, rb = find(a), find(b)
            if ra == rb:
                continue
            if ra > rb:
                ra, rb = rb, ra
            parent[rb] = ra
            live_count -= 1
            for c in range(ncols):
                nb_b = table[rb][c]
                if nb_b != -1:
                    table[rb][c] = -1
                    nb_b = find(nb_b)
                    if table[nb_b][inv_col[c]] == rb:
                        table[nb_b][inv_col[c]] = -1
                    # Link `ra --c--> nb_b`
                    cur_a = table[ra][c]
                    if cur_a != -1:
                        merge(cur_a, nb_b)
                    else:
                        cur_inv = table[nb_b][inv_col[c]]
                        if cur_inv != -1:
                            merge(ra, cur_inv)
                        else:
                            table[ra][c] = nb_b
                            table[nb_b][inv_col[c]] = ra

    def define(x: int, c: int) -> bool:
        nonlocal live_count
        rx = find(x)
        if table[rx][c] != -1:
            return True
        if len(table) >= max_cosets * 4 or live_count >= max_cosets:
            return False
        ny = len(table)
        table.append([-1] * ncols)
        parent.append(ny)
        live_count += 1
        table[rx][c] = ny
        table[ny][inv_col[c]] = rx
        return True

    rel_cols = [
        tuple((2 * (lit - 1)) if lit > 0 else (2 * (-lit - 1) + 1) for lit in r)
        for r in relators
        if r
    ]

    i = 0
    while i < len(table):
        if find(i) != i:
            i += 1
            continue
        # Scan all relators from coset `i`
        for rc in rel_cols:
            if find(i) != i:
                break
            f_pos, b_pos = 0, len(rc) - 1
            f_cos, b_cos = i, i
            while f_pos <= b_pos:
                f_cos = find(f_cos)
                nxt = table[f_cos][rc[f_pos]]
                if nxt != -1:
                    f_cos = nxt
                    f_pos += 1
                    continue
                b_cos = find(b_cos)
                prv = table[b_cos][inv_col[rc[b_pos]]]
                if prv != -1:
                    b_cos = prv
                    b_pos -= 1
                    continue
                if f_pos == b_pos:
                    # Deduction `f_cos --rc[f_pos]--> b_cos`!
                    col_ded = rc[f_pos]
                    if table[f_cos][col_ded] != -1:
                        merge(table[f_cos][col_ded], b_cos)
                    elif table[b_cos][inv_col[col_ded]] != -1:
                        merge(f_cos, table[b_cos][inv_col[col_ded]])
                    else:
                        table[f_cos][col_ded] = b_cos
                        table[b_cos][inv_col[col_ded]] = f_cos
                    process_coincidences()
                    break
                if not define(f_cos, rc[f_pos]):
                    return None
                process_coincidences()
            else:
                if find(f_cos) != find(b_cos):
                    merge(f_cos, b_cos)
                    process_coincidences()

        if find(i) != i:
            i += 1
            continue
        for c in range(ncols):
            if find(i) != i:
                break
            if table[i][c] == -1:
                if not define(i, c):
                    return None
                process_coincidences()
        i += 1

    return live_count


# ============================================================================
# 3. Kamada-Oshiro Octahedral Symmetric Quandle (Q_6, *, rho) & 3-Cocycles
# ============================================================================


def _q6_rho(x: int) -> int:
    return x ^ 1


def _q6_op(x: int, y: int) -> int:
    return x if (x >> 1) == (y >> 1) else (x ^ 1)


def _precompute_q6_symmetric_3cocycles_mod2() -> tuple[tuple[tuple[tuple[int, ...], ...], ...], ...]:
    """Compute canonical reduced-echelon basis of `H^3_{Q, rho}(Q_6; Z_2)` (`dim = 12`)."""

    def idx3(x: int, y: int, z: int) -> int:
        return 36 * x + 6 * y + z

    rows_set: set[int] = set()

    def add_eq(vars_list: list[int]) -> None:
        mask = 0
        for v in vars_list:
            mask ^= 1 << v
        if mask:
            rows_set.add(mask)

    for x in range(6):
        for y in range(6):
            add_eq([idx3(x, x, y)])
            add_eq([idx3(x, y, y)])
            for z in range(6):
                add_eq([idx3(_q6_rho(x), y, z), idx3(x, y, z)])
                add_eq([idx3(x, _q6_rho(y), z), idx3(_q6_op(x, y), y, z)])
                add_eq([idx3(x, y, _q6_rho(z)), idx3(_q6_op(x, z), _q6_op(y, z), z)])
                for w in range(6):
                    add_eq(
                        [
                            idx3(x, y, z),
                            idx3(x, z, w),
                            idx3(x, y, w),
                            idx3(_q6_op(x, z), _q6_op(y, z), w),
                            idx3(_q6_op(x, y), z, w),
                            idx3(_q6_op(x, w), _q6_op(y, w), _q6_op(z, w)),
                        ]
                    )

    piv: dict[int, int] = {}
    for mask in sorted(rows_set, reverse=True):
        m = mask
        while m:
            hb = m.bit_length() - 1
            if hb in piv:
                m ^= piv[hb]
            else:
                piv[hb] = m
                break
    for hb in sorted(piv.keys()):
        for hb2 in sorted(piv.keys()):
            if hb2 > hb and ((piv[hb2] >> hb) & 1):
                piv[hb2] ^= piv[hb]

    free_vars = [v for v in range(216) if v not in piv]
    z3_basis: list[int] = []
    for fv in free_vars:
        vec = 1 << fv
        for hb, rmask in piv.items():
            if (rmask >> fv) & 1:
                vec |= 1 << hb
        z3_basis.append(vec)

    def idx2(x: int, y: int) -> int:
        return 6 * x + y

    c2_piv: dict[int, int] = {}
    for x in range(6):
        for y in range(6):
            for eq in (
                [idx2(x, x)],
                [idx2(_q6_rho(x), y), idx2(x, y)],
                [idx2(x, _q6_rho(y)), idx2(_q6_op(x, y), y)],
            ):
                m = 0
                for v in eq:
                    m ^= 1 << v
                while m:
                    hb = m.bit_length() - 1
                    if hb in c2_piv:
                        m ^= c2_piv[hb]
                    else:
                        c2_piv[hb] = m
                        break
    for hb in sorted(c2_piv.keys()):
        for hb2 in sorted(c2_piv.keys()):
            if hb2 > hb and ((c2_piv[hb2] >> hb) & 1):
                c2_piv[hb2] ^= c2_piv[hb]

    c2_free = [v for v in range(36) if v not in c2_piv]
    b3_piv: dict[int, int] = {}
    for fv in c2_free:
        phi = 1 << fv
        for hb, rmask in c2_piv.items():
            if (rmask >> fv) & 1:
                phi |= 1 << hb
        dphi = 0
        for x in range(6):
            for y in range(6):
                for z in range(6):
                    val = (
                        ((phi >> idx2(x, y)) & 1)
                        ^ ((phi >> idx2(_q6_op(x, z), _q6_op(y, z))) & 1)
                        ^ ((phi >> idx2(x, z)) & 1)
                        ^ ((phi >> idx2(_q6_op(x, y), z)) & 1)
                    )
                    if val:
                        dphi |= 1 << idx3(x, y, z)
        m = dphi
        while m:
            hb = m.bit_length() - 1
            if hb in b3_piv:
                m ^= b3_piv[hb]
            else:
                b3_piv[hb] = m
                break

    h3_vecs: list[int] = []
    cur_piv = dict(b3_piv)
    for zvec in z3_basis:
        m = zvec
        while m:
            hb = m.bit_length() - 1
            if hb in cur_piv:
                m ^= cur_piv[hb]
            else:
                cur_piv[hb] = m
                h3_vecs.append(zvec)
                break

    tensors = tuple(
        tuple(
            tuple(tuple((vec >> idx3(x, y, z)) & 1 for z in range(6)) for y in range(6))
            for x in range(6)
        )
        for vec in h3_vecs
    )
    return tensors


_Q6_H3_BASIS: tuple[tuple[tuple[tuple[int, ...], ...], ...], ...] = (
    _precompute_q6_symmetric_3cocycles_mod2()
)


@dataclass(frozen=True)
class SymmetricQuandleResult:
    """Kamada-Oshiro octahedral symmetric quandle `(Q_6, *, rho)` coloring & 3-cocycle invariants."""

    col_q6: int
    cocycle_weight_distribution: tuple[tuple[int, int], ...]
    trivial_weight_count: int
    nontrivial_weight_count: int
    certifies_positive_triple_points: bool

    def summary_str(self) -> str:
        if self.col_q6 == 0:
            return "0"
        if self.nontrivial_weight_count == 0:
            return str(self.col_q6)
        return f"{self.trivial_weight_count} + {self.nontrivial_weight_count}*u"


def _enumerate_q6_graph_colorings(d: MarkedDiagram) -> list[list[int]]:
    n = d.num_vertices
    if n == 0:
        return [[]]

    comp_id = [-1] * (4 * n)
    comp_flip = [0] * (4 * n)
    adj: list[list[tuple[int, int]]] = [[] for _ in range(4 * n)]

    def add_rel(a: int, b: int, flip: int) -> None:
        adj[a].append((b, flip))
        adj[b].append((a, flip))

    for hd in range(4 * n):
        add_rel(hd, d.opp[hd], 1)

    classical_vs: list[int] = []
    for v in range(n):
        vt = d.vtype[v]
        if vt in ("Y", "Z"):
            add_rel(4 * v + 0, 4 * v + 1, 1)
            add_rel(4 * v + 0, 4 * v + 2, 0)
            add_rel(4 * v + 0, 4 * v + 3, 1)
        else:
            classical_vs.append(v)
            if vt == "X0":
                add_rel(4 * v + 1, 4 * v + 3, 1)
            else:
                add_rel(4 * v + 0, 4 * v + 2, 1)

    num_sheets = 0
    for start in range(4 * n):
        if comp_id[start] != -1:
            continue
        cid = num_sheets
        num_sheets += 1
        comp_id[start] = cid
        comp_flip[start] = 0
        st = [start]
        while st:
            u = st.pop()
            for w, fl in adj[u]:
                want_fl = comp_flip[u] ^ fl
                if comp_id[w] == -1:
                    comp_id[w] = cid
                    comp_flip[w] = want_fl
                    st.append(w)
                elif comp_id[w] != cid or comp_flip[w] != want_fl:
                    return []

    checks_at: list[list[tuple[int, int, int, int, int, int]]] = [[] for _ in range(num_sheets)]
    for v in classical_vs:
        vt = d.vtype[v]
        if vt == "X0":
            du0, du2, do0 = 4 * v + 0, 4 * v + 2, 4 * v + 1
        else:
            du0, du2, do0 = 4 * v + 1, 4 * v + 3, 4 * v + 0
        chk = (
            comp_id[du0],
            comp_flip[du0],
            comp_id[du2],
            comp_flip[du2],
            comp_id[do0],
            comp_flip[do0],
        )
        checks_at[max(chk[0], chk[2], chk[4])].append(chk)

    assign = [-1] * num_sheets
    out: list[list[int]] = []

    def dfs(idx: int) -> None:
        if idx == num_sheets:
            out.append([assign[comp_id[hd]] ^ comp_flip[hd] for hd in range(4 * n)])
            return
        for val in range(6):
            assign[idx] = val
            ok = True
            for cu0, fu0, cu2, fu2, co0, fo0 in checks_at[idx]:
                if (assign[cu2] ^ fu2) != _q6_op(_q6_rho(assign[cu0] ^ fu0), assign[co0] ^ fo0):
                    ok = False
                    break
            if ok:
                dfs(idx + 1)
        assign[idx] = -1

    dfs(0)
    return out


def _eval_q6_movie_bitmask(movie: list, init_col: list[int]) -> int:
    """Propagate a `(Q_6, *, rho)` coloring along `movie` and return the 12-bit mod-2 cocycle vector."""
    col = list(init_col)
    total_mask = 0

    for step in movie:
        kind = step.kind
        aft = step.after
        meta = step.meta

        if kind in ("R1-", "R2-"):
            old_to_new = meta["old_to_new"]
            m = aft.num_crossings
            new_col = [0] * (4 * m)
            for vo, vn in old_to_new.items():
                for s in range(4):
                    new_col[4 * vn + s] = col[4 * vo + s]
            col = new_col

        elif kind == "R2+":
            da, oda, db, odb, vt = meta["da"], meta["oda"], meta["db"], meta["odb"], meta["vt"]
            n_prev = step.before.num_crossings
            u, w = n_prev, n_prev + 1
            new_col = list(col) + [0] * 8
            new_col[4 * u + 0] = _q6_rho(col[da])
            new_col[4 * w + 0] = _q6_rho(col[oda])
            new_col[4 * w + 1] = _q6_rho(col[db])
            new_col[4 * u + 3] = _q6_rho(col[odb])

            if vt == "X1":
                new_col[4 * u + 2] = _q6_rho(new_col[4 * u + 0])
                new_col[4 * w + 2] = _q6_rho(new_col[4 * w + 0])
                new_col[4 * u + 1] = _q6_op(_q6_rho(new_col[4 * u + 3]), new_col[4 * u + 0])
                new_col[4 * w + 3] = _q6_rho(new_col[4 * u + 1])
            else:
                new_col[4 * u + 1] = _q6_rho(new_col[4 * u + 3])
                new_col[4 * w + 3] = _q6_rho(new_col[4 * w + 1])
                new_col[4 * u + 2] = _q6_op(_q6_rho(new_col[4 * u + 0]), new_col[4 * u + 1])
                new_col[4 * w + 2] = _q6_rho(new_col[4 * u + 2])
            col = new_col

        elif kind == "R3":
            vs = meta["vs"]
            ss = meta["ss"]
            heights = meta["heights"]
            ports = meta["ports"]
            port_pairs = meta["port_pairs"]

            idx_t = heights.index(2)
            idx_m = heights.index(1)
            idx_b = heights.index(0)

            d_tri = [(vs[i] << 2) | ss[i] for i in range(3)]
            s_in_prev = [step.before.opp[d_tri[(i - 1) % 3]] & 3 for i in range(3)]
            tri_on_left = [(s_in_prev[i] == ((ss[i] + 1) & 3)) for i in range(3)]
            c_out_tri = [
                _q6_rho(col[d_tri[i]]) if tri_on_left[i] else col[d_tri[i]]
                for i in range(3)
            ]
            x = c_out_tri[idx_b]
            y = c_out_tri[idx_m]
            z = c_out_tri[idx_t]

            step_mask = 0
            for b_idx, th in enumerate(_Q6_H3_BASIS):
                if th[x][y][z]:
                    step_mask |= 1 << b_idx
            total_mask ^= step_mask

            new_col = list(col)
            for i in range(3):
                wi = vs[i]
                pi_a, pi_b = port_pairs[i]
                new_col[4 * wi + 2] = col[ports[pi_a]]
                new_col[4 * wi + 3] = col[ports[pi_b]]
                vt_wi = aft.vtype[wi]
                if cl_is_over(vt_wi, 2):
                    new_col[4 * wi + 0] = _q6_rho(new_col[4 * wi + 2])
                    new_col[4 * wi + 1] = _q6_op(_q6_rho(new_col[4 * wi + 3]), new_col[4 * wi + 2])
                else:
                    new_col[4 * wi + 1] = _q6_rho(new_col[4 * wi + 3])
                    new_col[4 * wi + 0] = _q6_op(_q6_rho(new_col[4 * wi + 2]), new_col[4 * wi + 3])
            col = new_col

        else:
            raise ValueError(f"Unsupported Reidemeister movie step {kind!r}")

        for hd in range(4 * aft.num_crossings):
            od = aft.opp[hd]
            assert col[od] == _q6_rho(col[hd]), f"Q6 symmetric edge mismatch after {kind}!"

    return total_mask


def octahedral_symmetric_quandle_invariants(d: MarkedDiagram) -> SymmetricQuandleResult:
    """Compute `(Q_6, *, rho)` coloring count and `H^3_{Q, rho}(Q_6; Z_2)` state-sum invariant."""
    cols = _enumerate_q6_graph_colorings(d)
    free_mult = 6 ** d.free_loops
    n_cols = len(cols) * free_mult
    if n_cols == 0:
        return SymmetricQuandleResult(
            col_q6=0,
            cocycle_weight_distribution=(),
            trivial_weight_count=0,
            nontrivial_weight_count=0,
            certifies_positive_triple_points=False,
        )

    classical_vs = [v for v in range(d.num_vertices) if d.vtype[v] in ("X0", "X1")]
    if not classical_vs:
        return SymmetricQuandleResult(
            col_q6=n_cols,
            cocycle_weight_distribution=((0, n_cols),),
            trivial_weight_count=n_cols,
            nontrivial_weight_count=0,
            certifies_positive_triple_points=False,
        )

    mp = find_unlink_movie(d.resolve(1))
    mm = find_unlink_movie(d.resolve(-1))
    if mp is None or mm is None:
        raise ValueError("Cannot evaluate symmetric quandle 3-cocycle on non-admissible diagram")

    mask_counts: Counter[int] = Counter()
    for col_d in cols:
        init_c = [col_d[4 * vo + s] for vo in classical_vs for s in range(4)]
        net_mask = _eval_q6_movie_bitmask(mp, init_c) ^ _eval_q6_movie_bitmask(mm, init_c)
        mask_counts[net_mask] += free_mult

    triv = mask_counts.get(0, 0)
    nontriv = n_cols - triv
    return SymmetricQuandleResult(
        col_q6=n_cols,
        cocycle_weight_distribution=tuple(sorted(mask_counts.items())),
        trivial_weight_count=triv,
        nontrivial_weight_count=nontriv,
        certifies_positive_triple_points=(nontriv > 0),
    )


# ============================================================================
# 4. 4D Symmetry Profile (Chirality, Invertibility, Amphicheirality)
# ============================================================================


@dataclass(frozen=True)
class SymmetryProfile:
    """Rigorous 4D symmetry classification of `F(D)` under mirror `F*` and orientation `-F`.

    Note: `invariant_symmetry_flags` depends strictly on topological invariants (safe for census
    stratification keys), while `constructive_so4_amphicheiral` records whether bounded Yoshikawa
    search explicitly connected `D` to `D.mirror_3d()` in `SO(4)`.
    """

    normal_euler_tuple: tuple[int, ...]
    chiral_by_normal_euler: bool
    oriented_chiral_by_cocycle: bool
    non_invertible_by_cocycle: bool
    constructive_so4_amphicheiral: bool
    symmetry_summary: str

    def invariant_flags(self) -> tuple[bool, bool, bool]:
        return (
            self.chiral_by_normal_euler,
            self.oriented_chiral_by_cocycle,
            self.non_invertible_by_cocycle,
        )


def _are_so4_move_connected(d1: MarkedDiagram, d2: MarkedDiagram, max_states: int = 250) -> bool:
    target_k = canonical_marked_key(d2, up_to_mirror_and_reverse=False)
    start_k = canonical_marked_key(d1, up_to_mirror_and_reverse=False)
    if start_k == target_k:
        return True
    vis = {start_k}
    q = deque([d1])
    while q and len(vis) < max_states:
        cur = q.popleft()
        for nb in list_level_yoshikawa_neighbors(cur):
            nb_r = greedy_yoshikawa_reduce(nb)
            k = canonical_marked_key(nb_r, up_to_mirror_and_reverse=False)
            if k == target_k:
                return True
            if k not in vis:
                vis.add(k)
                q.append(nb_r)
    return False


def symmetry_profile(d: MarkedDiagram, max_bfs_states: int = 250) -> SymmetryProfile:
    """Compute algebraic chirality/invertibility obstructions and constructive amphicheirality."""
    topo = d.surface_topology()
    e_pos = tuple(sorted(c.normal_euler for c in topo.components))
    e_mir = tuple(sorted(-c.normal_euler for c in topo.components))
    chiral_euler = e_pos != e_mir

    ori_chiral_cocycle = False
    non_inv_cocycle = False
    if topo.orientability_label == "Ori" and topo.num_components >= 1:
        primes_to_test = (3, 5) if topo.num_components <= 2 else (3,)
        for p in primes_to_test:
            tup = mochizuki_3cocycle_invariant(d, p=p, up_to_mirror=False)
            if tup is None:
                continue
            if len(set(tup)) > 1:
                non_inv_cocycle = True
            # Check single-orientation mirror conjugate (`u^k -> u^{p-k}`):
            # If any orientation state-sum `s` in `tup` is not equal to its own conjugate `bar(s)`,
            # then `(+F)*` has `Phi((+F)*) = bar(Phi(+F)) != Phi(+F)`, so oriented `+F` is chiral!
            # And if the multiset `tup` != `tuple(sorted(bar(s) for s in tup))`, then `F` is EVEN
            # UNORIENTED CHIRAL (`F != F*` ignoring surface orientation)!
            for s_val in tup:
                if "u" in s_val:
                    # Check if `s_val` is invariant under `u^k <-> u^{p-k}`
                    parts = [pt.strip() for pt in s_val.split("+")]
                    coeff_map: dict[int, int] = {}
                    for pt in parts:
                        if "u" not in pt:
                            coeff_map[0] = int(pt)
                        else:
                            if "*" in pt:
                                c_str, u_part = pt.split("*", 1)
                                c_int = int(c_str)
                            else:
                                c_int, u_part = 1, pt
                            exp_int = 1 if u_part == "u" else int(u_part.split("^")[1])
                            coeff_map[exp_int] = c_int
                    conj_map = {((-k) % p): v for k, v in coeff_map.items()}
                    if conj_map != coeff_map:
                        ori_chiral_cocycle = True

    so4_amp = False
    if not chiral_euler and max_bfs_states > 0:
        so4_amp = _are_so4_move_connected(d, d.mirror_3d(), max_states=max_bfs_states)

    if chiral_euler:
        summary = "Chiral (normal Euler e != -e)"
    elif ori_chiral_cocycle and non_inv_cocycle:
        if so4_amp:
            summary = "Minus-amphicheiral, non-invertible (+F != -F, (+F)* ~= -F)"
        else:
            summary = "Oriented-chiral & non-invertible (Mochizuki 3-cocycle)"
    elif so4_amp:
        summary = "Unoriented-amphicheiral (constructive SO(4) Yoshikawa path)"
    else:
        summary = "Symmetry open / no obstruction detected"

    return SymmetryProfile(
        normal_euler_tuple=e_pos,
        chiral_by_normal_euler=chiral_euler,
        oriented_chiral_by_cocycle=ori_chiral_cocycle,
        non_invertible_by_cocycle=non_inv_cocycle,
        constructive_so4_amphicheiral=so4_amp,
        symmetry_summary=summary,
    )
