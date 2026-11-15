"""Kei, Universal Linear Bikei, Oriented Alexander Biquandle, and 3-Cocycle Invariants.

Implements:
1. `kei_module_invariants(d)`: Universal involutory quandle (kei) coloring module invariant factors
   over `Z`, giving exact dihedral coloring counts `|Col_{R_n}(F)|` for all `n >= 2`.
2. `bikei_module_invariants(d)`: Universal linear involutory biquandle (bikei) invariant over
   `Z[s] / (s^2 - 1)` directly from `MarkedDiagram` semi-edges (works on ALL surface-links,
   including non-orientable and mixed multi-component surface-links, with zero orientation choice!).
3. `alexander_biquandle_invariants(d)`: Oriented Alexander biquandle module over `Z[t^{+-1}, s^{+-1}]`
   (Ashihara 2012; Kim-Joung-Lee 2015) across all relative component orientations, extracting both:
   - exact specialized Smith Normal Form invariant factors over `Z[t, s] / (t^a - 1, s^b - 1)` for
     `(a, b) in {(1, 2), (1, 3), (2, 2), (2, 3), (3, 2)}` (separating surface-links with identical
     `pi_1` and classical Alexander ideals such as `8_1^{1,1}` vs `10_1^{1,1}`!), and
   - exact symbolic first elementary ideal `E_1^{BQ}(t, s)` when `t = 1` (`s`-monodromy ideal) and
     when `t = s`.
4. `mochizuki_3cocycle_invariant(d, p)`: Quandle 3-cocycle state-sum invariant `Phi_{theta_p}(F)`
   over dihedral quandle `R_p` (`p` odd prime) computed directly from constructive Reidemeister
   unlinking movies of `L_+(D)` and `L_-(D)` (Carter-Jelsovsky-Kamada-Langford-Saito 2003;
   Kamada-Kim-Lee 2015), detecting 4D chirality (e.g. `10_2`, `10_3`) and bounding triple-point
   number `t(F)` from below!
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import product
from math import gcd

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

from surfaceknots.diagrams import ClassicalLinkDiagram, MarkedDiagram
from surfaceknots.invariants.fundamental_group import GroupPresentation, wirtinger_presentation


@dataclass(frozen=True)
class BiquandleResult:
    """Summary of kei, universal linear bikei, and Alexander biquandle invariants."""

    kei_free_rank: int
    kei_torsion: tuple[int, ...]
    col_r3: int
    col_r4: int
    col_r5: int
    bikei_free_rank: int
    bikei_torsion: tuple[int, ...]
    oriented_abq_snf_profile: tuple | None
    oriented_s_monodromy_poly: tuple[str, ...] | None
    cocycle_r3: tuple[str, ...] | None = None
    cocycle_r5: tuple[str, ...] | None = None


def _snf_rank_and_torsion(rows: list[list[int]], num_cols: int) -> tuple[int, tuple[int, ...]]:
    """Compute `(corank, sorted_nontrivial_torsion_factors)` of an integer relation matrix."""
    if num_cols == 0:
        return (0, ())
    sp_rows: list[dict[int, int]] = []
    for r in rows:
        d = {j: int(v) for j, v in enumerate(r) if v != 0}
        if d:
            sp_rows.append(d)
    if not sp_rows:
        return (num_cols, ())

    rem_cols = num_cols
    # Sparse Markowitz `+-1` unit-pivot Gaussian elimination over Z:
    while sp_rows and rem_cols > 0:
        best_ri = -1
        best_col = -1
        best_len = 1 << 30
        for ri, rd in enumerate(sp_rows):
            rlen = len(rd)
            if rlen >= best_len:
                continue
            for cj, val in rd.items():
                if val in (1, -1) and max(abs(x) for x in rd.values()) <= 1000:
                    best_ri = ri
                    best_col = cj
                    best_len = rlen
                    break
            if best_len <= 2:
                break
        if best_ri == -1:
            break
        piv_row = sp_rows[best_ri]
        piv_sign = piv_row[best_col]  # +1 or -1: X_{best_col} = -piv_sign * sum_{k != best_col} piv_row[k] X_k
        repl = {k: -piv_sign * v for k, v in piv_row.items() if k != best_col}
        nxt_rows: list[dict[int, int]] = []
        for ri, rd in enumerate(sp_rows):
            if ri == best_ri:
                continue
            c_val = rd.get(best_col, 0)
            if c_val == 0:
                nxt_rows.append(rd)
            else:
                nd = {k: v for k, v in rd.items() if k != best_col}
                for rk, rv in repl.items():
                    nv = nd.get(rk, 0) + c_val * rv
                    if nv:
                        nd[rk] = nv
                    else:
                        nd.pop(rk, None)
                if nd:
                    nxt_rows.append(nd)
        sp_rows = nxt_rows
        rem_cols -= 1

    if not sp_rows:
        return (rem_cols, ())

    active_cols = sorted({cj for rd in sp_rows for cj in rd})
    unused_cols = rem_cols - len(active_cols)
    col_map = {cj: idx for idx, cj in enumerate(active_cols)}
    dense_rows = []
    for rd in sp_rows:
        dr = [0] * len(active_cols)
        for cj, val in rd.items():
            dr[col_map[cj]] = val
        dense_rows.append(dr)

    mat = sp.Matrix(dense_rows)
    snf = smith_normal_form(mat, domain=sp.ZZ)
    diag = [abs(int(snf[i, i])) for i in range(min(snf.rows, snf.cols))]
    nonzero = [x for x in diag if x > 0]
    free_rank = unused_cols + (len(active_cols) - len(nonzero))
    torsion = tuple(sorted(x for x in nonzero if x > 1))
    return (free_rank, torsion)


def kei_module_invariants(
    d: MarkedDiagram,
    pres: GroupPresentation | None = None,
) -> tuple[int, tuple[int, ...]]:
    """Return `(free_rank, torsion)` of the universal Fox/kei coloring module (`x |> y = 2y - x`).

    Note: for any dihedral kei `R_n` (`n >= 2`), the exact number of `R_n`-colorings is:
    `|Col_{R_n}(F)| = n^{free_rank} * prod_{t in torsion} gcd(n, t)`.
    """
    if pres is None:
        pres = wirtinger_presentation(d, from_resolution=-1, simplify=True)
    g = pres.num_generators
    if g == 0:
        return (0, ())
    # Evaluate Fox Jacobian at `psi(x_i) = -1`!
    #Wait: in `Hom(pi_1, Z_n rtimes_{-1} Z_2)`, mapping every meridian generator `x_i` to a
    # reflection `(c_i, -1)` satisfies `x_i^2 = 1` automatically! Let's compute directly on
    # over-arcs of `L_-(D)` (where every meridian generator IS a reflection `(c_i, -1)`!):
    # Don't use Tietze-simplified `pres` unless all generators are meridians (Tietze preserves
    # meridian reflections ONLY if every generator maps to `-1` in `{+1, -1}`!).
    # Instead, build directly from un-simplified `wirtinger_presentation(d, simplify=False)` where
    # EVERY generator IS a genuine meridian over-arc!
    raw_pres = wirtinger_presentation(d, from_resolution=-1, simplify=False)
    ng = raw_pres.num_generators
    if ng == 0:
        return (0, ())
    rows: list[list[int]] = []
    for rel in raw_pres.relators:
        # Map generator `j` -> `(c_j, -1)` in `Z rtimes_{-1} Z_2`:
        # `(a, s1) * (b, s2) = (a + s1 * b, s1 * s2)`
        # `(c_j, -1)^{-1} = (c_j, -1)`!
        # Wait: every relator in `raw_pres` is either:
        # - Wirtinger `(+-z, x, -+z, -y)` (length 4 -> sign = +1, linear part `2 c_z - c_x - c_y = 0`), or
        # - Saddle `(+-a1, -+a2)` (length 2 -> sign = +1, linear part `c_{a1} - c_{a2} = 0`!!)!
        # Check `(c_1, -1) * (c_2, -1) = (c_1 - c_2, +1)` -> independent of `lit` signs!
        vec = [0] * ng
        sign = 1
        for lit in rel:
            j = abs(lit) - 1
            vec[j] += sign
            sign = -sign
        # Note: if `sign == -1` (impossible on admissible diagrams since every relator has even length!),
        # no reflection coloring exists; on admissible diagrams `sign == +1` always!
        assert sign == 1
        if any(vec):
            rows.append(vec)
    return _snf_rank_and_torsion(rows, ng)


def dihedral_coloring_count(free_rank: int, torsion: tuple[int, ...], n: int) -> int:
    """Exact number of `R_n`-colorings `|Col_{R_n}(F)|` from `kei_module_invariants`."""
    ans = n**free_rank
    for t in torsion:
        ans *= gcd(n, t)
    return ans


def _diagram_edge_ids(d: MarkedDiagram) -> tuple[list[int], int]:
    """Assign undirected edge IDs `0..num_edges-1` to darts `0..4n-1` (`num_edges = 2n`)."""
    n = d.num_vertices
    dart_edge = [-1] * (4 * n)
    eid = 0
    for hd in range(4 * n):
        if dart_edge[hd] >= 0:
            continue
        od = d.opp[hd]
        dart_edge[hd] = eid
        dart_edge[od] = eid
        eid += 1
    return dart_edge, eid


def bikei_module_invariants(d: MarkedDiagram) -> tuple[int, tuple[int, ...]]:
    """Universal linear involutory biquandle (bikei) invariant over `R = Z[s] / (s^2 - 1)`.

    Works on ALL surface-links (orientable, non-orientable, or mixed), requiring NO orientation!
    Each undirected graph edge `e` (plus each `free_loop`) carries `a_e + b_e * s` (`2` Z-coords):
    - At marked vertex `v`: all 4 incident edges agree (`x_{e0} = x_{e1} = x_{e2} = x_{e3}`).
    - At classical vertex `u` with under-edges `(u1, u2)` and over-edges `(o1, o2)`:
      `o2 - s * o1 = 0` and `u1 + u2 - (1 + s) * o1 = 0`.
    Returns `(free_Z_rank, torsion_factors)` of this `Z`-module.
    """
    n = d.num_vertices
    if n == 0:
        return (2 * d.free_loops, ())
    dart_edge, num_edges = _diagram_edge_ids(d)
    # Unimodular s-monomial elimination (`s^2 = 1`): `X_e = s^{parity[e]} * Y_{root[e]}`
    total_edges = num_edges + d.free_loops
    parent = list(range(total_edges))
    s_par = [0] * total_edges  # 0 => 1, 1 => s

    def find_edge(x: int) -> tuple[int, int]:
        p = 0
        while parent[x] != x:
            p ^= s_par[x]
            x = parent[x]
        return x, p

    deferred_eqs: list[list[tuple[int, int, int]]] = []

    def equate_edges(ea: int, eb: int, rel_s: int) -> None:
        """Enforce `X_{ea} = s^{rel_s} * X_{eb}` (`rel_s in {0, 1}`)."""
        ra, pa = find_edge(ea)
        rb, pb = find_edge(eb)
        need = pa ^ pb ^ rel_s
        if ra != rb:
            parent[ra] = rb
            s_par[ra] = need
        elif need == 1:
            deferred_eqs.append([(rb, 1, -1)])

    for v in range(n):
        vt = d.vtype[v]
        e0 = dart_edge[4 * v + 0]
        e1 = dart_edge[4 * v + 1]
        e2 = dart_edge[4 * v + 2]
        e3 = dart_edge[4 * v + 3]
        if vt in ("Y", "Z"):
            for ek in (e1, e2, e3):
                if ek != e0:
                    equate_edges(ek, e0, 0)
        else:
            if vt == "X0":
                u1, u2, o1, o2 = e0, e2, e1, e3
            else:
                u1, u2, o1, o2 = e1, e3, e0, e2
            equate_edges(o2, o1, 1)
            deferred_eqs.append([(u1, 1, 0), (u2, 1, 0), (o1, -1, -1)])

    roots = [i for i in range(total_edges) if parent[i] == i]
    root_idx = {r: idx for idx, r in enumerate(roots)}
    total_vars = 2 * len(roots)
    rows: list[list[int]] = []

    for terms in deferred_eqs:
        r1 = [0] * total_vars
        rs = [0] * total_vars
        for e_orig, c1, cs in terms:
            rt, pe = find_edge(e_orig)
            ri = root_idx[rt]
            # Multiply `(c1 + cs*s)` by `s^{pe}`: if `pe == 1`, swap `(c1, cs)`!
            nc1, ncs = (cs, c1) if pe else (c1, cs)
            r1[2 * ri] += nc1
            r1[2 * ri + 1] += ncs
            rs[2 * ri] += ncs
            rs[2 * ri + 1] += nc1
        if any(r1):
            rows.append(r1)
        if any(rs):
            rows.append(rs)

    return _snf_rank_and_torsion(rows, total_vars)


def _orientable_edge_orientations(d: MarkedDiagram) -> list[list[int]] | None:
    """If `F(D)` is orientable, return up to `2^{c-1}` coherent dart orientations `dart_in_out`
    (`+1` = incoming to vertex, `-1` = outgoing from vertex) modulo global reversal!
    Returns `None` if `F(D)` is non-orientable.
    """
    n = d.num_vertices
    if n == 0:
        return [[]]
    dart_scomp, total_scomps = d._dart_surface_components()
    graph_scomps = total_scomps - d.free_loops

    # Check bipartite orientation consistency per surface component:
    # Constraints:
    # - Across graph edge: `dir[opp[hd]] == -dir[hd]`
    # - At classical `v`: `dir[4v + (s^2)] == -dir[4v + s]`
    # - At marked `v`: `dir[4v + ((s+1)&3)] == -dir[4v + s]` (so `0, 2` same, `1, 3` opposite!)
    base_dirs_per_comp: list[dict[int, int]] = []
    comp_ids_seen = sorted({dart_scomp[hd] for hd in range(4 * n)})

    for sc in comp_ids_seen:
        start_hd = next(hd for hd in range(4 * n) if dart_scomp[hd] == sc)
        cmap: dict[int, int] = {start_hd: 1}
        st = [start_hd]
        while st:
            hd = st.pop()
            val = cmap[hd]
            v, s = hd >> 2, hd & 3
            vt = d.vtype[v]
            # Edge Constraint
            od = d.opp[hd]
            if od in cmap:
                if cmap[od] != -val:
                    return None
            else:
                cmap[od] = -val
                st.append(od)
            # Vertex Constraint
            if vt in ("X0", "X1"):
                vd = (v << 2) | (s ^ 2)
                if vd in cmap:
                    if cmap[vd] != -val:
                        return None
                else:
                    cmap[vd] = -val
                    st.append(vd)
            else:
                for adj_s in ((s + 1) & 3, (s + 3) & 3):
                    vd = (v << 2) | adj_s
                    if vd in cmap:
                        if cmap[vd] != -val:
                            return None
                    else:
                        cmap[vd] = -val
                        st.append(vd)
        base_dirs_per_comp.append(cmap)

    if not base_dirs_per_comp:
        return [[]]

    # Fix first component sign to `+1`, vary signs of remaining components in `{+1, -1}`
    k = len(base_dirs_per_comp)
    out: list[list[int]] = []
    for signs in product((1, -1), repeat=k - 1):
        full_signs = (1,) + signs
        arr = [0] * (4 * n)
        for cmap, mult in zip(base_dirs_per_comp, full_signs):
            for hd, dv in cmap.items():
                arr[hd] = mult * dv
        out.append(arr)
    return out


def alexander_biquandle_invariants(
    d: MarkedDiagram,
) -> tuple[tuple, tuple[str, ...]] | None:
    """Compute oriented Alexander biquandle invariants over `Z[t^{+-1}, s^{+-1}]` if `F(D)` is orientable.

    Returns `(snf_profile, s_monodromy_polynomials)` or `None` if `F(D)` is non-orientable.
    - `snf_profile` is the canonical sorted tuple (across relative component orientations) of exact
      Smith Normal Form `(free_rank, torsion)` over the finite group rings `Z[t, s] / (t^a - 1, s^b - 1)`
      for `(a, b) in ((1, 2), (1, 3), (2, 2), (2, 3), (3, 2))`.
      Note: because `t -> t^{-1}, s -> s^{-1}` is an automorphism of `Z[t,s]/(t^a-1, s^b-1)`, every
      SNF tuple is **100% invariant under global orientation reversal AND ambient isotopy**!
    - `s_monodromy_polynomials`: symbolic generators of `E_1^{BQ}(t=1, s)` in `Z[s]`, capturing pure
      over-sheet `s`-twisting even when `pi_1` is abelian (`8_1^{1,1}` vs `10_1^{1,1}`)!
    """
    orientations = _orientable_edge_orientations(d)
    if orientations is None:
        return None
    n = d.num_vertices
    if n == 0:
        c = d.free_loops
        specs = ((1, 2), (1, 3), (2, 2), (2, 3), (3, 2))
        prof = tuple(((a * b * c, ()) for a, b in specs))
        return ((prof,) * (1 << max(0, c - 1)), ("1" if c <= 1 else "0",))

    dart_edge, num_edges = _diagram_edge_ids(d)
    total_edges = num_edges + d.free_loops
    specs = ((1, 2), (1, 3), (2, 2), (2, 3), (3, 2))

    all_orient_profiles = []
    s_sym = sp.Symbol("s")
    s_monodromy_strs: set[str] = set()

    for dart_dir in orientations:
        # Collect symbolic equations in `(e_idx, t_exp, s_exp, coeff)` where sum == 0:
        # Positive crossing (`diff == 3`):
        #   o_out - s * o_in = 0
        #   u_out - t * u_in - (1 - t*s) * o_in = 0
        # Negative crossing (`diff == 1`):
        #   o_out - s^{-1} * o_in = 0   <=>  s * o_out - o_in = 0
        #   u_out - t^{-1} * u_in - (1 - t^{-1}*s^{-1}) * o_in = 0
        #     <=>  t*s * u_out - s * u_in - (t*s - 1) * o_in = 0
        eqs: list[list[tuple[int, int, int, int]]] = []
        for v in range(n):
            vt = d.vtype[v]
            if vt in ("Y", "Z"):
                e0 = dart_edge[4 * v + 0]
                for s_slot in (1, 2, 3):
                    ek = dart_edge[4 * v + s_slot]
                    if ek != e0:
                        eqs.append([(e0, 0, 0, 1), (ek, 0, 0, -1)])
            else:
                under_slots = (0, 2) if vt == "X0" else (1, 3)
                over_slots = (1, 3) if vt == "X0" else (0, 2)
                u_in_s = under_slots[0] if dart_dir[4 * v + under_slots[0]] == 1 else under_slots[1]
                u_out_s = u_in_s ^ 2
                o_in_s = over_slots[0] if dart_dir[4 * v + over_slots[0]] == 1 else over_slots[1]
                o_out_s = o_in_s ^ 2

                eu_in = dart_edge[4 * v + u_in_s]
                eu_out = dart_edge[4 * v + u_out_s]
                eo_in = dart_edge[4 * v + o_in_s]
                eo_out = dart_edge[4 * v + o_out_s]

                diff = (o_in_s - u_in_s) & 3
                if diff == 3:
                    # Positive crossing
                    eqs.append([(eo_out, 0, 0, 1), (eo_in, 0, 1, -1)])
                    eqs.append(
                        [
                            (eu_out, 0, 0, 1),
                            (eu_in, 1, 0, -1),
                            (eo_in, 0, 0, -1),
                            (eo_in, 1, 1, 1),
                        ]
                    )
                else:
                    # Negative crossing
                    eqs.append([(eo_out, 0, 1, 1), (eo_in, 0, 0, -1)])
                    eqs.append(
                        [
                            (eu_out, 1, 1, 1),
                            (eu_in, 0, 1, -1),
                            (eo_in, 1, 1, -1),
                            (eo_in, 0, 0, 1),
                        ]
                    )

        # Unimodular Laurent-monomial elimination over `Z[t^{+-1}, s^{+-1}]`:
        # Any 2-term equation `t^{te1} s^{se1} X_{e1} - t^{te2} s^{se2} X_{e2} = 0` with `root(e1) != root(e2)`
        # substitutes `X_{r1} = t^{dt} s^{ds} X_{r2}` isomorphically across all quotient rings!
        uf_p = list(range(total_edges))
        uf_dt = [0] * total_edges
        uf_ds = [0] * total_edges

        def uf_find(x: int) -> tuple[int, int, int]:
            dt_acc = 0
            ds_acc = 0
            while uf_p[x] != x:
                dt_acc += uf_dt[x]
                ds_acc += uf_ds[x]
                x = uf_p[x]
            return x, dt_acc, ds_acc

        rem_eqs_raw: list[list[tuple[int, int, int, int]]] = []
        for eq in eqs:
            if len(eq) == 2 and sorted((eq[0][3], eq[1][3])) == [-1, 1]:
                (e_pos, te_pos, se_pos, c_pos), (e_neg, te_neg, se_neg, _) = (
                    (eq[0], eq[1]) if eq[0][3] == 1 else (eq[1], eq[0])
                )
                assert c_pos == 1
                r1, dt1, ds1 = uf_find(e_pos)
                r2, dt2, ds2 = uf_find(e_neg)
                if r1 != r2:
                    # `t^{te_pos + dt1} s^{se_pos + ds1} X_{r1} = t^{te_neg + dt2} s^{se_neg + ds2} X_{r2}`
                    uf_p[r1] = r2
                    uf_dt[r1] = (te_neg + dt2) - (te_pos + dt1)
                    uf_ds[r1] = (se_neg + ds2) - (se_pos + ds1)
                    continue
            rem_eqs_raw.append(eq)

        roots_abq = [i for i in range(total_edges) if uf_p[i] == i]
        root_map_abq = {r: idx for idx, r in enumerate(roots_abq)}
        num_red_edges = len(roots_abq)
        red_eqs: list[list[tuple[int, int, int, int]]] = []
        for eq in rem_eqs_raw:
            acc: dict[tuple[int, int, int], int] = {}
            for e_idx, te, se, coeff in eq:
                rt, dte, dse = uf_find(e_idx)
                key = (root_map_abq[rt], te + dte, se + dse)
                acc[key] = acc.get(key, 0) + coeff
            cleaned = [(ri, t_e, s_e, cf) for (ri, t_e, s_e), cf in acc.items() if cf != 0]
            if cleaned:
                red_eqs.append(cleaned)

        # 1. Compute SNF over `Z[t, s] / (t^a - 1, s^b - 1)` for each `(a, b)` in `specs`!
        cur_prof = []
        for a_mod, b_mod in specs:
            dim = a_mod * b_mod
            num_cols = num_red_edges * dim
            z_rows: list[list[int]] = []
            for eq in red_eqs:
                for bt in range(a_mod):
                    for bs in range(b_mod):
                        row = [0] * num_cols
                        for e_idx, te, se, coeff in eq:
                            nt = (bt + te) % a_mod
                            ns = (bs + se) % b_mod
                            col = e_idx * dim + (nt * b_mod + ns)
                            row[col] += coeff
                        if any(row):
                            z_rows.append(row)
            cur_prof.append(_snf_rank_and_torsion(z_rows, num_cols))
        all_orient_profiles.extend([tuple(cur_prof)] * (1 << d.free_loops))

        # 2. Symbolic `t = 1` specialization in `Z[s]` (pure `s`-monodromy polynomial!):
        sym_rows: list[dict[int, sp.Expr]] = []
        for eq in red_eqs:
            rd: dict[int, sp.Expr] = {}
            for e_idx, _te, se, coeff in eq:
                term = coeff * (s_sym**se)
                rd[e_idx] = sp.expand(rd.get(e_idx, sp.Integer(0)) + term)
            rd = {k: v for k, v in rd.items() if v != 0}
            if rd:
                sym_rows.append(rd)

        active_vars = list(range(num_red_edges))

        def _laurent_unit_inv(expr: sp.Expr) -> sp.Expr | None:
            for sign, e_abs in ((1, expr), (-1, -expr)):
                if e_abs == 1:
                    return sp.Integer(sign)
                if e_abs == s_sym:
                    return sign * (s_sym**-1)
                if isinstance(e_abs, sp.Pow) and e_abs.base == s_sym and e_abs.exp.is_Integer:
                    return sign * (s_sym ** (-e_abs.exp))
            return None

        # Eliminate equations where a variable has a Laurent-monomial unit coefficient `+- s^k`!
        changed = True
        while changed and len(active_vars) > 1:
            changed = False
            for r_i, rd in enumerate(sym_rows):
                piv_var = None
                piv_inv = sp.Integer(0)
                for k, v_expr in rd.items():
                    inv_ex = _laurent_unit_inv(v_expr)
                    if inv_ex is not None:
                        piv_var = k
                        piv_inv = inv_ex
                        break
                if piv_var is not None:
                    # `v_expr * X_{piv} + rest = 0` => `X_{piv} = -piv_inv * rest`
                    repl = {
                        k: sp.expand(-piv_inv * v_expr)
                        for k, v_expr in rd.items()
                        if k != piv_var
                    }
                    new_sym_rows: list[dict[int, sp.Expr]] = []
                    for r_j, rd2 in enumerate(sym_rows):
                        if r_j == r_i:
                            continue
                        if piv_var not in rd2:
                            new_sym_rows.append(rd2)
                        else:
                            c_piv = rd2[piv_var]
                            nd = {k: v for k, v in rd2.items() if k != piv_var}
                            for rk, rv in repl.items():
                                nd[rk] = sp.expand(nd.get(rk, sp.Integer(0)) + c_piv * rv)
                            nd = {k: v for k, v in nd.items() if v != 0}
                            if nd:
                                new_sym_rows.append(nd)
                    sym_rows = new_sym_rows
                    active_vars.remove(piv_var)
                    changed = True
                    break

        def _to_s_poly(expr: sp.Expr) -> sp.Poly:
            num, _den = sp.fraction(sp.together(expr))
            return sp.Poly(sp.expand(num), s_sym, domain=sp.ZZ)

        rem_g = len(active_vars)
        if rem_g <= 1:
            if not sym_rows:
                s_monodromy_strs.add("1")
            else:
                v0 = active_vars[0]
                polys = [rd.get(v0, sp.Integer(0)) for rd in sym_rows if rd.get(v0, sp.Integer(0)) != 0]
                cur_g = _to_s_poly(polys[0])
                for p_ex in polys[1:]:
                    cur_g = sp.gcd(cur_g, _to_s_poly(p_ex))
                s_monodromy_strs.add(_canonicalize_s_poly(cur_g.as_expr(), s_sym))
        else:
            minor_sz = rem_g - 1
            if len(sym_rows) < minor_sz:
                s_monodromy_strs.add("0")
            elif rem_g <= 8 and len(sym_rows) <= 14:
                from itertools import combinations

                v_map = {v: idx for idx, v in enumerate(active_vars)}
                mat = sp.zeros(len(sym_rows), rem_g)
                for ri, rd in enumerate(sym_rows):
                    for vk, vv in rd.items():
                        num_vv, den_vv = sp.fraction(sp.together(vv))
                        mat[ri, v_map[vk]] = sp.expand(num_vv) if den_vv == 1 else sp.together(vv)
                # Multiply each row by highest needed power of `s_sym` so `mat` is strictly polynomial!
                for ri in range(len(sym_rows)):
                    min_p = 0
                    for ci in range(rem_g):
                        entry = mat[ri, ci]
                        if entry == 0:
                            continue
                        num_e, den_e = sp.fraction(sp.together(entry))
                        if den_e != 1:
                            d_poly = sp.Poly(sp.expand(den_e), s_sym, domain=sp.ZZ)
                            min_p = max(min_p, d_poly.degree())
                    if min_p > 0:
                        for ci in range(rem_g):
                            if mat[ri, ci] != 0:
                                mat[ri, ci] = sp.expand(sp.simplify(mat[ri, ci] * (s_sym**min_p)))
                cur_g = None
                for rc in combinations(range(len(sym_rows)), minor_sz):
                    for cc in combinations(range(rem_g), minor_sz):
                        det = sp.expand(mat.extract(list(rc), list(cc)).det(method="berkowitz"))
                        if det != 0:
                            dp = sp.Poly(det, s_sym, domain=sp.ZZ)
                            # Strip `s^k` factor before gcd so Laurent units don't mask `1`!
                            terms = [(m[0], int(c)) for m, c in dp.terms() if c != 0]
                            m0 = min(e for e, _ in terms)
                            if m0 > 0:
                                dp = sp.Poly(
                                    sum(sp.Integer(c) * (s_sym ** (e - m0)) for e, c in terms),
                                    s_sym,
                                    domain=sp.ZZ,
                                )
                            cur_g = dp if cur_g is None else sp.gcd(cur_g, dp)
                            if cur_g.as_expr() in (1, -1):
                                break
                    if cur_g is not None and cur_g.as_expr() in (1, -1):
                        break
                s_monodromy_strs.add(
                    "0" if cur_g is None else _canonicalize_s_poly(cur_g.as_expr(), s_sym)
                )

    all_orient_profiles.sort()
    return (tuple(all_orient_profiles), tuple(sorted(s_monodromy_strs)))


def _canonicalize_s_poly(expr: sp.Expr, s_sym: sp.Symbol) -> str:
    """Strip monomial `s^k` unit factors and pick canonical representative between `f(s)` and `f(s^{-1})`."""
    if expr == 0:
        return "0"
    p = sp.Poly(sp.expand(expr), s_sym, domain=sp.ZZ)
    terms = [(mon[0], int(c)) for mon, c in p.terms() if c != 0]
    if not terms:
        return "0"
    min_e = min(e for e, _ in terms)
    max_e = max(e for e, _ in terms)
    deg = max_e - min_e
    fwd = sorted(((e - min_e, c) for e, c in terms), reverse=True)
    if fwd[0][1] < 0:
        fwd = [(e, -c) for e, c in fwd]
    rev = sorted(((deg - (e - min_e), c) for e, c in terms), reverse=True)
    if rev[0][1] < 0:
        rev = [(e, -c) for e, c in rev]
    best = min(fwd, rev)
    out = sum(sp.Integer(c) * (s_sym**e) for e, c in best)
    return str(sp.expand(out))


def compute_biquandle_summary(
    d: MarkedDiagram,
    pres: GroupPresentation | None = None,
) -> BiquandleResult:
    """Compute combined kei, universal linear bikei, Alexander biquandle, and 3-cocycle invariants."""
    k_rk, k_tor = kei_module_invariants(d, pres=pres)
    bk_rk, bk_tor = bikei_module_invariants(d)
    abq = alexander_biquandle_invariants(d)
    c_r3 = dihedral_coloring_count(k_rk, k_tor, 3)
    c_r5 = dihedral_coloring_count(k_rk, k_tor, 5)
    coc_3 = mochizuki_3cocycle_invariant(d, p=3, known_col_count=c_r3)
    coc_5 = mochizuki_3cocycle_invariant(d, p=5, known_col_count=c_r5)
    return BiquandleResult(
        kei_free_rank=k_rk,
        kei_torsion=k_tor,
        col_r3=c_r3,
        col_r4=dihedral_coloring_count(k_rk, k_tor, 4),
        col_r5=c_r5,
        bikei_free_rank=bk_rk,
        bikei_torsion=bk_tor,
        oriented_abq_snf_profile=abq[0] if abq is not None else None,
        oriented_s_monodromy_poly=abq[1] if abq is not None else None,
        cocycle_r3=coc_3,
        cocycle_r5=coc_5,
    )


def _mochizuki_theta(x: int, y: int, z: int, p: int) -> int:
    """Evaluate Mochizuki's 3-cocycle `theta_p(x, y, z) in Z_p` on dihedral quandle `R_p` (`p` odd prime):
    `theta_p(x, y, z) = (x - y) * (((2z - y) mod p)^p + y^p - 2 * z^p) / p (mod p)`.
    """
    if x == y or y == z:
        return 0
    w = (2 * z - y) % p
    num = (w**p + y**p - 2 * (z**p)) // p
    return ((x - y) * num) % p


def _format_zp_state_sum(counts: Counter[int], p: int) -> str:
    """Format state-sum `sum_{k=0}^{p-1} c_k * u^k` in `Z[Z_p]` as `'3 + 6*u'`."""
    terms = []
    if counts.get(0, 0):
        terms.append(str(counts[0]))
    for k in range(1, p):
        c = counts.get(k, 0)
        if c:
            u_part = "u" if k == 1 else f"u^{k}"
            terms.append(u_part if c == 1 else f"{c}*{u_part}")
    return " + ".join(terms) if terms else "0"


def _evaluate_classical_movie_cocycle(
    movie: list,
    init_col: list[int],
    init_dir: list[int],
    p: int,
) -> int:
    """Transport `(col, dart_dir)` along `movie` (`L -> O_k`) and return signed sum of `theta_p` mod `p`."""
    from surfaceknots.unlink import _is_over as cl_is_over

    col = list(init_col)
    ddir = list(init_dir)
    total_w = 0

    for step in movie:
        kind = step.kind
        aft = step.after
        meta = step.meta

        if kind in ("R1-", "R2-"):
            old_to_new: dict[int, int] = meta["old_to_new"]
            m = aft.num_crossings
            new_col = [0] * (4 * m)
            new_dir = [0] * (4 * m)
            for vo, vn in old_to_new.items():
                for s in range(4):
                    new_col[4 * vn + s] = col[4 * vo + s]
                    new_dir[4 * vn + s] = ddir[4 * vo + s]
            col, ddir = new_col, new_dir

        elif kind == "R2+":
            da, oda, db, odb, vt = meta["da"], meta["oda"], meta["db"], meta["odb"], meta["vt"]
            n_prev = step.before.num_crossings
            u, w = n_prev, n_prev + 1
            c_a = col[da]
            c_b = col[db]
            dir_da = ddir[da]
            dir_db = ddir[db]
            new_col = list(col) + [0] * 8
            new_dir = list(ddir) + [0] * 8
            # u: 0=da, 1=w.3, 2=w.2, 3=odb; w: 0=oda, 1=db, 2=u.2, 3=u.1
            new_dir[4 * u + 0] = -dir_da
            new_dir[4 * u + 2] = dir_da
            new_dir[4 * w + 2] = -dir_da
            new_dir[4 * w + 0] = dir_da

            new_dir[4 * w + 1] = -dir_db
            new_dir[4 * w + 3] = dir_db
            new_dir[4 * u + 1] = -dir_db
            new_dir[4 * u + 3] = dir_db

            if vt == "X1":  # A=(0,2) over B=(1,3)
                c_b_mid = (2 * c_a - c_b) % p
                new_col[4 * u + 0] = new_col[4 * u + 2] = c_a
                new_col[4 * w + 0] = new_col[4 * w + 2] = c_a
                new_col[4 * u + 3] = new_col[4 * w + 1] = c_b
                new_col[4 * u + 1] = new_col[4 * w + 3] = c_b_mid
            else:  # X0: B=(1,3) over A=(0,2)
                c_a_mid = (2 * c_b - c_a) % p
                new_col[4 * u + 1] = new_col[4 * u + 3] = c_b
                new_col[4 * w + 1] = new_col[4 * w + 3] = c_b
                new_col[4 * u + 0] = new_col[4 * w + 0] = c_a
                new_col[4 * u + 2] = new_col[4 * w + 2] = c_a_mid
            col, ddir = new_col, new_dir

        elif kind == "R3":
            vs = meta["vs"]
            ss = meta["ss"]
            heights = meta["heights"]
            ports = meta["ports"]
            port_pairs = meta["port_pairs"]

            d_tri = [(vs[i] << 2) | ss[i] for i in range(3)]
            sigma = [-ddir[d_tri[i]] for i in range(3)]
            c_tri = [col[d_tri[i]] for i in range(3)]

            def ext_port_across(strand_i: int, other_j: int) -> int:
                if other_j == (strand_i - 1) % 3:
                    return ports[2 * strand_i]
                return ports[(2 * strand_i + 3) % 6]

            idx_T = heights.index(2)
            idx_M = heights.index(1)
            idx_B = heights.index(0)

            z = c_tri[idx_T]
            sig_T, sig_M, sig_B = sigma[idx_T], sigma[idx_M], sigma[idx_B]
            y = c_tri[idx_M] if sig_T == -1 else col[ext_port_across(idx_M, idx_T)]

            if (sig_T, sig_M) == (-1, -1):
                x = c_tri[idx_B]
            elif (sig_T, sig_M) == (1, -1):
                x = col[ext_port_across(idx_B, idx_T)]
            elif (sig_T, sig_M) == (-1, 1):
                x = col[ext_port_across(idx_B, idx_M)]
            else:
                x = (2 * y - col[ext_port_across(idx_B, idx_T)]) % p

            perm_sign = 1 if (idx_T, idx_M, idx_B) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)) else -1
            eps = perm_sign * sig_T * sig_M * sig_B
            total_w = (total_w + eps * _mochizuki_theta(x, y, z, p)) % p

            # Update `(col, ddir)` on `after`
            new_col = list(col)
            new_dir = list(ddir)
            for i in range(3):
                wi = vs[i]
                pi_a, pi_b = port_pairs[i]
                new_col[4 * wi + 2] = col[ports[pi_a]]
                new_dir[4 * wi + 2] = ddir[ports[pi_a]]
                new_col[4 * wi + 3] = col[ports[pi_b]]
                new_dir[4 * wi + 3] = ddir[ports[pi_b]]
                new_dir[4 * wi + 0] = -new_dir[4 * wi + 2]
                new_dir[4 * wi + 1] = -new_dir[4 * wi + 3]

                vt_wi = aft.vtype[wi]
                if cl_is_over(vt_wi, 2):
                    new_col[4 * wi + 0] = new_col[4 * wi + 2]
                else:
                    new_col[4 * wi + 0] = (2 * new_col[4 * wi + 3] - new_col[4 * wi + 2]) % p

                if cl_is_over(vt_wi, 3):
                    new_col[4 * wi + 1] = new_col[4 * wi + 3]
                else:
                    new_col[4 * wi + 1] = (2 * new_col[4 * wi + 2] - new_col[4 * wi + 3]) % p

            col, ddir = new_col, new_dir

        else:
            raise ValueError(f"Unsupported movie step kind {kind!r}")

        # Self-check edge color and orientation consistency at every movie step!
        for hd in range(4 * aft.num_crossings):
            od = aft.opp[hd]
            assert col[od] == col[hd], f"Color mismatch after {kind}!"
            assert ddir[od] == -ddir[hd], f"Orientation mismatch after {kind}!"

    return total_w


def mochizuki_3cocycle_invariant(
    d: MarkedDiagram,
    p: int = 3,
    known_col_count: int | None = None,
    up_to_mirror: bool = True,
) -> tuple[str, ...] | None:
    """Compute Mochizuki's dihedral quandle 3-cocycle state-sum invariant `Phi_{theta_p}(F)` in `Z[Z_p]`
    across all `2^c` surface orientations of `F(D)` (returns `None` if `F(D)` is non-orientable).

    - When `up_to_mirror=False`, returns the exact sorted `2^c`-tuple of state-sums in fixed oriented
      `R^4` (evaluating both `+o` and `-o` directly on the unlinking movies). Comparing `d` against
      `d.mirror_3d()` certifies **4D chirality** (`F != F*`), while distinct entries within `d`
      certify **orientation sensitivity / non-invertibility** (`+F != -F`).
    - When `up_to_mirror=True` (default for Yoshikawa unoriented/unmirrored census comparison),
      canonicalizes between `Phi` and its mirror-conjugate `bar(Phi)` (`u^k -> u^{-k mod p}`).
    - Any state-sum with nontrivial `u^k` (`k != 0`) terms certifies `t(F) > 0` (non-ribbon!).
    """
    from surfaceknots.unlink import find_unlink_movie

    base_oris = _orientable_edge_orientations(d)
    if base_oris is None:
        return None

    if known_col_count is None:
        rk, tor = kei_module_invariants(d)
        known_col_count = dihedral_coloring_count(rk, tor, p)

    num_scomps = d.surface_topology().num_components
    if d.num_vertices == 0 or known_col_count == p:
        s_str = str(known_col_count)
        return tuple([s_str] * (1 << num_scomps))

    if known_col_count > 729:
        return None

    n = d.num_vertices
    classical_vs = [v for v in range(n) if d.vtype[v] in ("X0", "X1")]
    l_plus = d.resolve(1)
    l_minus = d.resolve(-1)

    mov_plus = find_unlink_movie(l_plus)
    mov_minus = find_unlink_movie(l_minus)
    if mov_plus is None or mov_minus is None:
        return None

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

    for hd in range(4 * n):
        union(hd, d.opp[hd])
    for v in range(n):
        vt = d.vtype[v]
        if vt in ("Y", "Z"):
            union(4 * v + 0, 4 * v + 1)
            union(4 * v + 1, 4 * v + 2)
            union(4 * v + 2, 4 * v + 3)
        elif vt == "X0":
            union(4 * v + 1, 4 * v + 3)
        else:
            union(4 * v + 0, 4 * v + 2)

    roots = sorted({find(hd) for hd in range(4 * n)})
    r_to_id = {r: idx for idx, r in enumerate(roots)}
    dart_var = [r_to_id[find(hd)] for hd in range(4 * n)]
    num_vars = len(roots) + d.free_loops

    eq_rows: list[list[int]] = []
    for v in classical_vs:
        vt = d.vtype[v]
        if vt == "X0":
            u1, u2, ov = dart_var[4 * v + 0], dart_var[4 * v + 2], dart_var[4 * v + 1]
        else:
            u1, u2, ov = dart_var[4 * v + 1], dart_var[4 * v + 3], dart_var[4 * v + 0]
        row = [0] * num_vars
        row[u1] = (row[u1] + 1) % p
        row[u2] = (row[u2] + 1) % p
        row[ov] = (row[ov] - 2) % p
        if any(row):
            eq_rows.append(row)

    piv_cols: list[int] = []
    rref = [r[:] for r in eq_rows]
    pr = 0
    for col_i in range(num_vars):
        found = -1
        for rr in range(pr, len(rref)):
            if rref[rr][col_i] % p != 0:
                found = rr
                break
        if found == -1:
            continue
        rref[pr], rref[found] = rref[found], rref[pr]
        inv_p = pow(rref[pr][col_i], p - 2, p)
        rref[pr] = [(val * inv_p) % p for val in rref[pr]]
        for rr in range(len(rref)):
            if rr != pr and rref[rr][col_i] != 0:
                fac = rref[rr][col_i]
                rref[rr] = [(rref[rr][cc] - fac * rref[pr][cc]) % p for cc in range(num_vars)]
        piv_cols.append(col_i)
        pr += 1

    piv_set = set(piv_cols)
    free_cols = [c_i for c_i in range(num_vars) if c_i not in piv_set]
    colorings_dart: list[list[int]] = []
    for free_vals in product(range(p), repeat=len(free_cols)):
        vec = [0] * num_vars
        for fc, fv in zip(free_cols, free_vals):
            vec[fc] = fv
        for rr, pc in enumerate(piv_cols):
            s_val = 0
            for fc in free_cols:
                if rref[rr][fc]:
                    s_val += rref[rr][fc] * vec[fc]
            vec[pc] = (-s_val) % p
        colorings_dart.append([vec[dart_var[hd]] for hd in range(4 * n)])

    assert len(colorings_dart) == known_col_count, (
        f"Coloring count mismatch: {len(colorings_dart)} vs {known_col_count}"
    )

    state_sums_pos: list[str] = []
    state_sums_mir: list[str] = []
    free_ori_mult = 1 << d.free_loops
    for base_ddir in base_oris:
        for sign_mult in (1, -1):
            init_dir_cl = [sign_mult * base_ddir[4 * v_old + s] for v_old in classical_vs for s in range(4)]
            counts: Counter[int] = Counter()
            counts_mir: Counter[int] = Counter()
            for col_d in colorings_dart:
                init_col_cl = [col_d[4 * v_old + s] for v_old in classical_vs for s in range(4)]
                wp = _evaluate_classical_movie_cocycle(mov_plus, init_col_cl, init_dir_cl, p)
                wm = _evaluate_classical_movie_cocycle(mov_minus, init_col_cl, init_dir_cl, p)
                w_net = (wp - wm) % p
                counts[w_net] += 1
                counts_mir[(-w_net) % p] += 1
            state_sums_pos.extend([_format_zp_state_sum(counts, p)] * free_ori_mult)
            state_sums_mir.extend([_format_zp_state_sum(counts_mir, p)] * free_ori_mult)

    tup_pos = tuple(sorted(state_sums_pos))
    if not up_to_mirror:
        return tup_pos
    tup_mir = tuple(sorted(state_sums_mir))
    return min(tup_pos, tup_mir)


