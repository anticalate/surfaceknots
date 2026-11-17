"""Independent mathematical audit of flagship ch=12 2-knot 12_2 (E_1 = (3, t^2+1))."""
from __future__ import annotations

import json
from collections import deque
from pathlib import Path
import sympy as sp

REPO = Path("/afs/csail.mit.edu/u/e/ejli/.gemini/antigravity-cli/worktrees/1dadd8f0-f304-455a-b682-76a831212ce6/subagent-4D-Topology---Mystery-Knot-Analyst-self-69be3150")
import sys
sys.path.insert(0, str(REPO / "code" / "src"))

from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.unlink import certify_classical_unlink


# ---------------------------------------------------------------------------
# 1. Independent from-scratch edge Wirtinger presentation directly from EPD
# ---------------------------------------------------------------------------
def independent_edge_wirtinger(d: MarkedDiagram):
    """Build exact edge-based Wirtinger presentation directly from opp & vtype.
    No reliance on surfaceknots.invariants!
    """
    n = d.num_vertices
    eid = {}
    k = 0
    for x in range(4 * n):
        if x not in eid:
            eid[x] = (k, 1)
            eid[d.opp[x]] = (k, -1)
            k += 1
    # Each generator g_0 .. g_{k-1} corresponds to meridian of directed edge;
    # determine consistent orientation signs in H_1 = Z first!
    rels_words = []  # words as tuples of signed 1-based ints (+(e+1) or -(e+1))

    def lit(dart: int) -> int:
        e, s = eid[dart]
        return (e + 1) if s == 1 else -(e + 1)

    for v in range(n):
        vt = d.vtype[v]
        b = 4 * v
        if vt in ("Y", "Z"):
            # Saddle: around vertex, signs alternate, all four arms have identical meridian up to orientation:
            # g(b+s) * g(b+s+1) = 1
            for s in range(3):
                rels_words.append((lit(b + s), lit(b + s + 1)))
        elif vt in ("X0", "X1"):
            u = 0 if vt == "X0" else 1
            o = 1 if vt == "X0" else 0
            # Over-strand passes straight through: g(b+o) * g(b+o+2) = 1
            rels_words.append((lit(b + o), lit(b + o + 2)))
            # Under-strand conjugation: g(b+u+2) * O * g(b+u) * O^-1 = 1
            eps = 1 if ((o - u) % 4) == 1 else -1
            O = lit(b + o) if eps > 0 else -lit(b + o)
            rels_words.append((lit(b + u + 2), O, lit(b + u), -O))
        else:
            raise ValueError(f"Unknown vertex type {vt}")
    return k, rels_words


# ---------------------------------------------------------------------------
# 2. Fast Union-Find HLT Todd-Coxeter Coset Enumerator (pure Python, exact)
# ---------------------------------------------------------------------------
def todd_coxeter_order(num_gens: int, relators: list[tuple[int, ...]], sub_gens: list[tuple[int, ...]] = (), max_cosets: int = 200000) -> int:
    """Exact HLT Todd-Coxeter coset enumeration with union-find coincidence processing."""
    # Columns: 0..2*num_gens-1 where 2*i is +(i+1) and 2*i+1 is -(i+1)
    ncols = 2 * num_gens

    def col_of(lit: int) -> int:
        g = abs(lit) - 1
        return 2 * g if lit > 0 else 2 * g + 1

    inv_col = [c ^ 1 for c in range(ncols)]
    rel_cols = [tuple(col_of(x) for x in r) for r in relators if r]
    sub_cols = [tuple(col_of(x) for x in w) for w in sub_gens if w]

    parent = [0]
    table = [[-1] * ncols]
    ded_q = deque()

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def new_coset() -> int:
        c = len(parent)
        if c >= max_cosets:
            raise RuntimeError(f"Exceeded max_cosets={max_cosets}")
        parent.append(c)
        table.append([-1] * ncols)
        return c

    def merge(a: int, b: int) -> None:
        q = deque([(a, b)])
        while q:
            u, v = q.popleft()
            u = find(u)
            v = find(v)
            if u == v:
                continue
            if u > v:
                u, v = v, u
            parent[v] = u
            row_u = table[u]
            row_v = table[v]
            for col in range(ncols):
                tv = row_v[col]
                if tv == -1:
                    continue
                row_v[col] = -1
                tv = find(tv)
                # Fix inverse pointer from tv
                ic = inv_col[col]
                if table[tv][ic] == v:
                    table[tv][ic] = -1
                tu = row_u[col]
                if tu == -1:
                    row_u[col] = tv
                    if table[tv][ic] == -1:
                        table[tv][ic] = u
                    else:
                        q.append((table[tv][ic], u))
                else:
                    q.append((tu, tv))

    def define(c: int, col: int) -> int:
        c = find(c)
        if table[c][col] != -1:
            return find(table[c][col])
        d = new_coset()
        table[c][col] = d
        table[d][inv_col[col]] = c
        return d

    def trace_word(start_c: int, w_cols: tuple[int, ...]) -> None:
        n = len(w_cols)
        if n == 0:
            return
        i = 0
        cur_f = find(start_c)
        while i < n:
            nxt = table[cur_f][w_cols[i]]
            if nxt == -1:
                break
            cur_f = find(nxt)
            i += 1
        if i == n:
            if cur_f != find(start_c):
                merge(cur_f, start_c)
            return
        j = n - 1
        cur_b = find(start_c)
        while j > i:
            nxt = table[cur_b][inv_col[w_cols[j]]]
            if nxt == -1:
                break
            cur_b = find(nxt)
            j -= 1
        while i < j:
            cur_f = define(cur_f, w_cols[i])
            i += 1
        # Now i == j: connect cur_f --w_cols[i]--> cur_b
        col = w_cols[i]
        ic = inv_col[col]
        cur_f = find(cur_f)
        cur_b = find(cur_b)
        if table[cur_f][col] != -1:
            merge(table[cur_f][col], cur_b)
        elif table[cur_b][ic] != -1:
            merge(cur_f, table[cur_b][ic])
        else:
            table[cur_f][col] = cur_b
            table[cur_b][ic] = cur_f

    # Subgroup generators stabilize coset 0
    for w in sub_cols:
        trace_word(0, w)

    c_idx = 0
    while c_idx < len(parent):
        if find(c_idx) == c_idx:
            for r in rel_cols:
                if find(c_idx) != c_idx:
                    break
                trace_word(c_idx, r)
            if find(c_idx) == c_idx:
                for col in range(ncols):
                    if find(c_idx) != c_idx:
                        break
                    if table[c_idx][col] == -1:
                        define(c_idx, col)
        c_idx += 1

    live = [i for i in range(len(parent)) if find(i) == i]
    return len(live)


# ---------------------------------------------------------------------------
# 3. Independent Fox Calculus over Z[t, t^-1]
# ---------------------------------------------------------------------------
def independent_fox_e1(num_gens: int, rels: list[tuple[int, ...]]):
    # Find abelianization map G -> Z (t^+-1) via Smith Normal Form
    rows = []
    for r in rels:
        row = [0] * num_gens
        for lit in r:
            row[abs(lit) - 1] += 1 if lit > 0 else -1
        if any(row):
            rows.append(row)
    M = sp.Matrix(rows)
    ns = M.nullspace()
    assert len(ns) == 1, f"Expected H_1 rank 1, got {len(ns)}"
    vec = ns[0]
    lcm_den = sp.lcm([x.q for x in vec])
    ivec = [int(x * lcm_den) for x in vec]
    g_all = sp.gcd(ivec)
    ivec = [int(x // g_all) for x in ivec]
    if ivec[0] < 0:
        ivec = [-x for x in ivec]
    assert all(abs(x) == 1 for x in ivec), f"Expected all meridians +-1, got {ivec}"

    t = sp.Symbol("t")
    # Compute Fox Jacobian J_{r, j} in Z[t] (cleared of negative powers of t)
    jac_rows = []
    for r in rels:
        derivs = [sp.Integer(0)] * num_gens
        pre = sp.Integer(1)
        for lit in r:
            g = abs(lit) - 1
            img_g = t ** ivec[g]
            if lit > 0:
                derivs[g] += pre
                pre = sp.expand(pre * img_g)
            else:
                pre = sp.expand(pre / img_g)
                derivs[g] -= pre
        jac_rows.append(derivs)
    # Drop column 0 (since x_1 maps to t^1, E_1 is generated by (n-1)x(n-1) minors of remaining columns!)
    J_sub = sp.Matrix([[jac_rows[i][j] for j in range(1, num_gens)] for i in range(len(rels))])
    return ivec, J_sub


def main():
    ch12 = json.loads((REPO / "code" / "data" / "expanded_ch12.json").read_text())
    e_12_2 = next(e for e in ch12["entries"] if e["label"] == "12_2")
    epds = [e_12_2["epd"]] + list(e_12_2["alternate_epds"])
    print(f"Found {len(epds)} diagrams in 12_2 packet.")

    for idx, epd in enumerate(epds):
        d = MarkedDiagram.from_epd(epd)
        c_pos = certify_classical_unlink(d.resolve(1))
        c_neg = certify_classical_unlink(d.resolve(-1))
        topo = d.surface_topology()
        print(f"\n--- Diagram 12_2.{chr(ord('a')+idx)} ---")
        print(f"EPD: {epd}")
        print(f"ch={d.ch_index}, (c,m)=({d.num_classical},{d.num_marked}), L_+ unlink={c_pos.is_unlink}, L_- unlink={c_neg.is_unlink}")
        print(f"Topology: {[c.homeomorphism_type for c in topo.components]}, chi={topo.euler_characteristic}, e={[c.normal_euler for c in topo.components]}")

        k, rels = independent_edge_wirtinger(d)
        ivec, J_sub = independent_fox_e1(k, rels)
        print(f"Raw edge Wirtinger: {k} generators, {len(rels)} relators")
        # Simplify J_sub over Q(t) and Z/pZ[t] via row operations? Better: reduce J_sub column by column over Z[t, t^-1]
        #Wait: let's also compute from 2-generator simplified presentation and check J_sub SNF mod primes!
        t = sp.Symbol("t")
        # Clear denominators in J_sub (multiply each entry by t^10 and takePoly over ZZ)
        J_poly = J_sub.applyfunc(lambda ex: sp.expand(ex * t**10))
        for p in (2, 3, 5, 7):
            # Rank / invariant factors over GF(p)[t]!
            # Let us do Euclidean column/row reduction over GF(p)[t]:
            M_p = [[sp.Poly(J_poly[r, c], t, modulus=p) for c in range(J_poly.cols)] for r in range(J_poly.rows)]
            # Simple Smith Normal Form over GF(p)[t]:
            nrows, ncols = len(M_p), len(M_p[0])
            diag_p = []
            r0 = 0
            for c0 in range(ncols):
                # Pivot search + Euclidean loop
                while True:
                    nonzeros = [(r, c, M_p[r][c].degree()) for r in range(r0, nrows) for c in range(r0, ncols) if not M_p[r][c].is_zero]
                    if not nonzeros:
                        break
                    pr, pc, _ = min(nonzeros, key=lambda x: x[2])
                    if pr != r0:
                        M_p[r0], M_p[pr] = M_p[pr], M_p[r0]
                    if pc != r0:
                        for rr in range(nrows):
                            M_p[rr][r0], M_p[rr][pc] = M_p[rr][pc], M_p[rr][r0]
                    piv = M_p[r0][r0]
                    changed = False
                    for r in range(r0 + 1, nrows):
                        if not M_p[r][r0].is_zero:
                            q, rem = sp.div(M_p[r][r0], piv, t, domain=sp.GF(p))
                            for cc in range(r0, ncols):
                                M_p[r][cc] = sp.Poly(M_p[r][cc].as_expr() - q.as_expr() * M_p[r0][cc].as_expr(), t, modulus=p)
                            if not M_p[r][r0].is_zero:
                                changed = True
                    for c in range(r0 + 1, ncols):
                        if not M_p[r0][c].is_zero:
                            q, rem = sp.div(M_p[r0][c], piv, t, domain=sp.GF(p))
                            for rr in range(r0, nrows):
                                M_p[rr][c] = sp.Poly(M_p[rr][c].as_expr() - q.as_expr() * M_p[rr][r0].as_expr(), t, modulus=p)
                            if not M_p[r0][c].is_zero:
                                changed = True
                    if not changed:
                        # Check divisibility of subblock
                        bad = None
                        for r in range(r0 + 1, nrows):
                            for c in range(r0 + 1, ncols):
                                _, rem = sp.div(M_p[r][c], piv, t, domain=sp.GF(p))
                                if not rem.is_zero:
                                    bad = (r, c)
                                    break
                            if bad:
                                break
                        if bad:
                            br, bc = bad
                            for cc in range(r0, ncols):
                                M_p[r0][cc] = sp.Poly(M_p[r0][cc].as_expr() + M_p[br][cc].as_expr(), t, modulus=p)
                        else:
                            # Make monic and strip t powers
                            expr = piv.monic().as_expr()
                            while expr != 0 and sp.rem(expr, t, domain=sp.GF(p)) == 0:
                                expr = sp.quo(expr, t, domain=sp.GF(p))
                            diag_p.append(expr)
                            r0 += 1
                            break
                    if r0 >= min(nrows, ncols):
                        break
                if r0 >= min(nrows, ncols):
                    break
            print(f"  SNF over GF({p})[t^+-1]: len={len(diag_p)}/{k-1}, non-units={[str(x) for x in diag_p if x != 1]}")

        # Collapse 2-letter arc equalities (u, +-v) deterministically (passing from 24 darts/edges to unbroken sheets)
        # while keeping generator 1 fixed as our base meridian!
        arc_parent = {g: (g, 1) for g in range(1, k + 1)}

        def find_arc(g: int) -> tuple[int, int]:
            cur, sgn = g, 1
            while arc_parent[cur][0] != cur:
                p, ps = arc_parent[cur]
                sgn *= ps
                cur = p
            return cur, sgn

        crossing_rels = []
        for r in rels:
            if len(r) == 2:
                u, v = r
                ru, su = find_arc(abs(u))
                if u < 0:
                    su = -su
                rv, sv = find_arc(abs(v))
                if v < 0:
                    sv = -sv
                # relation (su * ru) * (sv * rv) = 1 => rv = (-su*sv) * ru
                rel_sgn = -su * sv
                if ru != rv:
                    if ru > rv:
                        ru, rv = rv, ru
                    arc_parent[rv] = (ru, rel_sgn)
                else:
                    assert rel_sgn == 1
            else:
                crossing_rels.append(r)

        reps = sorted({find_arc(g)[0] for g in range(1, k + 1)})
        assert reps[0] == 1
        rep_idx = {r_g: idx_g + 1 for idx_g, r_g in enumerate(reps)}
        sheet_rels = []
        for r in crossing_rels:
            w = []
            for lit in r:
                rg, sg = find_arc(abs(lit))
                if lit < 0:
                    sg = -sg
                new_g = rep_idx[rg]
                w.append(new_g if sg > 0 else -new_g)
            sheet_rels.append(tuple(w))
        print(f"  Exact sheet Wirtinger: {len(reps)} sheet generators, {len(sheet_rels)} crossing relators: {sheet_rels}")

        # Verify raw 24-gen matches sheet Wirtinger for m=2..7, and compute m=2..8 on exact sheet Wirtinger!
        raw_orders_2_7 = {m: todd_coxeter_order(k, rels + [(1,) * m], max_cosets=500000) for m in range(2, 8)}
        orders = {m: todd_coxeter_order(len(reps), sheet_rels + [(1,) * m], max_cosets=500000) for m in range(2, 9)}
        assert all(raw_orders_2_7[m] == orders[m] for m in range(2, 8))
        print(f"  |G / <<a^m>>| (m=2..7) RAW 24-gen verified == sheet Wirtinger: {raw_orders_2_7}")
        print(f"  |G / <<a^m>>| (m=2..8) exact sheet Wirtinger: {orders}")

    # Now compare against G_{4_1 / 3} = <x, z | w x = z w, (z x^-1)^3 = 1> where w = x z^-1 x^-1 z!
    print("\n--- Comparing against Figure-Eight Quotient <x, z | w_41 x = z w_41, (z x^-1)^3 = 1> ---")
    rel_41 = (2, -1, 2, 1, -2, -1, 2, -1, -2, 1)  # z x^-1 ... wait let us write exact z^-1 x z^-1 x^-1 z x z^-1 x z x^-1
    # x = 1, z = 2:
    rel_41_exact = (-2, 1, -2, -1, 2, 1, -2, 1, 2, -1)
    rel_cube = (2, -1, 2, -1, 2, -1)
    orders_41 = {m: todd_coxeter_order(2, [rel_41_exact, rel_cube, (1,) * m]) for m in range(2, 9)}
    print(f"  |G_(4_1/3) / <<x^m>>| (m=2..8): {orders_41}")

    # And compare against (Z_3 * Z_3) \rtimes_\phi Z:
    # Gens: 1 = a, 2 = p0, 3 = p1
    # Rels: p0^3 = 1, p1^3 = 1, a p0 a^-1 = p1, a p1 a^-1 = p1^-1 p0^-1 p1
    rels_semidirect = [
        (2, 2, 2),
        (3, 3, 3),
        (1, 2, -1, -3),
        (1, 3, -1, -3, 2, 3),
    ]
    orders_semi = {m: todd_coxeter_order(3, rels_semidirect + [(1,) * m]) for m in range(2, 9)}
    print(f"  |((Z_3 * Z_3) rtimes Z) / <<a^m>>| (m=2..8): {orders_semi}")


if __name__ == "__main__":
    main()
