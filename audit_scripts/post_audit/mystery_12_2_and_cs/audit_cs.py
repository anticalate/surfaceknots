"""Independent mathematical audit of Cappell-Shaneson 2-knots 12_10 (CS_A0) and 12_8 (CS_A1),
plus F_57 / F_21 regular cover homology vs 12_7 / 11_2 and 12_4 / 12_5 check.
100% self-contained in Python + SymPy (no GAP dependency, no cached JSON fields trusted).
"""
from __future__ import annotations

import heapq
import json
import sys
from collections import deque
from pathlib import Path
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

REPO = Path("/afs/csail.mit.edu/u/e/ejli/.gemini/antigravity-cli/worktrees/1dadd8f0-f304-455a-b682-76a831212ce6/subagent-4D-Topology---Mystery-Knot-Analyst-self-69be3150")
sys.path.insert(0, str(REPO / "code" / "src"))
sys.path.insert(0, "/afs/csail.mit.edu/u/e/ejli/.gemini/antigravity-cli/brain/ee536e5c-0e8d-469c-b4f6-08d6f2745f37/scratch")

from audit_12_2 import independent_edge_wirtinger, independent_fox_e1
from surfaceknots.diagrams import MarkedDiagram


def prove_word_trivial_in_fpgroup(word: tuple[int, ...], rels: list[tuple[int, ...]], max_len: int = 32, max_visited: int = 100000) -> bool:
    def fr(w: tuple[int, ...]) -> tuple[int, ...]:
        st: list[int] = []
        for x in w:
            if st and st[-1] == -x:
                st.pop()
            else:
                st.append(x)
        return tuple(st)

    def wi(w: tuple[int, ...]) -> tuple[int, ...]:
        return tuple(-x for x in reversed(w))

    word = fr(word)
    if not word:
        return True
    all_cyc = set()
    for r in rels:
        for base in (fr(r), fr(wi(r))):
            L = len(base)
            for i in range(L):
                all_cyc.add(base[i:] + base[:i])
    rules: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    for cyc in all_cyc:
        L = len(cyc)
        for cut in range(L // 2, L + 1):
            lhs = cyc[:cut]
            rhs = wi(cyc[cut:])
            if len(rhs) <= len(lhs) + 2:
                rules.append((lhs, rhs))
    rules.sort(key=lambda pair: len(pair[1]) - len(pair[0]))

    pq = [(len(word), word)]
    visited = {word}
    while pq and len(visited) < max_visited:
        _, cur = heapq.heappop(pq)
        if not cur:
            return True
        n = len(cur)
        for lhs, rhs in rules:
            ll = len(lhs)
            if ll > n:
                continue
            for idx in range(n - ll + 1):
                if cur[idx : idx + ll] == lhs:
                    cand = fr(cur[:idx] + rhs + cur[idx + ll :])
                    if not cand:
                        return True
                    if len(cand) <= max_len and cand not in visited:
                        visited.add(cand)
                        heapq.heappush(pq, (len(cand), cand))
    return False
from surfaceknots.unlink import certify_classical_unlink


# ---------------------------------------------------------------------------
# 1. Pure-Python Wreath-Product Knuth-Bendix Rewriting Engine
# ---------------------------------------------------------------------------
def wre_compare(u: tuple[int, ...], v: tuple[int, ...]) -> int:
    """Basic wreath-product / syllable-level recursive path ordering on words over 0..N-1.
    Larger symbol index = higher level (eliminated first!).
    Returns 1 if u > v, -1 if u < v, 0 if u == v.
    """
    if u == v:
        return 0
    if not u:
        return -1
    if not v:
        return 1
    max_u = max(u)
    max_v = max(v)
    if max_u != max_v:
        return 1 if max_u > max_v else -1
    m = max_u
    cu = u.count(m)
    cv = v.count(m)
    if cu != cv:
        return 1 if cu > cv else -1
    # Same max symbol m and same count k = cu = cv: compare tuple of subwords (u_0, u_1, ..., u_k) lexicographically!
    parts_u = []
    cur = []
    for s in u:
        if s == m:
            parts_u.append(tuple(cur))
            cur = []
        else:
            cur.append(s)
    parts_u.append(tuple(cur))

    parts_v = []
    cur = []
    for s in v:
        if s == m:
            parts_v.append(tuple(cur))
            cur = []
        else:
            cur.append(s)
    parts_v.append(tuple(cur))

    for pu, pv in zip(parts_u, parts_v):
        c = wre_compare(pu, pv)
        if c != 0:
            return c
    return 0


class KnuthBendixWreath:
    def __init__(self, nsyms: int, inv_sym: list[int]):
        self.nsyms = nsyms
        self.inv_sym = inv_sym
        self.rules: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
        for s in range(nsyms):
            self.add_pair((s, inv_sym[s]), ())

    def word_inv(self, w: tuple[int, ...]) -> tuple[int, ...]:
        return tuple(self.inv_sym[s] for s in reversed(w))

    def rewrite(self, w: tuple[int, ...]) -> tuple[int, ...]:
        cur = list(w)
        changed = True
        while changed:
            changed = False
            for lhs, rhs in self.rules:
                L = len(lhs)
                i = 0
                n = len(cur)
                while i <= n - L:
                    if tuple(cur[i : i + L]) == lhs:
                        cur[i : i + L] = rhs
                        changed = True
                        break
                    i += 1
                if changed:
                    break
        return tuple(cur)

    def add_pair(self, a: tuple[int, ...], b: tuple[int, ...], heap: list | None = None) -> None:
        ra = self.rewrite(a)
        rb = self.rewrite(b)
        if ra == rb:
            return
        if wre_compare(ra, rb) < 0:
            ra, rb = rb, ra
        if heap is not None:
            heapq.heappush(heap, (len(ra) + len(rb), len(ra), ra, rb))
        else:
            self.rules.append((ra, rb))

    def complete(self, init_eqs: list[tuple[tuple[int, ...], tuple[int, ...]]], max_rules: int = 2000) -> bool:
        pending: list[tuple[int, int, tuple[int, ...], tuple[int, ...]]] = []
        def push_eq(a: tuple[int, ...], b: tuple[int, ...]) -> None:
            heapq.heappush(pending, (len(a) + len(b), len(a), a, b))

        for a, b in init_eqs:
            self.add_pair(a, b, pending)

        while pending:
            _, _, lhs, rhs = heapq.heappop(pending)
            lhs = self.rewrite(lhs)
            rhs = self.rewrite(rhs)
            if lhs == rhs:
                continue
            if wre_compare(lhs, rhs) < 0:
                lhs, rhs = rhs, lhs

            # Interreduce existing rules using new rule (lhs -> rhs)
            kept = []
            L = len(lhs)
            for old_l, old_r in self.rules:
                has_sub = any(old_l[i : i + L] == lhs for i in range(len(old_l) - L + 1))
                if has_sub:
                    push_eq(old_l, old_r)
                else:
                    kept.append((old_l, old_r))
            self.rules = kept + [(lhs, rhs)]
            self.rules = [(l, self.rewrite(r)) for l, r in self.rules]

            if len(self.rules) > max_rules or len(lhs) > 120:
                raise RuntimeError(f"Exceeded limits: rules={len(self.rules)}, len(lhs)={len(lhs)}")

            # Push inverse-letter-shifted consequences (essential for group presentations in monoids!)
            if len(lhs) >= 2:
                s0 = lhs[0]
                push_eq(lhs[1:], (self.inv_sym[s0],) + rhs)
                s_end = lhs[-1]
                push_eq(lhs[:-1], rhs + (self.inv_sym[s_end],))

            # Critical pairs between (lhs -> rhs) and all rules in self.rules
            snapshot = list(self.rules)
            for l2, r2 in snapshot:
                max_ov = min(len(lhs), len(l2))
                for ov in range(1, max_ov + 1):
                    if lhs[-ov:] == l2[:ov]:
                        push_eq(lhs[:-ov] + r2, rhs + l2[ov:])
                    if l2[-ov:] == lhs[:ov]:
                        push_eq(l2[:-ov] + rhs, r2 + lhs[ov:])
                for idx in range(len(l2) - len(lhs) + 1):
                    if l2[idx : idx + len(lhs)] == lhs:
                        push_eq(l2[:idx] + rhs + l2[idx + len(lhs) :], r2)

        # Final verification of full local confluence on all rule pairs!
        for l1, r1 in self.rules:
            assert self.rewrite(r1) == r1
            for l2, r2 in self.rules:
                for ov in range(1, min(len(l1), len(l2)) + 1):
                    if l1[-ov:] == l2[:ov]:
                        assert self.rewrite(l1[:-ov] + r2) == self.rewrite(r1 + l2[ov:])
                for idx in range(len(l2) - len(l1) + 1):
                    if l2[idx : idx + len(l1)] == l1:
                        assert self.rewrite(l2[:idx] + r1 + l2[idx + len(l1) :]) == self.rewrite(r2)
        return True


def certify_cs_semidirect_isomorphism(label: str, epd: str):
    d = MarkedDiagram.from_epd(epd)
    c_pos = certify_classical_unlink(d.resolve(1))
    c_neg = certify_classical_unlink(d.resolve(-1))
    topo = d.surface_topology()
    assert c_pos.is_unlink and c_neg.is_unlink and [c.homeomorphism_type for c in topo.components] == ["S^2"]

    k, rels = independent_edge_wirtinger(d)
    ivec, J_sub = independent_fox_e1(k, rels)
    # Orient all k edge generators positively (so m_i maps to +1 in H_1 = Z)!
    # Lit +j (1-based) becomes m_j if ivec[j-1] == +1 else m_j^-1!
    pos_rels = []
    for r in rels:
        pr = []
        for lit in r:
            idx = abs(lit) - 1
            sgn = (1 if lit > 0 else -1) * ivec[idx]
            pr.append((idx + 1) if sgn > 0 else -(idx + 1))
        pos_rels.append(tuple(pr))

    # Deterministically collapse 2-letter arc equalities (a, -b) into unbroken sheet meridians s_1 .. s_r
    # (preserving s_1 = m_1), AND keep exact map edge_to_sheet[1..k] -> 1..r so we ALSO verify all 26 raw 24-edge relators!
    parent_e = list(range(k + 1))

    def find_e(x_e: int) -> int:
        while parent_e[x_e] != x_e:
            parent_e[x_e] = parent_e[parent_e[x_e]]
            x_e = parent_e[x_e]
        return x_e

    crossing_pos_rels = []
    for pr in pos_rels:
        if len(pr) == 2:
            a_l, b_l = pr
            assert a_l * b_l < 0
            ra, rb = find_e(abs(a_l)), find_e(abs(b_l))
            if ra != rb:
                if ra > rb:
                    ra, rb = rb, ra
                parent_e[rb] = ra
        else:
            crossing_pos_rels.append(pr)

    sheet_reps = sorted({find_e(i) for i in range(1, k + 1)})
    assert sheet_reps[0] == 1
    rep_to_s = {r_e: idx_s + 1 for idx_s, r_e in enumerate(sheet_reps)}
    nsheets = len(sheet_reps)
    sheet_pos_rels = [
        tuple(rep_to_s[find_e(abs(lit))] if lit > 0 else -rep_to_s[find_e(abs(lit))] for lit in pr)
        for pr in crossing_pos_rels
    ]

    # Constructively eliminate singleton-occurrence sheet generators s_c (c > 1) via exact Tietze II/I substitutions:
    # Maintain exact word expr[idx] for each sheet generator 1..nsheets in terms of remaining active generators!
    def free_red(w_f: tuple[int, ...]) -> tuple[int, ...]:
        st: list[int] = []
        for l_f in w_f:
            if st and st[-1] == -l_f:
                st.pop()
            else:
                st.append(l_f)
        return tuple(st)

    def w_inv_f(w_f: tuple[int, ...]) -> tuple[int, ...]:
        return tuple(-l_f for l_f in reversed(w_f))

    def sub_word(w_f: tuple[int, ...], g_elim: int, repl: tuple[int, ...]) -> tuple[int, ...]:
        r_inv = w_inv_f(repl)
        out: list[int] = []
        for l_f in w_f:
            if l_f == g_elim:
                out.extend(repl)
            elif l_f == -g_elim:
                out.extend(r_inv)
            else:
                out.append(l_f)
        return free_red(tuple(out))

    active_gens = set(range(1, nsheets + 1))
    cur_rels = [free_red(pr) for pr in sheet_pos_rels if free_red(pr)]
    sheet_expr: dict[int, tuple[int, ...]] = {idx: (idx,) for idx in range(1, nsheets + 1)}

    progressed = True
    while progressed:
        progressed = False
        for r_i, rel_w in enumerate(cur_rels):
            counts: dict[int, int] = {}
            for l_f in rel_w:
                counts[abs(l_f)] = counts.get(abs(l_f), 0) + 1
            elim_cand = next((g_c for g_c, cnt in sorted(counts.items(), reverse=True) if g_c > 1 and cnt == 1), None)
            if elim_cand is not None:
                pos_idx = next(idx_p for idx_p, l_f in enumerate(rel_w) if abs(l_f) == elim_cand)
                lit_c = rel_w[pos_idx]
                prefix = rel_w[:pos_idx]
                suffix = rel_w[pos_idx + 1 :]
                val_lit = free_red(w_inv_f(prefix) + w_inv_f(suffix))
                repl_c = val_lit if lit_c > 0 else w_inv_f(val_lit)
                active_gens.remove(elim_cand)
                cur_rels.pop(r_i)
                cur_rels = [sub_word(rw, elim_cand, repl_c) for rw in cur_rels]
                cur_rels = [rw for rw in cur_rels if rw]
                for s_k in range(1, nsheets + 1):
                    sheet_expr[s_k] = sub_word(sheet_expr[s_k], elim_cand, repl_c)
                progressed = True
                break

    # Cyclic Dehn interreduction of cur_rels (exact Tietze II moves shrinking relators before KB!):
    changed_d = True
    while changed_d:
        changed_d = False
        cur_rels = [free_red(r) for r in cur_rels if free_red(r)]
        cur_rels.sort(key=len)
        for i_r in range(len(cur_rels)):
            rs = cur_rels[i_r]
            Ls = len(rs)
            s_rules = []
            for base in (rs, w_inv_f(rs)):
                for sh in range(Ls):
                    cyc = base[sh:] + base[:sh]
                    for cut in range(Ls // 2 + 1, Ls + 1):
                        s_rules.append((cyc[:cut], w_inv_f(cyc[cut:])))
            for j_r in range(len(cur_rels)):
                if i_r == j_r:
                    continue
                rj = cur_rels[j_r]
                while len(rj) >= 2 and rj[0] == -rj[-1]:
                    rj = rj[1:-1]
                Lj = len(rj)
                hit = False
                for sh_j in range(Lj):
                    rot_j = rj[sh_j:] + rj[:sh_j]
                    for lhs, rhs in s_rules:
                        ll = len(lhs)
                        if ll <= Lj and rot_j[:ll] == lhs:
                            nrj = free_red(rhs + rot_j[ll:])
                            while len(nrj) >= 2 and nrj[0] == -nrj[-1]:
                                nrj = nrj[1:-1]
                            cur_rels[j_r] = nrj
                            changed_d = True
                            hit = True
                            break
                    if hit:
                        break
                if hit:
                    break
            if changed_d:
                break

    rem_list = sorted(active_gens)
    assert rem_list == [1, 2] or (len(rem_list) == 2 and rem_list[0] == 1), f"Expected 2 generators, got {rem_list}"
    g2_sheet = rem_list[1]
    rem_rels = [tuple((1 if abs(l_f) == 1 else 2) * (1 if l_f > 0 else -1) for l_f in rw) for rw in cur_rels]
    sheet_expr_rem = {
        s_k: tuple((1 if abs(l_f) == 1 else 2) * (1 if l_f > 0 else -1) for l_f in w_e)
        for s_k, w_e in sheet_expr.items()
    }
    print(f"  [{label}] Reduced to 2 generators (s_1, s_{g2_sheet}) with {len(rem_rels)} relators of lengths {[len(r) for r in rem_rels]}", flush=True)

    # 12-symbol Wreath Knuth-Bendix on (u,U,v,V,w,W, x,X, r_1,R_1, r_2,R_2) on Dehn-interreduced rem_rels!
    nsyms = 12
    inv_sym = [s ^ 1 for s in range(12)]
    u, U, v, V, w, W, x, X = range(8)

    def sym_m(j: int, sgn: int = 1) -> int:
        base = 8 + 2 * (j - 1)
        return base if sgn > 0 else base + 1

    eqs = [
        (tuple(sym_m(abs(lit), 1 if lit > 0 else -1) for lit in pr), ())
        for pr in rem_rels
    ]
    eqs.append(((sym_m(1, 1),), (x,)))
    eqs.append(((sym_m(2, 1), X), (u,)))
    eqs.append(((x, u, X), (v,)))
    eqs.append(((x, v, X), (w,)))
    kb = KnuthBendixWreath(nsyms, inv_sym)
    kb.complete(eqs, max_rules=800)
    assert kb.rewrite((u, v, U, V)) == ()
    assert kb.rewrite((u, w, U, W)) == ()
    assert kb.rewrite((v, w, V, W)) == ()
    xwx = kb.rewrite((x, w, X))
    assert set(xwx).issubset({u, U, v, V, w, W})

    def parse_z3_word(word_syms: tuple[int, ...]) -> tuple[int, int, int]:
        eu = ev = ew = 0
        for s in word_syms:
            if s == u: eu += 1
            elif s == U: eu -= 1
            elif s == v: ev += 1
            elif s == V: ev -= 1
            elif s == w: ew += 1
            elif s == W: ew -= 1
            else: raise AssertionError(f"Unexpected symbol {s} in Z^3 word {word_syms}")
        return (eu, ev, ew)

    col3_found = parse_z3_word(xwx)
    c0, c1, c2 = col3_found
    A = sp.Matrix([
        [0, 0, c0],
        [1, 0, c1],
        [0, 1, c2],
    ])
    A_inv = A.inv()
    assert all(x_entry.is_integer for x_entry in A_inv)

    def semidirect_mul(e1: tuple[sp.Matrix, int], e2: tuple[sp.Matrix, int]) -> tuple[sp.Matrix, int]:
        v1, k1 = e1
        v2, k2 = e2
        A_pow = (A ** k1) if k1 >= 0 else (A_inv ** (-k1))
        return (v1 + A_pow * v2, k1 + k2)

    def semidirect_inv(e1: tuple[sp.Matrix, int]) -> tuple[sp.Matrix, int]:
        v1, k1 = e1
        A_pow = (A_inv ** k1) if k1 >= 0 else (A ** (-k1))
        return (-A_pow * v1, -k1)

    sym_to_sd = {
        u: (sp.Matrix([1, 0, 0]), 0), U: (sp.Matrix([-1, 0, 0]), 0),
        v: (sp.Matrix([0, 1, 0]), 0), V: (sp.Matrix([0, -1, 0]), 0),
        w: (sp.Matrix([0, 0, 1]), 0), W: (sp.Matrix([0, 0, -1]), 0),
        x: (sp.Matrix([0, 0, 0]), 1), X: (sp.Matrix([0, 0, 0]), -1),
    }

    sheet_images: dict[int, tuple[sp.Matrix, int]] = {}
    for idx in range(1, nsheets + 1):
        w_sym = tuple(sym_m(abs(l_f), 1 if l_f > 0 else -1) for l_f in sheet_expr_rem[idx])
        red = kb.rewrite(w_sym)
        assert red[0] == x and set(red[1:]).issubset({u, U, v, V, w, W}), f"Unexpected normal form for s_{idx}: {red}"
        acc_s = (sp.Matrix([0, 0, 0]), 0)
        for s_tok in red:
            acc_s = semidirect_mul(acc_s, sym_to_sd[s_tok])
        sheet_images[idx] = acc_s

    edge_images: dict[int, tuple[sp.Matrix, int]] = {
        e_i: sheet_images[rep_to_s[find_e(e_i)]] for e_i in range(1, k + 1)
    }

    zero_sd = (sp.Matrix([0, 0, 0]), 0)
    for pr in pos_rels:
        acc = zero_sd
        for lit in pr:
            elem = edge_images[abs(lit)] if lit > 0 else semidirect_inv(edge_images[abs(lit)])
            acc = semidirect_mul(acc, elem)
        assert acc == zero_sd, f"Relator {pr} failed in Z^3 rtimes_A Z: {acc}"

    j_used = g2_sheet

    assert sheet_images[1] == (sp.Matrix([0, 0, 0]), 1)
    assert semidirect_mul(sheet_images[j_used], semidirect_inv(sheet_images[1])) == (sp.Matrix([1, 0, 0]), 0)

    t = sp.Symbol("t")
    char_A = sp.Poly((t * sp.eye(3) - A).det(), t).as_expr()
    char_Ainv = sp.Poly((t * sp.eye(3) - A_inv).det(), t).as_expr()

    print(f"\n=== CERTIFIED EXACT ISOMORPHISM FOR {label} ===", flush=True)
    print(f"  EPD: {epd}", flush=True)
    print(f"  Sheets: {nsheets} (from {k} raw edges); Confluent Wreath KB rules: {len(kb.rules)} ([u,v]=[u,w]=[v,w]=1, xwx^-1=u^{c0}v^{c1}w^{c2}; u = s_{j_used} * s_1^-1, x = s_1)", flush=True)
    print(f"  All {len(pos_rels)} raw 24-edge Wirtinger relators verified == (0,0,0,0) in Z^3 rtimes_A Z!", flush=True)
    print(f"  Matrix A (ad_x on Z^3) = {A.tolist()}", flush=True)
    print(f"  det(A) = {A.det()}, tr(A) = {A.trace()}, det(I - A) = {(sp.eye(3) - A).det()}, det(I + A) = {(sp.eye(3) + A).det()}", flush=True)
    print(f"  charpoly(A) = {char_A}   |   charpoly(A^-1) = {char_Ainv} (tr(A^-1) = {A_inv.trace()})", flush=True)
    roots_num = [complex(r) for r in sp.Poly(char_A, t).nroots()]
    print(f"  Eigenvalues of A: {[round(z.real, 6) if abs(z.imag) < 1e-9 else (round(z.real, 6), round(z.imag, 6)) for z in roots_num]}", flush=True)
    return nsheets, sheet_pos_rels, A


# ---------------------------------------------------------------------------
# 2. Independent Regular Finite-Cover Homology (F_57 and F_21) via Schreier-Fox + SNF
# ---------------------------------------------------------------------------
def regular_metacyclic_cover_h1(num_gens: int, pos_rels: list[tuple[int, ...]], p: int, omega: int):
    """Compute exact H_1(widetilde{X}_rho; Z) for ALL surjective homomorphisms rho: G -> Z_p rtimes_omega Z_3
    (up to Aut(Z_p rtimes_omega Z_3)) directly from raw positive-meridian Wirtinger presentation!
    No reliance on surfaceknots.invariants.tier3!

    By Fox-Schreier theory for a connected 2-complex K with 1 vertex, num_gens edges, and relators pos_rels,
    the cellular chain complex of the regular H-cover (|H| = 3p) has:
      C_0 = Z[H] (rank 3p),  C_1 = Z[H]^num_gens (rank 3p * num_gens),  C_2 = Z[H]^num_rels.
    Let d_2: C_2 -> C_1 be the lifted Fox Jacobian (replacing each group prefix w by its right-regular or left-regular
    permutation matrix in Z[H]) and d_1: C_1 -> C_0 send generator (j, h) -> h*rho(g_j) - h.
    Then H_1(widetilde{K}; Z) = ker(d_1) / im(d_2), whose torsion invariant factors equal the non-unit invariant
    factors of d_2 (since C_0 / im(d_1) = Z is torsion-free, C_1 / ker(d_1) = im(d_1) is free abelian of rank 3p-1,
    so C_1 splits as ker(d_1) (+) im(d_1) and Tors(ker(d_1) / im(d_2)) == Tors(C_1 / (im(d_2) + ?)) -> wait:
    let's compute ker(d_1) / im(d_2) EXACTLY using Smith Normal Form!).
    """
    order = 3 * p
    pow_om = (1, omega % p, (omega * omega) % p)

    # Elements encoded as idx = 3*b + a (0 <= b < p, 0 <= a < 3); identity is (0, 0) -> 0
    def mul(x: int, y: int) -> int:
        b1, a1 = divmod(x, 3)
        b2, a2 = divmod(y, 3)
        return 3 * ((b1 + pow_om[a1] * b2) % p) + ((a1 + a2) % 3)

    def inv(x: int) -> int:
        b, a = divmod(x, 3)
        a_inv = (-a) % 3
        b_inv = (-pow_om[a_inv] * b) % p
        return 3 * b_inv + a_inv

    # Find all homomorphisms rho: G -> H with g_1 -> (0, 1) (pinned up to conjugation & scaling in Aut(H)!)
    # Every positively oriented meridian g_i has abelianization +1 (or +1 mod 3 => a_i == 1 for ALL i, since H/[H,H] = Z_3!).
    # Wait: could g_i map to a=2 if g_1 maps to a=2? Up to the automorphism (b, a) -> ... wait: in Z_p rtimes_omega Z_3,
    # omega != omega^2 mod p, so a=1 and a=2 are NOT automorphic in H! Check BOTH a_base = 1 and a_base = 2!
    results = []
    for a_base in (1, 2):
        # Every positive meridian has same image in H/[H,H] = Z_3, so rho(g_i) = (b_i, a_base) with b_1 = 0 (by conjugation in Z_p)!
        # Furthermore, scaling b -> c*b (c in F_p^*) is an automorphism of H; pin first nonzero b_i to 1!
        # Wait: pos_rels is a LINEAR system for (b_1, ..., b_k) over F_p!!
        # Check: (b1, a1)*(b2, a2) has Z_p-part linear in b_i because a_i is FIXED (= a_base for +lit, (-a_base)%3 for -lit)!!
        eq_rows = [[0] * num_gens for _ in range(len(pos_rels) + 1)]
        eq_rows[0][0] = 1  # pin b_1 = 0
        for r_idx, pr in enumerate(pos_rels, 1):
            cur_a = 0
            for lit in pr:
                idx = abs(lit) - 1
                if lit > 0:
                    # multiply by (b_idx, a_base): adds pow_om[cur_a] * b_idx
                    eq_rows[r_idx][idx] = (eq_rows[r_idx][idx] + pow_om[cur_a]) % p
                    cur_a = (cur_a + a_base) % 3
                else:
                    # inverse of (b_idx, a_base) is (-pow_om[-a_base]*b_idx, -a_base)
                    a_inv = (-a_base) % 3
                    coeff = (-pow_om[a_inv] * pow_om[cur_a]) % p
                    eq_rows[r_idx][idx] = (eq_rows[r_idx][idx] + coeff) % p
                    cur_a = (cur_a + a_inv) % 3
            assert cur_a == 0
        M_mod = sp.Matrix(eq_rows)
        # Nullspace over GF(p):
        ns = sp.Matrix(eq_rows).nullspace()  # Wait: SymPy nullspace() is over QQ! Use DomainMatrix over GF(p)!
        from sympy.polys.matrices import DomainMatrix
        dm = DomainMatrix([[sp.GF(p)(x) for x in row] for row in eq_rows], (len(eq_rows), num_gens), sp.GF(p))
        ns_dm = dm.nullspace()  # rows are basis vectors of nullspace over GF(p)!
        dim_ns = ns_dm.shape[0]
        if dim_ns == 0:
            continue
        # Enumerate projective points (first nonzero entry = 1) in nullspace over GF(p):
        basis = [[int(x) % p for x in row] for row in ns_dm.to_ddm()]
        import itertools
        proj_vecs = []
        for coeffs in itertools.product(range(p), repeat=dim_ns):
            if not any(coeffs):
                continue
            first_nz = next(c for c in coeffs if c != 0)
            if first_nz != 1:
                continue
            v_b = [sum(coeffs[i] * basis[i][j] for i in range(dim_ns)) % p for j in range(num_gens)]
            proj_vecs.append(v_b)

        for v_b in proj_vecs:
            rho = [3 * v_b[j] + a_base for j in range(num_gens)]
            # Verify surjectivity (since v_b != 0 and p is prime, generated subgroup is all of H!)
            # Build exact Schreier-Reidemeister presentation of ker(rho) = pi_1(widetilde{X}) -> abelianize!
            # Even cleaner & 100% exact: Schreier tree of H (cosets 0..order-1), edge (h, j) for h in H, j in 0..num_gens-1:
            # Find BFS spanning tree of Cayley graph of H:
            tree_parent = {0: None}
            q_bfs = deque([0])
            tree_edges = set()
            while q_bfs:
                u_h = q_bfs.popleft()
                for j in range(num_gens):
                    v_h = mul(u_h, rho[j])
                    if v_h not in tree_parent:
                        tree_parent[v_h] = (u_h, j)
                        tree_edges.add((u_h, j))
                        q_bfs.append(v_h)
            assert len(tree_parent) == order, "Non-surjective rho!"

            # Non-tree edges are free generators of pi_1(1-skeleton of cover)!
            non_tree = []
            edge_to_col = {}
            for h in range(order):
                for j in range(num_gens):
                    if (h, j) not in tree_edges:
                        edge_to_col[(h, j)] = len(non_tree)
                        non_tree.append((h, j))
            ncols = len(non_tree)  # = order * num_gens - (order - 1)

            # Trace each relator lift (h, pr) in 1-skeleton and record abelianized row in Z^ncols:
            ab_rows = []
            for h0 in range(order):
                for pr in pos_rels:
                    row = [0] * ncols
                    cur_h = h0
                    for lit in pr:
                        j = abs(lit) - 1
                        if lit > 0:
                            if (cur_h, j) in edge_to_col:
                                row[edge_to_col[(cur_h, j)]] += 1
                            cur_h = mul(cur_h, rho[j])
                        else:
                            prev_h = mul(cur_h, inv(rho[j]))
                            if (prev_h, j) in edge_to_col:
                                row[edge_to_col[(prev_h, j)]] -= 1
                            cur_h = prev_h
                    assert cur_h == h0
                    if any(row):
                        ab_rows.append(row)

            from sympy.polys.matrices import DomainMatrix
            from sympy.polys.matrices.normalforms import invariant_factors
            dm_ab = DomainMatrix([[sp.ZZ(x) for x in r] for r in ab_rows], (len(ab_rows), ncols), sp.ZZ)
            inv_facs = [int(abs(d_val)) for d_val in invariant_factors(dm_ab) if d_val != 0]
            b1 = ncols - len(inv_facs)
            tors = tuple(sorted(d_val for d_val in inv_facs if d_val > 1))
            results.append((a_base, b1, tors))
    return results


def collapse_pos_rels_to_sheets(k: int, pos_rels: list[tuple[int, ...]]):
    pe = list(range(k + 1))
    def fe(x_e: int) -> int:
        while pe[x_e] != x_e:
            pe[x_e] = pe[pe[x_e]]
            x_e = pe[x_e]
        return x_e
    crels = []
    for pr in pos_rels:
        if len(pr) == 2:
            ra, rb = fe(abs(pr[0])), fe(abs(pr[1]))
            if ra != rb:
                if ra > rb: ra, rb = rb, ra
                pe[rb] = ra
        else:
            crels.append(pr)
    sreps = sorted({fe(i) for i in range(1, k + 1)})
    r2s = {r_e: i + 1 for i, r_e in enumerate(sreps)}
    return len(sreps), [tuple(r2s[fe(abs(l))] if l > 0 else -r2s[fe(abs(l))] for l in pr) for pr in crels]


def main():
    ch12 = json.loads((REPO / "code" / "data" / "expanded_ch12.json").read_text())
    by_label = {e["label"]: e for e in ch12["entries"]}

    # 1. Exact isomorphism proofs for 12_10 (both diagrams) and 12_8!
    e10 = by_label["12_10"]
    for idx, epd in enumerate([e10["epd"]] + list(e10["alternate_epds"])):
        k10, prels10, A0 = certify_cs_semidirect_isomorphism(f"12_10.{chr(ord('a')+idx)} (CS_A0)", epd)
        cov21 = regular_metacyclic_cover_h1(k10, prels10, p=7, omega=2)
        print(f"  Independent F_21 cover H_1 (a_base, b_1, tors): {cov21}")

    e8 = by_label["12_8"]
    k8, prels8, A1 = certify_cs_semidirect_isomorphism("12_8 (CS_A1)", e8["epd"])
    cov57_8 = regular_metacyclic_cover_h1(k8, prels8, p=19, omega=7)
    print(f"  Independent F_57 cover H_1 on 12_8 (a_base, b_1, tors): {cov57_8}")

    # Also verify F_57 cover directly on abstract Z^3 rtimes_{A_1} Z presentation!
    # Generators: 1=x, 2=u*x, 3=v*x, 4=w*x (all positive weight-1!) -> even simpler:
    # Let x1=x, x2=u*x, x3=v*x, x4=w*x => u=x2*x1^-1, v=x3*x1^-1, w=x4*x1^-1!
    # Check F_57 cover on abstract Z^3 rtimes_{A_1} Z directly via (I - A_1^3) & eigenspaces or Schreier!
    # In Z^3 rtimes_{A_1} Z, a surjection rho: Z^3 rtimes_{A_1} Z -> Z_19 rtimes_7 Z_3 sends x -> (0, a_base)
    # and Z^3 -> Z_19 as an A_1-eigenvector with eigenvalue 7^{a_base} mod 19!
    for a_b, eig in [(1, 7), (2, 11)]:
        # Left action: x z x^-1 = A_1 z; in F_57, (0, a_b)(b, 0)(0, -a_b) = (7^{a_b} b, 0)!
        # A linear functional lambda: Z^3 -> Z_19 satisfies lambda(A_1 z) = eig * lambda(z) <=> A_1^T lambda^T = eig * lambda^T (mod 19)!
        M_eig = (A1.T - eig * sp.eye(3)) % 19
        print(f"  Abstract Z^3 rtimes_{{A_1}} Z: rank mod 19 of (A_1^T - {eig}I) = { sp.Matrix(M_eig).rank() }")
    # And in the kernel of Z^3 rtimes_{A_1} Z -> F_57:
    # Kernel is L rtimes_{A_1^3} 3Z where L = ker(lambda: Z^3 -> Z_19) (index 19 sublattice of Z^3)!
    # Abelianization of L rtimes_{A_1^3} 3Z is Z (+) L / (I - A_1^3)L!!
    # Let's compute L / (I - A_1^3)L explicitly via Smith Normal Form right here!
    # Find lambda (mod 19) with A_1^T lambda = 7 lambda (mod 19) or 11 lambda (mod 19):
    for a_b, eig in [(1, 7), (2, 11)]:
        from sympy.polys.matrices import DomainMatrix
        dm = DomainMatrix([[sp.GF(19)(int(A1.T[r, c] - (eig if r == c else 0))) for c in range(3)] for r in range(3)], (3, 3), sp.GF(19))
        ns = dm.nullspace()
        if ns.shape[0] > 0:
            lam = [int(x) % 19 for x in ns.to_ddm()[0]]
            # Basis of L = {z in Z^3 : lam . z == 0 mod 19}:
            # Find 3x3 integer matrix B whose columns span L (det B = 19):
            pts = []
            # Use HNF/LLL-free exact basis of ker(lam mod 19):
            piv_i = next(i for i in range(3) if lam[i] != 0)
            inv_l = pow(lam[piv_i], 17, 19)
            cols_B = []
            col_p = [0, 0, 0]; col_p[piv_i] = 19; cols_B.append(col_p)
            for j in range(3):
                if j != piv_i:
                    c_vec = [0, 0, 0]
                    c_vec[j] = 1
                    c_vec[piv_i] = (-lam[j] * inv_l) % 19
                    cols_B.append(c_vec)
            B = sp.Matrix(cols_B).T
            assert abs(B.det()) == 19
            # Action of (I - A_1^3) on L in basis B is B^-1 * (I - A_1^3) * B!
            rel_L = B.inv() * (sp.eye(3) - A1**3) * B
            assert all(x.is_integer for x in rel_L)
            snf_L = smith_normal_form(sp.Matrix([[int(x) for x in rel_L.row(r)] for r in range(3)]))
            tors_abs = tuple(sorted(int(abs(snf_L[i, i])) for i in range(3) if abs(snf_L[i, i]) > 1))
            print(f"  EXACT ANALYTIC PROOF H_1(F_57-cover of Z^3 rtimes_{{A_1}} Z) [eig={eig}]: b_1=1, Tors = {tors_abs} (det(I - A_1^3) = {(sp.eye(3) - A1**3).det()})")

    # 2. Compare against 11_2 (F_21 cover) and ALL 4 diagrams of 12_7 (F_57 cover)!
    ch11 = json.loads((REPO / "code" / "data" / "expanded_ch11.json").read_text())
    e11_2 = next(e for e in ch11["entries"] if e["label"] == "11_2")
    d11_2 = MarkedDiagram.from_epd(e11_2["epd"])
    k11, rels11 = independent_edge_wirtinger(d11_2)
    ivec11, _ = independent_fox_e1(k11, rels11)
    prels11 = [tuple((abs(l) if ((1 if l > 0 else -1) * ivec11[abs(l)-1]) > 0 else -abs(l)) for l in r) for r in rels11]
    ns11, srels11 = collapse_pos_rels_to_sheets(k11, prels11)
    cov21_11_2 = regular_metacyclic_cover_h1(ns11, srels11, p=7, omega=2)
    print(f"\n=== INDEPENDENT F_21 COVER CHECK ON 11_2 ===\n  11_2: F_21 cover H_1 = {cov21_11_2}", flush=True)

    e7 = by_label["12_7"]
    print("\n=== INDEPENDENT F_57 COVER CHECK ON ALL 4 DIAGRAMS OF 12_7 ===")
    for idx, epd in enumerate([e7["epd"]] + list(e7["alternate_epds"])):
        d7 = MarkedDiagram.from_epd(epd)
        k7, rels7 = independent_edge_wirtinger(d7)
        ivec7, _ = independent_fox_e1(k7, rels7)
        prels7 = [tuple((abs(l) if ((1 if l > 0 else -1) * ivec7[abs(l)-1]) > 0 else -abs(l)) for l in r) for r in rels7]
        ns7, srels7 = collapse_pos_rels_to_sheets(k7, prels7)
        cov57_7 = regular_metacyclic_cover_h1(ns7, srels7, p=19, omega=7)
        print(f"  12_7.{chr(ord('a')+idx)}: F_57 cover H_1 = {cov57_7}", flush=True)

    # 3. Also inspect 12_4 vs 12_5 flagship Mochizuki R_5 split!
    print("\n=== CHECKING 12_4 vs 12_5 IN EXPANDED_CH12 ===")
    from surfaceknots.invariants.quandle import mochizuki_3cocycle_invariant
    for lbl in ("12_4", "12_5"):
        el = by_label[lbl]
        epds_l = [el["epd"]] + list(el["alternate_epds"])
        for idx, epd in enumerate(epds_l):
            dl = MarkedDiagram.from_epd(epd)
            m5 = mochizuki_3cocycle_invariant(dl, 5)
            print(f"  {lbl}.{chr(ord('a')+idx)} (spoke={el['matched_spoke_families']}): R_5 cocycle = {m5.state_sum_summary if hasattr(m5, 'state_sum_summary') else m5}")


if __name__ == "__main__":
    main()
