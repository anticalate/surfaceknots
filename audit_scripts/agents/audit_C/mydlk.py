"""Independent 'q_i(ker lambda_ij)' : Rokhlin quadratic form of torus component F_i evaluated on the class
dual to the mod-2 linking functional lambda_ij (loops on F_i  |-> lk with F_j).
q(C) for a cross-section circle C (any mixed resolution of the saddles) = sum_{C' != C on F_i} lk(C,C') mod 2
(membrane = Seifert surface of C in the slice; framing obstruction 0)."""
from itertools import combinations

def analyse(d, ci, cj, max_marked=10):
    n = d.num_vertices
    dc, _ = d._dart_surface_components()
    darts = [h for h in range(4*n) if dc[h] == ci]
    # graph edges (undirected) with ids
    eidx = {}
    def E(a, b, tag):
        key = (min(a, b), max(a, b), tag)
        if key not in eidx: eidx[key] = len(eidx)
        return eidx[key]
    lam = {}
    def is_under(v, s):
        return (s & 1) == 0 if d.vtype[v] == "X0" else (s & 1) == 1
    for h in darts:
        E(h, d.opp[h], "e")
    for v in range(n):
        vt = d.vtype[v]
        if vt in ("X0", "X1"):
            for s in (0, 1):
                if dc[4*v+s] == ci:
                    e = E(4*v+s, 4*v+s+2, "c")
                    lam[e] = 1 if (is_under(v, s) and dc[4*v+((s+1) & 3)] == cj) else 0
        elif dc[4*v] == ci:
            for s in range(4):
                E(4*v+s, 4*v+((s+1) & 3), "m")
    ne = len(eidx)
    marked = [v for v in range(n) if d.vtype[v] in ("Y", "Z") and dc[4*v] == ci]
    def slice_cycles(S):
        seen = set(); cycs = []
        for h in darts:
            if h in seen: continue
            cyc = []; cur = h
            while cur not in seen:
                seen.add(cur)
                v, s = cur >> 2, cur & 3
                vt = d.vtype[v]
                if vt in ("X0", "X1"): ex = 4*v + (s ^ 2)
                else:
                    sgn = 1 if v in S else -1
                    p03 = (vt == "Y" and sgn > 0) or (vt == "Z" and sgn < 0)
                    ex = 4*v + ((3 - s) if p03 else (s ^ 1))
                seen.add(ex); cyc.append((cur, ex)); cur = d.opp[ex]
            cycs.append(cyc)
        return cycs
    def mask(cyc):
        m = 0
        for a, b in cyc:
            v = a >> 2
            tag = "c" if d.vtype[v] in ("X0", "X1") else "m"
            m ^= 1 << E(a, b, tag)
            m ^= 1 << E(b, d.opp[b], "e")
        return m
    assert len(eidx) == ne
    # boundaries: L- cycles, L+ cycles, vertex squares
    B = []
    for S in (set(), set(marked)):
        for cyc in slice_cycles(S): B.append(mask(cyc))
    for v in marked:
        m = 0
        for s in range(4): m ^= 1 << E(4*v+s, 4*v+((s+1) & 3), "m")
        B.append(m)
    basis = {}
    def red(x):
        while x:
            hb = x.bit_length() - 1
            if hb in basis: x ^= basis[hb]
            else: break
        return x
    for b in B:
        r = red(b)
        if r: basis[r.bit_length() - 1] = r
    # fully reduce helper (canonical form)
    def canon(x):
        for hb in sorted(basis, reverse=True):
            if (x >> hb) & 1: x ^= basis[hb]
        return x
    def lam_of(m):
        t = 0
        for e, w in lam.items():
            if (m >> e) & 1: t ^= w
        return t
    for b in B:
        assert lam_of(b) == 0, "lambda not a cocycle"
    qvals = {}   # class canon -> set of q values
    lamvals = {}
    subsets = []
    if len(marked) <= max_marked:
        for r in range(len(marked) + 1):
            for comb in combinations(marked, r): subsets.append(set(comb))
    else:
        import random
        rng = random.Random(1)
        for _ in range(600): subsets.append({v for v in marked if rng.random() < .5})
    self_vs = [v for v in range(n) if d.vtype[v] in ("X0", "X1") and dc[4*v] == ci and dc[4*v+1] == ci]
    for S in subsets:
        cycs = slice_cycles(S)
        which = {}; indart = set()
        for k, cyc in enumerate(cycs):
            for a, b in cyc:
                which[a] = k; which[b] = k; indart.add(a)
        lk2 = {}
        for v in self_vs:
            ka, kb = which[4*v], which[4*v+1]
            if ka == kb: continue
            vt = d.vtype[v]
            us = (0, 2) if vt == "X0" else (1, 3); os_ = (1, 3) if vt == "X0" else (0, 2)
            u_in = us[0] if 4*v+us[0] in indart else us[1]
            o_in = os_[0] if 4*v+os_[0] in indart else os_[1]
            sgn = 1 if ((o_in - u_in) & 3) == 3 else -1
            key = (min(ka, kb), max(ka, kb))
            lk2[key] = lk2.get(key, 0) + sgn
        for k, cyc in enumerate(cycs):
            cl = canon(mask(cyc))
            if cl == 0: continue
            q = 0
            for key, val in lk2.items():
                if k in key:
                    assert val % 2 == 0
                    q ^= (val // 2) & 1
            qvals.setdefault(cl, set()).add(q)
            lamvals.setdefault(cl, set()).add(lam_of(mask(cyc)))
    return qvals, lamvals

def dlk_true(d):
    topo = d.surface_topology()
    out = []
    c = topo.num_components
    for i in range(c):
        ct = topo.components[i]
        if not (ct.orientable and ct.genus_symbol == 1): continue
        for j in range(c):
            if j == i: continue
            qv, lv = analyse(d, i, j)
            consistent = all(len(s) == 1 for s in qv.values()) and all(len(s) == 1 for s in lv.values())
            ker = [cl for cl in qv if lv[cl] == {0}]
            nonker = [cl for cl in qv if lv[cl] == {1}]
            val = None
            if consistent and len(ker) == 1 and nonker:
                val = next(iter(qv[ker[0]]))
            elif consistent and not ker and len(nonker) == 2:
                # q(a+b) = q(a)+q(b)+1
                a, b = nonker
                val = next(iter(qv[a])) ^ next(iter(qv[b])) ^ 1
            cj = topo.components[j]
            out.append((ct.homeomorphism_type, f"{cj.homeomorphism_type}(e={cj.normal_euler})", val, consistent,
                        len(qv), sorted((sorted(qv[c_]), sorted(lv[c_])) for c_ in qv)))
    return sorted(out, key=lambda x: repr(x))
