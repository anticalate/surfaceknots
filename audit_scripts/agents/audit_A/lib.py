"""Independent helpers for the audit (planarity, Jones polynomial, colorings, PD parsing)."""
from __future__ import annotations

from surfaceknots.diagrams import ClassicalLinkDiagram, MarkedDiagram


def valid_involution(opp, n):
    if len(opp) != 4 * n:
        return False
    for d in range(4 * n):
        o = opp[d]
        if not (0 <= o < 4 * n) or o == d or opp[o] != d:
            return False
    return True


def vertex_components(opp, n):
    comp = [-1] * n
    c = 0
    for v in range(n):
        if comp[v] != -1:
            continue
        st = [v]
        comp[v] = c
        while st:
            x = st.pop()
            for s in range(4):
                y = opp[4 * x + s] >> 2
                if comp[y] == -1:
                    comp[y] = c
                    st.append(y)
        c += 1
    return comp, c


def is_planar(opp, n):
    """Genus 0 on every connected component of the 4-valent map (darts CCW)."""
    if n == 0:
        return True
    if not valid_involution(opp, n):
        return False
    comp, nc = vertex_components(opp, n)
    verts = [0] * nc
    for v in range(n):
        verts[comp[v]] += 1
    faces = [0] * nc
    seen = [False] * (4 * n)
    for d in range(4 * n):
        if seen[d]:
            continue
        faces[comp[d >> 2]] += 1
        cur = d
        while not seen[cur]:
            seen[cur] = True
            od = opp[cur]
            cur = (od & ~3) | ((od - 1) & 3)
    # V - E + F = 2 with E = 2V  ->  F = V + 2
    return all(faces[i] == verts[i] + 2 for i in range(nc))


def over_slots(vt):
    return (1, 3) if vt == "X0" else (0, 2)


def jones(d: ClassicalLinkDiagram):
    """Normalized Kauffman bracket f(A) = (-A^3)^{-w} <D>, with <O> = 1, as a sorted tuple of (exp, coeff)."""
    n = d.num_crossings
    opp = d.opp
    # orientation / writhe computed independently
    seen = [False] * (4 * n)
    direction = [0] * (4 * n)
    ncirc = d.free_loops
    for st in range(4 * n):
        if seen[st]:
            continue
        ncirc += 1
        cur = st
        while not seen[cur]:
            seen[cur] = True
            seen[cur ^ 2] = True
            direction[cur] = 1
            direction[cur ^ 2] = -1
            cur = opp[cur ^ 2]
    w = 0
    for v in range(n):
        o = over_slots(d.vtype[v])
        u = (o[0] ^ 1, o[1] ^ 1) if d.vtype[v] == "X1" else (0, 2)
        if d.vtype[v] == "X0":
            u = (0, 2)
        else:
            u = (1, 3)
        oin = o[0] if direction[4 * v + o[0]] == 1 else o[1]
        uin = u[0] if direction[4 * v + u[0]] == 1 else u[1]
        # over direction points to slot oin+2, under direction to uin+2; positive iff under dir is 90deg CCW of over dir
        w += 1 if ((uin + 2) - (oin + 2)) % 4 == 1 else -1
    bracket = {}
    for mask in range(1 << n):
        apow = 0
        pair = [0] * (4 * n)
        for v in range(n):
            o = over_slots(d.vtype[v])[0]
            b = (mask >> v) & 1
            if b == 0:
                # A smoothing: pairs (o, o-1) and (o+1, o+2)
                apow += 1
                prs = ((o, (o - 1) & 3), ((o + 1) & 3, (o + 2) & 3))
            else:
                apow -= 1
                prs = ((o, (o + 1) & 3), ((o + 2) & 3, (o + 3) & 3))
            for a, bb in prs:
                pair[4 * v + a] = 4 * v + bb
                pair[4 * v + bb] = 4 * v + a
        sn = [False] * (4 * n)
        loops = d.free_loops
        for st in range(4 * n):
            if sn[st]:
                continue
            loops += 1
            cur = st
            while not sn[cur]:
                sn[cur] = True
                p = pair[cur]
                sn[p] = True
                cur = opp[p]
        m = loops - 1
        # (-A^2 - A^-2)^m
        from math import comb

        for k in range(m + 1):
            e = apow + 2 * (m - 2 * k)
            bracket[e] = bracket.get(e, 0) + ((-1) ** m) * comb(m, k)
    res = {}
    for e, c in bracket.items():
        if c:
            res[e - 3 * w] = c * ((-1) ** (w % 2))
    return tuple(sorted(res.items())), ncirc


def unlink_jones(m):
    from math import comb

    res = {}
    mm = m - 1
    for k in range(mm + 1):
        e = 2 * (mm - 2 * k)
        res[e] = res.get(e, 0) + ((-1) ** mm) * comb(mm, k)
    return tuple(sorted(res.items()))


def colorings(d: ClassicalLinkDiagram, p: int) -> int:
    """Number of Fox p-colorings (p prime) -- p^(nullity)."""
    n = d.num_crossings
    opp = d.opp
    # arcs: union darts along edges and through over-strands
    parent = list(range(4 * n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for x in range(4 * n):
        a, b = find(x), find(opp[x])
        if a != b:
            parent[a] = b
    for v in range(n):
        o = over_slots(d.vtype[v])
        a, b = find(4 * v + o[0]), find(4 * v + o[1])
        if a != b:
            parent[a] = b
    roots = sorted({find(x) for x in range(4 * n)})
    idx = {r: i for i, r in enumerate(roots)}
    na = len(roots)
    rows = []
    for v in range(n):
        o = over_slots(d.vtype[v])
        u = (o[0] ^ 1, o[1] ^ 1)
        row = [0] * na
        row[idx[find(4 * v + o[0])]] += 2
        row[idx[find(4 * v + u[0])]] -= 1
        row[idx[find(4 * v + u[1])]] -= 1
        rows.append([x % p for x in row])
    # rank mod p
    rank = 0
    rows = [r[:] for r in rows]
    for col in range(na):
        piv = None
        for r in range(rank, len(rows)):
            if rows[r][col] % p:
                piv = r
                break
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        inv = pow(rows[rank][col], p - 2, p)
        rows[rank] = [(x * inv) % p for x in rows[rank]]
        for r in range(len(rows)):
            if r != rank and rows[r][col] % p:
                f = rows[r][col]
                rows[r] = [(x - f * y) % p for x, y in zip(rows[r], rows[rank])]
        rank += 1
    return p ** (na - rank + d.free_loops)


def classical_from_pd(pd):
    """pd: list of (a,b,c,d) with a = incoming under, CCW (KnotInfo convention). All X0."""
    n = len(pd)
    occ = {}
    for v, q in enumerate(pd):
        for s, lbl in enumerate(q):
            occ.setdefault(lbl, []).append(4 * v + s)
    opp = [0] * (4 * n)
    for lbl, ds in occ.items():
        assert len(ds) == 2, (lbl, ds)
        opp[ds[0]] = ds[1]
        opp[ds[1]] = ds[0]
    return ClassicalLinkDiagram(num_crossings=n, opp=tuple(opp), vtype=tuple(["X0"] * n))


def mirror_classical(d):
    return ClassicalLinkDiagram(
        num_crossings=d.num_crossings,
        opp=d.opp,
        vtype=tuple("X1" if t == "X0" else "X0" for t in d.vtype),
        free_loops=d.free_loops,
    )


def disjoint_union_classical(a, b):
    na = a.num_crossings
    opp = list(a.opp) + [x + 4 * na for x in b.opp]
    return ClassicalLinkDiagram(
        num_crossings=na + b.num_crossings, opp=tuple(opp), vtype=a.vtype + b.vtype, free_loops=a.free_loops + b.free_loops
    )


def connected_sum_classical(a, b, da=0, db=0):
    """Cut edge at dart da of a and dart db of b and splice (planar connected sum)."""
    na = a.num_crossings
    opp = list(a.opp) + [x + 4 * na for x in b.opp]
    a1, a2 = da, a.opp[da]
    b1, b2 = db + 4 * na, b.opp[db] + 4 * na
    opp[a1] = b2
    opp[b2] = a1
    opp[a2] = b1
    opp[b1] = a2
    return ClassicalLinkDiagram(num_crossings=na + b.num_crossings, opp=tuple(opp), vtype=a.vtype + b.vtype)


PD = {
    "trefoil": [(1, 5, 2, 4), (3, 1, 4, 6), (5, 3, 6, 2)],
    "figure8": [(4, 2, 5, 1), (8, 6, 1, 5), (6, 3, 7, 4), (2, 7, 3, 8)],
    "hopf": [(4, 1, 3, 2), (2, 3, 1, 4)],
    "5_1": [(2, 8, 3, 7), (4, 10, 5, 9), (6, 2, 7, 1), (8, 4, 9, 3), (10, 6, 1, 5)],
    "5_2": [(1, 5, 2, 4), (3, 9, 4, 8), (5, 1, 6, 10), (7, 3, 8, 2), (9, 7, 10, 6)],
    "6_1": [(1, 7, 2, 6), (3, 10, 4, 11), (5, 3, 6, 2), (7, 1, 8, 12), (9, 4, 10, 5), (11, 9, 12, 8)],
    "whitehead": [(6, 1, 7, 2), (10, 7, 5, 8), (4, 5, 1, 6), (2, 10, 3, 9), (8, 4, 9, 3)],
    "borromean": [(6, 1, 7, 2), (12, 8, 9, 7), (4, 12, 1, 11), (10, 5, 11, 6), (8, 4, 5, 3), (2, 9, 3, 10)],
}


# ---------------------------------------------------------------------------------------------
# Independent pi_1(S^4 \ F(D)) homomorphism counter (written from scratch; uses L_- + band relations)
# ---------------------------------------------------------------------------------------------
def _perm_group(gens):
    n = len(gens[0])
    ident = tuple(range(n))
    elems = {ident}
    frontier = [ident]
    while frontier:
        nf = []
        for e in frontier:
            for g in gens:
                h = tuple(g[e[i]] for i in range(n))
                if h not in elems:
                    elems.add(h)
                    nf.append(h)
        frontier = nf
    elems = sorted(elems)
    idx = {e: i for i, e in enumerate(elems)}
    m = len(elems)
    mul = [[idx[tuple(elems[a][elems[b][i]] for i in range(n))] for b in range(m)] for a in range(m)]
    inv = [0] * m
    e0 = idx[ident]
    for a in range(m):
        for b in range(m):
            if mul[a][b] == e0:
                inv[a] = b
    return mul, inv


_Q8_GENS = [(1, 2, 3, 0, 5, 6, 7, 4), (4, 7, 6, 5, 2, 1, 0, 3)]  # regular rep of Q8: i and j
GROUPS = {
    "S3": _perm_group([(1, 2, 0), (1, 0, 2)]),
    "A4": _perm_group([(1, 2, 0, 3), (0, 2, 3, 1)]),
    "Q8": _perm_group(_Q8_GENS),
    "S4": _perm_group([(1, 2, 3, 0), (1, 0, 2, 3)]),
}
assert len(GROUPS["Q8"][0]) == 8 and len(GROUPS["A4"][0]) == 12 and len(GROUPS["S4"][0]) == 24


def indep_presentation(d: MarkedDiagram):
    """Return (num_gens, relations) with relations as lists of signed ints (1-based) meaning product == 1."""
    n = d.num_vertices
    opp = d.opp
    vt = d.vtype

    def partner(x):
        v, s = x >> 2, x & 3
        t = vt[v]
        if t in ("X0", "X1"):
            return (v << 2) | (s ^ 2)
        if t == "Y":  # L-: (0,1),(2,3)
            return (v << 2) | (s ^ 1)
        return (v << 2) | (3 - s)  # Z L-: (0,3),(1,2)

    direction = [0] * (4 * n)
    for st in range(4 * n):
        if direction[st]:
            continue
        cur = st
        while not direction[cur]:
            direction[cur] = 1
            ex = partner(cur)
            direction[ex] = -1
            cur = opp[ex]
    parent = list(range(4 * n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        a, b = find(a), find(b)
        if a != b:
            parent[a] = b

    for x in range(4 * n):
        union(x, opp[x])
    for v in range(n):
        t = vt[v]
        if t == "X0":
            union(4 * v + 1, 4 * v + 3)
        elif t == "X1":
            union(4 * v, 4 * v + 2)
        else:
            for s in range(4):
                union(4 * v + s, partner(4 * v + s))
    roots = sorted({find(x) for x in range(4 * n)})
    gid = {r: i + 1 for i, r in enumerate(roots)}
    ng = len(roots) + d.free_loops
    rels = []
    for v in range(n):
        t = vt[v]
        if t in ("X0", "X1"):
            o = (1, 3) if t == "X0" else (0, 2)
            u = (0, 2) if t == "X0" else (1, 3)
            oin = o[0] if direction[4 * v + o[0]] == 1 else o[1]
            uin = u[0] if direction[4 * v + u[0]] == 1 else u[1]
            uout = uin ^ 2
            positive = (uin - oin) % 4 == 1
            go = gid[find(4 * v + oin)]
            ga = gid[find(4 * v + uin)]
            gb = gid[find(4 * v + uout)]
            if positive:  # b = o^-1 a o
                rels.append([-gb, -go, ga, go])
            else:  # b = o a o^-1
                rels.append([-gb, go, ga, -go])
        else:
            p = 0
            q = partner(4 * v + p) & 3
            # the other pair
            rest = [s for s in range(4) if s not in (p, q)]
            # strand A enters at p or q; B's slot opposite to A's entering slot
            a_in = p if direction[4 * v + p] == 1 else q
            b_opp = a_in ^ 2  # slot opposite to a_in belongs to strand B
            assert b_opp in rest
            ga = gid[find(4 * v + p)]
            gb = gid[find(4 * v + rest[0])]
            # the surface turns over like an inverted U across the band, so a meridian transported across the
            # band reverses its planar direction: antiparallel strands -> a = b ; parallel strands -> a = b^-1
            if direction[4 * v + b_opp] == 1:  # antiparallel
                rels.append([ga, -gb])
            else:
                rels.append([ga, gb])
    return ng, rels


def indep_hom_count(d: MarkedDiagram, gname: str) -> int:
    mul, inv = GROUPS[gname]
    m = len(mul)
    ng, rels = indep_presentation(d)
    if ng == 0:
        return 1
    # occurrences
    by_gen = [[] for _ in range(ng + 1)]
    for ri, r in enumerate(rels):
        for g in set(abs(x) for x in r):
            by_gen[g].append(ri)
    assign = [None] * (ng + 1)

    def eval_partial(r):
        """Return ('ok'|'bad'|('force', gen, val)|'open')."""
        unknown = [abs(x) for x in r if assign[abs(x)] is None]
        uset = set(unknown)
        if not uset:
            acc = None
            for x in r:
                e = assign[abs(x)] if x > 0 else inv[assign[abs(x)]]
                acc = e if acc is None else mul[acc][e]
            ident = mul[acc][inv[acc]]
            return "ok" if acc == ident else "bad"
        if len(uset) == 1 and len(unknown) == 1:
            # rotate relation so unknown is first: x * rest == 1 -> x = rest^-1
            k = next(i for i, x in enumerate(r) if assign[abs(x)] is None)
            rot = r[k:] + r[:k]
            acc = None
            for x in rot[1:]:
                e = assign[abs(x)] if x > 0 else inv[assign[abs(x)]]
                acc = e if acc is None else mul[acc][e]
            val = inv[acc]
            if rot[0] < 0:
                val = inv[val]
            return ("force", abs(rot[0]), val)
        return "open"

    def propagate(start_gens, trail):
        queue = list(start_gens)
        while queue:
            g = queue.pop()
            for ri in by_gen[g]:
                res = eval_partial(rels[ri])
                if res == "bad":
                    return False
                if isinstance(res, tuple):
                    _, gg, val = res
                    assign[gg] = val
                    trail.append(gg)
                    queue.append(gg)
        return True

    def rec():
        g = next((i for i in range(1, ng + 1) if assign[i] is None), None)
        if g is None:
            for r in rels:
                if eval_partial(r) != "ok":
                    return 0
            return 1
        total = 0
        for val in range(m):
            trail = [g]
            assign[g] = val
            if propagate([g], trail):
                total += rec()
            for t in trail:
                assign[t] = None
        return total

    return rec()


def surf_inv(d: MarkedDiagram):
    st = d.surface_topology()
    return tuple(sorted((c.orientable, c.euler_characteristic, c.normal_euler) for c in st.components))


def marked_ok(d: MarkedDiagram):
    return valid_involution(d.opp, d.num_vertices) and is_planar(d.opp, d.num_vertices)
