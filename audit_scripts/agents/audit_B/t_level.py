"""Differential test of C level moves: planarity, admissibility, invariants, and comparison with Python moves."""
import sys, time, ctypes, random, itertools
from multiprocessing import Pool
import aud
from surfaceknots.fast_enum import _get_lib, _CandidateStruct, decode_candidate
from surfaceknots.quadrangulations import generate_prime_shadows
from surfaceknots.moves import list_level_yoshikawa_neighbors, canonical_marked_key

def cands(n, sh, prune, mmin, mmax):
    lib = _get_lib()
    cap = 1 << 21
    buf = (_CandidateStruct * cap)()
    oa = (ctypes.c_int * (4*n))(*sh.opp)
    cnt = lib.enumerate_shadow_candidates(n, oa, mmin, mmax, prune, buf, cap)
    assert cnt <= cap
    return [(buf[i].m_mask, buf[i].yz_mask, buf[i].sign_mask) for i in range(cnt)]

def py_cut(d):
    n = d.num_vertices; opp = d.opp
    if n <= 1: return 0
    def comp(banned):
        vis = {0}; st = [0]
        while st:
            u = st.pop()
            for s in range(4):
                dd = 4*u+s
                if dd in banned: continue
                w = opp[dd]//4
                if w not in vis: vis.add(w); st.append(w)
        return len(vis)
    if comp(()) < n: return 1
    if n < 4: return 0
    edges = [(a, opp[a]) for a in range(4*n) if a < opp[a]]
    for (a,b),(c,e) in itertools.combinations(edges, 2):
        k = comp({a,b,c,e})
        if 2 <= k <= n-2: return 1
    return 0

def work(args):
    n, shopp, triples = args
    from surfaceknots.quadrangulations import PlanarShadow
    sh = PlanarShadow(num_vertices=n, opp=tuple(shopp))
    st = dict(diag=0, nbrs=0, nonplanar=0, nonadm=0, inv_mismatch=0, c_only=0, py_only=0, cut_mismatch=0, flips=0)
    bad = []
    for t in triples:
        d = decode_candidate(sh, *t)
        st["diag"] += 1
        i0 = aud.inv(d)
        cn = aud.level_neighbors(d)
        pn = list_level_yoshikawa_neighbors(d)
        ck = {}
        for x in cn:
            ck.setdefault(canonical_marked_key(x), x)
        pk = set(canonical_marked_key(x) for x in pn)
        st["c_only"] += len(set(ck) - pk); st["py_only"] += len(pk - set(ck))
        if set(ck) - pk and len(bad) < 3: bad.append(("c_only", d.to_epd()))
        for k, x in ck.items():
            st["nbrs"] += 1
            if not aud.is_planar(x):
                st["nonplanar"] += 1; bad.append(("nonplanar", d.to_epd(), x.to_epd())); continue
            if aud.has_cut(x) != py_cut(x):
                st["cut_mismatch"] += 1; bad.append(("cut", x.to_epd()))
            if not aud.admissible(x):
                st["nonadm"] += 1; bad.append(("nonadm", d.to_epd(), x.to_epd())); continue
            if aud.inv(x) != i0:
                st["inv_mismatch"] += 1; bad.append(("inv", d.to_epd(), x.to_epd()))
    return st, bad[:5]

if __name__ == "__main__":
    n = int(sys.argv[1]); per_shadow = int(sys.argv[2]); procs = int(sys.argv[3])
    random.seed(12345)
    tasks = []
    for sh in generate_prime_shadows(n):
        c = cands(n, sh, 0, 1, n)
        if len(c) > per_shadow: c = random.sample(c, per_shadow)
        if c: tasks.append((n, sh.opp, c))
    t = time.time()
    with Pool(procs) as p:
        res = p.map(work, tasks, chunksize=1)
    tot = {}
    for st, bad in res:
        for k, v in st.items(): tot[k] = tot.get(k, 0) + v
    print(f"n={n}", tot, "time %.1fs" % (time.time()-t))
    shown = 0
    for st, bad in res:
        for b in bad:
            if shown < 8: print("  BAD", b); shown += 1
