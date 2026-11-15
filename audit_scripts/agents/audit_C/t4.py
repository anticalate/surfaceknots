import sys, random, json
from multiprocessing import Pool
from mydlk import dlk_true
from variants import *
from surfaceknots import YOSHIKAWA_BY_LABEL, MarkedDiagram
from surfaceknots.invariants.linking import double_linking_mod2
LBL = ["6_1^0,1","8_1^1,1","10_1^1,1","9_1^1,-2","10_1^0,0,1"]
def work(a):
    lbl, seed, lvl = a
    rng = random.Random(seed)
    d0 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL[lbl].canonical_epd)
    d, path = random_variant(d0, rng, nsteps=rng.randint(1,5), use_level=lvl, max_v=16)
    try:
        r = double_linking_mod2(d)
    except Exception as ex:
        r = "EXC " + repr(ex)[:80]
    m = dlk_true(d)
    return lbl, seed, lvl, path, d.to_epd(), r, [(x[0],x[1],x[2]) for x in m], all(x[3] for x in m)
if __name__ == "__main__":
    N = int(sys.argv[1])
    base = {}
    for lbl in LBL:
        d0 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL[lbl].canonical_epd)
        base[lbl] = (double_linking_mod2(d0), [(x[0],x[1],x[2]) for x in dlk_true(d0)])
    jobs = [(lbl, i, lvl) for lbl in LBL for i in range(N) for lvl in (False, True)]
    with Pool(7) as p: res = p.map(work, jobs, chunksize=4)
    from collections import Counter
    c = Counter(); ex = {}
    for lbl, seed, lvl, path, epd, r, m, cons in res:
        k1 = (lbl, lvl, "repo_same" if r == base[lbl][0] else "repo_DIFF")
        k2 = (lbl, lvl, "mine_same" if m == base[lbl][1] else "mine_DIFF")
        k3 = (lbl, lvl, "mine_consistent" if cons else "mine_INCONSISTENT")
        for k in (k1,k2,k3):
            c[k]+=1
            if "DIFF" in k[2] or "INCONS" in k[2]:
                if k not in ex or len(path) < len(ex[k][0]): ex[k] = (path, epd, r, m)
    for k in sorted(c): print(k, c[k])
    for k,v in ex.items(): print("EX", k, v)
