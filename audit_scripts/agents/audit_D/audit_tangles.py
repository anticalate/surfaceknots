"""Audit of tangle_3d_key dedup: is it the isotopy-rel-boundary classification?  Also arch-mask soundness."""
import sys, time, itertools, hashlib, pickle
from collections import defaultdict
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, __file__.rsplit('/',1)[0])
from common import *
from surfaceknots.triplane_census import enumerate_minimal_braided_tangles, tangle_3d_key
from surfaceknots.triplane import tangle_endpoint_pairing

# S3 tables
perms = list(itertools.permutations(range(3)))
pid = {p: i for i, p in enumerate(perms)}
MUL = np.zeros((6, 6), dtype=np.int8)
INV = np.zeros(6, dtype=np.int8)
for i, p in enumerate(perms):
    for j, q in enumerate(perms):
        MUL[i, j] = pid[tuple(p[q[k]] for k in range(3))]
for i in range(6):
    for j in range(6):
        if MUL[i, j] == pid[(0, 1, 2)]:
            INV[i] = j
IDN = pid[(0, 1, 2)]

def hom_hash(b, tw):
    N = 6 ** b
    idx = np.arange(N)
    img = np.stack([(idx // 6 ** r) % 6 for r in range(b)]).astype(np.int8)
    imginv = INV[img]
    code = np.zeros(N, dtype=np.int64)
    for w in tw:
        acc = np.full(N, IDN, dtype=np.int8)
        for lit in w:
            acc = MUL[acc, img[lit - 1] if lit > 0 else imginv[-lit - 1]]
        code = code * 6 + acc
    code.sort()
    return hashlib.md5(code.tobytes()).hexdigest()

def work(args):
    b, m, bw = args
    tw = top_words(b, m, bw)
    return hom_hash(b, tw)

if __name__ == "__main__":
    b = int(sys.argv[1]); maxc = int(sys.argv[2])
    t0 = time.time()
    by_c, key_to_loc = enumerate_minimal_braided_tangles(b, maxc)
    print(f"b={b} repo counts by c:", [len(x) for x in by_c], f"({time.time()-t0:.1f}s)", flush=True)
    allt = [(c, t) for c in range(maxc + 1) for t in by_c[c]]
    # sanity: phi_T kills own cap words; repo top words == mine
    from surfaceknots.triplane_census import _tangle_top_words, _tangle_cap_port_words
    import random
    random.seed(1)
    for c, t in random.sample(allt, min(2000, len(allt))):
        tw = top_words(b, t.bottom_matching, t.braid_word)
        assert tuple(tw) == tuple(_tangle_top_words(t))
        assert isotopic(cap_words(b, t.bottom_matching, t.braid_word), tw)
    print("self-consistency ok", flush=True)
    # arch mask soundness: tau[u]==u+1 and len(w_u)==1  =>  w_u w_{u+1} == 1 ?
    bad_arch = 0; n_arch = 0; n_pair_nontriv = 0
    for c, t in allt:
        tw = top_words(b, t.bottom_matching, t.braid_word)
        tau = tangle_endpoint_pairing(t)
        for u in range(2 * b - 1):
            if tau[u] == u + 1:
                triv = len(red(tw[u] + tw[u + 1])) == 0
                if len(tw[u]) == 1:
                    n_arch += 1
                    if not triv:
                        bad_arch += 1
                        if bad_arch <= 5: print("  BAD ARCH", t, u, tw[u], tw[u+1])
                elif triv:
                    n_pair_nontriv += 1
    print(f"arch-mask bits: {n_arch}, of which NOT a trivial cap (w_u w_u+1 != 1): {bad_arch}; trivial caps missed by len==1 test: {n_pair_nontriv}", flush=True)
    # true isotopy classes
    with Pool(10) as pool:
        hs = pool.map(work, [(b, t.bottom_matching, t.braid_word) for c, t in allt], chunksize=200)
    print(f"hashes done ({time.time()-t0:.1f}s)", flush=True)
    buckets = defaultdict(list)
    for i, (c, t) in enumerate(allt):
        buckets[(tangle_endpoint_pairing(t), hs[i])].append(i)
    cls_min_c = []   # per true class: min crossing level
    merged_examples = []
    n_tests = 0
    for key, idxs in buckets.items():
        reps = []  # (top words, min c, idx)
        for i in idxs:
            c, t = allt[i]
            caps = cap_words(b, t.bottom_matching, t.braid_word)
            found = False
            for r in reps:
                n_tests += 1
                if isotopic(caps, r[0]):
                    found = True
                    if len(merged_examples) < 8:
                        merged_examples.append((allt[r[2]], allt[i]))
                    break
            if not found:
                reps.append((top_words(b, t.bottom_matching, t.braid_word), c, i))
        cls_min_c.extend(r[1] for r in reps)
    true_counts = [sum(1 for x in cls_min_c if x == c) for c in range(maxc + 1)]
    print(f"TRUE isotopy-class counts by min crossing level: {true_counts}  total {len(cls_min_c)} vs repo total {len(allt)}; exact kernel tests: {n_tests} ({time.time()-t0:.1f}s)")
    for a, bb in merged_examples:
        print("  isotopic but distinct keys:", a[0], a[1].bottom_matching, a[1].braid_word, " ~ ", bb[0], bb[1].bottom_matching, bb[1].braid_word)
    pickle.dump({"true_counts": true_counts}, open(__file__.rsplit("/",1)[0] + f"/tangle_classes_b{b}.pkl", "wb"))
