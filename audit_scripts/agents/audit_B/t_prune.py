"""Every candidate removed by prune_level 3 (vs 0) must have an explicit witness:
has_imm(D), or a C level-neighbour (validated separately) with has_imm or a 0/2-edge cut.
Also: has_1step implies such a witness (no stronger)."""
import sys, time, ctypes
from multiprocessing import Pool
import aud
from surfaceknots.fast_enum import _get_lib, _CandidateStruct, decode_candidate
from surfaceknots.quadrangulations import generate_prime_shadows, PlanarShadow

def cands(n, opp, prune, mmin, mmax):
    lib = _get_lib(); cap = 1 << 21
    buf = (_CandidateStruct * cap)()
    oa = (ctypes.c_int * (4*n))(*opp)
    cnt = lib.enumerate_shadow_candidates(n, oa, mmin, mmax, prune, buf, cap)
    assert cnt <= cap
    return set((buf[i].m_mask, buf[i].yz_mask, buf[i].sign_mask) for i in range(cnt))

def work(a):
    n, opp = a
    sh = PlanarShadow(num_vertices=n, opp=tuple(opp))
    s0 = cands(n, opp, 0, 0, n); s1 = cands(n, opp, 1, 0, n); s2 = cands(n, opp, 2, 0, n); s3 = cands(n, opp, 3, 0, n)
    assert s3 <= s2 <= s1 <= s0
    nowit = []
    for t in s0 - s3:
        d = decode_candidate(sh, *t)
        if aud.has_imm(d): continue
        ok = False
        for x in aud.level_neighbors(d):
            if aud.has_imm(x) or aud.has_cut(x): ok = True; break
        if not ok: nowit.append(t)
    nowit2 = []
    for t in s0 - s2:
        d = decode_candidate(sh, *t)
        if aud.has_imm(d): continue
        if not any(aud.has_imm(x) for x in aud.level_neighbors(d)): nowit2.append(t)
    return len(s0), len(s1), len(s2), len(s3), len(nowit), len(nowit2), nowit[:2]

if __name__ == "__main__":
    n = int(sys.argv[1]); procs = int(sys.argv[2])
    sh = generate_prime_shadows(n)
    t = time.time()
    with Pool(procs) as p: res = p.map(work, [(n, s.opp) for s in sh], chunksize=1)
    print(f"n={n} prune0={sum(r[0] for r in res)} prune1={sum(r[1] for r in res)} prune2={sum(r[2] for r in res)} prune3={sum(r[3] for r in res)} removed_without_witness(level3)={sum(r[4] for r in res)} (level2)={sum(r[5] for r in res)} time={time.time()-t:.1f}s")
