"""Brute force: every Python-certified admissible assignment must be emitted by C enumerate_shadow_candidates(prune_level=0)."""
import sys, time, ctypes, itertools
from multiprocessing import Pool
from surfaceknots.fast_enum import _get_lib, _CandidateStruct, decode_candidate
from surfaceknots.quadrangulations import generate_prime_shadows
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.unlink import certify_classical_unlink

def work(args):
    n, opp = args
    lib = _get_lib()
    cap = 1 << 20
    buf = (_CandidateStruct * cap)()
    oa = (ctypes.c_int * (4*n))(*opp)
    cnt = lib.enumerate_shadow_candidates(n, oa, 0, n, 0, buf, cap)
    assert cnt <= cap
    cset = set((buf[i].m_mask, buf[i].yz_mask, buf[i].sign_mask) for i in range(cnt))
    # brute force normalized assignments
    missing = []; n_adm = 0; n_c_not_cert = 0; reasons = {}
    seen_c = set()
    for m_mask in range(1 << n):
        mv = [v for v in range(n) if (m_mask >> v) & 1]
        cv = [v for v in range(n) if not (m_mask >> v) & 1]
        m = len(mv); c = len(cv)
        for yz in range(0, 1 << m, 2):
            for sm in range(0, 1 << c, 2):
                vt = [None]*n
                for i, v in enumerate(mv): vt[v] = "Z" if (yz >> i) & 1 else "Y"
                for i, v in enumerate(cv): vt[v] = "X1" if (sm >> i) & 1 else "X0"
                d = MarkedDiagram(num_vertices=n, opp=tuple(opp), vtype=tuple(vt), free_loops=0)
                cp = certify_classical_unlink(d.resolve(1))
                ok = cp.is_unlink
                cm = None
                if ok:
                    cm = certify_classical_unlink(d.resolve(-1)); ok = cm.is_unlink
                inC = (m_mask, yz, sm) in cset
                if ok:
                    n_adm += 1
                    if not inC: missing.append((m_mask, yz, sm))
                elif inC:
                    n_c_not_cert += 1
                    r = (cp.reason if not cp.is_unlink else cm.reason)
                    reasons[r] = reasons.get(r, 0) + 1
    return (cnt, n_adm, missing, n_c_not_cert, reasons)

if __name__ == "__main__":
    n = int(sys.argv[1]); procs = int(sys.argv[2])
    sh = generate_prime_shadows(n)
    t = time.time()
    with Pool(procs) as p:
        res = p.map(work, [(n, s.opp) for s in sh], chunksize=1)
    tot_c = sum(r[0] for r in res); tot_adm = sum(r[1] for r in res)
    miss = sum(len(r[2]) for r in res); cnc = sum(r[3] for r in res)
    reasons = {}
    for r in res:
        for k, v in r[4].items(): reasons[k] = reasons.get(k, 0) + v
    print(f"n={n} shadows={len(sh)} C_bracket_survivors={tot_c} py_certified_admissible={tot_adm} admissible_missing_from_C={miss} C_survivors_not_certified={cnc} reasons={reasons} time={time.time()-t:.1f}s")
    for i, r in enumerate(res):
        if r[2]: print(" shadow", i, "missing", r[2][:5])
