"""Quantify how pairwise unlink decisions are made: C R-move reduction (1), C rigorous rejection (0), or Python fallback (2)."""
import sys, time, ctypes, pickle
from collections import Counter
from surfaceknots.triplane_census import enumerate_minimal_braided_tangles
from surfaceknots.fast_enum import _get_lib, pack_tangles_c_buffer
from surfaceknots.triplane import union_two_tangles_to_classical_link
from surfaceknots.unlink import certify_classical_unlink, greedy_r1_r2, is_trivial_jones_polynomial
lib = _get_lib()
b = int(sys.argv[1])
by_c, _ = enumerate_minimal_braided_tangles(b, 5)
bufs = [pack_tangles_c_buffer(x) for x in by_c]
levels = [(0,c) for c in range(6)] + [(1,c) for c in range(1,6)] + [(2,2),(2,3),(2,4),(3,3)]
tot = Counter(); reasons = Counter(); hard = []
for (c1, c2) in levels:
    n1, n2 = len(by_c[c1]), len(by_c[c2])
    out = (ctypes.c_uint8 * (n1 * n2))()
    lib.batch_check_tangle_pairs_c(n1, bufs[c1], n2, bufs[c2], out)
    raw = bytes(out)
    cnt = Counter(raw)
    lv = Counter()
    idx = raw.find(2)
    while idx != -1:
        i, j = divmod(idx, n2)
        l = union_two_tangles_to_classical_link(by_c[c1][i], by_c[c2][j])
        cert = certify_classical_unlink(l)
        if cert.is_unlink:
            key = "py_certified(r3=%d,%s)" % (min(cert.r3_count, 9), "stab" if any(m in ("R1+", "R2+") for m in cert.moves) else "mono")
        else:
            redd, _ = greedy_r1_r2(l)
            key = "py_rejected:%s(jones_trivial=%s,c_red=%d)" % (cert.reason, is_trivial_jones_polynomial(redd), redd.num_crossings)
            if len(hard) < 2000:
                hard.append((c1, i, c2, j, cert.reason))
        lv[key] += 1
        idx = raw.find(2, idx + 1)
    print(f"b={b} level ({c1},{c2}) pairs={n1*n2}: C0(rigorous reject)={cnt.get(0,0)} C1(C R-move unlink)={cnt.get(1,0)} C2(to python)={cnt.get(2,0)} :: {dict(lv)}", flush=True)
    for k, v in cnt.items(): tot[k] += v
    reasons.update(lv)
print("TOTAL", dict(tot), dict(reasons))
pickle.dump(hard, open(__file__.rsplit("/",1)[0] + f"/unlink_hard_b{b}.pkl", "wb"))
