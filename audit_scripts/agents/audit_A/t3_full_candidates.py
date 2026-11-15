"""Q1: run the Python certifier on EVERY Kauffman/lk-surviving candidate emitted by the C kernel at a given ch
(prime shadows), log any candidate whose L+ or L- is rejected.  Usage: t3 ch nshards shard"""
import ctypes
import json
import sys
import time
from collections import Counter

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import surfaceknots.fast_enum as fe
from surfaceknots.quadrangulations import generate_prime_shadows
from surfaceknots.unlink import certify_classical_unlink

ch, nsh, sh_i = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
m_min = 2
m_max = max(2, ch - 4)
lib = fe._get_lib()
shadows = generate_prime_shadows(ch)[sh_i::nsh]
buf_cap = 1 << 20
buf = (fe._CandidateStruct * buf_cap)()
t0 = time.time()
total = 0
rej = []
reasons = Counter()
for si, sh in enumerate(shadows):
    opp_arr = (ctypes.c_int * (4 * ch))(*sh.opp)
    cnt = lib.enumerate_shadow_candidates(ch, opp_arr, m_min, m_max, 2, buf, buf_cap)
    if cnt > buf_cap:
        buf_cap = cnt + 1024
        buf = (fe._CandidateStruct * buf_cap)()
        cnt = lib.enumerate_shadow_candidates(ch, opp_arr, m_min, m_max, 2, buf, buf_cap)
    for idx in range(cnt):
        c = buf[idx]
        d = fe.decode_candidate(sh, int(c.m_mask), int(c.yz_mask), int(c.sign_mask))
        total += 1
        for sg in (1, -1):
            cert = certify_classical_unlink(d.resolve(sg))
            if not cert.is_unlink:
                reasons[cert.reason] += 1
                rej.append({"epd": d.to_epd(), "sign": sg, "reason": cert.reason})
                break
    if si % 50 == 0:
        print(f"shard {sh_i}: shadow {si}/{len(shadows)} total {total} rejected {len(rej)} t={time.time()-t0:.0f}s", flush=True)
print(f"DONE shard {sh_i} ch {ch}: candidates {total} rejected {len(rej)} reasons {dict(reasons)} time {time.time()-t0:.0f}s", flush=True)
json.dump(rej, open(f"{__file__.rsplit('/', 1)[0]}/t3_rej_ch{ch}_{sh_i}.json", "w"))
