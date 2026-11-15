"""Cross-check C pair classification against the independent Python certifier + full Jones polynomial on random samples."""
import random, ctypes
from collections import Counter
from surfaceknots.triplane_census import enumerate_minimal_braided_tangles
from surfaceknots.fast_enum import _get_lib, pack_tangles_c_buffer
from surfaceknots.triplane import union_two_tangles_to_classical_link
from surfaceknots.unlink import certify_classical_unlink, is_trivial_jones_polynomial
random.seed(3)
lib = _get_lib()
by_c, _ = enumerate_minimal_braided_tangles(4, 5)
res = Counter()
for (c1, c2, n) in [(1, 2, 1500), (2, 3, 2500), (1, 4, 2500), (3, 3, 3000), (2, 4, 3000), (1, 5, 3000), (0, 5, 1500)]:
    for _ in range(n):
        ta = random.choice(by_c[c1]); tb = random.choice(by_c[c2])
        out = (ctypes.c_uint8 * 1)()
        lib.batch_check_tangle_pairs_c(1, pack_tangles_c_buffer([ta]), 1, pack_tangles_c_buffer([tb]), out)
        st = out[0]
        l = union_two_tangles_to_classical_link(ta, tb)
        cert = certify_classical_unlink(l)
        lk = any(v != 0 for v in l.pairwise_linking_numbers_doubled().values())
        jt = is_trivial_jones_polynomial(l) if l.num_crossings > 0 else True
        res[(st, cert.is_unlink, "lk!=0" if lk else ("jones_trivial" if jt else "jones_nontrivial"))] += 1
for k, v in sorted(res.items()): print("C status %d | python certified unlink=%s | %s : %d" % (k[0], k[1], k[2], v))
