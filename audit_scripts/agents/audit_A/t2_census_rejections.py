"""Q1: do the Python certifier's 'not certified' verdicts silently drop Jones-trivial (very likely admissible)
diagrams in the census pipeline?  Monkeypatch fast_enum.certify_classical_unlink to log every rejection."""
import json
import sys
import time

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from lib import *  # noqa
import surfaceknots.fast_enum as fe
import surfaceknots.unlink as ul

ch = int(sys.argv[1])
mode = sys.argv[2] if len(sys.argv) > 2 else "census"
log = []
orig = ul.certify_classical_unlink
calls = [0]


def patched(d, max_r3_states=512):
    calls[0] += 1
    c = orig(d, max_r3_states=max_r3_states)
    if not c.is_unlink:
        j, nc = jones(d)
        log.append({"reason": c.reason, "c": d.num_crossings, "comps": nc, "jones_trivial": j == unlink_jones(nc),
                    "opp": list(d.opp), "vtype": list(d.vtype), "free": d.free_loops})
    return c


fe.certify_classical_unlink = patched
t0 = time.time()
if mode == "census":
    out, stats = fe.enumerate_prime_ch_census(ch, num_workers=4)
    print("ch", ch, "kept", len(out), stats)
else:
    out = fe.enumerate_admissible_ch_diagrams(ch)
    print("ch", ch, "admissible certified", len(out))
print("time", time.time() - t0, "certify calls", calls[0], "rejections", len(log))
from collections import Counter

print(Counter((r["reason"], r["c"], r["jones_trivial"]) for r in log))
json.dump(log, open(f"{__file__.rsplit('/', 1)[0]}/rej_{mode}_{ch}.json", "w"))
