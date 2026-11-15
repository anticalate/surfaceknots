"""Q2: empirical soundness of every move generator in moves.py.
Usage: t4 mode seed   (mode = onestep | walk)"""
import json
import random
import sys
import time
from collections import Counter

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from lib import *  # noqa
from surfaceknots import YOSHIKAWA_BY_LABEL, MarkedDiagram
from surfaceknots import moves as M
from surfaceknots.invariants.fundamental_group import wirtinger_presentation
from surfaceknots.unlink import certify_classical_unlink

mode = sys.argv[1]
seed = int(sys.argv[2])
rng = random.Random(seed)
DATA = "/tmp/claude-29349/-afs-csail-mit-edu-u-e-ejli-surfaceknots/cbe399a3-0fb9-5aba-8e99-ce46cd6219aa/scratchpad/sk/code/data/expanded_ch11.json"

base = [(lab, MarkedDiagram.from_epd(e.canonical_epd)) for lab, e in YOSHIKAWA_BY_LABEL.items() if lab != "0_1"]
j = json.load(open(DATA))
for e in j["entries"] + j["residual_collisions"]:
    base.append((e["label"], MarkedDiagram.from_epd(e["epd"])))

_cache = {}


def inv(d):
    k = (d.opp, d.vtype, d.free_loops)
    if k in _cache:
        return _cache[k]
    pres = wirtinger_presentation(d)
    r = (
        surf_inv(d),
        tuple(indep_hom_count(d, g) for g in ("S3", "Q8", "A4")),
        pres.abelianization_invariant_factors(),
    )
    if len(_cache) > 200000:
        _cache.clear()
    _cache[k] = r
    return r


def admissible(d):
    return certify_classical_unlink(d.resolve(1)).is_unlink and certify_classical_unlink(d.resolve(-1)).is_unlink


GENS = {
    "O3/O4": M._list_omega3_and_4_neighbors,
    "O5": M._list_omega5_neighbors,
    "O7": M._list_omega7_neighbors,
    "O8": M._list_omega8_neighbors,
    "tubeflip": M._list_marked_tube_flip_neighbors,
    "flype": M.list_4d_tait_flypes,
    "O1/6+": M.list_omega1_6_stabilizations,
    "O2+": M.list_omega2_pokes,
    "O1/6-": lambda d: [x for x in [M.find_omega1_or_6_reduction(d)] if x is not None],
    "O2-": lambda d: [x for x in [M.find_omega2_reduction(d)] if x is not None],
}
counts = Counter()
fails = []


def check(lab, tag, d0, d1, i0=None):
    counts[tag] += 1
    problems = []
    if not marked_ok(d1):
        problems.append("nonplanar_or_invalid")
    else:
        if i0 is None:
            i0 = inv(d0)
        i1 = inv(d1)
        if i0 != i1:
            problems.append(f"invariant_changed {i0} -> {i1}")
        if not admissible(d1):
            problems.append("not_certified_admissible")
    if problems:
        fails.append({"label": lab, "move": tag, "before": d0.to_epd(), "before_vt": list(d0.vtype), "before_opp": list(d0.opp),
                      "after": d1.to_epd(), "after_vt": list(d1.vtype), "after_opp": list(d1.opp), "problems": problems})
        print("FAIL", lab, tag, problems, flush=True)
        return False
    return True


t0 = time.time()
if mode == "onestep":
    for lab, d in base:
        assert marked_ok(d) and admissible(d), lab
        i0 = inv(d)
        for tag, fn in GENS.items():
            nb = fn(d)
            if tag in ("O1/6+", "O2+") and len(nb) > 60:
                nb = rng.sample(nb, 60)
            for x in nb:
                check(lab, tag, d, x, i0)
    print("onestep done", dict(counts), "fails", len(fails), f"{time.time()-t0:.0f}s")
else:
    nwalks = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    steps = int(sys.argv[4]) if len(sys.argv) > 4 else 60
    maxn = 15
    for lab, d in base:
        for w in range(nwalks):
            cur = d
            for st in range(steps):
                avail = []
                for tag, fn in GENS.items():
                    if tag in ("O1/6+", "O2+") and cur.num_vertices >= maxn - (1 if tag == "O1/6+" else 2) + 1:
                        continue
                    nb = fn(cur)
                    if nb:
                        avail.append((tag, nb))
                if not avail:
                    break
                # weight: prefer rarer level moves
                wts = [3.0 if t in ("O8", "O7", "O5", "tubeflip", "flype") else (1.5 if t in ("O3/O4",) else 1.0) for t, _ in avail]
                tag, nb = rng.choices(avail, weights=wts)[0]
                nxt = rng.choice(nb)
                ok = check(lab, tag, cur, nxt)
                if not ok:
                    break
                cur = nxt
        print(lab, "walks done", dict(counts), "fails", len(fails), f"{time.time()-t0:.0f}s", flush=True)
json.dump(fails, open(f"{__file__.rsplit('/', 1)[0]}/t4_fails_{mode}_{seed}.json", "w"))
print("FINAL", mode, seed, dict(counts), "fails", len(fails))
