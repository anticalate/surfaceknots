"""Q1: soundness of classical Reidemeister move implementations in unlink.py + verdicts on known links."""
import random
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from lib import *  # noqa
from surfaceknots.unlink import (
    certify_classical_unlink,
    find_r1_reduction,
    find_r2_reduction,
    find_unlink_movie,
    greedy_r1_r2,
    is_trivial_jones_polynomial,
    list_r1_plus_stabilizations,
    list_r2_plus_pokes,
    list_r3_neighbors,
)

rng = random.Random(1)

base = {k: classical_from_pd(v) for k, v in PD.items()}
for k, d in base.items():
    j, nc = jones(d)
    print(k, "planar", is_planar(d.opp, d.num_crossings), "comps", nc, "col3", colorings(d, 3), "col5", colorings(d, 5),
          "jones_trivial(mine)", j == unlink_jones(nc), "theirs", is_trivial_jones_polynomial(d),
          "lk2", d.pairwise_linking_numbers_doubled())

# extra diagrams
base["trefoil#trefoil"] = connected_sum_classical(base["trefoil"], base["trefoil"])
base["trefoil#mirror"] = connected_sum_classical(base["trefoil"], mirror_classical(base["trefoil"]))
base["fig8#fig8"] = connected_sum_classical(base["figure8"], base["figure8"])
base["trefoil+split+fig8"] = disjoint_union_classical(base["trefoil"], base["figure8"])
base["hopf#hopf"] = connected_sum_classical(base["hopf"], base["hopf"])
# simple unknot/unlink diagrams
kink = ClassicalLinkDiagram(1, (1, 0, 3, 2), ("X0",))
base["kink"] = kink
# two overlapping circles (R2 unlink): crossings u,w
# circle A through u(0,2), w(0,2); circle B through u(1,3), w(1,3). bigon u.1-w.? build via poke on two kinks instead
base["kink#kink"] = connected_sum_classical(kink, kink, 2, 2)


def check(name, d0, d1, tag):
    n = d1.num_crossings
    assert valid_involution(d1.opp, n), (name, tag, "invalid involution")
    assert is_planar(d1.opp, n), (name, tag, "non-planar")
    j0, c0 = jones(d0)
    j1, c1 = jones(d1)
    assert c0 == c1, (name, tag, "component count changed", c0, c1)
    if c0 == 1:
        assert j0 == j1, (name, tag, "Jones changed")
    else:
        def unit(j):
            if not j:
                return j
            e0 = j[0][0]
            sg = 1 if j[0][1] > 0 else -1
            return tuple((e - e0, c * sg) for e, c in j)
        assert unit(j0) == unit(j1), (name, tag, "Jones changed (up to unit)")
    for p in (3, 5):
        assert colorings(d0, p) == colorings(d1, p), (name, tag, "colorings changed", p)


stats = {}
for name, d in base.items():
    for trial in range(40):
        cur = d
        for step in range(25):
            opts = []
            r3 = list_r3_neighbors(cur)
            opts += [("R3", x) for x in r3]
            if cur.num_crossings <= 9:
                r1p = list_r1_plus_stabilizations(cur)
                r2p = list_r2_plus_pokes(cur)
                opts += [("R1+", x) for x in rng.sample(r1p, min(3, len(r1p)))]
                opts += [("R2+", x) for x in rng.sample(r2p, min(6, len(r2p)))]
            x = find_r1_reduction(cur)
            if x is not None and rng.random() < 0.3:
                opts.append(("R1-", x))
            x = find_r2_reduction(cur)
            if x is not None and rng.random() < 0.3:
                opts.append(("R2-", x))
            if not opts:
                break
            tag, nxt = rng.choice(opts)
            check(name, cur, nxt, tag)
            stats[tag] = stats.get(tag, 0) + 1
            cur = nxt
print("random-walk move applications all preserved planarity/components/Jones/colorings:", stats)

# exhaustive check of ALL neighbors for small diagrams (each generator)
cnt = 0
for name, d in base.items():
    for tag, fn in (("R3", list_r3_neighbors), ("R1+", list_r1_plus_stabilizations), ("R2+", list_r2_plus_pokes)):
        for nxt in fn(d):
            check(name, d, nxt, tag)
            cnt += 1
            if tag == "R2+":
                for n2 in list_r3_neighbors(nxt):
                    check(name, nxt, n2, "R3")
                    cnt += 1
                g, _ = greedy_r1_r2(nxt)
                check(name, nxt, g, "greedy")
print("exhaustive one-step neighbor checks passed:", cnt)

print("\nVERDICTS")
for name, d in base.items():
    for lab, dd in ((name, d), (name + "(mirror)", mirror_classical(d))):
        c = certify_classical_unlink(dd)
        mv = find_unlink_movie(dd)
        print(f"{lab:28s} crossings={dd.num_crossings} comps={c.num_components} is_unlink={c.is_unlink} reason={c.reason!r} nmoves={len(c.moves)} movie={'None' if mv is None else len(mv)}")
