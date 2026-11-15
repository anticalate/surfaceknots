"""Exercise the >12-marker branch of compute_torus_dlk_pair (Carbon's tree).
For each table / ch=11 link with a torus component, add Omega_6 marked kinks (and a few Omega_1 kinks and
Omega_2 pokes) at random until the diagram has many markers, then compare double_linking_mod2 with the
value on the original diagram.  Omega moves are isotopies, so the value must not change and no pair may
disappear.  Run from code/ with PYTHONPATH=src."""
import json, random, sys, collections, inspect
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots import moves as mv
from surfaceknots.invariants import linking as lk
from surfaceknots.invariants.linking import double_linking_mod2
import itertools
_orig_comb = itertools.combinations
BR = collections.Counter() if False else None
from surfaceknots.yoshikawa_data import YOSHIKAWA_TABLE
rng = random.Random(int(sys.argv[1])); target = int(sys.argv[2]); reps = int(sys.argv[3])
c11 = json.load(open("data/expanded_ch11.json"))
pool = [(e.label, MarkedDiagram.from_epd(e.canonical_epd)) for e in YOSHIKAWA_TABLE if e.canonical_epd]
pool += [(r["label"], MarkedDiagram.from_epd(r["epd"])) for r in c11["entries"]]
nm = lambda d: sum(1 for t in d.vtype if t in ("Y", "Z"))
res = collections.Counter(); bad = []
src = inspect.getsource(lk.compute_torus_dlk_pair)
assert "if len(marked_vs) <= 12:" in src
src = src.replace("if len(marked_vs) <= 12:", "_BR[len(marked_vs) > 12] += 1\n    if len(marked_vs) <= 12:", 1)
_BR = collections.Counter(); lk.__dict__["_BR"] = _BR
exec(compile(src, "<patched>", "exec"), lk.__dict__)
for lbl, d0 in pool:
    base = double_linking_mod2(d0)
    topo0 = d0.surface_topology()
    if topo0.num_components < 2 or not any(c.orientable and c.genus_symbol == 1 for c in topo0.components): continue
    for rep in range(reps):
        d = d0
        while nm(d) < target:
            L = mv.list_omega1_6_stabilizations(d)
            Lm = [x for x in L if nm(x) > nm(d)]
            d = rng.choice(Lm if (Lm and rng.random() < 0.85) else L)
            if rng.random() < 0.1:
                P = mv.list_omega2_pokes(d)
                if P: d = rng.choice(P)
        topo = d.surface_topology()
        val = double_linking_mod2(d)
        res[(lbl, val == base)] += 1
        if val != base: bad.append((lbl, nm(d), d.num_vertices, base, val))
    print(lbl, "base", base, {k[1]: v for k, v in res.items() if k[0] == lbl}, flush=True)
print("times the branch choice was reached: <=12 markers", _BR[False], "; >12 markers", _BR[True])
print("total", sum(res.values()), "unchanged", sum(v for k, v in res.items() if k[1]))
for b in bad[:10]: print("CHANGED", b)
