"""Is the invariant key unchanged under the two mirror operations the census quotients by?
  M1: switch every crossing (X0 <-> X1)        -- reflection of R^3
  M2: switch every marker (Y <-> Z)            -- reversal of the time direction
The enumeration fixes the type of the first crossing and of the first marker, i.e. it keeps one diagram
from each orbit of {1, M1, M2, M1M2}.  A key that changes under M1 or M2 could make the mirror image of a
table entry look like a new class.  Run from code/ with PYTHONPATH=src."""
import json, collections
from surfaceknots import MarkedDiagram, YOSHIKAWA_TABLE
from surfaceknots.yoshikawa_data import compute_full_invariant_record, invariant_filtration_layers
c11 = json.load(open("data/expanded_ch11.json"))
pool = [(e.label, e.canonical_epd) for e in YOSHIKAWA_TABLE if e.canonical_epd] + [(r["label"], r["epd"]) for r in c11["entries"]]
def layers(d):
    rec = compute_full_invariant_record(d, "x", include_sublinks=False, compute_so4_symmetry_bfs=False)
    L = invariant_filtration_layers(rec)
    return {k: repr(v) for k, v in (L.items() if isinstance(L, dict) else L)}
def op(d, cross, mark):
    m = {"X0": "X1" if cross else "X0", "X1": "X0" if cross else "X1", "Y": "Z" if mark else "Y", "Z": "Y" if mark else "Z"}
    return MarkedDiagram(d.num_vertices, d.opp, tuple(m[t] for t in d.vtype), d.free_loops)
bad = collections.defaultdict(list); n = 0
for lbl, epd in pool:
    d = MarkedDiagram.from_epd(epd)
    if d.num_vertices == 0: continue
    base = layers(d)
    for name, (c, m) in {"M1": (1, 0), "M2": (0, 1), "M1M2": (1, 1)}.items():
        try: L = layers(op(d, c, m))
        except Exception as ex:
            bad["EXCEPTION"].append((lbl, name, repr(ex)[:80])); continue
        n += 1
        for k in base:
            if L[k] != base[k]: bad[k].append((lbl, name))
print("comparisons", n)
for k, v in bad.items(): print("LAYER", k, "changed on", len(v), ":", v[:12])
if not bad: print("no layer changed under M1, M2 or M1M2")
