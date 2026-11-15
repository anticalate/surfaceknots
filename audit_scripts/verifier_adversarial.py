"""Adversarial test of move_search.verify_move_path_independent (Carbon's tree).

1. The known false path Tube(3_1) -> 8_1 # T^2 must be rejected.
2. Near misses: for each table / ch=11 diagram d, and for d after one random move, build d' by
   flipping one crossing (X0<->X1) or one marker (Y<->Z).  If d' is admissible, ask the single-step
   classifier and the full verifier about d -> d'.  Any acceptance is compared against the full
   invariant key of the two diagrams; an acceptance with different keys is a false certificate.
3. True moves: every neighbour produced by the repo's own generators must be accepted.
"""
import json, random, sys, itertools, collections
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots import move_search as ms
from surfaceknots.unlink import certify_classical_unlink
from surfaceknots.yoshikawa_data import YOSHIKAWA_TABLE, compute_full_invariant_record
from surfaceknots.yoshikawa_data import invariant_filtration_layers

random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
rest = "X[6,7,8,9], X[9,10,11,2], X[12,13,14,3], X[15,16,5,14], Z[16,17,18,7], X[18,19,10,8], Y[11,19,20,12], X[13,20,17,15]"
a = f"X[2,21,21,22], Z[1,22,3,4], Y[4,5,6,1], {rest}"
b = f"Y[2,21,22,23], Z[23,22,21,3], X[5,24,24,6], {rest}"
v = ms.verify_move_path_independent([a, b], allow_o4=True)
print("known false path:", v.is_valid, getattr(v, "step_kinds", None), str(getattr(v, "reason", getattr(v, "message", "")))[:160])
print("  template check alone:", ms._indep_is_single_local_move(MarkedDiagram.from_epd(a), MarkedDiagram.from_epd(b), True))

def key(d):
    r = compute_full_invariant_record(d, "x", include_sublinks=False, compute_so4_symmetry_bfs=False)
    L = dict(invariant_filtration_layers(r)) if not isinstance(invariant_filtration_layers(r), dict) else invariant_filtration_layers(r)
    return {k: repr(x) for k, x in L.items()}

def admissible(d):
    try:
        return certify_classical_unlink(d.resolve(1)).is_unlink and certify_classical_unlink(d.resolve(-1)).is_unlink
    except Exception:
        return False

def flips(d):
    for vtx in range(d.num_vertices):
        vt = list(d.vtype)
        vt[vtx] = {"X0": "X1", "X1": "X0", "Y": "Z", "Z": "Y"}[vt[vtx]]
        try:
            yield vtx, MarkedDiagram(d.num_vertices, d.opp, tuple(vt), d.free_loops)
        except Exception:
            try:
                yield vtx, MarkedDiagram(num_vertices=d.num_vertices, opp=d.opp, vtype=tuple(vt), free_loops=d.free_loops)
            except Exception as e:
                print("cannot construct", e); return

c11 = json.load(open("data/expanded_ch11.json"))
pool = [(e.label, MarkedDiagram.from_epd(e.canonical_epd)) for e in YOSHIKAWA_TABLE if e.canonical_epd]
pool += [(r["label"], MarkedDiagram.from_epd(r["epd"])) for r in c11["entries"]]
pool = [(l, d) for l, d in pool if d.num_vertices > 0]
print("pool", len(pool))

stats = collections.Counter(); bad = []
true_moves = collections.Counter(); true_rejected = []
for lbl, d0 in pool:
    n = d0.num_vertices
    nbs = ms.list_labeled_single_moves(d0, n + 2)
    bytag = collections.defaultdict(list)
    for tag, nb in nbs: bytag[tag].append(nb)
    sample = []
    for tag, L in bytag.items():
        random.shuffle(L)
        sample += [(tag, x) for x in L[:2]]
    k0 = None
    for tag, nb in sample:
        k = ms._indep_is_single_local_move(d0, nb, True)
        true_moves[(tag, k is not None)] += 1
        if k is None: true_rejected.append((lbl, tag))
        for vtx, d2 in flips(nb):
            stats["flips"] += 1
            if not admissible(d2): continue
            stats["admissible"] += 1
            k = ms._indep_is_single_local_move(d0, d2, True)
            if k is None: continue
            stats["template_accept:" + k] += 1
            if k0 is None: k0 = key(d0)
            k2 = key(d2)
            diff = [m for m in k0 if k0[m] != k2.get(m)]
            if diff:
                vv = ms.verify_move_path_independent([d0.to_epd(), d2.to_epd()], allow_o4=True)
                stats["template_accept_key_differs"] += 1
                if vv.is_valid: stats["FULL_VERIFIER_FALSE_CERT"] += 1
                bad.append((lbl, tag, vtx, nb.vtype[vtx], k, diff, vv.is_valid))
    print(lbl, dict(stats), flush=True)
print("true moves (tag, accepted):", dict(true_moves))
print("true moves rejected:", true_rejected[:30])
print("near misses accepted by the template check with a different key:")
for x in bad: print("  ", x)
