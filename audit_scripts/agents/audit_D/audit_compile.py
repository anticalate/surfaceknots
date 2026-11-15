"""Sanity check: compiled MarkedDiagram vs tri-plane invariants on a sample of compilable concentrated triples."""
import random, sys, json
from collections import Counter
from surfaceknots.triplane import (TangleDiagram, TriPlaneDiagram, compute_triplane_invariant_record,
    concentrated_triplane_to_marked_diagram, triplane_o4_key, generate_catalan_matchings)
from surfaceknots.triplane_census import enumerate_minimal_braided_tangles
from surfaceknots.yoshikawa_data import compute_full_invariant_record
random.seed(7)
def cmp(tp):
    md = concentrated_triplane_to_marked_diagram(tp)
    if md is None: return None
    r_tp = compute_triplane_invariant_record(tp)
    try:
        r_md = compute_full_invariant_record(md, label="c", include_sublinks=False)
    except AssertionError as ex:
        return ("INADMISSIBLE", str(ex)[:60])
    fields = ["num_components", "component_types", "orientability", "euler_characteristic", "abs_normal_euler_multiset", "h1_free_rank", "h1_torsion",
              "e1_fingerprint_hash", "hom_counts", "sigma2_h1_torsion", "sigma2_b1", "sigma2_b2", "sigma2_abs_signature"]
    diffs = []
    for f in fields:
        a = r_tp[f]; b_ = r_md.get(f)
        if f == "abs_normal_euler_multiset": b_ = sorted(abs(x) for x in r_md["normal_euler_numbers"])
        if f == "component_types": a, b_ = sorted(a), sorted(b_)
        if a != b_: diffs.append((f, a, b_))
    ch_ok = md.ch_index == tp.crossing_number + tp.bridge_number - (md.num_marked and (tp.bridge_number - md.num_marked) or tp.bridge_number)
    return ("OK" if not diffs else "MISMATCH", diffs, tuple(sorted(r_tp["component_types"])), md.num_classical, md.num_marked, triplane_o4_key(r_tp) == triplane_o4_key(r_md))
res = Counter(); classes = Counter(); bad = []
n_tested = 0
for b in (2, 3, 4):
    by_c, _ = enumerate_minimal_braided_tangles(b, 4)
    cats = [TangleDiagram(b, m, ()) for m in generate_catalan_matchings(b)]
    pool = [t for c in range(0, 5) for t in by_c[c]]
    tries = 0
    target = {2: 150, 3: 300, 4: 500}[b]
    got = 0
    while got < target and tries < 400000:
        tries += 1
        t2 = random.choice(pool); m1 = random.choice(cats); m3 = random.choice(cats)
        tp = TriPlaneDiagram(m1, t2, m3)
        if not tp.is_certified_valid(): continue
        out = cmp(tp)
        if out is None: res["not_compilable"] += 1; continue
        got += 1; n_tested += 1
        res[out[0]] += 1
        if out[0] == "OK":
            classes[(b, out[2])] += 1
            assert out[3] == t2.num_crossings
            if not out[5]: res["o4_key_differs"] += 1
        else:
            bad.append((tp, out))
# census witnesses (incl. knotted classes) from the committed JSON
d = json.load(open("data/complexity/triplane_b4_census.json"))
def mk(s): return TangleDiagram(s["num_bridges"], tuple(tuple(p) for p in s["bottom_matching"]), tuple(s["braid_word"]))
for e in d["entries"]:
    for key in ("witness_triplane", "witness_concentrated"):
        w = e.get(key)
        if not w: continue
        tp = TriPlaneDiagram(mk(w["t1"]), mk(w["t2"]), mk(w["t3"]))
        out = cmp(tp)
        if out is None: res["json_witness_not_compilable"] += 1; continue
        res["json_" + out[0]] += 1
        if out[0] != "OK": bad.append((tp, out))
    if e["best_compiled_epd"]:
        from surfaceknots.diagrams import MarkedDiagram
        try:
            md = MarkedDiagram.from_epd(e["best_compiled_epd"])
        except ValueError as ex:
            res["json_best_epd_unparseable"] += 1; print("  unparseable best_compiled_epd:", e["best_compiled_epd"], e["catalog_matches"]); continue
        r_md = compute_full_invariant_record(md, label="c", include_sublinks=False)
        same = triplane_o4_key(r_md) == triplane_o4_key(e["invariants"])
        res["json_best_epd_o4_match" if same else "json_best_epd_o4_MISMATCH"] += 1
        if not same: print("  best_epd mismatch", e["catalog_matches"], r_md["component_types"], r_md["normal_euler_numbers"], e["invariants"]["normal_euler_numbers"])
print(dict(res)); print(sorted(classes.items(), key=str))
for tp, out in bad[:10]: print("BAD", tp, out)
