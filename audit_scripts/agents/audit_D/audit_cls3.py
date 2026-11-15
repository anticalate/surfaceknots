import json
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.yoshikawa_data import compute_full_invariant_record
d = json.load(open("data/complexity/targeted_ch12_13_census.json"))
for e in d["ch12_m2_s2_knots"]["new_prime_classes"]:
    if e["num_move_components"] < 2: continue
    recs = [compute_full_invariant_record(MarkedDiagram.from_epd(x), label="x", include_sublinks=False) for x in e["all_epds"]]
    skip = {"label", "epd", "pi1_presentation", "symmetry", "ch_index_diagram", "bridge_bounds", "sl2f3_peripheral_sample", "pi1_generators", "pi1_relators", "num_classical", "num_marked"}
    diff = {k for k in recs[0] if k not in skip and any(str(r[k]) != str(recs[0][k]) for r in recs)}
    print(e["class_id"], "components:", len(recs), "| fields differing between move-components:", sorted(diff), "| pi1 presentations:", len({r["pi1_presentation"] for r in recs}), "| symmetry summaries:", {str(r["symmetry"].get("summary")) for r in recs})
    for k in sorted(diff):
        print("    ", k, [str(r[k])[:60] for r in recs])
