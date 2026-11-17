import json
from collections import Counter, defaultdict

def categorize_match(m: str) -> str:
    # m is a reference label from build_reference_catalog / build_def5_1_high_composite_matches
    if " u " in m:
        return "split_union"
    if " # " in m:
        return "connected_sum"
    # Prime or unknot:
    # Unknots: 0_1, 2_1^1, 2_1^-1 (or P+, P-)
    if m in ("0_1", "2_1^1", "2_1^-1", "P+", "P-"):
        return "unknot"
    # Check ch prefix
    head = m.split("_")[0]
    if head.isdigit():
        ch_val = int(head)
        if ch_val <= 7:
            # 6_1^0,1, 7_1^0,-2 are prime links!
            return f"prime_ch<=10"
        elif ch_val <= 10:
            return f"prime_ch<=10"
        elif ch_val == 11:
            return "prime_ch=11"
    return f"other({m})"

for ch, fn in [(11, "code/data/expanded_ch11.json"), (12, "code/data/expanded_ch12.json")]:
    with open(fn) as f:
        data = json.load(f)
    print(f"\n==========================================")
    print(f"CH = {ch} SUMMARY")
    print(f"==========================================")
    print("census_stats:", json.dumps(data["census_stats"], indent=2))
    print("post_filter_stats:", json.dumps(data["post_filter_stats"], indent=2))
    entries = data["entries"]
    res = data["residual_collisions"]
    print(f"entries: {len(entries)} (singletons={data['exact_singletons_count']}, packets={data['new_packets_count']}, total_diagrams={data['certified_new_prime_diagrams_total']})")
    print(f"residual_collisions: {len(res)} across {data['ambiguous_shadow_cells_count']} cells")

    print(f"\n--- NEW_PACKET entries at ch={ch} ---")
    for e in entries:
        if e["cell_status"] == "NEW_PACKET" or e["num_move_components"] > 1:
            epds = e.get("move_component_epds") or ([e["epd"]] + e.get("alternate_epds", []))
            from surfaceknots.diagrams import MarkedDiagram
            cm_list = []
            for ep in epds:
                d_obj = MarkedDiagram.from_epd(ep)
                cm_list.append((d_obj.num_classical, d_obj.num_marked))
            print(f"  {e['label']:16s} | comps={e['num_move_components']:2d} | topo={' u '.join(e['component_types']):16s} | e={e['normal_euler_numbers']} | (c,m)={Counter(cm_list)}")

    # Group residuals by cell (using invariant cell or exact collision_matches)
    # Wait: how does ch11/ch12_census count ambiguous_shadow_cells_count?
    # By distinct invariant cell key! Let's check collision_matches per diagram and per cell.
    by_matches = defaultdict(list)
    cat_diag_counts = Counter()
    exact_match_counts = Counter()
    for r in res:
        ms = tuple(r["collision_matches"])
        by_matches[ms].append(r)
        cats = tuple(sorted({categorize_match(m) for m in ms}))
        cat_diag_counts[cats] += 1
        exact_match_counts[ms] += 1

    print(f"\n--- Residual breakdown by category tuple (total={len(res)} diagrams) ---")
    for cats, cnt in cat_diag_counts.most_common():
        print(f"  {str(cats):45s}: {cnt:4d} diagrams")

    print(f"\n--- Residual breakdown by exact collision_matches ({len(exact_match_counts)} distinct match-tuples) ---")
    for ms, cnt in exact_match_counts.most_common():
        rs = by_matches[ms]
        topos = Counter(" u ".join(x["component_types"]) for x in rs)
        cms = Counter((x["num_classical"], x["num_marked"]) for x in rs)
        print(f"  matches={list(ms)} | count={cnt:3d} | topos={dict(topos)} | (c,m)={dict(cms)}")
