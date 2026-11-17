import json
from collections import Counter, defaultdict
from surfaceknots.yoshikawa_data import invariant_discrimination_key

UNKNOT_PRIMES = {"ch<=10:0_1", "ch<=10:2_1^1", "ch<=10:2_1^-1"}

def parse_single_tag(tag: str) -> str:
    if tag.startswith("unknot:") or tag in UNKNOT_PRIMES:
        return "UNKNOT"
    if tag.startswith("ch<=10:"):
        return "PRIME_CH<=10"
    if tag.startswith("ch=11:"):
        return "PRIME_CH=11"
    if tag.startswith("composite:"):
        return "COMPOSITE"
    if tag.startswith("split:"):
        return "SPLIT"
    return f"OTHER({tag})"

def is_viro_equivalent_to_unknot(matches: list[str]) -> bool:
    has_unk = any(m.startswith("unknot:") or m in UNKNOT_PRIMES for m in matches)
    if not has_unk:
        return False
    # Check if every composite in matches is 10_3 summed with non-orientable unknots
    for m in matches:
        if m.startswith("composite:"):
            # Viro: 10_3 # RP^2 ~= RP^2, hence 10_3 # (#k RP^2) ~= #k RP^2 for any k >= 1
            if "10_3" not in m:
                return False
    return True

for ch, fn in [(11, "code/data/expanded_ch11.json"), (12, "code/data/expanded_ch12.json")]:
    with open(fn) as f:
        data = json.load(f)
    res = data["residual_collisions"]
    # Group by invariant_discrimination_key to verify exact cell count == ambiguous_shadow_cells_count!
    cells = defaultdict(list)
    for r in res:
        ik = invariant_discrimination_key(r)
        cells[ik].append(r)
    print(f"=== CH = {ch}: {len(res)} residual diagrams in {len(cells)} invariant cells (JSON ambiguous_shadow_cells_count = {data['ambiguous_shadow_cells_count']}) ===")

    # High-level mutually exclusive classification of cells / diagrams
    bucket_cells = Counter()
    bucket_diags = Counter()
    sub_unknot_diags = Counter()
    sub_prime10_diags = Counter()
    sub_prime11_diags = Counter()
    sub_comp_diags = Counter()
    sub_split_diags = Counter()

    rows = []
    for ik, rlist in sorted(cells.items(), key=lambda kv: (-len(kv[1]), str(kv[1][0]["collision_matches"]))):
        ms = rlist[0]["collision_matches"]
        assert all(x["collision_matches"] == ms for x in rlist), "Mismatch inside cell!"
        cats = sorted({parse_single_tag(m) for m in ms})
        if cats == ["UNKNOT"]:
            primary = "1_UNKNOT_PURE"
        elif is_viro_equivalent_to_unknot(ms):
            primary = "1b_UNKNOT_PLUS_VIRO_10_3_COLLAPSE"
        elif cats == ["PRIME_CH<=10"]:
            primary = "2_PRIME_CH<=10"
        elif cats == ["PRIME_CH=11"]:
            primary = "3_PRIME_CH=11"
        elif cats == ["COMPOSITE"]:
            primary = "4_CONNECTED_SUM"
        elif cats == ["SPLIT"]:
            primary = "5_SPLIT_UNION"
        else:
            primary = "6_MIXED(" + "+".join(cats) + ")"

        bucket_cells[primary] += 1
        bucket_diags[primary] += len(rlist)

        if primary.startswith("1"):
            unk_tags = [m for m in ms if m.startswith("unknot:") or m in UNKNOT_PRIMES]
            sub_unknot_diags[tuple(unk_tags)] += len(rlist)
        elif primary == "2_PRIME_CH<=10":
            sub_prime10_diags[tuple(ms)] += len(rlist)
        elif primary == "3_PRIME_CH=11":
            sub_prime11_diags[tuple(ms)] += len(rlist)
        elif primary == "4_CONNECTED_SUM":
            sub_comp_diags[tuple(ms[:3]) + ((f"...+{len(ms)-3} more",) if len(ms) > 3 else ())] += len(rlist)
        elif primary == "5_SPLIT_UNION":
            sub_split_diags[tuple(ms)] += len(rlist)

        cms = Counter((x["num_classical"], x["num_marked"]) for x in rlist)
        topo = " u ".join(rlist[0]["component_types"])
        e_nums = rlist[0]["normal_euler_numbers"]
        short_ms = ms if len(ms) <= 4 else (ms[:3] + [f"...({len(ms)} total)"])
        rows.append((len(rlist), primary, topo, e_nums, dict(sorted(cms.items())), short_ms))

    print("\nSummary by Primary Bucket:")
    for k in sorted(bucket_diags.keys()):
        print(f"  {k:38s}: {bucket_cells[k]:3d} cells | {bucket_diags[k]:4d} diagrams")
    print(f"  {'TOTAL':38s}: {sum(bucket_cells.values()):3d} cells | {sum(bucket_diags.values()):4d} diagrams")

    print("\nDetailed Cell Table:")
    for idx, (cnt, prim, topo, e_nums, cms, sms) in enumerate(rows, 1):
        print(f"  Cell {idx:2d} | diags={cnt:3d} | {prim:34s} | topo={topo:14s} e={str(e_nums):10s} | (c,m)={cms} | matches={sms}")
    print()
