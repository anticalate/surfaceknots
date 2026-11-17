import json
from pathlib import Path
from collections import Counter, defaultdict
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.move_search import verify_move_path_independent
from surfaceknots.yoshikawa_data import invariant_discrimination_key

DATA_DIR = Path("code/data")
SCRATCH = Path("/afs/csail.mit.edu/u/e/ejli/.gemini/antigravity-cli/brain/c34f4119-775d-42cd-b8c1-b2dd16600960/scratch")

ch11_data = json.load((DATA_DIR / "expanded_ch11.json").open())
ch12_data = json.load((DATA_DIR / "expanded_ch12.json").open())
pkt_res = json.load((SCRATCH / "wide_packet_results.json").open())

print("=== DETAILED BRIDGE INSPECTION FOR MERGED PACKETS ===")
for pr in pkt_res:
    if pr["bridges"]:
        print(f"\nPacket ch={pr['ch']} {pr['label']} ({pr['num_owners']} -> {pr['num_distinct_roots']}):")
        for b in pr["bridges"]:
            path = b["epd_path"]
            ver = verify_move_path_independent(path, allow_o4=True)
            cm_trace = []
            for epd in path:
                d = MarkedDiagram.from_epd(epd)
                cm_trace.append((d.num_vertices, d.num_classical, d.num_marked))
            print(f"  Bridge {b['owner_a']} <-> {b['owner_b']}: wave=+{b['wave_h']}, max_v={b['max_ch_on_path']}, steps={b['num_steps']}, indep_valid={ver.is_valid}")
            print(f"    step_kinds: {Counter(ver.step_kinds)}")
            print(f"    (v, c, m) profile: {cm_trace}")

for ch, raw_data in [(11, ch11_data), (12, ch12_data)]:
    st1 = json.load((SCRATCH / f"residual_stage1_ch{ch}.json").open())
    st2 = json.load((SCRATCH / f"wide_residual_results_ch{ch}.json").open())

    # Count original cells by category
    def classify_cell(matches):
        has_unknot = any(m.startswith("unknot:") or m in ("ch<=10:0_1", "ch<=10:2_1^1", "ch<=10:2_1^-1") for m in matches)
        has_ch10 = any(m.startswith("ch<=10:") and m not in ("ch<=10:0_1", "ch<=10:2_1^1", "ch<=10:2_1^-1") for m in matches)
        has_ch11 = any(m.startswith("ch=11:") for m in matches)
        has_comp = any(m.startswith("composite") for m in matches)
        has_split = any(m.startswith("split:") for m in matches)
        if has_unknot:
            return "Unknot (incl. Viro 10_3 collapse)" if has_comp else "Pure Unknot"
        if has_ch10:
            return "ch <= 10 Prime"
        if has_ch11:
            return "ch = 11 Prime"
        if has_split:
            return "Split Union"
        if has_comp:
            return "Connected Sum K_1 # K_2"
        return "Other"

    orig_by_dk = defaultdict(list)
    cat_by_dk = {}
    tag_by_dk = {}
    for idx, rec in enumerate(raw_data["residual_collisions"]):
        dk = invariant_discrimination_key(rec)
        orig_by_dk[dk].append(idx)
        cat_by_dk[dk] = classify_cell(rec.get("collision_matches", []))
        tag_by_dk[dk] = " | ".join(rec.get("collision_matches", [])[:2])

    # Map Stage 1 survivals by dk
    st1_by_idx = {x["idx"]: x for x in st1}
    # Map Stage 2 results: note st2 has tag and input_diags, let's match by exact dk order!
    surv_indices = [item["idx"] for item in st1 if item["status"] == "survives"]
    by_cell_s1 = defaultdict(list)
    for idx in surv_indices:
        rec = raw_data["residual_collisions"][idx]
        dk = invariant_discrimination_key(rec)
        by_cell_s1[dk].append(rec["epd"])
    ordered_s1_dks = [dk for dk, _ in sorted(by_cell_s1.items(), key=lambda kv: (len(kv[1]), tag_by_dk[kv[0]]))]
    assert len(ordered_s1_dks) == len(st2)
    st2_by_dk = dict(zip(ordered_s1_dks, st2))

    cat_summary = defaultdict(lambda: [0, 0, 0, 0, 0, 0])  # [orig_cells, orig_diags, s1_surv, s2_surv_diags, s2_surv_comps, s2_surv_cells]
    wave_counts = Counter()
    kind_counts = Counter()
    ver_col_cnt = Counter()
    for cr in st2:
        for cd in cr["collapses_detail"]:
            wave_counts[cd["wave_h"]] += 1
            kind_counts[cd["kind"]] += 1
            ver_col_cnt[str(cd["verified"])] += 1

    for dk, idxs in orig_by_dk.items():
        cat = cat_by_dk[dk]
        cat_summary[cat][0] += 1
        cat_summary[cat][1] += len(idxs)
        s1_alive = sum(1 for i in idxs if st1_by_idx[i]["status"] == "survives")
        cat_summary[cat][2] += s1_alive
        if dk in st2_by_dk:
            cr = st2_by_dk[dk]
            cat_summary[cat][3] += cr["surviving_diags"]
            cat_summary[cat][4] += cr["surviving_components"]
            if cr["surviving_diags"] > 0:
                cat_summary[cat][5] += 1

    print(f"\n=== FULL CATEGORY SUMMARY ch={ch} ===")
    tot = [0] * 6
    for cat, vals in sorted(cat_summary.items(), key=lambda kv: -kv[1][1]):
        for j in range(6):
            tot[j] += vals[j]
        print(
            f"  {cat:34s}: orig={vals[1]:3d} diags ({vals[0]:2d} cells) -> "
            f"after_S1={vals[2]:3d} -> after_S2={vals[3]:3d} diags / {vals[4]:3d} comps ({vals[5]:2d} cells)"
        )
    print(
        f"  {'TOTAL':34s}: orig={tot[1]:3d} diags ({tot[0]:2d} cells) -> "
        f"after_S1={tot[2]:3d} -> after_S2={tot[3]:3d} diags / {tot[4]:3d} comps ({tot[5]:2d} cells)"
    )
    print(f"  Stage-2 C mountain-pass collapse waves (+h): {dict(sorted(wave_counts.items()))}")
    print(f"  Stage-2 C mountain-pass collapse kinds (1=<n, 2=split, 3=2cut@n, 4=2cut@n+1, 5=2cut@n+2): {dict(sorted(kind_counts.items()))}")
    print(f"  Stage-2 C collapse path independent verification: {dict(ver_col_cnt)}")
