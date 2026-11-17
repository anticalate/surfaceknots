import ctypes
import json
import os
import time
from pathlib import Path
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import multiprocessing

os.environ["OMP_NUM_THREADS"] = "1"

from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.ch11_census import _is_admissible
from surfaceknots.moves import (
    canonical_marked_key,
    connected_diagram_components,
    find_planar_connected_sum_decomposition,
    greedy_yoshikawa_reduce,
    list_level_yoshikawa_neighbors,
    list_4d_tait_flypes,
    list_omega1_6_stabilizations,
    list_omega2_pokes,
)
from surfaceknots.move_search import (
    _greedy_path_from,
    _indep_is_single_local_move,
    verify_move_path_independent,
)
from surfaceknots.yoshikawa_data import compute_full_invariant_record, invariant_discrimination_key

MAX_N = 20
VT_TO_INT = {"X0": 0, "X1": 1, "Y": 2, "Z": 3}
INT_TO_VT = ("X0", "X1", "Y", "Z")

SO_PATH = "/afs/csail.mit.edu/u/e/ejli/.gemini/antigravity-cli/brain/c34f4119-775d-42cd-b8c1-b2dd16600960/scratch/_wide_search.so"
_lib = ctypes.CDLL(SO_PATH)


class WideBridgeRecord(ctypes.Structure):
    _fields_ = [
        ("owner_a", ctypes.c_int32),
        ("parent_node_a", ctypes.c_int32),
        ("mid_nv", ctypes.c_int32),
        ("node_b", ctypes.c_int32),
        ("wave_h", ctypes.c_int32),
        ("mid_opp", ctypes.c_int32 * (4 * MAX_N)),
        ("mid_vt", ctypes.c_int32 * MAX_N),
    ]


class WideCollapseRecord(ctypes.Structure):
    _fields_ = [
        ("collapsed", ctypes.c_int32),
        ("parent_node", ctypes.c_int32),
        ("term_nv", ctypes.c_int32),
        ("wave_h", ctypes.c_int32),
        ("term_opp", ctypes.c_int32 * (4 * MAX_N)),
        ("term_vt", ctypes.c_int32 * MAX_N),
    ]


_lib.wide_bidirectional_search_c.argtypes = [
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(WideCollapseRecord),
    ctypes.POINTER(WideBridgeRecord),
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_int,
]
_lib.wide_bidirectional_search_c.restype = ctypes.c_int


def diag_from_c(nv: int, opp_seq, vt_seq) -> MarkedDiagram:
    if nv == 0:
        return MarkedDiagram.unknot_sphere()
    return MarkedDiagram(
        num_vertices=nv,
        opp=tuple(int(opp_seq[i]) for i in range(4 * nv)),
        vtype=tuple(INT_TO_VT[int(vt_seq[i])] for i in range(nv)),
        free_loops=0,
    )


def expand_level_flype_with_paths(root_epd: str, max_states: int = 1024) -> dict[tuple, tuple[MarkedDiagram, list[str]]]:
    """Return {canon_sym_key: (diagram, atomic_epd_path_from_root)} across level + admissible 4D Tait flypes."""
    d0 = MarkedDiagram.from_epd(root_epd)
    r0, steps0 = _greedy_path_from(d0)
    p0 = [d0.to_epd()] + [e for _, e in steps0]
    k0 = canonical_marked_key(r0, up_to_mirror_and_reverse=True)
    vis = {k0: (r0, p0)}
    q = [k0]
    qi = 0
    while qi < len(q) and len(vis) < max_states:
        ck = q[qi]
        qi += 1
        cur_d, cur_path = vis[ck]
        cands = list(list_level_yoshikawa_neighbors(cur_d))
        for fl in list_4d_tait_flypes(cur_d):
            if _is_admissible(fl):
                cands.append(fl)
        for nb in cands:
            if len(connected_diagram_components(nb)) >= 2:
                continue
            red, rsteps = _greedy_path_from(nb)
            if len(connected_diagram_components(red)) >= 2 or red.num_vertices != d0.num_vertices:
                continue
            rk = canonical_marked_key(red, up_to_mirror_and_reverse=True)
            if rk not in vis:
                npath = cur_path + [nb.to_epd()] + [e for _, e in rsteps]
                vis[rk] = (red, npath)
                if len(vis) >= max_states:
                    break
                q.append(rk)
    return vis


def unpack_chains(chains_buf, total_chains: int) -> list[list[MarkedDiagram]]:
    pos = 0
    cnt = int(chains_buf[pos])
    pos += 1
    assert cnt == total_chains
    out = []
    for _ in range(cnt):
        ln = int(chains_buf[pos])
        pos += 1
        if ln <= 0:
            out.append([])
            continue
        ch = []
        for _k in range(ln):
            nv = int(chains_buf[pos])
            pos += 1
            opp = [int(chains_buf[pos + j]) for j in range(4 * nv)]
            pos += 4 * nv
            vt = [int(chains_buf[pos + j]) for j in range(nv)]
            pos += nv
            ch.append(diag_from_c(nv, opp, vt))
        out.append(ch)
    return out


def expand_c_chain_to_atomic_epds(seed_map: dict[tuple, tuple[MarkedDiagram, list[str]]], chain: list[MarkedDiagram]) -> list[str]:
    """Convert a C parent chain `[nd_0, nd_1, ..., nd_k]` into a verified atomic EPD path starting at `root_epd`."""
    assert len(chain) >= 1
    k0 = canonical_marked_key(chain[0], up_to_mirror_and_reverse=True)
    seed_d, path = seed_map[k0]
    epds = list(path)
    if chain[0].to_epd() != epds[-1]:
        epds.append(chain[0].to_epd())

    cur = chain[0]
    for nxt in chain[1:]:
        # Check if cur -> nxt is already a single local move (e.g. Omega_1/6+ or Omega_2+ or pure level move):
        if _indep_is_single_local_move(cur, nxt, allow_o4=True) is not None:
            epds.append(nxt.to_epd())
            cur = nxt
            continue
        # Otherwise cur took 1 level move -> nb, followed by _greedy_path_from(nb) -> red ~ nxt:
        target_ck = canonical_marked_key(nxt, up_to_mirror_and_reverse=True)
        matched = False
        for nb in list_level_yoshikawa_neighbors(cur):
            red, rsteps = _greedy_path_from(nb)
            if red.num_vertices == nxt.num_vertices and canonical_marked_key(red, up_to_mirror_and_reverse=True) == target_ck:
                epds.append(nb.to_epd())
                for _, e_str in rsteps:
                    epds.append(e_str)
                if nxt.to_epd() != epds[-1]:
                    epds.append(nxt.to_epd())
                cur = nxt
                matched = True
                break
        if not matched:
            raise RuntimeError(f"Failed to reconstruct step from {cur.to_epd()} to {nxt.to_epd()}")
    return epds


def run_c_search_on_diagrams(
    n_base: int,
    epd_list: list[str],
    max_headroom: int = 4,
    wave_budget: int = 250000,
    up1_per_owner: int = 512,
    up2_per_owner: int = 384,
    stop_when_merged: int = 1,
    seed_bfs_cap: int = 1024,
):
    num_owners = len(epd_list)
    seed_maps = []
    flat_opp = []
    flat_vt = []
    seed_owner = []

    for o, epd in enumerate(epd_list):
        smap = expand_level_flype_with_paths(epd, max_states=seed_bfs_cap)
        seed_maps.append(smap)
        for sd, _ in smap.values():
            if sd.num_vertices != n_base:
                continue
            flat_opp.extend(sd.opp)
            flat_vt.extend(VT_TO_INT[v] for v in sd.vtype)
            seed_owner.append(o)

    num_seeds = len(seed_owner)
    if num_seeds == 0:
        return None

    c_opp = (ctypes.c_int * len(flat_opp))(*flat_opp)
    c_vt = (ctypes.c_int * len(flat_vt))(*flat_vt)
    c_own = (ctypes.c_int * num_seeds)(*seed_owner)

    out_roots = (ctypes.c_int * num_owners)()
    out_cols = (WideCollapseRecord * num_owners)()
    max_br = 64
    out_brs = (WideBridgeRecord * max_br)()
    out_num_br = ctypes.c_int(0)
    out_tot_nodes = ctypes.c_int(0)
    max_chains_ints = 262144
    out_chains = (ctypes.c_int * max_chains_ints)()

    rc = _lib.wide_bidirectional_search_c(
        n_base,
        num_owners,
        num_seeds,
        c_opp,
        c_vt,
        c_own,
        max_headroom,
        wave_budget,
        up1_per_owner,
        up2_per_owner,
        stop_when_merged,
        out_roots,
        out_cols,
        out_brs,
        max_br,
        ctypes.byref(out_num_br),
        ctypes.byref(out_tot_nodes),
        out_chains,
        max_chains_ints,
    )
    assert rc == 0, f"C search failed with rc={rc}"
    nbr = out_num_br.value
    chains = unpack_chains(out_chains, num_owners + 2 * nbr)

    # Verify bridges independently!
    verified_bridges = []
    for b in range(nbr):
        br = out_brs[b]
        chain_a = chains[num_owners + 2 * b]
        chain_b = chains[num_owners + 2 * b + 1]
        mid_d = diag_from_c(br.mid_nv, br.mid_opp, br.mid_vt)
        own_b = None
        k_b0 = canonical_marked_key(chain_b[0], up_to_mirror_and_reverse=True)
        for o_idx, sm in enumerate(seed_maps):
            if k_b0 in sm:
                own_b = o_idx
                break

        if br.parent_node_a == -1:
            # Bridged directly during seed insertion (level + 4D flype)!
            ka = canonical_marked_key(mid_d, up_to_mirror_and_reverse=True)
            epds_a = list(seed_maps[br.owner_a][ka][1])
        else:
            epds_a = expand_c_chain_to_atomic_epds(seed_maps[br.owner_a], chain_a + [mid_d])

        epds_b = expand_c_chain_to_atomic_epds(seed_maps[own_b], chain_b)
        # Join epds_a -> epds_b[::-1]
        full_path = list(epds_a)
        for rev_e in reversed(epds_b):
            if rev_e != full_path[-1]:
                full_path.append(rev_e)
        ver = verify_move_path_independent(full_path, allow_o4=True)
        max_v_on_path = max(MarkedDiagram.from_epd(x).num_vertices for x in full_path)
        verified_bridges.append({
            "owner_a": int(br.owner_a),
            "owner_b": int(own_b),
            "wave_h": int(br.wave_h),
            "max_ch_on_path": max_v_on_path,
            "num_steps": len(full_path) - 1,
            "verified": bool(ver.is_valid),
            "error": ver.error,
            "epd_path": full_path,
        })

    collapses = []
    for o in range(num_owners):
        col = out_cols[o]
        if col.collapsed != 0:
            term_d = diag_from_c(col.term_nv, col.term_opp, col.term_vt)
            ch_o = chains[o]
            ver_ok = None
            steps_cnt = 0
            if ch_o and col.collapsed in (1, 3, 4, 5) and col.term_nv > 0:
                try:
                    p_o = expand_c_chain_to_atomic_epds(seed_maps[o], ch_o + [term_d])
                    v_res = verify_move_path_independent(p_o, allow_o4=True)
                    ver_ok = bool(v_res.is_valid)
                    steps_cnt = len(p_o) - 1
                except Exception as exc:
                    ver_ok = f"err:{exc}"
            collapses.append({
                "owner": o,
                "kind": int(col.collapsed),  # 1=<n, 2=split, 3=2cut@n, 4=2cut@(n+1), 5=2cut@(n+2)
                "wave_h": int(col.wave_h),
                "term_nv": int(col.term_nv),
                "verified": ver_ok,
                "num_steps": steps_cnt,
            })

    roots = [int(out_roots[o]) for o in range(num_owners)]
    return {
        "num_owners": num_owners,
        "num_seeds": num_seeds,
        "total_nodes_explored": int(out_tot_nodes.value),
        "dsu_roots": roots,
        "num_distinct_roots": len(set(roots)),
        "collapses": collapses,
        "bridges": verified_bridges,
    }


def _packet_worker(args):
    ch, label, epds = args
    t0 = time.perf_counter()
    res = run_c_search_on_diagrams(
        n_base=ch,
        epd_list=epds,
        max_headroom=4,          # up to ch+4 (15 for ch=11, 16 for ch=12)!
        wave_budget=300000,      # up to 1.2M states across waves!
        up1_per_owner=512,
        up2_per_owner=384,
        stop_when_merged=1,
        seed_bfs_cap=2048,
    )
    res["ch"] = ch
    res["label"] = label
    res["runtime_sec"] = round(time.perf_counter() - t0, 3)
    return res


def _residual_cell_worker(args):
    ch, cell_idx, tag, epds = args
    t0 = time.perf_counter()
    # Cap headroom=4 for cells with <= 8 diagrams (includes ALL non-unknot cells!), headroom=3 for huge unknot cells
    mh = 4 if len(epds) <= 8 else 3
    wb = 180000 if len(epds) <= 8 else 120000
    # Run each diagram individually OR in batches of <= 8 so no large unknot cell dilutes per-diagram budget!
    collapsed_owners = {}
    all_roots_before = len(epds)
    # First run cell-together if len(epds) <= 12 (catches both collapses AND intra-cell merges!):
    if len(epds) <= 12:
        res = run_c_search_on_diagrams(
            n_base=ch,
            epd_list=epds,
            max_headroom=mh,
            wave_budget=wb,
            up1_per_owner=384,
            up2_per_owner=256,
            stop_when_merged=0,  # keep going to collapse all owners if possible!
            seed_bfs_cap=512,
        )
        for c in res["collapses"]:
            collapsed_owners[c["owner"]] = c
        # Propagate collapse across merged DSU components!
        roots = res["dsu_roots"]
        col_roots = {roots[o] for o in collapsed_owners}
        alive_roots = {roots[o] for o in range(len(epds)) if roots[o] not in col_roots}
        alive_diags = sum(1 for o in range(len(epds)) if roots[o] not in col_roots)
        return {
            "ch": ch,
            "cell_idx": cell_idx,
            "tag": tag,
            "input_diags": len(epds),
            "collapsed_diags": len(epds) - alive_diags,
            "surviving_diags": alive_diags,
            "surviving_components": len(alive_roots),
            "collapses_detail": res["collapses"],
            "bridges_count": len(res["bridges"]),
            "bridges_all_verified": all(b["verified"] is True for b in res["bridges"]),
            "runtime_sec": round(time.perf_counter() - t0, 3),
        }
    else:
        # Large unknot cell: test in batches of 8 so every diagram gets rich mountain-pass budget
        surv_diags = 0
        surv_comps = 0
        col_details = []
        br_cnt = 0
        br_ver = True
        for start in range(0, len(epds), 8):
            chunk = epds[start : start + 8]
            res = run_c_search_on_diagrams(
                n_base=ch,
                epd_list=chunk,
                max_headroom=3,
                wave_budget=140000,
                up1_per_owner=256,
                up2_per_owner=192,
                stop_when_merged=0,
                seed_bfs_cap=384,
            )
            c_owns = {c["owner"] for c in res["collapses"]}
            col_details.extend(res["collapses"])
            roots = res["dsu_roots"]
            c_roots = {roots[o] for o in c_owns}
            a_roots = {roots[o] for o in range(len(chunk)) if roots[o] not in c_roots}
            a_diags = sum(1 for o in range(len(chunk)) if roots[o] not in c_roots)
            surv_diags += a_diags
            surv_comps += len(a_roots)
            br_cnt += len(res["bridges"])
            if not all(b["verified"] is True for b in res["bridges"]):
                br_ver = False
        return {
            "ch": ch,
            "cell_idx": cell_idx,
            "tag": tag,
            "input_diags": len(epds),
            "collapsed_diags": len(epds) - surv_diags,
            "surviving_diags": surv_diags,
            "surviving_components": surv_comps,
            "collapses_detail": col_details,
            "bridges_count": br_cnt,
            "bridges_all_verified": br_ver,
            "runtime_sec": round(time.perf_counter() - t0, 3),
        }


if __name__ == "__main__":
    DATA_DIR = Path("code/data")
    SCRATCH_DIR = Path("/afs/csail.mit.edu/u/e/ejli/.gemini/antigravity-cli/brain/c34f4119-775d-42cd-b8c1-b2dd16600960/scratch")
    ch11_data = json.load((DATA_DIR / "expanded_ch11.json").open())
    ch12_data = json.load((DATA_DIR / "expanded_ch12.json").open())

    # =========================================================================
    # PART 1: ALL 15 NEW_PACKET MULTI-DIAGRAM CLASSES (up to +4 headroom!)
    # =========================================================================
    packet_tasks = []
    for e in ch11_data["entries"]:
        if e["cell_status"] == "NEW_PACKET":
            packet_tasks.append((11, e["label"], e["move_component_epds"]))
    for e in ch12_data["entries"]:
        if e["cell_status"] == "NEW_PACKET":
            packet_tasks.append((12, e["label"], e["move_component_epds"]))

    print(f"Running wide +1..+4 bidirectional C mountain-pass search on ALL {len(packet_tasks)} NEW_PACKET classes...")
    t_p0 = time.perf_counter()
    with ProcessPoolExecutor(max_workers=min(16, len(packet_tasks)), mp_context=multiprocessing.get_context("forkserver")) as pool:
        pkt_results = list(pool.map(_packet_worker, packet_tasks))

    print(f"\n=== RESULTS ON ALL {len(pkt_results)} MULTI-DIAGRAM NEW_PACKET CLASSES ({time.perf_counter()-t_p0:.1f}s) ===")
    total_before = 0
    total_after = 0
    completely_merged_packets = 0
    partially_merged_packets = 0
    for pr in pkt_results:
        bef = pr["num_owners"]
        aft = pr["num_distinct_roots"]
        total_before += bef
        total_after += aft
        if aft == 1:
            completely_merged_packets += 1
        elif aft < bef:
            partially_merged_packets += 1
        br_summary = [
            f"({b['owner_a']}<->{b['owner_b']} @h=+{b['wave_h']} max_ch={b['max_ch_on_path']} steps={b['num_steps']} ver={b['verified']})"
            for b in pr["bridges"]
        ]
        col_summary = [f"(own={c['owner']} kind={c['kind']} @h=+{c['wave_h']})" for c in pr["collapses"]]
        print(
            f"  ch={pr['ch']} {pr['label']:10s}: {bef} -> {aft} comps "
            f"[nodes={pr['total_nodes_explored']:7d}, {pr['runtime_sec']:5.2f}s] "
            f"bridges={br_summary} collapses={col_summary}"
        )
    print(
        f"  SUMMARY PACKETS: {len(pkt_results)} packets ({total_before} -> {total_after} move components; "
        f"fully merged to singleton: {completely_merged_packets}, partially merged: {partially_merged_packets})"
    )
    (SCRATCH_DIR / "wide_packet_results.json").write_text(json.dumps(pkt_results, indent=2))

    # =========================================================================
    # PART 2: ALL SURVIVING ch=11 AND ch=12 RESIDUALS (grouped by invariant cell!)
    # =========================================================================
    for ch, raw_data in [(11, ch11_data), (12, ch12_data)]:
        st1 = json.load((SCRATCH_DIR / f"residual_stage1_ch{ch}.json").open())
        surv_indices = [item["idx"] for item in st1 if item["status"] == "survives"]
        by_cell = defaultdict(list)
        cell_tag = {}
        for idx in surv_indices:
            rec = raw_data["residual_collisions"][idx]
            dk = invariant_discrimination_key(rec)
            by_cell[dk].append(rec["epd"])
            cell_tag[dk] = " | ".join(rec.get("collision_matches", [])[:2])

        cell_tasks = [
            (ch, c_i, cell_tag[dk], epds)
            for c_i, (dk, epds) in enumerate(sorted(by_cell.items(), key=lambda kv: (len(kv[1]), cell_tag[kv[0]])))
        ]
        print(f"\nRunning wide mountain-pass search on ch={ch} Stage-1 surviving residuals ({len(surv_indices)} diags across {len(cell_tasks)} cells)...")
        t_r0 = time.perf_counter()
        with ProcessPoolExecutor(max_workers=16, mp_context=multiprocessing.get_context("forkserver")) as pool:
            cell_results = list(pool.map(_residual_cell_worker, cell_tasks))

        in_tot = sum(cr["input_diags"] for cr in cell_results)
        col_tot = sum(cr["collapsed_diags"] for cr in cell_results)
        surv_tot = sum(cr["surviving_diags"] for cr in cell_results)
        comp_tot = sum(cr["surviving_components"] for cr in cell_results)
        cells_eliminated = sum(1 for cr in cell_results if cr["surviving_diags"] == 0)
        print(
            f"=== WIDE MOUNTAIN-PASS RESIDUAL RESULTS ch={ch} ({time.perf_counter()-t_r0:.1f}s) ===\n"
            f"  Stage-1 surviving input : {in_tot} diagrams ({len(cell_results)} cells)\n"
            f"  Collapsed by wide +1..+4: {col_tot} diagrams ({cells_eliminated} entire cells eliminated!)\n"
            f"  Remaining surviving     : {surv_tot} diagrams ({comp_tot} unique move components across {len(cell_results)-cells_eliminated} cells)\n"
            f"  All discovered bridges independently verified: {all(cr['bridges_all_verified'] for cr in cell_results)}"
        )
        for cr in sorted(cell_results, key=lambda x: (-x["collapsed_diags"], -x["input_diags"])):
            if cr["collapsed_diags"] > 0 or cr["surviving_components"] < cr["surviving_diags"] or cr["input_diags"] <= 5:
                print(
                    f"    in={cr['input_diags']:3d} -> surv={cr['surviving_diags']:3d} (comps={cr['surviving_components']:3d}, "
                    f"collapsed={cr['collapsed_diags']:3d}, bridges={cr['bridges_count']:2d}) | {cr['tag']}"
                )
        (SCRATCH_DIR / f"wide_residual_results_ch{ch}.json").write_text(json.dumps(cell_results, indent=2))
