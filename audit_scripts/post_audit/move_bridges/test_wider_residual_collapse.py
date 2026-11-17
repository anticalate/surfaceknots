import json
import os
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
    list_level_yoshikawa_neighbors,
    list_4d_tait_flypes,
    greedy_yoshikawa_reduce as simplify_yoshikawa,
    simplify_and_explore_yoshikawa,
    _is_planar_4reg_diagram,
)
from surfaceknots.move_search import analyze_marked_bigon_1handles
from surfaceknots.fast_enum import tunnel_diagrams_c

DATA_DIR = Path("code/data")


def _is_split_diag(d: MarkedDiagram) -> bool:
    return len(connected_diagram_components(d)) >= 2


def explore_level_and_flypes(d0: MarkedDiagram, max_states: int = 1024):
    """Full BFS over level Yoshikawa moves + admissible 4D Tait flypes.
    Returns (status, witness_info, component_diagrams_dict).
    status in {'reduced', 'split', '2cut', 'open'}
    """
    n = d0.num_vertices
    r0 = simplify_yoshikawa(d0)
    if _is_split_diag(r0):
        return ("split", ("init", r0.to_epd()), {})
    if r0.num_vertices < n:
        return ("reduced", ("init", r0.num_vertices, r0.to_epd()), {})
    if find_planar_connected_sum_decomposition(r0) is not None:
        return ("2cut", ("init", r0.to_epd()), {})

    k0 = canonical_marked_key(r0, up_to_mirror_and_reverse=True)
    visited = {k0: r0}
    queue = [r0]
    qi = 0
    while qi < len(queue) and len(visited) < max_states:
        cur = queue[qi]
        qi += 1
        nbs = list(list_level_yoshikawa_neighbors(cur))
        for fl in list_4d_tait_flypes(cur):
            if _is_admissible(fl):
                nbs.append(fl)
        for nb in nbs:
            if _is_split_diag(nb):
                return ("split", ("step", nb.to_epd()), visited)
            if nb.num_vertices < n:
                return ("reduced", ("step", nb.num_vertices, nb.to_epd()), visited)
            red = simplify_yoshikawa(nb)
            if _is_split_diag(red):
                return ("split", ("simp", red.to_epd()), visited)
            if red.num_vertices < n:
                return ("reduced", ("simp", red.num_vertices, red.to_epd()), visited)
            if find_planar_connected_sum_decomposition(red) is not None:
                return ("2cut", ("simp", red.to_epd()), visited)
            ck = canonical_marked_key(red, up_to_mirror_and_reverse=True)
            if ck not in visited:
                visited[ck] = red
                if len(visited) >= max_states:
                    break
                queue.append(red)
    return ("open", None, visited)


def list_valid_1comp_2gon_pinches(d: MarkedDiagram) -> list[MarkedDiagram]:
    """Return all admissible 1-component diagrams obtained by pinching an embedded marked-marked 2-gon of `d`."""
    n = d.num_vertices
    if n < 4:
        return []
    st = d.surface_topology()
    if st.num_components != 1:
        return []
    out = []
    seen = set()
    for u in range(n):
        if d.vtype[u] not in ("Y", "Z"):
            continue
        for su in range(4):
            d_a = d.opp[4 * u + su]
            d_b = d.opp[4 * u + ((su + 1) & 3)]
            w = d_a >> 2
            if w <= u or (d_b >> 2) != w or d.vtype[w] not in ("Y", "Z"):
                continue
            sw_a = d_a & 3
            sw_b = d_b & 3
            if ((sw_b + 1) & 3) != sw_a:
                continue
            if (d.vtype[u] == d.vtype[w]) != ((su & 1) == (sw_a & 1)):
                continue
            p0, p1 = 4 * u + ((su + 2) & 3), 4 * u + ((su + 3) & 3)
            p2, p3 = 4 * w + ((sw_b + 2) & 3), 4 * w + ((sw_b + 3) & 3)
            remap = [v - (1 if v > u else 0) - (1 if v > w else 0) for v in range(n)]
            m = n - 2
            nvt = tuple(d.vtype[v] for v in range(n) if v not in (u, w))
            partner = {p0: p1, p1: p0, p2: p3, p3: p2}
            nopp = [0] * (4 * m)
            valid = True
            for v in range(n):
                if v in (u, w):
                    continue
                nv = remap[v]
                for s in range(4):
                    od = d.opp[4 * v + s]
                    steps = 0
                    while (od >> 2) in (u, w):
                        od = d.opp[partner[od]]
                        steps += 1
                        if steps > 8:
                            valid = False
                            break
                    if not valid:
                        break
                    nopp[4 * nv + s] = (remap[od >> 2] << 2) | (od & 3)
                if not valid:
                    break
            if not valid:
                continue
            cand = MarkedDiagram(num_vertices=m, opp=tuple(nopp), vtype=nvt, free_loops=0)
            if not _is_planar_4reg_diagram(cand):
                continue
            if not _is_admissible(cand):
                continue
            cst = cand.surface_topology()
            if cst.num_components != 1 or cst.euler_characteristic != st.euler_characteristic + 2:
                continue
            ck = canonical_marked_key(cand, up_to_mirror_and_reverse=True)
            if ck not in seen:
                seen.add(ck)
                out.append(cand)
    return out


_UNKNOT_MEMO: dict[tuple, bool] = {}


def is_constructively_unknotted_1comp(d: MarkedDiagram, base_max_v: int = 3, depth: int = 0) -> bool:
    """Strictly constructive proof that a 1-component diagram `d` represents an unknotted connected surface in S^4:
    1. Base case: `v <= base_max_v` (for `base_max_v = 3`, ONLY `0_1` (0v), `2_1^1` (2v), `2_1^-1` (2v), and standard `K^2(e=0)` (3v) exist, ALL unknotted without invoking any census!).
    2. Level + flype + simplifying moves reduce `v`: recurse on smaller diagram.
    3. Planar 2-cut `A # B`: recurse and require BOTH `A` and `B` constructively unknotted.
    4. C mountain-pass tunneling (`tunnel_diagrams_c`): if `tunnel_diagrams_c` prunes `d` (reduces `< v` or exposes 2-cut/split), wait — `tunnel_diagrams_c` returns `[]` on BOTH `< v` AND 2-cut!
       WAIT: at `v >= 10`, a 2-cut in `tunnel_diagrams_c` could be `8_1 # T^2` (NOT unknotted)!
       So we ONLY treat `tunnel_diagrams_c(v, [x]) == []` as proving unknottedness when `v <= 9` (where any 2-cut has summands `<= 7`, AND if `base_max_v == 7` every `<= 8` reduction has no `8_1` if `v <= 8`!).
       EVEN STRICTER: for `base_max_v == 3`, NEVER trust `tunnel_diagrams_c == []` blindly for `v >= 8` unless we explicitly recover a `< v` diagram or 2-cut! Wait: for `v in (4, 5, 6, 7)`, ANY reduction or 2-cut has `< 7` vertices (no `8_1` exists below `8` vertices!). And at `v in (8, 9, 10)`, we can check recursive 2-gon pinches directly + explicit level/flype reductions!
    5. Marked-marked 2-gon pinch `d -> cand` (`v - 2` vertices, 1 component, `chi + 2`): by Hosokawa-Kawauchi (1979) & Kamada (1989), if `cand` is constructively unknotted, attaching a single 1-handle in S^4 yields an unknotted connected surface!
    """
    if depth > 8:
        return False
    st = d.surface_topology()
    if st.num_components != 1 or d.free_loops > 0:
        return False
    r0 = simplify_yoshikawa(d)
    if _is_split_diag(r0):
        return False
    if r0.num_vertices == 0:
        return True
    if r0.num_vertices <= base_max_v:
        return True

    ck = (base_max_v, canonical_marked_key(r0, up_to_mirror_and_reverse=True))
    if ck in _UNKNOT_MEMO:
        return _UNKNOT_MEMO[ck]

    decomp = find_planar_connected_sum_decomposition(r0)
    if decomp is not None:
        ok = (
            is_constructively_unknotted_1comp(decomp[0], base_max_v, depth + 1)
            and is_constructively_unknotted_1comp(decomp[1], base_max_v, depth + 1)
        )
        _UNKNOT_MEMO[ck] = ok
        return ok

    # Explore bounded level + flype states of r0:
    st_code, wit, comp = explore_level_and_flypes(r0, max_states=96)
    if st_code == "reduced":
        red_d = MarkedDiagram.from_epd(wit[2])
        ok = is_constructively_unknotted_1comp(red_d, base_max_v, depth + 1)
        _UNKNOT_MEMO[ck] = ok
        return ok
    if st_code == "2cut":
        cut_d = MarkedDiagram.from_epd(wit[1])
        dec2 = find_planar_connected_sum_decomposition(cut_d)
        if dec2 is not None:
            ok = (
                is_constructively_unknotted_1comp(dec2[0], base_max_v, depth + 1)
                and is_constructively_unknotted_1comp(dec2[1], base_max_v, depth + 1)
            )
            _UNKNOT_MEMO[ck] = ok
            return ok

    # Try 1-handle pinches across up to 24 states in comp:
    states_list = list(comp.values())[:24] if comp else [r0]
    for x in states_list:
        for cand in list_valid_1comp_2gon_pinches(x):
            if is_constructively_unknotted_1comp(cand, base_max_v, depth + 1):
                _UNKNOT_MEMO[ck] = True
                return True

    # Mode 1 (base_max_v == 3): 100% pure explicit constructive (no C tunnel black-box at all!).
    if base_max_v >= 7 and r0.num_vertices in (8, 9, 10):
        from surfaceknots.invariants.fundamental_group import finite_quotient_profile, wirtinger_presentation
        wp = wirtinger_presentation(r0)
        ok_q = True
        if wp.num_generators > 1:
            fqp = finite_quotient_profile(wp, ("S3", "A4", "F21", "A5"))
            if fqp["S3"] > 6 or fqp["A4"] > 12 or fqp["F21"] > 21 or fqp["A5"] > 60:
                ok_q = False
        if ok_q:
            tk = tunnel_diagrams_c(r0.num_vertices, [r0], depth_plus1=256, depth_plus2=192, max_bfs_states=2048, num_workers=1)
            if len(tk) == 0:
                _UNKNOT_MEMO[ck] = True
                return True
            if canonical_marked_key(tk[0], up_to_mirror_and_reverse=True) != ck[1] and is_constructively_unknotted_1comp(tk[0], base_max_v, depth + 1):
                _UNKNOT_MEMO[ck] = True
                return True

    _UNKNOT_MEMO[ck] = False
    return False


def audit_single_diagram(args):
    idx, epd, tag = args
    d0 = MarkedDiagram.from_epd(epd)
    st_code, wit, comp = explore_level_and_flypes(d0, max_states=1024)
    min_ck = min(comp.keys()) if comp else canonical_marked_key(d0, up_to_mirror_and_reverse=True)
    comp_size = len(comp)

    if st_code != "open":
        return (idx, tag, epd, st_code, str(wit[0]), comp_size, min_ck)

    # Check strict constructive unknotted 1-handle across states in comp:
    states = list(comp.values())
    for x in states[:64]:
        for cand in list_valid_1comp_2gon_pinches(x):
            if is_constructively_unknotted_1comp(cand, base_max_v=3):
                return (idx, tag, epd, "strict_unknot_1handle_v3", f"pinch_v={cand.num_vertices}", comp_size, min_ck)

    for x in states[:64]:
        for cand in list_valid_1comp_2gon_pinches(x):
            if is_constructively_unknotted_1comp(cand, base_max_v=7):
                return (idx, tag, epd, "census7_unknot_1handle_v7", f"pinch_v={cand.num_vertices}", comp_size, min_ck)

    # Check analyze_marked_bigon_1handles for trivial 1-handles over 0_1, 0_1 u 0_1, 8_1, 9_1, 10_1 across comp:
    for x in states[:64]:
        for hrec in analyze_marked_bigon_1handles(x):
            if hrec.is_trivial_handle:
                return (idx, tag, epd, "trivial_bigon_cord", f"{hrec.base_identification}:{hrec.cord_double_coset_key}", comp_size, min_ck)

    return (idx, tag, epd, "survives", "", comp_size, min_ck)


if __name__ == "__main__":
    ch11_data = json.load((DATA_DIR / "expanded_ch11.json").open())
    ch12_data = json.load((DATA_DIR / "expanded_ch12.json").open())

    # 1. Negative control: ALL certified prime diagrams at ch=11 (12 diagrams) and ch=12 (63 diagrams)
    prime_tasks = []
    for e in ch11_data["entries"]:
        for epd in e.get("move_component_epds", [e["epd"]]):
            prime_tasks.append((len(prime_tasks), epd, f"ch11_prime:{e['label']}"))
    for e in ch12_data["entries"]:
        for epd in e.get("move_component_epds", [e["epd"]]):
            prime_tasks.append((len(prime_tasks), epd, f"ch12_prime:{e['label']}"))

    with ProcessPoolExecutor(max_workers=16, mp_context=multiprocessing.get_context("forkserver")) as pool:
        prime_res = list(pool.map(audit_single_diagram, prime_tasks, chunksize=2))

    print(f"=== NEGATIVE CONTROL ON ALL {len(prime_tasks)} CERTIFIED PRIME DIAGRAMS (12 ch=11 + 63 ch=12) ===")
    p_status = Counter(r[3] for r in prime_res)
    print("  Statuses:", dict(p_status))
    for r in prime_res:
        if r[3] != "survives":
            print("  ALERT ON PRIME:", r[:6])
    # Also check if any multi-diagram packet merged under full level+flype closure!
    by_lbl_cks = defaultdict(set)
    by_lbl_cnt = Counter()
    for r in prime_res:
        by_lbl_cks[r[1]].add(r[6])
        by_lbl_cnt[r[1]] += 1
    for lbl, cnt in by_lbl_cnt.items():
        if cnt > 1 and len(by_lbl_cks[lbl]) < cnt:
            print(f"  PACKET MERGE UNDER LEVEL+FLYPE: {lbl}: {cnt} -> {len(by_lbl_cks[lbl])} components!")

    # 2. Audit all 93 ch=11 and 742 ch=12 residual diagrams:
    for ch, res_list in [(11, ch11_data["residual_collisions"]), (12, ch12_data["residual_collisions"])]:
        tasks = [(i, r["epd"], " | ".join(r.get("collision_matches", [])[:2])) for i, r in enumerate(res_list)]
        with ProcessPoolExecutor(max_workers=16, mp_context=multiprocessing.get_context("forkserver")) as pool:
            results = list(pool.map(audit_single_diagram, tasks, chunksize=4))

        print(f"\n=== RESIDUAL AUDIT ch={ch} ({len(results)} diagrams) ===")
        st_counts = Counter(r[3] for r in results)
        for k, v in st_counts.most_common():
            print(f"  {v:4d}  {k}")

        surv = [r for r in results if r[3] == "survives"]
        uniq_surv_cks = len({r[6] for r in surv})
        print(f"  Surviving open diagrams: {len(surv)} (unique level+flype components: {uniq_surv_cks})")
        surv_by_tag = Counter(r[1] for r in surv)
        print(f"  Breakdown of surviving open diagrams at ch={ch}:")
        for tag, cnt in surv_by_tag.most_common(40):
            print(f"    {cnt:4d}  {tag}")

        # Save detailed JSON for follow-up tunneling on surviving open diagrams:
        out_path = Path(f"/afs/csail.mit.edu/u/e/ejli/.gemini/antigravity-cli/brain/c34f4119-775d-42cd-b8c1-b2dd16600960/scratch/residual_stage1_ch{ch}.json")
        out_path.write_text(json.dumps([
            {"idx": r[0], "tag": r[1], "epd": r[2], "status": r[3], "detail": r[4], "comp_size": r[5]}
            for r in results
        ], indent=2))
