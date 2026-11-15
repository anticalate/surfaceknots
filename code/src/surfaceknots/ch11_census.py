"""Stratification, composite/split reference catalog, 4D flype merge, and dataset export for the
expanded `ch = 11` surface-link census.
"""

from __future__ import annotations

from collections import defaultdict
import csv
import json
from pathlib import Path

from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.moves import (
    canonical_marked_key,
    find_planar_connected_sum_decomposition,
    has_certified_cyclic_unknotted_tube,
    list_4d_tait_flypes,
    simplify_and_explore_yoshikawa,
)
from surfaceknots.unlink import certify_classical_unlink
from surfaceknots.yoshikawa_data import (
    YOSHIKAWA_BY_LABEL,
    YOSHIKAWA_TABLE,
    assign_nearest_separating_witnesses,
    compute_full_invariant_record,
    invariant_discrimination_key,
)


def planar_connected_sum(
    d1: MarkedDiagram,
    dart1: int,
    d2: MarkedDiagram,
    dart2: int,
    flip_splice: bool = False,
) -> MarkedDiagram:
    """Form the planar connected sum `d1 # d2` by splicing edge `(dart1, opp[dart1])` with `(dart2, opp[dart2])`."""
    n1, n2 = d1.num_vertices, d2.num_vertices
    assert n1 > 0 and n2 > 0
    b1 = d1.opp[dart1]
    b2 = d2.opp[dart2]
    off = 4 * n1

    new_opp = list(d1.opp) + [x + off for x in d2.opp]
    if not flip_splice:
        new_opp[dart1] = b2 + off
        new_opp[b2 + off] = dart1
        new_opp[b1] = dart2 + off
        new_opp[dart2 + off] = b1
    else:
        new_opp[dart1] = dart2 + off
        new_opp[dart2 + off] = dart1
        new_opp[b1] = b2 + off
        new_opp[b2 + off] = b1

    res = MarkedDiagram(
        num_vertices=n1 + n2,
        opp=tuple(new_opp),
        vtype=d1.vtype + d2.vtype,
        free_loops=d1.free_loops + d2.free_loops,
    )
    return res


def disjoint_union(d1: MarkedDiagram, d2: MarkedDiagram) -> MarkedDiagram:
    """Form the split union `d1 u d2` of two `MarkedDiagram`s."""
    n1, n2 = d1.num_vertices, d2.num_vertices
    off = 4 * n1
    new_opp = list(d1.opp) + [x + off for x in d2.opp]
    return MarkedDiagram(
        num_vertices=n1 + n2,
        opp=tuple(new_opp),
        vtype=d1.vtype + d2.vtype,
        free_loops=d1.free_loops + d2.free_loops,
    )


def _is_admissible(d: MarkedDiagram) -> bool:
    return (
        certify_classical_unlink(d.resolve(1)).is_unlink
        and certify_classical_unlink(d.resolve(-1)).is_unlink
    )


def _compute_record_worker(args: tuple[MarkedDiagram, str, bool]) -> dict[str, object]:
    d, tag, inc_sub = args
    return compute_full_invariant_record(d, label=tag, include_sublinks=inc_sub)


def _oriented_diagram_variants(lbl: str, d: MarkedDiagram) -> list[tuple[str, MarkedDiagram]]:
    """Return up to 4 distinct oriented-diagram representatives (`id`, `mirror_3d`, `mirror_time`,
    `reverse_orientation_zt`) of `d`, deduplicated up to planar oriented diagram isomorphism.
    """
    if d.num_vertices == 0:
        return [(lbl, d)]
    seen: set[tuple] = set()
    out: list[tuple[str, MarkedDiagram]] = []
    for suffix, cand in (
        ("", d),
        ("*", d.mirror_3d()),
        ("*t", d.mirror_time()),
        ("*zt", d.reverse_orientation_zt()),
    ):
        ck = canonical_marked_key(cand, up_to_mirror_and_reverse=False)
        if ck not in seen:
            seen.add(ck)
            out.append((f"{lbl}{suffix}", cand))
    return out


def build_reference_catalog(max_composite_ch: int = 12, num_workers: int = 16) -> dict[tuple, list[str]]:
    """Build an exhaustive invariant lookup table of:
    1. All 23 `ch <= 10` Yoshikawa prime surface-links (`ch<=10:{label}`),
    2. Minimal-`ch` unknotted connected surfaces up to `ch <= 12` (`unknot:{type}(e={e})`,
       including the 3-vertex `K^2(e=0)` diagram),
    3. Split unions of unknots and `ch <= 10` entries up to `ch <= max(11, max_composite_ch)` (`split:{A} u {B}`),
    4. Nontrivial connected sums `A # B` across all 4 orientation/mirror variants and both splice
       orientations up to `ch <= max_composite_ch` (`composite:{A} # {B}`).
    """
    from surfaceknots.complexity import construct_minimal_unknot_diagram

    catalog: dict[tuple, list[str]] = defaultdict(list)
    queued: list[tuple[MarkedDiagram, str, bool]] = []
    tags_by_ck: dict[tuple, list[str]] = defaultdict(list)

    def register(d: MarkedDiagram, tag: str) -> None:
        if not _is_admissible(d):
            return
        ck = canonical_marked_key(d, up_to_mirror_and_reverse=True)
        if ck not in tags_by_ck:
            queued.append((d, tag, False))
        if tag not in tags_by_ck[ck]:
            tags_by_ck[ck].append(tag)

    # 1. All 23 Yoshikawa entries (`ch <= 10`):
    yosh_diags: list[tuple[str, MarkedDiagram]] = []
    for entry in YOSHIKAWA_TABLE:
        d = MarkedDiagram.from_epd(entry.canonical_epd)
        register(d, f"ch<=10:{entry.label}")
        yosh_diags.append((entry.label, d))

    # 2. Build all minimal-`ch` connected unknotted surfaces up to ch = 12 (including 3-vertex K^2(e=0)):
    unknots_by_topo: dict[tuple[int, int], tuple[str, MarkedDiagram]] = {
        (0, 0): ("0_1", MarkedDiagram.unknot_sphere()),
    }
    for g in range(1, 6):
        ud = construct_minimal_unknot_diagram(g, 0)
        ct = ud.surface_topology().components[0]
        lbl = "2_1^1" if g == 1 else f"unknot({ct.homeomorphism_type},e=0)"
        unknots_by_topo[(g, 0)] = (lbl, ud)
        if g >= 2:
            register(ud, f"unknot:{ct.homeomorphism_type}(e=0)")
    for k in range(1, 11):
        for e in range(-2 * k, 2 * k + 1):
            if (e - 2 * k) % 4 != 0:
                continue
            if max(1, abs(e) // 2) + k > 11:
                continue
            ud = construct_minimal_unknot_diagram(-k, e)
            ct = ud.surface_topology().components[0]
            if k == 1:
                lbl = "P+" if e == 2 else "P-"
            else:
                lbl = f"unknot({ct.homeomorphism_type},e={e})"
                register(ud, f"unknot:{ct.homeomorphism_type}(e={e})")
            unknots_by_topo[(-k, e)] = (lbl, ud)

    # 3. Nontrivial connected sums `A # U` and `A # B` with `ch(A) + ch(B) <= max_composite_ch`:
    summands: list[tuple[str, MarkedDiagram]] = []
    for entry in YOSHIKAWA_TABLE:
        if entry.ch_index == 0 or entry.ch_index > max_composite_ch - 2:
            continue
        d = MarkedDiagram.from_epd(entry.canonical_epd)
        summands.append((entry.label, d))
    for (gs, es), (ulbl, ud) in unknots_by_topo.items():
        if 3 <= ud.ch_index <= 6:
            summands.append((ulbl, ud))

    def _is_pure_unknot_label(lbl: str) -> bool:
        base = lbl.split("*")[0]
        return base in ("0_1", "2_1^1", "2_1^-1", "P+", "P-") or base.startswith("unknot")

    for i, (lbl_a, da) in enumerate(summands):
        if da.num_vertices == 0:
            continue
        dc_a, nc_a = da._dart_surface_components()
        rep_dart_a = {}
        for hd in range(4 * da.num_vertices):
            rep_dart_a.setdefault(dc_a[hd], hd)

        for j in range(i, len(summands)):
            lbl_b, db = summands[j]
            if db.num_vertices == 0 or da.num_vertices + db.num_vertices > max_composite_ch:
                continue
            # Skip unknot # unknot (already registered as `unknot:*`!)
            if _is_pure_unknot_label(lbl_a) and _is_pure_unknot_label(lbl_b):
                continue
            for lbl_b_var, db_var in _oriented_diagram_variants(lbl_b, db):
                dc_b, _ = db_var._dart_surface_components()
                rep_dart_b = {}
                for hd in range(4 * db_var.num_vertices):
                    rep_dart_b.setdefault(dc_b[hd], hd)

                for ca, hda in rep_dart_a.items():
                    for cb, hdb in rep_dart_b.items():
                        for fl in (False, True):
                            csum = planar_connected_sum(da, hda, db_var, hdb, flip_splice=fl)
                            if _is_admissible(csum):
                                register(csum, f"composite:{lbl_a}#{lbl_b_var}(c{ca}+c{cb})")

    # 4. Split unions `A u B` (and `A u 0_1 u 0_1`, `unknot u unknot`) up to ch <= max(11, max_composite_ch):
    max_split_ch = max(11, max_composite_ch)
    split_blocks: list[tuple[str, MarkedDiagram]] = [("0_1", MarkedDiagram.unknot_sphere())]
    for (gs, es), (ulbl, ud) in unknots_by_topo.items():
        if 2 <= ud.ch_index <= max(8, max_split_ch - 2):
            split_blocks.append((ulbl, ud))
    for entry in YOSHIKAWA_TABLE:
        if 6 <= entry.ch_index <= 10:
            split_blocks.append((entry.label, MarkedDiagram.from_epd(entry.canonical_epd)))

    for i, (la, da) in enumerate(split_blocks):
        for j in range(i, len(split_blocks)):
            lb, db = split_blocks[j]
            if da.ch_index + db.ch_index > max_split_ch:
                continue
            sp2 = disjoint_union(da, db)
            register(sp2, f"split:{la}_u_{lb}")
            # Also 3-component splits with 0_1:
            if (
                sp2.ch_index <= max(9, max_split_ch - 2)
                and _is_admissible(sp2)
                and sp2.surface_topology().num_components == 2
            ):
                register(disjoint_union(sp2, MarkedDiagram.unknot_sphere()), f"split:{la}_u_{lb}_u_0_1")

    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing

    with ProcessPoolExecutor(
        max_workers=max(1, num_workers), mp_context=multiprocessing.get_context("forkserver")
    ) as pool:
        for (d_rep, _, _), rec in zip(queued, pool.map(_compute_record_worker, queued, chunksize=8)):
            k = invariant_discrimination_key(rec)
            ck_rep = canonical_marked_key(d_rep, up_to_mirror_and_reverse=True)
            for tag in tags_by_ck[ck_rep]:
                if tag not in catalog[k]:
                    catalog[k].append(tag)

    return dict(catalog)


def post_filter_and_merge_ch11_survivors(
    raw_epds: list[str],
) -> tuple[list[MarkedDiagram], dict[str, int]]:
    """Apply 4D Tait flype simplification, post-flype C mountain-pass tunneling, composite-cut check,
    and cyclic-tube pruning on `ch = 11` survivors.
    """
    from surfaceknots.fast_enum import tunnel_diagrams_c

    kept_by_key: dict[tuple, MarkedDiagram] = {}
    pruned_flype_reduce = 0
    pruned_flype_cut_or_split = 0
    pruned_cyclic_tube = 0

    for epd in raw_epds:
        d0 = MarkedDiagram.from_epd(epd)
        best, is_sp, is_cp = simplify_and_explore_yoshikawa(d0, max_states=64)
        if is_sp or is_cp:
            pruned_flype_cut_or_split += 1
            continue
        if best.num_vertices < 11:
            pruned_flype_reduce += 1
            continue

        # Probe 4D Tait flypes from `best`:
        reduced_via_flype = False
        for fl_d in list_4d_tait_flypes(best):
            if not _is_admissible(fl_d):
                continue
            b2, sp2, cp2 = simplify_and_explore_yoshikawa(fl_d, max_states=64)
            if sp2 or cp2:
                pruned_flype_cut_or_split += 1
                reduced_via_flype = True
                break
            if b2.num_vertices < 11:
                pruned_flype_reduce += 1
                reduced_via_flype = True
                break
            k2 = canonical_marked_key(b2, up_to_mirror_and_reverse=True)
            k_cur = canonical_marked_key(best, up_to_mirror_and_reverse=True)
            if k2 < k_cur:
                best = b2
        if reduced_via_flype:
            continue

        if find_planar_connected_sum_decomposition(best) is not None or has_certified_cyclic_unknotted_tube(best):
            pruned_cyclic_tube += 1
            continue

        ckey = canonical_marked_key(best, up_to_mirror_and_reverse=True)
        if ckey not in kept_by_key:
            kept_by_key[ckey] = best

    pre_tunnel_survivors = [kept_by_key[k] for k in sorted(kept_by_key)]
    # Post-flype C mountain-pass reduction (4D flypes unlock additional <= 13-vertex mountain-pass collapses):
    tunneled = tunnel_diagrams_c(
        11,
        pre_tunnel_survivors,
        depth_plus1=384,
        depth_plus2=256,
        max_bfs_states=4096,
    )
    final_by_key: dict[tuple, MarkedDiagram] = {}
    for td in tunneled:
        if find_planar_connected_sum_decomposition(td) is not None or has_certified_cyclic_unknotted_tube(td):
            pruned_cyclic_tube += 1
            continue
        tk = canonical_marked_key(td, up_to_mirror_and_reverse=True)
        if tk not in final_by_key:
            final_by_key[tk] = td

    survivors = [final_by_key[k] for k in sorted(final_by_key)]
    pruned_post_flype_c_tunnel = len(pre_tunnel_survivors) - len(tunneled)
    return survivors, {
        "input_survivors": len(raw_epds),
        "pruned_flype_reduce": pruned_flype_reduce,
        "pruned_flype_cut_or_split": pruned_flype_cut_or_split,
        "pruned_cyclic_tube": pruned_cyclic_tube,
        "pre_post_flype_tunnel_diagrams": len(pre_tunnel_survivors),
        "pruned_post_flype_c_tunnel": pruned_post_flype_c_tunnel,
        "post_flype_unique_diagrams": len(survivors),
    }


def enumerate_ch11_checkpointed(
    out_dir: Path,
    num_workers: int = 32,
    depth_plus1: int = 256,
    depth_plus2: int = 256,
    max_bfs_states: int = 2048,
    force_recompute: bool = False,
) -> tuple[list[str], dict[str, object]]:
    """Enumerate all `ch = 11` prime-shadow diagrams across `m = 1..10` (`m = 0` and `m = 11` proved
    unknotted) with clean per-`m` telemetry saved to `out_dir / 'ch11_raw_survivors.json'`.
    """
    import time
    from surfaceknots.fast_enum import (
        _filter_cyclic_tubes_parallel,
        cascade_tunnel_diagrams_c,
        enumerate_targeted_census_level_candidates,
        merge_diagrams_c,
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    ckpt_path = out_dir / "ch11_raw_survivors.json"
    if ckpt_path.exists() and not force_recompute:
        ckpt = json.loads(ckpt_path.read_text())
        if ckpt.get("status") == "COMPLETE" and "survivor_epds" in ckpt:
            return list(ckpt["survivor_epds"]), dict(ckpt["census_stats"])

    t0 = time.perf_counter()
    per_m_stats: dict[str, dict[str, object]] = {}
    per_m_diags: list[MarkedDiagram] = []
    for m_val in range(1, 11):
        tm0 = time.perf_counter()
        diags_m, st_m = enumerate_targeted_census_level_candidates(
            ch=11,
            m_min=m_val,
            m_max=m_val,
            max_bfs_states=max_bfs_states,
            num_workers=num_workers,
            _certify_unlink=False,
        )
        tun_m = cascade_tunnel_diagrams_c(
            ch=11,
            diagrams=diags_m,
            depth_plus1=depth_plus1,
            depth_plus2=depth_plus2,
            max_bfs_states=max_bfs_states,
            num_workers=num_workers,
        )
        dtm = round(time.perf_counter() - tm0, 3)
        st_rec = dict(st_m)
        st_rec["intra_m_tunneled_kept"] = len(tun_m)
        st_rec["runtime_sec"] = dtm
        per_m_stats[str(m_val)] = st_rec
        per_m_diags.extend(tun_m)

    merged_all = merge_diagrams_c(
        11,
        per_m_diags,
        max_bfs_states=max(4096, max_bfs_states * 2),
        num_workers=num_workers,
        _certify_unlink=False,
    )
    tunneled_all = cascade_tunnel_diagrams_c(
        ch=11,
        diagrams=merged_all,
        depth_plus1=depth_plus1,
        depth_plus2=depth_plus2,
        max_bfs_states=max_bfs_states,
        num_workers=num_workers,
    )
    filtered_all = _filter_cyclic_tubes_parallel(tunneled_all, num_workers=num_workers)

    m2_7 = range(2, 8)
    m_all = range(1, 11)
    census_stats: dict[str, object] = {
        "num_shadows": 3378,
        "marker_range_scanned": [1, 10],
        "marker_extremes_note": "m=0 (0 saddles => unlink of S^2) and m=11 (0 crossings => planar surface in R^3) are unknotted",
        "num_workers": num_workers,
        "raw_cands": sum(int(per_m_stats[str(m)]["raw_cands"]) for m in m_all),
        "raw_cands_m2_to_7": sum(int(per_m_stats[str(m)]["raw_cands"]) for m in m2_7),
        "components": sum(int(per_m_stats[str(m)]["components"]) for m in m_all),
        "components_m2_to_7": sum(int(per_m_stats[str(m)]["components"]) for m in m2_7),
        "discarded": sum(int(per_m_stats[str(m)]["discarded"]) for m in m_all),
        "level_kept": sum(int(per_m_stats[str(m)]["kept"]) for m in m_all),
        "level_kept_m2_to_7": sum(int(per_m_stats[str(m)]["kept"]) for m in m2_7),
        "tunneled_kept": len(tunneled_all),
        "cyclic_tube_pruned": len(tunneled_all) - len(filtered_all),
        "cyclic_tube_filtered_kept": len(filtered_all),
        "total_enumeration_wall_sec": round(time.perf_counter() - t0, 3),
        "per_m_stats": per_m_stats,
    }
    survivor_epds = [d.to_epd() for d in filtered_all]
    ckpt_path.write_text(
        json.dumps(
            {
                "status": "COMPLETE",
                "ch_index": 11,
                "census_stats": census_stats,
                "survivor_epds": survivor_epds,
            },
            indent=2,
        )
        + "\n"
    )
    return survivor_epds, census_stats


def stratify_and_export_ch11(
    survivor_epds: list[str],
    out_dir: Path,
    census_stats: dict[str, object] | None = None,
    num_workers: int = 16,
) -> dict[str, object]:
    """Stratify `ch = 11` survivors against the exhaustive `ch <= 10` / composite / split / unknot catalog
    (including Yoshikawa Def. 5.1 weak-primality verification up to 20 vertices) and export
    `expanded_ch11.json` and `expanded_ch11.csv`.
    """
    from concurrent.futures import ProcessPoolExecutor
    from dataclasses import asdict
    from surfaceknots.ch12_census import (
        build_def5_1_high_composite_matches,
        compute_indecomposability_certificate,
    )
    from surfaceknots.move_search import analyze_marked_bigon_1handles
    from surfaceknots.yoshikawa_data import export_yoshikawa_datasets

    out_dir.mkdir(parents=True, exist_ok=True)
    yosh_json = out_dir / "yoshikawa_ch10.json"
    if yosh_json.exists():
        yosh_records: list[dict[str, object]] = json.loads(yosh_json.read_text())
        if not yosh_records or "sublink_invariant_profile" not in yosh_records[0]:
            yosh_records = export_yoshikawa_datasets(out_dir)
    else:
        yosh_records = export_yoshikawa_datasets(out_dir)

    survivors, post_stats = post_filter_and_merge_ch11_survivors(survivor_epds)
    ref_catalog = build_reference_catalog(max_composite_ch=13, num_workers=num_workers)

    # Compute full invariants on every surviving `ch = 11` diagram in parallel:
    grouped_new: dict[tuple, list[dict[str, object]]] = defaultdict(list)
    grouped_residual: dict[tuple, list[dict[str, object]]] = defaultdict(list)
    residual_collisions: list[dict[str, object]] = []

    surv_tasks = [(d, f"ch11_cand_{idx}", True) for idx, d in enumerate(survivors)]
    import multiprocessing

    with ProcessPoolExecutor(
        max_workers=max(1, num_workers), mp_context=multiprocessing.get_context("forkserver")
    ) as pool:
        for d_obj, rec in zip(survivors, pool.map(_compute_record_worker, surv_tasks, chunksize=4)):
            rec["indecomposability"] = compute_indecomposability_certificate(d_obj, rec, 20)
            dkey = invariant_discrimination_key(rec)
            if dkey in ref_catalog:
                rec["cell_status"] = "AMBIGUOUS_SHADOW"
                rec["collision_matches"] = ref_catalog[dkey]
                rec["marked_bigon_1handles"] = [asdict(h) for h in analyze_marked_bigon_1handles(d_obj)]
                residual_collisions.append(rec)
                grouped_residual[dkey].append(rec)
            else:
                grouped_new[dkey].append(rec)

    # Verify Yoshikawa Definition 5.1 weak primality against all pairs (F1, F2) with ch(F1) <= 10, ch(F2) <= 10:
    prelim_reps = [grp[0] for grp in grouped_new.values()]
    def5_1_extra = build_def5_1_high_composite_matches(11, prelim_reps, out_dir, num_workers=num_workers)
    for dkey, tags in def5_1_extra.items():
        if dkey in grouped_new:
            moved = grouped_new.pop(dkey)
            for rec in moved:
                rec["cell_status"] = "AMBIGUOUS_SHADOW"
                rec["collision_matches"] = tags
                residual_collisions.append(rec)
                grouped_residual[dkey].append(rec)

    # Sort certified new `ch = 11` invariant classes in Yoshikawa order:
    ori_rank = {"Ori": 0, "Mixed": 1, "Non": 2}
    sorted_groups = sorted(
        grouped_new.values(),
        key=lambda grp: (
            grp[0]["num_components"],
            ori_rank.get(str(grp[0]["orientability"]), 9),
            tuple(grp[0]["yoshikawa_genera"]),  # type: ignore[arg-type]
            tuple(grp[0]["normal_euler_numbers"]),  # type: ignore[arg-type]
            tuple(grp[0]["e1_generators"]),  # type: ignore[arg-type]
            tuple(sorted(grp[0]["hom_counts"].items())),  # type: ignore[union-attr]
        ),
    )

    # Assign Yoshikawa-style labels `11_k^{g_1,...,g_c}` (`11_k` for connected `S^2` knots):
    genera_counters: dict[tuple[int, ...], int] = defaultdict(int)
    certified_new_entries: list[dict[str, object]] = []

    for grp in sorted_groups:
        grp_sorted = sorted(
            grp,
            key=lambda r: (int(r["pi1_relators"]), int(r["pi1_generators"]), str(r["epd"])),  # type: ignore[arg-type]
        )
        rep = dict(grp_sorted[0])
        gen_tup = tuple(int(x) for x in rep["yoshikawa_genera"])  # type: ignore[union-attr]
        genera_counters[gen_tup] += 1
        k_idx = genera_counters[gen_tup]
        if gen_tup == (0,):
            y_label = f"11_{k_idx}"
        else:
            sup = ",".join(str(g) for g in gen_tup)
            y_label = f"11_{k_idx}^{sup}"
        m_comps = len(grp_sorted)
        rep["label"] = y_label
        rep["cell_status"] = "EXACT_SINGLETON" if m_comps == 1 else "NEW_PACKET"
        rep["num_move_components"] = m_comps
        rep["move_component_ids"] = (
            [y_label] if m_comps == 1 else [f"{y_label}.{chr(ord('a') + i)}" for i in range(m_comps)]
        )
        rep["move_component_epds"] = [str(r["epd"]) for r in grp_sorted]
        rep["num_diagram_representatives"] = m_comps
        rep["alternate_epds"] = [str(r["epd"]) for r in grp_sorted[1:]]
        certified_new_entries.append(rep)

    # Compute nearest separating witnesses across the full 32-entry (`ch <= 10` + `ch = 11`) prime catalog:
    assign_nearest_separating_witnesses(
        certified_new_entries,
        reference_pool=yosh_records + certified_new_entries,
    )

    exact_singletons_count = sum(1 for r in certified_new_entries if r["cell_status"] == "EXACT_SINGLETON")
    new_packets_count = sum(1 for r in certified_new_entries if r["cell_status"] == "NEW_PACKET")

    summary_payload = {
        "ch_index": 11,
        "census_stats": census_stats or {},
        "post_filter_stats": post_stats,
        "certified_new_prime_invariant_classes": len(certified_new_entries),
        "exact_singletons_count": exact_singletons_count,
        "new_packets_count": new_packets_count,
        "ambiguous_shadow_cells_count": len(grouped_residual),
        "certified_new_prime_diagrams_total": sum(
            int(r["num_diagram_representatives"]) for r in certified_new_entries  # type: ignore[arg-type]
        ),
        "residual_colliding_diagrams": len(residual_collisions),
        "upper_bound_total_diagram_classes": len(certified_new_entries) + len(residual_collisions),
        "entries": certified_new_entries,
        "residual_collisions": residual_collisions,
    }

    json_path = out_dir / "expanded_ch11.json"
    json_path.write_text(json.dumps(summary_payload, indent=2) + "\n")

    csv_path = out_dir / "expanded_ch11.csv"
    csv_fields = [
        "label",
        "cell_status",
        "num_move_components",
        "ch_index_diagram",
        "num_diagram_representatives",
        "num_components",
        "component_types",
        "orientability",
        "euler_characteristic",
        "normal_euler_numbers",
        "pi1_presentation",
        "e1_generators",
        "sigma2_b1",
        "sigma2_h1_torsion",
        "sigma2_b2",
        "sigma2_abs_signature",
        "certified_pi2_z_rank",
        "hom_S3",
        "hom_A4",
        "hom_S4",
        "hom_A5",
        "hom_Q8",
        "hom_Dic16",
        "sl2f3_hom_count",
        "sl2f3_non_meridional_peripheral_count",
        "double_linking_mod2",
        "mochizuki_c3_mod3",
        "sym_quandle_q6_col",
        "sym_quandle_q6_nontrivial_weights",
        "symmetry_summary",
        "nearest_separating_label",
        "nearest_separating_layer",
        "epd",
    ]
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_fields)
        writer.writeheader()
        for r in certified_new_entries:
            hc = r["hom_counts"]  # type: ignore[assignment]
            wit = r["nearest_separating_witness"]  # type: ignore[assignment]
            sym_d = r["symmetry"]  # type: ignore[assignment]
            writer.writerow(
                {
                    "label": r["label"],
                    "cell_status": r["cell_status"],
                    "num_move_components": r["num_move_components"],
                    "ch_index_diagram": r["ch_index_diagram"],
                    "num_diagram_representatives": r["num_diagram_representatives"],
                    "num_components": r["num_components"],
                    "component_types": " u ".join(r["component_types"]),  # type: ignore[arg-type]
                    "orientability": r["orientability"],
                    "euler_characteristic": r["euler_characteristic"],
                    "normal_euler_numbers": json.dumps(r["normal_euler_numbers"]),
                    "pi1_presentation": r["pi1_presentation"],
                    "e1_generators": json.dumps(r["e1_generators"]),
                    "sigma2_b1": r["sigma2_b1"],
                    "sigma2_h1_torsion": json.dumps(r["sigma2_h1_torsion"]),
                    "sigma2_b2": r["sigma2_b2"],
                    "sigma2_abs_signature": r["sigma2_abs_signature"],
                    "certified_pi2_z_rank": r["certified_pi2_z_rank"],
                    "hom_S3": hc["S3"],  # type: ignore[index]
                    "hom_A4": hc["A4"],  # type: ignore[index]
                    "hom_S4": hc["S4"],  # type: ignore[index]
                    "hom_A5": hc["A5"],  # type: ignore[index]
                    "hom_Q8": hc["Q8"],  # type: ignore[index]
                    "hom_Dic16": hc["Dic16"],  # type: ignore[index]
                    "sl2f3_hom_count": r["sl2f3_hom_count"],
                    "sl2f3_non_meridional_peripheral_count": r["sl2f3_non_meridional_peripheral_count"],
                    "double_linking_mod2": json.dumps(r["double_linking_mod2"]),
                    "mochizuki_c3_mod3": json.dumps(r["mochizuki_c3_mod3"]),
                    "sym_quandle_q6_col": r["sym_quandle_q6_col"],
                    "sym_quandle_q6_nontrivial_weights": r["sym_quandle_q6_nontrivial_weights"],
                    "symmetry_summary": sym_d["summary"],  # type: ignore[index]
                    "nearest_separating_label": wit["nearest_label"],  # type: ignore[index]
                    "nearest_separating_layer": wit["separating_layer"],  # type: ignore[index]
                    "epd": r["epd"],
                }
            )

    return summary_payload
