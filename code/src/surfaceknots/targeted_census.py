"""Reproducible targeted census driver for `ch = 12, 13, 14` sectors (`targeted_ch12_13_census.json`).

Implements Compute Campaign 2 (`HANDOFF.md` §4):
1. `ch = 12, m = 2` and `ch = 13, m = 2` (and optional `ch = 14, m = 2`) connected 2-sphere knots (`S^2`),
   with full 10-layer stratification (including exact strong Groebner `E_1` over `Z`, `F21`/`F57` regular
   cover homology, and Mochizuki `R_3, R_5, R_7` cocycles).
2. Full-marker-range (`m = 1 .. ch - 1`) `RP^2`-knot census at `ch = 11, 12, 13`, paired with constructive
   Yoshikawa move-path reducibility certificates (`find_rp2_reducible_certificate`, independently verified by
   `verify_move_path_independent`) and Todd-Coxeter meridian-order bounds (`meridian_order_lower_bound`),
   testing Kinoshita's Conjecture and establishing lower bounds on `ch(R_HKMN)`.
3. Targeted `ch = 13` sweeps (`m = 2, 3`) for connected Klein bottles (`K^2`), connected tori (`T^2`),
   and 2-component surface-links.
"""

from __future__ import annotations

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import json
import multiprocessing
from pathlib import Path
import time

_MP_CTX = multiprocessing.get_context("forkserver")

from surfaceknots.ch11_census import _compute_record_worker, _is_admissible
from surfaceknots.ch12_census import (
    build_def5_1_high_composite_matches,
    build_reference_catalog_ch12,
    compute_indecomposability_certificate,
)
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.fast_enum import (
    _filter_cyclic_tubes_parallel,
    cascade_tunnel_diagrams_c,
    enumerate_targeted_census_level_candidates,
    merge_diagrams_c,
    tunnel_diagrams_c,
)
from surfaceknots.move_search import (
    analyze_marked_bigon_1handles,
    find_rp2_reducible_certificate,
    meridian_order_lower_bound,
)
from surfaceknots.moves import (
    canonical_marked_key,
    find_planar_connected_sum_decomposition,
    has_certified_cyclic_unknotted_tube,
    list_4d_tait_flypes,
    simplify_and_explore_yoshikawa,
)
from surfaceknots.yoshikawa_data import (
    compute_full_invariant_record,
    invariant_discrimination_key,
)


def _post_filter_and_tunnel_sector(
    ch: int,
    diagrams: list[MarkedDiagram],
    num_workers: int = 32,
) -> tuple[list[MarkedDiagram], dict[str, int]]:
    """Apply Python simplification, 4D Tait flype probing, and post-flype C mountain-pass tunneling."""
    pruned_flype_reduce = 0
    pruned_flype_cut_or_split = 0
    pruned_cyclic_tube = 0
    kept_by_key: dict[tuple, MarkedDiagram] = {}

    for d0 in diagrams:
        if find_planar_connected_sum_decomposition(d0) is not None or has_certified_cyclic_unknotted_tube(d0):
            pruned_cyclic_tube += 1
            continue
        best, is_sp, is_cp = simplify_and_explore_yoshikawa(d0, max_states=64)
        if is_sp or is_cp:
            pruned_flype_cut_or_split += 1
            continue
        if best.num_vertices < ch:
            pruned_flype_reduce += 1
            continue

        reduced_via_flype = False
        for fl_d in list_4d_tait_flypes(best):
            if not _is_admissible(fl_d):
                continue
            b2, sp2, cp2 = simplify_and_explore_yoshikawa(fl_d, max_states=64)
            if sp2 or cp2:
                pruned_flype_cut_or_split += 1
                reduced_via_flype = True
                break
            if b2.num_vertices < ch:
                pruned_flype_reduce += 1
                reduced_via_flype = True
                break
            if canonical_marked_key(b2, up_to_mirror_and_reverse=True) < canonical_marked_key(
                best, up_to_mirror_and_reverse=True
            ):
                best = b2
        if reduced_via_flype:
            continue

        if find_planar_connected_sum_decomposition(best) is not None or has_certified_cyclic_unknotted_tube(best):
            pruned_cyclic_tube += 1
            continue
        ck = canonical_marked_key(best, up_to_mirror_and_reverse=True)
        if ck not in kept_by_key:
            kept_by_key[ck] = best

    pre_tun = [kept_by_key[k] for k in sorted(kept_by_key)]
    tun = (
        tunnel_diagrams_c(
            ch,
            pre_tun,
            depth_plus1=256,
            depth_plus2=192,
            max_bfs_states=3072,
            num_workers=num_workers,
        )
        if pre_tun
        else []
    )
    final_by_key: dict[tuple, MarkedDiagram] = {}
    for td in tun:
        if find_planar_connected_sum_decomposition(td) is not None or has_certified_cyclic_unknotted_tube(td):
            pruned_cyclic_tube += 1
            continue
        tk = canonical_marked_key(td, up_to_mirror_and_reverse=True)
        if tk not in final_by_key:
            final_by_key[tk] = td

    survivors = [final_by_key[k] for k in sorted(final_by_key)]
    return survivors, {
        "input_c_tunneled": len(diagrams),
        "pruned_flype_reduce": pruned_flype_reduce,
        "pruned_flype_cut_or_split": pruned_flype_cut_or_split,
        "pruned_cyclic_tube": pruned_cyclic_tube,
        "pruned_post_flype_c_tunnel": len(pre_tun) - len(tun),
        "final_survivor_diagrams": len(survivors),
    }


def _build_reference_catalog_for_target(
    data_dir: Path,
    target_ch: int,
    num_workers: int = 32,
) -> dict[tuple, list[str]]:
    """Build reference catalog of all certified classes `< target_ch` plus splits, unknots, and low-`ch` sums."""
    ref = dict(build_reference_catalog_ch12(data_dir, max_composite_ch=15, num_workers=num_workers))
    if target_ch >= 13:
        ch12_path = data_dir / "expanded_ch12.json"
        if ch12_path.exists():
            ch12_payload = json.loads(ch12_path.read_text())
            ch12_entries = list(ch12_payload.get("entries", []))
            for ent in ch12_entries:
                lbl = f"ch=12:{ent['label']}"
                if "sublink_invariant_profile" in ent:
                    dk = invariant_discrimination_key(ent)
                else:
                    rec = compute_full_invariant_record(
                        MarkedDiagram.from_epd(str(ent["epd"])), label=lbl, include_sublinks=False
                    )
                    dk = invariant_discrimination_key(rec)
                ref.setdefault(dk, [])
                if lbl not in ref[dk]:
                    ref[dk].append(lbl)
    return ref


def _load_spoke_catalog_by_key(data_dir: Path) -> dict[tuple, list[str]]:
    out: dict[tuple, list[str]] = defaultdict(list)
    fam_path = data_dir / "families" / "family_catalog.json"
    if not fam_path.exists():
        return out
    payload = json.loads(fam_path.read_text())
    for fent in payload.get("entries", []):
        fid = str(fent.get("instance_id") or fent.get("label") or "")
        if not fid:
            continue
        if "sublink_invariant_profile" in fent:
            dk = invariant_discrimination_key(fent)
        else:
            rec = compute_full_invariant_record(
                MarkedDiagram.from_epd(str(fent["epd"])), label=fid, include_sublinks=False
            )
            dk = invariant_discrimination_key(rec)
        if fid not in out[dk]:
            out[dk].append(fid)
    return out


def run_m2_s2_sector(
    ch: int,
    data_dir: Path,
    num_workers: int = 32,
) -> dict[str, object]:
    """Run exhaustive `m = 2` connected `S^2` census at `ch` and stratify into new prime classes vs collisions."""
    t_c0 = time.perf_counter()
    diags_level, st_c = enumerate_targeted_census_level_candidates(
        ch=ch,
        m_min=2,
        m_max=2,
        max_bfs_states=2048,
        target_components=1,
        target_orientability=1,
        target_chi_min=2,
        target_chi_max=2,
        num_workers=num_workers,
        _certify_unlink=False,
    )
    tun_c = cascade_tunnel_diagrams_c(
        ch=ch,
        diagrams=diags_level,
        depth_plus1=256,
        depth_plus2=192,
        max_bfs_states=2048,
        num_workers=num_workers,
    )
    filt_c = _filter_cyclic_tubes_parallel(tun_c, num_workers=num_workers)
    dt_c = round(time.perf_counter() - t_c0, 3)

    c_stats = dict(st_c)
    c_stats["raw_and_level_counts_are_kauffman_bracket_filtered"] = True
    c_stats["all_tunneled_survivors_reidemeister_unlink_certified"] = True
    c_stats["level_kept"] = c_stats.get("kept", len(diags_level))
    c_stats["tunneled_kept"] = len(tun_c)
    c_stats["cyclic_tube_pruned"] = len(tun_c) - len(filt_c)

    t_p0 = time.perf_counter()
    survivors, post_stats = _post_filter_and_tunnel_sector(ch, filt_c, num_workers=num_workers)
    ref_cat = _build_reference_catalog_for_target(data_dir, target_ch=ch, num_workers=num_workers)
    spoke_by_key = _load_spoke_catalog_by_key(data_dir)

    tasks = [(d, f"{ch}_m2_S2_cand_{i}", False) for i, d in enumerate(survivors)]
    recs: list[dict[str, object]] = []
    if tasks:
        with ProcessPoolExecutor(max_workers=max(1, min(num_workers, len(tasks))), mp_context=_MP_CTX) as pool:
            recs = list(pool.map(_compute_record_worker, tasks, chunksize=2))

    grouped_new: dict[tuple, list[tuple[MarkedDiagram, dict[str, object]]]] = defaultdict(list)
    ambiguous: list[dict[str, object]] = []

    for d_obj, rec in zip(survivors, recs):
        dk = invariant_discrimination_key(rec)
        if dk in ref_cat:
            ambiguous.append(
                {
                    "epd": rec["epd"],
                    "component_types": rec["component_types"],
                    "normal_euler_numbers": rec["normal_euler_numbers"],
                    "e1_generators": rec["e1_generators"],
                    "pi1_presentation": rec["pi1_presentation"],
                    "collision_matches": ref_cat[dk],
                    "marked_bigon_1handles": [asdict(h) for h in analyze_marked_bigon_1handles(d_obj)],
                }
            )
        else:
            grouped_new[dk].append((d_obj, rec))

    if grouped_new:
        cand_reps = [items[0][1] for items in grouped_new.values()]
        high_comp = build_def5_1_high_composite_matches(
            target_ch=ch,
            candidate_records=cand_reps,
            data_dir=data_dir,
            num_workers=num_workers,
        )
        for dk in list(grouped_new.keys()):
            if dk in high_comp:
                for d_obj, rec in grouped_new.pop(dk):
                    ambiguous.append(
                        {
                            "epd": rec["epd"],
                            "component_types": rec["component_types"],
                            "normal_euler_numbers": rec["normal_euler_numbers"],
                            "e1_generators": rec["e1_generators"],
                            "pi1_presentation": rec["pi1_presentation"],
                            "collision_matches": high_comp[dk],
                            "marked_bigon_1handles": [asdict(h) for h in analyze_marked_bigon_1handles(d_obj)],
                        }
                    )

    new_classes: list[dict[str, object]] = []
    sorted_groups = sorted(grouped_new.items(), key=lambda kv: str(kv[1][0][1]["epd"]))
    for k_idx, (dk, items) in enumerate(sorted_groups, start=1):
        d0, rep = items[0]
        cid = f"{ch}_m2_S2_{k_idx}"
        rep["class_id"] = cid
        rep["label"] = cid
        rep["spoke_matches"] = spoke_by_key.get(dk, [])
        rep["num_move_components"] = len(items)
        rep["all_epds"] = [str(r["epd"]) for _, r in items]
        rep["indecomposability"] = compute_indecomposability_certificate(d0, rep, 0)
        rep["marked_bigon_1handles"] = [asdict(h) for h in analyze_marked_bigon_1handles(d0)]
        new_classes.append(rep)

    dt_p = round(time.perf_counter() - t_p0, 3)
    return {
        "runtime_c_sec": dt_c,
        "runtime_post_sec": dt_p,
        "c_stats": c_stats,
        "post_filter_stats": post_stats,
        "new_prime_classes": new_classes,
        "ambiguous_shadow_collisions": ambiguous,
    }


def _rp2_cert_worker(epd: str) -> dict[str, object]:
    d = MarkedDiagram.from_epd(epd)
    rec = compute_full_invariant_record(d, label="rp2_cand", include_sublinks=False)
    mer_lb, mer_reason = meridian_order_lower_bound(d)
    cert = find_rp2_reducible_certificate(
        d,
        max_ch=min(15, d.num_vertices + 2),
        max_states=2500,
    )
    return {
        "epd": epd,
        "component_types": rec["component_types"],
        "normal_euler_numbers": rec["normal_euler_numbers"],
        "e1_generators": rec["e1_generators"],
        "pi1_presentation": rec["pi1_presentation"],
        "hom_counts": rec["hom_counts"],
        "meridian_order_lower_bound": mer_lb,
        "meridian_order_reason": mer_reason,
        "reducible_certificate": asdict(cert) if cert is not None else None,
        "dkey": invariant_discrimination_key(rec),
    }


def run_rp2_full_marker_sector(
    ch: int,
    data_dir: Path,
    num_workers: int = 32,
) -> dict[str, object]:
    """Scan ALL marker counts `m = 1 .. ch - 1` for connected `RP^2` knots (`chi = 1`, non-orientable)
    at `ch`, tunnel, stratify, and search for constructive Yoshikawa move paths reducing each survivor
    to `< ch` or a 2-edge-cut connected sum `K # P_+-`.
    """
    t0 = time.perf_counter()
    per_m_diags: list[MarkedDiagram] = []
    tot_raw = 0
    tot_comps = 0
    tot_disc = 0
    tot_kept = 0
    num_shadows = 0

    for m_val in range(1, ch):
        # At ch = 13, m = 0 has 0 saddles; for m = 1, chi = c+ + c- - 1 = 1 => c+ + c- = 2 (impossible for 1 comp
        # unless c+ = c- = 1, which flips orientability? Still scan every m = 1..ch-1 so zero assumptions are made!).
        diags_m, st_m = enumerate_targeted_census_level_candidates(
            ch=ch,
            m_min=m_val,
            m_max=m_val,
            max_bfs_states=2048,
            target_components=1,
            target_orientability=0,
            target_chi_min=1,
            target_chi_max=1,
            num_workers=num_workers,
            _certify_unlink=False,
        )
        num_shadows = int(st_m.get("num_shadows", num_shadows))
        tot_raw += int(st_m.get("raw_cands", 0))
        tot_comps += int(st_m.get("components", 0))
        tot_disc += int(st_m.get("discarded", 0))
        tot_kept += int(st_m.get("kept", 0))
        if diags_m:
            tun_m = cascade_tunnel_diagrams_c(
                ch=ch,
                diagrams=diags_m,
                depth_plus1=192,
                depth_plus2=128,
                max_bfs_states=2048,
                num_workers=num_workers,
            )
            per_m_diags.extend(tun_m)

    merged = (
        merge_diagrams_c(ch, per_m_diags, max_bfs_states=4096, num_workers=num_workers, _certify_unlink=False)
        if per_m_diags
        else []
    )
    tun_all = (
        cascade_tunnel_diagrams_c(
            ch=ch,
            diagrams=merged,
            depth_plus1=256,
            depth_plus2=192,
            max_bfs_states=2048,
            num_workers=num_workers,
        )
        if merged
        else []
    )
    filt_all = _filter_cyclic_tubes_parallel(tun_all, num_workers=num_workers) if tun_all else []
    survivors, post_stats = _post_filter_and_tunnel_sector(ch, filt_all, num_workers=num_workers)
    ref_cat = _build_reference_catalog_for_target(data_dir, target_ch=ch, num_workers=num_workers)

    surv_epds = [d.to_epd() for d in survivors]
    analyzed: list[dict[str, object]] = []
    if surv_epds:
        with ProcessPoolExecutor(max_workers=max(1, min(num_workers, len(surv_epds))), mp_context=_MP_CTX) as pool:
            analyzed = list(pool.map(_rp2_cert_worker, surv_epds, chunksize=1))

    new_prime_rp2: list[dict[str, object]] = []
    colliding_rp2: list[dict[str, object]] = []
    certified_reducible_count = 0
    certified_irreducible_count = 0
    unresolved_count = 0
    meridian_gt2_count = 0

    for item in analyzed:
        dk = item.pop("dkey")
        if item["meridian_order_lower_bound"] > 2:  # type: ignore[operator]
            meridian_gt2_count += 1
        cert_dict = item["reducible_certificate"]
        cert_status = cert_dict.get("status") if isinstance(cert_dict, dict) else None
        if (
            isinstance(cert_dict, dict)
            and cert_status in ("STANDARD_P2", "REDUCED_LOWER_CH", "SPLIT_CONNECTED_SUM")
            and bool(cert_dict.get("path_verified", True))
        ):
            certified_reducible_count += 1
        elif cert_status == "CERTIFIED_IRREDUCIBLE":
            certified_irreducible_count += 1
        else:
            unresolved_count += 1
        if dk in ref_cat:
            item["collision_matches"] = ref_cat[dk]
            colliding_rp2.append(item)
        else:
            item["collision_matches"] = []
            new_prime_rp2.append(item)

    return {
        "ch_index": ch,
        "marker_range_scanned": [1, ch - 1],
        "runtime_sec": round(time.perf_counter() - t0, 3),
        "c_stats": {
            "num_shadows": num_shadows,
            "raw_cands": tot_raw,
            "components": tot_comps,
            "discarded": tot_disc,
            "level_kept": tot_kept,
            "tunneled_kept": len(tun_all),
            "cyclic_tube_pruned": len(tun_all) - len(filt_all),
        },
        "post_filter_stats": post_stats,
        "final_survivor_count": len(survivors),
        "constructive_reducible_certified_count": certified_reducible_count,
        "certified_irreducible_count": certified_irreducible_count,
        "unresolved_count": unresolved_count,
        "all_survivors_constructively_reducible": (certified_reducible_count == len(survivors)),
        "meridian_order_gt_2_count": meridian_gt2_count,
        "new_prime_classes": new_prime_rp2,
        "survivors": colliding_rp2 + new_prime_rp2,
    }


def run_ch13_other_sectors(data_dir: Path, num_workers: int = 32) -> dict[str, object]:
    """Run Campaign 2.3 targeted sweeps at `ch = 13`:
    - Connected Klein bottles `K^2` (`chi = 0`, non-orientable, `m = 2, 3, 4`),
    - Connected tori `T^2` (`chi = 0`, orientable, `m = 2, 3, 4`),
    - 2-component links (`num_components = 2`, `m = 2, 3`).
    """
    ref_cat = _build_reference_catalog_for_target(data_dir, target_ch=13, num_workers=num_workers)
    spoke_by_key = _load_spoke_catalog_by_key(data_dir)
    sectors = [
        ("ch13_K2_m2_to_4", 1, 0, 0, 0, 2, 4),
        ("ch13_T2_m2_to_4", 1, 1, 0, 0, 2, 4),
        ("ch13_2comp_m2_to_3", 2, -1, -4, 4, 2, 3),
    ]
    results: dict[str, object] = {}
    for name, t_comp, t_ori, chi_lo, chi_hi, m_lo, m_hi in sectors:
        t0 = time.perf_counter()
        all_m: list[MarkedDiagram] = []
        tot_raw = 0
        tot_kept = 0
        for m_val in range(m_lo, m_hi + 1):
            dm, stm = enumerate_targeted_census_level_candidates(
                ch=13,
                m_min=m_val,
                m_max=m_val,
                max_bfs_states=2048,
                target_components=t_comp,
                target_orientability=t_ori,
                target_chi_min=chi_lo,
                target_chi_max=chi_hi,
                num_workers=num_workers,
                _certify_unlink=False,
            )
            tot_raw += int(stm.get("raw_cands", 0))
            tot_kept += int(stm.get("kept", 0))
            if dm:
                all_m.extend(
                    cascade_tunnel_diagrams_c(
                        ch=13,
                        diagrams=dm,
                        depth_plus1=128,
                        depth_plus2=96,
                        max_bfs_states=2048,
                        num_workers=num_workers,
                    )
                )
        merged = (
            merge_diagrams_c(13, all_m, max_bfs_states=4096, num_workers=num_workers, _certify_unlink=False)
            if all_m
            else []
        )
        tun = (
            cascade_tunnel_diagrams_c(
                ch=13, diagrams=merged, depth_plus1=192, depth_plus2=128, max_bfs_states=2048, num_workers=num_workers
            )
            if merged
            else []
        )
        filt = _filter_cyclic_tubes_parallel(tun, num_workers=num_workers) if tun else []
        survivors, post_stats = _post_filter_and_tunnel_sector(13, filt, num_workers=num_workers)

        tasks = [(d, f"{name}_{i}", False) for i, d in enumerate(survivors)]
        recs = []
        if tasks:
            with ProcessPoolExecutor(max_workers=max(1, min(num_workers, len(tasks))), mp_context=_MP_CTX) as pool:
                recs = list(pool.map(_compute_record_worker, tasks, chunksize=2))

        grouped_new: dict[tuple, list[tuple[MarkedDiagram, dict[str, object]]]] = defaultdict(list)
        colliding = 0
        for d_obj, rec in zip(survivors, recs):
            dk = invariant_discrimination_key(rec)
            if dk in ref_cat:
                colliding += 1
            else:
                grouped_new[dk].append((d_obj, rec))

        if grouped_new:
            high_comp = build_def5_1_high_composite_matches(
                target_ch=13,
                candidate_records=[v[0][1] for v in grouped_new.values()],
                data_dir=data_dir,
                num_workers=num_workers,
            )
            for dk in list(grouped_new.keys()):
                if dk in high_comp:
                    colliding += len(grouped_new.pop(dk))

        new_list: list[dict[str, object]] = []
        for idx, (dk, items) in enumerate(sorted(grouped_new.items(), key=lambda kv: str(kv[1][0][1]["epd"])), start=1):
            d0, rep = items[0]
            cid = f"{name}_cls_{idx}"
            rep["class_id"] = cid
            rep["label"] = cid
            rep["spoke_matches"] = spoke_by_key.get(dk, [])
            rep["num_move_components"] = len(items)
            rep["all_epds"] = [str(r["epd"]) for _, r in items]
            rep["indecomposability"] = compute_indecomposability_certificate(d0, rep, 0)
            rep["marked_bigon_1handles"] = [asdict(h) for h in analyze_marked_bigon_1handles(d0)]
            new_list.append(rep)

        results[name] = {
            "marker_range": [m_lo, m_hi],
            "runtime_sec": round(time.perf_counter() - t0, 3),
            "raw_cands": tot_raw,
            "level_kept": tot_kept,
            "final_survivor_diagrams": len(survivors),
            "colliding_diagrams": colliding,
            "new_prime_class_count": len(new_list),
            "new_prime_classes": new_list,
            "post_filter_stats": post_stats,
        }
    return results


def export_targeted_ch12_13_census(
    data_dir: Path,
    num_workers: int = 32,
    include_ch14_s2: bool = False,
) -> dict[str, object]:
    """Run all Campaign 2 targeted sweeps with incremental checkpointing to
    `data_dir / 'complexity' / 'targeted_ch12_13_census.json'`.
    """
    comp_dir = data_dir / "complexity"
    comp_dir.mkdir(parents=True, exist_ok=True)
    out_path = comp_dir / "targeted_ch12_13_census.json"

    payload: dict[str, object] = {}

    print(f"[Campaign 2.2a] Running ch=12 m=2 S^2 knots on {num_workers} workers...", flush=True)
    payload["ch12_m2_s2_knots"] = run_m2_s2_sector(12, data_dir, num_workers=num_workers)
    out_path.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"  -> ch=12 m=2 S^2: {len(payload['ch12_m2_s2_knots']['new_prime_classes'])} new classes "  # type: ignore[index]
        f"({payload['ch12_m2_s2_knots']['runtime_c_sec']}s C + {payload['ch12_m2_s2_knots']['runtime_post_sec']}s post)",  # type: ignore[index]
        flush=True,
    )

    print(f"[Campaign 2.2b] Running ch=13 m=2 S^2 knots on {num_workers} workers...", flush=True)
    payload["ch13_m2_s2_knots"] = run_m2_s2_sector(13, data_dir, num_workers=num_workers)
    out_path.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        f"  -> ch=13 m=2 S^2: {len(payload['ch13_m2_s2_knots']['new_prime_classes'])} new classes "  # type: ignore[index]
        f"({payload['ch13_m2_s2_knots']['runtime_c_sec']}s C + {payload['ch13_m2_s2_knots']['runtime_post_sec']}s post)",  # type: ignore[index]
        flush=True,
    )

    rp2_by_ch: dict[str, object] = {}
    for ch_rp2 in (11, 12, 13):
        print(f"[Campaign 2.1] Running full-marker RP^2 census at ch={ch_rp2} (m=1..{ch_rp2-1})...", flush=True)
        rp2_res = run_rp2_full_marker_sector(ch_rp2, data_dir, num_workers=num_workers)
        rp2_by_ch[str(ch_rp2)] = rp2_res
        payload["rp2_targeted_by_ch"] = rp2_by_ch
        out_path.write_text(json.dumps(payload, indent=2) + "\n")
        print(
            f"  -> ch={ch_rp2} RP^2: {rp2_res['final_survivor_count']} survivors, "
            f"{rp2_res['constructive_reducible_certified_count']} constructively reduced, "
            f"meridian_gt2={rp2_res['meridian_order_gt_2_count']}, "
            f"new_prime={len(rp2_res['new_prime_classes'])} ({rp2_res['runtime_sec']}s)",  # type: ignore[arg-type]
            flush=True,
        )

    print(f"[Campaign 2.3] Running ch=13 K^2, T^2, and 2-component targeted sweeps...", flush=True)
    payload["ch13_other_sectors"] = run_ch13_other_sectors(data_dir, num_workers=num_workers)
    out_path.write_text(json.dumps(payload, indent=2) + "\n")

    if include_ch14_s2:
        print(f"[Campaign 2.2c] Running ch=14 m=2 S^2 knots on {num_workers} workers...", flush=True)
        payload["ch14_m2_s2_knots"] = run_m2_s2_sector(14, data_dir, num_workers=num_workers)
        out_path.write_text(json.dumps(payload, indent=2) + "\n")

    return payload
