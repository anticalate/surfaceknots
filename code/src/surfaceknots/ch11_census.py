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


def build_reference_catalog(max_composite_ch: int = 12) -> dict[tuple, list[str]]:
    """Build an exhaustive invariant lookup table of:
    1. All 23 `ch <= 10` Yoshikawa prime surface-links (`ch<=10:{label}`),
    2. Standard unknotted connected surfaces up to `ch <= 12` (`unknot:{type}(e={e})`),
    3. Split unions of unknots and `ch <= 10` entries up to `ch <= 11` (`split:{A} u {B}`),
    4. Nontrivial connected sums `A # B` up to `ch <= max_composite_ch` (`composite:{A} # {B}`).
    """
    catalog: dict[tuple, list[str]] = defaultdict(list)
    queued: list[tuple[MarkedDiagram, str, bool]] = []

    def register(d: MarkedDiagram, tag: str) -> None:
        if not _is_admissible(d):
            return
        queued.append((d, tag, False))

    # 1. All 23 Yoshikawa entries (`ch <= 10`):
    yosh_diags: list[tuple[str, MarkedDiagram]] = []
    for entry in YOSHIKAWA_TABLE:
        d = MarkedDiagram.from_epd(entry.canonical_epd)
        register(d, f"ch<=10:{entry.label}")
        yosh_diags.append((entry.label, d))

    # 2. Build all connected unknotted surfaces up to ch = 12 by iterated connected sums of `2_1^1`, `P_+`, `P_-`:
    t2 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^1"].canonical_epd)
    rp2_pos = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^-1"].canonical_epd)
    rp2_neg = rp2_pos.mirror_3d()

    unknot_atoms: list[tuple[str, MarkedDiagram]] = [
        ("T^2", t2),
        ("P+", rp2_pos),
        ("P-", rp2_neg),
    ]
    # Iteratively build connected unknots up to ch = 10 (`g <= 5` orientable, `k <= 5` non-orientable):
    unknots_by_topo: dict[tuple[int, int], tuple[str, MarkedDiagram]] = {
        (0, 0): ("0_1", MarkedDiagram.unknot_sphere()),
        (1, 0): ("2_1^1", t2),
        (-1, 2): ("P+", rp2_pos),
        (-1, -2): ("P-", rp2_neg),
    }
    frontier = [(1, 0), (-1, 2), (-1, -2)]
    while frontier:
        g_cur, e_cur = frontier.pop(0)
        lbl_cur, d_cur = unknots_by_topo[(g_cur, e_cur)]
        if d_cur.ch_index + 2 > 10:
            continue
        for atom_lbl, d_atom in unknot_atoms:
            cand = planar_connected_sum(d_cur, 0, d_atom, 0, flip_splice=False)
            if not _is_admissible(cand):
                cand = planar_connected_sum(d_cur, 0, d_atom, 0, flip_splice=True)
            if not _is_admissible(cand):
                continue
            ct = cand.surface_topology().components[0]
            sig = (ct.genus_symbol, ct.normal_euler)
            if sig not in unknots_by_topo:
                new_lbl = f"unknot({ct.homeomorphism_type},e={ct.normal_euler})"
                unknots_by_topo[sig] = (new_lbl, cand)
                register(cand, f"unknot:{ct.homeomorphism_type}(e={ct.normal_euler})")
                frontier.append(sig)

    # 3. Nontrivial connected sums `A # U` and `A # B` with `ch(A) + ch(B) <= max_composite_ch`:
    # Summands: all nontrivial Yoshikawa entries (`ch >= 2`) + unknots (`2 <= ch <= 6`)!
    summands: list[tuple[str, MarkedDiagram]] = []
    for entry in YOSHIKAWA_TABLE:
        if entry.ch_index == 0 or entry.ch_index > max_composite_ch - 2:
            continue
        d = MarkedDiagram.from_epd(entry.canonical_epd)
        summands.append((entry.label, d))
        if entry.label in ("2_1^-1", "10_2"):
            summands.append((f"{entry.label}*", d.mirror_3d()))
    for (gs, es), (ulbl, ud) in unknots_by_topo.items():
        if 4 <= ud.ch_index <= 6:
            summands.append((ulbl, ud))

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
            if (lbl_a.startswith("2_1") or lbl_a.startswith("unknot")) and (
                lbl_b.startswith("2_1") or lbl_b.startswith("unknot")
            ):
                continue
            dc_b, _ = db._dart_surface_components()
            rep_dart_b = {}
            for hd in range(4 * db.num_vertices):
                rep_dart_b.setdefault(dc_b[hd], hd)

            for ca, hda in rep_dart_a.items():
                for cb, hdb in rep_dart_b.items():
                    for fl in (False, True):
                        csum = planar_connected_sum(da, hda, db, hdb, flip_splice=fl)
                        if _is_admissible(csum):
                            register(csum, f"composite:{lbl_a}#{lbl_b}(c{ca}+c{cb})")
                            break

    # 4. Split unions `A u B` (and `A u 0_1 u 0_1`, `unknot u unknot`) up to ch <= 11:
    split_blocks: list[tuple[str, MarkedDiagram]] = [("0_1", MarkedDiagram.unknot_sphere())]
    for (gs, es), (ulbl, ud) in unknots_by_topo.items():
        if 2 <= ud.ch_index <= 8:
            split_blocks.append((ulbl, ud))
    for entry in YOSHIKAWA_TABLE:
        if 6 <= entry.ch_index <= 10:
            split_blocks.append((entry.label, MarkedDiagram.from_epd(entry.canonical_epd)))

    for i, (la, da) in enumerate(split_blocks):
        for j in range(i, len(split_blocks)):
            lb, db = split_blocks[j]
            if da.ch_index + db.ch_index > 11:
                continue
            sp2 = disjoint_union(da, db)
            register(sp2, f"split:{la}_u_{lb}")
            # Also 3-component splits with 0_1 up to ch <= 10:
            if sp2.ch_index <= 9 and _is_admissible(sp2) and sp2.surface_topology().num_components == 2:
                register(disjoint_union(sp2, MarkedDiagram.unknot_sphere()), f"split:{la}_u_{lb}_u_0_1")

    from concurrent.futures import ProcessPoolExecutor

    with ProcessPoolExecutor(max_workers=2) as pool:
        for (_, tag, _), rec in zip(queued, pool.map(_compute_record_worker, queued, chunksize=8)):
            k = invariant_discrimination_key(rec)
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


def stratify_and_export_ch11(
    survivor_epds: list[str],
    out_dir: Path,
    census_stats: dict[str, object] | None = None,
) -> dict[str, object]:
    """Stratify `ch = 11` survivors against the exhaustive `ch <= 10` / composite / split / unknot catalog
    and export `expanded_ch11.json` and `expanded_ch11.csv`.
    """
    from concurrent.futures import ProcessPoolExecutor
    from surfaceknots.yoshikawa_data import export_yoshikawa_datasets

    out_dir.mkdir(parents=True, exist_ok=True)
    yosh_json = out_dir / "yoshikawa_ch10.json"
    if yosh_json.exists():
        yosh_records: list[dict[str, object]] = json.loads(yosh_json.read_text())
        if not yosh_records or "sigma2_euler" not in yosh_records[0]:
            yosh_records = export_yoshikawa_datasets(out_dir)
    else:
        yosh_records = export_yoshikawa_datasets(out_dir)

    survivors, post_stats = post_filter_and_merge_ch11_survivors(survivor_epds)
    ref_catalog = build_reference_catalog(max_composite_ch=12)

    # Compute full invariants on every surviving `ch = 11` diagram in parallel:
    grouped_new: dict[tuple, list[dict[str, object]]] = defaultdict(list)
    grouped_residual: dict[tuple, list[dict[str, object]]] = defaultdict(list)
    residual_collisions: list[dict[str, object]] = []

    surv_tasks = [(d, f"ch11_cand_{idx}", True) for idx, d in enumerate(survivors)]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for rec in pool.map(_compute_record_worker, surv_tasks, chunksize=4):
            dkey = invariant_discrimination_key(rec)
            if dkey in ref_catalog:
                rec["cell_status"] = "AMBIGUOUS_SHADOW"
                rec["collision_matches"] = ref_catalog[dkey]
                residual_collisions.append(rec)
                grouped_residual[dkey].append(rec)
            else:
                grouped_new[dkey].append(rec)

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
