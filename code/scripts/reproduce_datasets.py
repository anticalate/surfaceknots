"""Reproducible CLI pipeline to regenerate the SurfaceKnots Hub, Spoke, Complexity, Gallery, and Web Atlas datasets.

Usage from `code/`:
    PYTHONPATH=src .venv/bin/python scripts/reproduce_datasets.py --yoshikawa
    PYTHONPATH=src .venv/bin/python scripts/reproduce_datasets.py --ch11-enumerate --ch11-export
    PYTHONPATH=src .venv/bin/python scripts/reproduce_datasets.py --ch12-enumerate --ch12-export
    PYTHONPATH=src .venv/bin/python scripts/reproduce_datasets.py --families --gallery --extended-gallery
    PYTHONPATH=src .venv/bin/python scripts/reproduce_datasets.py --triplane --complexity --atlas
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from surfaceknots.ch11_census import stratify_and_export_ch11
from surfaceknots.triplane import (
    TangleDiagram,
    TriPlaneDiagram,
    generate_catalan_matchings,
    standard_rp2_triplane,
    standard_torus_triplane,
    standard_unknot_s2_triplane,
)
from surfaceknots.yoshikawa_data import export_yoshikawa_datasets


DEFAULT_DATA_DIR = Path(__file__).resolve().parents[1] / "data"


def run_b_le_4_triplane_checks() -> dict[str, object]:
    """Run exhaustive concentrated low-crossing tri-plane checks for b <= 3 and representative b = 4 models."""
    out: dict[str, object] = {}
    c0_by_b: dict[int, list[tuple[tuple[int, int, int], int]]] = {}
    for b in (1, 2, 3):
        mats = generate_catalan_matchings(b)
        found = set()
        for m1 in mats:
            t1 = TangleDiagram(b, m1, ())
            for m2 in mats:
                t2 = TangleDiagram(b, m2, ())
                for m3 in mats:
                    t3 = TangleDiagram(b, m3, ())
                    tp = TriPlaneDiagram(t1, t2, t3)
                    if tp.is_certified_valid():
                        found.add((tuple(sorted(tp.patch_counts())), tp.euler_characteristic()))
        c0_by_b[b] = sorted(found)
    out["crossingless_b1_to_3"] = {str(k): v for k, v in c0_by_b.items()}

    conc_by_b: dict[int, list[tuple[int, tuple[int, int, int], int]]] = {}
    for b in (2, 3):
        mats = generate_catalan_matchings(b)
        gens = [g for k in range(1, 2 * b) for g in (k, -k)]
        braids: list[tuple[int, ...]] = [()]
        for g1 in gens:
            braids.append((g1,))
            for g2 in gens:
                if g1 != -g2:
                    braids.append((g1, g2))
        found_c = set()
        for bw in braids:
            for m1 in mats:
                t1 = TangleDiagram(b, m1, bw)
                for m2 in mats:
                    t2 = TangleDiagram(b, m2, ())
                    for m3 in mats:
                        t3 = TangleDiagram(b, m3, ())
                        tp = TriPlaneDiagram(t1, t2, t3)
                        if tp.is_certified_valid():
                            found_c.add((len(bw), tuple(sorted(tp.patch_counts())), tp.euler_characteristic()))
        conc_by_b[b] = sorted(found_c)
    out["concentrated_c_le_2_b2_b3"] = {str(k): v for k, v in conc_by_b.items()}

    m4 = generate_catalan_matchings(4)
    b4_valid_6 = 0
    for m_a in m4:
        for m_b in m4:
            for bw in ((2, 2, 2, -6, -6, -6), (2, -1, 2, -6, 7, -6)):
                tp = TriPlaneDiagram(TangleDiagram(4, m_a, bw), TangleDiagram(4, m_a, ()), TangleDiagram(4, m_b, ()))
                if tp.patch_counts() == (2, 2, 2) and tp.is_certified_valid():
                    b4_valid_6 += 1
    out["b4_spun_trefoil_4_222_valid_models"] = b4_valid_6
    out["standard_models_valid"] = (
        standard_unknot_s2_triplane().is_certified_valid()
        and standard_rp2_triplane(1).is_certified_valid()
        and standard_torus_triplane().is_certified_valid()
    )
    return out


def collect_census_stats(data_dir: Path) -> dict[str, object]:
    sh_path = data_dir / "_ch11_shards.json"
    if sh_path.exists():
        sh_state = json.loads(sh_path.read_text())
        tot_raw = sum(v[0] for v in sh_state["stats"].values())
        tot_comps = sum(v[1] for v in sh_state["stats"].values())
        tot_disc = sum(v[2] for v in sh_state["stats"].values())
    else:
        tot_raw, tot_comps, tot_disc = -1, -1, -1

    def _cnt(name: str) -> int:
        p = data_dir / name
        return len(json.loads(p.read_text())) if p.exists() else -1

    return {
        "prime_shadows_oeis_a113201": 3378,
        "raw_kauffman_1step_survivors": tot_raw,
        "level_components_explored": tot_comps,
        "level_components_discarded_in_shards": tot_disc,
        "level_irreducible_global_merged": _cnt("_level_kept_ch11.json"),
        "passA_plus1_d1_kept": _cnt("_stageA_ch11.json"),
        "passB_plus12_d2_kept": _cnt("_stageB_ch11.json"),
        "passC1_plus12_96_48_kept": _cnt("_stageC1_ch11.json"),
        "passC2_plus12_256_192_kept": _cnt("_census_ch11_all_tunneled.json"),
        "cyclic_tube_filtered_kept": _cnt("_census_ch11.json"),
        "triplane_b_le_4_crosscheck": run_b_le_4_triplane_checks(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Regenerate surfaceknots census datasets and diagram gallery.")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--pics-dir", type=Path, default=DEFAULT_DATA_DIR.parents[1] / "pics")
    parser.add_argument(
        "--atlas-out",
        type=Path,
        default=None,
        help="Output path for interactive Web Atlas HTML (default: expo/surfaceknots_atlas.html)",
    )
    parser.add_argument("--workers", type=int, default=32, help="Worker thread/process count for parallel stages")
    parser.add_argument("--yoshikawa", action="store_true", help="Export yoshikawa_ch10.{json,csv}")
    parser.add_argument("--ch11-enumerate", action="store_true", help="Run checkpointed ch=11 enumeration/tunneling")
    parser.add_argument("--ch11-export", action="store_true", help="Stratify & export expanded_ch11.{json,csv}")
    parser.add_argument("--ch12-enumerate", action="store_true", help="Run checkpointed ch=12 enumeration/tunneling")
    parser.add_argument("--ch12-export", action="store_true", help="Stratify & export expanded_ch12.{json,csv}")
    parser.add_argument("--families", action="store_true", help="Build & export parametric family catalogs and Hub join")
    parser.add_argument("--gallery", action="store_true", help="Render SVG & TikZ gallery into pics/{svg,tikz}")
    parser.add_argument(
        "--extended-gallery",
        action="store_true",
        help="Render splitting movie triptychs & showcase Spoke family SVG/TikZ into pics/",
    )
    parser.add_argument(
        "--targeted-census",
        action="store_true",
        help="Run Campaign 2 targeted ch=12..13 S^2, RP^2, and T^2/link sectors",
    )
    parser.add_argument(
        "--move-search",
        action="store_true",
        help="Run Campaign 3 bidirectional Yoshikawa move searches (Jablonowski, Viro, Klein tube)",
    )
    parser.add_argument(
        "--triplane",
        action="store_true",
        help="Run Campaign 4/6.1 tri-plane censuses (braided b<=4 up to c_tri<=7, concentrated b=5 c_tri<=4, general planar b=4 c<=3)",
    )
    parser.add_argument(
        "--complexity",
        action="store_true",
        help="Export complexity_comparison.{json,csv} across Hub, ch=12, Spoke, and tri-plane catalogs",
    )
    parser.add_argument(
        "--atlas",
        action="store_true",
        help="Export self-contained interactive Web Atlas to expo/surfaceknots_atlas.html",
    )
    args = parser.parse_args()

    stage_flags = (
        args.yoshikawa,
        args.ch11_enumerate,
        args.ch11_export,
        args.ch12_enumerate,
        args.ch12_export,
        args.families,
        args.gallery,
        args.extended_gallery,
        args.targeted_census,
        args.move_search,
        args.triplane,
        args.complexity,
        args.atlas,
    )
    if not any(stage_flags):
        parser.print_help()
        return

    if args.yoshikawa:
        export_yoshikawa_datasets(args.data_dir)
        print(f"Exported yoshikawa_ch10.{{json,csv}} to {args.data_dir}")
    if args.ch11_enumerate:
        from surfaceknots.ch11_census import enumerate_ch11_checkpointed

        epds11, st11 = enumerate_ch11_checkpointed(args.data_dir, num_workers=args.workers)
        print(
            f"Completed ch=11 enumeration (m=1..10): {len(epds11)} post-tube survivors "
            f"(level_kept={st11['level_kept']}, wall={st11['total_enumeration_wall_sec']}s)."
        )
    if args.ch11_export:
        from surfaceknots.ch11_census import enumerate_ch11_checkpointed

        ckpt11 = args.data_dir / "ch11_raw_survivors.json"
        census_cache = args.data_dir / "_census_ch11.json"
        existing_json = args.data_dir / "expanded_ch11.json"
        if ckpt11.exists():
            c11 = json.loads(ckpt11.read_text())
            epds = list(c11["survivor_epds"])
            c_stats = dict(c11.get("census_stats", {}))
            c_stats["triplane_b_le_4_crosscheck"] = run_b_le_4_triplane_checks()
        elif census_cache.exists():
            epds = json.loads(census_cache.read_text())
            c_stats = collect_census_stats(args.data_dir)
        elif existing_json.exists():
            prev = json.loads(existing_json.read_text())
            epds = (
                [e["epd"] for e in prev["entries"]]
                + [a for e in prev["entries"] for a in e.get("alternate_epds", [])]
                + [r["epd"] for r in prev["residual_collisions"]]
            )
            c_stats = prev.get("census_stats") or collect_census_stats(args.data_dir)
        else:
            epds, c_stats = enumerate_ch11_checkpointed(args.data_dir, num_workers=args.workers)
            c_stats["triplane_b_le_4_crosscheck"] = run_b_le_4_triplane_checks()
        res = stratify_and_export_ch11(epds, args.data_dir, census_stats=c_stats, num_workers=args.workers)
        print(
            f"Exported expanded_ch11.{{json,csv}}: "
            f"{res['certified_new_prime_invariant_classes']} certified new prime classes, "
            f"{res['residual_colliding_diagrams']} residual collisions."
        )
    if args.ch12_enumerate:
        from surfaceknots.ch12_census import enumerate_ch12_checkpointed

        epds12, st12 = enumerate_ch12_checkpointed(args.data_dir, num_workers=args.workers)
        print(
            f"Completed ch=12 enumeration: {len(epds12)} post-tube survivors "
            f"(level_kept={st12['level_kept']}, wall={st12['total_enumeration_wall_sec']}s)."
        )
    if args.ch12_export:
        from surfaceknots.ch12_census import enumerate_ch12_checkpointed, stratify_and_export_ch12

        ckpt12 = args.data_dir / "ch12_raw_survivors.json"
        if ckpt12.exists():
            c12 = json.loads(ckpt12.read_text())
            epds12 = list(c12["survivor_epds"])
            st12 = dict(c12.get("census_stats", {}))
        else:
            epds12, st12 = enumerate_ch12_checkpointed(args.data_dir, num_workers=args.workers)
        res12 = stratify_and_export_ch12(
            epds12,
            args.data_dir,
            census_stats=st12,
            num_workers=args.workers,
            max_composite_ch=15,
        )
        print(
            f"Exported expanded_ch12.{{json,csv}}: "
            f"{res12['certified_new_prime_invariant_classes']} certified new classes "
            f"({res12['exact_singletons_count']} singletons, {res12['new_packets_count']} packets), "
            f"{res12['residual_colliding_diagrams']} residual collisions."
        )
    if args.families:
        from surfaceknots.families import export_family_datasets

        fam_res = export_family_datasets(args.data_dir)
        ex12 = fam_res.get("exact_ch_12_certified_instances", [])
        ex13 = fam_res.get("exact_ch_13_certified_instances", [])
        print(
            f"Exported family_catalog.{{json,csv}} & hub_realizations.json: "
            f"{fam_res['total_family_instances']} instances, "
            f"{fam_res['distinct_hub_entries_realized']} Hub entries realized (hub_max_ch={fam_res.get('hub_max_ch', 11)}), "
            f"exact_ch12={len(ex12)} {ex12}, exact_ch13={len(ex13)} {ex13}."
        )
    if args.gallery:
        from surfaceknots.drawing import export_canonical_gallery

        manifest = export_canonical_gallery(args.data_dir, args.pics_dir)
        print(f"Exported {len(manifest)} SVG + TikZ diagram pairs to {args.pics_dir}")
    if args.extended_gallery:
        from surfaceknots.drawing import export_extended_gallery

        ext = export_extended_gallery(args.data_dir, args.pics_dir)
        print(
            f"Exported extended gallery: {len(ext['canonical'])} canonical, "
            f"{len(ext['families'])} family SVG+TikZ, {len(ext['movies'])} splitting movies."
        )
    if args.targeted_census:
        from surfaceknots.targeted_census import export_targeted_ch12_13_census

        tc_res = export_targeted_ch12_13_census(args.data_dir, num_workers=args.workers)
        print(
            f"Exported targeted_ch12_13_census.json: "
            f"ch12_m2={len(tc_res['ch12_m2_s2_knots']['new_prime_classes'])} classes, "
            f"ch13_m2={len(tc_res['ch13_m2_s2_knots']['new_prime_classes'])} classes."
        )
    if args.move_search:
        from surfaceknots.move_search import export_move_search_results

        ms_res = export_move_search_results(args.data_dir, num_workers=args.workers)
        print(
            f"Exported move_search_results.json in {ms_res['total_wall_sec']}s: "
            + ", ".join(
                f"{k}={v['connected']}" for k, v in ms_res["targets"].items()  # type: ignore[union-attr]
            )
        )
    if args.triplane:
        from surfaceknots.triplane_census import run_triplane_b4_census

        tp_res = run_triplane_b4_census(args.data_dir, max_c_tri=6, num_workers=args.workers)
        print(
            f"Exported triplane_b4_census.json: {tp_res['distinct_o4_classes_count']} O(4) classes "
            f"({tp_res['runtime_enum_sec']}s enum + {tp_res['runtime_invariants_sec']}s inv)."
        )
    if args.complexity:
        from surfaceknots.complexity import export_complexity_comparison_datasets

        comp_res = export_complexity_comparison_datasets(args.data_dir)
        print(
            f"Exported complexity_comparison.{{json,csv}} ({comp_res['total_distinct_classes_compared']} classes)."
        )
    if args.atlas:
        from surfaceknots.atlas import export_interactive_atlas

        out_atlas = args.atlas_out or (args.data_dir.parents[1] / "expo" / "surfaceknots_atlas.html")
        counts = export_interactive_atlas(args.data_dir, out_atlas)
        print(f"Exported interactive Web Atlas ({counts['total_distinct_classes']} classes) to {out_atlas}")


if __name__ == "__main__":
    main()
