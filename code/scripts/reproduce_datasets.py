"""Reproducible CLI pipeline to regenerate the Yoshikawa (ch <= 10) and expanded (ch = 11) datasets.

Usage from `code/`:
    PYTHONPATH=src .venv/bin/python scripts/reproduce_datasets.py --yoshikawa
    PYTHONPATH=src .venv/bin/python scripts/reproduce_datasets.py --ch11-enumerate --max-seconds 3600
    PYTHONPATH=src .venv/bin/python scripts/reproduce_datasets.py --ch11-export
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import time

from surfaceknots.ch11_census import stratify_and_export_ch11
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.fast_enum import (
    _get_lib,
    _run_c_census_shard,
    merge_diagrams_c,
    tunnel_diagrams_c,
)
from surfaceknots.moves import has_certified_cyclic_unknotted_tube
from surfaceknots.quadrangulations import generate_prime_shadows
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


def run_ch11_enumeration(data_dir: Path, max_seconds: float = 3600.0) -> bool:
    """Run checkpointed ch = 11 shadow enumeration and multi-pass tunneling.

    Returns True when `_census_ch11.json` is complete.
    """
    data_dir.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    _get_lib()
    final_path = data_dir / "_census_ch11.json"
    if final_path.exists():
        epds = json.loads(final_path.read_text())
        print(f"[DONE] {final_path.name} already complete ({len(epds)} survivors).")
        return True

    sh_path = data_dir / "_shadows_ch11.json"
    if not sh_path.exists():
        t0 = time.time()
        opp_list = [list(s.opp) for s in generate_prime_shadows(11)]
        sh_path.write_text(json.dumps(opp_list))
        print(f"[Step 1] Generated {len(opp_list)} ch=11 prime shadows in {time.time() - t0:.2f}s")
    else:
        opp_list = json.loads(sh_path.read_text())

    shards_path = data_dir / "_ch11_shards.json"
    num_batches = 12
    batch_size = (len(opp_list) + num_batches - 1) // num_batches
    shard_state = json.loads(shards_path.read_text()) if shards_path.exists() else {"completed": {}, "stats": {}}

    for b_idx in range(num_batches):
        key = str(b_idx)
        if key in shard_state["completed"]:
            continue
        if time.time() - t_start > max_seconds:
            print(f"[PAUSE] Time budget reached ({len(shard_state['completed'])}/{num_batches} batches done).")
            return False
        tb = time.time()
        sub = [tuple(x) for x in opp_list[b_idx * batch_size : (b_idx + 1) * batch_size]]
        micro_sz = max(1, (len(sub) + 5) // 6)
        parts = [sub[i : i + micro_sz] for i in range(0, len(sub), micro_sz)]
        batch_diags: list[MarkedDiagram] = []
        raw_c, comp_c, disc_c = 0, 0, 0
        with ThreadPoolExecutor(max_workers=2) as ex:
            for d_list, st in ex.map(_run_c_census_shard, [(11, p, 2, 7, 384, False) for p in parts if p]):
                batch_diags.extend(d_list)
                raw_c += st[0]
                comp_c += st[1]
                disc_c += st[2]
        merged_b = merge_diagrams_c(11, batch_diags, max_bfs_states=512)
        shard_state["completed"][key] = [d.to_epd() for d in merged_b]
        shard_state["stats"][key] = [raw_c, comp_c, disc_c, len(batch_diags), len(merged_b)]
        shards_path.write_text(json.dumps(shard_state))
        print(
            f"[Step 2] Batch {b_idx + 1}/{num_batches}: raw={raw_c}, comps={comp_c}, "
            f"merged={len(merged_b)} in {time.time() - tb:.2f}s"
        )

    level_path = data_dir / "_level_kept_ch11.json"
    if not level_path.exists():
        if time.time() - t_start > max_seconds:
            return False
        all_shard_diags = [
            MarkedDiagram.from_epd(s)
            for k in sorted(shard_state["completed"], key=int)
            for s in shard_state["completed"][k]
        ]
        if len(all_shard_diags) > 12000:
            q_sz = (len(all_shard_diags) + 3) // 4
            tier1: list[MarkedDiagram] = []
            for qi in range(0, len(all_shard_diags), q_sz):
                tier1.extend(merge_diagrams_c(11, all_shard_diags[qi : qi + q_sz], max_bfs_states=1024))
            all_shard_diags = tier1
        if len(all_shard_diags) > 6000:
            mid = len(all_shard_diags) // 2
            all_shard_diags = merge_diagrams_c(11, all_shard_diags[:mid], max_bfs_states=1536) + merge_diagrams_c(
                11, all_shard_diags[mid:], max_bfs_states=1536
            )
        level_kept = merge_diagrams_c(11, all_shard_diags, max_bfs_states=2048)
        level_path.write_text(json.dumps([d.to_epd() for d in level_kept]))
        print(f"[Step 3] Global ch=11 level merge: {len(level_kept)} level-irreducible diagrams")

    def _run_parallel_chunks(
        src_epds: list[str],
        chunk_sz: int,
        part_path: Path,
        d1: int,
        d2: int,
        mbfs: int,
        step_tag: str,
    ) -> dict[str, list[str]] | None:
        state = json.loads(part_path.read_text()) if part_path.exists() else {}
        num_chunks = max(1, (len(src_epds) + chunk_sz - 1) // chunk_sz)
        pending = [ci for ci in range(num_chunks) if str(ci) not in state]
        for pi in range(0, len(pending), 2):
            if time.time() - t_start > max_seconds:
                return None
            pair = pending[pi : pi + 2]

            def _do_one(c_idx: int) -> tuple[int, int, list[str]]:
                sub_d = [MarkedDiagram.from_epd(s) for s in src_epds[c_idx * chunk_sz : (c_idx + 1) * chunk_sz]]
                tun_c = tunnel_diagrams_c(11, sub_d, depth_plus1=d1, depth_plus2=d2, max_bfs_states=mbfs)
                return c_idx, len(sub_d), [d.to_epd() for d in tun_c]

            with ThreadPoolExecutor(max_workers=len(pair)) as ex:
                results = list(ex.map(_do_one, pair))
            for c_idx, n_in, out_list in results:
                state[str(c_idx)] = out_list
                print(f"[{step_tag}] Chunk {c_idx + 1}/{num_chunks}: {n_in} -> {len(out_list)}")
            part_path.write_text(json.dumps(state))
        return state

    stage_a_path = data_dir / "_stageA_ch11.json"
    if not stage_a_path.exists():
        level_epds = json.loads(level_path.read_text())
        a_state = _run_parallel_chunks(level_epds, 1800, data_dir / "_stageA_chunks_ch11.json", 1, 0, 64, "Pass A")
        if a_state is None or time.time() - t_start > max_seconds:
            return False
        combined_a = [MarkedDiagram.from_epd(s) for k in sorted(a_state, key=int) for s in a_state[k]]
        if len(combined_a) > 3600:
            mid = len(combined_a) // 2
            with ThreadPoolExecutor(max_workers=2) as ex:
                halves = list(
                    ex.map(
                        lambda sd: tunnel_diagrams_c(11, sd, depth_plus1=1, depth_plus2=0, max_bfs_states=96),
                        [combined_a[:mid], combined_a[mid:]],
                    )
                )
            combined_a = halves[0] + halves[1]
        final_a = tunnel_diagrams_c(11, combined_a, depth_plus1=1, depth_plus2=0, max_bfs_states=128)
        stage_a_path.write_text(json.dumps([d.to_epd() for d in final_a]))

    stage_b_path = data_dir / "_stageB_ch11.json"
    if not stage_b_path.exists():
        a_epds = json.loads(stage_a_path.read_text())
        b_state = _run_parallel_chunks(a_epds, 450, data_dir / "_stageB_chunks_ch11.json", 2, 2, 256, "Pass B")
        if b_state is None or time.time() - t_start > max_seconds:
            return False
        combined_b = [MarkedDiagram.from_epd(s) for k in sorted(b_state, key=int) for s in b_state[k]]
        final_b = tunnel_diagrams_c(11, combined_b, depth_plus1=2, depth_plus2=2, max_bfs_states=384)
        stage_b_path.write_text(json.dumps([d.to_epd() for d in final_b]))

    stage_c1_path = data_dir / "_stageC1_ch11.json"
    if not stage_c1_path.exists():
        b_epds = json.loads(stage_b_path.read_text())
        c1_state = _run_parallel_chunks(b_epds, 120, data_dir / "_stageC1_chunks_ch11.json", 96, 48, 1536, "Pass C1")
        if c1_state is None or time.time() - t_start > max_seconds:
            return False
        combined_c1 = [MarkedDiagram.from_epd(s) for k in sorted(c1_state, key=int) for s in c1_state[k]]
        final_c1 = tunnel_diagrams_c(11, combined_c1, depth_plus1=96, depth_plus2=48, max_bfs_states=1536)
        stage_c1_path.write_text(json.dumps([d.to_epd() for d in final_c1]))

    c2_chunks = _run_parallel_chunks(
        json.loads(stage_c1_path.read_text()), 75, data_dir / "_stageC2_chunks_ch11.json", 256, 192, 2048, "Pass C2"
    )
    if c2_chunks is None or time.time() - t_start > max_seconds:
        return False
    c2_combined = [MarkedDiagram.from_epd(s) for k in sorted(c2_chunks, key=int) for s in c2_chunks[k]]
    tun_final = tunnel_diagrams_c(11, c2_combined, depth_plus1=256, depth_plus2=192, max_bfs_states=2048)
    (data_dir / "_census_ch11_all_tunneled.json").write_text(json.dumps([d.to_epd() for d in tun_final]))
    filtered = [d for d in tun_final if not has_certified_cyclic_unknotted_tube(d)]
    final_path.write_text(json.dumps([d.to_epd() for d in filtered]))
    print(f"[DONE] Saved {len(filtered)} final ch=11 survivors to {final_path.name}.")
    return True


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
    parser = argparse.ArgumentParser(description="Regenerate surfaceknots census datasets.")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--yoshikawa", action="store_true", help="Export yoshikawa_ch10.{json,csv}")
    parser.add_argument("--ch11-enumerate", action="store_true", help="Run checkpointed ch=11 enumeration/tunneling")
    parser.add_argument("--ch11-export", action="store_true", help="Stratify & export expanded_ch11.{json,csv}")
    parser.add_argument("--max-seconds", type=float, default=3600.0)
    args = parser.parse_args()

    if args.yoshikawa:
        export_yoshikawa_datasets(args.data_dir)
        print(f"Exported yoshikawa_ch10.{{json,csv}} to {args.data_dir}")
    if args.ch11_enumerate:
        run_ch11_enumeration(args.data_dir, max_seconds=args.max_seconds)
    if args.ch11_export:
        census_cache = args.data_dir / "_census_ch11.json"
        existing_json = args.data_dir / "expanded_ch11.json"
        if census_cache.exists():
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
            raise FileNotFoundError(f"Neither {census_cache} nor {existing_json} exists!")
        res = stratify_and_export_ch11(epds, args.data_dir, census_stats=c_stats)
        print(
            f"Exported expanded_ch11.{{json,csv}}: "
            f"{res['certified_new_prime_invariant_classes']} certified new prime classes, "
            f"{res['residual_colliding_diagrams']} residual collisions."
        )


if __name__ == "__main__":
    main()
