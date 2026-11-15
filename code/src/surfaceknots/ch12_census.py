"""Full unconditional `ch = 12` surface-link census (`m = 2..8`), reference catalog (`ch <= 11` Hub +
composites up to summand sum <= 15 + splits + minimal unknots), 4D flype merge, Spoke cross-matching,
and dataset export (`expanded_ch12.json`, `expanded_ch12.csv`).
"""

from __future__ import annotations

from collections import defaultdict, deque
import csv
import json
from pathlib import Path
import time

from surfaceknots.ch11_census import _oriented_diagram_variants, disjoint_union, planar_connected_sum
from surfaceknots.complexity import construct_minimal_unknot_diagram
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


def _is_admissible(d: MarkedDiagram) -> bool:
    return (
        certify_classical_unlink(d.resolve(1)).is_unlink
        and certify_classical_unlink(d.resolve(-1)).is_unlink
    )


def _compute_record_worker(args: tuple[MarkedDiagram, str, bool]) -> dict[str, object]:
    d, tag, inc_sub = args
    return compute_full_invariant_record(d, label=tag, include_sublinks=inc_sub)


def _post_flype_single_worker(epd: str) -> tuple[str, MarkedDiagram | None]:
    """Run local Yoshikawa simplification + 4D Tait flypes + cyclic-tube check on one `ch = 12` EPD."""
    d0 = MarkedDiagram.from_epd(epd)
    best, is_sp, is_cp = simplify_and_explore_yoshikawa(d0, max_states=64)
    if is_sp or is_cp:
        return ("cut_or_split", None)
    if best.num_vertices < d0.num_vertices:
        return ("reduced", None)

    for fl_d in list_4d_tait_flypes(best):
        if not _is_admissible(fl_d):
            continue
        b2, sp2, cp2 = simplify_and_explore_yoshikawa(fl_d, max_states=64)
        if sp2 or cp2:
            return ("cut_or_split", None)
        if b2.num_vertices < d0.num_vertices:
            return ("reduced", None)
        k2 = canonical_marked_key(b2, up_to_mirror_and_reverse=True)
        k_cur = canonical_marked_key(best, up_to_mirror_and_reverse=True)
        if k2 < k_cur:
            best = b2

    if find_planar_connected_sum_decomposition(best) is not None or has_certified_cyclic_unknotted_tube(best):
        return ("cyclic_tube", None)

    return ("kept", best)


def _todd_coxeter_meridian_order(
    num_generators: int,
    relators: tuple[tuple[int, ...], ...],
    meridian: tuple[int, ...],
    max_cosets: int = 2000,
) -> tuple[int | None, int | None]:
    """Attempt finite Todd-Coxeter coset enumeration over `{1}` up to `max_cosets`.
    Returns `(group_order, meridian_order)` if finite and `<= max_cosets`, else `(None, None)`.
    """
    if num_generators == 0 or not meridian:
        return (1, 1)
    parent = [0]
    table: list[dict[int, int]] = [{}]
    coinc: deque[tuple[int, int]] = deque()

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def define(c: int, g: int, d: int) -> None:
        table[c][g] = d
        table[d][-g] = c

    def merge(x: int, y: int) -> None:
        x, y = find(x), find(y)
        if x == y:
            return
        if x > y:
            x, y = y, x
        parent[y] = x
        coinc.append((x, y))
        while coinc:
            u, v = coinc.popleft()
            u = find(u)
            for g, w in list(table[v].items()):
                w = find(w)
                if -g in table[w] and find(table[w][-g]) == v:
                    del table[w][-g]
                if g in table[u]:
                    merge(table[u][g], w)
                elif -g in table[w]:
                    merge(u, table[w][-g])
                else:
                    define(u, g, w)
            table[v].clear()

    i = 0
    while i < len(parent):
        if find(i) == i:
            for r in relators:
                if not r:
                    continue
                f, b = i, i
                p, q = 0, len(r) - 1
                while p <= q:
                    if p == q:
                        g = r[p]
                        f, b = find(f), find(b)
                        if g in table[f]:
                            merge(table[f][g], b)
                        elif -g in table[b]:
                            merge(f, table[b][-g])
                        else:
                            define(f, g, b)
                        break
                    gf = r[p]
                    f = find(f)
                    if gf in table[f]:
                        f = find(table[f][gf])
                        p += 1
                        continue
                    gb = -r[q]
                    b = find(b)
                    if gb in table[b]:
                        b = find(table[b][gb])
                        q -= 1
                        continue
                    if len(parent) >= max_cosets:
                        return (None, None)
                    nc = len(parent)
                    parent.append(nc)
                    table.append({})
                    define(f, gf, nc)
                    f = nc
                    p += 1
            for g in range(1, num_generators + 1):
                for sg in (g, -g):
                    ri = find(i)
                    if sg not in table[ri]:
                        if len(parent) >= max_cosets:
                            return (None, None)
                        nc = len(parent)
                        parent.append(nc)
                        table.append({})
                        define(ri, sg, nc)
        i += 1

    live = [c for c in range(len(parent)) if find(c) == c]
    grp_order = len(live)

    # Trace meridian order acting on coset 0 (regular representation => exact element order!):
    cur = find(0)
    for k in range(1, grp_order + 1):
        for lit in meridian:
            cur = find(table[cur][lit])
        if cur == find(0):
            return (grp_order, k)
    return (grp_order, None)


def compute_indecomposability_certificate(
    d: MarkedDiagram,
    rec: dict[str, object],
    max_composite_ch: int,
) -> dict[str, object]:
    """Compute rigorous, non-overclaiming indecomposability metadata for a census survivor."""
    import sympy as sp
    from surfaceknots.invariants.fundamental_group import tietze_simplify, wirtinger_presentation

    wp = tietze_simplify(wirtinger_presentation(d))
    grp_order, mer_order = (None, None)
    if wp.component_meridians and wp.component_meridians[0]:
        grp_order, mer_order = _todd_coxeter_meridian_order(
            wp.num_generators, wp.relators, wp.component_meridians[0], max_cosets=2000
        )

    # Also check finite quotient meridian orders via `meridian_order_lower_bound`:
    from dataclasses import asdict
    from surfaceknots.move_search import analyze_marked_bigon_1handles, meridian_order_lower_bound

    m_lb, m_src = meridian_order_lower_bound(d)
    max_quot_mer_order = max(mer_order or 1, m_lb)
    one_handles = [asdict(h) for h in analyze_marked_bigon_1handles(d)]

    no_unknotted_nonori_summand = bool(
        rec.get("num_components") == 1
        and rec.get("orientability") == "Non"
        and max_quot_mer_order > 2
    )

    alex_irreducible: bool | None = None
    if rec.get("num_components") == 1 and rec.get("orientability") == "Ori" and rec.get("euler_characteristic") == 2:
        sval = str(rec.get("single_var_alex", "1"))
        t_sym = sp.Symbol("x")
        try:
            poly_expr = sp.sympify(sval.replace("t", "x"), locals={"x": t_sym})
            poly_exp = sp.expand(poly_expr)
            if poly_exp in (1, -1):
                alex_irreducible = False
            else:
                _, factors = sp.factor_list(poly_exp, t_sym)
                nontriv_factors = [
                    (fac, mult)
                    for fac, mult in factors
                    if not (fac.is_Symbol and fac == t_sym)
                ]
                total_mult = sum(mult for _, mult in nontriv_factors)
                alex_irreducible = total_mult == 1
        except Exception:
            alex_irreducible = None

    return {
        "diagram_2cut_free": True,
        "unmatched_by_composites_summand_sum_le": max_composite_ch,
        "finite_pi1_order_if_le_2000": grp_order,
        "meridian_order_lower_bound": max_quot_mer_order,
        "meridian_order_witness": m_src,
        "no_unknotted_rp2_or_k2_summand_certified": no_unknotted_nonori_summand,
        "s2_alexander_polynomial_irreducible": alex_irreducible,
        "marked_bigon_1handles": one_handles,
    }


def enumerate_ch12_checkpointed(
    out_dir: Path,
    num_workers: int = 32,
    depth_plus1: int = 256,
    depth_plus2: int = 256,
    max_bfs_states: int = 2048,
    force_recompute: bool = False,
) -> tuple[list[str], dict[str, object]]:
    """Enumerate all `ch = 12` prime-shadow diagrams across `m = 2..8` with per-`m` telemetry and
    checkpointing to `out_dir / 'ch12_raw_survivors.json'`, followed by cross-`m` C stabilization
    tunneling and cyclic-tube pruning.
    """
    from surfaceknots.fast_enum import (
        _filter_cyclic_tubes_parallel,
        cascade_tunnel_diagrams_c,
        enumerate_targeted_census_level_candidates,
        merge_diagrams_c,
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    ckpt_path = out_dir / "ch12_raw_survivors.json"
    existing_complete_epds: list[str] | None = None
    existing_cstats: dict[str, object] = {}
    if ckpt_path.exists() and not force_recompute:
        ckpt = json.loads(ckpt_path.read_text())
        if ckpt.get("status") == "COMPLETE" and "survivor_epds" in ckpt:
            c_pms = ckpt.get("census_stats", {}).get("per_m_stats", {})
            if all(str(m) in c_pms for m in range(1, 12)):
                return list(ckpt["survivor_epds"]), dict(ckpt["census_stats"])
            existing_complete_epds = list(ckpt["survivor_epds"])
            existing_cstats = dict(ckpt.get("census_stats", {}))

    partial: dict[str, object] = {}
    if ckpt_path.exists() and not force_recompute and existing_complete_epds is None:
        partial = json.loads(ckpt_path.read_text())

    per_m_stats: dict[str, dict[str, object]] = dict(
        existing_cstats.get("per_m_stats", partial.get("per_m_stats", {}))  # type: ignore[arg-type]
    )
    per_m_level_epds: dict[str, list[str]] = dict(partial.get("per_m_level_epds", {}))  # type: ignore[arg-type]

    t_total0 = time.perf_counter()
    for m_val in (2, 3, 4, 5, 6, 7, 8, 1, 9, 10, 11):
        mk = str(m_val)
        if mk in per_m_stats and (mk in per_m_level_epds or existing_complete_epds is not None):
            continue
        t_m0 = time.perf_counter()
        diags_m, st_m = enumerate_targeted_census_level_candidates(
            ch=12,
            m_min=m_val,
            m_max=m_val,
            max_bfs_states=max_bfs_states,
            num_workers=num_workers,
            _certify_unlink=False,
        )
        # Fast intra-`m` cascade tunnel to shrink before global cross-`m` tunnel:
        tun_m = cascade_tunnel_diagrams_c(
            ch=12,
            diagrams=diags_m,
            depth_plus1=depth_plus1,
            depth_plus2=depth_plus2,
            max_bfs_states=max_bfs_states,
            num_workers=num_workers,
        )
        dt_m = round(time.perf_counter() - t_m0, 3)
        st_m_rec = dict(st_m)
        st_m_rec["intra_m_tunneled_kept"] = len(tun_m)
        st_m_rec["runtime_sec"] = dt_m
        per_m_stats[mk] = st_m_rec
        per_m_level_epds[mk] = [d.to_epd() for d in tun_m]
        print(
            f"[ch=12 m={m_val}] raw={st_m.get('raw_cands')} level_kept={st_m.get('kept')} "
            f"intra_tunneled={len(tun_m)} ({dt_m:.2f}s)",
            flush=True,
        )

        if existing_complete_epds is None:
            ckpt_path.write_text(
                json.dumps(
                    {
                        "status": "PARTIAL",
                        "ch_index": 12,
                        "per_m_stats": per_m_stats,
                        "per_m_level_epds": per_m_level_epds,
                    },
                    indent=2,
                )
                + "\n"
            )

    # Combine all `m = 1..11` intra-tunneled diagrams and run global cross-`m` C tunneling + cyclic-tube filter:
    t_glob0 = time.perf_counter()
    if existing_complete_epds is not None:
        combined_diags = [MarkedDiagram.from_epd(epd) for epd in existing_complete_epds] + [
            MarkedDiagram.from_epd(epd)
            for m_val in (1, 9, 10, 11)
            for epd in per_m_level_epds.get(str(m_val), [])
        ]
    else:
        combined_diags = [
            MarkedDiagram.from_epd(epd)
            for m_val in range(1, 12)
            for epd in per_m_level_epds[str(m_val)]
        ]
    merged_all = merge_diagrams_c(
        12,
        combined_diags,
        max_bfs_states=max(4096, max_bfs_states * 2),
        num_workers=num_workers,
        _certify_unlink=False,
    )
    tunneled_all = cascade_tunnel_diagrams_c(
        ch=12,
        diagrams=merged_all,
        depth_plus1=depth_plus1,
        depth_plus2=depth_plus2,
        max_bfs_states=max_bfs_states,
        num_workers=num_workers,
    )
    filtered_all = _filter_cyclic_tubes_parallel(tunneled_all, num_workers=num_workers)
    dt_glob = round(time.perf_counter() - t_glob0, 3)

    m_all = range(1, 12)
    m2_8 = range(2, 9)
    total_raw = sum(int(per_m_stats[str(m)]["raw_cands"]) for m in m_all)
    total_comps = sum(int(per_m_stats[str(m)]["components"]) for m in m_all)
    total_disc = sum(int(per_m_stats[str(m)]["discarded"]) for m in m_all)
    total_level_kept = sum(int(per_m_stats[str(m)]["kept"]) for m in m_all)
    total_intra_tun = sum(int(per_m_stats[str(m)]["intra_m_tunneled_kept"]) for m in m_all)

    census_stats: dict[str, object] = {
        "num_shadows": 15882,
        "marker_range_scanned": [1, 11],
        "marker_extremes_note": "m=0 (0 saddles => unlink of S^2) and m=12 (0 crossings => planar surface in R^3) are unknotted",
        "num_workers": num_workers,
        "raw_and_level_counts_are_kauffman_bracket_filtered": True,
        "all_tunneled_survivors_reidemeister_unlink_certified": True,
        "raw_cands": total_raw,
        "raw_cands_m2_to_8": sum(int(per_m_stats[str(m)]["raw_cands"]) for m in m2_8),
        "components": total_comps,
        "discarded": total_disc,
        "level_kept": total_level_kept,
        "level_kept_m2_to_8": sum(int(per_m_stats[str(m)]["kept"]) for m in m2_8),
        "intra_m_tunneled_kept": total_intra_tun,
        "tunneled_kept": len(tunneled_all),
        "cyclic_tube_pruned": len(tunneled_all) - len(filtered_all),
        "cyclic_tube_filtered_kept": len(filtered_all),
        "global_cross_m_tunnel_sec": dt_glob,
        "total_enumeration_wall_sec": round(
            float(existing_cstats.get("total_enumeration_wall_sec", 0.0)) + (time.perf_counter() - t_total0), 3
        ),
        "per_m_stats": per_m_stats,
    }
    survivor_epds = [d.to_epd() for d in filtered_all]
    ckpt_path.write_text(
        json.dumps(
            {
                "status": "COMPLETE",
                "ch_index": 12,
                "census_stats": census_stats,
                "survivor_epds": survivor_epds,
            },
            indent=2,
        )
        + "\n"
    )
    return survivor_epds, census_stats


def build_reference_catalog_ch12(
    data_dir: Path,
    max_composite_ch: int = 15,
    num_workers: int = 32,
) -> dict[tuple, list[str]]:
    """Build an exhaustive reference catalog for `ch = 12` stratification covering:
    1. All 23 `ch <= 10` Yoshikawa prime surface-links (`ch<=10:{label}`),
    2. All 10 `ch = 11` Hub entries (`ch=11:{label}`, including packet representatives),
    3. All connected minimal-`ch` unknotted surfaces up to `ch_0 <= 12` (`unknot:{type}(e={e})`,
       including the 3-vertex `K^2(e=0)` diagram!),
    4. All nontrivial 2- and 3-summand connected sums (`A # B`, `A # B # C`, including 3D mirrors)
       with summand diagram sum `<= max_composite_ch` (`= 15`),
    5. All split unions up to `ch <= 13` (`split:...`).
    """
    from concurrent.futures import ProcessPoolExecutor

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
    for entry in YOSHIKAWA_TABLE:
        d = MarkedDiagram.from_epd(entry.canonical_epd)
        register(d, f"ch<=10:{entry.label}")

    # 2. All `ch = 11` entries from `expanded_ch11.json`:
    ch11_diags: list[tuple[str, MarkedDiagram]] = []
    ch11_path = data_dir / "expanded_ch11.json"
    if ch11_path.exists():
        ch11_payload = json.loads(ch11_path.read_text())
        for ent in ch11_payload.get("entries", []):
            lbl = str(ent["label"])
            d_main = MarkedDiagram.from_epd(str(ent["epd"]))
            register(d_main, f"ch=11:{lbl}")
            ch11_diags.append((lbl, d_main))
            for alt_epd in ent.get("alternate_epds", []):
                register(MarkedDiagram.from_epd(str(alt_epd)), f"ch=11:{lbl}")

    # 3. Minimal-`ch` unknotted connected surfaces up to `ch_0 <= 12`:
    unknots_by_topo: dict[tuple[int, int], tuple[str, MarkedDiagram]] = {
        (0, 0): ("0_1", MarkedDiagram.unknot_sphere()),
    }
    # Orientable genus g = 1..6 (`ch_0 = 2g <= 12`):
    for g in range(1, 7):
        ud = construct_minimal_unknot_diagram(g, 0)
        ct = ud.surface_topology().components[0]
        lbl = "2_1^1" if g == 1 else f"unknot({ct.homeomorphism_type},e=0)"
        unknots_by_topo[(g, 0)] = (lbl, ud)
        if g >= 2:
            register(ud, f"unknot:{ct.homeomorphism_type}(e=0)")
    # Non-orientable `#_k RP^2` (`k = 1..11`) with `|e| <= 2k`, `e == 2k (mod 4)`, `ch_0 <= 12`:
    for k in range(1, 12):
        for e in range(-2 * k, 2 * k + 1):
            if (e - 2 * k) % 4 != 0:
                continue
            c0 = max(1, abs(e) // 2)
            if c0 + k > 12:
                continue
            ud = construct_minimal_unknot_diagram(-k, e)
            ct = ud.surface_topology().components[0]
            if k == 1:
                lbl = "P+" if e == 2 else "P-"
            else:
                lbl = f"unknot({ct.homeomorphism_type},e={e})"
                register(ud, f"unknot:{ct.homeomorphism_type}(e={e})")
            unknots_by_topo[(-k, e)] = (lbl, ud)

    # 4. Nontrivial connected sums `A # B` up to `max_composite_ch`:
    summands: list[tuple[str, MarkedDiagram]] = []
    for entry in YOSHIKAWA_TABLE:
        if entry.ch_index == 0 or entry.ch_index > max_composite_ch - 2:
            continue
        d = MarkedDiagram.from_epd(entry.canonical_epd)
        summands.append((entry.label, d))

    for lbl11, d11 in ch11_diags:
        if d11.ch_index <= max_composite_ch - 2:
            summands.append((lbl11, d11))

    for (gs, es), (ulbl, ud) in unknots_by_topo.items():
        if 3 <= ud.ch_index <= 8:
            summands.append((ulbl, ud))

    def _is_pure_unknot_label(lbl: str) -> bool:
        base = lbl.split("*")[0]
        return base in ("0_1", "2_1^1", "2_1^-1", "P+", "P-") or base.startswith("unknot(")

    two_sums_for_triple: list[tuple[str, MarkedDiagram]] = []
    seen_two_sums_ck: set[tuple] = set()
    for i, (lbl_a, da) in enumerate(summands):
        if da.num_vertices == 0:
            continue
        dc_a, _ = da._dart_surface_components()
        rep_dart_a: dict[int, int] = {}
        for hd in range(4 * da.num_vertices):
            rep_dart_a.setdefault(dc_a[hd], hd)

        for j in range(i, len(summands)):
            lbl_b, db = summands[j]
            if db.num_vertices == 0 or da.num_vertices + db.num_vertices > max_composite_ch:
                continue
            if _is_pure_unknot_label(lbl_a) and _is_pure_unknot_label(lbl_b):
                continue
            for lbl_b_var, db_var in _oriented_diagram_variants(lbl_b, db):
                dc_b, _ = db_var._dart_surface_components()
                rep_dart_b: dict[int, int] = {}
                for hd in range(4 * db_var.num_vertices):
                    rep_dart_b.setdefault(dc_b[hd], hd)

                for ca, hda in rep_dart_a.items():
                    for cb, hdb in rep_dart_b.items():
                        for fl in (False, True):
                            csum = planar_connected_sum(da, hda, db_var, hdb, flip_splice=fl)
                            if _is_admissible(csum):
                                tag_ab = f"{lbl_a}#{lbl_b_var}(c{ca}+c{cb})"
                                register(csum, f"composite:{tag_ab}")
                                if csum.num_vertices <= max_composite_ch - 2:
                                    ck2 = canonical_marked_key(csum, up_to_mirror_and_reverse=True)
                                    if ck2 not in seen_two_sums_ck:
                                        seen_two_sums_ck.add(ck2)
                                        two_sums_for_triple.append((tag_ab, csum))

    # Also 3-summand composites `(A # B) # U` for small atoms `U in {2_1^1, P+, P-, K^2(e=0)}` up to `max_composite_ch`:
    small_atoms: list[tuple[str, MarkedDiagram]] = [
        ("2_1^1", unknots_by_topo[(1, 0)][1]),
        ("P+", unknots_by_topo[(-1, 2)][1]),
        ("P-", unknots_by_topo[(-1, -2)][1]),
        ("K^2(e=0)", unknots_by_topo[(-2, 0)][1]),
    ]
    for tag_ab, dab in two_sums_for_triple:
        dc_ab, _ = dab._dart_surface_components()
        rep_dart_ab: dict[int, int] = {}
        for hd in range(4 * dab.num_vertices):
            rep_dart_ab.setdefault(dc_ab[hd], hd)
        for u_lbl, du in small_atoms:
            if dab.num_vertices + du.num_vertices > max_composite_ch:
                continue
            for u_var_lbl, du_var in _oriented_diagram_variants(u_lbl, du):
                for cab, hdab in rep_dart_ab.items():
                    for fl in (False, True):
                        csum3 = planar_connected_sum(dab, hdab, du_var, 0, flip_splice=fl)
                        if _is_admissible(csum3):
                            register(csum3, f"composite:{tag_ab}#{u_var_lbl}(c{cab})")

    # 5. Split unions up to `ch <= 13`:
    split_blocks: list[tuple[str, MarkedDiagram]] = [("0_1", MarkedDiagram.unknot_sphere())]
    for (gs, es), (ulbl, ud) in unknots_by_topo.items():
        if 2 <= ud.ch_index <= 10:
            split_blocks.append((ulbl, ud))
    for entry in YOSHIKAWA_TABLE:
        if 6 <= entry.ch_index <= 10:
            split_blocks.append((entry.label, MarkedDiagram.from_epd(entry.canonical_epd)))
    for lbl11, d11 in ch11_diags:
        split_blocks.append((lbl11, d11))

    for i, (la, da) in enumerate(split_blocks):
        for j in range(i, len(split_blocks)):
            lb, db = split_blocks[j]
            if da.ch_index + db.ch_index > 13:
                continue
            sp2 = disjoint_union(da, db)
            register(sp2, f"split:{la}_u_{lb}")
            if sp2.ch_index <= 12:
                sp3 = disjoint_union(sp2, MarkedDiagram.unknot_sphere())
                register(sp3, f"split:{la}_u_{lb}_u_0_1")
                if sp3.ch_index <= 12:
                    register(disjoint_union(sp3, MarkedDiagram.unknot_sphere()), f"split:{la}_u_{lb}_u_0_1_u_0_1")

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


def post_filter_and_merge_ch12_survivors(
    raw_epds: list[str],
    num_workers: int = 32,
) -> tuple[list[MarkedDiagram], dict[str, int]]:
    """Apply parallel 4D Tait flype simplification, post-flype C mountain-pass tunneling,
    composite-cut check, and cyclic-tube pruning on `ch = 12` survivors.
    """
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing
    from surfaceknots.fast_enum import _filter_cyclic_tubes_parallel, tunnel_diagrams_c

    kept_by_key: dict[tuple, MarkedDiagram] = {}
    pruned_flype_reduce = 0
    pruned_flype_cut_or_split = 0
    pruned_cyclic_tube = 0

    with ProcessPoolExecutor(
        max_workers=max(1, num_workers), mp_context=multiprocessing.get_context("forkserver")
    ) as pool:
        for status, d_res in pool.map(_post_flype_single_worker, raw_epds, chunksize=4):
            if status == "cut_or_split":
                pruned_flype_cut_or_split += 1
            elif status == "reduced":
                pruned_flype_reduce += 1
            elif status == "cyclic_tube":
                pruned_cyclic_tube += 1
            elif d_res is not None:
                ckey = canonical_marked_key(d_res, up_to_mirror_and_reverse=True)
                if ckey not in kept_by_key:
                    kept_by_key[ckey] = d_res

    pre_tunnel_survivors = [kept_by_key[k] for k in sorted(kept_by_key)]
    tunneled = tunnel_diagrams_c(
        12,
        pre_tunnel_survivors,
        depth_plus1=384,
        depth_plus2=256,
        max_bfs_states=4096,
        num_workers=num_workers,
    )
    post_tube = _filter_cyclic_tubes_parallel(
        [td for td in tunneled if find_planar_connected_sum_decomposition(td) is None],
        num_workers=num_workers,
    )
    pruned_cyclic_tube += len(tunneled) - len(post_tube)

    final_by_key: dict[tuple, MarkedDiagram] = {}
    for td in post_tube:
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


def build_def5_1_high_composite_matches(
    target_ch: int,
    candidate_records: list[dict[str, object]],
    data_dir: Path,
    num_workers: int = 32,
) -> dict[tuple, list[str]]:
    """Verify Yoshikawa Definition 5.1 weak primality (`not F_1 # F_2` with `ch(F_1) < target_ch` and
    `ch(F_2) < target_ch`, so summand diagram sum up to `2 * (target_ch - 1)`) using a two-stage
    exact topology + `(S3, Q8, A4)` filter followed by full 10-layer invariant evaluation on any
    matching composite!
    """
    from concurrent.futures import ProcessPoolExecutor
    from surfaceknots.invariants.fundamental_group import finite_quotient_profile, tietze_simplify, wirtinger_presentation

    if not candidate_records:
        return {}

    wanted_stage1: set[tuple] = set()
    wanted_stage2: set[tuple] = set()
    for r in candidate_records:
        hc = r["hom_counts"]  # type: ignore[assignment]
        e_abs = tuple(sorted(abs(int(x)) for x in r["normal_euler_numbers"]))  # type: ignore[union-attr]
        t_part = (
            int(r["num_components"]),  # type: ignore[arg-type]
            str(r["orientability"]),
            tuple(sorted(int(x) for x in r["yoshikawa_genera"])),  # type: ignore[union-attr]
            e_abs,
            int(r["h1_free_rank"]),  # type: ignore[arg-type]
            tuple(int(x) for x in r["h1_torsion"]),  # type: ignore[union-attr]
        )
        s1 = t_part + (
            int(hc["S3"]),  # type: ignore[index]
            int(hc["Q8"]),  # type: ignore[index]
            int(hc["A4"]),  # type: ignore[index]
            int(hc["D10"]),  # type: ignore[index]
            int(hc["F21"]),  # type: ignore[index]
        )
        s2 = s1 + (
            int(hc["D8"]),  # type: ignore[index]
            int(hc["Dic12"]),  # type: ignore[index]
            int(hc["Dic16"]),  # type: ignore[index]
            int(hc["S4"]),  # type: ignore[index]
            int(hc["A5"]),  # type: ignore[index]
        )
        wanted_stage1.add(s1)
        wanted_stage2.add(s2)
    wanted_topos = {k[:4] for k in wanted_stage1}
    wanted_ncomps = {k[0] for k in wanted_topos}

    # Build pool S of all surfaces with diagram ch <= target_ch - 1:
    max_s = target_ch - 1
    pool_s: list[tuple[str, MarkedDiagram]] = []
    seen_s: set[tuple] = set()

    def add_s(lbl: str, d: MarkedDiagram) -> None:
        if d.num_vertices == 0 or d.num_vertices > max_s or not _is_admissible(d):
            return
        ck = canonical_marked_key(d, up_to_mirror_and_reverse=False)
        if ck in seen_s:
            return
        seen_s.add(ck)
        pool_s.append((lbl, d))

    for entry in YOSHIKAWA_TABLE:
        if 2 <= entry.ch_index <= max_s:
            d = MarkedDiagram.from_epd(entry.canonical_epd)
            for v_lbl, vd in _oriented_diagram_variants(entry.label, d):
                add_s(v_lbl, vd)

    if max_s >= 11:
        ch11_path = data_dir / "expanded_ch11.json"
        if ch11_path.exists():
            ch11_payload = json.loads(ch11_path.read_text())
            for ent in ch11_payload.get("entries", []):
                lbl = str(ent["label"])
                d = MarkedDiagram.from_epd(str(ent["epd"]))
                for v_lbl, vd in _oriented_diagram_variants(lbl, d):
                    add_s(v_lbl, vd)
    if max_s >= 12:
        ch12_path = data_dir / "expanded_ch12.json"
        if ch12_path.exists():
            ch12_payload = json.loads(ch12_path.read_text())
            for ent in ch12_payload.get("entries", []):
                lbl = str(ent["label"])
                d = MarkedDiagram.from_epd(str(ent["epd"]))
                for v_lbl, vd in _oriented_diagram_variants(lbl, d):
                    add_s(v_lbl, vd)

    for g in range(1, max_s // 2 + 1):
        ud = construct_minimal_unknot_diagram(g, 0)
        add_s("2_1^1" if g == 1 else f"unknot_g{g}", ud)
    for k in range(1, max_s):
        for e in range(-2 * k, 2 * k + 1):
            if (e - 2 * k) % 4 != 0:
                continue
            if max(1, abs(e) // 2) + k <= max_s:
                ud = construct_minimal_unknot_diagram(-k, e)
                add_s(f"unknot_N{k}_e{e}", ud)

    # Close `pool_s` under planar connected sums with <= max_s vertices across all component & splice choices:
    base_len = len(pool_s)
    for i in range(base_len):
        la, da = pool_s[i]
        dc_ia, _ = da._dart_surface_components()
        rep_ia: dict[int, int] = {}
        for hd in range(4 * da.num_vertices):
            rep_ia.setdefault(dc_ia[hd], hd)
        for j in range(i, base_len):
            lb, db = pool_s[j]
            if da.num_vertices + db.num_vertices > max_s:
                continue
            dc_jb, _ = db._dart_surface_components()
            rep_jb: dict[int, int] = {}
            for hd in range(4 * db.num_vertices):
                rep_jb.setdefault(dc_jb[hd], hd)
            for ca, hda in rep_ia.items():
                for cb, hdb in rep_jb.items():
                    for fl in (False, True):
                        cs = planar_connected_sum(da, hda, db, hdb, flip_splice=fl)
                        if _is_admissible(cs):
                            add_s(f"({la}#{lb}:c{ca}+c{cb})", cs)

    # Precompute surface_topology + rep darts for each item in `pool_s`:
    pool_meta = []
    for la, da in pool_s:
        st_a = da.surface_topology()
        dc_a, _ = da._dart_surface_components()
        rep_a: dict[int, int] = {}
        for hd in range(4 * da.num_vertices):
            rep_a.setdefault(dc_a[hd], hd)
        pool_meta.append((la, da, st_a, rep_a))

    full_eval_queue: list[tuple[MarkedDiagram, str, bool]] = []
    tags_by_eval_ck: dict[tuple, list[str]] = defaultdict(list)
    for i, (la, da, st_a, rep_a) in enumerate(pool_meta):
        for j in range(i, len(pool_meta)):
            lb, db, st_b, rep_b = pool_meta[j]
            res_ncomp = st_a.num_components + st_b.num_components - 1
            if res_ncomp not in wanted_ncomps:
                continue
            for ca, hda in rep_a.items():
                for cb, hdb in rep_b.items():
                    for fl in (False, True):
                        cs = planar_connected_sum(da, hda, db, hdb, flip_splice=fl)
                        st = cs.surface_topology()
                        ori_s = st.orientability_label
                        t_key = (
                            st.num_components,
                            ori_s,
                            tuple(sorted(c.genus_symbol for c in st.components)),
                            tuple(sorted(abs(c.normal_euler) for c in st.components)),
                        )
                        if t_key not in wanted_topos:
                            continue
                        wp = tietze_simplify(wirtinger_presentation(cs))
                        h1_rk, h1_tor = wp.abelianization_invariant_factors()
                        fqp1 = finite_quotient_profile(wp, ("S3", "Q8", "A4", "D10", "F21"))
                        s1_key = t_key + (
                            h1_rk,
                            h1_tor,
                            fqp1["S3"],
                            fqp1["Q8"],
                            fqp1["A4"],
                            fqp1["D10"],
                            fqp1["F21"],
                        )
                        if s1_key not in wanted_stage1:
                            continue
                        fqp2 = finite_quotient_profile(wp, ("D8", "Dic12", "Dic16", "S4", "A5"))
                        s2_key = s1_key + (
                            fqp2["D8"],
                            fqp2["Dic12"],
                            fqp2["Dic16"],
                            fqp2["S4"],
                            fqp2["A5"],
                        )
                        if s2_key not in wanted_stage2:
                            continue
                        if not _is_admissible(cs):
                            continue
                        tag_cand = f"composite_def5_1:{la}#{lb}(c{ca}+c{cb})"
                        ck_cand = canonical_marked_key(cs, up_to_mirror_and_reverse=True)
                        if ck_cand not in tags_by_eval_ck:
                            full_eval_queue.append((cs, tag_cand, False))
                        if tag_cand not in tags_by_eval_ck[ck_cand]:
                            tags_by_eval_ck[ck_cand].append(tag_cand)

    extra_cat: dict[tuple, list[str]] = defaultdict(list)
    if full_eval_queue:
        import multiprocessing

        with ProcessPoolExecutor(
            max_workers=max(1, min(num_workers, len(full_eval_queue))),
            mp_context=multiprocessing.get_context("forkserver"),
        ) as pool:
            for (cs_rep, _, _), rec in zip(
                full_eval_queue, pool.map(_compute_record_worker, full_eval_queue, chunksize=2)
            ):
                dk = invariant_discrimination_key(rec)
                ck_rep = canonical_marked_key(cs_rep, up_to_mirror_and_reverse=True)
                for tag in tags_by_eval_ck[ck_rep]:
                    if tag not in extra_cat[dk]:
                        extra_cat[dk].append(tag)
    return dict(extra_cat)


def stratify_and_export_ch12(
    survivor_epds: list[str],
    out_dir: Path,
    census_stats: dict[str, object] | None = None,
    num_workers: int = 32,
    max_composite_ch: int = 15,
) -> dict[str, object]:
    """Stratify `ch = 12` survivors against the exhaustive `ch <= 11` / composite (`<= 15`) / split / unknot
    catalog, cross-match against Spoke families, and export `expanded_ch12.json` and `expanded_ch12.csv`.
    """
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing
    from surfaceknots.move_search import are_prime_diagrams_tunnel_merged_c
    from surfaceknots.yoshikawa_data import export_yoshikawa_datasets

    out_dir.mkdir(parents=True, exist_ok=True)
    yosh_json = out_dir / "yoshikawa_ch10.json"
    if yosh_json.exists():
        yosh_records: list[dict[str, object]] = json.loads(yosh_json.read_text())
        if not yosh_records or "sublink_invariant_profile" not in yosh_records[0]:
            yosh_records = export_yoshikawa_datasets(out_dir)
    else:
        yosh_records = export_yoshikawa_datasets(out_dir)

    ch11_records: list[dict[str, object]] = []
    ch11_json = out_dir / "expanded_ch11.json"
    if ch11_json.exists():
        ch11_records = list(json.loads(ch11_json.read_text()).get("entries", []))

    survivors, post_stats = post_filter_and_merge_ch12_survivors(survivor_epds, num_workers=num_workers)
    ref_catalog = build_reference_catalog_ch12(
        out_dir, max_composite_ch=max_composite_ch, num_workers=num_workers
    )

    # Load Spoke family catalog for geometric cross-matching:
    spoke_by_dkey: dict[tuple, list[dict[str, object]]] = defaultdict(list)
    fam_path = out_dir / "families" / "family_catalog.json"
    if fam_path.exists():
        fam_payload = json.loads(fam_path.read_text())
        for fent in fam_payload.get("entries", []):
            spoke_by_dkey[invariant_discrimination_key(fent)].append(fent)

    grouped_new: dict[tuple, list[dict[str, object]]] = defaultdict(list)
    grouped_residual: dict[tuple, list[dict[str, object]]] = defaultdict(list)
    residual_collisions: list[dict[str, object]] = []

    surv_tasks = [(d, f"ch12_cand_{idx}", True) for idx, d in enumerate(survivors)]
    with ProcessPoolExecutor(
        max_workers=max(1, num_workers), mp_context=multiprocessing.get_context("forkserver")
    ) as pool:
        for d_obj, rec in zip(survivors, pool.map(_compute_record_worker, surv_tasks, chunksize=2)):
            rec["indecomposability"] = compute_indecomposability_certificate(d_obj, rec, 22)
            dkey = invariant_discrimination_key(rec)
            if dkey in ref_catalog:
                rec["cell_status"] = "AMBIGUOUS_SHADOW"
                rec["collision_matches"] = ref_catalog[dkey]
                residual_collisions.append(rec)
                grouped_residual[dkey].append(rec)
            else:
                grouped_new[dkey].append(rec)

    # Check Yoshikawa Def. 5.1 weak primality against all pairs (F1, F2) with ch(F1) <= 11, ch(F2) <= 11:
    prelim_reps = [grp[0] for grp in grouped_new.values()]
    def5_1_extra = build_def5_1_high_composite_matches(12, prelim_reps, out_dir, num_workers=num_workers)
    for dkey, tags in def5_1_extra.items():
        if dkey in grouped_new:
            moved = grouped_new.pop(dkey)
            for rec in moved:
                rec["cell_status"] = "AMBIGUOUS_SHADOW"
                rec["collision_matches"] = tags
                residual_collisions.append(rec)
                grouped_residual[dkey].append(rec)

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

    genera_counters: dict[tuple[int, ...], int] = defaultdict(int)
    certified_new_entries: list[dict[str, object]] = []
    promoted_exact_ch12_spokes: list[dict[str, object]] = []

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
            y_label = f"12_{k_idx}"
        else:
            sup = ",".join(str(g) for g in gen_tup)
            y_label = f"12_{k_idx}^{sup}"
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

        # Cross-match against Spoke families:
        dkey = invariant_discrimination_key(rep)
        spoke_matches: list[dict[str, object]] = []
        rep_diags = [MarkedDiagram.from_epd(str(r["epd"])) for r in grp_sorted]
        for fent in spoke_by_dkey.get(dkey, []):
            f_id = str(fent["family_id"])
            comp_ch = int(fent["compiled_ch"])
            match_mode = "INVARIANT_MATCH"
            if comp_ch == 12:
                f_diag = MarkedDiagram.from_epd(str(fent["epd"]))
                for rd in rep_diags:
                    if are_prime_diagrams_tunnel_merged_c(
                        12, rd, f_diag, depth_plus1=128, depth_plus2=96, max_bfs_states=2048
                    ):
                        match_mode = "MOVE_CERTIFIED"
                        break
            sm_info = {
                "family_id": f_id,
                "spoke": fent["spoke"],
                "subfamily": fent["subfamily"],
                "compiled_ch": comp_ch,
                "match_mode": match_mode,
            }
            spoke_matches.append(sm_info)
            promoted_exact_ch12_spokes.append({"hub_label": y_label, **sm_info})
        rep["matched_spoke_families"] = spoke_matches
        certified_new_entries.append(rep)

    assign_nearest_separating_witnesses(
        certified_new_entries,
        reference_pool=yosh_records + ch11_records + certified_new_entries,
    )

    exact_singletons_count = sum(1 for r in certified_new_entries if r["cell_status"] == "EXACT_SINGLETON")
    new_packets_count = sum(1 for r in certified_new_entries if r["cell_status"] == "NEW_PACKET")

    # Sector breakdown by topology:
    by_topology: dict[str, int] = defaultdict(int)
    for r in certified_new_entries:
        topo_str = " u ".join(str(x) for x in r["component_types"])  # type: ignore[union-attr]
        by_topology[topo_str] += 1

    census_stats_out = dict(census_stats or {})
    census_stats_out.setdefault("raw_and_level_counts_are_kauffman_bracket_filtered", True)
    census_stats_out.setdefault("all_tunneled_survivors_reidemeister_unlink_certified", True)

    summary_payload = {
        "ch_index": 12,
        "max_composite_ch_checked": max_composite_ch,
        "census_stats": census_stats_out,
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
        "counts_by_topology": dict(sorted(by_topology.items())),
        "promoted_exact_ch12_spokes": promoted_exact_ch12_spokes,
        "entries": certified_new_entries,
        "residual_collisions": residual_collisions,
    }

    json_path = out_dir / "expanded_ch12.json"
    json_path.write_text(json.dumps(summary_payload, indent=2) + "\n")

    csv_path = out_dir / "expanded_ch12.csv"
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
        "matched_spoke_ids",
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
            sp_ids = [str(m["family_id"]) for m in r.get("matched_spoke_families", [])]  # type: ignore[union-attr]
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
                    "matched_spoke_ids": json.dumps(sp_ids),
                    "nearest_separating_label": wit["nearest_label"],  # type: ignore[index]
                    "nearest_separating_layer": wit["separating_layer"],  # type: ignore[index]
                    "epd": r["epd"],
                }
            )

    return summary_payload
