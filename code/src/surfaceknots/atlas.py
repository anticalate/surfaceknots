"""Interactive KnotInfo-style HTML Web Atlas generator for the SurfaceKnots repository.

Unifies all canonical prime Hub classes (`ch <= 10` Yoshikawa + `ch = 11, 12` Expanded Censuses)
with all pairwise-distinct remaining Spoke family classes into a single self-contained explorer.
Includes an exact client-side JavaScript port of the IRLS `p`-Laplacian Tutte + block-cut pocket
layout + cubic Bezier hyperbolic splitting and LaTeX/TikZ engine, with optional deflate payload
compression (`compress=True`) and explicit epistemic rigor badges across Hub and Spoke tiers.
"""

from __future__ import annotations

import base64
import hashlib
import json
from pathlib import Path
import re
from typing import Any
import zlib

from surfaceknots.diagrams import MarkedDiagram


_EULER_TAG_RE = re.compile(r"e=(-?\d+)")


def _short_hash(obj: Any) -> str:
    raw = json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _o4_canonicalize_linking(items: Any) -> list[str] | None:
    """Canonicalize mod-2 double/triple linking rows up to global `O(4)` mirror `e -> -e`."""
    if items is None:
        return None
    if not items:
        return []

    def _flip_e(obj: Any) -> Any:
        if isinstance(obj, str):
            return _EULER_TAG_RE.sub(lambda m: f"e={-int(m.group(1))}", obj)
        if isinstance(obj, (list, tuple)):
            return [_flip_e(x) for x in obj]
        return obj

    pos = sorted(json.dumps(x, sort_keys=True, separators=(",", ":")) for x in items)
    neg = sorted(json.dumps(_flip_e(x), sort_keys=True, separators=(",", ":")) for x in items)
    return min(pos, neg)


def _fmt_torsion(tors: list[int]) -> str:
    return " \u2295 ".join(f"Z/{t}" for t in tors)


def _realization_rigor_badge(mode: str) -> str:
    if mode == "MOVE_CERTIFIED":
        return "PROVABLE (Move-Certified / Exact Singleton)"
    return "HEURISTIC (10-Layer Invariant Match, Unverified Move Path)"


def _normalize_entry(
    rec: dict[str, Any],
    *,
    tier_code: str,
    tier_label: str,
    hub_reals: list[dict[str, Any]],
) -> dict[str, Any]:
    ch_diag = int(rec["ch_index_diagram"])
    ch_lo = int(rec.get("ch_lower", ch_diag))
    ch_hi = int(rec.get("ch_upper", ch_diag))
    is_hub = tier_code in ("yoshikawa_ch10", "expanded_ch11", "expanded_ch12")
    exact_ch = bool(rec.get("exact_ch_certified", is_hub))

    if is_hub:
        cell_status = str(rec.get("cell_status", "EXACT_SINGLETON"))
        cell_reason = str(rec.get("cell_status_reason", "Canonical Hub prime census class."))
    else:
        cell_status = f"EXACT_CH{ch_lo}_SPOKE" if exact_ch else "SPOKE_BOUNDED"
        cell_reason = str(rec.get("complexity_certification_reason", ""))

    # Explicit epistemic rigor classification for ch(F) / Hub identification:
    if rec.get("hub_match_mode") == "INVARIANT_MATCH":
        rigor_code = "heuristic"
        rigor_badge = "HEURISTIC (10-Layer Invariant Match, Unverified Move Path)"
    elif is_hub:
        if cell_status == "NEW_PACKET":
            rigor_code = "packet"
            rigor_badge = "PACKET (Unmerged Diagrams, Same Invariants)"
        else:
            rigor_code = "provable"
            rigor_badge = "PROVABLE (Move-Certified / Exact Singleton)"
    elif exact_ch:
        rigor_code = "provable"
        rigor_badge = "PROVABLE (Move-Certified / Exact Singleton)"
    else:
        rigor_code = "bounded"
        rigor_badge = f"BOUNDED ([{ch_lo}, {ch_hi}])"

    has_heuristic_spoke_match = any(r.get("match_mode") != "MOVE_CERTIFIED" for r in hub_reals)
    spoke_match_badge = (
        "HEURISTIC (10-Layer Invariant Match, Unverified Move Path)"
        if has_heuristic_spoke_match
        else ("PROVABLE (Move-Certified / Exact Singleton)" if hub_reals else "")
    )

    sublinks = rec.get("sublinks") or {}
    sym = rec.get("symmetry") or {}

    # Layer 1: strict unoriented O(4) intrinsic topology & normal Euler multisets up to global e -> -e
    topo = MarkedDiagram.from_epd(rec["epd"]).surface_topology()
    c_pos = sorted([bool(c.orientable), int(c.genus_symbol), int(c.normal_euler)] for c in topo.components)
    c_neg = sorted([bool(c.orientable), int(c.genus_symbol), -int(c.normal_euler)] for c in topo.components)
    topo_o4 = min(c_pos, c_neg)

    c3_raw = rec["mochizuki_c3_mod3"]
    c3_mir = rec["mochizuki_c3_mirror_mod3"]
    c3_unord = sorted([c3_raw, c3_mir]) if (c3_raw is not None and c3_mir is not None) else None
    c5_raw = rec.get("mochizuki_c5_mod5")
    c5_mir = rec.get("mochizuki_c5_mirror_mod5")
    c5_unord = sorted([c5_raw, c5_mir]) if (c5_raw is not None and c5_mir is not None) else None
    c7_raw = rec.get("mochizuki_c7_mod7")
    c7_mir = rec.get("mochizuki_c7_mirror_mod7")
    c7_unord = sorted([c7_raw, c7_mir]) if (c7_raw is not None and c7_mir is not None) else None

    # Strict O(4) topological invariants ONLY across all 10 filtration layers:
    # NO presentation/orientation-dependent generators, greedy min_ch_found, or bounded-BFS booleans.
    raw_layers = [
        (rec["orientability"], int(rec["euler_characteristic"]), topo_o4),
        (int(rec["h1_free_rank"]), [int(t) for t in rec["h1_torsion"]]),
        (str(rec["e1_fingerprint_hash"]), [str(x) for x in (rec.get("twisted_s3_delta0") or [])]),
        (
            rec["hom_counts"],
            int(rec["sl2f3_hom_count"]),
            int(rec["sl2f3_non_meridional_peripheral_count"]),
        ),
        (
            rec["kei_invariant_factors"],
            rec["bikei_invariant_factors"],
            rec["alexander_biquandle_t"],
            rec["alexander_biquandle_s"],
        ),
        (
            c3_unord,
            c5_unord,
            c7_unord,
            int(rec["sym_quandle_q6_col"]),
            int(rec["sym_quandle_q6_trivial_weights"]),
            int(rec["sym_quandle_q6_nontrivial_weights"]),
        ),
        (
            _o4_canonicalize_linking(rec.get("double_linking_mod2")),
            _o4_canonicalize_linking(rec.get("triple_linking_mod2")),
            rec.get("sublink_invariant_profile", []),
        ),
        (
            int(rec["sigma2_euler"]),
            int(rec["sigma2_b1"]),
            [int(t) for t in rec["sigma2_h1_torsion"]],
            int(rec["sigma2_b2"]),
            rec["sigma2_abs_signature"],
        ),
        (rec.get("cover_betti_profile") or {},),
        (
            bool(sym.get("chiral_by_normal_euler")),
            bool(sym.get("oriented_chiral_by_cocycle")),
            bool(sym.get("non_invertible_by_cocycle")),
        ),
    ]
    layer_hashes = [_short_hash(x) for x in raw_layers]

    h1_str = f"Z^{rec['h1_free_rank']}" + (
        f" \u2295 {_fmt_torsion(rec['h1_torsion'])}" if rec["h1_torsion"] else ""
    )
    hom = rec["hom_counts"]
    hom_str = (
        f"S3:{hom['S3']} \u00b7 D8:{hom['D8']} \u00b7 A4:{hom['A4']} \u00b7 "
        f"S4:{hom['S4']} \u00b7 A5:{hom['A5']} \u00b7 SL(2,F3):{rec['sl2f3_hom_count']} "
        f"(periph:{rec['sl2f3_non_meridional_peripheral_count']})"
    )
    kei_r, kei_t = rec["kei_invariant_factors"]
    bik_r, bik_t = rec["bikei_invariant_factors"]
    kei_str = f"Kei=Z^{kei_r}" + (f"\u2295{kei_t}" if kei_t else "")
    bik_str = f"Bikei=Z^{bik_r}" + (f"\u2295{bik_t}" if bik_t else "")
    c3_str = f"[{', '.join(c3_raw)}]" if c3_raw is not None else "N/A (non-ori)"
    sig_val = rec["sigma2_abs_signature"]
    sig_str = str(sig_val) if sig_val is not None else "\u2014"
    sigma2_str = (
        f"\u03c7={rec['sigma2_euler']}, b1={rec['sigma2_b1']}, "
        f"tors={rec['sigma2_h1_torsion']}, b2={rec['sigma2_b2']}, |\u03c3|={sig_str}"
    )
    cov_keys = list((rec.get("cover_betti_profile") or {}).keys())
    cov_part = f"covers:[{','.join(cov_keys)}]" if cov_keys else "covers:none"
    pi2_rk = rec.get("certified_pi2_z_rank")
    pi2_reason = rec.get("pi2_certification_reason")
    pi2_str = f"rk={pi2_rk} ({pi2_reason}) | {cov_part}" if pi2_rk is not None else cov_part

    tw_s3 = rec.get("twisted_s3_delta0") or []
    tw_s3_part = f" | \u03940(S3)=[{', '.join(str(x) for x in tw_s3)}]" if tw_s3 else ""

    cup_sep = " \u222a "
    comp_joined = cup_sep.join(rec["component_types"])
    layer_displays = [
        f"{comp_joined} ({rec['orientability']}), \u03c7={rec['euler_characteristic']}, e={rec['normal_euler_numbers']}",
        h1_str,
        f"E1=({', '.join(rec['e1_generators'])}){tw_s3_part} | \u03941={rec['single_var_alex']}",
        hom_str,
        f"{kei_str} | {bik_str}",
        f"C3={c3_str} | Q6: col={rec['sym_quandle_q6_col']} (nontriv={rec['sym_quandle_q6_nontrivial_weights']})",
        f"dlk2={rec['double_linking_mod2']}, tlk2={rec['triple_linking_mod2']}, sublinks={len(rec.get('sublink_invariant_profile', []))}",
        sigma2_str,
        pi2_str,
        str(sym.get("summary", "")),
    ]

    nw = rec.get("nearest_separating_witness")
    nw_str = ""
    if isinstance(nw, dict):
        if nw.get("nearest_label"):
            nw_str = (
                f"Nearest: {nw['nearest_label']} "
                f"(shared {nw.get('shared_prefix_layers', 0)}/10 layers; "
                f"separates at L{nw.get('separating_layer_index', '?')} {nw.get('separating_layer', '')})"
            )
        elif nw.get("nearest_entry"):
            nw_str = (
                f"Nearest: {nw['nearest_entry']} (dist={nw.get('num_separating_layers', '?')}/10 via "
                f"{', '.join(nw.get('separating_layers', [])[:3])})"
            )

    return {
        "label": rec["label"],
        "tier_code": tier_code,
        "tier_label": tier_label,
        "spoke": rec.get("spoke", ""),
        "subfamily": rec.get("subfamily", ""),
        "desc": rec.get("description", ""),
        "ch_d": ch_diag,
        "ch_lo": ch_lo,
        "ch_hi": ch_hi,
        "ex_ch": exact_ch,
        "rigor_code": rigor_code,
        "rigor_badge": rigor_badge,
        "spoke_match_badge": spoke_match_badge,
        "c": int(rec["num_classical"]),
        "m": int(rec["num_marked"]),
        "cm": int(rec["c_minus"]),
        "cp": int(rec["c_plus"]),
        "nc": int(rec["num_components"]),
        "types": rec["component_types"],
        "ori": rec["orientability"],
        "chi": int(rec["euler_characteristic"]),
        "eul": rec["normal_euler_numbers"],
        "rib": rec.get("ribbon_status", "Unknown"),
        "cell": cell_status,
        "reason": cell_reason,
        "nw": nw_str,
        "err": bool(rec.get("has_yoshikawa_table_ii_erratum", False)),
        "err_note": rec.get("erratum_note", ""),
        "err_print": rec.get("yoshikawa_printed_e1", []),
        "pi1": rec["pi1_presentation"],
        "e1": rec["e1_generators"],
        "d1": rec["single_var_alex"],
        "brunn": bool(sublinks.get("is_brunnian", False)),
        "chiral": bool(sym.get("chiral_by_normal_euler") or sym.get("oriented_chiral_by_cocycle")),
        "q6_pos": int(rec["sym_quandle_q6_nontrivial_weights"]) > 0,
        "q_hs4": bool(rec["sigma2_is_q_homology_4sphere"]),
        "epd": rec["epd"],
        "reals": [
            {
                "id": r["family_id"],
                "mode": r["match_mode"],
                "ch": r["compiled_ch"],
                "rigor": _realization_rigor_badge(str(r["match_mode"])),
            }
            for r in hub_reals
        ],
        "lh": layer_hashes,
        "ld": layer_displays,
    }


def build_atlas_payload(data_dir: Path) -> dict[str, Any]:
    """Assemble the unified Hub + Spoke payload for the interactive Web Atlas."""
    from surfaceknots.yoshikawa_data import invariant_discrimination_key

    yosh_list = json.loads((data_dir / "yoshikawa_ch10.json").read_text())
    ch11_payload = json.loads((data_dir / "expanded_ch11.json").read_text())
    ch12_path = data_dir / "expanded_ch12.json"
    ch12_payload = json.loads(ch12_path.read_text()) if ch12_path.exists() else {}
    ch12_entries = ch12_payload.get("entries", [])
    fam_payload = json.loads((data_dir / "families" / "family_catalog.json").read_text())

    hub_reals_map: dict[str, list[dict[str, Any]]] = fam_payload.get("hub_realizations", {})
    promoted_spoke_ids: set[str] = {
        str(p["family_id"]) for p in ch12_payload.get("promoted_exact_ch12_spokes", []) if p.get("family_id")
    }
    entries: list[dict[str, Any]] = []
    seen_hub_iks: set[tuple] = set()

    for r in yosh_list:
        seen_hub_iks.add(invariant_discrimination_key(r))
        entries.append(
            _normalize_entry(
                r,
                tier_code="yoshikawa_ch10",
                tier_label="Yoshikawa (ch \u2264 10)",
                hub_reals=hub_reals_map.get(r["label"], []),
            )
        )

    for r in ch11_payload["entries"]:
        seen_hub_iks.add(invariant_discrimination_key(r))
        entries.append(
            _normalize_entry(
                r,
                tier_code="expanded_ch11",
                tier_label="Expanded Census (ch = 11)",
                hub_reals=hub_reals_map.get(r["label"], []),
            )
        )

    for r in ch12_entries:
        seen_hub_iks.add(invariant_discrimination_key(r))
        reals12 = hub_reals_map.get(r["label"], [])
        if not reals12 and r.get("matched_spoke_families"):
            reals12 = r["matched_spoke_families"]
        for sm in reals12:
            if sm.get("family_id"):
                promoted_spoke_ids.add(str(sm["family_id"]))
        entries.append(
            _normalize_entry(
                r,
                tier_code="expanded_ch12",
                tier_label="Expanded Census (ch = 12)",
                hub_reals=reals12,
            )
        )

    seen_spoke_iks: set[tuple] = set()
    for r in fam_payload["entries"]:
        fam_id = str(r.get("family_id") or r.get("label") or "")
        if r.get("hub_match") is not None or r.get("reference_collisions") or fam_id in promoted_spoke_ids:
            continue
        ik = invariant_discrimination_key(r)
        if ik in seen_hub_iks or ik in seen_spoke_iks:
            continue
        seen_spoke_iks.add(ik)
        ch_lo = int(r["ch_lower"])
        is_ex = bool(r.get("exact_ch_certified")) and ch_lo >= 12
        entries.append(
            _normalize_entry(
                r,
                tier_code=f"exact_ch{ch_lo}" if is_ex else "spoke_bounded",
                tier_label=f"Exact ch = {ch_lo} (Spoke)" if is_ex else f"Spoke Catalog (ch \u2265 {ch_lo})",
                hub_reals=[],
            )
        )

    # Convert per-layer hashes into exact collision-free integer equivalence IDs 0..k
    for layer_i in range(10):
        seen_ids: dict[str, int] = {}
        for ent in entries:
            h_val = ent["lh"][layer_i]
            if h_val not in seen_ids:
                seen_ids[h_val] = len(seen_ids)
            ent["lh"][layer_i] = seen_ids[h_val]

    hub_total = len(yosh_list) + len(ch11_payload["entries"]) + len(ch12_entries)
    spoke_total = len(entries) - hub_total
    exact_total = sum(1 for e in entries if e["ex_ch"])
    hub_realized = sum(
        1
        for e in entries
        if e["tier_code"] in ("yoshikawa_ch10", "expanded_ch11", "expanded_ch12") and len(e["reals"]) > 0
    )

    return {
        "version": "0.6.0",
        "counts": {
            "total_distinct_classes": len(entries),
            "hub_total": hub_total,
            "yoshikawa_ch10": len(yosh_list),
            "expanded_ch11": len(ch11_payload["entries"]),
            "expanded_ch12": len(ch12_entries),
            "exact_ch12": sum(1 for e in entries if e["ex_ch"] and e["ch_lo"] == 12),
            "exact_ch_total": exact_total,
            "spoke_total": spoke_total,
            "hub_realized_count": hub_realized,
            "compiled_family_instances": int(fam_payload["total_family_instances"]),
        },
        "entries": entries,
    }


_ATLAS_HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>SurfaceKnots Atlas — Knotted Surfaces &amp; Surface-Links in R^4</title>
  <style>
    :root {
      --background: #f8fafc;
      --foreground: #0f172a;
      --card: #ffffff;
      --card-border: #e2e8f0;
      --border: #cbd5e1;
      --secondary: #e2e8f0;
      --secondary-foreground: #0f172a;
      --muted-foreground: #475569;
      --primary: #1e3a8a;
      --primary-foreground: #ffffff;
      --link: #1d4ed8;
      --warning: #b45309;
      --error: #dc2626;
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    }
    @media (prefers-color-scheme: dark) {
      :root:not(.light) {
        --background: #0b1120;
        --foreground: #f1f5f9;
        --card: #111827;
        --card-border: #1e293b;
        --border: #334155;
        --secondary: #1e293b;
        --secondary-foreground: #f1f5f9;
        --muted-foreground: #94a3b8;
        --primary: #3b82f6;
        --primary-foreground: #ffffff;
        --link: #60a5fa;
        --warning: #f59e0b;
        --error: #f87171;
      }
    }
    *, *::before, *::after { box-sizing: border-box; }
    body { background: var(--background); color: var(--foreground); margin: 0; padding: 1rem; line-height: 1.45; -webkit-font-smoothing: antialiased; }
    button, input, select, textarea { font: inherit; color: inherit; }
    button { cursor: pointer; background: transparent; }
    code, pre, .font-mono { font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace; }
    .hidden { display: none !important; }
    .block { display: block; }
    .inline-block { display: inline-block; }
    .inline-flex { display: inline-flex; }
    .flex { display: flex; }
    .flex-col { flex-direction: column; }
    .flex-wrap { flex-wrap: wrap; }
    .items-center { align-items: center; }
    .items-start { align-items: flex-start; }
    .justify-between { justify-content: space-between; }
    .justify-end { justify-content: flex-end; }
    .shrink-0 { flex-shrink: 0; }
    .min-w-0 { min-width: 0; }
    .w-full { width: 100%; }
    .w-2\.5 { width: 0.625rem; }
    .h-2\.5 { height: 0.625rem; }
    .max-w-7xl { max-width: 82rem; }
    .max-w-\[225px\] { max-width: 225px; }
    .max-h-\[650px\] { max-height: 650px; }
    .max-h-72 { max-height: 18rem; }
    .max-h-56 { max-height: 14rem; }
    .max-h-52 { max-height: 13rem; }
    .mx-auto { margin-left: auto; margin-right: auto; }
    .m-0 { margin: 0; }
    .mt-0\.5 { margin-top: 0.125rem; }
    .mt-1 { margin-top: 0.25rem; }
    .mt-2 { margin-top: 0.5rem; }
    .mb-1 { margin-bottom: 0.25rem; }
    .p-2 { padding: 0.5rem; }
    .p-2\.5 { padding: 0.625rem; }
    .p-3 { padding: 0.75rem; }
    .p-3\.5 { padding: 0.875rem; }
    .p-4 { padding: 1rem; }
    .p-6 { padding: 1.5rem; }
    .px-1 { padding-left: 0.25rem; padding-right: 0.25rem; }
    .px-1\.5 { padding-left: 0.375rem; padding-right: 0.375rem; }
    .px-2 { padding-left: 0.5rem; padding-right: 0.5rem; }
    .px-2\.5 { padding-left: 0.625rem; padding-right: 0.625rem; }
    .px-3 { padding-left: 0.75rem; padding-right: 0.75rem; }
    .px-4 { padding-left: 1rem; padding-right: 1rem; }
    .py-0\.5 { padding-top: 0.125rem; padding-bottom: 0.125rem; }
    .py-1 { padding-top: 0.25rem; padding-bottom: 0.25rem; }
    .py-1\.5 { padding-top: 0.375rem; padding-bottom: 0.375rem; }
    .py-2 { padding-top: 0.5rem; padding-bottom: 0.5rem; }
    .py-2\.5 { padding-top: 0.625rem; padding-bottom: 0.625rem; }
    .pt-1 { padding-top: 0.25rem; }
    .pb-1 { padding-bottom: 0.25rem; }
    .pb-2 { padding-bottom: 0.5rem; }
    .pb-3 { padding-bottom: 0.75rem; }
    .gap-1\.5 { gap: 0.375rem; }
    .gap-2 { gap: 0.5rem; }
    .gap-2\.5 { gap: 0.625rem; }
    .gap-3 { gap: 0.75rem; }
    .gap-4 { gap: 1rem; }
    .space-y-0\.5 > * + * { margin-top: 0.125rem; }
    .space-y-1 > * + * { margin-top: 0.25rem; }
    .space-y-1\.5 > * + * { margin-top: 0.375rem; }
    .space-y-2 > * + * { margin-top: 0.5rem; }
    .space-y-2\.5 > * + * { margin-top: 0.625rem; }
    .space-y-3 > * + * { margin-top: 0.75rem; }
    .space-y-4 > * + * { margin-top: 1rem; }
    .divide-y > * + * { border-top: 1px solid var(--card-border); }
    .grid { display: grid; }
    .grid-cols-1 { grid-template-columns: repeat(1, minmax(0, 1fr)); }
    .grid-cols-2 { grid-template-columns: repeat(2, minmax(0, 1fr)); }
    @media (min-width: 640px) {
      .sm\:grid-cols-2 { grid-template-columns: repeat(2, minmax(0, 1fr)); }
      .sm\:grid-cols-3 { grid-template-columns: repeat(3, minmax(0, 1fr)); }
    }
    @media (min-width: 768px) {
      body { padding: 1.25rem; }
      .md\:grid-cols-2 { grid-template-columns: repeat(2, minmax(0, 1fr)); }
      .md\:text-lg { font-size: 1.125rem; }
    }
    @media (min-width: 1024px) {
      .lg\:flex-row { flex-direction: row; }
      .lg\:items-center { align-items: center; }
      .lg\:justify-between { justify-content: space-between; }
      .lg\:grid-cols-2 { grid-template-columns: repeat(2, minmax(0, 1fr)); }
      .lg\:grid-cols-4 { grid-template-columns: repeat(4, minmax(0, 1fr)); }
      .lg\:grid-cols-5 { grid-template-columns: repeat(5, minmax(0, 1fr)); }
      .lg\:grid-cols-12 { grid-template-columns: repeat(12, minmax(0, 1fr)); }
      .lg\:col-span-5 { grid-column: span 5 / span 5; }
      .lg\:col-span-7 { grid-column: span 7 / span 7; }
    }
    .rounded { border-radius: 0.25rem; }
    .rounded-md { border-radius: 0.375rem; }
    .rounded-lg { border-radius: 0.5rem; }
    .rounded-xl { border-radius: 0.75rem; }
    .rounded-full { border-radius: 9999px; }
    .border { border: 1px solid var(--border); }
    .border-t { border-top: 1px solid var(--card-border); }
    .border-b { border-bottom: 1px solid var(--card-border); }
    .border-transparent { border-color: transparent; }
    .border-\[var\(--border\)\] { border-color: var(--border); }
    .border-\[var\(--card-border\)\] { border-color: var(--card-border); }
    .border-\[var\(--warning\)\] { border-color: var(--warning); }
    .bg-\[var\(--background\)\] { background-color: var(--background); }
    .bg-\[var\(--card\)\] { background-color: var(--card); }
    .bg-\[var\(--secondary\)\] { background-color: var(--secondary); }
    .bg-\[var\(--primary\)\] { background-color: var(--primary); }
    .text-\[var\(--foreground\)\] { color: var(--foreground); }
    .text-\[var\(--muted-foreground\)\] { color: var(--muted-foreground); }
    .text-\[var\(--secondary-foreground\)\] { color: var(--secondary-foreground); }
    .text-\[var\(--primary-foreground\)\] { color: var(--primary-foreground); }
    .text-\[var\(--link\)\] { color: var(--link); }
    .text-\[var\(--warning\)\] { color: var(--warning); }
    .text-\[var\(--error\)\] { color: var(--error); }
    .hover\:bg-\[var\(--background\)\]:hover { background-color: var(--background); }
    .hover\:bg-\[var\(--secondary\)\]:hover { background-color: var(--secondary); }
    .hover\:underline:hover { text-decoration: underline; }
    .focus\:outline-none:focus { outline: none; }
    .text-xs { font-size: 0.75rem; }
    .text-sm { font-size: 0.875rem; }
    .text-base { font-size: 1rem; }
    .text-xl { font-size: 1.25rem; }
    .font-medium { font-weight: 500; }
    .font-semibold { font-weight: 600; }
    .font-bold { font-weight: 700; }
    .text-left { text-align: left; }
    .text-center { text-align: center; }
    .text-right { text-align: right; }
    .truncate { overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
    .break-words { overflow-wrap: break-word; }
    .break-all { word-break: break-all; }
    .whitespace-pre-wrap { white-space: pre-wrap; }
    .overflow-hidden { overflow: hidden; }
    .overflow-y-auto { overflow-y: auto; }
    .overflow-x-auto { overflow-x: auto; }
    .border-collapse { border-collapse: collapse; }
    .transition { transition: background-color 0.15s ease, border-color 0.15s ease, color 0.15s ease; }
    .cursor-pointer { cursor: pointer; }
  </style>
</head>
<body class="bg-[var(--background)] text-[var(--foreground)] antialiased m-0 p-3 md:p-5">
  <div class="max-w-7xl mx-auto space-y-4">
    <!-- Header Banner -->
    <header class="bg-[var(--card)] border border-[var(--border)] rounded-xl p-4 flex flex-col lg:flex-row lg:items-center lg:justify-between gap-3">
      <div class="space-y-1">
        <div class="flex flex-wrap items-center gap-2">
          <h1 class="text-base md:text-lg font-bold tracking-tight text-[var(--foreground)]">
            SurfaceKnots Atlas
          </h1>
          <span id="hdr-badge" class="px-2 py-0.5 rounded-md text-xs font-mono bg-[var(--secondary)] text-[var(--secondary-foreground)] border border-[var(--border)]">
            v0.6.0 &middot; Certified Classes in R&#8308;
          </span>
        </div>
        <p class="text-xs text-[var(--muted-foreground)]">
          Interactive Explorer: Yoshikawa (ch &le; 10) &middot; Prime Censuses (ch = 11, 12) &middot; Parametric Families
        </p>
      </div>
      <!-- Primary Mode Tabs -->
      <nav class="flex flex-wrap items-center gap-1.5" role="tablist" aria-label="Atlas views">
        <button id="tab-explorer" type="button" data-main-tab="explorer"
          class="px-3 py-1.5 rounded-lg text-xs font-medium transition border">
          Atlas Explorer
        </button>
        <button id="tab-compare" type="button" data-main-tab="compare"
          class="px-3 py-1.5 rounded-lg text-xs font-medium transition border">
          Pairwise Separator
        </button>
        <button id="tab-sandbox" type="button" data-main-tab="sandbox"
          class="px-3 py-1.5 rounded-lg text-xs font-medium transition border">
          Live EPD Sandbox
        </button>
        <button id="tab-overview" type="button" data-main-tab="overview"
          class="px-3 py-1.5 rounded-lg text-xs font-medium transition border">
          Census Summary
        </button>
      </nav>
    </header>

    <!-- VIEW 1: ATLAS EXPLORER -->
    <section id="view-explorer" class="space-y-3">
      <div class="bg-[var(--card)] border border-[var(--border)] rounded-xl p-3 grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-5 gap-2.5">
        <div>
          <label for="f-search" class="block text-xs font-medium text-[var(--muted-foreground)] mb-1">Search label / invariant</label>
          <input id="f-search" type="text" placeholder="e.g. 11_1, Klein, Borromean, 2t-1..."
            class="w-full px-2.5 py-1.5 rounded-lg text-xs bg-[var(--background)] text-[var(--foreground)] border border-[var(--border)] focus:outline-none" />
        </div>
        <div>
          <label for="f-tier" class="block text-xs font-medium text-[var(--muted-foreground)] mb-1">Catalog Tier</label>
          <select id="f-tier" class="w-full px-2.5 py-1.5 rounded-lg text-xs bg-[var(--background)] text-[var(--foreground)] border border-[var(--border)]">
            <option value="all">All Tiers</option>
            <option value="hub">Canonical Hub (ch &le; 12)</option>
            <option value="yoshikawa_ch10">Yoshikawa Table (ch &le; 10)</option>
            <option value="expanded_ch11">Expanded Census (ch = 11)</option>
            <option value="expanded_ch12">Expanded Census (ch = 12)</option>
            <option value="spoke">Parametric Spoke Catalog</option>
          </select>
        </div>
        <div>
          <label for="f-topo" class="block text-xs font-medium text-[var(--muted-foreground)] mb-1">Surface Topology</label>
          <select id="f-topo" class="w-full px-2.5 py-1.5 rounded-lg text-xs bg-[var(--background)] text-[var(--foreground)] border border-[var(--border)]">
            <option value="all">All Topologies</option>
            <option value="s2_knot">2-Knots (S&#178;)</option>
            <option value="t2_knot">Torus Knots (T&#178;)</option>
            <option value="nonori_knot">Non-Orientable Knots (RP&#178;, K&#178;)</option>
            <option value="link_2">2-Component Surface-Links</option>
            <option value="link_3">3-Component Surface-Links</option>
            <option value="brunnian">Brunnian Surface-Links</option>
          </select>
        </div>
        <div>
          <label for="f-ori" class="block text-xs font-medium text-[var(--muted-foreground)] mb-1">Orientability</label>
          <select id="f-ori" class="w-full px-2.5 py-1.5 rounded-lg text-xs bg-[var(--background)] text-[var(--foreground)] border border-[var(--border)]">
            <option value="all">All (Ori / Non / Mixed)</option>
            <option value="Ori">Orientable (Ori)</option>
            <option value="Non">Non-Orientable (Non)</option>
            <option value="Mixed">Mixed Orientability</option>
          </select>
        </div>
        <div>
          <label for="f-prop" class="block text-xs font-medium text-[var(--muted-foreground)] mb-1">Certification / Property</label>
          <select id="f-prop" class="w-full px-2.5 py-1.5 rounded-lg text-xs bg-[var(--background)] text-[var(--foreground)] border border-[var(--border)]">
            <option value="all">All Entries</option>
            <option value="exact_ch">Exact ch(F) Proven</option>
            <option value="rigor_provable">Rigor: PROVABLE (Move / Singleton)</option>
            <option value="rigor_packet">Rigor: PACKET (Unmerged Same-Invariant)</option>
            <option value="rigor_heuristic">Rigor: HEURISTIC (10-Layer Match)</option>
            <option value="rigor_bounded">Rigor: BOUNDED ([ch_min, ch_max])</option>
            <option value="hub_realized">Hub Entries w/ Spoke Model</option>
            <option value="chiral">Certified Chiral</option>
            <option value="q6_triple">Positive Triple Points (Q&#8326;)</option>
            <option value="erratum">Yoshikawa Table II Erratum</option>
          </select>
        </div>
      </div>

      <div class="grid grid-cols-1 lg:grid-cols-12 gap-4 items-start">
        <div class="lg:col-span-5 bg-[var(--card)] border border-[var(--border)] rounded-xl overflow-hidden">
          <div class="px-3 py-2 border-b border-[var(--card-border)] flex items-center justify-between text-xs">
            <span id="filter-count" class="font-medium text-[var(--foreground)]">Entries</span>
            <button id="btn-reset-filters" type="button" class="text-xs text-[var(--link)] hover:underline">
              Reset filters
            </button>
          </div>
          <div class="max-h-[650px] overflow-y-auto divide-y divide-[var(--card-border)]" id="catalog-list"></div>
        </div>
        <div class="lg:col-span-7 bg-[var(--card)] border border-[var(--border)] rounded-xl p-4 space-y-4" id="detail-panel"></div>
      </div>
    </section>

    <!-- VIEW 2: PAIRWISE INVARIANT SEPARATOR -->
    <section id="view-compare" class="hidden space-y-4">
      <div class="bg-[var(--card)] border border-[var(--border)] rounded-xl p-4 space-y-3">
        <div class="flex flex-wrap items-center justify-between gap-2">
          <div>
            <h2 class="text-sm font-bold text-[var(--foreground)]">10-Layer Pairwise Invariant Separator</h2>
            <p class="text-xs text-[var(--muted-foreground)]">
              Compare any two surface-links side-by-side across all 10 layers of the exact topological invariant filtration.
            </p>
          </div>
          <div class="flex flex-wrap gap-1.5" id="compare-presets"></div>
        </div>
        <div class="grid grid-cols-1 md:grid-cols-2 gap-3 pt-1">
          <div>
            <label for="cmp-a" class="block text-xs font-medium text-[var(--muted-foreground)] mb-1">Surface A</label>
            <select id="cmp-a" class="w-full px-2.5 py-1.5 rounded-lg text-xs font-mono bg-[var(--background)] text-[var(--foreground)] border border-[var(--border)]"></select>
          </div>
          <div>
            <label for="cmp-b" class="block text-xs font-medium text-[var(--muted-foreground)] mb-1">Surface B</label>
            <select id="cmp-b" class="w-full px-2.5 py-1.5 rounded-lg text-xs font-mono bg-[var(--background)] text-[var(--foreground)] border border-[var(--border)]"></select>
          </div>
        </div>
      </div>
      <div id="compare-output" class="space-y-4"></div>
    </section>

    <!-- VIEW 3: LIVE EPD SANDBOX -->
    <section id="view-sandbox" class="hidden space-y-4">
      <div class="bg-[var(--card)] border border-[var(--border)] rounded-xl p-4 space-y-3">
        <div class="flex flex-wrap items-center justify-between gap-2">
          <div>
            <h2 class="text-sm font-bold text-[var(--foreground)]">Live Planar EPD &amp; Splitting Movie Renderer</h2>
            <p class="text-xs text-[var(--muted-foreground)]">
              Paste any planar marked graph EPD (<code class="font-mono">X[a,b,c,d], Y[...], Z[...]</code>) to compute its crossing-free IRLS p-Laplacian Tutte layout, hyperbolic splitting movie, and LaTeX/TikZ source in real time.
            </p>
          </div>
          <button id="btn-render-epd" type="button"
            class="px-3 py-1.5 rounded-lg text-xs font-medium bg-[var(--primary)] text-[var(--primary-foreground)]">
            Render EPD
          </button>
        </div>
        <textarea id="sandbox-epd" rows="3"
          class="w-full p-2.5 rounded-lg text-xs font-mono bg-[var(--background)] text-[var(--foreground)] border border-[var(--border)]">Y[1,6,2,5], Z[2,6,3,7], X[3,12,4,11], X[4,9,5,10], X[1,8,12,9], X[7,11,8,10]</textarea>
      </div>
      <div id="sandbox-output" class="bg-[var(--card)] border border-[var(--border)] rounded-xl p-4"></div>
    </section>

    <!-- VIEW 4: HUB & SPOKE ARCHITECTURE SUMMARY -->
    <section id="view-overview" class="hidden space-y-4"></section>
  </div>

  __ATLAS_SCRIPT_TAG__
  <script>
  (async function () {
    const dataEl = document.getElementById('atlas-data');
    let DATA;
    if (dataEl.getAttribute('data-encoding') === 'deflate-b64') {
      const bin = atob(dataEl.textContent.trim());
      const bytes = new Uint8Array(bin.length);
      for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
      const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('deflate'));
      DATA = JSON.parse(await new Response(stream).text());
    } else {
      DATA = JSON.parse(dataEl.textContent);
    }
    const ENTRIES = DATA.entries;
    const BY_LABEL = Object.fromEntries(ENTRIES.map((e) => [e.label, e]));

    const PALETTE_LIGHT = ['#1E3A8A', '#C2410C', '#0F766E', '#7E22CE', '#B45309', '#BE123C'];
    const PALETTE_DARK  = ['#60A5FA', '#FB923C', '#2DD4BF', '#C084FC', '#FBBF24', '#FB7185'];
    const TIKZ_RGB = [[30,58,138],[194,65,12],[15,118,110],[126,34,206],[180,83,9],[190,18,60]];

    const layoutCache = new Map();

    /* ---------- Exact IRLS p-Laplacian Tutte & Splitting Engine ---------- */
    function norm2(vx, vy) {
      const d = Math.hypot(vx, vy);
      return d < 1e-12 ? [1.0, 0.0] : [vx / d, vy / d];
    }

    function solveLinearSystem(aMat, rhsX, rhsY) {
      const n = aMat.length;
      if (n === 0) return [[], []];
      const aug = aMat.map((row, i) => [...row, rhsX[i], rhsY[i]]);
      for (let col = 0; col < n; col++) {
        let piv = col;
        let best = Math.abs(aug[col][col]);
        for (let r = col + 1; r < n; r++) {
          const v = Math.abs(aug[r][col]);
          if (v > best) { best = v; piv = r; }
        }
        if (best < 1e-12) throw new Error('Non-planar or singular Tutte block');
        if (piv !== col) { const tmp = aug[col]; aug[col] = aug[piv]; aug[piv] = tmp; }
        const inv = 1.0 / aug[col][col];
        const rc = aug[col];
        for (let j = col; j < n + 2; j++) rc[j] *= inv;
        for (let r = 0; r < n; r++) {
          if (r === col) continue;
          const f = aug[r][col];
          if (Math.abs(f) > 1e-15) {
            const rr = aug[r];
            for (let j = col; j < n + 2; j++) rr[j] -= f * rc[j];
          }
        }
      }
      return [aug.map((r) => r[n]), aug.map((r) => r[n + 1])];
    }

    function parseEpd(epdStr) {
      const s = (epdStr || '').trim();
      if (!s || s === 'O[1]' || s === '0_1') {
        return { n: 0, opp: [], vtype: [], freeLoops: 1 };
      }
      let freeLoops = 0;
      const oMatches = s.match(/O\s*\[\s*1\s*\]/g);
      if (oMatches) freeLoops = oMatches.length;
      const re = /([XYZ][01]?)\s*\[\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\]/g;
      const items = [];
      let m;
      while ((m = re.exec(s)) !== null) {
        items.push([m[1], Number(m[2]), Number(m[3]), Number(m[4]), Number(m[5])]);
      }
      if (items.length === 0 && freeLoops > 0) {
        return { n: 0, opp: [], vtype: [], freeLoops };
      }
      if (items.length === 0) throw new Error('Could not parse any X/Y/Z[a,b,c,d] quadruples from EPD.');
      const n = items.length;
      const vtype = [];
      const edgeOcc = new Map();
      for (let v = 0; v < n; v++) {
        const [k, a, b, c, d] = items[v];
        vtype.push(k === 'X' ? 'X0' : k);
        [a, b, c, d].forEach((lbl, slot) => {
          if (!edgeOcc.has(lbl)) edgeOcc.set(lbl, []);
          edgeOcc.get(lbl).push(4 * v + slot);
        });
      }
      const opp = new Array(4 * n).fill(0);
      for (const [lbl, darts] of edgeOcc.entries()) {
        if (darts.length !== 2) throw new Error(`Edge label ${lbl} appears ${darts.length} times (expected 2).`);
        opp[darts[0]] = darts[1];
        opp[darts[1]] = darts[0];
      }
      return { n, opp, vtype, freeLoops };
    }

    function dartSurfaceComponents(diag) {
      const { n, opp, vtype, freeLoops } = diag;
      const parent = Array.from({ length: 4 * n }, (_, i) => i);
      const find = (x) => { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
      const union = (a, b) => { const ra = find(a), rb = find(b); if (ra !== rb) parent[rb] = ra; };
      for (let d = 0; d < 4 * n; d++) union(d, opp[d]);
      for (let v = 0; v < n; v++) {
        if (vtype[v] === 'X0' || vtype[v] === 'X1') {
          union(4 * v, 4 * v + 2); union(4 * v + 1, 4 * v + 3);
        } else {
          union(4 * v, 4 * v + 1); union(4 * v + 1, 4 * v + 2); union(4 * v + 2, 4 * v + 3);
        }
      }
      const rootMap = new Map();
      const comp = new Array(4 * n).fill(0);
      for (let d = 0; d < 4 * n; d++) {
        const r = find(d);
        if (!rootMap.has(r)) rootMap.set(r, rootMap.size);
        comp[d] = rootMap.get(r);
      }
      return { dartComp: comp, graphComps: rootMap.size, totalComps: rootMap.size + freeLoops };
    }

    function countResolutionCircles(diag, sign) {
      const { n, opp, vtype, freeLoops } = diag;
      const vis = new Array(4 * n).fill(false);
      let circles = freeLoops;
      for (let start = 0; start < 4 * n; start++) {
        if (vis[start]) continue;
        circles++;
        let cur = start;
        while (!vis[cur]) {
          vis[cur] = true;
          const v = cur >> 2, s = cur & 3, vt = vtype[v];
          let partner;
          if (vt === 'X0' || vt === 'X1') {
            partner = (v << 2) | (s ^ 2);
          } else {
            const p03 = (vt === 'Y' && sign === 1) || (vt === 'Z' && sign === -1);
            const map = p03 ? [3, 2, 1, 0] : [1, 0, 3, 2];
            partner = (v << 2) | map[s];
          }
          vis[partner] = true;
          cur = opp[partner];
        }
      }
      return circles;
    }

    const TWO_PI = 2.0 * Math.PI;
    function posMod2Pi(x) {
      const r = x % TWO_PI;
      return r < 0 ? r + TWO_PI : r;
    }

    function verifyCcwCyclicOrder(px, py, qs) {
      const angs = qs.map(([qx, qy]) => Math.atan2(qy - py, qx - px));
      let total = 0.0;
      for (let s = 0; s < 4; s++) {
        const diff = posMod2Pi(angs[(s + 1) & 3] - angs[s]);
        if (diff < 1e-5 || diff > TWO_PI - 1e-5) return false;
        total += diff;
      }
      return Math.abs(total - TWO_PI) < 1e-3;
    }

    function findSelfLoop(m, localOpp) {
      if (m < 2) return null;
      for (let v = 0; v < m; v++) {
        for (let s0 = 0; s0 < 4; s0++) {
          if (localOpp[4 * v + s0] === 4 * v + ((s0 + 1) & 3)) return [v, s0];
        }
      }
      return null;
    }

    function findCutSplit(m, localOpp) {
      if (m < 3) return null;
      for (let v = 0; v < m; v++) {
        const seen = new Set([v]);
        const comps = [];
        for (let st = 0; st < m; st++) {
          if (seen.has(st)) continue;
          const comp = [];
          const stack = [st];
          seen.add(st);
          while (stack.length) {
            const u = stack.pop();
            comp.push(u);
            for (let s = 0; s < 4; s++) {
              const w = localOpp[4 * u + s] >> 2;
              if (!seen.has(w)) { seen.add(w); stack.push(w); }
            }
          }
          comp.sort((a, b) => a - b);
          comps.push(comp);
        }
        if (comps.length <= 1) continue;
        const vOfSlot = [0, 1, 2, 3].map((s) => localOpp[4 * v + s] >> 2);
        const compIdOf = new Map();
        comps.forEach((cverts, cid) => cverts.forEach((u) => compIdOf.set(u, cid)));
        for (let s0 = 0; s0 < 4; s0++) {
          const s1 = (s0 + 1) & 3, s2 = (s0 + 2) & 3, s3 = (s0 + 3) & 3;
          const ca = new Set([compIdOf.get(vOfSlot[s0]), compIdOf.get(vOfSlot[s1])]);
          const cb = new Set([compIdOf.get(vOfSlot[s2]), compIdOf.get(vOfSlot[s3])]);
          let disj = true;
          for (const x of ca) if (cb.has(x)) disj = false;
          if (ca.size && cb.size && disj) {
            const vertsA = [], vertsB = [];
            comps.forEach((cv, cid) => { if (ca.has(cid)) vertsA.push(...cv); else vertsB.push(...cv); });
            vertsA.sort((a, b) => a - b);
            vertsB.sort((a, b) => a - b);
            return vertsA.length >= vertsB.length ? [v, s0, vertsA, vertsB] : [v, s2, vertsB, vertsA];
          }
        }
      }
      return null;
    }

    function layout2ConnectedLoopFreeBlock(m, localOpp, localVtype, forcedOuterDart = null) {
      if (m === 1) {
        if (localOpp[0] === 2) throw new Error('Non-planar 1-vertex crossing loop.');
        const sub = [[0, 0]];
        for (let s = 0; s < 4; s++) sub.push([0.65 * Math.cos(0.5 * Math.PI * s), 0.65 * Math.sin(0.5 * Math.PI * s)]);
        return sub;
      }
      const visited = new Array(4 * m).fill(false);
      const faces = [];
      for (let st = 0; st < 4 * m; st++) {
        if (visited[st]) continue;
        const cyc = [];
        let cur = st;
        while (!visited[cur]) {
          visited[cur] = true;
          cyc.push(cur);
          const od = localOpp[cur];
          cur = (od & ~3) | (((od & 3) + 3) & 3);
        }
        faces.push(cyc);
      }
      if (faces.length !== m + 2) {
        throw new Error(`Non-planar shadow: V=${m}, E=${2 * m}, F=${faces.length} (expected ${m + 2}).`);
      }
      const faceWalks = faces.map((cyc) => {
        const w = [];
        for (const hd of cyc) w.push(hd >> 2, m + hd, m + localOpp[hd]);
        return w;
      });
      let outerIdx = 0;
      let bestKey = [-1, -1, -1, -1, 1];
      for (let fi = 0; fi < faces.length; fi++) {
        const wf = faceWalks[fi];
        const hasForced = (forcedOuterDart !== null && faces[fi].includes(forcedOuterDart)) ? 1 : 0;
        const simple = new Set(wf).size === wf.length ? 1 : 0;
        const mkCnt = faces[fi].reduce((acc, hd) => acc + ((localVtype[hd >> 2] === 'Y' || localVtype[hd >> 2] === 'Z') ? 1 : 0), 0);
        const key = [hasForced, simple, faces[fi].length, mkCnt, -fi];
        for (let k = 0; k < 5; k++) {
          if (key[k] !== bestKey[k]) { if (key[k] > bestKey[k]) { bestKey = key; outerIdx = fi; } break; }
        }
      }
      const wOut = [];
      const seenOut = new Set();
      for (const nd of faceWalks[outerIdx]) { if (!seenOut.has(nd)) { seenOut.add(nd); wOut.push(nd); } }
      const interior = [];
      for (let fi = 0; fi < faces.length; fi++) if (fi !== outerIdx) interior.push(fi);
      const totalAug = 5 * m + interior.length;

      const fixedPos = new Map();
      const L = wOut.length;
      for (let j = 0; j < L; j++) {
        const th = Math.PI - (2.0 * Math.PI * j) / L;
        fixedPos.set(wOut[j], [Math.cos(th), Math.sin(th)]);
      }

      const augEdges = [];
      for (let hd = 0; hd < 4 * m; hd++) {
        augEdges.push([hd >> 2, m + hd, 1.35, true]);
        const od = localOpp[hd];
        if (hd < od) augEdges.push([m + hd, m + od, 1.0, true]);
      }
      interior.forEach((fi, idxIn) => {
        const cf = 5 * m + idxIn;
        const deg = faces[fi].length;
        const uniq = [...new Set(faceWalks[fi])];
        for (const nd of uniq) {
          const wr = deg <= 2 ? (nd < m ? 0.85 : 0.65) : (nd < m ? 0.42 : 0.28);
          augEdges.push([cf, nd, wr, false]);
        }
      });

      const freeNodes = [];
      const freeIdx = new Map();
      for (let nd = 0; nd < totalAug; nd++) {
        if (!fixedPos.has(nd)) { freeIdx.set(nd, freeNodes.length); freeNodes.push(nd); }
      }
      const nFree = freeNodes.length;
      const curX = new Array(totalAug).fill(0);
      const curY = new Array(totalAug).fill(0);
      for (const [nd, [fx, fy]] of fixedPos.entries()) { curX[nd] = fx; curY[nd] = fy; }

      for (let it = 0; it < 10; it++) {
        const mat = Array.from({ length: nFree }, () => new Array(nFree).fill(0));
        const rx = new Array(nFree).fill(0);
        const ry = new Array(nFree).fill(0);
        for (const [u, v, bw, isPrim] of augEdges) {
          let w = bw;
          if (it > 0) {
            const dist = Math.hypot(curX[u] - curX[v], curY[u] - curY[v]);
            w = bw * Math.pow(dist + 0.035, isPrim ? 1.35 : 1.05);
          }
          const iu = freeIdx.get(u);
          const iv = freeIdx.get(v);
          if (iu !== undefined) {
            mat[iu][iu] += w;
            if (iv !== undefined) mat[iu][iv] -= w;
            else { rx[iu] += w * curX[v]; ry[iu] += w * curY[v]; }
          }
          if (iv !== undefined) {
            mat[iv][iv] += w;
            if (iu !== undefined) mat[iv][iu] -= w;
            else { rx[iv] += w * curX[u]; ry[iv] += w * curY[u]; }
          }
        }
        const [sx, sy] = solveLinearSystem(mat, rx, ry);
        for (let i = 0; i < nFree; i++) { curX[freeNodes[i]] = sx[i]; curY[freeNodes[i]] = sy[i]; }
      }

      const cosStep = Math.max(Math.cos((2.0 * Math.PI) / Math.max(L, 6)), 0.68);
      const rPuff = Math.min(1.58, 1.14 / cosStep);
      for (const hd of faces[outerIdx]) {
        for (const nd of [m + hd, m + localOpp[hd]]) {
          if (fixedPos.has(nd)) { curX[nd] *= rPuff; curY[nd] *= rPuff; }
        }
      }
      const subPos = Array.from({ length: 5 * m }, (_, i) => [curX[i], curY[i]]);
      for (let v = 0; v < m; v++) {
        const [vx, vy] = subPos[v];
        const qs = [0, 1, 2, 3].map((s) => subPos[m + 4 * v + s]);
        if (!verifyCcwCyclicOrder(vx, vy, qs)) throw new Error(`CCW cyclic dart order violated at vertex ${v}`);
      }
      return subPos;
    }

    function layoutConnectedBlock(m, localOpp, localVtype, forcedOuterDart = null) {
      if (m >= 2) {
        const visChk = new Array(4 * m).fill(false);
        let fCnt = 0;
        for (let st = 0; st < 4 * m; st++) {
          if (visChk[st]) continue;
          fCnt++;
          let cur = st;
          while (!visChk[cur]) {
            visChk[cur] = true;
            const od = localOpp[cur];
            cur = (od & ~3) | (((od & 3) + 3) & 3);
          }
        }
        if (fCnt !== m + 2) {
          throw new Error(`Non-planar shadow: V=${m}, E=${2 * m}, F=${fCnt} (expected ${m + 2}).`);
        }
      }

      const sl = findSelfLoop(m, localOpp);
      if (sl !== null) {
        const [vSl, s0] = sl;
        const s1 = (s0 + 1) & 3, s2 = (s0 + 2) & 3, s3 = (s0 + 3) & 3;
        const od2 = localOpp[4 * vSl + s2], od3 = localOpp[4 * vSl + s3];
        const remNodes = [];
        const idxR = new Map();
        for (let u = 0; u < m; u++) if (u !== vSl) { idxR.set(u, remNodes.length); remNodes.push(u); }
        const mr = m - 1;
        const oppR = new Array(4 * mr).fill(0);
        const vtR = remNodes.map((u) => localVtype[u]);
        const od2R = 4 * idxR.get(od2 >> 2) + (od2 & 3);
        const od3R = 4 * idxR.get(od3 >> 2) + (od3 & 3);
        for (const u of remNodes) {
          const iu = idxR.get(u);
          for (let s = 0; s < 4; s++) {
            const od = localOpp[4 * u + s];
            if (od === 4 * vSl + s2) oppR[4 * iu + s] = od3R;
            else if (od === 4 * vSl + s3) oppR[4 * iu + s] = od2R;
            else oppR[4 * iu + s] = 4 * idxR.get(od >> 2) + (od & 3);
          }
        }
        let forcedR = null;
        if (forcedOuterDart !== null) {
          const fv = forcedOuterDart >> 2, fs = forcedOuterDart & 3;
          forcedR = fv !== vSl ? (4 * idxR.get(fv) + fs) : ((fs === s1 || fs === s3) ? od2R : od3R);
        }
        const posR = layoutConnectedBlock(mr, oppR, vtR, forcedR);
        let [p2x, p2y] = posR[mr + od2R];
        const [p3x, p3y] = posR[mr + od3R];
        if (Math.hypot(p2x - p3x, p2y - p3y) < 1e-6) {
          const [v2x, v2y] = posR[od2R >> 2];
          p2x = 0.5 * (v2x + p2x); p2y = 0.5 * (v2y + p2y);
        }
        const pvx = 0.5 * (p2x + p3x), pvy = 0.5 * (p2y + p3y);
        const segHalf = 0.5 * Math.hypot(p2x - p3x, p2y - p3y);
        const theta2 = Math.atan2(p2y - pvy, p2x - pvx);
        const thetaN = theta2 - 0.5 * Math.PI;
        const nx = Math.cos(thetaN), ny = Math.sin(thetaN);
        let rLoop = Math.min(0.85 * segHalf, 0.42);
        const prx = pvx + rLoop * nx, pry = pvy + rLoop * ny;
        for (let idxPt = 0; idxPt < 5 * mr; idxPt++) {
          if (idxPt === mr + od2R || idxPt === mr + od3R) continue;
          const [qx, qy] = posR[idxPt];
          rLoop = Math.min(rLoop, 0.45 * Math.hypot(qx - pvx, qy - pvy), 0.6 * Math.hypot(qx - prx, qy - pry));
        }
        rLoop = Math.max(rLoop, 0.12 * segHalf);
        const outV = new Array(m), outU = new Array(4 * m);
        for (const u of remNodes) {
          const iu = idxR.get(u);
          outV[u] = posR[iu];
          for (let s = 0; s < 4; s++) outU[4 * u + s] = posR[mr + 4 * iu + s];
        }
        const v2Pos = posR[od2R >> 2], v3Pos = posR[od3R >> 2];
        outU[od2] = [0.35 * v2Pos[0] + 0.65 * p2x, 0.35 * v2Pos[1] + 0.65 * p2y];
        outU[od3] = [0.35 * v3Pos[0] + 0.65 * p3x, 0.35 * v3Pos[1] + 0.65 * p3y];
        outV[vSl] = [pvx, pvy];
        outU[4 * vSl + s2] = [pvx + 0.45 * (p2x - pvx), pvy + 0.45 * (p2y - pvy)];
        outU[4 * vSl + s3] = [pvx + 0.45 * (p3x - pvx), pvy + 0.45 * (p3y - pvy)];
        const ang0 = theta2 - (2.0 * Math.PI / 3.0), ang1 = theta2 - (1.0 * Math.PI / 3.0);
        outU[4 * vSl + s0] = [pvx + rLoop * Math.cos(ang0), pvy + rLoop * Math.sin(ang0)];
        outU[4 * vSl + s1] = [pvx + rLoop * Math.cos(ang1), pvy + rLoop * Math.sin(ang1)];
        const subPos = [...outV, ...outU];
        for (let v = 0; v < m; v++) {
          const [vx, vy] = subPos[v];
          if (!verifyCcwCyclicOrder(vx, vy, [0, 1, 2, 3].map((s) => subPos[m + 4 * v + s]))) {
            throw new Error(`CCW cyclic dart order violated at vertex ${v}`);
          }
        }
        return subPos;
      }

      const split = findCutSplit(m, localOpp);
      if (split !== null) {
        let [vCut, s0, vertsA, vertsB] = split;
        if (forcedOuterDart !== null) {
          const fv = forcedOuterDart >> 2, fs = forcedOuterDart & 3;
          if (vertsB.includes(fv) || (fv === vCut && (fs === ((s0 + 2) & 3) || fs === ((s0 + 3) & 3)))) {
            s0 = (s0 + 2) & 3;
            [vertsA, vertsB] = [vertsB, vertsA];
          }
        }
        const s1 = (s0 + 1) & 3, s2 = (s0 + 2) & 3, s3 = (s0 + 3) & 3;
        const nodesA = [vCut, ...vertsA];
        const idxA = new Map(nodesA.map((u, i) => [u, i]));
        const ma = nodesA.length;
        const oppA = new Array(4 * ma).fill(0);
        for (const u of nodesA) {
          const iu = idxA.get(u);
          for (let s = 0; s < 4; s++) {
            if (u === vCut && (s === s2 || s === s3)) oppA[4 * iu + s] = 4 * iu + (s === s2 ? s3 : s2);
            else { const od = localOpp[4 * u + s]; oppA[4 * iu + s] = 4 * idxA.get(od >> 2) + (od & 3); }
          }
        }
        const forcedA = forcedOuterDart !== null ? (4 * idxA.get(forcedOuterDart >> 2) + (forcedOuterDart & 3)) : null;
        const posA = layoutConnectedBlock(ma, oppA, nodesA.map((u) => localVtype[u]), forcedA);

        const nodesB = [vCut, ...vertsB];
        const idxB = new Map(nodesB.map((u, i) => [u, i]));
        const mb = nodesB.length;
        const oppB = new Array(4 * mb).fill(0);
        for (const u of nodesB) {
          const iu = idxB.get(u);
          for (let s = 0; s < 4; s++) {
            if (u === vCut && (s === s0 || s === s1)) oppB[4 * iu + s] = 4 * iu + (s === s0 ? s1 : s0);
            else { const od = localOpp[4 * u + s]; oppB[4 * iu + s] = 4 * idxB.get(od >> 2) + (od & 3); }
          }
        }
        const posB = layoutConnectedBlock(mb, oppB, nodesB.map((u) => localVtype[u]), s1);

        const [vaX, vaY] = posA[0];
        const [qa2X, qa2Y] = posA[ma + s2], [qa3X, qa3Y] = posA[ma + s3];
        let pocketR = Math.min(Math.hypot(qa2X - vaX, qa2Y - vaY), Math.hypot(qa3X - vaX, qa3Y - vaY));
        const midPx = 0.5 * (qa2X + qa3X), midPy = 0.5 * (qa2Y + qa3Y);
        for (let idxPt = 0; idxPt < 5 * ma; idxPt++) {
          if (idxPt === 0 || idxPt === ma + s2 || idxPt === ma + s3) continue;
          const [px, py] = posA[idxPt];
          pocketR = Math.min(pocketR, 0.75 * Math.hypot(px - midPx, py - midPy));
        }
        const ang2A = Math.atan2(qa2Y - vaY, qa2X - vaX);
        const angTarget = ang2A + 0.5 * posMod2Pi(Math.atan2(qa3Y - vaY, qa3X - vaX) - ang2A);

        const [vbX, vbY] = posB[0];
        const [qb2X, qb2Y] = posB[mb + s2], [qb3X, qb3Y] = posB[mb + s3];
        const ang2B = Math.atan2(qb2Y - vbY, qb2X - vbX);
        const angSrc = ang2B + 0.5 * posMod2Pi(Math.atan2(qb3Y - vbY, qb3X - vbX) - ang2B);

        let maxDistB = 1e-6;
        for (let i = 1; i < mb; i++) maxDistB = Math.max(maxDistB, Math.hypot(posB[i][0] - vbX, posB[i][1] - vbY));
        for (let iu = 0; iu < mb; iu++) {
          for (let s = 0; s < 4; s++) {
            if (iu === 0 && (s === s0 || s === s1)) continue;
            const [px, py] = posB[mb + 4 * iu + s];
            maxDistB = Math.max(maxDistB, Math.hypot(px - vbX, py - vbY));
          }
        }
        const scale = (0.9 * pocketR) / maxDistB;
        const rot = angTarget - angSrc;
        const cosR = Math.cos(rot), sinR = Math.sin(rot);
        const xformB = ([px, py]) => {
          const dx = px - vbX, dy = py - vbY;
          return [vaX + scale * (dx * cosR - dy * sinR), vaY + scale * (dx * sinR + dy * cosR)];
        };

        const outV = new Array(m), outU = new Array(4 * m);
        for (const u of nodesA) {
          const iu = idxA.get(u);
          outV[u] = posA[iu];
          for (let s = 0; s < 4; s++) {
            if (!(u === vCut && (s === s2 || s === s3))) outU[4 * u + s] = posA[ma + 4 * iu + s];
          }
        }
        for (const u of nodesB) {
          const iu = idxB.get(u);
          if (u !== vCut) outV[u] = xformB(posB[iu]);
          for (let s = 0; s < 4; s++) {
            if (!(u === vCut && (s === s0 || s === s1))) outU[4 * u + s] = xformB(posB[mb + 4 * iu + s]);
          }
        }
        const subPos = [...outV, ...outU];
        for (let v = 0; v < m; v++) {
          const [vx, vy] = subPos[v];
          if (!verifyCcwCyclicOrder(vx, vy, [0, 1, 2, 3].map((s) => subPos[m + 4 * v + s]))) {
            throw new Error(`CCW cyclic dart order violated at vertex ${v}`);
          }
        }
        return subPos;
      }

      return layout2ConnectedLoopFreeBlock(m, localOpp, localVtype, forcedOuterDart);
    }

    function buildSplines(m, localOpp, localVtype, localScomp, subPos, resolution) {
      const vPos = subPos.slice(0, m);
      const uPos = subPos.slice(m, 5 * m);
      const xs = subPos.map((p) => p[0]);
      const ys = subPos.map((p) => p[1]);
      const diagSpan = Math.max(Math.max(...xs) - Math.min(...xs), Math.max(...ys) - Math.min(...ys), 1.0);

      const vAngles = [];
      const vMinArm = [];
      for (let v = 0; v < m; v++) {
        const [vx, vy] = vPos[v];
        let csum = 0, ssum = 0, minD = Infinity;
        for (let s = 0; s < 4; s++) {
          const [ux, uy] = uPos[4 * v + s];
          const dx = ux - vx, dy = uy - vy;
          const d = Math.hypot(dx, dy);
          if (d < minD) minD = d;
          const ang = Math.atan2(dy, dx) - s * 0.5 * Math.PI;
          csum += Math.cos(ang); ssum += Math.sin(ang);
        }
        vAngles.push(Math.atan2(ssum, csum));
        vMinArm.push(minD);
      }

      function slotTangent(cv, cs) {
        const [px, py] = vPos[cv];
        const [qx, qy] = uPos[4 * cv + cs];
        const dx = qx - px, dy = qy - py;
        const natD = Math.hypot(dx, dy);
        const natU = norm2(dx, dy);
        const phi = vAngles[cv] + cs * 0.5 * Math.PI;
        const orthU = [Math.cos(phi), Math.sin(phi)];
        let uv = orthU;
        if (orthU[0] * natU[0] + orthU[1] * natU[1] < 0.65) {
          uv = norm2(0.45 * orthU[0] + 0.55 * natU[0], 0.45 * orthU[1] + 0.55 * natU[1]);
          const oq = uPos[4 * cv + (cs ^ 2)];
          const oppU = norm2(oq[0] - px, oq[1] - py);
          const axU = norm2(natU[0] - oppU[0], natU[1] - oppU[1]);
          if (axU[0] * natU[0] + axU[1] * natU[1] > 0.35) uv = axU;
        }
        return [uv, natD];
      }

      const bzs = [];
      const mks = [];
      for (let hd = 0; hd < 4 * m; hd++) {
        const od = localOpp[hd];
        if (hd >= od) continue;
        const qh = uPos[hd], qo = uPos[od];
        const mid = [0.5 * (qh[0] + qo[0]), 0.5 * (qh[1] + qo[1])];
        const scomp = localScomp[hd];
        for (const [cv, cs, cq] of [[hd >> 2, hd & 3, qh], [od >> 2, od & 3, qo]]) {
          const vt = localVtype[cv];
          const isUnder = (vt === 'X0' && (cs & 1) === 0) || (vt === 'X1' && (cs & 1) === 1);
          const isSm = resolution !== 0 && (vt === 'Y' || vt === 'Z');
          const pc = vPos[cv];
          const [uv, natD] = slotTangent(cv, cs);
          const gapR = Math.min(0.32 * vMinArm[cv], 0.056 * diagSpan);
          let p0, layer, l1;
          if (isUnder) {
            p0 = [pc[0] + gapR * uv[0], pc[1] + gapR * uv[1]];
            layer = 0;
            l1 = Math.min(0.34 * natD, 0.36 * vMinArm[cv]);
          } else if (isSm) {
            const gapS = Math.min(0.42 * vMinArm[cv], 0.072 * diagSpan);
            p0 = [pc[0] + gapS * uv[0], pc[1] + gapS * uv[1]];
            layer = 1;
            l1 = Math.min(0.28 * Math.max(natD - gapS, 0.22 * natD), 0.30 * vMinArm[cv]);
          } else {
            p0 = pc;
            layer = 1;
            l1 = Math.min(0.34 * natD, 0.36 * vMinArm[cv]);
          }
          const p1 = [p0[0] + l1 * uv[0], p0[1] + l1 * uv[1]];
          bzs.push([p0, p1, cq, mid, scomp, layer]);
        }
      }

      for (let v = 0; v < m; v++) {
        const vt = localVtype[v];
        if (vt !== 'Y' && vt !== 'Z') continue;
        const [vx, vy] = vPos[v];
        const [sa, sb] = vt === 'Y' ? [0, 1] : [1, 2];
        const ua = norm2(uPos[4 * v + sa][0] - vx, uPos[4 * v + sa][1] - vy);
        const ub = norm2(uPos[4 * v + sb][0] - vx, uPos[4 * v + sb][1] - vy);
        const bis = norm2(ua[0] + ub[0], ua[1] + ub[1]);
        const br = Math.min(0.32 * vMinArm[v], 0.058 * diagSpan);
        mks.push([[vx - br * bis[0], vy - br * bis[1]], [vx + br * bis[0], vy + br * bis[1]], [vx, vy], vt, localScomp[4 * v]]);
        if (resolution !== 0) {
          const gapS = Math.min(0.42 * vMinArm[v], 0.072 * diagSpan);
          const uvs = [0, 1, 2, 3].map((si) => slotTangent(v, si)[0]);
          const pairs = ((vt === 'Y' && resolution === -1) || (vt === 'Z' && resolution === 1))
            ? [[0, 1], [2, 3]] : [[0, 3], [1, 2]];
          for (const [ia, ib] of pairs) {
            const pa = [vx + gapS * uvs[ia][0], vy + gapS * uvs[ia][1]];
            const pb = [vx + gapS * uvs[ib][0], vy + gapS * uvs[ib][1]];
            const ca = [pa[0] - 0.56 * gapS * uvs[ia][0], pa[1] - 0.56 * gapS * uvs[ia][1]];
            const cb = [pb[0] - 0.56 * gapS * uvs[ib][0], pb[1] - 0.56 * gapS * uvs[ib][1]];
            bzs.push([pa, ca, cb, pb, localScomp[4 * v + ia], 1]);
          }
        }
      }
      return { vPos, bzs, mks };
    }

    function getSplittingBundle(epdStr) {
      if (layoutCache.has(epdStr)) return layoutCache.get(epdStr);
      const diag = parseEpd(epdStr);
      const { n, opp, vtype, freeLoops } = diag;
      const { dartComp, graphComps, totalComps } = dartSurfaceComponents(diag);

      const compOfV = new Array(n).fill(-1);
      const blocks = [];
      for (let sv = 0; sv < n; sv++) {
        if (compOfV[sv] !== -1) continue;
        const bvs = [sv];
        compOfV[sv] = blocks.length;
        const q = [sv];
        while (q.length) {
          const u = q.pop();
          for (let s = 0; s < 4; s++) {
            const w = opp[4 * u + s] >> 2;
            if (compOfV[w] === -1) { compOfV[w] = blocks.length; bvs.push(w); q.push(w); }
          }
        }
        bvs.sort((a, b) => a - b);
        blocks.push(bvs);
      }

      const solved = blocks.map((bvs) => {
        const m = bvs.length;
        const g2l = new Map(bvs.map((gv, lv) => [gv, lv]));
        const lOpp = new Array(4 * m);
        const lVt = bvs.map((gv) => vtype[gv]);
        const lSc = new Array(4 * m);
        for (let lv = 0; lv < m; lv++) {
          const gv = bvs[lv];
          for (let s = 0; s < 4; s++) {
            const od = opp[4 * gv + s];
            lOpp[4 * lv + s] = 4 * g2l.get(od >> 2) + (od & 3);
            lSc[4 * lv + s] = dartComp[4 * gv + s];
          }
        }
        const subPos = layoutConnectedBlock(m, lOpp, lVt);
        return { m, lOpp, lVt, lSc, subPos };
      });

      function assemble(resMode, refBbox) {
        let cursorX = 0;
        const gap = 0.45;
        const allBz = [], allMk = [], allCc = [];
        for (const blk of solved) {
          const sp0 = buildSplines(blk.m, blk.lOpp, blk.lVt, blk.lSc, blk.subPos, 0);
          const sp = resMode === 0 ? sp0 : buildSplines(blk.m, blk.lOpp, blk.lVt, blk.lSc, blk.subPos, resMode);
          const px = [];
          for (const b of sp0.bzs) for (let k = 0; k < 4; k++) px.push(b[k][0]);
          const xmin = Math.min(...px), xmax = Math.max(...px);
          const sh = cursorX - xmin;
          for (const b of sp.bzs) {
            allBz.push([[b[0][0] + sh, b[0][1]], [b[1][0] + sh, b[1][1]], [b[2][0] + sh, b[2][1]], [b[3][0] + sh, b[3][1]], b[4], b[5]]);
          }
          for (const mk of sp.mks) {
            allMk.push([[mk[0][0] + sh, mk[0][1]], [mk[1][0] + sh, mk[1][1]], [mk[2][0] + sh, mk[2][1]], mk[3], mk[4]]);
          }
          cursorX += (xmax - xmin) + gap;
        }
        for (let i = 0; i < freeLoops; i++) {
          const r = 0.72, cx = cursorX + r;
          allCc.push([[cx, 0.0], r, graphComps + i]);
          cursorX += 2 * r + gap;
        }
        let bbox = refBbox;
        if (!bbox) {
          const xs = [], ys = [];
          for (const b of allBz) for (let k = 0; k < 4; k++) { xs.push(b[k][0]); ys.push(b[k][1]); }
          for (const mk of allMk) { xs.push(mk[0][0], mk[1][0]); ys.push(mk[0][1], mk[1][1]); }
          for (const cc of allCc) { xs.push(cc[0][0] - cc[1], cc[0][0] + cc[1]); ys.push(cc[0][1] - cc[1], cc[0][1] + cc[1]); }
          bbox = xs.length ? [Math.min(...xs), Math.min(...ys), Math.max(...xs), Math.max(...ys)] : [-1, -1, 1, 1];
        }
        return { bzs: allBz, mks: allMk, ccs: allCc, bbox };
      }

      const zero = assemble(0, null);
      const minus = assemble(-1, zero.bbox);
      const plus = assemble(1, zero.bbox);
      const cCnt = vtype.filter((t) => t === 'X0' || t === 'X1').length;
      const mCnt = vtype.filter((t) => t === 'Y' || t === 'Z').length;
      const cm = countResolutionCircles(diag, -1);
      const cp = countResolutionCircles(diag, 1);
      const bundle = { minus, zero, plus, bbox: zero.bbox, totalComps, n, cCnt, mCnt, cm, cp, chi: cm - mCnt + cp };
      layoutCache.set(epdStr, bundle);
      return bundle;
    }

    function renderTikzFromZero(zero, label) {
      const [xmin, ymin, xmax, ymax] = zero.bbox;
      const span = Math.max(xmax - xmin, ymax - ymin, 1e-3);
      const cx = 0.5 * (xmin + xmax), cy = 0.5 * (ymin + ymax);
      const sc = 5.4 / span;
      const tx = (p) => [((p[0] - cx) * sc).toFixed(4), ((p[1] - cy) * sc).toFixed(4)];
      const used = new Set([...zero.bzs.map((b) => b[4]), ...zero.ccs.map((c) => c[2])]);
      const lines = [`% Canonical ch-diagram for ${label}`, `\\begin{tikzpicture}[line join=round]`];
      [...used].sort((a, b) => a - b).forEach((ci) => {
        const [r, g, b] = TIKZ_RGB[ci % TIKZ_RGB.length];
        lines.push(`  \\definecolor{skcomp${ci}}{RGB}{${r},${g},${b}}`);
      });
      lines.push(`  \\definecolor{skmark}{RGB}{15,23,42}`);
      for (const cc of zero.ccs) {
        const [x, y] = tx(cc[0]);
        lines.push(`  \\draw[skcomp${cc[2]}, line width=1.15pt] (${x},${y}) circle (${(cc[1] * sc).toFixed(4)});`);
      }
      for (const [layerIdx, cap] of [[0, 'butt'], [1, 'round']]) {
        for (const b of zero.bzs) {
          if (b[5] !== layerIdx) continue;
          const p0 = tx(b[0]), p1 = tx(b[1]), p2 = tx(b[2]), p3 = tx(b[3]);
          lines.push(`  \\draw[skcomp${b[4]}, line width=1.15pt, line cap=${cap}] (${p0[0]},${p0[1]}) .. controls (${p1[0]},${p1[1]}) and (${p2[0]},${p2[1]}) .. (${p3[0]},${p3[1]});`);
        }
      }
      for (const mk of zero.mks) {
        const pa = tx(mk[0]), pb = tx(mk[1]), pc = tx(mk[2]);
        lines.push(`  \\draw[white, line width=2.8pt, line cap=round] (${pa[0]},${pa[1]}) -- (${pb[0]},${pb[1]});`);
        lines.push(`  \\draw[skmark, line width=1.85pt, line cap=round] (${pa[0]},${pa[1]}) -- (${pb[0]},${pb[1]});`);
        lines.push(`  \\fill[skmark] (${pc[0]},${pc[1]}) circle (1.2pt);`);
      }
      lines.push(`\\end{tikzpicture}`);
      return lines.join('\n');
    }

    /* ---------- Theme & SVG Rendering ---------- */
    function isLightTheme() {
      const root = document.documentElement;
      if (root.classList.contains('light')) return true;
      if (root.classList.contains('dark')) return false;
      const bg = getComputedStyle(document.body).backgroundColor || '';
      const m = bg.match(/(\d+)\s*,\s*(\d+)\s*,\s*(\d+)/);
      if (!m) return true;
      return (0.299 * Number(m[1]) + 0.587 * Number(m[2]) + 0.114 * Number(m[3])) / 255 > 0.5;
    }

    function getPalette() {
      return isLightTheme() ? PALETTE_LIGHT : PALETTE_DARK;
    }

    function esc(s) {
      return String(s ?? '')
        .replace(/&/g, '&amp;')
        .replace(/</g, '&lt;')
        .replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;');
    }

    function fmtCh(e) {
      return (e.ex_ch || e.ch_lo === e.ch_hi) ? String(e.ch_lo) : `[${e.ch_lo}, ${e.ch_hi}]`;
    }

    function renderSliceSvgInner(slice, bbox, resMode, strokeW) {
      const [xmin, ymin, xmax, ymax] = bbox;
      const span = Math.max(xmax - xmin, ymax - ymin, 1e-3);
      const cx = 0.5 * (xmin + xmax), cy = 0.5 * (ymin + ymax);
      const scale = (300 - 44) / span;
      const tr = (p) => [(150 + (p[0] - cx) * scale).toFixed(2), (150 - (p[1] - cy) * scale).toFixed(2)];
      const pal = getPalette();
      const parts = [];

      if (resMode !== 0) {
        for (const mk of slice.mks) {
          const pa = tr(mk[0]), pb = tr(mk[1]);
          parts.push(`<line x1="${pa[0]}" y1="${pa[1]}" x2="${pb[0]}" y2="${pb[1]}" stroke="var(--muted-foreground)" stroke-opacity="0.55" stroke-width="1.4" stroke-dasharray="2.5,2.5" stroke-linecap="round"/>`);
        }
      }
      for (const cc of slice.ccs) {
        const pc = tr(cc[0]);
        const col = pal[cc[2] % pal.length];
        parts.push(`<circle cx="${pc[0]}" cy="${pc[1]}" r="${(cc[1] * scale).toFixed(2)}" fill="none" stroke="${col}" stroke-width="${strokeW}"/>`);
      }
      for (const [layerIdx, cap] of [[0, 'butt'], [1, 'round']]) {
        for (const b of slice.bzs) {
          if (b[5] !== layerIdx) continue;
          const p0 = tr(b[0]), p1 = tr(b[1]), p2 = tr(b[2]), p3 = tr(b[3]);
          const col = pal[b[4] % pal.length];
          parts.push(`<path d="M ${p0[0]} ${p0[1]} C ${p1[0]} ${p1[1]}, ${p2[0]} ${p2[1]}, ${p3[0]} ${p3[1]}" fill="none" stroke="${col}" stroke-width="${strokeW}" stroke-linecap="${cap}"/>`);
        }
      }
      if (resMode === 0) {
        for (const mk of slice.mks) {
          const pa = tr(mk[0]), pb = tr(mk[1]), pc = tr(mk[2]);
          parts.push(`<line x1="${pa[0]}" y1="${pa[1]}" x2="${pb[0]}" y2="${pb[1]}" stroke="var(--background)" stroke-width="${(strokeW * 2.45).toFixed(2)}" stroke-linecap="round"/>`);
          parts.push(`<line x1="${pa[0]}" y1="${pa[1]}" x2="${pb[0]}" y2="${pb[1]}" stroke="var(--foreground)" stroke-width="${(strokeW * 1.65).toFixed(2)}" stroke-linecap="round"/>`);
          parts.push(`<circle cx="${pc[0]}" cy="${pc[1]}" r="${(strokeW * 0.72).toFixed(2)}" fill="var(--foreground)"/>`);
        }
      }
      return parts.join('');
    }

    function renderSingleSvgFromEpd(epd, label, maxPx) {
      const b = getSplittingBundle(epd);
      return `<svg viewBox="0 0 300 300" class="w-full mx-auto block" style="max-width:${maxPx}px" role="img" aria-label="${esc(label)} diagram">${renderSliceSvgInner(b.zero, b.bbox, 0, 2.85)}</svg>`;
    }

    function renderMovieTriptychFromEpd(epd) {
      const b = getSplittingBundle(epd);
      const panels = [
        { slice: b.minus, mode: -1, title: `t = -\u03b5 : L\u208b(D)`, meta: `${b.cm} circle${b.cm === 1 ? '' : 's'}` },
        { slice: b.zero,  mode: 0,  title: `t = 0 : Marked D`,      meta: `c=${b.cCnt}, m=${b.mCnt}` },
        { slice: b.plus,  mode: 1,  title: `t = +\u03b5 : L\u208a(D)`, meta: `${b.cp} circle${b.cp === 1 ? '' : 's'}` },
      ];
      const cards = panels.map((p) => (
        `<div class="bg-[var(--background)] border border-[var(--border)] rounded-lg p-2 flex flex-col items-center">` +
          `<div class="w-full flex items-center justify-between text-xs font-mono px-1 pb-1 border-b border-[var(--card-border)]">` +
            `<span class="font-semibold text-[var(--foreground)]">${esc(p.title)}</span>` +
            `<span class="text-[var(--muted-foreground)]">${esc(p.meta)}</span>` +
          `</div>` +
          `<svg viewBox="0 0 300 300" class="w-full max-w-[225px] block mt-1">${renderSliceSvgInner(p.slice, b.bbox, p.mode, 2.75)}</svg>` +
        `</div>`
      )).join('');
      return (
        `<div class="space-y-2">` +
          `<div class="grid grid-cols-1 sm:grid-cols-3 gap-2.5">${cards}</div>` +
          `<div class="text-center text-xs font-mono text-[var(--muted-foreground)] bg-[var(--background)] border border-[var(--border)] rounded-lg py-1.5 px-2">` +
            `Hyperbolic splitting Euler balance: \u03c7(F) = c\u208b &minus; m + c\u208a = ${b.cm} &minus; ${b.mCnt} + ${b.cp} = ${b.chi}` +
          `</div>` +
        `</div>`
      );
    }

    /* ---------- UI State & Views ---------- */
    function resolveLabel(lbl) {
      if (BY_LABEL[lbl]) return lbl;
      const found = ENTRIES.find((e) => e.reals.some((r) => r.id === lbl));
      return found ? found.label : (ENTRIES[0] ? ENTRIES[0].label : lbl);
    }

    function rigorPillHtml(badgeText, compact = false) {
      if (!badgeText) return '';
      const isWarn = badgeText.startsWith('HEURISTIC') || badgeText.startsWith('PACKET');
      const cls = isWarn
        ? 'border-[var(--warning)] text-[var(--warning)] bg-[var(--background)]'
        : badgeText.startsWith('PROVABLE')
          ? 'border-[var(--border)] text-[var(--foreground)] bg-[var(--secondary)]'
          : 'border-[var(--border)] text-[var(--muted-foreground)] bg-[var(--background)]';
      const pad = compact ? 'px-1.5 py-0.5' : 'px-2 py-0.5';
      return `<span class="${pad} rounded text-xs font-mono border ${cls}">${esc(badgeText)}</span>`;
    }

    let state = {
      mainTab: 'explorer',
      selectedLabel: '10_1^0,1',
      diagramMode: 'zero',
      cmpA: resolveLabel('R(b^2*a)'),
      cmpB: resolveLabel('R(b^2*a^-1)'),
    };

    function matchesFilters(e) {
      const q = document.getElementById('f-search').value.trim().toLowerCase();
      const tier = document.getElementById('f-tier').value;
      const topo = document.getElementById('f-topo').value;
      const ori = document.getElementById('f-ori').value;
      const prop = document.getElementById('f-prop').value;

      const isHub = ['yoshikawa_ch10', 'expanded_ch11', 'expanded_ch12'].includes(e.tier_code);
      if (tier === 'hub' && !isHub) return false;
      if (tier === 'spoke' && isHub) return false;
      if (['yoshikawa_ch10', 'expanded_ch11', 'expanded_ch12'].includes(tier) && e.tier_code !== tier) return false;

      if (topo === 's2_knot' && !(e.nc === 1 && e.types[0] === 'S^2')) return false;
      if (topo === 't2_knot' && !(e.nc === 1 && e.types[0] === 'T^2')) return false;
      if (topo === 'nonori_knot' && !(e.nc === 1 && e.ori === 'Non')) return false;
      if (topo === 'link_2' && e.nc !== 2) return false;
      if (topo === 'link_3' && e.nc !== 3) return false;
      if (topo === 'brunnian' && !e.brunn) return false;

      if (ori !== 'all' && e.ori !== ori) return false;

      if (prop === 'exact_ch' && !e.ex_ch) return false;
      if (prop === 'rigor_provable' && e.rigor_code !== 'provable') return false;
      if (prop === 'rigor_packet' && e.rigor_code !== 'packet') return false;
      if (prop === 'rigor_heuristic' && !(e.rigor_code === 'heuristic' || (e.spoke_match_badge && e.spoke_match_badge.startsWith('HEURISTIC')))) return false;
      if (prop === 'rigor_bounded' && e.rigor_code !== 'bounded') return false;
      if (prop === 'hub_realized' && e.reals.length === 0) return false;
      if (prop === 'chiral' && !e.chiral) return false;
      if (prop === 'q6_triple' && !e.q6_pos) return false;
      if (prop === 'erratum' && !e.err) return false;

      if (q) {
        const hay = [
          e.label, e.desc, e.types.join(' '), e.spoke, e.subfamily,
          e.rigor_badge, e.spoke_match_badge, e.nw,
          e.e1.join(' '), e.d1, e.pi1, e.reals.map((r) => `${r.id} ${r.rigor || ''}`).join(' '),
        ].join(' ').toLowerCase();
        if (!hay.includes(q)) return false;
      }
      return true;
    }

    function renderCatalogList() {
      const filtered = ENTRIES.filter(matchesFilters);
      document.getElementById('filter-count').textContent = `${filtered.length} of ${ENTRIES.length} surface-links`;
      if (filtered.length > 0 && !filtered.some((e) => e.label === state.selectedLabel)) {
        state.selectedLabel = filtered[0].label;
      }
      const listEl = document.getElementById('catalog-list');
      if (filtered.length === 0) {
        listEl.innerHTML = `<div class="p-6 text-center text-xs text-[var(--muted-foreground)]">No surface-links match the active filters.</div>`;
        document.getElementById('detail-panel').innerHTML = `<div class="p-6 text-center text-xs text-[var(--muted-foreground)]">Reset filters to inspect a surface-link.</div>`;
        return;
      }
      listEl.innerHTML = filtered.map((e) => {
        const active = e.label === state.selectedLabel;
        const rowCls = active ? 'bg-[var(--secondary)]' : 'hover:bg-[var(--background)]';
        const badge = e.ex_ch
          ? `<span class="px-1.5 py-0.5 rounded text-xs font-mono bg-[var(--background)] border border-[var(--border)]">ch=${fmtCh(e)}</span>`
          : `<span class="px-1.5 py-0.5 rounded text-xs font-mono text-[var(--muted-foreground)] border border-[var(--border)]">ch&isin;${fmtCh(e)}</span>`;
        const spokeWarn = (e.spoke_match_badge && e.spoke_match_badge.startsWith('HEURISTIC'))
          ? ` ${rigorPillHtml('Spoke Join: ' + e.spoke_match_badge, true)}`
          : '';
        return (
          `<button type="button" data-select-label="${esc(e.label)}" class="w-full text-left px-3 py-2.5 transition flex items-start justify-between gap-2 ${rowCls}">` +
            `<div class="min-w-0">` +
              `<div class="flex items-center gap-2">` +
                `<span class="font-mono text-xs font-bold text-[var(--foreground)] truncate">${esc(e.label)}</span>` +
                `<span class="text-xs text-[var(--muted-foreground)] truncate">${esc(e.types.join(' \u222a '))}</span>` +
              `</div>` +
              `<div class="text-xs text-[var(--muted-foreground)] truncate mt-0.5">${esc(e.desc || e.tier_label)}</div>` +
              `<div class="flex flex-wrap gap-1.5 mt-1">${rigorPillHtml(e.rigor_badge, true)}${spokeWarn}</div>` +
            `</div>` +
            `<div class="shrink-0">${badge}</div>` +
          `</button>`
        );
      }).join('');
      listEl.querySelectorAll('[data-select-label]').forEach((btn) => {
        btn.addEventListener('click', () => {
          state.selectedLabel = btn.getAttribute('data-select-label');
          renderCatalogList();
        });
      });
      renderDetailPanel(BY_LABEL[state.selectedLabel]);
    }

    function renderDetailPanel(e) {
      if (!e) return;
      const pal = getPalette();
      const compDots = e.types.map((ct, idx) => {
        const col = pal[idx % pal.length];
        const eul = e.eul[idx];
        return (
          `<span class="inline-flex items-center gap-1.5 px-2 py-0.5 rounded-md text-xs font-mono bg-[var(--background)] border border-[var(--border)]">` +
            `<span class="w-2.5 h-2.5 rounded-full inline-block" style="background:${col}"></span>` +
            `${esc(ct)} (e=${eul >= 0 ? '+' + eul : eul})` +
          `</span>`
        );
      }).join('');

      let stageHtml = '';
      if (state.diagramMode === 'zero') {
        stageHtml = `<div class="bg-[var(--background)] border border-[var(--border)] rounded-xl p-3">${renderSingleSvgFromEpd(e.epd, e.label, 320)}</div>`;
      } else if (state.diagramMode === 'movie') {
        stageHtml = renderMovieTriptychFromEpd(e.epd);
      } else {
        const tikzCode = renderTikzFromZero(getSplittingBundle(e.epd).zero, e.label);
        stageHtml = (
          `<div class="space-y-2.5 text-xs">` +
            `<div>` +
              `<div class="flex items-center justify-between mb-1">` +
                `<span class="font-semibold text-[var(--muted-foreground)]">Canonical Planar EPD (${e.ch_d} vertices)</span>` +
                `<button type="button" data-copy-text="${esc(e.epd)}" class="text-xs text-[var(--link)] hover:underline">Copy EPD</button>` +
              `</div>` +
              `<pre class="p-2.5 rounded-lg bg-[var(--background)] border border-[var(--border)] font-mono text-xs overflow-x-auto whitespace-pre-wrap break-all">${esc(e.epd)}</pre>` +
            `</div>` +
            `<div>` +
              `<div class="flex items-center justify-between mb-1">` +
                `<span class="font-semibold text-[var(--muted-foreground)]">LaTeX / TikZ Snippet</span>` +
                `<button type="button" data-copy-text="${esc(tikzCode)}" class="text-xs text-[var(--link)] hover:underline">Copy TikZ</button>` +
              `</div>` +
              `<pre class="p-2.5 rounded-lg bg-[var(--background)] border border-[var(--border)] font-mono text-xs max-h-56 overflow-y-auto">${esc(tikzCode)}</pre>` +
            `</div>` +
          `</div>`
        );
      }

      const errBanner = e.err ? (
        `<div class="p-2.5 rounded-lg border border-[var(--warning)] bg-[var(--background)] text-xs space-y-0.5">` +
          `<div class="font-bold text-[var(--warning)]">Verified Yoshikawa Table II Erratum</div>` +
          `<div class="text-[var(--foreground)]">${esc(e.err_note)} (Printed: <code class="font-mono">${esc(e.err_print.join(', '))}</code>)</div>` +
        `</div>`
      ) : '';

      const hubRealHtml = e.reals.length > 0 ? (
        `<div class="p-2.5 rounded-lg bg-[var(--background)] border border-[var(--border)] text-xs space-y-1.5">` +
          `<div class="font-semibold text-[var(--foreground)]">Constructive Parametric Spoke Realizations (${e.reals.length})</div>` +
          `<div class="flex flex-wrap gap-1.5">` +
            e.reals.map((r) => {
              const rWarn = (r.rigor || '').startsWith('HEURISTIC');
              const rCls = rWarn ? 'border-[var(--warning)]' : 'border-[var(--border)]';
              return (
                `<span class="px-2 py-0.5 rounded font-mono text-xs bg-[var(--card)] border ${rCls}">` +
                  `${esc(r.id)} &middot; <span class="${rWarn ? 'text-[var(--warning)]' : 'text-[var(--muted-foreground)]'}">${esc(r.rigor || r.mode)} (diag ch=${r.ch})</span>` +
                `</span>`
              );
            }).join('') +
          `</div>` +
        `</div>`
      ) : '';

      const modeBtn = (id, lbl) => {
        const cls = state.diagramMode === id
          ? 'bg-[var(--primary)] text-[var(--primary-foreground)] border-transparent'
          : 'bg-[var(--background)] text-[var(--foreground)] border-[var(--border)]';
        return `<button type="button" data-diag-mode="${id}" class="px-2.5 py-1 rounded-md text-xs font-medium border transition ${cls}">${lbl}</button>`;
      };

      const el = document.getElementById('detail-panel');
      el.innerHTML = (
        `<div class="flex flex-wrap items-start justify-between gap-2 pb-3 border-b border-[var(--card-border)]">` +
          `<div class="space-y-1">` +
            `<div class="flex flex-wrap items-center gap-2">` +
              `<h2 class="text-base font-bold font-mono text-[var(--foreground)]">${esc(e.label)}</h2>` +
              `<span class="px-2 py-0.5 rounded text-xs bg-[var(--secondary)] text-[var(--secondary-foreground)]">${esc(e.tier_label)}</span>` +
              `${rigorPillHtml(e.rigor_badge)}` +
            `</div>` +
            `<p class="text-xs text-[var(--muted-foreground)]">${esc(e.desc)}</p>` +
            `<div class="flex flex-wrap gap-1.5 pt-1">${compDots}</div>` +
          `</div>` +
          `<div class="flex flex-wrap gap-1.5">` +
            modeBtn('zero', 'Marked Diagram D') +
            modeBtn('movie', '3-Slice Movie') +
            modeBtn('code', 'EPD & TikZ') +
          `</div>` +
        `</div>` +
        errBanner +
        stageHtml +
        hubRealHtml +
        `<div class="grid grid-cols-1 md:grid-cols-2 gap-2.5 text-xs">` +
          `<div class="p-2.5 rounded-lg bg-[var(--background)] border border-[var(--border)] space-y-1">` +
            `<div class="font-semibold text-[var(--muted-foreground)]">Complexity &amp; Epistemic Status</div>` +
            `<div><span class="text-[var(--muted-foreground)]">ch(F):</span> <span class="font-mono font-semibold">${e.ex_ch ? e.ch_lo + ' (EXACT)' : '[' + e.ch_lo + ', ' + e.ch_hi + ']'}</span> (diag=${e.ch_d}: c=${e.c}, m=${e.m})</div>` +
            `<div><span class="text-[var(--muted-foreground)]">Epistemic Rigor:</span> <span class="font-mono">${esc(e.rigor_badge)}</span></div>` +
            `<div><span class="text-[var(--muted-foreground)]">Orientability / &chi;(F):</span> <span class="font-mono">${esc(e.ori)} / &chi;=${e.chi}</span> &middot; Ribbon: <span class="font-mono">${esc(e.rib)}</span></div>` +
            `<div><span class="text-[var(--muted-foreground)]">Notes:</span> ${esc(e.reason)}${e.nw ? ' &middot; ' + esc(e.nw) : ''}</div>` +
          `</div>` +
          `<div class="p-2.5 rounded-lg bg-[var(--background)] border border-[var(--border)] space-y-1">` +
            `<div class="font-semibold text-[var(--muted-foreground)]">Knot Group &amp; Alexander Ideals</div>` +
            `<div class="truncate" title="${esc(e.pi1)}"><span class="text-[var(--muted-foreground)]">&pi;&#8321;:</span> <code class="font-mono">${esc(e.pi1)}</code></div>` +
            `<div><span class="text-[var(--muted-foreground)]">H&#8321;:</span> <span class="font-mono">${esc(e.ld[1])}</span></div>` +
            `<div><span class="text-[var(--muted-foreground)]">Ideals:</span> <code class="font-mono">${esc(e.ld[2])}</code></div>` +
          `</div>` +
          `<div class="p-2.5 rounded-lg bg-[var(--background)] border border-[var(--border)] space-y-1">` +
            `<div class="font-semibold text-[var(--muted-foreground)]">Tier-2 Quandle, Cocycle &amp; Quotients</div>` +
            `<div class="font-mono text-xs break-words">${esc(e.ld[3])}</div>` +
            `<div><code class="font-mono">${esc(e.ld[4])}</code></div>` +
            `<div><code class="font-mono">${esc(e.ld[5])}</code></div>` +
          `</div>` +
          `<div class="p-2.5 rounded-lg bg-[var(--background)] border border-[var(--border)] space-y-1">` +
            `<div class="font-semibold text-[var(--muted-foreground)]">Tier-3 Branched Cover &Sigma;&#8322;, &pi;&#8322; &amp; Symmetry</div>` +
            `<div><span class="text-[var(--muted-foreground)]">&Sigma;&#8322;(F):</span> <code class="font-mono">${esc(e.ld[7])}</code>${e.q_hs4 ? ' <span class="px-1 rounded bg-[var(--secondary)] text-xs">Q-HS&#8308;</span>' : ''}</div>` +
            `<div><span class="text-[var(--muted-foreground)]">&pi;&#8322; / Covers:</span> <code class="font-mono">${esc(e.ld[8])}</code></div>` +
            `<div><span class="text-[var(--muted-foreground)]">Symmetry:</span> ${esc(e.ld[9])}</div>` +
          `</div>` +
        `</div>` +
        `<div class="flex flex-wrap items-center justify-end gap-2 pt-1">` +
          `<button type="button" id="btn-load-sandbox" class="px-2.5 py-1 rounded text-xs border border-[var(--border)] bg-[var(--background)] hover:bg-[var(--secondary)]">Open in EPD Sandbox</button>` +
          `<button type="button" id="btn-send-cmp-a" class="px-2.5 py-1 rounded text-xs border border-[var(--border)] bg-[var(--background)] hover:bg-[var(--secondary)]">Set as Compare A</button>` +
          `<button type="button" id="btn-send-cmp-b" class="px-2.5 py-1 rounded text-xs border border-[var(--border)] bg-[var(--background)] hover:bg-[var(--secondary)]">Set as Compare B</button>` +
        `</div>`
      );

      el.querySelectorAll('[data-diag-mode]').forEach((b) => {
        b.addEventListener('click', () => { state.diagramMode = b.getAttribute('data-diag-mode'); renderDetailPanel(e); });
      });
      el.querySelectorAll('[data-copy-text]').forEach((b) => {
        b.addEventListener('click', () => {
          const txt = b.getAttribute('data-copy-text');
          if (navigator.clipboard) navigator.clipboard.writeText(txt);
          const prev = b.textContent;
          b.textContent = 'Copied!';
          setTimeout(() => { b.textContent = prev; }, 1200);
        });
      });
      document.getElementById('btn-load-sandbox').addEventListener('click', () => {
        document.getElementById('sandbox-epd').value = e.epd;
        setMainTab('sandbox');
      });
      document.getElementById('btn-send-cmp-a').addEventListener('click', () => {
        state.cmpA = e.label; document.getElementById('cmp-a').value = e.label; setMainTab('compare');
      });
      document.getElementById('btn-send-cmp-b').addEventListener('click', () => {
        state.cmpB = e.label; document.getElementById('cmp-b').value = e.label; setMainTab('compare');
      });
    }

    const LAYER_NAMES = [
      '1. Intrinsic O(4) Topology & Normal Euler',
      '2. Abelianized H\u2081(S\u2074 \\ F)',
      '3. Alexander Ideal E\u2081 & Twisted S\u2083',
      '4. Finite & SL(2,F\u2083) Quotients',
      '5. Kei, Bikei & Alexander Biquandle',
      '6. Mochizuki 3-Cocycle & Q\u2086 Triple-Points',
      '7. O(4) Mod-2 Double/Triple Linking & Sublinks',
      '8. 4D Double Branched Cover \u03a3\u2082(F)',
      '9. Finite Cover Homology Profile',
      '10. Exact O(4) Chirality & Invertibility Obstructions',
    ];

    function renderCompareView() {
      const a = BY_LABEL[state.cmpA] || ENTRIES[0];
      const b = BY_LABEL[state.cmpB] || ENTRIES[1];
      let sepCount = 0;

      const rowsHtml = LAYER_NAMES.map((name, idx) => {
        const diff = a.lh[idx] !== b.lh[idx];
        if (diff) sepCount += 1;
        const badge = diff
          ? `<span class="px-2 py-0.5 rounded text-xs font-semibold bg-[var(--primary)] text-[var(--primary-foreground)]">SEPARATES</span>`
          : `<span class="px-2 py-0.5 rounded text-xs text-[var(--muted-foreground)] border border-[var(--border)]">MATCH</span>`;
        return (
          `<tr class="border-t border-[var(--card-border)] ${diff ? 'bg-[var(--secondary)]' : ''}">` +
            `<td class="px-3 py-2 font-medium text-xs">${esc(name)}</td>` +
            `<td class="px-3 py-2 font-mono text-xs break-words">${esc(a.ld[idx])}</td>` +
            `<td class="px-3 py-2 font-mono text-xs break-words">${esc(b.ld[idx])}</td>` +
            `<td class="px-3 py-2 text-right">${badge}</td>` +
          `</tr>`
        );
      }).join('');

      const verdict = a.label === b.label
        ? `Identical entry selected (${esc(a.label)}).`
        : sepCount > 0
          ? `Strictly separated by ${sepCount} of 10 exact invariant layers (${esc(a.label)} \u2260 ${esc(b.label)} in R\u2074).`
          : `All 10 invariant layers agree on ${esc(a.label)} and ${esc(b.label)}.`;

      document.getElementById('compare-output').innerHTML = (
        `<div class="grid grid-cols-1 md:grid-cols-2 gap-4">` +
          [a, b].map((x) => (
            `<div class="bg-[var(--card)] border border-[var(--border)] rounded-xl p-3 space-y-2">` +
              `<div class="flex items-center justify-between text-xs border-b border-[var(--card-border)] pb-2">` +
                `<span class="font-mono font-bold text-sm">${esc(x.label)}</span>` +
                `<span class="text-[var(--muted-foreground)]">${esc(x.types.join(' \u222a '))} &middot; ch=${fmtCh(x)}</span>` +
              `</div>` +
              `<div class="bg-[var(--background)] border border-[var(--border)] rounded-lg p-2">${renderSingleSvgFromEpd(x.epd, x.label, 230)}</div>` +
              `<div class="text-xs text-[var(--muted-foreground)] truncate">${esc(x.desc)}</div>` +
            `</div>`
          )).join('') +
        `</div>` +
        `<div class="bg-[var(--card)] border border-[var(--border)] rounded-xl overflow-hidden">` +
          `<div class="px-4 py-2.5 border-b border-[var(--card-border)] flex items-center justify-between text-xs">` +
            `<span class="font-bold text-[var(--foreground)]">${esc(verdict)}</span>` +
            `<span class="font-mono text-[var(--muted-foreground)]">${sepCount}/10 layers separating</span>` +
          `</div>` +
          `<div class="overflow-x-auto">` +
            `<table class="w-full text-left border-collapse">` +
              `<thead>` +
                `<tr class="text-xs text-[var(--muted-foreground)] bg-[var(--background)]">` +
                  `<th class="px-3 py-2">Filtration Layer</th>` +
                  `<th class="px-3 py-2 font-mono">${esc(a.label)}</th>` +
                  `<th class="px-3 py-2 font-mono">${esc(b.label)}</th>` +
                  `<th class="px-3 py-2 text-right">Status</th>` +
                `</tr>` +
              `</thead>` +
              `<tbody>${rowsHtml}</tbody>` +
            `</table>` +
          `</div>` +
        `</div>`
      );
    }

    function renderSandboxView() {
      const out = document.getElementById('sandbox-output');
      const raw = document.getElementById('sandbox-epd').value;
      try {
        const b = getSplittingBundle(raw);
        const tikz = renderTikzFromZero(b.zero, 'custom_epd');
        out.innerHTML = (
          `<div class="space-y-4">` +
            `<div class="flex flex-wrap items-center justify-between gap-2 text-xs pb-2 border-b border-[var(--card-border)]">` +
              `<span class="font-semibold">Planar Diagram Verified: ch_diag = ${b.n} (c=${b.cCnt} classical, m=${b.mCnt} marked) &middot; Surface components = ${b.totalComps} &middot; &chi;(F) = ${b.chi}</span>` +
            `</div>` +
            renderMovieTriptychFromEpd(raw) +
            `<details class="text-xs">` +
              `<summary class="cursor-pointer font-semibold text-[var(--link)]">Show Generated LaTeX / TikZ Source</summary>` +
              `<pre class="mt-2 p-2.5 rounded-lg bg-[var(--background)] border border-[var(--border)] font-mono text-xs max-h-52 overflow-y-auto">${esc(tikz)}</pre>` +
            `</details>` +
          `</div>`
        );
      } catch (err) {
        out.innerHTML = `<div class="text-xs text-[var(--error)] font-mono">EPD Parse / Planarity Error: ${esc(err.message)}</div>`;
      }
    }

    function renderOverviewView() {
      const c = DATA.counts;
      const ch12List = ENTRIES.filter((e) => e.tier_code === 'expanded_ch12' || e.tier_code === 'exact_ch12');
      const hubRealized = ENTRIES.filter((e) => e.reals.length > 0);
      document.getElementById('view-overview').innerHTML = (
        `<div class="grid grid-cols-2 lg:grid-cols-4 gap-3">` +
          `<div class="bg-[var(--card)] border border-[var(--border)] rounded-xl p-3.5">` +
            `<div class="text-xs text-[var(--muted-foreground)]">Unified Distinct Classes</div>` +
            `<div class="text-xl font-bold font-mono mt-1">${c.total_distinct_classes}</div>` +
            `<div class="text-xs text-[var(--muted-foreground)] mt-1">${c.hub_total} Hub + ${c.spoke_total} Spoke</div>` +
          `</div>` +
          `<div class="bg-[var(--card)] border border-[var(--border)] rounded-xl p-3.5">` +
            `<div class="text-xs text-[var(--muted-foreground)]">Exact ch(F) Certified</div>` +
            `<div class="text-xl font-bold font-mono mt-1">${c.exact_ch_total}</div>` +
            `<div class="text-xs text-[var(--muted-foreground)] mt-1">${c.yoshikawa_ch10} (ch&le;10) + ${c.expanded_ch11} (ch=11) + ${c.expanded_ch12 || c.exact_ch12} (ch=12)</div>` +
          `</div>` +
          `<div class="bg-[var(--card)] border border-[var(--border)] rounded-xl p-3.5">` +
            `<div class="text-xs text-[var(--muted-foreground)]">Hub Classes Realized by Spokes</div>` +
            `<div class="text-xl font-bold font-mono mt-1">${c.hub_realized_count} / ${c.hub_total}</div>` +
            `<div class="text-xs text-[var(--muted-foreground)] mt-1">Spoke parametric realizations</div>` +
          `</div>` +
          `<div class="bg-[var(--card)] border border-[var(--border)] rounded-xl p-3.5">` +
            `<div class="text-xs text-[var(--muted-foreground)]">Compiled Family Diagrams</div>` +
            `<div class="text-xl font-bold font-mono mt-1">${c.compiled_family_instances}</div>` +
            `<div class="text-xs text-[var(--muted-foreground)] mt-1">Spokes 1–3 constructive compilers</div>` +
          `</div>` +
        `</div>` +
        `<div class="grid grid-cols-1 lg:grid-cols-2 gap-4">` +
          `<div class="bg-[var(--card)] border border-[var(--border)] rounded-xl p-4 space-y-2.5">` +
            `<h3 class="text-sm font-bold">${ch12List.length} Certified Exact ch(F) = 12 Surface-Links</h3>` +
            `<p class="text-xs text-[var(--muted-foreground)]">Exhaustive ch = 12 prime census classes proven disjoint from the complete ch &le; 11 prime + composite + split reference catalog:</p>` +
            `<div class="max-h-72 overflow-y-auto divide-y divide-[var(--card-border)]">` +
              ch12List.map((e) => (
                `<button type="button" data-jump="${esc(e.label)}" class="w-full py-2 text-left flex items-center justify-between hover:bg-[var(--background)] px-2 rounded text-xs">` +
                  `<div class="min-w-0 truncate"><span class="font-mono font-bold">${esc(e.label)}</span> &middot; <span class="text-[var(--muted-foreground)]">${esc(e.desc || e.ld[2])}</span></div>` +
                  `<span class="font-mono text-xs px-1.5 py-0.5 rounded border border-[var(--border)] shrink-0">${esc(e.types.join(' \u222a '))}</span>` +
                `</button>`
              )).join('') +
            `</div>` +
          `</div>` +
          `<div class="bg-[var(--card)] border border-[var(--border)] rounded-xl p-4 space-y-2.5">` +
            `<h3 class="text-sm font-bold">${hubRealized.length} Canonical Hub Entries Realized by Parametric Spokes</h3>` +
            `<div class="max-h-72 overflow-y-auto divide-y divide-[var(--card-border)]">` +
              hubRealized.map((e) => (
                `<button type="button" data-jump="${esc(e.label)}" class="w-full py-2 text-left flex items-center justify-between gap-2 hover:bg-[var(--background)] px-2 rounded text-xs">` +
                  `<span class="font-mono font-bold shrink-0">${esc(e.label)}</span>` +
                  `<span class="font-mono text-xs text-[var(--muted-foreground)] truncate">${esc(e.reals.map((r) => r.id + ' [' + r.mode + ']').join(', '))}</span>` +
                `</button>`
              )).join('') +
            `</div>` +
          `</div>` +
        `</div>`
      );
      document.getElementById('view-overview').querySelectorAll('[data-jump]').forEach((b) => {
        b.addEventListener('click', () => {
          state.selectedLabel = b.getAttribute('data-jump');
          document.getElementById('btn-reset-filters').click();
          setMainTab('explorer');
        });
      });
    }

    document.getElementById('hdr-badge').textContent = `v${DATA.version} \u00b7 ${DATA.counts.total_distinct_classes} Certified Classes in R\u2074`;
    document.getElementById('tab-explorer').textContent = `Atlas Explorer (${DATA.counts.total_distinct_classes})`;

    function setMainTab(tab) {
      state.mainTab = tab;
      ['explorer', 'compare', 'sandbox', 'overview'].forEach((t) => {
        const btn = document.getElementById('tab-' + t);
        const sec = document.getElementById('view-' + t);
        if (t === tab) {
          btn.className = 'px-3 py-1.5 rounded-lg text-xs font-medium transition border bg-[var(--primary)] text-[var(--primary-foreground)] border-transparent';
          sec.classList.remove('hidden');
        } else {
          btn.className = 'px-3 py-1.5 rounded-lg text-xs font-medium transition border bg-[var(--background)] text-[var(--foreground)] border-[var(--border)]';
          sec.classList.add('hidden');
        }
      });
      if (tab === 'explorer') renderCatalogList();
      if (tab === 'compare') renderCompareView();
      if (tab === 'sandbox') renderSandboxView();
      if (tab === 'overview') renderOverviewView();
    }

    document.querySelectorAll('[data-main-tab]').forEach((b) => {
      b.addEventListener('click', () => setMainTab(b.getAttribute('data-main-tab')));
    });
    ['f-search', 'f-tier', 'f-topo', 'f-ori', 'f-prop'].forEach((id) => {
      const el = document.getElementById(id);
      el.addEventListener(id === 'f-search' ? 'input' : 'change', renderCatalogList);
    });
    document.getElementById('btn-reset-filters').addEventListener('click', () => {
      document.getElementById('f-search').value = '';
      document.getElementById('f-tier').value = 'all';
      document.getElementById('f-topo').value = 'all';
      document.getElementById('f-ori').value = 'all';
      document.getElementById('f-prop').value = 'all';
      renderCatalogList();
    });
    document.getElementById('btn-render-epd').addEventListener('click', renderSandboxView);

    const optsHtml = ENTRIES.map((e) =>
      `<option value="${esc(e.label)}">${esc(e.label)} — ${esc(e.types.join(' u '))} (ch=${fmtCh(e)})</option>`
    ).join('');
    const selA = document.getElementById('cmp-a');
    const selB = document.getElementById('cmp-b');
    selA.innerHTML = optsHtml;
    selB.innerHTML = optsHtml;
    selA.value = state.cmpA;
    selB.value = state.cmpB;
    selA.addEventListener('change', () => { state.cmpA = selA.value; renderCompareView(); });
    selB.addEventListener('change', () => { state.cmpB = selB.value; renderCompareView(); });

    const PRESETS = [
      ['8_1', '9_1', 'Spun vs Ribbon 3_1'],
      ['11_1', '11_2', 'ch=11 2-Knots'],
      ['R(b^2*a)', 'R(b^2*a^-1)', 'Non-Reciprocal Z-HS4 Pair'],
      ['10_2', '10_3', '2- vs 3-Twist Spun'],
      ['10_1^0,0,1', '11_1^0,0,-2', '3-Comp Links'],
    ];
    const presetWrap = document.getElementById('compare-presets');
    presetWrap.innerHTML = PRESETS.map(([pa, pb, lbl], idx) =>
      `<button type="button" data-preset="${idx}" class="px-2 py-1 rounded-md text-xs bg-[var(--background)] border border-[var(--border)] hover:bg-[var(--secondary)]">${esc(lbl)}</button>`
    ).join('');
    presetWrap.querySelectorAll('[data-preset]').forEach((b) => {
      b.addEventListener('click', () => {
        const [rawA, rawB] = PRESETS[Number(b.getAttribute('data-preset'))];
        const pa = resolveLabel(rawA), pb = resolveLabel(rawB);
        state.cmpA = pa; state.cmpB = pb;
        selA.value = pa; selB.value = pb;
        renderCompareView();
      });
    });

    window.addEventListener('themechange', () => setMainTab(state.mainTab));
    setMainTab('explorer');
  })();
  </script>
</body>
</html>
"""


def render_interactive_atlas_html(data_dir: Path, *, compress: bool = False) -> str:
    """Build the self-contained interactive SurfaceKnots Web Atlas HTML string."""
    payload = build_atlas_payload(data_dir)
    raw_json = json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
    if compress:
        b64 = base64.b64encode(zlib.compress(raw_json.encode("utf-8"), level=9)).decode("ascii")
        tag = f'<script id="atlas-data" type="text/plain" data-encoding="deflate-b64">{b64}</script>'
    else:
        safe_json = raw_json.replace("</", "<\\/")
        tag = f'<script id="atlas-data" type="application/json">{safe_json}</script>'
    return _ATLAS_HTML_TEMPLATE.replace("__ATLAS_SCRIPT_TAG__", tag)


def export_interactive_atlas(data_dir: Path, out_html: Path) -> dict[str, int]:
    """Generate and write the self-contained interactive HTML Web Atlas to `out_html`."""
    html_str = render_interactive_atlas_html(data_dir)
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html_str)
    payload = build_atlas_payload(data_dir)
    return dict(payload["counts"])
