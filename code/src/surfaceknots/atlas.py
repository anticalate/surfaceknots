"""Interactive KnotInfo-style HTML Web Atlas generator for the SurfaceKnots repository.

Unifies all 33 canonical prime Hub classes (`ch <= 10` Yoshikawa + `ch = 11` Expanded Census)
with all 41 pairwise-distinct `ch >= 12` Spoke family classes (74 surface-links total).
Includes an exact client-side JavaScript port of the IRLS `p`-Laplacian Tutte + cubic Bezier
hyperbolic splitting layout and LaTeX/TikZ engine, keeping the entire self-contained HTML
under ~100 KB while unlocking a live interactive EPD diagram & splitting-movie sandbox.
"""

from __future__ import annotations

import base64
import hashlib
import json
from pathlib import Path
from typing import Any
import zlib


def _short_hash(obj: Any) -> str:
    raw = json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _fmt_torsion(tors: list[int]) -> str:
    return " \u2295 ".join(f"Z/{t}" for t in tors)


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
    exact_ch = bool(rec.get("exact_ch_certified", tier_code in ("yoshikawa_ch10", "expanded_ch11")))

    if tier_code in ("yoshikawa_ch10", "expanded_ch11"):
        cell_status = str(rec.get("cell_status", "EXACT_SINGLETON"))
        cell_reason = str(rec.get("cell_status_reason", "Canonical Hub prime census class."))
    else:
        cell_status = "EXACT_CH12_SPOKE" if (exact_ch and ch_lo == 12) else "SPOKE_BOUNDED"
        cell_reason = str(rec.get("complexity_certification_reason", ""))

    sublinks = rec.get("sublinks") or {}
    sym = rec.get("symmetry") or {}

    c3_raw = rec["mochizuki_c3_mod3"]
    c3_mir = rec["mochizuki_c3_mirror_mod3"]
    c3_unord = sorted([c3_raw, c3_mir]) if (c3_raw is not None and c3_mir is not None) else None

    # Exact canonical objects for each of the 10 invariant filtration layers
    raw_layers = [
        (rec["component_types"], rec["orientability"], rec["euler_characteristic"], rec["normal_euler_numbers"]),
        (rec["h1_free_rank"], rec["h1_torsion"]),
        (rec["e1_fingerprint_hash"], rec["e2_generators"], rec["single_var_alex"]),
        (
            rec["hom_counts"],
            rec["sl2f3_hom_count"],
            rec["sl2f3_non_meridional_peripheral_count"],
            rec["twisted_s3_delta0"],
        ),
        (
            rec["kei_invariant_factors"],
            rec["bikei_invariant_factors"],
            rec["alexander_biquandle_t"],
            rec["alexander_biquandle_s"],
        ),
        (
            c3_unord,
            rec["sym_quandle_q6_col"],
            rec["sym_quandle_q6_trivial_weights"],
            rec["sym_quandle_q6_nontrivial_weights"],
        ),
        (
            rec["double_linking_mod2"],
            rec["triple_linking_mod2"],
            bool(sublinks.get("is_brunnian", False)),
            sublinks.get("proper_sublinks", []),
        ),
        (
            rec["sigma2_euler"],
            rec["sigma2_b1"],
            rec["sigma2_h1_torsion"],
            rec["sigma2_b2"],
            rec["sigma2_abs_signature"],
        ),
        (rec["cover_betti_profile"], rec["certified_pi2_z_rank"]),
        (
            bool(sym.get("chiral_by_normal_euler", False)),
            bool(sym.get("oriented_chiral_by_cocycle", False)),
            bool(sym.get("non_invertible_by_cocycle", False)),
            bool(sym.get("constructive_so4_amphicheiral", False)),
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
    cov_keys = list(rec["cover_betti_profile"].keys())
    pi2_rk = rec["certified_pi2_z_rank"]
    pi2_str = (
        f"rk\u2265{pi2_rk if pi2_rk is not None else '?'} ({rec['pi2_certification_reason']}) "
        f"| covers:[{','.join(cov_keys)}]"
    )

    cup_sep = " \u222a "
    comp_joined = cup_sep.join(rec["component_types"])
    layer_displays = [
        f"{comp_joined} ({rec['orientability']}), \u03c7={rec['euler_characteristic']}, e={rec['normal_euler_numbers']}",
        h1_str,
        f"E1=({', '.join(rec['e1_generators'])}) | E2=({', '.join(rec['e2_generators'])}) | \u03941={rec['single_var_alex']}",
        hom_str,
        f"{kei_str} | {bik_str}",
        f"C3={c3_str} | Q6: col={rec['sym_quandle_q6_col']} (nontriv={rec['sym_quandle_q6_nontrivial_weights']})",
        f"dlk2={rec['double_linking_mod2']}, tlk2={rec['triple_linking_mod2']}, Brunnian={bool(sublinks.get('is_brunnian', False))}",
        sigma2_str,
        pi2_str,
        str(sym.get("summary", "")),
    ]

    nw = rec.get("nearest_separating_witness")
    nw_str = ""
    if isinstance(nw, dict) and nw.get("nearest_entry"):
        nw_str = (
            f"Nearest: {nw['nearest_entry']} (dist={nw['num_separating_layers']}/10 via "
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
            {"id": r["family_id"], "mode": r["match_mode"], "ch": r["compiled_ch"]}
            for r in hub_reals
        ],
        "lh": layer_hashes,
        "ld": layer_displays,
    }


def build_atlas_payload(data_dir: Path) -> dict[str, Any]:
    """Assemble the unified 74-entry Hub + Spoke payload for the interactive Web Atlas."""
    yosh_list = json.loads((data_dir / "yoshikawa_ch10.json").read_text())
    ch11_payload = json.loads((data_dir / "expanded_ch11.json").read_text())
    fam_payload = json.loads((data_dir / "families" / "family_catalog.json").read_text())

    hub_reals_map: dict[str, list[dict[str, Any]]] = fam_payload.get("hub_realizations", {})
    entries: list[dict[str, Any]] = []

    for r in yosh_list:
        entries.append(
            _normalize_entry(
                r,
                tier_code="yoshikawa_ch10",
                tier_label="Yoshikawa (ch \u2264 10)",
                hub_reals=hub_reals_map.get(r["label"], []),
            )
        )

    for r in ch11_payload["entries"]:
        entries.append(
            _normalize_entry(
                r,
                tier_code="expanded_ch11",
                tier_label="Expanded Census (ch = 11)",
                hub_reals=hub_reals_map.get(r["label"], []),
            )
        )

    for r in fam_payload["entries"]:
        if r.get("hub_match") is not None:
            continue
        is_ex12 = bool(r.get("exact_ch_certified")) and int(r["ch_lower"]) == 12
        entries.append(
            _normalize_entry(
                r,
                tier_code="exact_ch12" if is_ex12 else "spoke_ch_ge_12",
                tier_label="Exact ch = 12 (Spoke)" if is_ex12 else "Spoke Catalog (ch \u2265 12)",
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

    return {
        "version": "0.5.1",
        "counts": {
            "total_distinct_classes": len(entries),
            "yoshikawa_ch10": len(yosh_list),
            "expanded_ch11": len(ch11_payload["entries"]),
            "exact_ch12": sum(1 for e in entries if e["tier_code"] == "exact_ch12"),
            "spoke_total": sum(
                1 for e in entries if e["tier_code"] in ("exact_ch12", "spoke_ch_ge_12")
            ),
            "hub_realized_count": int(fam_payload["distinct_hub_entries_realized"]),
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
  <script src="https://www.gstatic.com/antigravity/web/dev/tailwindcss.min.js"></script>
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
          <span class="px-2 py-0.5 rounded-md text-xs font-mono bg-[var(--secondary)] text-[var(--secondary-foreground)] border border-[var(--border)]">
            v0.5.1 &middot; 74 Certified Classes in R&#8308;
          </span>
        </div>
        <p class="text-xs text-[var(--muted-foreground)]">
          Unified Hub-and-Spoke Explorer: Yoshikawa (ch &le; 10) &middot; Expanded Prime Census (ch = 11) &middot; Certified Exact ch = 12 &amp; Parametric Spoke Catalogs
        </p>
      </div>
      <!-- Primary Mode Tabs -->
      <nav class="flex flex-wrap items-center gap-1.5" role="tablist" aria-label="Atlas views">
        <button id="tab-explorer" type="button" data-main-tab="explorer"
          class="px-3 py-1.5 rounded-lg text-xs font-medium transition border">
          Atlas Explorer (74)
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
          Hub &amp; Spoke Summary
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
            <option value="all">All Tiers (74)</option>
            <option value="hub">Canonical Hub (ch &le; 11) (33)</option>
            <option value="yoshikawa_ch10">Yoshikawa Table (ch &le; 10) (23)</option>
            <option value="expanded_ch11">Expanded Census (ch = 11) (10)</option>
            <option value="exact_ch12">Certified Exact ch = 12 (5)</option>
            <option value="spoke">All Spoke ch &ge; 12 Classes (41)</option>
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
            <option value="exact_ch">Exact ch(F) Proven (38)</option>
            <option value="hub_realized">Hub Entries w/ Spoke Model (15)</option>
            <option value="chiral">Certified Chiral</option>
            <option value="q6_triple">Positive Triple Points (Q&#8326;)</option>
            <option value="erratum">Yoshikawa Table II Erratum</option>
          </select>
        </div>
      </div>

      <div class="grid grid-cols-1 lg:grid-cols-12 gap-4 items-start">
        <div class="lg:col-span-5 bg-[var(--card)] border border-[var(--border)] rounded-xl overflow-hidden">
          <div class="px-3 py-2 border-b border-[var(--card-border)] flex items-center justify-between text-xs">
            <span id="filter-count" class="font-medium text-[var(--foreground)]">74 entries</span>
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
      const re = /([XYZ])\s*\[\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\]/g;
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

    function layoutConnectedBlock(m, localOpp, localVtype) {
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
      let bestKey = [-1, -1, -1, 1];
      for (let fi = 0; fi < faces.length; fi++) {
        const wf = faceWalks[fi];
        const simple = new Set(wf).size === wf.length ? 1 : 0;
        const mkCnt = faces[fi].reduce((acc, hd) => acc + ((localVtype[hd >> 2] === 'Y' || localVtype[hd >> 2] === 'Z') ? 1 : 0), 0);
        const key = [simple, faces[fi].length, mkCnt, -fi];
        for (let k = 0; k < 4; k++) {
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
      return Array.from({ length: 5 * m }, (_, i) => [curX[i], curY[i]]);
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
    let state = {
      mainTab: 'explorer',
      selectedLabel: '10_1^0,1',
      diagramMode: 'zero',
      cmpA: 'R(b^2*a)',
      cmpB: 'R(b^2*a^-1)',
    };

    function matchesFilters(e) {
      const q = document.getElementById('f-search').value.trim().toLowerCase();
      const tier = document.getElementById('f-tier').value;
      const topo = document.getElementById('f-topo').value;
      const ori = document.getElementById('f-ori').value;
      const prop = document.getElementById('f-prop').value;

      if (tier === 'hub' && !(e.tier_code === 'yoshikawa_ch10' || e.tier_code === 'expanded_ch11')) return false;
      if (tier === 'spoke' && !(e.tier_code === 'exact_ch12' || e.tier_code === 'spoke_ch_ge_12')) return false;
      if (['yoshikawa_ch10', 'expanded_ch11', 'exact_ch12'].includes(tier) && e.tier_code !== tier) return false;

      if (topo === 's2_knot' && !(e.nc === 1 && e.types[0] === 'S^2')) return false;
      if (topo === 't2_knot' && !(e.nc === 1 && e.types[0] === 'T^2')) return false;
      if (topo === 'nonori_knot' && !(e.nc === 1 && e.ori === 'Non')) return false;
      if (topo === 'link_2' && e.nc !== 2) return false;
      if (topo === 'link_3' && e.nc !== 3) return false;
      if (topo === 'brunnian' && !e.brunn) return false;

      if (ori !== 'all' && e.ori !== ori) return false;

      if (prop === 'exact_ch' && !e.ex_ch) return false;
      if (prop === 'hub_realized' && e.reals.length === 0) return false;
      if (prop === 'chiral' && !e.chiral) return false;
      if (prop === 'q6_triple' && !e.q6_pos) return false;
      if (prop === 'erratum' && !e.err) return false;

      if (q) {
        const hay = [
          e.label, e.desc, e.types.join(' '), e.spoke, e.subfamily,
          e.e1.join(' '), e.d1, e.pi1, e.reals.map((r) => r.id).join(' '),
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
        return (
          `<button type="button" data-select-label="${esc(e.label)}" class="w-full text-left px-3 py-2.5 transition flex items-center justify-between gap-2 ${rowCls}">` +
            `<div class="min-w-0">` +
              `<div class="flex items-center gap-2">` +
                `<span class="font-mono text-xs font-bold text-[var(--foreground)] truncate">${esc(e.label)}</span>` +
                `<span class="text-xs text-[var(--muted-foreground)] truncate">${esc(e.types.join(' \u222a '))}</span>` +
              `</div>` +
              `<div class="text-xs text-[var(--muted-foreground)] truncate mt-0.5">${esc(e.desc || e.tier_label)}</div>` +
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
            e.reals.map((r) => (
              `<span class="px-2 py-0.5 rounded font-mono text-xs bg-[var(--card)] border border-[var(--border)]">` +
                `${esc(r.id)} &middot; <span class="text-[var(--muted-foreground)]">${esc(r.mode)} (ch=${r.ch})</span>` +
              `</span>`
            )).join('') +
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
              `<span class="px-2 py-0.5 rounded text-xs font-mono border border-[var(--border)]">${esc(e.cell)}</span>` +
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
            `<div class="font-semibold text-[var(--muted-foreground)]">Complexity &amp; Equivalence Cell</div>` +
            `<div><span class="text-[var(--muted-foreground)]">ch(F):</span> <span class="font-mono font-semibold">${e.ex_ch ? e.ch_lo + ' (EXACT)' : '[' + e.ch_lo + ', ' + e.ch_hi + ']'}</span> (diag=${e.ch_d}: c=${e.c}, m=${e.m})</div>` +
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
            `<div><span class="text-[var(--muted-foreground)]">&pi;&#8322;:</span> <code class="font-mono">${esc(e.ld[8])}</code></div>` +
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
      '1. Intrinsic Topology & Normal Euler',
      '2. Abelianized H\u2081(S\u2074 \\ F)',
      '3. Alexander Ideals (E\u2081, E\u2082, \u0394\u2081)',
      '4. Finite & SL(2,F\u2083) Quotients',
      '5. Kei, Bikei & Alexander Biquandle',
      '6. Mochizuki 3-Cocycle & Q\u2086 Triple-Points',
      '7. Mod-2 Double/Triple Linking & Sublinks',
      '8. 4D Double Branched Cover \u03a3\u2082(F)',
      '9. Lomonaco \u03c0\u2082 & Finite Cover Homology',
      '10. O(4) Symmetry & Chirality Profile',
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
      const exact12 = ENTRIES.filter((e) => e.tier_code === 'exact_ch12');
      const hubRealized = ENTRIES.filter((e) => e.reals.length > 0);
      document.getElementById('view-overview').innerHTML = (
        `<div class="grid grid-cols-2 lg:grid-cols-4 gap-3">` +
          `<div class="bg-[var(--card)] border border-[var(--border)] rounded-xl p-3.5">` +
            `<div class="text-xs text-[var(--muted-foreground)]">Unified Distinct Classes</div>` +
            `<div class="text-xl font-bold font-mono mt-1">${c.total_distinct_classes}</div>` +
            `<div class="text-xs text-[var(--muted-foreground)] mt-1">33 Hub + 41 Spoke (ch &ge; 12)</div>` +
          `</div>` +
          `<div class="bg-[var(--card)] border border-[var(--border)] rounded-xl p-3.5">` +
            `<div class="text-xs text-[var(--muted-foreground)]">Exact ch(F) Certified</div>` +
            `<div class="text-xl font-bold font-mono mt-1">${c.yoshikawa_ch10 + c.expanded_ch11 + c.exact_ch12}</div>` +
            `<div class="text-xs text-[var(--muted-foreground)] mt-1">23 (ch&le;10) + 10 (ch=11) + 5 (ch=12)</div>` +
          `</div>` +
          `<div class="bg-[var(--card)] border border-[var(--border)] rounded-xl p-3.5">` +
            `<div class="text-xs text-[var(--muted-foreground)]">Hub Classes Realized by Spokes</div>` +
            `<div class="text-xl font-bold font-mono mt-1">${c.hub_realized_count} / 33</div>` +
            `<div class="text-xs text-[var(--muted-foreground)] mt-1">12 move-certified via C mountain-pass</div>` +
          `</div>` +
          `<div class="bg-[var(--card)] border border-[var(--border)] rounded-xl p-3.5">` +
            `<div class="text-xs text-[var(--muted-foreground)]">Compiled Family Diagrams</div>` +
            `<div class="text-xl font-bold font-mono mt-1">${c.compiled_family_instances}</div>` +
            `<div class="text-xs text-[var(--muted-foreground)] mt-1">Spokes 1–3 constructive compilers</div>` +
          `</div>` +
        `</div>` +
        `<div class="grid grid-cols-1 lg:grid-cols-2 gap-4">` +
          `<div class="bg-[var(--card)] border border-[var(--border)] rounded-xl p-4 space-y-2.5">` +
            `<h3 class="text-sm font-bold">5 Certified Exact ch(F) = 12 Surface-Links</h3>` +
            `<p class="text-xs text-[var(--muted-foreground)]">Compiled into 12-vertex planar ch-diagrams and proven disjoint from the complete ch &le; 11 prime + composite + split reference catalog:</p>` +
            `<div class="divide-y divide-[var(--card-border)]">` +
              exact12.map((e) => (
                `<button type="button" data-jump="${esc(e.label)}" class="w-full py-2 text-left flex items-center justify-between hover:bg-[var(--background)] px-2 rounded text-xs">` +
                  `<div><span class="font-mono font-bold">${esc(e.label)}</span> &middot; <span class="text-[var(--muted-foreground)]">${esc(e.desc)}</span></div>` +
                  `<span class="font-mono text-xs px-1.5 py-0.5 rounded border border-[var(--border)]">${esc(e.types.join(' \u222a '))}</span>` +
                `</button>`
              )).join('') +
            `</div>` +
          `</div>` +
          `<div class="bg-[var(--card)] border border-[var(--border)] rounded-xl p-4 space-y-2.5">` +
            `<h3 class="text-sm font-bold">15 Canonical Hub Entries Realized by Parametric Spokes</h3>` +
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
        const [pa, pb] = PRESETS[Number(b.getAttribute('data-preset'))];
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
