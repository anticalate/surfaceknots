"""Exhaustive braided tri-plane diagram census for `b <= 4`, `c_tri <= 6`.

Implements:
- `tangle_3d_key` & `enumerate_minimal_braided_tangles`: exact 3D-isotopy deduplication of braided
  `b`-bridge trivial tangles rel `dB^3` via Artin's faithful free-group representation on `F_b`.
- `enumerate_triplane_for_bridge`: exhaustive valid `(T_1, T_2, T_3)` enumeration up to `max_c_tri = 6`
  using C99 batch R1/R2 + linking-number pre-classification (`batch_check_tangle_pairs_unlink`) plus
  exact Python Reidemeister unlink certification (`TriPlaneDiagram.is_certified_valid`), including the
  exact concentrated `c_tri = 6` algebraic identity `(T_prev * g) u bar(M_i) == T_prev u bar(M_i * (-g))`.
- `run_triplane_b4_census`: cross-joins all discovered `b <= 4, c_tri <= 6` surface-links against the
  Hub (`ch <= 11`), Spoke (`ch >= 12`), Composite, and Targeted `ch = 12` catalogs, tests
  Gong-Lewis-Monkman-Osnes (2026, arXiv:2606.03799v2) Conjectures 3.6 & 3.7 on 4-bridge braided
  tri-plane diagrams, and measures the concentrated Meander Interleaving Obstruction.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
import json
from pathlib import Path
import time

from surfaceknots.ch11_census import _is_admissible, planar_connected_sum
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.fast_enum import batch_check_tangle_pairs_unlink, pack_tangles_c_buffer
from surfaceknots.invariants.fundamental_group import tietze_simplify
from surfaceknots.triplane import (
    TangleDiagram,
    TriPlaneDiagram,
    _propagate_braid_words,
    compute_triplane_invariant_record,
    concentrated_triplane_to_marked_diagram,
    generate_catalan_matchings,
    tangle_endpoint_pairing,
    triplane_fundamental_group,
    triplane_o4_key,
)
from surfaceknots.yoshikawa_data import YOSHIKAWA_BY_LABEL, compute_full_invariant_record


def tangle_3d_key(t: TangleDiagram) -> tuple[tuple[int, ...], tuple[tuple[int, ...], ...]]:
    """Canonical 3D-isotopy key `(endpoint_pairing, canonical_artin_meridian_words)` of `t` rel `dB^3`.

    By Artin's theorem and Waldhausen/Zieschang's compression-body uniqueness, two trivial `b`-bridge
    tangles in `B^3` with identical boundary arc pairing and identical boundary meridian kernel
    `ker(pi_1(dB^3 \\ dT) -> pi_1(B^3 \\ T))` are ambient isotopic in `B^3` rel `dB^3`.
    """
    b = t.num_bridges
    w_bot: list[tuple[int, ...]] = [()] * (2 * b)
    for r, (u, v) in enumerate(t.bottom_matching):
        u_lo, v_hi = (u, v) if u < v else (v, u)
        w_bot[u_lo] = (r + 1,)
        w_bot[v_hi] = (-(r + 1),)
    m_top = _propagate_braid_words(w_bot, t.braid_word)
    remap: dict[int, int] = {}
    nxt = 1
    canon_top: list[tuple[int, ...]] = []
    for w in m_top:
        cw: list[int] = []
        for lit in w:
            g = abs(lit)
            s = 1 if lit > 0 else -1
            if g not in remap:
                remap[g] = nxt * s
                nxt += 1
            cw.append(remap[g] * s)
        canon_top.append(tuple(cw))
    return (tangle_endpoint_pairing(t), tuple(canon_top))


def untouched_boundary_arches_mask(t: TangleDiagram) -> int:
    """Bitmask over `u in 0..2b-2` of adjacent bottom arches `(u, u+1)` untouched by `t.braid_word`."""
    touched = 0
    for g in t.braid_word:
        k = abs(g) - 1
        touched |= (1 << k) | (1 << (k + 1))
    mask = 0
    for u, v in t.bottom_matching:
        lo, hi = (u, v) if u < v else (v, u)
        if hi == lo + 1 and ((touched >> lo) & 3) == 0:
            mask |= 1 << lo
    return mask


def enumerate_minimal_braided_tangles(
    b: int,
    max_crossings: int = 5,
) -> tuple[list[list[TangleDiagram]], dict[tuple, tuple[int, int]]]:
    """Enumerate minimal-crossing 3D braided tangles `by_c[0..max_crossings]` at bridge number `b`."""
    by_c: list[list[TangleDiagram]] = [[] for _ in range(max_crossings + 1)]
    key_to_loc: dict[tuple, tuple[int, int]] = {}

    for m in generate_catalan_matchings(b):
        t = TangleDiagram(b, m, ())
        k = tangle_3d_key(t)
        if k not in key_to_loc:
            key_to_loc[k] = (0, len(by_c[0]))
            by_c[0].append(t)

    if b <= 1:
        return by_c, key_to_loc

    gens = [g for k in range(1, 2 * b) for g in (k, -k)]
    for c in range(1, max_crossings + 1):
        for t_prev in by_c[c - 1]:
            bw = t_prev.braid_word
            last_g = bw[-1] if bw else 0
            bm = t_prev.bottom_matching
            for g in gens:
                if last_g == -g:
                    continue
                t_new = TangleDiagram(b, bm, bw + (g,))
                k = tangle_3d_key(t_new)
                if k not in key_to_loc:
                    key_to_loc[k] = (c, len(by_c[c]))
                    by_c[c].append(t_new)

    return by_c, key_to_loc


def _serialize_tangle(t: TangleDiagram) -> dict[str, object]:
    return {
        "num_bridges": t.num_bridges,
        "bottom_matching": [list(p) for p in t.bottom_matching],
        "braid_word": list(t.braid_word),
    }


def _serialize_triplane(tp: TriPlaneDiagram) -> dict[str, object]:
    return {
        "bridge_number": tp.bridge_number,
        "c_tri": tp.crossing_number,
        "crossing_triple": list(tp.crossing_triple),
        "is_concentrated": tp.is_concentrated,
        "patch_counts": list(tp.patch_counts()),
        "t1": _serialize_tangle(tp.t1),
        "t2": _serialize_tangle(tp.t2),
        "t3": _serialize_tangle(tp.t3),
    }


@dataclass
class _FastGroupBucket:
    rep_tp: TriPlaneDiagram
    min_b: int
    min_c_tri: int
    min_conc_tp: TriPlaneDiagram | None = None
    min_unconc_tp: TriPlaneDiagram | None = None
    best_compiled_md: MarkedDiagram | None = None
    best_compiled_source_tp: TriPlaneDiagram | None = None
    crossing_triples_seen: set[tuple[int, int, int]] = field(default_factory=set)
    valid_triple_count: int = 0
    conc_compilable_count: int = 0
    conc_interleaved_count: int = 0


def _int_snf_small(rows: list[list[int]], ncols: int) -> tuple[int, tuple[int, ...]]:
    """Exact pure-Python integer Smith Normal Form `(nullity, torsion_factors)` for `ncols <= 6`."""
    a = [r[:] for r in rows if any(r)]
    m = len(a)
    k = 0
    while k < m and k < ncols:
        best_val = 0
        br = -1
        bc = -1
        for r in range(k, m):
            row_r = a[r]
            for c in range(k, ncols):
                v = abs(row_r[c])
                if v and (best_val == 0 or v < best_val):
                    best_val = v
                    br = r
                    bc = c
        if br == -1:
            break
        if br != k:
            a[k], a[br] = a[br], a[k]
        if bc != k:
            for r in range(m):
                a[r][k], a[r][bc] = a[r][bc], a[r][k]

        while True:
            piv = a[k][k]
            changed = False
            for r in range(k + 1, m):
                if a[r][k] != 0:
                    q = a[r][k] // piv
                    for c in range(k, ncols):
                        a[r][c] -= q * a[k][c]
                    if a[r][k] != 0:
                        a[k], a[r] = a[r], a[k]
                        piv = a[k][k]
                        changed = True
            for c in range(k + 1, ncols):
                if a[k][c] != 0:
                    q = a[k][c] // piv
                    for r in range(k, m):
                        a[r][c] -= q * a[r][k]
                    if a[k][c] != 0:
                        for r in range(k, m):
                            a[r][k], a[r][c] = a[r][c], a[r][k]
                        piv = a[k][k]
                        changed = True
            if changed:
                continue
            rem_found = False
            for r in range(k + 1, m):
                for c in range(k + 1, ncols):
                    if a[r][c] % piv != 0:
                        for cc in range(k, ncols):
                            a[k][cc] += a[r][cc]
                        rem_found = True
                        break
                if rem_found:
                    break
            if not rem_found:
                break
        if a[k][k] < 0:
            for c in range(k, ncols):
                a[k][c] = -a[k][c]
        k += 1

    diag = [abs(a[i][i]) for i in range(k) if a[i][i] != 0]
    r_null = ncols - len(diag)
    tors = tuple(sorted(d for d in diag if d > 1))
    return (r_null, tors)


def _meander_marker_count(
    b: int,
    m_minus_pairs: tuple[tuple[int, int], ...],
    m_plus_pairs: tuple[tuple[int, int], ...],
) -> int | None:
    """Return number of marked vertices `m` needed to compile `(m_minus, m_plus)` in `D_+` with zero
    extra classical crossings, or `None` if obstructed by meander interleaving on `S^1`.
    """
    w_count = 2 * b
    m_minus: dict[int, int] = {}
    for u, v in m_minus_pairs:
        m_minus[u] = v
        m_minus[v] = u
    m_plus: dict[int, int] = {}
    for u, v in m_plus_pairs:
        m_plus[u] = v
        m_plus[v] = u

    active = list(range(w_count))
    next_virt_id = w_count
    markers = 0

    while active:
        k_act = len(active)
        found_common = -1
        for idx in range(k_act):
            pa = active[idx]
            pb = active[(idx + 1) % k_act]
            if m_minus[pa] == pb and m_plus[pa] == pb:
                found_common = idx
                break
        if found_common != -1:
            pa = active[found_common]
            pb = active[(found_common + 1) % k_act]
            del m_minus[pa], m_minus[pb], m_plus[pa], m_plus[pb]
            if found_common + 1 < k_act:
                del active[found_common + 1]
                del active[found_common]
            else:
                active.pop()
                active.pop(0)
            continue

        if k_act < 4:
            return None
        found_triple: tuple[int, bool] | None = None
        for idx in range(k_act):
            p0 = active[idx]
            p1 = active[(idx + 1) % k_act]
            p2 = active[(idx + 2) % k_act]
            if m_minus[p0] == p1 and m_plus[p1] == p2:
                found_triple = (idx, True)
                break
            if m_plus[p0] == p1 and m_minus[p1] == p2:
                found_triple = (idx, False)
                break
        if found_triple is None:
            return None

        idx, is_y = found_triple
        p0 = active[idx]
        p1 = active[(idx + 1) % k_act]
        p2 = active[(idx + 2) % k_act]
        markers += 1
        merged_id = next_virt_id
        next_virt_id += 1
        if is_y:
            q_minus = m_minus[p2]
            q_plus = m_plus[p0]
        else:
            q_minus = m_minus[p0]
            q_plus = m_plus[p2]
        del m_minus[p0], m_minus[p1], m_minus[p2]
        del m_plus[p0], m_plus[p1], m_plus[p2]
        m_minus[merged_id] = q_minus
        m_minus[q_minus] = merged_id
        m_plus[merged_id] = q_plus
        m_plus[q_plus] = merged_id

        rem_set = {idx % k_act, (idx + 1) % k_act, (idx + 2) % k_act}
        new_active: list[int] = []
        for j_pos in range(k_act):
            if j_pos == idx:
                new_active.append(merged_id)
            elif j_pos not in rem_set:
                new_active.append(active[j_pos])
        active = new_active

    return markers


def _tangle_top_words(t: TangleDiagram) -> tuple[tuple[int, ...], ...]:
    b = t.num_bridges
    w_bot: list[tuple[int, ...]] = [()] * (2 * b)
    for r, (u, v) in enumerate(t.bottom_matching):
        u_lo, v_hi = (u, v) if u < v else (v, u)
        w_bot[u_lo] = (r + 1,)
        w_bot[v_hi] = (-(r + 1),)
    return tuple(_propagate_braid_words(w_bot, t.braid_word))


def _canon_relators_up_to_gen_perm(num_gen: int, relators: tuple[tuple[int, ...], ...]) -> tuple:
    """Canonicalize simplified relators up to `S_{num_gen}` generator permutation (`num_gen <= 4`)."""
    import itertools

    if not relators or num_gen <= 1 or num_gen > 4:
        return relators
    best: tuple | None = None
    for perm in itertools.permutations(range(1, num_gen + 1)):
        mapped = []
        for rel in relators:
            mr = tuple(perm[lit - 1] if lit > 0 else -perm[-lit - 1] for lit in rel)
            mr_inv = tuple(-x for x in reversed(mr))
            mapped.append(mr if mr <= mr_inv else mr_inv)
        cand = tuple(sorted(mapped))
        if best is None or cand < best:
            best = cand
    return best if best is not None else relators


_TIETZE_CACHE: dict[tuple, tuple[int, tuple]] = {}
_SNF_CACHE: dict[tuple, tuple[int, tuple[int, ...]]] = {}


def _tangle_cap_port_words(t: TangleDiagram) -> tuple[tuple[int, ...], ...]:
    """Precompute the `b` cyclically-reduced top-port relation words (`+-1..+-2b`) of `t`."""
    from surfaceknots.triplane import cyclically_reduce_word

    b = t.num_bridges
    w_ports: list[tuple[int, ...]] = [(p + 1,) for p in range(2 * b)]
    if t.braid_word:
        inv_bw = tuple(-g for g in reversed(t.braid_word))
        w_bot = _propagate_braid_words(w_ports, inv_bw)
    else:
        w_bot = w_ports
    out: list[tuple[int, ...]] = []
    for u, v in t.bottom_matching:
        u_lo, v_hi = (u, v) if u < v else (v, u)
        out.append(cyclically_reduce_word(w_bot[u_lo] + w_bot[v_hi]))
    return tuple(out)


def _pack_tangle_meta_array(
    tangles: list[TangleDiagram],
    taus_c: list[tuple[int, ...]],
    top_w_c: list[tuple[tuple[int, ...], ...]],
    include_caps: bool,
):
    """Pack `tangles` into a contiguous `(TangleFastMeta * max(1, n))` buffer."""
    import ctypes
    from surfaceknots.fast_enum import TangleFastMeta

    n = len(tangles)
    stride = ctypes.sizeof(TangleFastMeta)
    assert stride == 1448
    ba = bytearray(stride * max(1, n))

    for idx, t in enumerate(tangles):
        base = idx * stride
        b = t.num_bridges
        w_count = 2 * b
        tau = taus_c[idx]
        for p in range(w_count):
            ba[base + p] = tau[p]

        wire_orig = list(range(w_count))
        for g in t.braid_word:
            k = abs(g) - 1
            wire_orig[k], wire_orig[k + 1] = wire_orig[k + 1], wire_orig[k]
        orig_to_top = [0] * w_count
        for tx, bi in enumerate(wire_orig):
            orig_to_top[bi] = tx

        ie_mask = 0
        wire_port = [0] * w_count
        wire_up = [0] * w_count
        for u, v in t.bottom_matching:
            u_lo, v_hi = (u, v) if u < v else (v, u)
            pe = orig_to_top[u_lo]
            ie_mask |= 1 << pe
            wire_port[u_lo] = pe
            wire_up[u_lo] = -1
            wire_port[v_hi] = pe
            wire_up[v_hi] = 1
        ba[base + 8] = ie_mask
        nx = len(t.braid_word)
        assert nx <= 6
        ba[base + 9] = nx
        for x_i, g in enumerate(t.braid_word):
            k = abs(g) - 1
            sg = 1 if g > 0 else -1
            pL, uL = wire_port[k], wire_up[k]
            pR, uR = wire_port[k + 1], wire_up[k + 1]
            ba[base + 10 + x_i] = pL
            ba[base + 16 + x_i] = pR
            ba[base + 22 + x_i] = (sg * uL * uR) & 0xFF
            wire_port[k], wire_port[k + 1] = pR, pL
            wire_up[k], wire_up[k + 1] = uR, uL

        if include_caps:
            cwords = _tangle_cap_port_words(t)
            for r_i, cw in enumerate(cwords):
                cl = len(cw)
                assert cl <= 32
                ba[base + 28 + r_i] = cl
                off_c = base + 32 + 32 * r_i
                for q_i, lit in enumerate(cw):
                    ba[off_c + q_i] = lit & 0xFF

        w3 = top_w_c[idx]
        for p in range(w_count):
            wp = w3[p]
            wl = len(wp)
            assert wl <= 160, f"w3 word length {wl} > 160"
            ba[base + 160 + p] = wl
            off_w = base + 168 + 160 * p
            for q_i, lit in enumerate(wp):
                ba[off_w + q_i] = lit & 0xFF

    return (TangleFastMeta * max(1, n)).from_buffer(ba), ba


def _rows_to_u64_array(rows: list[int], words: int):
    import ctypes

    n_rows = len(rows)
    byte_stride = 8 * words
    ba = bytearray(max(8, n_rows * byte_stride))
    for r_i, msk in enumerate(rows):
        if msk:
            ba[r_i * byte_stride : (r_i + 1) * byte_stride] = msk.to_bytes(byte_stride, "little")
    return (ctypes.c_uint64 * max(1, n_rows * words)).from_buffer(ba), ba


def _fast_prekey(
    tp: TriPlaneDiagram,
    tau1: tuple[int, ...],
    tau2: tuple[int, ...],
    tau3: tuple[int, ...],
    w3: tuple[tuple[int, ...], ...],
) -> tuple:
    """Fast bridge-permutation-invariant `(topology, kei_snf, canon_pi1)` prekey with inline
    Tietze-II elimination and presentation memoization.
    """
    from surfaceknots.invariants.fundamental_group import GroupPresentation, reduce_word
    from surfaceknots.triplane import cyclically_reduce_word

    b = tp.bridge_number
    w_count = 2 * b

    port_scomp = [-1] * w_count
    port_sign = [0] * w_count
    orb_sizes: list[int] = []
    ori_flags: list[bool] = []
    for sp_start in range(w_count):
        if port_scomp[sp_start] != -1:
            continue
        sc = len(orb_sizes)
        cnt = 0
        is_ori = True
        port_scomp[sp_start] = sc
        port_sign[sp_start] = 1
        stk = [sp_start]
        while stk:
            p = stk.pop()
            cnt += 1
            exp_s = -port_sign[p]
            for tau in (tau1, tau2, tau3):
                q = tau[p]
                if port_scomp[q] == -1:
                    port_scomp[q] = sc
                    port_sign[q] = exp_s
                    stk.append(q)
                elif port_sign[q] != exp_s:
                    is_ori = False
        orb_sizes.append(cnt)
        ori_flags.append(is_ori)

    num_sc = len(orb_sizes)
    if all(ori_flags):
        p_sum = [0] * num_sc
        for ta, tb in ((tau1, tau2), (tau2, tau3), (tau3, tau1)):
            seen = [False] * w_count
            for p in range(w_count):
                if seen[p]:
                    continue
                p_sum[port_scomp[p]] += 1
                cur = p
                while not seen[cur]:
                    seen[cur] = True
                    nxt = ta[cur]
                    seen[nxt] = True
                    cur = tb[nxt]
        comp_sig = tuple(sorted((True, p_sum[sc] - (orb_sizes[sc] >> 1), 0) for sc in range(num_sc)))
        e_sum = 0
    else:
        topo = tp.surface_topology()
        comp_sig = tuple(
            sorted((c.orientable, c.euler_characteristic, abs(c.normal_euler)) for c in topo.components)
        )
        e_sum = abs(sum(c.normal_euler for c in topo.components))

    rels: list[tuple[int, ...]] = []
    for ti in (tp.t1, tp.t2):
        if not ti.braid_word:
            w_bot = w3
        else:
            inv_bw = tuple(-g for g in reversed(ti.braid_word))
            w_bot = _propagate_braid_words(list(w3), inv_bw)
        for u, v in ti.bottom_matching:
            u_lo, v_hi = (u, v) if u < v else (v, u)
            rel = cyclically_reduce_word(w_bot[u_lo] + w_bot[v_hi])
            if rel:
                rels.append(rel)

    kei_rows: list[tuple[int, ...]] = []
    for rel in rels:
        sgn = 1
        vec = [0] * b
        for lit in rel:
            vec[abs(lit) - 1] += sgn
            sgn = -sgn
        if any(vec):
            kei_rows.append(tuple(vec))

    snf_key = (b, tuple(sorted(kei_rows)))
    kei_snf = _SNF_CACHE.get(snf_key)
    if kei_snf is None:
        kei_snf = _int_snf_small([list(r) for r in kei_rows], b)
        _SNF_CACHE[snf_key] = kei_snf

    # Inline eliminate len == 1 and len == 2 (distinct gens) before invoking tietze_simplify:
    cur_rels = list(rels)
    alive_gens = set(range(1, b + 1))
    while True:
        cur_rels = [cyclically_reduce_word(r) for r in cur_rels if r]
        cur_rels = [r for r in cur_rels if r]
        elim_g = 0
        sub_lit = 0
        is_del = False
        for r in cur_rels:
            if len(r) == 1:
                elim_g = abs(r[0])
                is_del = True
                break
        if not is_del:
            for r in cur_rels:
                if len(r) == 2 and abs(r[0]) != abs(r[1]):
                    elim_g = abs(r[0])
                    sub_lit = -r[1] if r[0] > 0 else r[1]
                    break
        if elim_g == 0:
            break
        alive_gens.discard(elim_g)
        nxt_rels: list[tuple[int, ...]] = []
        for r in cur_rels:
            nr: list[int] = []
            for lit in r:
                if abs(lit) == elim_g:
                    if is_del:
                        continue
                    nr.append(sub_lit if lit > 0 else -sub_lit)
                else:
                    nr.append(lit)
            rw = reduce_word(tuple(nr))
            if rw:
                nxt_rels.append(rw)
        cur_rels = nxt_rels

    max_rel_len = max((len(r) for r in cur_rels), default=0)
    if max_rel_len <= 2:
        r_z2 = len({abs(r[0]) for r in cur_rels if len(r) == 2})
        return ("FREE_Z2", comp_sig, e_sum, kei_snf, len(alive_gens), r_z2)

    sorted_gens = sorted(alive_gens)
    gen_map = {g: idx + 1 for idx, g in enumerate(sorted_gens)}
    rem_b = len(sorted_gens)
    mapped_rels: list[tuple[int, ...]] = []
    for r in cur_rels:
        mr = tuple(gen_map[lit] if lit > 0 else -gen_map[-lit] for lit in r)
        mr_inv = tuple(-x for x in reversed(mr))
        mapped_rels.append(mr if mr <= mr_inv else mr_inv)
    t_key = (rem_b, tuple(sorted(set(mapped_rels))))
    cached = _TIETZE_CACHE.get(t_key)
    if cached is None:
        raw_pres = GroupPresentation(
            num_generators=rem_b,
            relators=t_key[1],
            component_meridians=((1,),) if rem_b >= 1 else ((),),
        )
        simp_pres = tietze_simplify(raw_pres)
        canon_rels = _canon_relators_up_to_gen_perm(simp_pres.num_generators, simp_pres.relators)
        cached = (simp_pres.num_generators, canon_rels)
        _TIETZE_CACHE[t_key] = cached

    return (
        comp_sig,
        e_sum,
        kei_snf,
        cached[0],
        cached[1],
    )


def _ingest_valid_triple(
    tp: TriPlaneDiagram,
    tau1: tuple[int, ...],
    tau2: tuple[int, ...],
    tau3: tuple[int, ...],
    w3: tuple[tuple[int, ...], ...],
    conc_markers: int | None,
    buckets: dict[tuple, _FastGroupBucket],
    stats: dict[str, int],
) -> None:
    stats["valid_triples_ingested"] += 1
    pk = _fast_prekey(tp, tau1, tau2, tau3, w3)
    ct = tuple(sorted(tp.crossing_triple))
    c_tri = tp.crossing_number
    b = tp.bridge_number
    is_conc = tp.is_concentrated

    bkt = buckets.get(pk)
    if bkt is None:
        assert tp.is_certified_valid(), f"Uncertified triplane triple passed C filter: {tp}"
        bkt = _FastGroupBucket(rep_tp=tp, min_b=b, min_c_tri=c_tri)
        buckets[pk] = bkt
    else:
        if (b, c_tri) < (bkt.min_b, bkt.min_c_tri):
            assert tp.is_certified_valid()
            bkt.rep_tp = tp
        if b < bkt.min_b:
            bkt.min_b = b
        if c_tri < bkt.min_c_tri:
            bkt.min_c_tri = c_tri

    bkt.valid_triple_count += 1
    bkt.crossing_triples_seen.add(ct)  # type: ignore[arg-type]

    if is_conc:
        if bkt.min_conc_tp is None or c_tri < bkt.min_conc_tp.crossing_number:
            assert tp.is_certified_valid()
            bkt.min_conc_tp = tp
        if conc_markers is not None:
            bkt.conc_compilable_count += 1
            stats["concentrated_compilable"] += 1
            cand_ch = c_tri + conc_markers
            if bkt.best_compiled_md is None or cand_ch < bkt.best_compiled_md.ch_index:
                md = concentrated_triplane_to_marked_diagram(tp)
                assert md is not None
                bkt.best_compiled_md = md
                bkt.best_compiled_source_tp = tp
        else:
            bkt.conc_interleaved_count += 1
            stats["concentrated_meander_obstructed"] += 1
    else:
        if bkt.min_unconc_tp is None or c_tri < bkt.min_unconc_tp.crossing_number:
            assert tp.is_certified_valid()
            bkt.min_unconc_tp = tp


def _ingest_c_resolved_bucket(
    bo,
    b: int,
    ct: tuple[int, int, int],
    c_tri: int,
    is_conc: bool,
    list1: list[TangleDiagram],
    list2: list[TangleDiagram],
    list3: list[TangleDiagram],
    buckets: dict[tuple, _FastGroupBucket],
    stats: dict[str, int],
) -> None:
    pk = ("C128", int(bo.key_hi), int(bo.key_lo))
    vcnt = int(bo.valid_count)
    bkt = buckets.get(pk)
    if bkt is None:
        rep_tp = TriPlaneDiagram(list1[bo.rep_i], list2[bo.rep_j], list3[bo.rep_k])
        assert rep_tp.is_certified_valid(), f"Uncertified C bucket rep: {rep_tp}"
        bkt = _FastGroupBucket(rep_tp=rep_tp, min_b=b, min_c_tri=c_tri)
        buckets[pk] = bkt
    else:
        if (b, c_tri) < (bkt.min_b, bkt.min_c_tri):
            rep_tp = TriPlaneDiagram(list1[bo.rep_i], list2[bo.rep_j], list3[bo.rep_k])
            assert rep_tp.is_certified_valid()
            bkt.rep_tp = rep_tp
        if b < bkt.min_b:
            bkt.min_b = b
        if c_tri < bkt.min_c_tri:
            bkt.min_c_tri = c_tri

    bkt.valid_triple_count += vcnt
    bkt.crossing_triples_seen.add(ct)

    if is_conc:
        if bkt.min_conc_tp is None or c_tri < bkt.min_conc_tp.crossing_number:
            rep_tp = TriPlaneDiagram(list1[bo.rep_i], list2[bo.rep_j], list3[bo.rep_k])
            assert rep_tp.is_certified_valid()
            bkt.min_conc_tp = tp_conc = rep_tp
        cc_cnt = int(bo.conc_compilable_count)
        ci_cnt = int(bo.conc_interleaved_count)
        bkt.conc_compilable_count += cc_cnt
        bkt.conc_interleaved_count += ci_cnt
        stats["concentrated_compilable"] += cc_cnt
        stats["concentrated_meander_obstructed"] += ci_cnt
        if bo.best_comp_markers >= 0:
            cand_ch = c_tri + int(bo.best_comp_markers)
            if bkt.best_compiled_md is None or cand_ch < bkt.best_compiled_md.ch_index:
                comp_tp = TriPlaneDiagram(
                    list1[bo.best_comp_i], list2[bo.best_comp_j], list3[bo.best_comp_k]
                )
                md = concentrated_triplane_to_marked_diagram(comp_tp)
                assert md is not None
                bkt.best_compiled_md = md
                bkt.best_compiled_source_tp = comp_tp
    else:
        if bkt.min_unconc_tp is None or c_tri < bkt.min_unconc_tp.crossing_number:
            rep_tp = TriPlaneDiagram(list1[bo.rep_i], list2[bo.rep_j], list3[bo.rep_k])
            assert rep_tp.is_certified_valid()
            bkt.min_unconc_tp = rep_tp


def _step_braid_one(
    tau: tuple[int, ...],
    m_top: tuple[tuple[int, ...], ...],
    g: int,
) -> tuple[tuple[tuple[int, ...], tuple[tuple[int, ...], ...]], tuple[tuple[int, ...], ...]]:
    k = abs(g) - 1
    ntau = list(tau)
    uk, uk1 = tau[k], tau[k + 1]
    if uk == k + 1:
        ntau[k] = k + 1
        ntau[k + 1] = k
    else:
        ntau[k] = uk1
        ntau[uk1] = k
        ntau[k + 1] = uk
        ntau[uk] = k + 1
    nm = tuple(_propagate_braid_words(list(m_top), (g,)))
    remap: dict[int, int] = {}
    nxt = 1
    canon: list[tuple[int, ...]] = []
    for w in nm:
        cw: list[int] = []
        for lit in w:
            ab = abs(lit)
            s = 1 if lit > 0 else -1
            v = remap.get(ab)
            if v is None:
                v = nxt * s
                remap[ab] = v
                nxt += 1
            cw.append(v * s)
        canon.append(tuple(cw))
    return (tuple(ntau), tuple(canon)), nm


def enumerate_triplane_up_to_b4(max_c_tri: int = 6) -> tuple[dict[tuple, _FastGroupBucket], dict[str, object]]:
    """Exhaustively enumerate braided tri-plane diagrams for `b = 1, 2, 3, 4` and `c_tri <= max_c_tri`."""
    import ctypes
    from surfaceknots.fast_enum import TriplaneBucketOut, TriplaneOverflowTriple, _get_lib

    assert max_c_tri <= 6
    lib = _get_lib()
    buckets: dict[tuple, _FastGroupBucket] = {}
    stats: dict[str, int] = {
        "valid_triples_ingested": 0,
        "c_resolved_triples": 0,
        "python_tietze_overflow_triples": 0,
        "pruned_shared_untouched_arch_b4": 0,
        "pruned_trivial_double_tangle": 0,
        "pruned_halftwist_reducible": 0,
        "concentrated_compilable": 0,
        "concentrated_meander_obstructed": 0,
    }
    per_bridge_tangle_counts: dict[str, list[int]] = {}

    max_bkt_cap = 4096
    max_ovf_cap = 1_000_000
    buckets_out = (TriplaneBucketOut * max_bkt_cap)()
    overflow_out = (TriplaneOverflowTriple * max_ovf_cap)()

    for b in (1, 2, 3, 4):
        max_single_c = min(5, max_c_tri) if b >= 2 else 0
        t0 = time.perf_counter()

        # Incremental tangle enumeration + exact `m_top` / `tau` / `key_to_loc`:
        by_c: list[list[TangleDiagram]] = [[] for _ in range(max_single_c + 1)]
        taus: list[list[tuple[int, ...]]] = [[] for _ in range(max_single_c + 1)]
        top_w: list[list[tuple[tuple[int, ...], ...]]] = [[] for _ in range(max_single_c + 1)]
        key_to_loc: dict[tuple, tuple[int, int]] = {}

        for m in generate_catalan_matchings(b):
            t = TangleDiagram(b, m, ())
            wb: list[tuple[int, ...]] = [()] * (2 * b)
            for r, (u, v) in enumerate(m):
                u_lo, v_hi = (u, v) if u < v else (v, u)
                wb[u_lo] = (r + 1,)
                wb[v_hi] = (-(r + 1),)
            tau = tangle_endpoint_pairing(t)
            remap: dict[int, int] = {}
            nxt = 1
            canon: list[tuple[int, ...]] = []
            for w in wb:
                lit = w[0]
                ab = abs(lit)
                s = 1 if lit > 0 else -1
                if ab not in remap:
                    remap[ab] = nxt * s
                    nxt += 1
                canon.append((remap[ab] * s,))
            k0 = (tau, tuple(canon))
            if k0 not in key_to_loc:
                key_to_loc[k0] = (0, len(by_c[0]))
                by_c[0].append(t)
                taus[0].append(tau)
                top_w[0].append(tuple(wb))

        gens = [g for k in range(1, 2 * b) for g in (k, -k)] if b >= 2 else []
        gen_to_bit = {g: gi for gi, g in enumerate(gens)}

        for c in range(1, max_single_c + 1):
            for idx_p, t_prev in enumerate(by_c[c - 1]):
                bw = t_prev.braid_word
                last_g = bw[-1] if bw else 0
                bm = t_prev.bottom_matching
                tau_p = taus[c - 1][idx_p]
                mt_p = top_w[c - 1][idx_p]
                for g in gens:
                    if last_g == -g:
                        continue
                    kg, nmg = _step_braid_one(tau_p, mt_p, g)
                    if kg not in key_to_loc:
                        key_to_loc[kg] = (c, len(by_c[c]))
                        by_c[c].append(TangleDiagram(b, bm, bw + (g,)))
                        taus[c].append(kg[0])
                        top_w[c].append(nmg)

        per_bridge_tangle_counts[str(b)] = [len(by_c[c]) for c in range(max_single_c + 1)]

        # 3D boundary-parallel cap mask: `tau[u] == u + 1` AND `len(top_w[u]) == 1`:
        arch_masks = [
            [
                sum(
                    (1 << u)
                    for u in range(2 * b - 1)
                    if taus[c][idx_t][u] == u + 1 and len(top_w[c][idx_t][u]) == 1
                )
                for idx_t in range(len(by_c[c]))
            ]
            for c in range(max_single_c + 1)
        ]
        pair_masks = [
            [
                sum((1 << u) for u in range(2 * b - 1) if taus[c][idx_t][u] == u + 1)
                for idx_t in range(len(by_c[c]))
            ]
            for c in range(max_single_c + 1)
        ]
        bufs = [pack_tangles_c_buffer(by_c[c]) if by_c[c] else None for c in range(max_single_c + 1)]

        n0 = len(by_c[0])
        conc_lut = [
            [
                _meander_marker_count(b, by_c[0][mi].bottom_matching, by_c[0][mj].bottom_matching)
                for mj in range(n0)
            ]
            for mi in range(n0)
        ]
        conc_lut_c = (ctypes.c_int8 * max(1, n0 * n0))(
            *(
                (conc_lut[mi][mj] if conc_lut[mi][mj] is not None else -1)
                for mi in range(n0)
                for mj in range(n0)
            )
        )

        # Full exact half-twist reduction masks across ALL `c = 0..max_single_c` (including `c = 5`!):
        dec_masks: list[list[int]] = [[0] * len(by_c[c]) for c in range(max_single_c + 1)]
        le_masks: list[list[int]] = [[0] * len(by_c[c]) for c in range(max_single_c + 1)]
        if b >= 2:
            for m_i in range(n0):
                le_m = 0
                pm_val = pair_masks[0][m_i]
                for gi, g in enumerate(gens):
                    if (pm_val >> (abs(g) - 1)) & 1:
                        le_m |= 1 << gi
                le_masks[0][m_i] = le_m
            for c in range(1, max_single_c + 1):
                for idx_t in range(len(by_c[c])):
                    dm = 0
                    lm = 0
                    tau_t = taus[c][idx_t]
                    mt_t = top_w[c][idx_t]
                    for gi, g in enumerate(gens):
                        kh, _ = _step_braid_one(tau_t, mt_t, -g)
                        loc = key_to_loc.get(kh)
                        if loc is not None:
                            if loc[0] < c:
                                dm |= 1 << gi
                            if loc[0] <= c:
                                lm |= 1 << gi
                    dec_masks[c][idx_t] = dm
                    le_masks[c][idx_t] = lm

        shifted_loc: dict[int, list[tuple[int, int]]] = {}
        shifted_t1_set: set[int] = set()
        if b >= 2 and max_c_tri >= 6 and max_single_c >= 5 and by_c[5]:
            for g in gens:
                locs_g: list[tuple[int, int]] = []
                for m_i in range(n0):
                    kh, _ = _step_braid_one(taus[0][m_i], top_w[0][m_i], -g)
                    c_loc, idx_loc = key_to_loc[kh]
                    assert c_loc in (0, 1)
                    locs_g.append((c_loc, idx_loc))
                    if c_loc == 1:
                        shifted_t1_set.add(idx_loc)
                shifted_loc[g] = locs_g

        adj_bits: dict[tuple[int, int], list[int]] = {}

        def build_adj(
            c1: int,
            c2: int,
            restrict_t0_compat: bool = False,
            row_subset: set[int] | None = None,
        ) -> list[int]:
            if (c1, c2) in adj_bits:
                return adj_bits[(c1, c2)]
            n1 = len(by_c[c1])
            n2 = len(by_c[c2])
            if n1 == 0 or n2 == 0:
                adj_bits[(c1, c2)] = [0] * n1
                return adj_bits[(c1, c2)]
            if restrict_t0_compat or row_subset is not None:
                t0_1 = adj_bits[(0, c1)]
                t0_2 = adj_bits[(0, c2)]
                if row_subset is not None:
                    ok1 = sorted(row_subset)
                else:
                    ok1 = [idx for idx in range(n1) if any((t0_1[m_i] >> idx) & 1 for m_i in range(n0))]
                if restrict_t0_compat:
                    ok2 = [idx for idx in range(n2) if any((t0_2[m_i] >> idx) & 1 for m_i in range(n0))]
                else:
                    ok2 = list(range(n2))
                sub1 = [by_c[c1][idx] for idx in ok1]
                sub2 = [by_c[c2][idx] for idx in ok2]
                raw = batch_check_tangle_pairs_unlink(sub1, sub2)
                rows = [0] * n1
                nb_sub = len(ok2)
                for r_pos, orig_i in enumerate(ok1):
                    msk = 0
                    off = r_pos * nb_sub
                    pos = raw.find(1, off, off + nb_sub)
                    while pos != -1:
                        msk |= 1 << ok2[pos - off]
                        pos = raw.find(1, pos + 1, off + nb_sub)
                    rows[orig_i] = msk
                adj_bits[(c1, c2)] = rows
                return rows

            raw = batch_check_tangle_pairs_unlink(by_c[c1], by_c[c2], bufs[c1], bufs[c2])
            rows = [0] * n1
            for i in range(n1):
                msk = 0
                off = i * n2
                pos = raw.find(1, off, off + n2)
                while pos != -1:
                    msk |= 1 << (pos - off)
                    pos = raw.find(1, pos + 1, off + n2)
                rows[i] = msk
            adj_bits[(c1, c2)] = rows
            return rows

        for c2 in range(max_single_c + 1):
            build_adj(0, c2)
        if max_single_c >= 1:
            for c2 in range(1, min(4, max_single_c) + 1):
                if 1 + c2 + c2 <= max_c_tri or 0 + 1 + c2 <= max_c_tri:
                    build_adj(1, c2)
            if max_single_c >= 5 and max_c_tri >= 6:
                t0_1_rows = {
                    idx for idx in range(len(by_c[1])) if any((adj_bits[(0, 1)][m_i] >> idx) & 1 for m_i in range(n0))
                }
                build_adj(1, 5, row_subset=(t0_1_rows | shifted_t1_set))
        if max_single_c >= 2:
            if max_c_tri >= 4:
                build_adj(2, 2)
            if max_single_c >= 3 and max_c_tri >= 5:
                build_adj(2, 3)
            if max_single_c >= 4 and max_c_tri >= 6:
                build_adj(2, 4, restrict_t0_compat=True)
        if max_single_c >= 3 and max_c_tri >= 6:
            build_adj(3, 3, restrict_t0_compat=True)

        # Pack C metadata & bitset arrays per `c in 0..max_single_c`:
        meta_c = []
        dm_c = []
        lm_c = []
        am_c = []
        pm_c = []
        db_c = []
        lb_c = []
        ab_c = []
        npb_c = []
        words_c = []
        for c in range(max_single_c + 1):
            nc = len(by_c[c])
            wc = max(1, (nc + 63) // 64)
            words_c.append(wc)
            meta_arr, _mba = _pack_tangle_meta_array(by_c[c], taus[c], top_w[c], include_caps=(c <= 3))
            meta_c.append((meta_arr, _mba))
            dm_c.append((ctypes.c_uint16 * max(1, nc))(*dec_masks[c]))
            lm_c.append((ctypes.c_uint16 * max(1, nc))(*le_masks[c]))
            am_c.append((ctypes.c_uint8 * max(1, nc))(*arch_masks[c]))
            pm_c.append((ctypes.c_uint8 * max(1, nc))(*pair_masks[c]))

            d_rows = [
                sum((1 << k) for k in range(nc) if (dec_masks[c][k] >> gi) & 1)
                for gi in range(14)
            ]
            l_rows = [
                sum((1 << k) for k in range(nc) if (le_masks[c][k] >> gi) & 1)
                for gi in range(14)
            ]
            a_rows = [
                sum((1 << k) for k in range(nc) if (arch_masks[c][k] >> u) & 1)
                for u in range(7)
            ]
            np_rows = [
                sum((1 << k) for k in range(nc) if not ((pair_masks[c][k] >> u) & 1))
                for u in range(7)
            ]
            db_c.append(_rows_to_u64_array(d_rows, wc))
            lb_c.append(_rows_to_u64_array(l_rows, wc))
            ab_c.append(_rows_to_u64_array(a_rows, wc))
            npb_c.append(_rows_to_u64_array(np_rows, wc))

        adj_u64: dict[tuple[int, int], tuple] = {}
        for (c1_k, c2_k), rows_k in adj_bits.items():
            adj_u64[(c1_k, c2_k)] = _rows_to_u64_array(rows_k, words_c[c2_k])

        out_stats_c = (ctypes.c_uint64 * 4)()
        num_bkt_c = ctypes.c_int(0)
        num_ovf_c = ctypes.c_int(0)

        for c1 in range(max_single_c + 1):
            for c2 in range(c1, max_single_c + 1):
                for c3 in range(c2, max_single_c + 1):
                    if c1 + c2 + c3 > max_c_tri:
                        break
                    n1, n2, n3 = len(by_c[c1]), len(by_c[c2]), len(by_c[c3])
                    if n1 == 0 or n2 == 0 or n3 == 0:
                        continue
                    is_conc_lvl = (c1 == 0 and c2 == 0)
                    lib.scan_triplane_level_c(
                        b,
                        c1,
                        c2,
                        c3,
                        n1,
                        n2,
                        n3,
                        meta_c[c1][0],
                        meta_c[c2][0],
                        meta_c[c3][0],
                        words_c[c2],
                        words_c[c3],
                        adj_u64[(c1, c2)][0],
                        adj_u64[(c1, c3)][0],
                        adj_u64[(c2, c3)][0],
                        dm_c[c1],
                        lm_c[c1],
                        dm_c[c2],
                        lm_c[c2],
                        am_c[c1],
                        pm_c[c1],
                        am_c[c2],
                        pm_c[c2],
                        db_c[c3][0],
                        lb_c[c3][0],
                        ab_c[c3][0],
                        npb_c[c3][0],
                        conc_lut_c if is_conc_lvl else None,
                        out_stats_c,
                        buckets_out,
                        max_bkt_cap,
                        ctypes.byref(num_bkt_c),
                        overflow_out,
                        max_ovf_cap,
                        ctypes.byref(num_ovf_c),
                    )
                    n_bkt = num_bkt_c.value
                    n_ovf = num_ovf_c.value
                    assert n_bkt < max_bkt_cap, f"Bucket cap overflow: {n_bkt}"
                    assert n_ovf <= max_ovf_cap, f"Overflow cap exceeded: {n_ovf} > {max_ovf_cap}"

                    stats["pruned_trivial_double_tangle"] += int(out_stats_c[1])
                    stats["pruned_halftwist_reducible"] += int(out_stats_c[2])
                    stats["pruned_shared_untouched_arch_b4"] += int(out_stats_c[3])

                    ct_lvl = (c1, c2, c3)
                    c_tri_lvl = c1 + c2 + c3
                    for bi in range(n_bkt):
                        bo = buckets_out[bi]
                        vcnt = int(bo.valid_count)
                        stats["valid_triples_ingested"] += vcnt
                        stats["c_resolved_triples"] += vcnt
                        _ingest_c_resolved_bucket(
                            bo, b, ct_lvl, c_tri_lvl, is_conc_lvl,
                            by_c[c1], by_c[c2], by_c[c3], buckets, stats,
                        )

                    stats["python_tietze_overflow_triples"] += n_ovf
                    tau_c1, tau_c2, tau_c3 = taus[c1], taus[c2], taus[c3]
                    w_c3 = top_w[c3]
                    for oi in range(n_ovf):
                        ov = overflow_out[oi]
                        ii, jj, kk = int(ov.i), int(ov.j), int(ov.k)
                        tp_cur = TriPlaneDiagram(by_c[c1][ii], by_c[c2][jj], by_c[c3][kk])
                        cm_val = (
                            (conc_lut[ii][kk] if c3 == 0 else conc_lut[ii][jj])
                            if is_conc_lvl
                            else None
                        )
                        _ingest_valid_triple(
                            tp_cur,
                            tau_c1[ii],
                            tau_c2[jj],
                            tau_c3[kk],
                            w_c3[kk],
                            cm_val,
                            buckets,
                            stats,
                        )

        # Concentrated `(0, 0, 6)` for `b >= 2` and `max_c_tri == 6` via exact algebraic join + C batch:
        if b >= 2 and max_c_tri >= 6 and max_single_c >= 5 and by_c[5]:
            n5 = len(by_c[5])
            t5_mask0 = [0] * n5
            for m_i, row_m in enumerate(adj_bits[(0, 5)]):
                tmp = row_m
                bit_m = 1 << m_i
                while tmp:
                    lb = tmp & -tmp
                    t5_mask0[lb.bit_length() - 1] |= bit_m
                    tmp ^= lb

            t5_mask1 = [0] * n5
            for t1_i, row_1 in enumerate(adj_bits[(1, 5)]):
                if row_1 == 0:
                    continue
                tmp = row_1
                bit_1 = 1 << t1_i
                while tmp:
                    lb = tmp & -tmp
                    t5_mask1[lb.bit_length() - 1] |= bit_1
                    tmp ^= lb

            am0 = arch_masks[0]
            pm0 = pair_masks[0]
            lm0 = le_masks[0]
            tau_c0 = taus[0]
            seen_t6_keys: set[tuple] = set()
            t6_list: list[TangleDiagram] = []
            t6_taus: list[tuple[int, ...]] = []
            t6_w: list[tuple[tuple[int, ...], ...]] = []
            c6_triples: list[tuple[int, int, int]] = []

            for idx5, t_prev in enumerate(by_c[5]):
                msk0 = t5_mask0[idx5]
                msk1 = t5_mask1[idx5]
                if msk0 == 0 and msk1 == 0:
                    continue
                last_g = t_prev.braid_word[-1]
                tau_prev = taus[5][idx5]
                mt_prev = top_w[5][idx5]
                for g in gens:
                    if g == -last_g:
                        continue
                    locs_g = shifted_loc[g]
                    compat_m: list[int] = []
                    for m_i in range(n0):
                        c_loc, idx_loc = locs_g[m_i]
                        if c_loc == 0:
                            if (msk0 >> idx_loc) & 1:
                                compat_m.append(m_i)
                        else:
                            if (msk1 >> idx_loc) & 1:
                                compat_m.append(m_i)
                    if len(compat_m) < 2:
                        if len(compat_m) == 1:
                            stats["pruned_trivial_double_tangle"] += 1
                        continue
                    k6, w6 = _step_braid_one(tau_prev, mt_prev, g)
                    if k6 in key_to_loc or k6 in seen_t6_keys:
                        continue
                    seen_t6_keys.add(k6)
                    tau6 = k6[0]
                    # Check full half-twist reduction mask on t6 against `key_to_loc` (c <= 5):
                    d6_mask = 0
                    for gi_chk, g_chk in enumerate(gens):
                        kh6, _ = _step_braid_one(tau6, w6, -g_chk)
                        if kh6 in key_to_loc:
                            d6_mask |= 1 << gi_chk
                    a6 = sum((1 << u) for u in range(2 * b - 1) if tau6[u] == u + 1 and len(w6[u]) == 1)
                    p6 = sum((1 << u) for u in range(2 * b - 1) if tau6[u] == u + 1)

                    pairs_kept: list[tuple[int, int]] = []
                    for p_a in range(len(compat_m)):
                        mi = compat_m[p_a]
                        ami = am0[mi]
                        pmi = pm0[mi]
                        lmi = lm0[mi]
                        for p_b in range(p_a + 1, len(compat_m)):
                            mj = compat_m[p_b]
                            if d6_mask & lmi & lm0[mj]:
                                stats["pruned_halftwist_reducible"] += 1
                                continue
                            if b == 4:
                                amj = am0[mj]
                                pmj = pm0[mj]
                                a_ij = ami & amj
                                if (
                                    (a_ij & a6)
                                    or (a_ij & ~p6)
                                    or ((ami & a6) & ~pmj)
                                    or ((amj & a6) & ~pmi)
                                ):
                                    stats["pruned_shared_untouched_arch_b4"] += 1
                                    continue
                            pairs_kept.append((mi, mj))

                    if pairs_kept:
                        k6_idx = len(t6_list)
                        t6_list.append(TangleDiagram(b, t_prev.bottom_matching, t_prev.braid_word + (g,)))
                        t6_taus.append(tau6)
                        t6_w.append(w6)
                        for mi, mj in pairs_kept:
                            c6_triples.append((mi, mj, k6_idx))

            if c6_triples:
                meta6_arr, _m6ba = _pack_tangle_meta_array(t6_list, t6_taus, t6_w, include_caps=False)
                n_tr6 = len(c6_triples)
                tr6_in = (TriplaneOverflowTriple * n_tr6)(
                    *(TriplaneOverflowTriple(mi, mj, k6_i) for mi, mj, k6_i in c6_triples)
                )
                ovf6_out = (TriplaneOverflowTriple * max(1, n_tr6))()
                lib.eval_triplane_triples_list_c(
                    b,
                    n0,
                    meta_c[0][0],
                    meta6_arr,
                    conc_lut_c,
                    n_tr6,
                    tr6_in,
                    buckets_out,
                    max_bkt_cap,
                    ctypes.byref(num_bkt_c),
                    ovf6_out,
                    n_tr6,
                    ctypes.byref(num_ovf_c),
                )
                n_bkt6 = num_bkt_c.value
                n_ovf6 = num_ovf_c.value
                for bi in range(n_bkt6):
                    bo = buckets_out[bi]
                    vcnt = int(bo.valid_count)
                    stats["valid_triples_ingested"] += vcnt
                    stats["c_resolved_triples"] += vcnt
                    _ingest_c_resolved_bucket(
                        bo, b, (0, 0, 6), 6, True,
                        by_c[0], by_c[0], t6_list, buckets, stats,
                    )
                stats["python_tietze_overflow_triples"] += n_ovf6
                for oi in range(n_ovf6):
                    ov = ovf6_out[oi]
                    mi, mj, k6_i = int(ov.i), int(ov.j), int(ov.k)
                    _ingest_valid_triple(
                        TriPlaneDiagram(by_c[0][mi], by_c[0][mj], t6_list[k6_i]),
                        tau_c0[mi],
                        tau_c0[mj],
                        t6_taus[k6_i],
                        t6_w[k6_i],
                        conc_lut[mi][mj],
                        buckets,
                        stats,
                    )

        stats[f"elapsed_b{b}_ms"] = int(round((time.perf_counter() - t0) * 1000))

    return buckets, {
        **stats,
        "tangle_counts_by_bridge": per_bridge_tangle_counts,
        "distinct_prekey_buckets": len(buckets),
    }


def _build_catalog_lookup(data_dir: Path) -> dict[tuple, list[str]]:
    """Build `triplane_o4_key -> list[catalog_label]` across Hub (`ch <= 11`), Spoke (`ch >= 12`),
    Targeted `ch = 12`, Unknots (`ch_0`), and Composites (`ch <= 12`).
    """
    from surfaceknots.complexity import MINIMAL_KLEIN_BOTTLE_E0_EPD
    from surfaceknots.ch11_census import disjoint_union

    lookup: dict[tuple, list[str]] = defaultdict(list)

    def add_rec(rec: dict[str, object], tag: str) -> None:
        k = triplane_o4_key(rec)
        if tag not in lookup[k]:
            lookup[k].append(tag)

    # 1. Yoshikawa `ch <= 10`
    for e in json.loads((data_dir / "yoshikawa_ch10.json").read_text()):
        add_rec(e, f"Hub:{e['label']}")

    # 2. Expanded `ch = 11`
    for e in json.loads((data_dir / "expanded_ch11.json").read_text())["entries"]:
        add_rec(e, f"Hub11:{e['label']}")

    # 3. Family Spoke catalog
    fam_path = data_dir / "families" / "family_catalog.json"
    if fam_path.exists():
        for e in json.loads(fam_path.read_text())["entries"]:
            if e["hub_match"] is None:
                add_rec(e, f"Spoke:{e['family_id']}")

    # 4. Targeted `ch = 12` S^2 census if available
    t12_path = data_dir / "complexity" / "targeted_ch12_13_census.json"
    if t12_path.exists():
        t12 = json.loads(t12_path.read_text())
        for e in t12.get("ch12_m2_s2_knots", {}).get("new_prime_classes", []):
            add_rec(e, f"Targeted12:{e['class_id']}")

    # 5. Standard small unknots, unlinks, Klein bottles, and `# RP^2` / `# T^2` stabilizations
    k2_0 = MarkedDiagram.from_epd(MINIMAL_KLEIN_BOTTLE_E0_EPD)
    rp_raw = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^-1"].canonical_epd)
    p_pos = rp_raw if rp_raw.surface_topology().components[0].normal_euler == 2 else rp_raw.mirror_3d()
    p_neg = p_pos.mirror_3d()
    t2 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^1"].canonical_epd)
    s2 = MarkedDiagram.unknot_sphere()

    def sum_ok(a: MarkedDiagram, b_d: MarkedDiagram) -> MarkedDiagram:
        for fl in (False, True):
            c = planar_connected_sum(a, 0, b_d, 0, flip_splice=fl)
            if _is_admissible(c):
                return c
        raise RuntimeError("Unadmissible connected sum")

    p2_e4 = sum_ok(p_pos, p_pos)
    p3_e6 = sum_ok(p2_e4, p_pos)
    p3_e2 = sum_ok(k2_0, p_pos)
    t2_g2 = sum_ok(t2, t2)

    named_extras: list[tuple[str, MarkedDiagram]] = [
        ("Unknot:K^2(e=0)", k2_0),
        ("Unknot:K^2(e=+-4)", p2_e4),
        ("Unknot:#3RP^2(e=+-2)", p3_e2),
        ("Unknot:#3RP^2(e=+-6)", p3_e6),
        ("Unknot:Sigma_2(g=2)", t2_g2),
        ("Unlink:0_1_u_0_1", MarkedDiagram(0, (), (), free_loops=2)),
        ("Unlink:0_1_u_0_1_u_0_1", MarkedDiagram(0, (), (), free_loops=3)),
        ("Unlink:0_1_u_0_1_u_0_1_u_0_1", MarkedDiagram(0, (), (), free_loops=4)),
        ("Split:0_1_u_2_1^-1", disjoint_union(s2, p_pos)),
        ("Split:0_1_u_2_1^1", disjoint_union(s2, t2)),
        ("Split:2_1^-1_u_2_1^-1(e=2,2)", disjoint_union(p_pos, p_pos)),
        ("Split:2_1^-1_u_2_1^-1*(e=2,-2)", disjoint_union(p_pos, p_neg)),
        ("Split:0_1_u_K^2(e=0)", disjoint_union(s2, k2_0)),
        ("Composite:8_1#2_1^-1", sum_ok(MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["8_1"].canonical_epd), p_pos)),
        ("Composite:8_1#2_1^1", sum_ok(MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["8_1"].canonical_epd), t2)),
        ("Split:0_1_u_8_1", disjoint_union(s2, MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["8_1"].canonical_epd))),
        ("Split:0_1_u_6_1^0,1", disjoint_union(s2, MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["6_1^0,1"].canonical_epd))),
    ]
    for tag, md in named_extras:
        add_rec(compute_full_invariant_record(md, label=tag, include_sublinks=False), tag)

    return dict(lookup)


def run_triplane_b4_census(data_dir: Path, max_c_tri: int = 6) -> dict[str, object]:
    """Run the exhaustive `b <= 4, c_tri <= max_c_tri` braided tri-plane census, cross-join with
    the `surfaceknots` catalog, and export `code/data/complexity/triplane_b4_census.json`.
    """
    t0 = time.perf_counter()
    buckets, enum_stats = enumerate_triplane_up_to_b4(max_c_tri=max_c_tri)
    dt_enum = time.perf_counter() - t0

    t1 = time.perf_counter()
    cat_lookup = _build_catalog_lookup(data_dir)

    # Merge `_FastGroupBucket`s by full `O(4)` invariant key `triplane_o4_key`:
    merged: dict[tuple, dict[str, object]] = {}
    for bkt in buckets.values():
        rec = compute_triplane_invariant_record(bkt.rep_tp)
        ok = triplane_o4_key(rec)
        cur = merged.get(ok)
        conc_c = bkt.min_conc_tp.crossing_number if bkt.min_conc_tp is not None else None
        unconc_c = bkt.min_unconc_tp.crossing_number if bkt.min_unconc_tp is not None else None
        comp_ch = bkt.best_compiled_md.ch_index if bkt.best_compiled_md is not None else None

        if cur is None:
            merged[ok] = {
                "o4_key_repr": str(ok[:4]),
                "catalog_matches": cat_lookup.get(ok, []),
                "min_bridge_in_census": bkt.min_b,
                "min_c_tri": bkt.min_c_tri,
                "min_concentrated_c_tri": conc_c,
                "min_unconcentrated_c_tri": unconc_c,
                "best_compiled_ch": comp_ch,
                "best_compiled_epd": bkt.best_compiled_md.to_epd() if bkt.best_compiled_md is not None else None,
                "crossing_triples_seen": set(bkt.crossing_triples_seen),
                "valid_triple_count": bkt.valid_triple_count,
                "conc_compilable_count": bkt.conc_compilable_count,
                "conc_interleaved_count": bkt.conc_interleaved_count,
                "witness_triplane": _serialize_triplane(bkt.rep_tp),
                "witness_unconcentrated": (
                    _serialize_triplane(bkt.min_unconc_tp) if bkt.min_unconc_tp is not None else None
                ),
                "witness_concentrated": (
                    _serialize_triplane(bkt.min_conc_tp) if bkt.min_conc_tp is not None else None
                ),
                "invariants": rec,
            }
        else:
            if (bkt.min_b, bkt.min_c_tri) < (int(cur["min_bridge_in_census"]), int(cur["min_c_tri"])):  # type: ignore[arg-type]
                cur["witness_triplane"] = _serialize_triplane(bkt.rep_tp)
                cur["invariants"] = rec
            cur["min_bridge_in_census"] = min(int(cur["min_bridge_in_census"]), bkt.min_b)  # type: ignore[arg-type]
            cur["min_c_tri"] = min(int(cur["min_c_tri"]), bkt.min_c_tri)  # type: ignore[arg-type]
            if conc_c is not None:
                if cur["min_concentrated_c_tri"] is None or conc_c < int(cur["min_concentrated_c_tri"]):  # type: ignore[arg-type]
                    cur["min_concentrated_c_tri"] = conc_c
                    cur["witness_concentrated"] = _serialize_triplane(bkt.min_conc_tp)  # type: ignore[arg-type]
            if unconc_c is not None:
                if cur["min_unconcentrated_c_tri"] is None or unconc_c < int(cur["min_unconcentrated_c_tri"]):  # type: ignore[arg-type]
                    cur["min_unconcentrated_c_tri"] = unconc_c
                    cur["witness_unconcentrated"] = _serialize_triplane(bkt.min_unconc_tp)  # type: ignore[arg-type]
            if comp_ch is not None:
                if cur["best_compiled_ch"] is None or comp_ch < int(cur["best_compiled_ch"]):  # type: ignore[arg-type]
                    cur["best_compiled_ch"] = comp_ch
                    cur["best_compiled_epd"] = bkt.best_compiled_md.to_epd()  # type: ignore[union-attr]
            cur["crossing_triples_seen"].update(bkt.crossing_triples_seen)  # type: ignore[union-attr]
            cur["valid_triple_count"] = int(cur["valid_triple_count"]) + bkt.valid_triple_count  # type: ignore[arg-type]
            cur["conc_compilable_count"] = int(cur["conc_compilable_count"]) + bkt.conc_compilable_count  # type: ignore[arg-type]
            cur["conc_interleaved_count"] = int(cur["conc_interleaved_count"]) + bkt.conc_interleaved_count  # type: ignore[arg-type]

    entries: list[dict[str, object]] = []
    for item in sorted(
        merged.values(),
        key=lambda r: (
            int(r["min_c_tri"]),
            int(r["min_bridge_in_census"]),
            int(r["invariants"]["num_components"]),  # type: ignore[index]
            str(r["invariants"]["orientability"]),  # type: ignore[index]
            -int(r["invariants"]["euler_characteristic"]),  # type: ignore[index]
            str(r["catalog_matches"]),
        ),
    ):
        item_out = dict(item)
        item_out["crossing_triples_seen"] = [list(ct) for ct in sorted(item["crossing_triples_seen"])]  # type: ignore[arg-type]
        entries.append(item_out)

    # Extract connected S^2 knots (`num_components == 1, orientability == 'Ori', chi == 2`):
    s2_knots = [
        e
        for e in entries
        if e["invariants"]["num_components"] == 1  # type: ignore[index]
        and e["invariants"]["orientability"] == "Ori"  # type: ignore[index]
        and e["invariants"]["euler_characteristic"] == 2  # type: ignore[index]
    ]
    s2_labels = [e["catalog_matches"] for e in s2_knots]

    dt_inv = time.perf_counter() - t1
    payload = {
        "max_bridge": 4,
        "max_c_tri": max_c_tri,
        "runtime_enum_sec": round(dt_enum, 3),
        "runtime_invariants_sec": round(dt_inv, 3),
        "enumeration_stats": enum_stats,
        "distinct_o4_classes_count": len(entries),
        "connected_s2_classes_up_to_c6": [
            {
                "catalog_matches": e["catalog_matches"],
                "min_bridge": e["min_bridge_in_census"],
                "min_c_tri": e["min_c_tri"],
                "min_concentrated_c_tri": e["min_concentrated_c_tri"],
                "min_unconcentrated_c_tri": e["min_unconcentrated_c_tri"],
                "best_compiled_ch": e["best_compiled_ch"],
                "crossing_triples_seen": e["crossing_triples_seen"],
                "e1_generators": e["invariants"]["e1_generators"],  # type: ignore[index]
                "sigma2_h1_torsion": e["invariants"]["sigma2_h1_torsion"],  # type: ignore[index]
            }
            for e in s2_knots
        ],
        "gong_lewis_monkman_osnes_2026_verification": {
            "paper": "Gong, Lewis-Monkman, Osnes (June 2026, arXiv:2606.03799v2)",
            "conjecture_3_6_spun_trefoil_c_tri_eq_6_at_b4": any(
                "Hub:8_1" in e["catalog_matches"] and e["min_c_tri"] == 6 for e in s2_knots  # type: ignore[operator]
            ),
            "corollary_3_5_two_twist_spun_trefoil_c_tri_eq_6_at_b4": any(
                "Hub:10_2" in e["catalog_matches"]
                and e["min_c_tri"] == 6
                and e["best_compiled_ch"] is None
                for e in s2_knots  # type: ignore[operator]
            ),
            "conjecture_3_7_only_0_1_8_1_10_2_for_s2_up_to_c6_at_b4": (
                len(s2_knots) == 3
                and all(
                    any(lbl in m for lbl in ("Hub:0_1", "Hub:8_1", "Hub:10_2"))
                    for m in s2_labels
                )
            ),
            "rank_inversion_9_1_absent_at_c6": not any(
                "Hub:9_1" in e["catalog_matches"] for e in s2_knots  # type: ignore[operator]
            ),
        },
        "entries": entries,
    }

    out_dir = data_dir / "complexity"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "triplane_b4_census.json"
    out_path.write_text(json.dumps(payload, indent=2) + "\n")
    return payload


if __name__ == "__main__":
    res = run_triplane_b4_census(Path("code/data"), max_c_tri=6)
    print(
        f"Done in {res['runtime_enum_sec']}s + {res['runtime_invariants_sec']}s! "
        f"Distinct O(4) classes: {res['distinct_o4_classes_count']}"
    )
    print("GLMO 2026 verification:", json.dumps(res["gong_lewis_monkman_osnes_2026_verification"], indent=2))
    for e in res["entries"]:
        inv = e["invariants"]
        print(
            f"  c_tri={e['min_c_tri']} (conc={e['min_concentrated_c_tri']}, unconc={e['min_unconcentrated_c_tri']}) "
            f"b={e['min_bridge_in_census']} comp_ch={e['best_compiled_ch']} "
            f"topo={inv['component_types']} e={inv['normal_euler_numbers']} "
            f"E1={inv['e1_generators']} matches={e['catalog_matches']}"
        )
