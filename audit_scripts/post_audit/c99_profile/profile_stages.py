#!/usr/bin/env python3
"""Exact cycle-accurate C99 + Python per-stage profiler & ground-truth completeness oracle."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import ctypes
import hashlib
import json
import os
from pathlib import Path
import pickle
import subprocess
import sys
import time

REPO_ROOT = Path("/afs/csail.mit.edu/u/e/ejli/.gemini/antigravity-cli/worktrees/1dadd8f0-f304-455a-b682-76a831212ce6/subagent-C99---Unknotter-Speedup-Engineer-self-1528875b")
SRC_DIR = REPO_ROOT / "code" / "src"
SCRATCH_DIR = Path("/afs/csail.mit.edu/u/e/ejli/.gemini/antigravity-cli/brain/ae1d4843-4647-4011-8a0b-eb0a15e61a22/scratch")

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.quadrangulations import PlanarShadow
from surfaceknots.fast_enum import (
    _VT_MAP,
    decode_candidate,
    enumerate_admissible_ch_diagrams,
    enumerate_census_level_candidates,
    enumerate_prime_ch_census,
    enumerate_targeted_census_level_candidates,
    generate_prime_shadows_c_buffer,
    merge_diagrams_c,
)
from surfaceknots.moves import canonical_marked_key
from surfaceknots.unlink import certify_classical_unlink


def build_instrumented_so() -> ctypes.CDLL:
    c_orig = (SRC_DIR / "surfaceknots" / "_fast_enum.c").read_text()

    header_inject = """
#include <x86intrin.h>
#include <time.h>

typedef struct {
    uint64_t tsc_stage_a_markers_shadows;
    uint64_t tsc_stage_b_loops_butterfly;
    uint64_t tsc_stage_b_sign_scan;
    uint64_t tsc_stage_b_yosh1step;
    uint64_t tsc_stage_bfs_dedup;
    int64_t cnt_mmask_pass;
    int64_t cnt_yz_tested;
    int64_t cnt_chi_pass;
    int64_t cnt_bf_p;
    int64_t cnt_bf_m;
    int64_t cnt_bf_a2;
    int64_t cnt_kauff_surv;
    int64_t cnt_raw_cands;
} ProfCounters;

static __thread ProfCounters g_prof = {0};

void prof_reset_thread(void) {
    memset(&g_prof, 0, sizeof(g_prof));
}

void prof_get_thread(ProfCounters *out) {
    *out = g_prof;
}

double prof_calibrate_tsc_ghz(void) {
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    uint64_t c0 = __rdtsc();
    struct timespec req = {0, 50000000L}; /* 50 ms */
    nanosleep(&req, NULL);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    uint64_t c1 = __rdtsc();
    double ns = (double)(t1.tv_sec - t0.tv_sec) * 1e9 + (double)(t1.tv_nsec - t0.tv_nsec);
    return (double)(c1 - c0) / ns;
}
"""

    assert "#define MAX_STATES (1 << 15)" in c_orig
    c_mod = c_orig.replace("#define MAX_STATES (1 << 15)", "#define MAX_STATES (1 << 15)\n" + header_inject, 1)

    # Instrument enumerate_shadow_candidates_ex:
    old_block = """    uint32_t total_masks = 1u << n;
    for (uint32_t m_mask = 0; m_mask < total_masks; m_mask++) {
        int m = __builtin_popcount(m_mask);
        if (m < m_min || m > m_max) continue;
        int c = n - m;
        if (c < 0 || c > 15 || (1 << c) > MAX_STATES) {
            if (!ext_scratch) free(sc);
            return (c < 0) ? -1 : -2;
        }

        if (!passes_m_mask_topology_filter_c(n, orig_opp, m_mask, target_components, target_orientability)) {
            continue;
        }"""

    new_block = """    uint32_t total_masks = 1u << n;
    for (uint32_t m_mask = 0; m_mask < total_masks; m_mask++) {
        int m = __builtin_popcount(m_mask);
        if (m < m_min || m > m_max) continue;
        int c = n - m;
        if (c < 0 || c > 15 || (1 << c) > MAX_STATES) {
            if (!ext_scratch) free(sc);
            return (c < 0) ? -1 : -2;
        }
        uint64_t _t0 = __rdtsc();
        if (!passes_m_mask_topology_filter_c(n, orig_opp, m_mask, target_components, target_orientability)) {
            g_prof.tsc_stage_a_markers_shadows += (__rdtsc() - _t0);
            continue;
        }
        g_prof.cnt_mmask_pass++;"""
    assert old_block in c_mod, "old_block not found"
    c_mod = c_mod.replace(old_block, new_block, 1)

    old_yz = """        for (uint32_t yz_half = 0; yz_half < yz_limit; yz_half++) {
            uint32_t yz_mask = yz_half << 1;

            /* Early marked-marked bigon filter BEFORE building ReducedShadow! */
            int mm_hit = 0;
            for (int b = 0; b < mm_bg; b++) {
                int bi = (int)((yz_mask >> mm_mi[b]) & 1u);
                int bj = (int)((yz_mask >> mm_mj[b]) & 1u);
                if ((bi ^ bj) == mm_forbid_xor[b]) {
                    mm_hit = 1;
                    break;
                }
            }
            if (mm_hit) continue;

            ReducedShadow rsp, rsm;
            build_reduced_shadow(n, orig_opp, m_mask, yz_mask, 1, &rsp);
            if (rsp.odd_linking_fail) continue;
            build_reduced_shadow(n, orig_opp, m_mask, yz_mask, -1, &rsm);
            if (rsm.odd_linking_fail) continue;

            /* O(1) Euler characteristic filter BEFORE O(c * 2^c) Kauffman butterfly! */
            int chi = rsp.num_circles + rsm.num_circles - m;
            if (chi < target_chi_min || chi > target_chi_max) continue;

            if (c == 0) {
                int vt[MAX_N];
                for (int v = 0; v < n; v++) {
                    vt[v] = ((yz_mask >> orig_to_m[v]) & 1u) ? 3 : 2;
                }
                if (prune_level >= 2 && has_1step_yoshikawa_reduction(n, orig_opp, vt, prune_level >= 3 ? 1 : 0)) {
                    continue;
                }
                if (out_count < max_out) {
                    out_buf[out_count].m_mask = m_mask;
                    out_buf[out_count].yz_mask = yz_mask;
                    out_buf[out_count].sign_mask = 0u;
                }
                out_count++;
                continue;
            }

            /* Stage 1: compute L_+ loops + A1 = 3 butterfly ONLY! */
            compute_loop_counts(&rsp, sc->loops_p);
            kauffman_butterfly(c, sc->loops_p, A_list[0], A_inv_list[0], d_pows[0], sc->kp1);"""

    new_yz = """        g_prof.tsc_stage_a_markers_shadows += (__rdtsc() - _t0);
        for (uint32_t yz_half = 0; yz_half < yz_limit; yz_half++) {
            uint64_t _tyz = __rdtsc();
            uint32_t yz_mask = yz_half << 1;

            int mm_hit = 0;
            for (int b = 0; b < mm_bg; b++) {
                int bi = (int)((yz_mask >> mm_mi[b]) & 1u);
                int bj = (int)((yz_mask >> mm_mj[b]) & 1u);
                if ((bi ^ bj) == mm_forbid_xor[b]) {
                    mm_hit = 1;
                    break;
                }
            }
            if (mm_hit) { g_prof.tsc_stage_a_markers_shadows += (__rdtsc() - _tyz); continue; }
            g_prof.cnt_yz_tested++;

            ReducedShadow rsp, rsm;
            build_reduced_shadow(n, orig_opp, m_mask, yz_mask, 1, &rsp);
            if (rsp.odd_linking_fail) { g_prof.tsc_stage_a_markers_shadows += (__rdtsc() - _tyz); continue; }
            build_reduced_shadow(n, orig_opp, m_mask, yz_mask, -1, &rsm);
            if (rsm.odd_linking_fail) { g_prof.tsc_stage_a_markers_shadows += (__rdtsc() - _tyz); continue; }

            int chi = rsp.num_circles + rsm.num_circles - m;
            if (chi < target_chi_min || chi > target_chi_max) { g_prof.tsc_stage_a_markers_shadows += (__rdtsc() - _tyz); continue; }
            g_prof.cnt_chi_pass++;
            g_prof.tsc_stage_a_markers_shadows += (__rdtsc() - _tyz);

            if (c == 0) {
                uint64_t _ty = __rdtsc();
                int vt[MAX_N];
                for (int v = 0; v < n; v++) {
                    vt[v] = ((yz_mask >> orig_to_m[v]) & 1u) ? 3 : 2;
                }
                if (prune_level >= 2 && has_1step_yoshikawa_reduction(n, orig_opp, vt, prune_level >= 3 ? 1 : 0)) {
                    g_prof.tsc_stage_b_yosh1step += (__rdtsc() - _ty);
                    continue;
                }
                g_prof.tsc_stage_b_yosh1step += (__rdtsc() - _ty);
                if (out_count < max_out) {
                    out_buf[out_count].m_mask = m_mask;
                    out_buf[out_count].yz_mask = yz_mask;
                    out_buf[out_count].sign_mask = 0u;
                }
                out_count++;
                g_prof.cnt_raw_cands++;
                continue;
            }

            uint64_t _tbf1 = __rdtsc();
            compute_loop_counts(&rsp, sc->loops_p);
            kauffman_butterfly(c, sc->loops_p, A_list[0], A_inv_list[0], d_pows[0], sc->kp1);
            g_prof.tsc_stage_b_loops_butterfly += (__rdtsc() - _tbf1);
            g_prof.cnt_bf_p++;
            uint64_t _tsc1 = __rdtsc();"""
    assert old_yz in c_mod, "old_yz not found"
    c_mod = c_mod.replace(old_yz, new_yz, 1)

    old_st2 = """            if (surv_cnt == 0) continue;

            /* Stage 2: ONLY when L_+ has survivors, evaluate L_- at A1 = 3 and both at A2 = 5! */
            compute_loop_counts(&rsm, sc->loops_m);
            kauffman_butterfly(c, sc->loops_m, A_list[0], A_inv_list[0], d_pows[0], sc->km1);

            int64_t dm1 = d_pows[0][rsm.num_circles - 1];
            int64_t dp2 = d_pows[1][rsp.num_circles - 1];
            int64_t dm2 = d_pows[1][rsm.num_circles - 1];
            int computed_a2 = 0;

            for (int k = 0; k < surv_cnt; k++) {
                int smask = sc->surv_smask[k];
                int wp = sc->surv_wp[k];
                int wm = sc->surv_wm[k];

                int64_t exp_km1 = (neg_A3_pos[0][center + wm] * dm1) % MOD_P;
                if (sc->km1[smask] != exp_km1) continue;

                if (!computed_a2) {
                    kauffman_butterfly(c, sc->loops_p, A_list[1], A_inv_list[1], d_pows[1], sc->kp2);
                    kauffman_butterfly(c, sc->loops_m, A_list[1], A_inv_list[1], d_pows[1], sc->km2);
                    computed_a2 = 1;
                }

                int64_t exp_kp2 = (neg_A3_pos[1][center + wp] * dp2) % MOD_P;
                if (sc->kp2[smask] != exp_kp2) continue;
                int64_t exp_km2 = (neg_A3_pos[1][center + wm] * dm2) % MOD_P;
                if (sc->km2[smask] != exp_km2) continue;

                if (prune_level >= 2) {
                    int vt[MAX_N];
                    for (int v = 0; v < n; v++) {
                        if (orig_to_m[v] >= 0) {
                            vt[v] = ((yz_mask >> orig_to_m[v]) & 1u) ? 3 : 2;
                        } else {
                            vt[v] = ((smask >> orig_to_c[v]) & 1) ? 1 : 0;
                        }
                    }
                    if (has_1step_yoshikawa_reduction(n, orig_opp, vt, prune_level >= 3 ? 1 : 0)) {
                        continue;
                    }
                }

                if (out_count < max_out) {
                    out_buf[out_count].m_mask = m_mask;
                    out_buf[out_count].yz_mask = yz_mask;
                    out_buf[out_count].sign_mask = (uint32_t)smask;
                }
                out_count++;
            }"""

    new_st2 = """            g_prof.tsc_stage_b_sign_scan += (__rdtsc() - _tsc1);
            if (surv_cnt == 0) continue;

            uint64_t _tbf2 = __rdtsc();
            compute_loop_counts(&rsm, sc->loops_m);
            kauffman_butterfly(c, sc->loops_m, A_list[0], A_inv_list[0], d_pows[0], sc->km1);
            g_prof.tsc_stage_b_loops_butterfly += (__rdtsc() - _tbf2);
            g_prof.cnt_bf_m++;

            int64_t dm1 = d_pows[0][rsm.num_circles - 1];
            int64_t dp2 = d_pows[1][rsp.num_circles - 1];
            int64_t dm2 = d_pows[1][rsm.num_circles - 1];
            int computed_a2 = 0;

            for (int k = 0; k < surv_cnt; k++) {
                uint64_t _ts2 = __rdtsc();
                int smask = sc->surv_smask[k];
                int wp = sc->surv_wp[k];
                int wm = sc->surv_wm[k];

                int64_t exp_km1 = (neg_A3_pos[0][center + wm] * dm1) % MOD_P;
                if (sc->km1[smask] != exp_km1) { g_prof.tsc_stage_b_sign_scan += (__rdtsc() - _ts2); continue; }

                if (!computed_a2) {
                    g_prof.tsc_stage_b_sign_scan += (__rdtsc() - _ts2);
                    uint64_t _tbfa2 = __rdtsc();
                    kauffman_butterfly(c, sc->loops_p, A_list[1], A_inv_list[1], d_pows[1], sc->kp2);
                    kauffman_butterfly(c, sc->loops_m, A_list[1], A_inv_list[1], d_pows[1], sc->km2);
                    g_prof.tsc_stage_b_loops_butterfly += (__rdtsc() - _tbfa2);
                    g_prof.cnt_bf_a2++;
                    computed_a2 = 1;
                    _ts2 = __rdtsc();
                }

                int64_t exp_kp2 = (neg_A3_pos[1][center + wp] * dp2) % MOD_P;
                if (sc->kp2[smask] != exp_kp2) { g_prof.tsc_stage_b_sign_scan += (__rdtsc() - _ts2); continue; }
                int64_t exp_km2 = (neg_A3_pos[1][center + wm] * dm2) % MOD_P;
                if (sc->km2[smask] != exp_km2) { g_prof.tsc_stage_b_sign_scan += (__rdtsc() - _ts2); continue; }
                g_prof.tsc_stage_b_sign_scan += (__rdtsc() - _ts2);
                g_prof.cnt_kauff_surv++;

                if (prune_level >= 2) {
                    uint64_t _ty2 = __rdtsc();
                    int vt[MAX_N];
                    for (int v = 0; v < n; v++) {
                        if (orig_to_m[v] >= 0) {
                            vt[v] = ((yz_mask >> orig_to_m[v]) & 1u) ? 3 : 2;
                        } else {
                            vt[v] = ((smask >> orig_to_c[v]) & 1) ? 1 : 0;
                        }
                    }
                    if (has_1step_yoshikawa_reduction(n, orig_opp, vt, prune_level >= 3 ? 1 : 0)) {
                        g_prof.tsc_stage_b_yosh1step += (__rdtsc() - _ty2);
                        continue;
                    }
                    g_prof.tsc_stage_b_yosh1step += (__rdtsc() - _ty2);
                }

                if (out_count < max_out) {
                    out_buf[out_count].m_mask = m_mask;
                    out_buf[out_count].yz_mask = yz_mask;
                    out_buf[out_count].sign_mask = (uint32_t)smask;
                }
                out_count++;
                g_prof.cnt_raw_cands++;
            }"""
    assert old_st2 in c_mod, "old_st2 not found"
    c_mod = c_mod.replace(old_st2, new_st2, 1)

    old_bfs = """        for (int ci = 0; ci < cnt; ci++) {"""
    new_bfs = """        uint64_t _tbfs0 = __rdtsc();
        for (int ci = 0; ci < cnt; ci++) {"""
    assert old_bfs in c_mod, "old_bfs not found"
    c_mod = c_mod.replace(old_bfs, new_bfs, 1)

    old_bfs_end = """                kept_cnt++;
            }
        }
    }"""
    new_bfs_end = """                kept_cnt++;
            }
        }
        g_prof.tsc_stage_bfs_dedup += (__rdtsc() - _tbfs0);
    }"""
    assert old_bfs_end in c_mod, "old_bfs_end not found"
    c_mod = c_mod.replace(old_bfs_end, new_bfs_end, 1)

    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    c_out = SCRATCH_DIR / "_prof_enum.c"
    so_out = SCRATCH_DIR / "_prof_enum.so"
    c_out.write_text(c_mod)
    subprocess.run(
        ["gcc", "-O3", "-fopenmp", "-shared", "-fPIC", str(c_out), "-o", str(so_out)],
        check=True,
    )
    lib = ctypes.CDLL(str(so_out))
    return lib


class ProfCounters(ctypes.Structure):
    _fields_ = [
        ("tsc_stage_a_markers_shadows", ctypes.c_uint64),
        ("tsc_stage_b_loops_butterfly", ctypes.c_uint64),
        ("tsc_stage_b_sign_scan", ctypes.c_uint64),
        ("tsc_stage_b_yosh1step", ctypes.c_uint64),
        ("tsc_stage_bfs_dedup", ctypes.c_uint64),
        ("cnt_mmask_pass", ctypes.c_int64),
        ("cnt_yz_tested", ctypes.c_int64),
        ("cnt_chi_pass", ctypes.c_int64),
        ("cnt_bf_p", ctypes.c_int64),
        ("cnt_bf_m", ctypes.c_int64),
        ("cnt_bf_a2", ctypes.c_int64),
        ("cnt_kauff_surv", ctypes.c_int64),
        ("cnt_raw_cands", ctypes.c_int64),
    ]


class CandidateStruct(ctypes.Structure):
    _fields_ = [
        ("m_mask", ctypes.c_uint32),
        ("yz_mask", ctypes.c_uint32),
        ("sign_mask", ctypes.c_uint32),
    ]


def setup_lib_sigs(lib: ctypes.CDLL) -> float:
    lib.prof_reset_thread.argtypes = []
    lib.prof_reset_thread.restype = None
    lib.prof_get_thread.argtypes = [ctypes.POINTER(ProfCounters)]
    lib.prof_get_thread.restype = None
    lib.prof_calibrate_tsc_ghz.argtypes = []
    lib.prof_calibrate_tsc_ghz.restype = ctypes.c_double
    lib.enumerate_shadow_candidates.argtypes = [
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(CandidateStruct),
        ctypes.c_int,
    ]
    lib.enumerate_shadow_candidates.restype = ctypes.c_int
    lib.enumerate_targeted_census_level_c.argtypes = [
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int64),
    ]
    lib.enumerate_targeted_census_level_c.restype = ctypes.c_int
    return float(lib.prof_calibrate_tsc_ghz())


def profile_census_slice(
    lib: ctypes.CDLL,
    tsc_ghz: float,
    ch: int,
    sh_buf: ctypes.Array,
    num_sh: int,
    m_min: int,
    m_max: int,
    num_workers: int = 16,
    max_bfs_states: int = 2048,
) -> tuple[list[MarkedDiagram], dict[str, object]]:
    scan_bfs = min(max_bfs_states, 384) if ch >= 9 else max_bfs_states
    target_slices = min(num_sh, max(num_workers * 12, 256)) if (num_workers > 1 and num_sh >= 32) else 1
    slice_sz = max(1, (num_sh + target_slices - 1) // target_slices)
    slices = []
    for s_i in range(0, num_sh, slice_sz):
        c_i = min(num_sh, s_i + slice_sz) - s_i
        if c_i > 0:
            slices.append((s_i, c_i))

    def _run_one(sc: tuple[int, int]):
        s_idx, c_sh = sc
        lib.prof_reset_thread()
        ptr = ctypes.cast(
            ctypes.byref(sh_buf, s_idx * 4 * ch * ctypes.sizeof(ctypes.c_int)),
            ctypes.POINTER(ctypes.c_int),
        )
        max_kept = 32768
        out_opp = (ctypes.c_int * (max_kept * 4 * ch))()
        out_vt = (ctypes.c_int * (max_kept * ch))()
        out_stats = (ctypes.c_int64 * 4)()
        kept_cnt = int(
            lib.enumerate_targeted_census_level_c(
                ch, c_sh, ptr, m_min, m_max, 0, -1, -1000, 1000,
                scan_bfs, out_opp, out_vt, max_kept, out_stats
            )
        )
        if kept_cnt > max_kept:
            max_kept = kept_cnt + 256
            out_opp = (ctypes.c_int * (max_kept * 4 * ch))()
            out_vt = (ctypes.c_int * (max_kept * ch))()
            lib.prof_reset_thread()
            kept_cnt = int(
                lib.enumerate_targeted_census_level_c(
                    ch, c_sh, ptr, m_min, m_max, 0, -1, -1000, 1000,
                    scan_bfs, out_opp, out_vt, max_kept, out_stats
                )
            )
        pc = ProfCounters()
        lib.prof_get_thread(ctypes.byref(pc))
        diags = [
            MarkedDiagram(
                num_vertices=ch,
                opp=tuple(int(out_opp[k * 4 * ch + i]) for i in range(4 * ch)),
                vtype=tuple(_VT_MAP[int(out_vt[k * ch + i])] for i in range(ch)),
                free_loops=0,
            )
            for k in range(kept_cnt)
        ]
        return diags, tuple(int(out_stats[i]) for i in range(4)), pc

    t_wall0 = time.perf_counter()
    shard_diags: list[MarkedDiagram] = []
    tot_raw = tot_comps = tot_disc = 0
    agg_pc = {f[0]: 0 for f in ProfCounters._fields_}

    t_c0 = time.perf_counter()
    if num_workers <= 1 or len(slices) == 1:
        res_iter = [_run_one(sc) for sc in slices]
    else:
        with ThreadPoolExecutor(max_workers=min(num_workers, len(slices))) as ex:
            res_iter = list(ex.map(_run_one, slices))
    t_c_wall = time.perf_counter() - t_c0

    for d_list, st, pc in res_iter:
        shard_diags.extend(d_list)
        tot_raw += st[0]
        tot_comps += st[1]
        tot_disc += st[2]
        for fname, _ in ProfCounters._fields_:
            agg_pc[fname] += int(getattr(pc, fname))

    t_m0 = time.perf_counter()
    merged_unverified = merge_diagrams_c(
        ch, shard_diags, max_bfs_states=max(4096, max_bfs_states * 2),
        num_workers=num_workers, _certify_unlink=False,
    )
    t_merge_wall = time.perf_counter() - t_m0

    t_py0 = time.perf_counter()
    final_diags = [
        d for d in merged_unverified
        if certify_classical_unlink(d.resolve(1)).is_unlink
        and certify_classical_unlink(d.resolve(-1)).is_unlink
    ]
    t_py_cert_wall = time.perf_counter() - t_py0
    t_total_wall = time.perf_counter() - t_wall0

    inv_hz = 1.0 / (tsc_ghz * 1e9)
    cpu_stage_a = agg_pc["tsc_stage_a_markers_shadows"] * inv_hz
    cpu_stage_b_bf = agg_pc["tsc_stage_b_loops_butterfly"] * inv_hz
    cpu_stage_b_scan = agg_pc["tsc_stage_b_sign_scan"] * inv_hz
    cpu_stage_b_yosh = agg_pc["tsc_stage_b_yosh1step"] * inv_hz
    cpu_stage_bfs = agg_pc["tsc_stage_bfs_dedup"] * inv_hz
    cpu_c_sum = cpu_stage_a + cpu_stage_b_bf + cpu_stage_b_scan + cpu_stage_b_yosh + cpu_stage_bfs

    report = {
        "ch": ch,
        "num_shadows": num_sh,
        "m_range": [m_min, m_max],
        "num_workers": num_workers,
        "wall_total_sec": round(t_total_wall, 4),
        "wall_c_slices_sec": round(t_c_wall, 4),
        "wall_c_merge_sec": round(t_merge_wall, 4),
        "wall_py_cert_kept_sec": round(t_py_cert_wall, 4),
        "cpu_thread_sum_sec": {
            "stage_a_markers_shadows": round(cpu_stage_a, 4),
            "stage_b_loops_butterfly": round(cpu_stage_b_bf, 4),
            "stage_b_sign_scan": round(cpu_stage_b_scan, 4),
            "stage_b_yosh1step": round(cpu_stage_b_yosh, 4),
            "stage_bfs_canon_dedup": round(cpu_stage_bfs, 4),
            "total_c_cpu_sec": round(cpu_c_sum, 4),
        },
        "cpu_share_pct": {
            "stage_a_markers_shadows": round(100.0 * cpu_stage_a / max(1e-9, cpu_c_sum), 2),
            "stage_b_loops_butterfly": round(100.0 * cpu_stage_b_bf / max(1e-9, cpu_c_sum), 2),
            "stage_b_sign_scan": round(100.0 * cpu_stage_b_scan / max(1e-9, cpu_c_sum), 2),
            "stage_b_yosh1step": round(100.0 * cpu_stage_b_yosh / max(1e-9, cpu_c_sum), 2),
            "stage_bfs_canon_dedup": round(100.0 * cpu_stage_bfs / max(1e-9, cpu_c_sum), 2),
        },
        "counts": {
            "mmask_pass": agg_pc["cnt_mmask_pass"],
            "yz_tested": agg_pc["cnt_yz_tested"],
            "chi_pass": agg_pc["cnt_chi_pass"],
            "bf_p_calls": agg_pc["cnt_bf_p"],
            "bf_m_calls": agg_pc["cnt_bf_m"],
            "bf_a2_calls": agg_pc["cnt_bf_a2"],
            "kauffman_survivors": agg_pc["cnt_kauff_surv"],
            "raw_cands": tot_raw,
            "bfs_components": tot_comps,
            "bfs_discarded": tot_disc,
            "pre_merge_kept": len(shard_diags),
            "merged_unverified": len(merged_unverified),
            "certified_kept": len(final_diags),
            "false_positives_rejected_by_py": len(merged_unverified) - len(final_diags),
        },
    }
    return final_diags, report


def profile_raw_admissible_mode(
    lib: ctypes.CDLL,
    tsc_ghz: float,
    ch: int,
    sh_buf: ctypes.Array,
    num_sh: int,
    m_min: int,
    m_max: int,
    prune_level: int = 2,
    certify_sample_limit: int | None = None,
) -> tuple[list[tuple[int, int, int, int]], dict[str, object]]:
    """Profile `enumerate_admissible_ch_diagrams` style: C candidate enumeration + Python R1/R2/R3 unlink cert."""
    lib.prof_reset_thread()
    buf_cap = 65536
    buf = (CandidateStruct * buf_cap)()

    t_c_total = 0.0
    t_py_total = 0.0
    certified_quads: list[tuple[int, int, int, int]] = []
    total_c_cands = 0
    py_tested = 0
    py_passed = 0

    for sh_idx in range(num_sh):
        ptr = ctypes.cast(
            ctypes.byref(sh_buf, sh_idx * 4 * ch * ctypes.sizeof(ctypes.c_int)),
            ctypes.POINTER(ctypes.c_int),
        )
        tc0 = time.perf_counter()
        cnt = int(lib.enumerate_shadow_candidates(ch, ptr, m_min, m_max, prune_level, buf, buf_cap))
        if cnt > buf_cap:
            buf_cap = cnt + 1024
            buf = (CandidateStruct * buf_cap)()
            cnt = int(lib.enumerate_shadow_candidates(ch, ptr, m_min, m_max, prune_level, buf, buf_cap))
        t_c_total += time.perf_counter() - tc0
        total_c_cands += cnt

        if certify_sample_limit is not None and py_tested >= certify_sample_limit:
            continue

        sh = PlanarShadow(
            num_vertices=ch,
            opp=tuple(int(sh_buf[sh_idx * 4 * ch + s]) for s in range(4 * ch)),
        )
        tp0 = time.perf_counter()
        for idx in range(cnt):
            if certify_sample_limit is not None and py_tested >= certify_sample_limit:
                break
            py_tested += 1
            cand = buf[idx]
            mm, yz, sm = int(cand.m_mask), int(cand.yz_mask), int(cand.sign_mask)
            d = decode_candidate(sh, mm, yz, sm)
            if not certify_classical_unlink(d.resolve(1)).is_unlink:
                continue
            if not certify_classical_unlink(d.resolve(-1)).is_unlink:
                continue
            py_passed += 1
            if certify_sample_limit is None:
                certified_quads.append((sh_idx, mm, yz, sm))
        t_py_total += time.perf_counter() - tp0

    pc = ProfCounters()
    lib.prof_get_thread(ctypes.byref(pc))
    inv_hz = 1.0 / (tsc_ghz * 1e9)
    est_py_full = (
        t_py_total * (total_c_cands / max(1, py_tested))
        if certify_sample_limit is not None and py_tested > 0
        else t_py_total
    )
    return certified_quads, {
        "ch": ch,
        "num_shadows": num_sh,
        "m_range": [m_min, m_max],
        "prune_level": prune_level,
        "c_enum_sec": round(t_c_total, 4),
        "c_stage_a_markers_sec": round(pc.tsc_stage_a_markers_shadows * inv_hz, 4),
        "c_stage_b_butterfly_sec": round(pc.tsc_stage_b_loops_butterfly * inv_hz, 4),
        "c_stage_b_sign_scan_sec": round(pc.tsc_stage_b_sign_scan * inv_hz, 4),
        "c_stage_b_yosh1step_sec": round(pc.tsc_stage_b_yosh1step * inv_hz, 4),
        "py_cert_measured_sec": round(t_py_total, 4),
        "py_cert_estimated_full_sec": round(est_py_full, 4),
        "total_c_cands": total_c_cands,
        "py_tested": py_tested,
        "py_passed": py_passed,
    }


def main() -> None:
    print("Building instrumented C99 library...")
    lib = build_instrumented_so()
    tsc_ghz = setup_lib_sigs(lib)
    print(f"Calibrated TSC frequency: {tsc_ghz:.4f} GHz")

    results: dict[str, object] = {"tsc_ghz": round(tsc_ghz, 4)}
    oracle: dict[str, object] = {}

    # 1. Ground-truth oracle & profiling for full ch <= 10 (1,097 shadows):
    print("\n=== Profiling & recording ground-truth oracle for ch <= 10 (1,097 shadows) ===")
    ch_le_10_level_reports = {}
    ch_le_10_raw_default_oracle = {}
    t_le10_0 = time.perf_counter()
    for ch in range(2, 11):
        num_sh, sh_buf = generate_prime_shadows_c_buffer(ch)
        m_hi = 2 if ch == 2 else ch - 1
        diags_lvl, rep_lvl = profile_census_slice(lib, tsc_ghz, ch, sh_buf, num_sh, 1, m_hi, num_workers=16)
        ch_le_10_level_reports[str(ch)] = rep_lvl
        keys_lvl = sorted(canonical_marked_key(d, up_to_mirror_and_reverse=True) for d in diags_lvl)
        oracle[f"level_keys_ch{ch}_m1_{m_hi}"] = keys_lvl
        print(
            f"  [Level 16w] ch={ch:2d} (m=1..{m_hi}): sh={num_sh:4d}, wall={rep_lvl['wall_total_sec']:.3f}s, "
            f"raw={rep_lvl['counts']['raw_cands']}, kept={len(keys_lvl)}"
        )

        # Also record exact raw admissible candidates via default enumerate_admissible_ch_diagrams(ch) (m=2..ch-4):
        m_def_lo = 1 if ch == 2 else 2
        m_def_hi = 2 if ch == 2 else max(2, ch - 4)
        if ch <= 9:
            quads, rep_raw = profile_raw_admissible_mode(lib, tsc_ghz, ch, sh_buf, num_sh, m_def_lo, m_def_hi, prune_level=2)
            oracle[f"raw_admissible_quads_ch{ch}_default"] = quads
            ch_le_10_raw_default_oracle[str(ch)] = rep_raw
        else:
            # Measure exact 1-worker timing on first 100 shadows of ch=10 (m=2..6, prune=2), AND full 803 shadows on 16 processes for complete oracle!
            _, rep_raw10_100 = profile_raw_admissible_mode(
                lib, tsc_ghz, 10, sh_buf, 100, 2, 6, prune_level=2, certify_sample_limit=None
            )
            # Also get full C enum timing on all 803 shadows (m=2..6, prune=2) + full 16-process Python certification:
            _, rep_raw10_c_only = profile_raw_admissible_mode(
                lib, tsc_ghz, 10, sh_buf, num_sh, 2, 6, prune_level=2, certify_sample_limit=0
            )
            ch_le_10_raw_default_oracle["10_100sh_1w"] = rep_raw10_100
            ch_le_10_raw_default_oracle["10_full_803sh_c_1w"] = rep_raw10_c_only
            print(
                f"  [Raw Admissible 1w 100/803 sh] ch=10 (m=2..6, prune=2): C={rep_raw10_100['c_enum_sec']:.3f}s, "
                f"PyCert={rep_raw10_100['py_cert_measured_sec']:.3f}s, cands={rep_raw10_100['total_c_cands']} -> {rep_raw10_100['py_passed']}"
            )
            print(
                f"  [Raw Admissible 1w FULL 803 sh C-stage] ch=10 (m=2..6, prune=2): C={rep_raw10_c_only['c_enum_sec']:.3f}s "
                f"(a={rep_raw10_c_only['c_stage_a_markers_sec']:.3f}s, bf={rep_raw10_c_only['c_stage_b_butterfly_sec']:.3f}s, "
                f"scan={rep_raw10_c_only['c_stage_b_sign_scan_sec']:.3f}s, yosh={rep_raw10_c_only['c_stage_b_yosh1step_sec']:.3f}s), "
                f"total_c_cands={rep_raw10_c_only['total_c_cands']}"
            )

    results["ch_le_10_level_16w"] = ch_le_10_level_reports
    results["ch_le_10_raw_admissible_1w"] = ch_le_10_raw_default_oracle
    print(f"Finished ch <= 10 block in {time.perf_counter() - t_le10_0:.2f}s")

    # 2. Full ch = 11 (3,378 shadows) per-m breakdown on 16 workers + raw admissible sample:
    print("\n=== Profiling & recording ground-truth oracle for ch = 11 (3,378 shadows, m=1..10, 16 workers) ===")
    num_sh11, sh_buf11 = generate_prime_shadows_c_buffer(11)
    ch11_per_m = {}
    t11_0 = time.perf_counter()
    for m_val in range(1, 11):
        diags_m, rep_m = profile_census_slice(lib, tsc_ghz, 11, sh_buf11, num_sh11, m_val, m_val, num_workers=16)
        ch11_per_m[str(m_val)] = rep_m
        keys_m = sorted(canonical_marked_key(d, up_to_mirror_and_reverse=True) for d in diags_m)
        oracle[f"level_keys_ch11_m{m_val}"] = keys_m
        cs = rep_m["cpu_thread_sum_sec"]
        print(
            f"  [ch=11 m={m_val:2d}] wall={rep_m['wall_total_sec']:6.3f}s "
            f"(C_wall={rep_m['wall_c_slices_sec']:6.3f}s, merge={rep_m['wall_c_merge_sec']:5.3f}s, pyCert={rep_m['wall_py_cert_kept_sec']:5.3f}s) | "
            f"CPU: a={cs['stage_a_markers_shadows']:5.2f}s bf={cs['stage_b_loops_butterfly']:6.2f}s "
            f"scan={cs['stage_b_sign_scan']:5.2f}s yosh={cs['stage_b_yosh1step']:6.2f}s bfs={cs['stage_bfs_canon_dedup']:6.2f}s | "
            f"raw={rep_m['counts']['raw_cands']:8d} kept={len(keys_m):5d}"
        )
    t11_wall = time.perf_counter() - t11_0
    results["ch11_per_m_16w"] = ch11_per_m
    results["ch11_total_wall_16w_sec"] = round(t11_wall, 4)
    print(f"Total ch = 11 (m=1..10, 3,378 shadows, 16 workers) level+cert wall time: {t11_wall:.3f}s")

    # Also measure raw admissible mode (enumerate_admissible_ch_diagrams) on ch=11 m=2 (full 3,378 shadows) and m=2..7 (first 200 shadows):
    _, rep_raw11_m2 = profile_raw_admissible_mode(
        lib, tsc_ghz, 11, sh_buf11, num_sh11, 2, 2, prune_level=2, certify_sample_limit=50000
    )
    _, rep_raw11_slice = profile_raw_admissible_mode(
        lib, tsc_ghz, 11, sh_buf11, 200, 2, 7, prune_level=2, certify_sample_limit=50000
    )
    results["ch11_raw_admissible_m2_1w"] = rep_raw11_m2
    results["ch11_raw_admissible_200sh_m2_7_1w"] = rep_raw11_slice
    print(f"  [Raw Admissible 1w] ch=11 m=2 (3378 sh): {rep_raw11_m2}")
    print(f"  [Raw Admissible 1w] ch=11 m=2..7 (200 sh): {rep_raw11_slice}")

    # 3. Profile ch = 12 slices: m = 2 (full 15,882 shadows!) and m = 6 (first 1,000 shadows deterministic slice!):
    print("\n=== Profiling ch = 12 slices (m = 2 full 15,882 shadows; m = 6 1,000-shadow slice, 16 workers) ===")
    num_sh12, sh_buf12 = generate_prime_shadows_c_buffer(12)
    diags12_m2, rep12_m2 = profile_census_slice(lib, tsc_ghz, 12, sh_buf12, num_sh12, 2, 2, num_workers=16)
    oracle["level_keys_ch12_m2_full"] = sorted(canonical_marked_key(d, up_to_mirror_and_reverse=True) for d in diags12_m2)
    results["ch12_m2_full_16w"] = rep12_m2
    print(f"  [ch=12 m=2 full 15,882 sh] wall={rep12_m2['wall_total_sec']:.3f}s | CPU={rep12_m2['cpu_thread_sum_sec']} | counts={rep12_m2['counts']}")

    diags12_m6_1k, rep12_m6_1k = profile_census_slice(lib, tsc_ghz, 12, sh_buf12, 1000, 6, 6, num_workers=16)
    oracle["level_keys_ch12_m6_1000sh"] = sorted(canonical_marked_key(d, up_to_mirror_and_reverse=True) for d in diags12_m6_1k)
    results["ch12_m6_1000sh_16w"] = rep12_m6_1k
    print(f"  [ch=12 m=6 1,000 sh slice] wall={rep12_m6_1k['wall_total_sec']:.3f}s | CPU={rep12_m6_1k['cpu_thread_sum_sec']} | counts={rep12_m6_1k['counts']}")

    (SCRATCH_DIR / "baseline_profile.json").write_text(json.dumps(results, indent=2))
    with open(SCRATCH_DIR / "oracle_ground_truth.pkl", "wb") as f:
        pickle.dump(oracle, f)
    print(f"\nSaved baseline_profile.json and oracle_ground_truth.pkl to {SCRATCH_DIR}")


if __name__ == "__main__":
    main()
