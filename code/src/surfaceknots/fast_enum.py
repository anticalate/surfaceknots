"""ctypes wrapper around `_fast_enum.c` + certified Reidemeister verification."""

from __future__ import annotations

import ctypes
from pathlib import Path
import subprocess

from surfaceknots.diagrams import MarkedDiagram, VertexKind
from surfaceknots.quadrangulations import PlanarShadow, generate_prime_shadows
from surfaceknots.unlink import certify_classical_unlink


class _CandidateStruct(ctypes.Structure):
    _fields_ = [
        ("m_mask", ctypes.c_uint16),
        ("yz_mask", ctypes.c_uint16),
        ("sign_mask", ctypes.c_uint16),
    ]


_LIB: ctypes.CDLL | None = None


def _get_lib() -> ctypes.CDLL:
    global _LIB
    if _LIB is not None:
        return _LIB
    here = Path(__file__).resolve().parent
    c_src = here / "_fast_enum.c"
    so_path = here / "_fast_enum.so"
    if not so_path.exists() or c_src.stat().st_mtime > so_path.stat().st_mtime:
        subprocess.run(
            [
                "gcc",
                "-O3",
                "-march=native",
                "-shared",
                "-fPIC",
                "-Wall",
                "-Werror",
                str(c_src),
                "-o",
                str(so_path),
            ],
            check=True,
        )
    lib = ctypes.CDLL(str(so_path))
    lib.enumerate_shadow_candidates.argtypes = [
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(_CandidateStruct),
        ctypes.c_int,
    ]
    lib.enumerate_shadow_candidates.restype = ctypes.c_int
    lib.enumerate_census_level_c.argtypes = [
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int64),
    ]
    lib.enumerate_census_level_c.restype = ctypes.c_int
    lib.merge_and_filter_diagrams_c.argtypes = [
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
    ]
    lib.merge_and_filter_diagrams_c.restype = ctypes.c_int
    lib.tunnel_and_merge_diagrams_c.argtypes = [
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
    ]
    lib.tunnel_and_merge_diagrams_c.restype = ctypes.c_int
    _LIB = lib
    return lib


def decode_candidate(shadow: PlanarShadow, m_mask: int, yz_mask: int, sign_mask: int) -> MarkedDiagram:
    """Decode `(m_mask, yz_mask, sign_mask)` on `shadow` into a `MarkedDiagram`."""
    n = shadow.num_vertices
    vtypes: list[VertexKind] = []
    m_idx = 0
    c_idx = 0
    for v in range(n):
        if (m_mask >> v) & 1:
            vtypes.append("Z" if ((yz_mask >> m_idx) & 1) else "Y")
            m_idx += 1
        else:
            vtypes.append("X1" if ((sign_mask >> c_idx) & 1) else "X0")
            c_idx += 1
    return MarkedDiagram(num_vertices=n, opp=shadow.opp, vtype=tuple(vtypes), free_loops=0)


def enumerate_admissible_ch_diagrams(
    ch: int,
    m_min: int | None = None,
    m_max: int | None = None,
    skip_classical_r2: bool = True,
    prune_level: int | None = None,
) -> list[MarkedDiagram]:
    """Enumerate certified admissible ch-diagrams on prime shadows of size `ch` (`ch >= 2`)."""
    if ch == 0:
        return [MarkedDiagram.unknot_sphere()]
    if ch < 2:
        return []

    if m_min is None:
        m_min = 1 if ch == 2 else 2
    if m_max is None:
        m_max = 2 if ch == 2 else max(2, ch - 4)
    if m_min > m_max:
        return []
    if prune_level is None:
        prune_level = 2 if skip_classical_r2 else 0

    lib = _get_lib()
    shadows = generate_prime_shadows(ch)
    buf_cap = 65536
    buf = (_CandidateStruct * buf_cap)()

    certified: list[MarkedDiagram] = []
    for sh in shadows:
        opp_arr = (ctypes.c_int * (4 * ch))(*sh.opp)
        cnt = lib.enumerate_shadow_candidates(
            ch,
            opp_arr,
            m_min,
            m_max,
            prune_level,
            buf,
            buf_cap,
        )
        if cnt > buf_cap:
            buf_cap = cnt + 1024
            buf = (_CandidateStruct * buf_cap)()
            cnt = lib.enumerate_shadow_candidates(
                ch,
                opp_arr,
                m_min,
                m_max,
                prune_level,
                buf,
                buf_cap,
            )

        for idx in range(cnt):
            cand = buf[idx]
            d = decode_candidate(sh, int(cand.m_mask), int(cand.yz_mask), int(cand.sign_mask))
            if not certify_classical_unlink(d.resolve(1)).is_unlink:
                continue
            if not certify_classical_unlink(d.resolve(-1)).is_unlink:
                continue
            certified.append(d)
    return certified


_VT_MAP: tuple[VertexKind, ...] = ("X0", "X1", "Y", "Z")
_VT_INV: dict[VertexKind, int] = {"X0": 0, "X1": 1, "Y": 2, "Z": 3}


def _run_c_census_shard(
    args: tuple[int, list[tuple[int, ...]], int, int, int, bool]
) -> tuple[list[MarkedDiagram], tuple[int, int, int, int]]:
    ch, opp_tuples, m_min, m_max, max_bfs_states, do_certify = args
    num_sh = len(opp_tuples)
    if num_sh == 0:
        return [], (0, 0, 0, 0)
    lib = _get_lib()
    flat_opp = [x for tup in opp_tuples for x in tup]
    opp_arr = (ctypes.c_int * (num_sh * 4 * ch))(*flat_opp)

    max_kept = 32768
    out_opp = (ctypes.c_int * (max_kept * 4 * ch))()
    out_vt = (ctypes.c_int * (max_kept * ch))()
    out_stats = (ctypes.c_int64 * 4)()

    kept_cnt = lib.enumerate_census_level_c(
        ch,
        num_sh,
        opp_arr,
        m_min,
        m_max,
        max_bfs_states,
        out_opp,
        out_vt,
        max_kept,
        out_stats,
    )
    if kept_cnt > max_kept:
        max_kept = kept_cnt + 256
        out_opp = (ctypes.c_int * (max_kept * 4 * ch))()
        out_vt = (ctypes.c_int * (max_kept * ch))()
        kept_cnt = lib.enumerate_census_level_c(
            ch,
            num_sh,
            opp_arr,
            m_min,
            m_max,
            max_bfs_states,
            out_opp,
            out_vt,
            max_kept,
            out_stats,
        )

    diags: list[MarkedDiagram] = []
    for k in range(kept_cnt):
        base_o = k * 4 * ch
        base_v = k * ch
        d = MarkedDiagram(
            num_vertices=ch,
            opp=tuple(int(out_opp[base_o + i]) for i in range(4 * ch)),
            vtype=tuple(_VT_MAP[int(out_vt[base_v + i])] for i in range(ch)),
            free_loops=0,
        )
        if do_certify:
            if not certify_classical_unlink(d.resolve(1)).is_unlink:
                continue
            if not certify_classical_unlink(d.resolve(-1)).is_unlink:
                continue
        diags.append(d)

    stats = (int(out_stats[0]), int(out_stats[1]), int(out_stats[2]), int(out_stats[3]))
    return diags, stats


def merge_diagrams_c(
    ch: int,
    diagrams: list[MarkedDiagram],
    max_bfs_states: int = 4096,
) -> list[MarkedDiagram]:
    """Deduplicate and filter `ch`-vertex diagrams in C with up to `max_bfs_states` states."""
    if not diagrams:
        return []
    lib = _get_lib()
    num_in = len(diagrams)
    flat_opp = [x for d in diagrams for x in d.opp]
    flat_vt = [_VT_INV[vt] for d in diagrams for vt in d.vtype]
    in_opp = (ctypes.c_int * (num_in * 4 * ch))(*flat_opp)
    in_vt = (ctypes.c_int * (num_in * ch))(*flat_vt)
    max_kept = num_in + 16
    out_opp = (ctypes.c_int * (max_kept * 4 * ch))()
    out_vt = (ctypes.c_int * (max_kept * ch))()

    kept_cnt = lib.merge_and_filter_diagrams_c(
        ch,
        num_in,
        in_opp,
        in_vt,
        max_bfs_states,
        out_opp,
        out_vt,
        max_kept,
    )
    out: list[MarkedDiagram] = []
    for k in range(kept_cnt):
        base_o = k * 4 * ch
        base_v = k * ch
        d = MarkedDiagram(
            num_vertices=ch,
            opp=tuple(int(out_opp[base_o + i]) for i in range(4 * ch)),
            vtype=tuple(_VT_MAP[int(out_vt[base_v + i])] for i in range(ch)),
            free_loops=0,
        )
        if not certify_classical_unlink(d.resolve(1)).is_unlink:
            continue
        if not certify_classical_unlink(d.resolve(-1)).is_unlink:
            continue
        out.append(d)
    return out


def tunnel_diagrams_c(
    ch: int,
    diagrams: list[MarkedDiagram],
    depth_plus1: int = 2,
    depth_plus2: int = 1,
    max_bfs_states: int = 2048,
) -> list[MarkedDiagram]:
    """Run `+1` (`Omega_1+, Omega_6+`) and `+2` (`Omega_2+`) mountain-pass tunneling in C."""
    if not diagrams:
        return []
    if ch <= 2:
        return merge_diagrams_c(ch, diagrams, max_bfs_states=max_bfs_states)
    lib = _get_lib()
    num_in = len(diagrams)
    flat_opp = [x for d in diagrams for x in d.opp]
    flat_vt = [_VT_INV[vt] for d in diagrams for vt in d.vtype]
    in_opp = (ctypes.c_int * (num_in * 4 * ch))(*flat_opp)
    in_vt = (ctypes.c_int * (num_in * ch))(*flat_vt)
    max_kept = num_in + 16
    out_opp = (ctypes.c_int * (max_kept * 4 * ch))()
    out_vt = (ctypes.c_int * (max_kept * ch))()

    kept_cnt = lib.tunnel_and_merge_diagrams_c(
        ch,
        num_in,
        in_opp,
        in_vt,
        depth_plus1,
        depth_plus2,
        max_bfs_states,
        out_opp,
        out_vt,
        max_kept,
    )
    out: list[MarkedDiagram] = []
    for k in range(kept_cnt):
        base_o = k * 4 * ch
        base_v = k * ch
        d = MarkedDiagram(
            num_vertices=ch,
            opp=tuple(int(out_opp[base_o + i]) for i in range(4 * ch)),
            vtype=tuple(_VT_MAP[int(out_vt[base_v + i])] for i in range(ch)),
            free_loops=0,
        )
        if not certify_classical_unlink(d.resolve(1)).is_unlink:
            continue
        if not certify_classical_unlink(d.resolve(-1)).is_unlink:
            continue
        out.append(d)
    return out


def enumerate_census_level_candidates(
    ch: int,
    m_min: int | None = None,
    m_max: int | None = None,
    max_bfs_states: int = 2048,
    num_workers: int = 2,
) -> tuple[list[MarkedDiagram], dict[str, int]]:
    """Enumerate certified level-irreducible, 2-cut-free component representatives at `ch`."""
    if ch == 0:
        return [MarkedDiagram.unknot_sphere()], {"raw_cands": 1, "components": 1, "discarded": 0, "kept": 1}
    if ch < 2:
        return [], {"raw_cands": 0, "components": 0, "discarded": 0, "kept": 0}

    if m_min is None:
        m_min = 1 if ch == 2 else 2
    if m_max is None:
        m_max = 2 if ch == 2 else max(2, ch - 4)
    if m_min > m_max:
        return [], {"raw_cands": 0, "components": 0, "discarded": 0, "kept": 0}

    _get_lib()
    shadows = generate_prime_shadows(ch)
    opp_tuples = [sh.opp for sh in shadows]
    scan_bfs_states = min(max_bfs_states, 384) if ch >= 9 else max_bfs_states

    if num_workers <= 1 or len(opp_tuples) < 32:
        diags, st = _run_c_census_shard((ch, opp_tuples, m_min, m_max, scan_bfs_states, False))
        final_diags = merge_diagrams_c(ch, diags, max_bfs_states=max(4096, max_bfs_states * 2))
        return final_diags, {
            "raw_cands": st[0],
            "components": st[1],
            "discarded": st[2],
            "kept": len(final_diags),
        }

    from concurrent.futures import ThreadPoolExecutor

    shards = [opp_tuples[w::num_workers] for w in range(num_workers)]
    merged_raw = 0
    merged_comps = 0
    merged_disc = 0
    shard_diags: list[MarkedDiagram] = []

    with ThreadPoolExecutor(max_workers=num_workers) as ex:
        for d_list, st in ex.map(
            _run_c_census_shard,
            [(ch, sh_sub, m_min, m_max, scan_bfs_states, False) for sh_sub in shards],
        ):
            shard_diags.extend(d_list)
            merged_raw += st[0]
            merged_comps += st[1]
            merged_disc += st[2]

    final_diags = merge_diagrams_c(ch, shard_diags, max_bfs_states=max(4096, max_bfs_states * 2))
    return final_diags, {
        "raw_cands": merged_raw,
        "components": merged_comps,
        "discarded": merged_disc,
        "kept": len(final_diags),
    }


def cascade_tunnel_diagrams_c(
    ch: int,
    diagrams: list[MarkedDiagram],
    depth_plus1: int = 256,
    depth_plus2: int = 192,
    max_bfs_states: int = 2048,
) -> list[MarkedDiagram]:
    """Run a progressive `+1/+2` mountain-pass tunneling cascade in C."""
    if not diagrams:
        return []
    if ch <= 2:
        return merge_diagrams_c(ch, diagrams, max_bfs_states=max_bfs_states)
    cur = diagrams
    if len(cur) > 64 and (depth_plus1 > 32 or depth_plus2 > 32):
        cur = tunnel_diagrams_c(
            ch=ch,
            diagrams=cur,
            depth_plus1=1,
            depth_plus2=0,
            max_bfs_states=min(max_bfs_states, 64),
        )
        if len(cur) > 64 and depth_plus2 > 0:
            cur = tunnel_diagrams_c(
                ch=ch,
                diagrams=cur,
                depth_plus1=2,
                depth_plus2=2,
                max_bfs_states=min(max_bfs_states, 256),
            )
        if len(cur) > 48 and (depth_plus1 > 96 or depth_plus2 > 48):
            cur = tunnel_diagrams_c(
                ch=ch,
                diagrams=cur,
                depth_plus1=min(depth_plus1, 96),
                depth_plus2=min(depth_plus2, 48),
                max_bfs_states=max_bfs_states,
            )
    return tunnel_diagrams_c(
        ch=ch,
        diagrams=cur,
        depth_plus1=depth_plus1,
        depth_plus2=depth_plus2,
        max_bfs_states=max_bfs_states,
    )


def enumerate_prime_ch_census(
    ch: int,
    depth_plus1: int = 256,
    depth_plus2: int = 256,
    max_bfs_states: int = 2048,
    num_workers: int = 2,
    filter_cyclic_tubes: bool = True,
) -> tuple[list[MarkedDiagram], dict[str, int]]:
    """Enumerate prime, non-splittable `ch`-minimal candidates via level scan + C stabilization tunneling."""
    level_diags, stats = enumerate_census_level_candidates(
        ch=ch,
        max_bfs_states=max_bfs_states,
        num_workers=num_workers,
    )
    stats = dict(stats)
    stats["level_kept"] = len(level_diags)
    if ch <= 2:
        stats["tunneled_kept"] = len(level_diags)
        stats["cyclic_tube_pruned"] = 0
        return level_diags, stats
    tunneled = cascade_tunnel_diagrams_c(
        ch=ch,
        diagrams=level_diags,
        depth_plus1=depth_plus1,
        depth_plus2=depth_plus2,
        max_bfs_states=max_bfs_states,
    )
    stats["tunneled_kept"] = len(tunneled)
    if filter_cyclic_tubes and ch >= 9:
        from .moves import has_certified_cyclic_unknotted_tube

        filtered = [d for d in tunneled if not has_certified_cyclic_unknotted_tube(d)]
        stats["cyclic_tube_pruned"] = len(tunneled) - len(filtered)
        return filtered, stats
    stats["cyclic_tube_pruned"] = 0
    return tunneled, stats

