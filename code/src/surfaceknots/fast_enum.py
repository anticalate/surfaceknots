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
        ("m_mask", ctypes.c_uint32),
        ("yz_mask", ctypes.c_uint32),
        ("sign_mask", ctypes.c_uint32),
    ]


class TangleFastMeta(ctypes.Structure):
    _fields_ = [
        ("tau", ctypes.c_uint8 * 8),
        ("is_entry_mask", ctypes.c_uint8),
        ("nxings", ctypes.c_uint8),
        ("xing_pL", ctypes.c_uint8 * 6),
        ("xing_pR", ctypes.c_uint8 * 6),
        ("xing_bs", ctypes.c_int8 * 6),
        ("cap_len", ctypes.c_uint8 * 4),
        ("cap_words", (ctypes.c_int8 * 32) * 4),
        ("w3_len", ctypes.c_uint8 * 8),
        ("w3_words", (ctypes.c_int8 * 160) * 8),
    ]


class TriplaneBucketOut(ctypes.Structure):
    _fields_ = [
        ("key_hi", ctypes.c_uint64),
        ("key_lo", ctypes.c_uint64),
        ("valid_count", ctypes.c_uint64),
        ("conc_compilable_count", ctypes.c_uint64),
        ("conc_interleaved_count", ctypes.c_uint64),
        ("rep_i", ctypes.c_int32),
        ("rep_j", ctypes.c_int32),
        ("rep_k", ctypes.c_int32),
        ("best_comp_markers", ctypes.c_int32),
        ("best_comp_i", ctypes.c_int32),
        ("best_comp_j", ctypes.c_int32),
        ("best_comp_k", ctypes.c_int32),
        ("pad0", ctypes.c_int32),
    ]


class TriplaneOverflowTriple(ctypes.Structure):
    _fields_ = [
        ("i", ctypes.c_int32),
        ("j", ctypes.c_int32),
        ("k", ctypes.c_int32),
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
    lib.generate_prime_shadows_c.argtypes = [
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_int,
    ]
    lib.generate_prime_shadows_c.restype = ctypes.c_int
    lib.batch_check_tangle_pairs_c.argtypes = [
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int8),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int8),
        ctypes.POINTER(ctypes.c_uint8),
    ]
    lib.batch_check_tangle_pairs_c.restype = None
    lib.scan_triplane_level_c.argtypes = [
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(TangleFastMeta),
        ctypes.POINTER(TangleFastMeta),
        ctypes.POINTER(TangleFastMeta),
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_uint64),
        ctypes.POINTER(ctypes.c_uint64),
        ctypes.POINTER(ctypes.c_uint64),
        ctypes.POINTER(ctypes.c_uint16),
        ctypes.POINTER(ctypes.c_uint16),
        ctypes.POINTER(ctypes.c_uint16),
        ctypes.POINTER(ctypes.c_uint16),
        ctypes.POINTER(ctypes.c_uint8),
        ctypes.POINTER(ctypes.c_uint8),
        ctypes.POINTER(ctypes.c_uint8),
        ctypes.POINTER(ctypes.c_uint8),
        ctypes.POINTER(ctypes.c_uint64),
        ctypes.POINTER(ctypes.c_uint64),
        ctypes.POINTER(ctypes.c_uint64),
        ctypes.POINTER(ctypes.c_uint64),
        ctypes.POINTER(ctypes.c_int8),
        ctypes.POINTER(ctypes.c_uint64),
        ctypes.POINTER(TriplaneBucketOut),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(TriplaneOverflowTriple),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.scan_triplane_level_c.restype = None
    lib.eval_triplane_triples_list_c.argtypes = [
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(TangleFastMeta),
        ctypes.POINTER(TangleFastMeta),
        ctypes.POINTER(ctypes.c_int8),
        ctypes.c_int,
        ctypes.POINTER(TriplaneOverflowTriple),
        ctypes.POINTER(TriplaneBucketOut),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(TriplaneOverflowTriple),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.eval_triplane_triples_list_c.restype = None
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


# Exact OEIS A113201 upper capacities up to n = 16 for 1-pass allocation:
_A113201_HINTS: dict[int, int] = {
    2: 1,
    3: 1,
    4: 2,
    5: 3,
    6: 9,
    7: 18,
    8: 62,
    9: 198,
    10: 803,
    11: 3378,
    12: 15882,
    13: 77185,
}


def generate_prime_shadows_c_buffer(num_vertices: int) -> tuple[int, ctypes.Array[ctypes.c_int]]:
    """Generate all `A113201(num_vertices)` prime 4-regular spherical shadows directly into a
    contiguous `c_int` buffer of length `count * 4 * num_vertices` via C99 quadrangulation dualization.
    """
    if num_vertices < 2:
        raise ValueError(f"num_vertices must be >= 2, got {num_vertices}")
    lib = _get_lib()
    cap = _A113201_HINTS.get(num_vertices)
    if cap is None:
        cap = int(lib.generate_prime_shadows_c(num_vertices, None, 0))
        if cap < 0:
            raise ValueError(f"num_vertices={num_vertices} exceeds C MAX_QV limit")
    buf = (ctypes.c_int * (cap * 4 * num_vertices))()
    cnt = int(lib.generate_prime_shadows_c(num_vertices, buf, cap))
    if cnt < 0:
        raise ValueError(f"num_vertices={num_vertices} exceeds C MAX_QV limit")
    if cnt > cap:
        buf = (ctypes.c_int * (cnt * 4 * num_vertices))()
        cnt = int(lib.generate_prime_shadows_c(num_vertices, buf, cnt))
    return cnt, buf


def generate_prime_shadows_c_list(num_vertices: int) -> tuple[PlanarShadow, ...]:
    """Return `tuple[PlanarShadow, ...]` generated via the C99 quadrangulation kernel."""
    cnt, buf = generate_prime_shadows_c_buffer(num_vertices)
    stride = 4 * num_vertices
    return tuple(
        PlanarShadow(num_vertices=num_vertices, opp=tuple(int(buf[i * stride + s]) for s in range(stride)))
        for i in range(cnt)
    )


def enumerate_targeted_census_level_candidates(
    ch: int,
    m_min: int | None = None,
    m_max: int | None = None,
    target_components: int = 0,
    target_orientability: int = -1,
    target_chi_min: int = -1000,
    target_chi_max: int = 1000,
    max_bfs_states: int = 2048,
    num_workers: int = 2,
) -> tuple[list[MarkedDiagram], dict[str, int]]:
    """Enumerate certified level-irreducible, 2-cut-free component representatives at `ch`
    filtered in O(n)/O(1) inside C99 BEFORE Kauffman butterfly transforms:
    - `target_components`: `0` = any, `> 0` = exact `|pi_0(F)|`,
    - `target_orientability`: `-1` = any, `1` = Ori, `0` = Non, `2` = Mixed,
    - `target_chi_min..target_chi_max`: inclusive range for `chi(F) = c_+ + c_- - m`.
    """
    if ch == 0:
        d0 = MarkedDiagram.unknot_sphere()
        ok = (target_components in (0, 1)) and (target_orientability in (-1, 1)) and (target_chi_min <= 2 <= target_chi_max)
        return ([d0] if ok else []), {"num_shadows": 0, "raw_cands": int(ok), "components": int(ok), "discarded": 0, "kept": int(ok)}
    if ch < 2:
        return [], {"num_shadows": 0, "raw_cands": 0, "components": 0, "discarded": 0, "kept": 0}

    if m_min is None:
        m_min = 1 if ch == 2 else 2
    if m_max is None:
        m_max = 2 if ch == 2 else max(2, ch - 4)
    if m_min > m_max:
        return [], {"num_shadows": 0, "raw_cands": 0, "components": 0, "discarded": 0, "kept": 0}

    lib = _get_lib()
    num_sh, sh_buf = generate_prime_shadows_c_buffer(ch)
    scan_bfs_states = min(max_bfs_states, 384) if ch >= 9 else max_bfs_states

    def _run_slice(start_idx: int, count_sh: int) -> tuple[list[MarkedDiagram], tuple[int, int, int, int]]:
        if count_sh <= 0:
            return [], (0, 0, 0, 0)
        ptr = ctypes.cast(
            ctypes.byref(sh_buf, start_idx * 4 * ch * ctypes.sizeof(ctypes.c_int)),
            ctypes.POINTER(ctypes.c_int),
        )
        max_kept = 32768
        out_opp = (ctypes.c_int * (max_kept * 4 * ch))()
        out_vt = (ctypes.c_int * (max_kept * ch))()
        out_stats = (ctypes.c_int64 * 4)()
        kept_cnt = int(
            lib.enumerate_targeted_census_level_c(
                ch,
                count_sh,
                ptr,
                m_min,
                m_max,
                target_components,
                target_orientability,
                target_chi_min,
                target_chi_max,
                scan_bfs_states,
                out_opp,
                out_vt,
                max_kept,
                out_stats,
            )
        )
        if kept_cnt > max_kept:
            max_kept = kept_cnt + 256
            out_opp = (ctypes.c_int * (max_kept * 4 * ch))()
            out_vt = (ctypes.c_int * (max_kept * ch))()
            kept_cnt = int(
                lib.enumerate_targeted_census_level_c(
                    ch,
                    count_sh,
                    ptr,
                    m_min,
                    m_max,
                    target_components,
                    target_orientability,
                    target_chi_min,
                    target_chi_max,
                    scan_bfs_states,
                    out_opp,
                    out_vt,
                    max_kept,
                    out_stats,
                )
            )
        diags: list[MarkedDiagram] = []
        for k in range(kept_cnt):
            bo = k * 4 * ch
            bv = k * ch
            diags.append(
                MarkedDiagram(
                    num_vertices=ch,
                    opp=tuple(int(out_opp[bo + i]) for i in range(4 * ch)),
                    vtype=tuple(_VT_MAP[int(out_vt[bv + i])] for i in range(ch)),
                    free_loops=0,
                )
            )
        return diags, (int(out_stats[0]), int(out_stats[1]), int(out_stats[2]), int(out_stats[3]))

    if num_workers <= 1 or num_sh < 32:
        diags, st = _run_slice(0, num_sh)
        final_diags = merge_diagrams_c(ch, diags, max_bfs_states=max(4096, max_bfs_states * 2))
        return final_diags, {
            "num_shadows": num_sh,
            "raw_cands": st[0],
            "components": st[1],
            "discarded": st[2],
            "kept": len(final_diags),
        }

    from concurrent.futures import ThreadPoolExecutor

    # Partition contiguous slices cleanly across `num_workers` threads:
    slices: list[tuple[int, int]] = []
    chunk = (num_sh + num_workers - 1) // num_workers
    for w in range(num_workers):
        s_i = w * chunk
        c_i = max(0, min(num_sh, s_i + chunk) - s_i)
        if c_i > 0:
            slices.append((s_i, c_i))

    merged_raw = 0
    merged_comps = 0
    merged_disc = 0
    shard_diags: list[MarkedDiagram] = []
    with ThreadPoolExecutor(max_workers=len(slices)) as ex:
        for d_list, st in ex.map(lambda sc: _run_slice(sc[0], sc[1]), slices):
            shard_diags.extend(d_list)
            merged_raw += st[0]
            merged_comps += st[1]
            merged_disc += st[2]

    final_diags = merge_diagrams_c(ch, shard_diags, max_bfs_states=max(4096, max_bfs_states * 2))
    return final_diags, {
        "num_shadows": num_sh,
        "raw_cands": merged_raw,
        "components": merged_comps,
        "discarded": merged_disc,
        "kept": len(final_diags),
    }


def enumerate_targeted_prime_ch_census(
    ch: int,
    m_min: int | None = None,
    m_max: int | None = None,
    target_components: int = 0,
    target_orientability: int = -1,
    target_chi_min: int = -1000,
    target_chi_max: int = 1000,
    depth_plus1: int = 256,
    depth_plus2: int = 192,
    max_bfs_states: int = 2048,
    num_workers: int = 2,
    filter_cyclic_tubes: bool = True,
) -> tuple[list[MarkedDiagram], dict[str, int]]:
    """Run targeted level scan + C mountain-pass tunneling + cyclic unknotted-tube pruning."""
    level_diags, stats = enumerate_targeted_census_level_candidates(
        ch=ch,
        m_min=m_min,
        m_max=m_max,
        target_components=target_components,
        target_orientability=target_orientability,
        target_chi_min=target_chi_min,
        target_chi_max=target_chi_max,
        max_bfs_states=max_bfs_states,
        num_workers=num_workers,
    )
    stats = dict(stats)
    stats["level_kept"] = len(level_diags)
    if ch <= 2 or not level_diags:
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


def pack_tangles_c_buffer(tangles: list) -> ctypes.Array:
    """Pack a sequence of `TangleDiagram`s (`b <= 6, c <= 10`) into a flat 24-byte/tangle `c_int8` array."""
    n = len(tangles)
    buf = (ctypes.c_int8 * (24 * max(1, n)))()
    for idx, t in enumerate(tangles):
        base = 24 * idx
        b = t.num_bridges
        c = len(t.braid_word)
        assert 1 <= b <= 6 and 0 <= c <= 10
        buf[base + 0] = b
        buf[base + 1] = c
        for u, v in t.bottom_matching:
            buf[base + 2 + u] = v
            buf[base + 2 + v] = u
        for k, g in enumerate(t.braid_word):
            buf[base + 14 + k] = g
    return buf


def batch_check_tangle_pairs_unlink(
    tangles_a: list,
    tangles_b: list,
    buf_a: ctypes.Array | None = None,
    buf_b: ctypes.Array | None = None,
) -> bytes:
    """Return `len(tangles_a) * len(tangles_b)` status bytes (`1` = certified R1/R2 unlink,
    `0` = certified non-unlink, `2` = needs full R1/R2/R3 `certify_classical_unlink` check)
    for `T_a union bar(T_b)`.
    """
    from surfaceknots.triplane import union_two_tangles_to_classical_link

    n_a = len(tangles_a)
    n_b = len(tangles_b)
    if n_a == 0 or n_b == 0:
        return b""
    lib = _get_lib()
    if buf_a is None:
        buf_a = pack_tangles_c_buffer(tangles_a)
    if buf_b is None:
        buf_b = pack_tangles_c_buffer(tangles_b)
    out = (ctypes.c_uint8 * (n_a * n_b))()
    lib.batch_check_tangle_pairs_c(n_a, buf_a, n_b, buf_b, out)
    raw = bytearray(out)
    # Upgrade any status == 2 via exact Python Reidemeister unlink certification:
    idx = raw.find(2)
    while idx != -1:
        i, j = divmod(idx, n_b)
        l_ij = union_two_tangles_to_classical_link(tangles_a[i], tangles_b[j])
        raw[idx] = 1 if certify_classical_unlink(l_ij).is_unlink else 0
        idx = raw.find(2, idx + 1)
    return bytes(raw)



