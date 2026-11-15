"""Tests for simple planar quadrangulations, prime 4-regular spherical link shadows (OEIS A113201), and C99 kernel guards."""

from __future__ import annotations

import ctypes
import pytest

from surfaceknots import fast_enum as fe
from surfaceknots.ch11_census import disjoint_union, planar_connected_sum
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.moves import find_planar_connected_sum_decomposition
from surfaceknots.quadrangulations import generate_prime_shadows, generate_quadrangulations
from surfaceknots.yoshikawa_data import YOSHIKAWA_BY_LABEL


OEIS_A113201 = {
    2: 1,
    3: 1,
    4: 2,
    5: 3,
    6: 9,
    7: 18,
    8: 62,
    9: 198,
}


def test_oeis_a113201_counts() -> None:
    for n, expected in OEIS_A113201.items():
        quads = generate_quadrangulations(n)
        shadows = generate_prime_shadows(n)
        c_shadows = fe.generate_prime_shadows_c_list(n)
        assert len(quads) == expected, f"Mismatch at n={n}: got {len(quads)}, expected {expected}"
        assert len(shadows) == expected
        assert len(c_shadows) == expected


def test_prime_shadows_have_no_kinks_or_2_edge_cuts() -> None:
    for n in range(2, 8):
        for sh in generate_prime_shadows(n):
            assert len(sh.opp) == 4 * n
            # Involution without fixed points or 1-gons
            for v in range(n):
                for s in range(4):
                    d = 4 * v + s
                    od = sh.opp[d]
                    assert od != d
                    assert sh.opp[od] == d
                    assert od != 4 * v + ((s + 1) & 3)
            if n >= 4:
                diag = MarkedDiagram(num_vertices=n, opp=sh.opp, vtype=("X0",) * n)
                assert find_planar_connected_sum_decomposition(diag) is None


def test_c_shadow_and_fast_enum_bounds_and_validation() -> None:
    with pytest.raises(ValueError):
        fe.generate_prime_shadows_c_list(1)
    with pytest.raises(ValueError):
        fe.generate_prime_shadows_c_list(21)

    for bad_ch in (-1, 21, 24):
        with pytest.raises(ValueError):
            fe.enumerate_admissible_ch_diagrams(bad_ch)
        with pytest.raises(ValueError):
            fe.merge_diagrams_c(bad_ch, [])
        with pytest.raises(ValueError):
            fe.tunnel_diagrams_c(bad_ch, [])
        with pytest.raises(ValueError):
            fe.enumerate_targeted_census_level_candidates(bad_ch)

    # c > 15 (ch=16, m=0 -> c=16) fails loudly instead of UB shift / range(-2):
    sh16 = [0] * 64
    for v in range(16):
        nv = (v + 1) % 16
        sh16[4 * v + 0] = 4 * nv + 2
        sh16[4 * nv + 2] = 4 * v + 0
        sh16[4 * v + 1] = 4 * nv + 3
        sh16[4 * nv + 3] = 4 * v + 1
    with pytest.raises(ValueError, match="exceeds C kernel limit"):
        fe.enumerate_admissible_ch_diagrams(16, m_min=0, m_max=0)
    with pytest.raises(RuntimeError, match="enumerate_census_level_c failed"):
        fe._run_c_census_shard((16, [tuple(sh16)], 0, 0, 64, False))

    # Shape mismatch and free_loops != 0 are rejected before flattening d.opp:
    t2 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^1"].canonical_epd)
    with pytest.raises(ValueError, match="num_vertices"):
        fe.merge_diagrams_c(6, [t2])
    with pytest.raises(ValueError, match="num_vertices"):
        fe.tunnel_diagrams_c(6, [t2])

    t2_fl = MarkedDiagram(t2.num_vertices, t2.opp, t2.vtype, free_loops=1)
    with pytest.raises(ValueError, match="free_loops"):
        fe.merge_diagrams_c(2, [t2_fl])
    with pytest.raises(ValueError, match="free_loops"):
        fe.tunnel_diagrams_c(2, [t2_fl])


def test_disconnected_and_step0_reducible_composite_filtering_c() -> None:
    t2 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["2_1^1"].canonical_epd)
    kink1 = MarkedDiagram(1, (1, 0, 3, 2), ("X0",))
    d_split2 = disjoint_union(kink1, kink1)
    d_split3 = disjoint_union(kink1, t2)
    d_split4 = disjoint_union(t2, t2)

    for d_sp in (d_split2, d_split3, d_split4):
        n_sp = d_sp.num_vertices
        assert fe.merge_diagrams_c(n_sp, [d_sp]) == []
        assert fe.tunnel_diagrams_c(n_sp, [d_sp], depth_plus1=16, depth_plus2=16) == []
        kept_shard, _ = fe._run_c_census_shard((n_sp, [d_sp.opp], 0, n_sp, 64, False))
        assert kept_shard == []

    # Connected diagram in a mixed batch with disconnected diagram survives cleanly:
    d_hopf = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["6_1^0,1"].canonical_epd)
    d_split6 = disjoint_union(t2, d_split4)
    assert len(fe.merge_diagrams_c(6, [d_split6, d_hopf])) == 1
    assert len(fe.tunnel_diagrams_c(6, [d_split6, d_hopf], depth_plus1=16, depth_plus2=16)) == 1

    # Step-0 filtering at ch=20 when nnbr == 0 (composite 10_3 # 10_3 and 1-gon kinked diagram):
    d10_3 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_3"].canonical_epd)
    d20_comp = planar_connected_sum(d10_3, 0, d10_3, 0)
    assert d20_comp.num_vertices == 20
    assert fe.tunnel_diagrams_c(20, [d20_comp], depth_plus1=96, depth_plus2=64, num_workers=2) == []

    d9_1 = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["9_1"].canonical_epd)
    d19 = planar_connected_sum(d10_3, 0, d9_1, 0)
    opp20 = list(d19.opp) + [0, 0, 0, 0]
    od0 = opp20[0]
    u = 19
    opp20[4 * u + 0] = 4 * u + 1
    opp20[4 * u + 1] = 4 * u + 0
    opp20[4 * u + 2] = od0
    opp20[od0] = 4 * u + 2
    opp20[4 * u + 3] = 0
    opp20[0] = 4 * u + 3
    d20_kink = MarkedDiagram(20, tuple(opp20), d19.vtype + ("X0",))
    assert fe.tunnel_diagrams_c(20, [d20_kink], depth_plus1=96, depth_plus2=64, num_workers=2) == []


def test_triplane_c_hashtable_bounds_overflow_and_normal_euler_guard() -> None:
    lib = fe._get_lib()
    nb = ctypes.c_int(0)
    novf = ctypes.c_int(0)
    stats = (ctypes.c_uint64 * 4)()
    b_small = (fe.TriplaneBucketOut * 2)()

    # 1. max_buckets > TP_HT_SIZE // 2 (524_288) fails loudly with -1; max_buckets = 500_000 succeeds:
    rc_bad = lib.scan_triplane_level_c(
        2, 0, 0, 0, 0, 0, 0, None, None, None, 0, 0,
        None, None, None, None, None, None, None, None, None, None, None,
        None, None, None, None, None, stats, b_small, 524_289,
        ctypes.byref(nb), None, 0, ctypes.byref(novf),
    )
    assert rc_bad == -1 and nb.value == 524_290

    rc_ok = lib.eval_triplane_triples_list_c(
        2, 0, None, None, None, 0, None, b_small, 500_000,
        ctypes.byref(nb), None, 0, ctypes.byref(novf),
    )
    assert rc_ok == 0 and nb.value == 0

    # 2. Construct non-orientable b=2 tangles testing |e|=32 overflow fallback vs |e|=30 resolution,
    #    and testing bucket-cap overflow reporting (num_buckets > max_buckets):
    meta0 = (fe.TangleFastMeta * 3)()  # [0]=m1(10 xings), [1]=m2(6 xings -> |e|=32), [2]=m2_30(5 xings -> |e|=30)
    meta6 = (fe.TangleFastMeta * 1)()  # [0]=m3(0 xings)
    for p, v in enumerate((1, 0, 3, 2)):
        meta0[0].tau[p] = v
    for p, v in enumerate((2, 3, 0, 1)):
        meta0[1].tau[p] = v
        meta0[2].tau[p] = v
    for p, v in enumerate((3, 2, 1, 0)):
        meta6[0].tau[p] = v
    meta0[0].is_entry_mask = 0b0101
    meta0[1].is_entry_mask = 0b0011
    meta0[2].is_entry_mask = 0b0011
    meta6[0].is_entry_mask = 0b0011

    meta0[0].nxings = 10
    for x in range(10):
        meta0[0].xing_pL[x] = 0
        meta0[0].xing_pR[x] = 2
        meta0[0].xing_bs[x] = -1
    meta0[1].nxings = 6
    for x in range(6):
        meta0[1].xing_pL[x] = 0
        meta0[1].xing_pR[x] = 1
        meta0[1].xing_bs[x] = 1
    meta0[2].nxings = 5
    for x in range(5):
        meta0[2].xing_pL[x] = 0
        meta0[2].xing_pR[x] = 1
        meta0[2].xing_bs[x] = 1

    conc_lut = (ctypes.c_int8 * 9)(*([-1] * 9))
    tr_e32 = (fe.TriplaneOverflowTriple * 1)(fe.TriplaneOverflowTriple(0, 1, 0))
    ovf_out = (fe.TriplaneOverflowTriple * 2)()

    rc_32 = lib.eval_triplane_triples_list_c(
        2, 3, meta0, meta6, conc_lut, 1, tr_e32,
        b_small, 2, ctypes.byref(nb), ovf_out, 2, ctypes.byref(novf),
    )
    assert rc_32 == 0
    assert nb.value == 0 and novf.value == 1

    # Compare two distinct C-resolvable triples ((0, 2, 0) with |e|=30 and (2, 2, 0) with different topology)
    # under max_buckets = 1 -> reports nb.value == 2 > max_buckets without OOB write!
    tr_two = (fe.TriplaneOverflowTriple * 2)(
        fe.TriplaneOverflowTriple(0, 2, 0),
        fe.TriplaneOverflowTriple(1, 2, 0),
    )
    b_cap1 = (fe.TriplaneBucketOut * 1)()
    rc_cap = lib.eval_triplane_triples_list_c(
        2, 3, meta0, meta6, conc_lut, 2, tr_two,
        b_cap1, 1, ctypes.byref(nb), ovf_out, 2, ctypes.byref(novf),
    )
    assert rc_cap == 0
    assert nb.value == 2 and novf.value == 0
    assert (b_cap1[0].key_hi & 31) == 30
