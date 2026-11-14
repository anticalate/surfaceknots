"""Tests for simple planar quadrangulations and prime 4-regular spherical link shadows (OEIS A113201)."""

from __future__ import annotations

from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.moves import find_planar_connected_sum_decomposition
from surfaceknots.quadrangulations import generate_prime_shadows, generate_quadrangulations


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
        assert len(quads) == expected, f"Mismatch at n={n}: got {len(quads)}, expected {expected}"
        assert len(shadows) == expected


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
