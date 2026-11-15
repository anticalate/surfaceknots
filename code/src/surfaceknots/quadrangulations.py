"""Exact generation of simple planar quadrangulations and prime 4-regular spherical link shadows.

By planar duality on S^2 (Brinkmann et al. 2005; Cantarella-Chapman-Mastin 2016; Jablonowski 2018):
- Isomorphism classes of prime (no connected sum, no loop, no 1-gon kink, no nugatory crossing)
  connected 4-regular spherical shadows with n >= 2 vertices (modulo reflection) are in exact
  bijection with isomorphism classes of simple planar quadrangulations Q on V = n + 2 vertices
  (OEIS A113201: 1, 1, 2, 3, 9, 18, 62, 198, 803, 3378, 15882 for n = 2..12).
- Every simple planar quadrangulation on V >= 5 vertices is obtained from C_4 (V = 4) by iterative
  vertex-splitting along pairs of incident edges at a vertex.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache


RotationSystem = tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class PlanarShadow:
    """A connected 4-regular planar/spherical shadow on n vertices.

    Each vertex v in 0..n-1 has 4 incident darts (half-edges) in counterclockwise (CCW)
    order: 4*v + 0, 4*v + 1, 4*v + 2, 4*v + 3.
    `opp[d]` gives the opposite dart paired with dart `d` along an undirected edge
    (`opp[opp[d]] == d` and `opp[d] != d`).
    """

    num_vertices: int
    opp: tuple[int, ...]

    def edges_as_labels(self) -> tuple[tuple[int, int, int, int], ...]:
        """Assign edge labels 1..2n so each vertex has CCW edge tuple (a, b, c, d).

        Edges are numbered consecutively following transverse circular components of the
        flat shadow when possible, matching standard PD/EPD conventions.
        """
        n = self.num_vertices
        if n == 0:
            return ()
        dart_edge = [0] * (4 * n)
        next_label = 1
        for start_d in range(4 * n):
            if dart_edge[start_d] != 0:
                continue
            d = start_d
            while dart_edge[d] == 0:
                od = self.opp[d]
                dart_edge[d] = next_label
                dart_edge[od] = next_label
                next_label += 1
                v = od >> 2
                s = od & 3
                # Continue straight across vertex v: slot (s + 2) mod 4
                d = (v << 2) | ((s + 2) & 3)
        return tuple(
            (dart_edge[4 * v], dart_edge[4 * v + 1], dart_edge[4 * v + 2], dart_edge[4 * v + 3])
            for v in range(n)
        )


def canonical_quadrangulation(adj: RotationSystem, allow_reflection: bool = True) -> RotationSystem:
    """Compute canonical BFS-rotation code of a simple embedded spherical map."""
    n = len(adj)
    best: RotationSystem | None = None
    orientations = (1, -1) if allow_reflection else (1,)
    for v0 in range(n):
        nbrs0 = adj[v0]
        for u0_idx, u0 in enumerate(nbrs0):
            for orient in orientations:
                old_to_new = [-1] * n
                new_to_old = [v0, u0]
                old_to_new[v0] = 0
                old_to_new[u0] = 1
                code: list[tuple[int, ...]] = []
                idx = 0
                while idx < len(new_to_old):
                    v = new_to_old[idx]
                    nbrs = adj[v]
                    d = len(nbrs)
                    if idx == 0:
                        start = u0_idx
                    else:
                        start = 0
                        best_val = 10**9
                        for k in range(d):
                            nv = old_to_new[nbrs[k]]
                            if 0 <= nv < best_val:
                                best_val = nv
                                start = k
                    seq = [0] * d
                    for step in range(d):
                        w = nbrs[(start + orient * step) % d]
                        nw = old_to_new[w]
                        if nw == -1:
                            nw = len(new_to_old)
                            old_to_new[w] = nw
                            new_to_old.append(w)
                        seq[step] = nw
                    code.append(tuple(seq))
                    idx += 1
                cand = tuple(code)
                if best is None or cand < best:
                    best = cand
    assert best is not None
    return best


def _expand_quadrangulation_once(adj: RotationSystem) -> dict[RotationSystem, RotationSystem]:
    """All simple quadrangulations on V+1 vertices obtained by splitting one vertex of `adj`."""
    v_count = len(adj)
    out: dict[RotationSystem, RotationSystem] = {}
    for b in range(v_count):
        nbrs = adj[b]
        d = len(nbrs)
        for i in range(d):
            a = nbrs[i]
            for j in range(i + 1, d):
                c = nbrs[j]
                if a == c:
                    continue
                new_adj = list(adj)
                b1_nbrs = (a,) + nbrs[i + 1 : j] + (c,)
                b2_nbrs = (c,) + nbrs[j + 1 :] + nbrs[:i] + (a,)
                new_adj[b] = b1_nbrs
                new_adj.append(b2_nbrs)
                for w in nbrs[j + 1 :] + nbrs[:i]:
                    new_adj[w] = tuple(v_count if x == b else x for x in new_adj[w])
                idx_a = new_adj[a].index(b)
                new_adj[a] = new_adj[a][:idx_a] + (b, v_count) + new_adj[a][idx_a + 1 :]
                idx_c = new_adj[c].index(b)
                new_adj[c] = new_adj[c][:idx_c] + (v_count, b) + new_adj[c][idx_c + 1 :]
                cand = tuple(new_adj)
                canon = canonical_quadrangulation(cand, allow_reflection=True)
                if canon not in out:
                    out[canon] = canon
    return out


@lru_cache(maxsize=16)
def generate_quadrangulations(num_faces: int) -> tuple[RotationSystem, ...]:
    """Return all simple spherical quadrangulations with `num_faces` faces (V = num_faces + 2).

    Valid for `num_faces >= 2`. Counts match OEIS A113201:
    n=2: 1, n=3: 1, n=4: 2, n=5: 3, n=6: 9, n=7: 18, n=8: 62, n=9: 198, n=10: 803, n=11: 3378.
    """
    if num_faces < 2:
        raise ValueError(f"num_faces must be >= 2, got {num_faces}")
    if num_faces == 2:
        c4: RotationSystem = ((1, 3), (2, 0), (3, 1), (0, 2))
        return (canonical_quadrangulation(c4, allow_reflection=True),)
    prev = generate_quadrangulations(num_faces - 1)
    merged: dict[RotationSystem, RotationSystem] = {}
    for g in prev:
        merged.update(_expand_quadrangulation_once(g))
    return tuple(sorted(merged.keys()))


def quadrangulation_to_shadow(adj: RotationSystem) -> PlanarShadow:
    """Dualize a simple planar quadrangulation Q (V = n + 2, F = n) into a 4-regular PlanarShadow."""
    v_count = len(adj)
    # Directed darts of Q: (u, k) meaning u -> adj[u][k]
    # Next dart on the LEFT face walk: from u -> w (where w = adj[u][k]), take w -> prev_CCW_at_w(u).
    pos_in_nbrs: list[dict[int, int]] = [
        {w: k for k, w in enumerate(adj[u])} for u in range(v_count)
    ]
    visited: dict[tuple[int, int], tuple[int, int]] = {}  # Q-dart -> (dual_vertex, slot 0..3)
    faces: list[list[tuple[int, int]]] = []

    for u in range(v_count):
        for k in range(len(adj[u])):
            start = (u, k)
            if start in visited:
                continue
            face_id = len(faces)
            walk: list[tuple[int, int]] = []
            cur = start
            while cur not in visited:
                visited[cur] = (face_id, len(walk))
                walk.append(cur)
                cu, ck = cur
                cw = adj[cu][ck]
                idx_back = pos_in_nbrs[cw][cu]
                next_k = (idx_back - 1) % len(adj[cw])
                cur = (cw, next_k)
            if len(walk) != 4:
                raise ValueError(f"Expected quadrangulation (face size 4), got face of size {len(walk)}")
            faces.append(walk)

    num_faces = len(faces)
    # Note: walk around a face on the LEFT of Q-darts goes CCW around the face interior!
    # Dual dart (face_id, s) crosses Q-dart (u -> w) from face_id to opposite face across (w -> u).
    opp = [0] * (4 * num_faces)
    for (u, k), (fid, s) in visited.items():
        w = adj[u][k]
        rev_k = pos_in_nbrs[w][u]
        ofid, os = visited[(w, rev_k)]
        opp[4 * fid + s] = 4 * ofid + os
    return PlanarShadow(num_vertices=num_faces, opp=tuple(opp))


@lru_cache(maxsize=16)
def generate_prime_shadows(num_vertices: int) -> tuple[PlanarShadow, ...]:
    """Return all prime 4-regular spherical link shadows with `num_vertices` crossings (n >= 2)."""
    if num_vertices >= 9:
        try:
            from surfaceknots.fast_enum import generate_prime_shadows_c_list

            return generate_prime_shadows_c_list(num_vertices)
        except Exception:
            pass
    return tuple(quadrangulation_to_shadow(q) for q in generate_quadrangulations(num_vertices))

