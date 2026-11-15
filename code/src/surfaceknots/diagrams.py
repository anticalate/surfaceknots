"""Marked graph diagrams (ch-diagrams), classical link resolutions, and surface topology.

Conventions (matching Yoshikawa 1994, Carter-Kamada-Saito 2004, and Jablonowski 2018):
- A 4-valent planar diagram `D` has `n` vertices `0..n-1` and `free_loops` vertex-free circles.
- Each vertex `v` has 4 incident darts `4*v + s` (`s in {0, 1, 2, 3}`) in counterclockwise (CCW)
  planar cyclic order around `v`. `opp[d]` gives the opposite dart of the edge containing `d`.
- Vertex types (`vtype[v]`):
  * `"X0"`: classical crossing where `(0, 2)` is the UNDER-strand and `(1, 3)` is the OVER-strand
            (matches `X[a, b, c, d]` in EPD notation).
  * `"X1"`: classical crossing where `(1, 3)` is the UNDER-strand and `(0, 2)` is the OVER-strand
            (equivalent to `X[b, c, d, a]` in EPD notation).
  * `"Y"`:  marked saddle vertex where:
            - positive resolution `L_+(D)` pairs `(0, 3)` and `(1, 2)`,
            - negative resolution `L_-(D)` pairs `(0, 1)` and `(2, 3)`.
  * `"Z"`:  marked saddle vertex where:
            - positive resolution `L_+(D)` pairs `(0, 1)` and `(2, 3)`,
            - negative resolution `L_-(D)` pairs `(0, 3)` and `(1, 2)`.
"""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Literal

VertexKind = Literal["X0", "X1", "Y", "Z"]


@dataclass(frozen=True)
class ClassicalLinkDiagram:
    """A planar classical link diagram represented via half-edges (darts) + free circles.

    Every vertex `v` (`0..num_crossings-1`) is a classical crossing of type `"X0"` (`(0,2)` under,
    `(1,3)` over) or `"X1"` (`(1,3)` under, `(0,2)` over), with CCW darts `4*v + 0..3`.
    `free_loops` counts trivial crossing-free circles disjoint from the crossing graph
    (or carrying component tags via `free_loop_tags`).
    """

    num_crossings: int
    opp: tuple[int, ...]
    vtype: tuple[Literal["X0", "X1"], ...]
    free_loops: int = 0
    dart_surface_comp: tuple[int, ...] | None = None
    free_loop_surface_comp: tuple[int, ...] | None = None

    def circle_components(self) -> tuple[list[list[int]], int]:
        """Trace transverse circle components.

        Returns `(cycles, total_circles)` where each entry in `cycles` is the ordered list of
        directed entry-darts traversed along one circle (each crossing traversed twice: entries
        at `d` exit straight across at `d ^ 2`, then jump via `opp[d ^ 2]` to the next entry dart),
        and `total_circles = len(cycles) + free_loops`.
        """
        n = self.num_crossings
        seen = [False] * (4 * n)
        cycles: list[list[int]] = []
        for start in range(4 * n):
            if seen[start]:
                continue
            cyc: list[int] = []
            cur = start
            while not seen[cur]:
                seen[cur] = True
                exit_d = cur ^ 2
                seen[exit_d] = True
                cyc.append(cur)
                cur = self.opp[exit_d]
            cycles.append(cyc)
        return cycles, len(cycles) + self.free_loops

    def num_circles(self) -> int:
        """Total number of closed circle components."""
        _, total = self.circle_components()
        return total

    def crossing_signs_and_circles(self) -> tuple[list[int], list[tuple[int, int]], int]:
        """Orient each circle arbitrarily and return `(signs, (over_comp, under_comp), num_circles)`.

        Convention for crossing sign at `v` with CCW slots `(0, 1, 2, 3)`:
        Let the oriented over-strand enter at slot `s_over_in` and exit at `s_over_out = s_over_in ^ 2`,
        and the oriented under-strand enter at `s_under_in` and exit at `s_under_out = s_under_in ^ 2`.
        In CCW order around `v`:
        - Positive (`+1`) right-hand crossing: if `s_under_in == 0` (going `0 -> 2`), over goes `3 -> 1`
          (`s_over_in == 3`, i.e., `(s_over_in - s_under_in) % 4 == 3`).
        - Negative (`-1`) crossing: `(s_over_in - s_under_in) % 4 == 1`.
        """
        cycles, total = self.circle_components()
        n = self.num_crossings
        dart_dir = [0] * (4 * n)  # +1 if circle enters at d, -1 if circle exits at d
        dart_comp = [-1] * (4 * n)
        for cid, cyc in enumerate(cycles):
            for d in cyc:
                dart_dir[d] = 1
                dart_dir[d ^ 2] = -1
                dart_comp[d] = cid
                dart_comp[d ^ 2] = cid

        signs: list[int] = []
        comp_pairs: list[tuple[int, int]] = []
        for v in range(n):
            vt = self.vtype[v]
            if vt == "X0":
                under_slots = (0, 2)
                over_slots = (1, 3)
            else:
                under_slots = (1, 3)
                over_slots = (0, 2)
            s_under_in = under_slots[0] if dart_dir[4 * v + under_slots[0]] == 1 else under_slots[1]
            s_over_in = over_slots[0] if dart_dir[4 * v + over_slots[0]] == 1 else over_slots[1]
            diff = (s_over_in - s_under_in) & 3
            sign = 1 if diff == 3 else -1
            signs.append(sign)
            comp_pairs.append((dart_comp[4 * v + s_over_in], dart_comp[4 * v + s_under_in]))
        return signs, comp_pairs, total

    def pairwise_linking_numbers_doubled(self) -> dict[tuple[int, int], int]:
        """Return `2 * lk(C_i, C_j)` (`i < j`) for all circle pairs that cross each other."""
        signs, pairs, _ = self.crossing_signs_and_circles()
        lk2: dict[tuple[int, int], int] = {}
        for sgn, (c1, c2) in zip(signs, pairs):
            if c1 != c2:
                key = (c1, c2) if c1 < c2 else (c2, c1)
                lk2[key] = lk2.get(key, 0) + sgn
        return lk2

    def self_writhe_by_surface_component(self, num_surface_comps: int) -> list[int]:
        """Return sum of self-crossing signs for each surface component `0..num_surface_comps-1`.

        Note: self-crossing signs of a circle `C` are invariant under flipping the orientation of `C`,
        and for an unlink resolution `L_+(D)` or `L_-(D)` any two distinct circles have `lk = 0` so
        cross-circle crossings cancel in pairs anyway; restricting to `c1 == c2` (or same surface
        component) makes this manifestly independent of circle orientations.
        """
        signs, pairs, _ = self.crossing_signs_and_circles()
        res = [0] * num_surface_comps
        if self.dart_surface_comp is None:
            return res
        for v, (sgn, (c1, c2)) in enumerate(zip(signs, pairs)):
            sc1 = self.dart_surface_comp[4 * v]
            sc2 = self.dart_surface_comp[4 * v + 1]
            if sc1 == sc2 and c1 == c2:
                res[sc1] += sgn
        return res


@dataclass(frozen=True)
class ComponentTopology:
    """Intrinsic topology and normal Euler number of one connected surface component."""

    orientable: bool
    euler_characteristic: int
    genus_symbol: int  # g >= 0 for orientable Sigma_g; -k < 0 for non-orientable #_k RP^2
    normal_euler: int  # e(F_i)

    @property
    def homeomorphism_type(self) -> str:
        if self.orientable:
            if self.genus_symbol == 0:
                return "S^2"
            if self.genus_symbol == 1:
                return "T^2"
            return f"Sigma_{self.genus_symbol}"
        if self.genus_symbol == -1:
            return "RP^2"
        if self.genus_symbol == -2:
            return "K^2"
        return f"#{abs(self.genus_symbol)}RP^2"


@dataclass(frozen=True)
class SurfaceTopology:
    """Intrinsic topology of a surface-link F(D) reconstructed from its hyperbolic splitting."""

    components: tuple[ComponentTopology, ...]
    num_circles_plus: int
    num_circles_minus: int

    @property
    def num_components(self) -> int:
        return len(self.components)

    @property
    def euler_characteristic(self) -> int:
        return sum(c.euler_characteristic for c in self.components)

    @property
    def orientability_label(self) -> str:
        ors = [c.orientable for c in self.components]
        if all(ors):
            return "Ori"
        if not any(ors):
            return "Non"
        return "Mixed"

    @property
    def yoshikawa_genera(self) -> tuple[int, ...]:
        """Sorted Yoshikawa superscript tuple (e.g. `(0, 1)`, `(0, -2)`, `(-1, -1)`)."""
        # Yoshikawa orders orientable components in increasing genus (0, 0, 1, ...)
        # followed by non-orientable components (-1, -2, ...).
        ori = sorted(c.genus_symbol for c in self.components if c.orientable)
        non = sorted((c.genus_symbol for c in self.components if not c.orientable), reverse=True)
        return tuple(ori + non)


@dataclass(frozen=True)
class MarkedDiagram:
    """A 4-valent planar marked graph diagram `D`."""

    num_vertices: int
    opp: tuple[int, ...]
    vtype: tuple[VertexKind, ...]
    free_loops: int = 0

    @staticmethod
    def unknot_sphere() -> MarkedDiagram:
        """Standard ch=0 diagram `0_1` of the unknotted 2-sphere."""
        return MarkedDiagram(num_vertices=0, opp=(), vtype=(), free_loops=1)

    @property
    def ch_index(self) -> int:
        return self.num_vertices

    @property
    def num_classical(self) -> int:
        return sum(1 for t in self.vtype if t in ("X0", "X1"))

    @property
    def num_marked(self) -> int:
        return sum(1 for t in self.vtype if t in ("Y", "Z"))

    @classmethod
    def from_epd(cls, epd: str | list[tuple[str, int, int, int, int]]) -> MarkedDiagram:
        """Parse an EPD code string or tuple list like `'X[1,5,2,4], Y[5,19,6,20], ...'`.

        Supports `X[a,b,c,d]`, `Y[a,b,c,d]`, `Z[a,b,c,d]`. Every edge label must appear
        exactly twice across all quadruples.
        """
        if isinstance(epd, str):
            s = epd.strip()
            if not s or s == "0_1":
                return cls.unknot_sphere()
            if s == "EMPTY":
                return cls(num_vertices=0, opp=(), vtype=(), free_loops=0)
            token_re = re.compile(
                r"([XYZ])\s*\[\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\]|O\s*\[\s*(\d+)\s*\]"
            )
            items = []
            free_loops = 0
            pos = 0
            n_str = len(s)
            while pos < n_str:
                while pos < n_str and s[pos] in " \t\r\n,;":
                    pos += 1
                if pos >= n_str:
                    break
                m = token_re.match(s, pos)
                if m is None:
                    raise ValueError(f"Malformed token at position {pos} in EPD string: {epd!r}")
                if m.group(1) is not None:
                    items.append(
                        (
                            m.group(1),
                            int(m.group(2)),
                            int(m.group(3)),
                            int(m.group(4)),
                            int(m.group(5)),
                        )
                    )
                else:
                    free_loops += 1
                pos = m.end()
            if not items:
                if free_loops <= 0:
                    raise ValueError(f"Could not parse EPD string: {epd!r}")
                return cls(num_vertices=0, opp=(), vtype=(), free_loops=free_loops)
        else:
            items = [(t, int(a), int(b), int(c), int(d)) for (t, a, b, c, d) in epd]
            free_loops = 0
            if not items:
                return cls.unknot_sphere()

        n = len(items)
        edge_occ: dict[int, list[int]] = {}
        vtypes: list[VertexKind] = []
        for v, (kind, a, b, c, d) in enumerate(items):
            if kind == "X":
                vtypes.append("X0")
            elif kind in ("Y", "Z"):
                vtypes.append(kind)  # type: ignore[arg-type]
            else:
                raise ValueError(f"Unknown EPD vertex kind {kind!r}")
            for s_slot, lbl in enumerate((a, b, c, d)):
                edge_occ.setdefault(lbl, []).append(4 * v + s_slot)

        opp = [0] * (4 * n)
        for lbl, darts in edge_occ.items():
            if len(darts) != 2:
                raise ValueError(f"Edge label {lbl} appears {len(darts)} times in EPD (expected 2)")
            d1, d2 = darts
            opp[d1] = d2
            opp[d2] = d1

        # Verify planar genus-0 ribbon embedding on each connected graph component (F_c - V_c == 2):
        vis_v = [False] * n
        comp_of_v = [-1] * n
        num_gcomps = 0
        comp_sizes: list[int] = []
        for sv in range(n):
            if not vis_v[sv]:
                cid = num_gcomps
                num_gcomps += 1
                vis_v[sv] = True
                comp_of_v[sv] = cid
                stk = [sv]
                cnt = 0
                while stk:
                    u = stk.pop()
                    cnt += 1
                    for sl in range(4):
                        w = opp[4 * u + sl] >> 2
                        if not vis_v[w]:
                            vis_v[w] = True
                            comp_of_v[w] = cid
                            stk.append(w)
                comp_sizes.append(cnt)
        face_vis = [False] * (4 * n)
        comp_faces = [0] * num_gcomps
        for hd in range(4 * n):
            if not face_vis[hd]:
                comp_faces[comp_of_v[hd >> 2]] += 1
                cur = hd
                while not face_vis[cur]:
                    face_vis[cur] = True
                    od = opp[cur]
                    cur = ((od >> 2) << 2) | (((od & 3) - 1) & 3)
        for cid in range(num_gcomps):
            if comp_faces[cid] != comp_sizes[cid] + 2:
                raise ValueError(
                    f"Non-planar 4-valent EPD component (V={comp_sizes[cid]}, F={comp_faces[cid]})"
                )

        return cls(num_vertices=n, opp=tuple(opp), vtype=tuple(vtypes), free_loops=free_loops)

    def to_epd(self) -> str:
        """Serialize to standard EPD string `X[a,b,c,d], Y[...], Z[...]`."""
        if self.num_vertices == 0:
            if self.free_loops == 0:
                return "EMPTY"
            return ", ".join(f"O[{i}]" for i in range(1, self.free_loops + 1))
        n = self.num_vertices
        dart_label = [0] * (4 * n)
        nxt = 1
        for d in range(4 * n):
            if dart_label[d] == 0:
                od = self.opp[d]
                dart_label[d] = nxt
                dart_label[od] = nxt
                nxt += 1
        parts: list[str] = []
        for v in range(n):
            vt = self.vtype[v]
            a, b, c, d = (dart_label[4 * v + s] for s in range(4))
            if vt == "X0":
                parts.append(f"X[{a},{b},{c},{d}]")
            elif vt == "X1":
                # Rotate by 1 CCW so (1,3) under-strand sits at slots (0,2) of X[...]
                parts.append(f"X[{b},{c},{d},{a}]")
            else:
                parts.append(f"{vt}[{a},{b},{c},{d}]")
        for i in range(self.free_loops):
            parts.append(f"O[{nxt + i}]")
        return ", ".join(parts)

    def mirror_3d(self) -> MarkedDiagram:
        """4D mirror image `F^*` via `(x, y, z, t) -> (x, y, -z, t)` (flips classical crossings)."""
        flip = {"X0": "X1", "X1": "X0", "Y": "Y", "Z": "Z"}
        return MarkedDiagram(
            num_vertices=self.num_vertices,
            opp=self.opp,
            vtype=tuple(flip[t] for t in self.vtype),  # type: ignore[arg-type]
            free_loops=self.free_loops,
        )

    def mirror_time(self) -> MarkedDiagram:
        """4D mirror image `F^*` via `(x, y, z, t) -> (x, y, z, -t)` (swaps `Y <-> Z`)."""
        flip = {"X0": "X0", "X1": "X1", "Y": "Z", "Z": "Y"}
        return MarkedDiagram(
            num_vertices=self.num_vertices,
            opp=self.opp,
            vtype=tuple(flip[t] for t in self.vtype),  # type: ignore[arg-type]
            free_loops=self.free_loops,
        )

    def reverse_orientation_zt(self) -> MarkedDiagram:
        """Ambient-isotopic `-F` via `(x, y, z, t) -> (x, y, -z, -t)` (`det = +1` in SO(4))."""
        flip = {"X0": "X1", "X1": "X0", "Y": "Z", "Z": "Y"}
        return MarkedDiagram(
            num_vertices=self.num_vertices,
            opp=self.opp,
            vtype=tuple(flip[t] for t in self.vtype),  # type: ignore[arg-type]
            free_loops=self.free_loops,
        )

    def _dart_surface_components(self) -> tuple[list[int], int]:
        """Partition darts `0..4n-1` into connected 2-manifold components of `F(D)`.

        Across a classical crossing `v` (`X0` or `X1`), strand `(0, 2)` is one sheet and strand
        `(1, 3)` is a (possibly distinct) sheet!
        Across a marked vertex `v` (`Y` or `Z`), all 4 incident darts meet on the connected saddle!
        """
        n = self.num_vertices
        parent = list(range(4 * n))

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x: int, y: int) -> None:
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[ry] = rx

        for d in range(4 * n):
            union(d, self.opp[d])
        for v in range(n):
            vt = self.vtype[v]
            if vt in ("X0", "X1"):
                union(4 * v + 0, 4 * v + 2)
                union(4 * v + 1, 4 * v + 3)
            else:
                union(4 * v + 0, 4 * v + 1)
                union(4 * v + 1, 4 * v + 2)
                union(4 * v + 2, 4 * v + 3)

        root_to_id: dict[int, int] = {}
        dart_comp = [0] * (4 * n)
        for d in range(4 * n):
            r = find(d)
            if r not in root_to_id:
                root_to_id[r] = len(root_to_id)
            dart_comp[d] = root_to_id[r]
        return dart_comp, len(root_to_id) + self.free_loops

    def resolve(self, sign: Literal[1, -1]) -> ClassicalLinkDiagram:
        """Compute positive (`sign = +1`, `L_+(D)`) or negative (`sign = -1`, `L_-(D)`) resolution."""
        n = self.num_vertices
        dart_scomp, _ = self._dart_surface_components()

        # Map classical vertices to 0..c-1
        classical_old_to_new: dict[int, int] = {}
        new_vtypes: list[Literal["X0", "X1"]] = []
        for v in range(n):
            vt = self.vtype[v]
            if vt in ("X0", "X1"):
                classical_old_to_new[v] = len(new_vtypes)
                new_vtypes.append(vt)
        c = len(new_vtypes)

        # At each marked vertex v, define internal splicing permutation `splice[d]` on its 4 darts:
        # Pair (0,3)&(1,2) when (vt=='Y' and sign==+1) or (vt=='Z' and sign==-1).
        # Pair (0,1)&(2,3) when (vt=='Y' and sign==-1) or (vt=='Z' and sign==+1).
        splice: dict[int, int] = {}
        for v in range(n):
            vt = self.vtype[v]
            if vt in ("Y", "Z"):
                pair_03_12 = (vt == "Y" and sign == 1) or (vt == "Z" and sign == -1)
                if pair_03_12:
                    splice[4 * v + 0] = 4 * v + 3
                    splice[4 * v + 3] = 4 * v + 0
                    splice[4 * v + 1] = 4 * v + 2
                    splice[4 * v + 2] = 4 * v + 1
                else:
                    splice[4 * v + 0] = 4 * v + 1
                    splice[4 * v + 1] = 4 * v + 0
                    splice[4 * v + 2] = 4 * v + 3
                    splice[4 * v + 3] = 4 * v + 2

        new_opp = [0] * (4 * c)
        new_dart_scomp = [0] * (4 * c)
        visited_marked_dart = [False] * (4 * n)

        # Trace paths starting from classical darts
        for v_old, v_new in classical_old_to_new.items():
            for s in range(4):
                start_d = 4 * v_old + s
                new_dart_scomp[4 * v_new + s] = dart_scomp[start_d]
                cur = self.opp[start_d]
                while cur in splice:
                    visited_marked_dart[cur] = True
                    nxt_in_v = splice[cur]
                    visited_marked_dart[nxt_in_v] = True
                    cur = self.opp[nxt_in_v]
                w_old = cur >> 2
                ws = cur & 3
                w_new = classical_old_to_new[w_old]
                new_opp[4 * v_new + s] = 4 * w_new + ws

        # Also count closed loops formed purely by marked vertices with 0 classical crossings!
        extra_free_scomps: list[int] = []
        for d in splice:
            if visited_marked_dart[d]:
                continue
            extra_free_scomps.append(dart_scomp[d])
            cur = d
            while not visited_marked_dart[cur]:
                visited_marked_dart[cur] = True
                nxt_in_v = splice[cur]
                visited_marked_dart[nxt_in_v] = True
                cur = self.opp[nxt_in_v]

        base_graph_scomps = max(dart_scomp, default=-1) + 1
        for i in range(self.free_loops):
            extra_free_scomps.append(base_graph_scomps + i)

        return ClassicalLinkDiagram(
            num_crossings=c,
            opp=tuple(new_opp),
            vtype=tuple(new_vtypes),
            free_loops=len(extra_free_scomps),
            dart_surface_comp=tuple(new_dart_scomp),
            free_loop_surface_comp=tuple(extra_free_scomps),
        )

    def surface_topology(self) -> SurfaceTopology:
        """Compute per-component homeomorphism type, orientability, chi, and normal Euler e(F_i).

        Intrinsic topological invariants (`orientable`, `euler_characteristic`, `genus_symbol`)
        are defined on the abstract saddle cobordism between `L_-(D)` and `L_+(D)` for any marked
        diagram `D`. By contrast, `normal_euler` (`e(F_i) = w(L_-(D)|_{F_i}) - w(L_+(D)|_{F_i})`)
        requires vanishing pairwise linking numbers in `L_+(D)` and `L_-(D)` (in particular,
        admissible ch-diagrams whose resolutions bound embedded caps in 4-space) and uses the
        time-first 4-space orientation `R_t x R^3_{x,y,z}` (`dt ^ dx ^ dy ^ dz`).
        """
        n = self.num_vertices
        dart_scomp, total_scomps = self._dart_surface_components()
        graph_scomps = total_scomps - self.free_loops

        l_plus = self.resolve(1)
        l_minus = self.resolve(-1)

        # Count circles of L_+ and L_- belonging to each surface component
        def count_circles_per_scomp(ld: ClassicalLinkDiagram) -> list[int]:
            counts = [0] * total_scomps
            cycles, _ = ld.circle_components()
            assert ld.dart_surface_comp is not None
            for cyc in cycles:
                sc = ld.dart_surface_comp[cyc[0]]
                counts[sc] += 1
            if ld.free_loop_surface_comp is not None:
                for sc in ld.free_loop_surface_comp:
                    counts[sc] += 1
            return counts

        c_plus = count_circles_per_scomp(l_plus)
        c_minus = count_circles_per_scomp(l_minus)

        # Count marked vertices per surface component
        m_counts = [0] * total_scomps
        for v in range(n):
            if self.vtype[v] in ("Y", "Z"):
                m_counts[dart_scomp[4 * v]] += 1

        # Per-component orientability via Proposition 13 (Jablonowski 2018 / Yoshikawa 1994):
        # Assign dir[d] in {+1 (Toward vertex), -1 (From vertex)} on each surface component:
        # - Along each edge: dir[opp[d]] == -dir[d].
        # - Across classical crossing strand: dir[4*v + (s ^ 2)] == -dir[4*v + s].
        # - Around marked vertex: alternating [T, F, T, F] or [F, T, F, T], i.e.:
        #   dir[4*v + 1] == -dir[4*v + 0], dir[4*v + 2] == dir[4*v + 0], dir[4*v + 3] == -dir[4*v + 0].
        comp_orientable = [True] * total_scomps
        dart_dir = [0] * (4 * n)
        for sc in range(graph_scomps):
            # BFS over all darts in surface component `sc`
            for start_d in range(4 * n):
                if dart_scomp[start_d] != sc or dart_dir[start_d] != 0:
                    continue
                dart_dir[start_d] = 1
                queue = [start_d]
                while queue:
                    d = queue.pop()
                    val = dart_dir[d]
                    # 1. Edge constraint: opp[d] has -val
                    od = self.opp[d]
                    if dart_dir[od] == 0:
                        dart_dir[od] = -val
                        queue.append(od)
                    elif dart_dir[od] != -val:
                        comp_orientable[sc] = False
                    # 2. Vertex constraint
                    v = d >> 2
                    s = d & 3
                    vt = self.vtype[v]
                    if vt in ("X0", "X1"):
                        across = (v << 2) | (s ^ 2)
                        if dart_dir[across] == 0:
                            dart_dir[across] = -val
                            queue.append(across)
                        elif dart_dir[across] != -val:
                            comp_orientable[sc] = False
                    else:
                        for t in range(4):
                            td = (v << 2) | t
                            expected = val if ((t - s) & 1) == 0 else -val
                            if dart_dir[td] == 0:
                                dart_dir[td] = expected
                                queue.append(td)
                            elif dart_dir[td] != expected:
                                comp_orientable[sc] = False

        # Per-component normal Euler number: e(F_i) = w(L_-(F_i)) - w(L_+(F_i))
        # (Convention check: for standard P_+ vs P_-, {+2, -2}; orientable always 0.)
        w_plus = l_plus.self_writhe_by_surface_component(total_scomps)
        w_minus = l_minus.self_writhe_by_surface_component(total_scomps)

        comps: list[ComponentTopology] = []
        for sc in range(total_scomps):
            chi = c_plus[sc] + c_minus[sc] - m_counts[sc]
            ori = comp_orientable[sc]
            e_num = w_minus[sc] - w_plus[sc]
            if ori:
                # Orientable closed connected surface: chi = 2 - 2g
                g_sym = (2 - chi) // 2
            else:
                # Non-orientable closed connected surface: chi = 2 - k (represented as -k)
                g_sym = -(2 - chi)
            comps.append(
                ComponentTopology(
                    orientable=ori,
                    euler_characteristic=chi,
                    genus_symbol=g_sym,
                    normal_euler=e_num,
                )
            )

        return SurfaceTopology(
            components=tuple(comps),
            num_circles_plus=sum(c_plus),
            num_circles_minus=sum(c_minus),
        )
