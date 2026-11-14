"""Automated planar layout and SVG / TikZ vector rendering for marked graph diagrams (`ch`-diagrams).

Implements a theorem-backed planar layout pipeline:
1. Connected 4-regular block decomposition + surface-component coloring (`_dart_surface_components`).
2. Doubly-subdivided primal/face-star triangulation `G_aug` + outer-face selection.
3. Iteratively reweighted `p`-Laplacian Tutte harmonic embedding on the unit disk (strictly
   positive spring weights guarantee convex interior faces and zero edge crossings), followed by
   exterior polar puffing of outer-face arcs into `r > 1` for rounded knot-table lobes.
4. Orthogonal `90`-degree crossing frames, `C^1`-continuous cubic Bezier splines, geometric
   under-crossing gaps, and Yoshikawa/Lomonaco `Y`/`Z` saddle marker bars.
5. Export to crisp `.svg` vector cards and standalone compilable LaTeX/TikZ `.tex` figures.
"""

from __future__ import annotations

from dataclasses import dataclass
import html
import json
import math
from pathlib import Path

from surfaceknots.diagrams import MarkedDiagram


# Colorblind-safe, high-contrast academic ink palette (Okabe-Ito / KnotInfo inspired)
COMPONENT_PALETTE: tuple[str, ...] = (
    "#1E3A8A",  # Component 0: Deep sapphire indigo
    "#C2410C",  # Component 1: Vermilion / burnt sienna
    "#0F766E",  # Component 2: Deep viridian teal
    "#7E22CE",  # Component 3: Royal violet
    "#B45309",  # Component 4: Warm ochre
    "#BE123C",  # Component 5: Crimson
)

TIKZ_PALETTE_RGB: tuple[tuple[int, int, int], ...] = (
    (30, 58, 138),
    (194, 65, 12),
    (15, 118, 110),
    (126, 34, 206),
    (180, 83, 9),
    (190, 18, 60),
)


@dataclass(frozen=True)
class CubicBezierSegment:
    """A single cubic Bezier curve `p0 -> p1 -> p2 -> p3` in mathematical `(x, y)` coordinates."""

    p0: tuple[float, float]
    p1: tuple[float, float]
    p2: tuple[float, float]
    p3: tuple[float, float]
    surface_component: int
    layer: int  # 0 = under-strand, 1 = over-strand or marked-strand


@dataclass(frozen=True)
class SaddleMarkerBar:
    """Orthogonal saddle marker bar at a `Y` or `Z` marked vertex."""

    center: tuple[float, float]
    a: tuple[float, float]
    b: tuple[float, float]
    kind: str  # "Y" or "Z"
    surface_component: int


@dataclass(frozen=True)
class FreeCirclePrimitive:
    """Vertex-free circle component (`O[1]`)."""

    center: tuple[float, float]
    radius: float
    surface_component: int


@dataclass(frozen=True)
class DiagramLayout:
    """Complete geometric layout of a `MarkedDiagram` in mathematical `(x, y)` coordinates."""

    vertex_positions: tuple[tuple[float, float], ...]
    vertex_angles: tuple[float, ...]
    beziers: tuple[CubicBezierSegment, ...]
    markers: tuple[SaddleMarkerBar, ...]
    circles: tuple[FreeCirclePrimitive, ...]
    bbox: tuple[float, float, float, float]  # (xmin, ymin, xmax, ymax)
    num_surface_components: int


def label_to_stem(label: str) -> str:
    """Convert a surface-link catalog label (e.g. `'10_1^0,0,1'`, `'7_1^0,-2'`, `'11_1.a'`)
    into a collision-free, filesystem- and TeX-safe stem (`[0-9A-Za-z_]+`).
    """
    s = label.strip()
    suffix = ""
    if "." in s:
        s, sub = s.split(".", 1)
        suffix = f"_{sub}"
    if "^" not in s:
        base_clean = "".join(ch if ch.isalnum() else "_" for ch in s)
        return f"{base_clean}{suffix}"
    base, sup = s.split("^", 1)
    base_clean = "".join(ch if ch.isalnum() else "_" for ch in base)
    parts: list[str] = []
    for tok in sup.split(","):
        tok = tok.strip()
        if tok.startswith("-"):
            parts.append(f"m{tok[1:]}")
        else:
            parts.append(tok)
    return f"{base_clean}_sup_{'_'.join(parts)}{suffix}"


def _solve_linear_system(a_mat: list[list[float]], rhs_x: list[float], rhs_y: list[float]) -> tuple[list[float], list[float]]:
    """Solve `A * [x, y] = [rhs_x, rhs_y]` via Gaussian elimination with partial pivoting."""
    n = len(a_mat)
    if n == 0:
        return [], []
    aug = [list(a_mat[i]) + [rhs_x[i], rhs_y[i]] for i in range(n)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[piv][col]) < 1e-13:
            raise ValueError("Singular Tutte harmonic system (disconnected or degenerate planar block)")
        if piv != col:
            aug[col], aug[piv] = aug[piv], aug[col]
        inv_diag = 1.0 / aug[col][col]
        row_c = aug[col]
        for j in range(col, n + 2):
            row_c[j] *= inv_diag
        for r in range(n):
            if r == col:
                continue
            factor = aug[r][col]
            if abs(factor) > 1e-16:
                row_r = aug[r]
                for j in range(col, n + 2):
                    row_r[j] -= factor * row_c[j]
    return [aug[i][n] for i in range(n)], [aug[i][n + 1] for i in range(n)]


def _norm2(vx: float, vy: float) -> tuple[float, float]:
    d = math.hypot(vx, vy)
    if d < 1e-12:
        return (1.0, 0.0)
    return (vx / d, vy / d)


def _verify_ccw_cyclic_order(
    px: float,
    py: float,
    qs: list[tuple[float, float]],
) -> bool:
    """Check that the 4 points `qs[0..3]` appear in strict CCW angular order around `(px, py)`."""
    angs = [math.atan2(qy - py, qx - px) for qx, qy in qs]
    # Sum of positive CCW angular increments `((angs[(s+1)%4] - angs[s]) mod 2*pi)` must equal `2*pi`!
    total = 0.0
    for s in range(4):
        diff = (angs[(s + 1) & 3] - angs[s]) % (2.0 * math.pi)
        if diff < 1e-5 or diff > 2.0 * math.pi - 1e-5:
            return False
        total += diff
    return abs(total - 2.0 * math.pi) < 1e-3


def _layout_connected_block(
    m: int,
    local_opp: list[int],
    local_vtype: list[str],
    local_dart_scomp: list[int],
) -> tuple[
    list[tuple[float, float]],
    list[float],
    list[CubicBezierSegment],
    list[SaddleMarkerBar],
]:
    """Compute crossing-free planar coordinates and smooth cubic Bezier curves for a connected
    `m`-vertex 4-regular planar block (`m >= 1`).
    """
    if m == 1:
        # 1-vertex figure-eight / double-loop block (`4` darts paired in `local_opp`)
        pos_v = [(0.0, 0.0)]
        phi = 0.0
        # Check whether pairing is adjacent `(0,1),(2,3)` or `(0,3),(1,2)` (planar!).
        if local_opp[0] == 2:
            raise ValueError("Non-planar 1-vertex crossing loop opp=(2,3,0,1)")
        pos_u: list[tuple[float, float]] = []
        for s in range(4):
            ang = phi + s * (0.5 * math.pi) + 0.25 * math.pi * 0.0
            pos_u.append((0.65 * math.cos(ang), 0.65 * math.sin(ang)))
        # Push loop midpoints out cleanly
        all_pos = pos_v + pos_u
        return _build_splines_from_positions(m, local_opp, local_vtype, local_dart_scomp, all_pos)

    # 1. Extract combinatorial left-faces (`next_dart(hd) = (opp[hd] & ~3) | (((opp[hd] & 3) - 1) & 3)`)
    visited = [False] * (4 * m)
    faces: list[list[int]] = []
    for start in range(4 * m):
        if visited[start]:
            continue
        cyc: list[int] = []
        cur = start
        while not visited[cur]:
            visited[cur] = True
            cyc.append(cur)
            od = local_opp[cur]
            cur = (od & ~3) | (((od & 3) - 1) & 3)
        faces.append(cyc)

    # Euler characteristic check on connected 4-regular spherical map: `V - E + F = m - 2m + F = F - m == 2`
    if len(faces) != m + 2:
        raise ValueError(
            f"Connected 4-regular block is non-planar: V={m}, E={2 * m}, F={len(faces)} (expected F={m + 2})"
        )

    # 2. Build doubly-subdivided boundary cycles `W_f` (`0..m-1` = original vertices; `m+hd` = dart nodes)
    face_walks: list[list[int]] = []
    for cyc in faces:
        w_f: list[int] = []
        for hd in cyc:
            w_f.append(hd >> 2)
            w_f.append(m + hd)
            w_f.append(m + local_opp[hd])
        face_walks.append(w_f)

    # 3. Select outer face maximizing `(is_simple_cycle, degree, marked_count, -f_idx)`
    def _face_score(f_idx: int) -> tuple[int, int, int, int]:
        w_f = face_walks[f_idx]
        is_simple = 1 if len(set(w_f)) == len(w_f) else 0
        marked_cnt = sum(1 for hd in faces[f_idx] if local_vtype[hd >> 2] in ("Y", "Z"))
        return (is_simple, len(faces[f_idx]), marked_cnt, -f_idx)

    outer_f_idx = max(range(len(faces)), key=_face_score)
    w_out_raw = face_walks[outer_f_idx]
    # Deduplicate while preserving order if cut-vertex ever repeats
    seen_out: set[int] = set()
    w_out: list[int] = []
    for nd in w_out_raw:
        if nd not in seen_out:
            seen_out.add(nd)
            w_out.append(nd)

    interior_faces = [fi for fi in range(len(faces)) if fi != outer_f_idx]
    total_aug_nodes = 5 * m + len(interior_faces)

    # Place outer walk CLOCKWISE on the unit circle (`theta_j = pi - 2*pi*j/L`)
    fixed_pos: dict[int, tuple[float, float]] = {}
    l_out = len(w_out)
    for j, nd in enumerate(w_out):
        theta = math.pi - (2.0 * math.pi * j) / l_out
        fixed_pos[nd] = (math.cos(theta), math.sin(theta))

    # Build base weighted edge list of `G_aug`
    # `(u, v, base_weight, is_primal)`
    aug_edges: list[tuple[int, int, float, bool]] = []
    for hd in range(4 * m):
        aug_edges.append((hd >> 2, m + hd, 1.35, True))
        od = local_opp[hd]
        if hd < od:
            aug_edges.append((m + hd, m + od, 1.0, True))

    for idx_in, fi in enumerate(interior_faces):
        cf = 5 * m + idx_in
        deg_f = len(faces[fi])
        for nd in dict.fromkeys(face_walks[fi]):
            if deg_f <= 2:
                w_rad = 0.85 if nd < m else 0.65
            else:
                w_rad = 0.42 if nd < m else 0.28
            aug_edges.append((cf, nd, w_rad, False))

    free_nodes = [nd for nd in range(total_aug_nodes) if nd not in fixed_pos]
    free_idx = {nd: i for i, nd in enumerate(free_nodes)}
    n_free = len(free_nodes)

    cur_x = [0.0] * total_aug_nodes
    cur_y = [0.0] * total_aug_nodes
    for nd, (fx, fy) in fixed_pos.items():
        cur_x[nd] = fx
        cur_y[nd] = fy

    # 4. IRLS `p`-Laplacian Tutte harmonic solver (`p = 3.3`: strictly positive weights -> 100% crossing-free!)
    for it in range(10):
        mat = [[0.0] * n_free for _ in range(n_free)]
        rx = [0.0] * n_free
        ry = [0.0] * n_free
        for u, v, bw, is_prim in aug_edges:
            if it == 0:
                w = bw
            else:
                dist = math.hypot(cur_x[u] - cur_x[v], cur_y[u] - cur_y[v])
                exp_p = 1.35 if is_prim else 1.05
                w = bw * ((dist + 0.035) ** exp_p)
            iu = free_idx.get(u)
            iv = free_idx.get(v)
            if iu is not None:
                mat[iu][iu] += w
                if iv is not None:
                    mat[iu][iv] -= w
                else:
                    rx[iu] += w * cur_x[v]
                    ry[iu] += w * cur_y[v]
            if iv is not None:
                mat[iv][iv] += w
                if iu is not None:
                    mat[iv][iu] -= w
                else:
                    rx[iv] += w * cur_x[u]
                    ry[iv] += w * cur_y[u]

        sx, sy = _solve_linear_system(mat, rx, ry)
        for nd, i in free_idx.items():
            cur_x[nd] = sx[i]
            cur_y[nd] = sy[i]

    # 5. Exterior polar puffing of outer-face subdivision nodes into `r > 1`:
    # All interior nodes satisfy `r < 1`; pushing outer-face subdivision nodes `m + hd`, `m + opp[hd]`
    # along their existing polar rays to `r_puff` opens natural ~90-degree X-crossings on the boundary
    # without intersecting any `r < 1` interior edge or disjoint polar sector!
    cos_step = max(math.cos((2.0 * math.pi) / max(l_out, 6)), 0.68)
    r_puff = min(1.58, 1.14 / cos_step)
    for hd in faces[outer_f_idx]:
        for nd in (m + hd, m + local_opp[hd]):
            if nd in fixed_pos:
                cur_x[nd] *= r_puff
                cur_y[nd] *= r_puff

    sub_pos = [(cur_x[i], cur_y[i]) for i in range(5 * m)]

    # Verify strict CCW dart order `(0, 1, 2, 3)` at every vertex!
    for v in range(m):
        vx, vy = sub_pos[v]
        qs = [sub_pos[m + 4 * v + s] for s in range(4)]
        if not _verify_ccw_cyclic_order(vx, vy, qs):
            raise RuntimeError(f"CCW cyclic dart order violated at vertex {v}")

    return _build_splines_from_positions(m, local_opp, local_vtype, local_dart_scomp, sub_pos)


def _build_splines_from_positions(
    m: int,
    local_opp: list[int],
    local_vtype: list[str],
    local_dart_scomp: list[int],
    sub_pos: list[tuple[float, float]],
) -> tuple[
    list[tuple[float, float]],
    list[float],
    list[CubicBezierSegment],
    list[SaddleMarkerBar],
]:
    """Fit orthogonal crossing frames, under-crossing gaps, cubic Bezier splines, and `Y`/`Z` bars."""
    v_pos = sub_pos[:m]
    u_pos = sub_pos[m : 5 * m]

    # Characteristic diagram radius
    xs = [p[0] for p in sub_pos]
    ys = [p[1] for p in sub_pos]
    diag_span = max(max(xs) - min(xs), max(ys) - min(ys), 1.0)

    v_angles: list[float] = []
    v_min_arm: list[float] = []
    for v in range(m):
        vx, vy = v_pos[v]
        csum, ssum = 0.0, 0.0
        min_d = float("inf")
        for s in range(4):
            ux, uy = u_pos[4 * v + s]
            dx, dy = ux - vx, uy - vy
            d = math.hypot(dx, dy)
            min_d = min(min_d, d)
            ang = math.atan2(dy, dx) - s * (0.5 * math.pi)
            csum += math.cos(ang)
            ssum += math.sin(ang)
        v_angles.append(math.atan2(ssum, csum))
        v_min_arm.append(min_d)

    beziers: list[CubicBezierSegment] = []
    markers: list[SaddleMarkerBar] = []

    # Build two `C^1`-joined Bezier segments for each undirected edge `hd < od`
    # Since `mid = 0.5 * (qh + qo)`, setting `p2 = qh` on the `hd` side and `p2 = qo` on the `od`
    # side gives `B'_hd(1) = 3*(mid - qh) = 1.5*(qo - qh) = -B'_od(1)`: exact `C^1` continuity
    # at `mid` while keeping every spline strictly inside its disjoint Tutte edge corridor!
    for hd in range(4 * m):
        od = local_opp[hd]
        if hd >= od:
            continue
        v = hd >> 2
        s = hd & 3
        w = od >> 2
        os = od & 3
        scomp = local_dart_scomp[hd]

        qh = u_pos[hd]
        qo = u_pos[od]
        mid = (0.5 * (qh[0] + qo[0]), 0.5 * (qh[1] + qo[1]))

        for cur_v, cur_s, cur_q in (
            (v, s, qh),
            (w, os, qo),
        ):
            vt = local_vtype[cur_v]
            is_under = (vt == "X0" and (cur_s & 1) == 0) or (vt == "X1" and (cur_s & 1) == 1)
            p_center = v_pos[cur_v]
            nat_dx = cur_q[0] - p_center[0]
            nat_dy = cur_q[1] - p_center[1]
            nat_d = math.hypot(nat_dx, nat_dy)
            nat_u = _norm2(nat_dx, nat_dy)

            phi_slot = v_angles[cur_v] + cur_s * (0.5 * math.pi)
            orth_u = (math.cos(phi_slot), math.sin(phi_slot))
            # If the orthogonal ray deviates by > 50 deg from the Tutte corridor ray `nat_u`
            # (e.g. on tight 2-vertex multi-edge bigons like `2_1^1`), clamp toward `nat_u`
            # so `p1` never leaves the private face sector of `cur_q`!
            dot_uo = orth_u[0] * nat_u[0] + orth_u[1] * nat_u[1]
            if dot_uo < 0.65:
                uv = _norm2(0.45 * orth_u[0] + 0.55 * nat_u[0], 0.45 * orth_u[1] + 0.55 * nat_u[1])
                # Keep opposite slots `(cur_s, cur_s ^ 2)` collinear at classical crossings:
                opp_q = u_pos[4 * cur_v + (cur_s ^ 2)]
                opp_u = _norm2(opp_q[0] - p_center[0], opp_q[1] - p_center[1])
                axis_u = _norm2(nat_u[0] - opp_u[0], nat_u[1] - opp_u[1])
                if axis_u[0] * nat_u[0] + axis_u[1] * nat_u[1] > 0.35:
                    uv = axis_u
            else:
                uv = orth_u

            gap_r = min(0.32 * v_min_arm[cur_v], 0.056 * diag_span)
            if is_under:
                p0 = (p_center[0] + gap_r * uv[0], p_center[1] + gap_r * uv[1])
                layer = 0
            else:
                p0 = p_center
                layer = 1

            l1 = min(0.34 * nat_d, 0.36 * v_min_arm[cur_v])
            p1 = (p0[0] + l1 * uv[0], p0[1] + l1 * uv[1])
            p2 = cur_q

            beziers.append(
                CubicBezierSegment(
                    p0=p0,
                    p1=p1,
                    p2=p2,
                    p3=mid,
                    surface_component=scomp,
                    layer=layer,
                )
            )

    # Build marked saddle bars (`Y`: `phi + pi/4` spanning `(0,1)-(2,3)`; `Z`: `phi + 3*pi/4` spanning `(1,2)-(3,0)`)
    for v in range(m):
        vt = local_vtype[v]
        if vt not in ("Y", "Z"):
            continue
        # Bisect actual quadrant `(0,1)` for `Y` or `(1,2)` for `Z` directly so bar is always centered between branches!
        vx, vy = v_pos[v]
        s_a, s_b = (0, 1) if vt == "Y" else (1, 2)
        ua = _norm2(u_pos[4 * v + s_a][0] - vx, u_pos[4 * v + s_a][1] - vy)
        ub = _norm2(u_pos[4 * v + s_b][0] - vx, u_pos[4 * v + s_b][1] - vy)
        bis = _norm2(ua[0] + ub[0], ua[1] + ub[1])
        bar_r = min(0.32 * v_min_arm[v], 0.058 * diag_span)
        dx, dy = bar_r * bis[0], bar_r * bis[1]
        markers.append(
            SaddleMarkerBar(
                center=(vx, vy),
                a=(vx - dx, vy - dy),
                b=(vx + dx, vy + dy),
                kind=vt,
                surface_component=local_dart_scomp[4 * v],
            )
        )

    return v_pos, v_angles, beziers, markers


def compute_diagram_layout(d: MarkedDiagram) -> DiagramLayout:
    """Compute a complete, crossing-free `DiagramLayout` for any planar `MarkedDiagram`."""
    n = d.num_vertices
    dart_scomp, total_scomps = d._dart_surface_components()
    base_graph_scomps = max(dart_scomp, default=-1) + 1

    # Partition vertices `0..n-1` into connected 4-regular graph components
    comp_of_v = [-1] * n
    blocks: list[list[int]] = []
    for start_v in range(n):
        if comp_of_v[start_v] != -1:
            continue
        cid = len(blocks)
        b_vs = [start_v]
        comp_of_v[start_v] = cid
        q = [start_v]
        while q:
            u = q.pop()
            for s in range(4):
                w = d.opp[4 * u + s] >> 2
                if comp_of_v[w] == -1:
                    comp_of_v[w] = cid
                    b_vs.append(w)
                    q.append(w)
        blocks.append(sorted(b_vs))

    global_v_pos: list[tuple[float, float]] = [(0.0, 0.0)] * n
    global_v_ang: list[float] = [0.0] * n
    all_beziers: list[CubicBezierSegment] = []
    all_markers: list[SaddleMarkerBar] = []
    all_circles: list[FreeCirclePrimitive] = []

    cursor_x = 0.0
    spacing_gap = 0.45

    for b_vs in blocks:
        m = len(b_vs)
        g2l = {gv: lv for lv, gv in enumerate(b_vs)}
        local_opp = [0] * (4 * m)
        local_vtype = [d.vtype[gv] for gv in b_vs]
        local_scomp = [0] * (4 * m)
        for lv, gv in enumerate(b_vs):
            for s in range(4):
                od = d.opp[4 * gv + s]
                local_opp[4 * lv + s] = 4 * g2l[od >> 2] + (od & 3)
                local_scomp[4 * lv + s] = dart_scomp[4 * gv + s]

        vp, va, bzs, mks = _layout_connected_block(m, local_opp, local_vtype, local_scomp)
        pts_x = [p[0] for bz in bzs for p in (bz.p0, bz.p1, bz.p2, bz.p3)] + [p[0] for p in vp]
        b_xmin, b_xmax = min(pts_x), max(pts_x)
        shift_x = cursor_x - b_xmin

        for lv, gv in enumerate(b_vs):
            global_v_pos[gv] = (vp[lv][0] + shift_x, vp[lv][1])
            global_v_ang[gv] = va[lv]
        for bz in bzs:
            all_beziers.append(
                CubicBezierSegment(
                    p0=(bz.p0[0] + shift_x, bz.p0[1]),
                    p1=(bz.p1[0] + shift_x, bz.p1[1]),
                    p2=(bz.p2[0] + shift_x, bz.p2[1]),
                    p3=(bz.p3[0] + shift_x, bz.p3[1]),
                    surface_component=bz.surface_component,
                    layer=bz.layer,
                )
            )
        for mk in mks:
            all_markers.append(
                SaddleMarkerBar(
                    center=(mk.center[0] + shift_x, mk.center[1]),
                    a=(mk.a[0] + shift_x, mk.a[1]),
                    b=(mk.b[0] + shift_x, mk.b[1]),
                    kind=mk.kind,
                    surface_component=mk.surface_component,
                )
            )
        cursor_x += (b_xmax - b_xmin) + spacing_gap

    for i in range(d.free_loops):
        r_circ = 0.72
        cx = cursor_x + r_circ
        all_circles.append(
            FreeCirclePrimitive(
                center=(cx, 0.0),
                radius=r_circ,
                surface_component=base_graph_scomps + i,
            )
        )
        cursor_x += 2.0 * r_circ + spacing_gap

    # Compute overall bounding box
    xs_all: list[float] = []
    ys_all: list[float] = []
    for bz in all_beziers:
        for px, py in (bz.p0, bz.p1, bz.p2, bz.p3):
            xs_all.append(px)
            ys_all.append(py)
    for mk in all_markers:
        for px, py in (mk.a, mk.b):
            xs_all.append(px)
            ys_all.append(py)
    for cc in all_circles:
        xs_all.extend([cc.center[0] - cc.radius, cc.center[0] + cc.radius])
        ys_all.extend([cc.center[1] - cc.radius, cc.center[1] + cc.radius])

    if not xs_all:
        bbox = (-1.0, -1.0, 1.0, 1.0)
    else:
        bbox = (min(xs_all), min(ys_all), max(xs_all), max(ys_all))

    return DiagramLayout(
        vertex_positions=tuple(global_v_pos),
        vertex_angles=tuple(global_v_ang),
        beziers=tuple(all_beziers),
        markers=tuple(all_markers),
        circles=tuple(all_circles),
        bbox=bbox,
        num_surface_components=total_scomps,
    )


def render_svg(
    d: MarkedDiagram,
    *,
    title: str | None = None,
    subtitle: str | None = None,
    size: int = 380,
    show_card: bool = True,
) -> str:
    """Render `d` to a standalone, crisp SVG string."""
    layout = compute_diagram_layout(d)
    xmin, ymin, xmax, ymax = layout.bbox
    w_math = max(xmax - xmin, 1e-3)
    h_math = max(ymax - ymin, 1e-3)
    span = max(w_math, h_math)
    cx_math = 0.5 * (xmin + xmax)
    cy_math = 0.5 * (ymin + ymax)

    header_h = 46 if (title or subtitle) else 0
    pad = 30.0
    draw_box = float(size) - 2.0 * pad
    scale = draw_box / span
    total_w = size
    total_h = size + header_h

    def to_svg(pt: tuple[float, float]) -> tuple[float, float]:
        sx = 0.5 * total_w + (pt[0] - cx_math) * scale
        # Flip `+y` up (math) -> `+y` down (SVG) so CCW cyclic dart order stays visually CCW!
        sy = header_h + 0.5 * size - (pt[1] - cy_math) * scale
        return (sx, sy)

    stroke_w = max(2.7, min(3.8, size * 0.0088))
    bar_w = stroke_w * 1.65
    halo_w = stroke_w * 2.55
    bg_color = "#FAFAFC"

    lines: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {total_w} {total_h}" '
        f'width="{total_w}" height="{total_h}" role="img" '
        f'aria-label="{html.escape(title or "Marked surface-link diagram")}">',
    ]
    if show_card:
        lines.append(
            f'  <rect x="1" y="1" width="{total_w - 2}" height="{total_h - 2}" rx="12" '
            f'fill="{bg_color}" stroke="#E2E8F0" stroke-width="1.5"/>'
        )

    if title or subtitle:
        if title:
            lines.append(
                f'  <text x="18" y="24" font-family="ui-sans-serif, -apple-system, Segoe UI, Helvetica, Arial, sans-serif" '
                f'font-size="15" font-weight="700" fill="#0F172A">{html.escape(title)}</text>'
            )
        if subtitle:
            lines.append(
                f'  <text x="18" y="41" font-family="ui-monospace, SFMono-Regular, Menlo, Consolas, monospace" '
                f'font-size="11" fill="#475569">{html.escape(subtitle)}</text>'
            )
        # Component color legend dots in top-right when c >= 2
        if layout.num_surface_components >= 2:
            for ci in range(layout.num_surface_components):
                col = COMPONENT_PALETTE[ci % len(COMPONENT_PALETTE)]
                dot_x = total_w - 18 - (layout.num_surface_components - 1 - ci) * 14
                lines.append(f'  <circle cx="{dot_x:.1f}" cy="20" r="4.5" fill="{col}"/>')

    # 1. Free circles (`O[1]`)
    for cc in layout.circles:
        scx, scy = to_svg(cc.center)
        sr = cc.radius * scale
        col = COMPONENT_PALETTE[cc.surface_component % len(COMPONENT_PALETTE)]
        lines.append(
            f'  <circle cx="{scx:.2f}" cy="{scy:.2f}" r="{sr:.2f}" fill="none" '
            f'stroke="{col}" stroke-width="{stroke_w:.2f}"/>'
        )

    # 2. Layer 0: Under-strand Bezier segments
    for bz in layout.beziers:
        if bz.layer != 0:
            continue
        p0, p1, p2, p3 = to_svg(bz.p0), to_svg(bz.p1), to_svg(bz.p2), to_svg(bz.p3)
        col = COMPONENT_PALETTE[bz.surface_component % len(COMPONENT_PALETTE)]
        lines.append(
            f'  <path d="M {p0[0]:.2f} {p0[1]:.2f} C {p1[0]:.2f} {p1[1]:.2f}, '
            f'{p2[0]:.2f} {p2[1]:.2f}, {p3[0]:.2f} {p3[1]:.2f}" '
            f'fill="none" stroke="{col}" stroke-width="{stroke_w:.2f}" stroke-linecap="butt"/>'
        )

    # 3. Layer 1: Over-strand & marked-strand Bezier segments
    for bz in layout.beziers:
        if bz.layer != 1:
            continue
        p0, p1, p2, p3 = to_svg(bz.p0), to_svg(bz.p1), to_svg(bz.p2), to_svg(bz.p3)
        col = COMPONENT_PALETTE[bz.surface_component % len(COMPONENT_PALETTE)]
        lines.append(
            f'  <path d="M {p0[0]:.2f} {p0[1]:.2f} C {p1[0]:.2f} {p1[1]:.2f}, '
            f'{p2[0]:.2f} {p2[1]:.2f}, {p3[0]:.2f} {p3[1]:.2f}" '
            f'fill="none" stroke="{col}" stroke-width="{stroke_w:.2f}" stroke-linecap="round"/>'
        )

    # 4. Layer 2: Marked saddle vertex orthogonal bars (`Y` and `Z`)
    for mk in layout.markers:
        sa, sb = to_svg(mk.a), to_svg(mk.b)
        sc = to_svg(mk.center)
        if show_card:
            lines.append(
                f'  <line x1="{sa[0]:.2f}" y1="{sa[1]:.2f}" x2="{sb[0]:.2f}" y2="{sb[1]:.2f}" '
                f'stroke="{bg_color}" stroke-width="{halo_w:.2f}" stroke-linecap="round"/>'
            )
        lines.append(
            f'  <line x1="{sa[0]:.2f}" y1="{sa[1]:.2f}" x2="{sb[0]:.2f}" y2="{sb[1]:.2f}" '
            f'stroke="#0F172A" stroke-width="{bar_w:.2f}" stroke-linecap="round"/>'
        )
        lines.append(f'  <circle cx="{sc[0]:.2f}" cy="{sc[1]:.2f}" r="{0.72 * stroke_w:.2f}" fill="#0F172A"/>')

    lines.append("</svg>\n")
    return "\n".join(lines)


def render_tikz(
    d: MarkedDiagram,
    *,
    label: str | None = None,
    width_cm: float = 5.4,
    standalone: bool = True,
) -> str:
    """Render `d` to a publication-ready LaTeX TikZ picture (standalone `.tex` if `standalone=True`)."""
    layout = compute_diagram_layout(d)
    xmin, ymin, xmax, ymax = layout.bbox
    span = max(xmax - xmin, ymax - ymin, 1e-3)
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    scale = width_cm / span

    def tx(pt: tuple[float, float]) -> tuple[float, float]:
        return ((pt[0] - cx) * scale, (pt[1] - cy) * scale)

    body: list[str] = []
    if label:
        safe_comment = label.replace("%", "")
        body.append(f"% Canonical ch-diagram for {safe_comment}")
    body.append(r"\begin{tikzpicture}[line join=round]")

    used_comps = sorted(
        {bz.surface_component for bz in layout.beziers}
        | {cc.surface_component for cc in layout.circles}
    )
    for ci in used_comps:
        r_v, g_v, b_v = TIKZ_PALETTE_RGB[ci % len(TIKZ_PALETTE_RGB)]
        body.append(f"  \\definecolor{{skcomp{ci}}}{{RGB}}{{{r_v},{g_v},{b_v}}}")
    body.append("  \\definecolor{skmark}{RGB}{15,23,42}")

    for cc in layout.circles:
        ccx, ccy = tx(cc.center)
        cr = cc.radius * scale
        body.append(
            f"  \\draw[skcomp{cc.surface_component}, line width=1.15pt] "
            f"({ccx:.4f},{ccy:.4f}) circle ({cr:.4f});"
        )

    for layer_idx, cap_style in ((0, "butt"), (1, "round")):
        for bz in layout.beziers:
            if bz.layer != layer_idx:
                continue
            p0, p1, p2, p3 = tx(bz.p0), tx(bz.p1), tx(bz.p2), tx(bz.p3)
            body.append(
                f"  \\draw[skcomp{bz.surface_component}, line width=1.15pt, line cap={cap_style}] "
                f"({p0[0]:.4f},{p0[1]:.4f}) .. controls ({p1[0]:.4f},{p1[1]:.4f}) "
                f"and ({p2[0]:.4f},{p2[1]:.4f}) .. ({p3[0]:.4f},{p3[1]:.4f});"
            )

    for mk in layout.markers:
        pa, pb = tx(mk.a), tx(mk.b)
        pc = tx(mk.center)
        body.append(
            f"  \\draw[white, line width=2.8pt, line cap=round] "
            f"({pa[0]:.4f},{pa[1]:.4f}) -- ({pb[0]:.4f},{pb[1]:.4f});"
        )
        body.append(
            f"  \\draw[skmark, line width=1.85pt, line cap=round] "
            f"({pa[0]:.4f},{pa[1]:.4f}) -- ({pb[0]:.4f},{pb[1]:.4f});"
        )
        body.append(f"  \\fill[skmark] ({pc[0]:.4f},{pc[1]:.4f}) circle (1.2pt);")

    body.append(r"\end{tikzpicture}")
    tikz_block = "\n".join(body) + "\n"
    if not standalone:
        return tikz_block

    return (
        "\\documentclass[tikz,border=6pt]{standalone}\n"
        "\\usepackage{xcolor}\n"
        "\\begin{document}\n"
        f"{tikz_block}"
        "\\end{document}\n"
    )


def _format_grid_section(items: list[dict[str, str]], cols: int = 3) -> str:
    rows: list[str] = []
    header = "| " + " | ".join(f"Diagram {i + 1}" for i in range(cols)) + " |"
    sep = "| " + " | ".join(":---:" for _ in range(cols)) + " |"
    rows.extend([header, sep])
    for i in range(0, len(items), cols):
        chunk = items[i : i + cols]
        img_cells: list[str] = []
        cap_cells: list[str] = []
        for it in chunk:
            lbl = it["label"]
            svg_r = it["svg_rel"]
            tex_r = it["tikz_rel"]
            img_cells.append(f"![{lbl}](./{svg_r})")
            cap_cells.append(
                f"**`{lbl}`** (`{it['components']}`)<br>[SVG](./{svg_r}) · [TikZ](./{tex_r})"
            )
        while len(img_cells) < cols:
            img_cells.append("")
            cap_cells.append("")
        rows.append("| " + " | ".join(img_cells) + " |")
        rows.append("| " + " | ".join(cap_cells) + " |")
    return "\n".join(rows)


def _build_gallery_readme(manifest: list[dict[str, str]]) -> str:
    yosh = [m for m in manifest if m["group"] == "yoshikawa_ch10"]
    ch11 = [m for m in manifest if m["group"] == "expanded_ch11"]
    pkt = [m for m in manifest if m["group"] == "expanded_ch11_packet"]

    return f"""# `/pics` — Diagrammatic Gallery & Planar Rendering Specification

This directory houses automated vector illustrations (`.svg` in [`./svg/`](./svg) for web/Markdown viewing and standalone compilable `.tex` TikZ figures in [`./tikz/`](./tikz) for LaTeX manuscripts) of canonical marked graph diagrams (`ch`-diagrams) for knotted surfaces and surface-links in $\\mathbb{{R}}^4$, generated by [`drawing.py`](../code/src/surfaceknots/drawing.py) (`python scripts/reproduce_datasets.py --gallery`).

---

## 1. Visual Conventions for Marked `ch`-Diagrams

Every surface-link $F \\subset \\mathbb{{R}}^4$ in hyperbolic splitting position is represented by a 4-regular planar shadow on $S^2$ decorated at each vertex (`MarkedDiagram` in [`diagrams.py`](../code/src/surfaceknots/diagrams.py), with CCW dart indices `0, 1, 2, 3`):

1. **Classical Crossings (`X0`, `X1`)**:
   - Drawn as two smooth transverse strands with a geometric under-strand gap (`(0, 2)` for `X0`, `(1, 3)` for `X1`) and a continuous stroke on the over-passing arc (`(1, 3)` for `X0`, `(0, 2)` for `X1`).
2. **Marked Saddle Vertices (`Y`, `Z`)**:
   - Drawn as a transverse 4-valent intersection carrying a **bold orthogonal marker bar** (Yoshikawa/Lomonaco convention) bisecting the two opposite quadrants smoothed together in $L_-(D)$:
     - For `Y`: $L_-(D)$ pairs `(0, 1)` and `(2, 3)`, while $L_+(D)$ pairs `(0, 3)` and `(1, 2)`; the marker bar spans the quadrant bisector separating `(0, 3),(1, 2)` and connecting the `(0, 1)` and `(2, 3)` lobes.
     - For `Z`: $L_-(D)$ pairs `(0, 3)` and `(1, 2)`, while $L_+(D)$ pairs `(0, 1)` and `(2, 3)`; the marker bar spans the dual bisector connecting the `(0, 3)` and `(1, 2)` lobes.
3. **Component Palette**:
   - Distinct connected components $F_1, \\dots, F_c$ of $F(D)$ are rendered in high-contrast, colorblind-safe strokes (deep sapphire indigo `#1E3A8A`, vermilion `#C2410C`, viridian teal `#0F766E`, royal violet `#7E22CE`).

---

## 2. Automated Planar Layout Pipeline (`drawing.py`)

1. **Doubly-Subdivided Face-Star Triangulation (`G_aug`)**:
   - Extracts combinatorial faces directly from the CCW rotation system (`opp`), checks Euler planarity ($V - E + F = 2$), subdivides every edge twice so parallel multi-edges (`2_1^1`, `2_1^-1`) and kinks become simple cycles, and inserts radial face-star vertices inside bounded faces.
2. **IRLS $p$-Laplacian Tutte Harmonic Embedding + Polar Lobe Puffing**:
   - Solves $10$ iteratively reweighted $p$-Laplacian Tutte linear systems ($p \\approx 3.35$) with strictly positive spring weights, guaranteeing convex interior faces and zero crossings while equalizing interior face areas, then puffs outer-face subdivision control points into $r > 1$ for rounded knot-table lobes.
3. **Convex-Corridor $C^1$ Cubic Bézier Splines & Dual Targets**:
   - Routes every half-edge inside its disjoint Tutte corridor with $C^1$-collinear tangents at edge midpoints and crossing centers, exporting both `./svg/<stem>.svg` and standalone compilable `./tikz/<stem>.tex`.

---

## 3. Yoshikawa's Canonical Prime Table ($\\operatorname{{ch}}(F) \\le 10$; 23 Entries)

{_format_grid_section(yosh, cols=3)}

---

## 4. Expanded Certified Prime Table ($\\operatorname{{ch}}(F) = 11$; 10 Entries)

{_format_grid_section(ch11, cols=3)}

---

## 5. Move-Component Packet Representatives (`11_1.a`, `11_1.b`, `11_1.c`)

{_format_grid_section(pkt, cols=3)}
"""


def export_canonical_gallery(
    data_dir: Path,
    pics_dir: Path,
    *,
    update_readme: bool = True,
) -> list[dict[str, str]]:
    """Generate `.svg` and standalone `.tex` TikZ files in `pics/svg/` and `pics/tikz/`
    for all 23 Yoshikawa `ch <= 10` entries and all 10 certified prime `ch = 11` entries
    (plus packet move-components `11_1.a, 11_1.b, 11_1.c`).
    """
    svg_dir = pics_dir / "svg"
    tikz_dir = pics_dir / "tikz"
    svg_dir.mkdir(parents=True, exist_ok=True)
    tikz_dir.mkdir(parents=True, exist_ok=True)

    yosh_entries = json.loads((data_dir / "yoshikawa_ch10.json").read_text())
    ch11_payload = json.loads((data_dir / "expanded_ch11.json").read_text())
    ch11_entries = ch11_payload["entries"]

    manifest: list[dict[str, str]] = []

    def _emit(
        lbl: str,
        epd: str,
        comp_str: str,
        ch_val: int,
        cell_status: str,
        group_tag: str,
    ) -> None:
        d = MarkedDiagram.from_epd(epd)
        stem = label_to_stem(lbl)
        sub = f"ch={ch_val} | {comp_str} | {cell_status}"
        svg_txt = render_svg(d, title=lbl, subtitle=sub)
        tex_txt = render_tikz(d, label=lbl, standalone=True)
        svg_path = svg_dir / f"{stem}.svg"
        tex_path = tikz_dir / f"{stem}.tex"
        svg_path.write_text(svg_txt)
        tex_path.write_text(tex_txt)
        manifest.append(
            {
                "label": lbl,
                "stem": stem,
                "group": group_tag,
                "components": comp_str,
                "ch_index": str(ch_val),
                "cell_status": cell_status,
                "svg_rel": f"svg/{stem}.svg",
                "tikz_rel": f"tikz/{stem}.tex",
            }
        )

    for rec in yosh_entries:
        comp_str = " u ".join(rec["component_types"])
        _emit(
            rec["label"],
            rec["epd"],
            comp_str,
            int(rec["ch_index_diagram"]),
            str(rec["cell_status"]),
            "yoshikawa_ch10",
        )

    for rec in ch11_entries:
        comp_str = " u ".join(rec["component_types"])
        _emit(
            rec["label"],
            rec["epd"],
            comp_str,
            int(rec["ch_index_diagram"]),
            str(rec["cell_status"]),
            "expanded_ch11",
        )
        if int(rec["num_move_components"]) > 1:
            for sub_id, sub_epd in zip(rec["move_component_ids"], rec["move_component_epds"]):
                _emit(
                    sub_id,
                    sub_epd,
                    comp_str,
                    int(rec["ch_index_diagram"]),
                    "NEW_PACKET_MEMBER",
                    "expanded_ch11_packet",
                )

    if update_readme:
        (pics_dir / "README.md").write_text(_build_gallery_readme(manifest))

    return manifest
