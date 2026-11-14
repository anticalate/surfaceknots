# `/pics` — Diagrammatic Gallery & Planar Rendering Specification

This directory houses vector illustrations (`.svg` for web/markdown viewing and standalone `.tex` TikZ figures for LaTeX manuscripts) of canonical marked graph diagrams (`ch`-diagrams) for knotted surfaces and surface-links in $\mathbb{R}^4$.

---

## 1. Visual Conventions for Marked `ch`-Diagrams

Every surface-link $F \subset \mathbb{R}^4$ in hyperbolic splitting position is represented by a 4-regular planar shadow on $S^2$ decorated at each vertex (`MarkedDiagram` in [`diagrams.py`](../code/src/surfaceknots/diagrams.py), with CCW dart indices `0, 1, 2, 3`):

1. **Classical Crossings (`X0`, `X1`)**:
   - Drawn as two smooth transverse strands with a clean **white halo / under-strand gap** on the under-passing arc (`(0, 2)` for `X0`, `(1, 3)` for `X1`) and a continuous stroke on the over-passing arc (`(1, 3)` for `X0`, `(0, 2)` for `X1`).
2. **Marked Saddle Vertices (`Y`, `Z`)**:
   - Drawn as a transverse 4-valent intersection carrying a **bold orthogonal marker bar** (Yoshikawa/Lomonaco convention) bisecting the two opposite quadrants that get smoothed together in $L_-(D)$:
     - For `Y`: $L_-(D)$ pairs `(0, 1)` and `(2, 3)`, while $L_+(D)$ pairs `(0, 3)` and `(1, 2)`; the marker bar spans the axis separating `(0, 1)` from `(2, 3)` (i.e., connecting the `(0, 1)` and `(2, 3)` lobes).
     - For `Z`: $L_-(D)$ pairs `(0, 3)` and `(1, 2)`, while $L_+(D)$ pairs `(0, 1)` and `(2, 3)`; the marker bar spans the dual axis connecting the `(0, 3)` and `(1, 2)` lobes.
3. **Component Palette**:
   - Distinct connected components $F_1, \dots, F_c$ of $F(D)$ are rendered in high-contrast, colorblind-safe strokes (deep indigo, vermilion, teal, amber) so multi-component links ($S^2 \sqcup T^2$, $S^2 \sqcup K^2$, $T^2 \sqcup T^2$, $\mathbb{R}P^2 \sqcup \mathbb{R}P^2$, $S^2 \sqcup S^2 \sqcup T^2$) are legible at a glance.

---

## 2. Automated Planar Layout Pipeline (`drawing.py`)

Because EPD codes specify purely combinatorial rotation systems (`opp` involution on CCW darts `0..4n-1` plus `free_loops`), generating aesthetic knot-table figures without hand-tuning coordinates requires a two-stage geometric solver:

1. **Outer-Face Selection & Subdivided Planar Dual/Primal Graph**:
   - Extract the planar faces of the connected components of $D$ directly from the combinatorial rotation system.
   - Subdivide every graph edge (including parallel multi-edges and self-loops/kinks!) with midpoint routing nodes so the subdivided graph is **simple and 3-connected** (or biconnected with virtual outer rim), preventing bigons (`2`-cycles) and R1 kinks (`1`-cycles) from collapsing onto collinear segments.
   - Choose the outer face maximizing `(face_degree, symmetry_score)` and place its boundary vertices evenly around the unit circle.
2. **Tutte Barycentric Initialization + Angular/Repulsive Spring Relaxation**:
   - Solve Tutte's linear harmonic system so all interior vertices and edge midpoints are guaranteed **strictly crossing-free** inside the unit disk.
   - Refine interior positions with damped spring + angular-orthogonality relaxation (encouraging $90^\circ$ separation between adjacent darts `0, 1, 2, 3` at every 4-valent vertex while preserving positive face areas so no spurious edge crossings ever appear).
3. **Smooth Cubic Bézier Splines & Dual Output Targets**:
   - Route each half-edge from vertex $v$ along its exact orthogonal slot ray $\theta_v + s \cdot \frac{\pi}{2}$ to the edge midpoint with $C^1$-continuous cubic Bézier tangents.
   - Export both `pics/svg/<id>.svg` and standalone `pics/tikz/<id>.tex` for all canonical `ch <= 10` (Yoshikawa's 23 entries) and `ch = 11` (expanded prime entries) diagrams.
