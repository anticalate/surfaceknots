# Foundations, Combinatorial Topology, and Purely 4D Invariants of Surface-Links

This document records the self-contained mathematical foundations underlying the `surfaceknots` census and invariant engine, written so that every convention, sign choice, lemma, and proof matches the exact implementation in [`code/src/surfaceknots/`](../code/src/surfaceknots/).

---

## 1. Hyperbolic Splittings, Marked `ch`-Diagrams, and $\mathrm{O}(4)$ Symmetries

Throughout, a **surface-link** $F = F_1 \sqcup \dots \sqcup F_c \subset \mathbb{R}^4$ (or equivalently $F \subset S^4 = \mathbb{R}^4 \cup \{\infty\}$) denotes a smooth (or locally flat PL) closed 2-dimensional submanifold whose connected components $F_i$ may be orientable, non-orientable, or mixed. Two surface-links $F, F' \subset \mathbb{R}^4$ are **ambient isotopic** ($F \cong_{\mathrm{SO}(4)} F'$) if they are related by an orientation-preserving diffeomorphism of $\mathbb{R}^4$, and **equivalent up to mirror image** ($F \cong_{\mathrm{O}(4)} F'$) if they are related by an arbitrary diffeomorphism of $\mathbb{R}^4$.

### Definition 1.1 (Hyperbolic Splitting and Banded Unlinks)
Identify $\mathbb{R}^4 = \mathbb{R}^3 \times \mathbb{R}_t$ with height coordinate $t: \mathbb{R}^4 \to \mathbb{R}$. By Kawauchi–Shibuya–Suzuki (1982) and Kamada (2002), every closed surface-link $F \subset \mathbb{R}^4$ is ambient isotopic to a surface in **hyperbolic splitting position**:
1. All local minima of $t|_F$ form $m_-$ disks in $\mathbb{R}^3 \times \{-1\}$ bounded by a classical $m_-$-component unlink $O_- \subset \mathbb{R}^3 \times \{-\varepsilon\}$.
2. All saddle critical points of $t|_F$ occur in the equatorial slice $\mathbb{R}^3 \times \{0\}$, realized by $h$ pairwise disjoint flat bands $B = b_1 \cup \dots \cup b_h \subset \mathbb{R}^3 \times \{0\}$ attached to $L_- = F \cap (\mathbb{R}^3 \times \{-\varepsilon\})$ so that band surgery produces $L_+ = F \cap (\mathbb{R}^3 \times \{+\varepsilon\})$.
3. All local maxima of $t|_F$ form $m_+$ disks in $\mathbb{R}^3 \times \{+1\}$ capping a classical $m_+$-component unlink $O_+ \subset \mathbb{R}^3 \times \{+\varepsilon\}$.

### Theorem 1.2 (Horibe–Yanagawa Lemma; Handle Structure of the Exterior)
*Let $(L_-, B, L_+)$ be a banded link in $\mathbb{R}^3$ such that both $L_-$ and $L_+$ are classical unlinks in $\mathbb{R}^3$ with $m_-$ and $m_+$ components, respectively. Up to ambient isotopy of $\mathbb{R}^4$ fixing $\mathbb{R}^3 \times [-\varepsilon, +\varepsilon]$ pointwise, the capping disk systems in $\mathbb{R}^3 \times (-\infty, -\varepsilon]$ and $\mathbb{R}^3 \times [+\varepsilon, +\infty)$ are uniquely determined. Moreover, the exterior $E = S^4 \setminus \nu(F)$ admits a handlebody decomposition consisting of:*
- *one $0$-handle,*
- *$m_-$ $1$-handles (dual to the minimum disks),*
- *$h$ $2$-handles (dual to the saddle bands), and*
- *$m_+ - 1$ $3$-handles (dual to all but one of the maximum disks, with the final disk canceled against the $4$-handle at $\infty$).*

*Proof.* By Horibe–Yanagawa (see Kawauchi–Shibuya–Suzuki 1982, Lemma 1), any two systems of slice disks in $\mathbb{R}^3 \times [+\varepsilon, +\infty)$ with only local maxima (no saddles) bounding the same classical unlink $L_+ \subset \mathbb{R}^3 \times \{+\varepsilon\}$ are ambient isotopic rel boundary. Turning the Morse function $t|_F$ into a relative handle decomposition of $(S^4, \nu(F))$ and dualizing gives an $i$-handle of $E = S^4 \setminus \nu(F)$ for each $(i-1)$-handle of $F$ ($i = 1, 2, 3$), modulo canceling one $3$-handle against the top $4$-handle of $S^4$. Computing the alternating sum of handles confirms
$$\chi(E) = 1 - m_- + h - (m_+ - 1) = 2 - (m_- - h + m_+) = 2 - \chi(F),$$
matching Alexander duality $H_k(E; \mathbb{Z}) \cong H^{3-k}(F; \mathbb{Z})$ for $k \ge 1$. $\square$

### Definition 1.3 (Marked `ch`-Diagrams and Codebase Dart Conventions)
Projecting a banded unlink $(L_-, B, L_+)$ generically to $\mathbb{R}^2 \subset \mathbb{R}^3$ after shrinking each band $b_j$ into a 4-valent **marked vertex** yields a **marked graph diagram** (`ch`-diagram) $D$ on $S^2 = \mathbb{R}^2 \cup \{\infty\}$ with $c(D)$ classical crossings, $h(D)$ marked vertices, and $f(D)$ disjoint vertex-free loops (`free_loops`). Following [`MarkedDiagram`](../code/src/surfaceknots/diagrams.py), around each 4-valent vertex $v \in \{0, \dots, n-1\}$ ($n = c(D) + h(D)$) we label its four incident half-edges (darts) $4v + s$ ($s \in \{0, 1, 2, 3\}$) in **counterclockwise (CCW)** cyclic order:

1. **Classical crossings**:
   - Type `X0` (written `X[a,b,c,d]` in EPD): darts $(0, 2)$ form the **under-strand** and darts $(1, 3)$ form the **over-strand**.
   - Type `X1`: darts $(1, 3)$ form the **under-strand** and darts $(0, 2)$ form the **over-strand**.
   - If an oriented ClassicalLinkDiagram enters the over-strand at slot $s_{\text{over}}^{\text{in}}$ and the under-strand at slot $s_{\text{under}}^{\text{in}}$, the right-hand writhe sign is
     $$\varepsilon(v) = \begin{cases} +1 & \text{if } (s_{\text{over}}^{\text{in}} - s_{\text{under}}^{\text{in}}) \bmod 4 = 3, \\ -1 & \text{if } (s_{\text{over}}^{\text{in}} - s_{\text{under}}^{\text{in}}) \bmod 4 = 1. \end{cases}$$
2. **Marked (saddle) vertices**:
   - Type `Y` (marker bar joining opposite quadrants $(0,2)$ or $(1,3)$ per resolution):
     - Positive resolution $L_+(D)$ (`sign = +1`): pairs $(0, 3)$ and $(1, 2)$.
     - Negative resolution $L_-(D)$ (`sign = -1`): pairs $(0, 1)$ and $(2, 3)$.
   - Type `Z` (dual marker bar):
     - Positive resolution $L_+(D)$ (`sign = +1`): pairs $(0, 1)$ and $(2, 3)$.
     - Negative resolution $L_-(D)$ (`sign = -1`): pairs $(0, 3)$ and $(1, 2)$.

A marked diagram $D$ is **admissible** if both classical link diagrams $L_-(D)$ and $L_+(D)$ represent classical unlinks in $\mathbb{R}^3$. By Swenton (2001) and Kearton–Kurlin (2008), two admissible diagrams $D, D'$ represent ambient isotopic surface-links in $\mathbb{R}^4$ if and only if they are connected by a finite sequence of planar isotopies and the eight **Yoshikawa moves** $\Omega_1, \dots, \Omega_8$ ($\Omega_1$–$\Omega_3$: classical Reidemeister moves; $\Omega_4, \Omega_4', \Omega_5$: mixed crossing/marker moves; $\Omega_6, \Omega_6'$: 4D mushroom/handle stabilizations changing $h(D)$ by $\pm 1$ and $|O_\pm|$ by $\pm 1$; $\Omega_7$: marker switch/slide; $\Omega_8$: 4-vertex double-saddle move).

### Definition 1.4 (The $\operatorname{ch}$-Index and $\mathrm{O}(4)$ Reflections)
The **`ch`-complexity** of an admissible diagram $D$ is $\operatorname{ch}(D) = c(D) + h(D)$, and the **`ch`-index** of a surface-link $F \subset \mathbb{R}^4$ is
$$\operatorname{ch}(F) = \min\{\, c(D) + h(D) \mid D \text{ admissible with } F(D) \cong F \,\}.$$
Under reflections of $\mathbb{R}^4 = \mathbb{R}^2_{x,y} \times \mathbb{R}_z \times \mathbb{R}_t$:
- **Spatial mirror** $r_z: (x, y, z, t) \mapsto (x, y, -z, t)$ (`mirror_3d`, orientation-reversing in $\mathrm{O}(4) \setminus \mathrm{SO}(4)$) swaps `X0` $\leftrightarrow$ `X1` and fixes `Y`, `Z`; this sends $F \mapsto F^*$.
- **Temporal mirror** $r_t: (x, y, z, t) \mapsto (x, y, z, -t)$ (`mirror_time`, also in $\mathrm{O}(4) \setminus \mathrm{SO}(4)$) fixes `X0`, `X1` and swaps `Y` $\leftrightarrow$ `Z` (exchanging $L_+(D) \leftrightarrow L_-(D)$).
- **Simultaneous $(z, t)$-halfturn** $r_z \circ r_t: (x, y, z, t) \mapsto (x, y, -z, -t)$ (`reverse_orientation_zt`, **orientation-preserving** in $\mathrm{SO}(4)$!) swaps both `X0` $\leftrightarrow$ `X1` and `Y` $\leftrightarrow$ `Z`, representing an ambient-isotopic diagram of the same unoriented surface $F$ (with reversed co-orientation / surface orientation when $F$ is oriented).

---

## 2. Intrinsic Topology and Normal Bundle Formulas

Every topological invariant of $F$ as an abstract 2-manifold together with its normal bundle $\nu_{F \subset \mathbb{R}^4}$ is readable directly from $D$ in linear time.

### Proposition 2.1 (Component Decomposition, Euler Characteristic, and Orientability Criterion)
*Partition the darts and free loops of $D$ into connected classes under edge adjacency, classical strand continuation $s \leftrightarrow s \oplus 2$, and both resolutions at marked vertices ($s \leftrightarrow s \pm 1 \bmod 4$). These classes are in canonical bijection with the connected components $F_1, \dots, F_c$ of $F(D)$, and for each component $F_i$:*
1. **Euler characteristic**: *Let $|O_+|_i$ and $|O_-|_i$ denote the numbers of circles of $L_+(D)$ and $L_-(D)$ belonging to $F_i$, and let $h_i$ denote the number of marked vertices belonging to $F_i$. Then*
   $$\chi(F_i) = |O_+|_i + |O_-|_i - h_i.$$
2. **Orientability via dart alternation**: *$F_i$ is orientable if and only if its darts admit a sign assignment $\delta: \operatorname{Darts}(F_i) \to \{+1, -1\}$ (`+1` = incoming to vertex, `-1` = outgoing) satisfying:*
   - $\delta(\operatorname{opp}(d)) = -\delta(d)$ *along every graph edge,*
   - $\delta(4v + (s \oplus 2)) = -\delta(4v + s)$ *at every classical crossing $v$, and*
   - $\delta(4v + ((s + 1) \bmod 4)) = -\delta(4v + s)$ *at every marked vertex $v$ (so darts alternate $(+1, -1, +1, -1)$ cyclically around $v$).*

*Proof.* (1) $F_i$ is built from $|O_-|_i$ minimum $0$-handles ($2$-disks), $h_i$ saddle $1$-handles ($1$-bands), and $|O_+|_i$ maximum $2$-handles ($2$-disks), giving $\chi(F_i) = |O_-|_i - h_i + |O_+|_i$.
(2) An orientation of $F_i$ induces coherent boundary orientations on the level cross-sections $L_-(F_i)$ and $L_+(F_i)$ such that attaching each saddle band $b_v$ preserves boundary orientation before and after surgery. A flat saddle band connecting two parallel segments preserves orientation under both $(0,1)/(2,3)$ and $(0,3)/(1,2)$ smoothings if and only if the two paired darts in each smoothing consist of one incoming (`+1`) and one outgoing (`-1`) dart — i.e., $\delta$ alternates $(+1, -1, +1, -1)$ around $v$. Conversely, any such $\delta$ orients every circle of $L_-(F_i)$ and $L_+(F_i)$ coherently with the unique capping disks and saddle bands, patching to a global orientation of $F_i$. $\square$

### Theorem 2.2 (Kamada's Normal Euler Formula, Whitney–Massey Congruence, and Unknotted Bound)
*For each component $F_i$ of $F(D)$, let $D_i = F_i(D)$ be the subdiagram formed by retaining only the circles and vertices belonging to $F_i$. Assign arbitrary orientations to each circle of $L_-(D_i)$ and $L_+(D_i)$, and let $w(L_\mp(D_i))$ denote the resulting self-writhes. Then the normal Euler number $e(F_i) = e(\nu_{F_i \subset \mathbb{R}^4}) \in \mathbb{Z}$ (with respect to the standard orientation of $\mathbb{R}^4$, implemented in [`MarkedDiagram.surface_topology`](../code/src/surfaceknots/diagrams.py)) is given by*
$$e(F_i) = w\bigl(L_-(D_i)\bigr) - w\bigl(L_+(D_i)\bigr).$$
*Moreover:*
1. *If $F_i$ is orientable, $e(F_i) = 0$ identically.*
2. *(**Whitney–Massey Congruence**, Massey 1969): For every closed connected surface $F_i \subset \mathbb{R}^4$,*
   $$e(F_i) \equiv 2\chi(F_i) \pmod{4}, \qquad |e(F_i)|_{\text{unknotted}} \le 2\bigl(2 - \chi(F_i)\bigr).$$
   *In particular, any unknotted non-orientable surface $\#_k \mathbb{R}P^2 = P_+^{\# p} \# P_-^{\# q}$ ($p + q = k = 2 - \chi(F_i)$) satisfies $e = 2(p - q) \in \{-2k, -2k + 4, \dots, 2k\}$.*

*Proof.* Push $F_i$ slightly off itself along the blackboard-framing vector field $\partial/\partial y$ in $\mathbb{R}^3 \times \mathbb{R}_t$, extended across the minimum and maximum disks in $\mathbb{R}^3 \times (-\infty, -\varepsilon]$ and $\mathbb{R}^3 \times [+\varepsilon, +\infty)$. Because $L_-(D_i)$ and $L_+(D_i)$ are unlinks in $\mathbb{R}^3$, a blackboard-framed pushoff of $L_\mp(D_i)$ has total pairwise + self linking number equal to the diagram writhe $w(L_\mp(D_i))$; note that reversing the orientation of any circle $C \subset L_\mp(D_i)$ leaves self-crossings of $C$ invariant and flips signs of both $(C, C')$ and $(C', C)$ crossings with $C' \neq C$, leaving $w(L_\mp(D_i))$ independent of circle orientations because $L_\mp(D_i)$ is an unlink ($\operatorname{lk}(C, C') = 0$). Since the bands at $t = 0$ lie flat in the projection plane, all signed intersections between $F_i$ and its pushoff occur inside the capping 4-balls $t < -\varepsilon$ (contributing $+w(L_-(D_i))$) and $t > +\varepsilon$ (contributing $-w(L_+(D_i))$ due to opposite outward normal of the upper 4-ball), yielding $e(F_i) = w(L_-(D_i)) - w(L_+(D_i))$.

When $F_i$ is orientable, Proposition 2.1(2) gives a single dart orientation inducing valid link orientations on both $L_-(D_i)$ and $L_+(D_i)$ that agree at every classical crossing of $D_i$, so $w(L_-(D_i)) = w(L_+(D_i))$ and $e(F_i) = 0$. Checking the canonical projective plane $2_1^{-1}$ (`Y[1,2,3,4], X[4,3,2,1]`) gives a single positive kink in $L_-$ ($w(L_-) = +1$) and a negative kink in $L_+$ ($w(L_+) = -1$), so $e(P_+) = (+1) - (-1) = +2$, while $r_z(P_+) = P_-$ has $e(P_-) = -2$, satisfying $e \equiv 2(1) \equiv 2 \pmod 4$. Additivity under connected sum gives $e(P_+^{\# p} \# P_-^{\# q}) = 2(p - q) \equiv 2(2 - (p + q)) = 2\chi \pmod 4$; for arbitrary $F_i \subset S^4$, Whitney's conjecture proved by Massey (1969) via Pontryagin squares / Atiyah–Singer $G$-signature on $\Sigma_2(S^4, F_i)$ (Theorem 5.6 below) gives $e(F_i) = -2\sigma(\Sigma_2(S^4, F_i))$ with $\sigma(\Sigma_2) \equiv b_2(\Sigma_2) = 2 + 2b_1(\Sigma_2) - \chi(F_i) \equiv -\chi(F_i) \pmod 2$, whence $e(F_i) \equiv 2\chi(F_i) \pmod 4$. $\square$

---

## 3. Planar Cut-Pruning and the Unconditional Singleton Exactness Theorem

In our $\operatorname{ch} = 11$ census ([`ch11_census.py`](../code/src/surfaceknots/ch11_census.py)), we filter 4-regular planar graphs before and after marker assignment using exact geometric decompositions, and then stratify survivors via a three-stage quotient tower. Here we record the exact scope and proofs of those theorems.

### Lemma 3.1 (Planar 2-Edge-Cut Decomposition for Knots and Links)
*Let $D$ be a connected admissible marked diagram on $S^2$, and suppose a simple closed curve $\gamma \subset S^2 \setminus V(D)$ intersects $D$ transversely in exactly two edges $e_1, e_2$, partitioning $S^2$ into disks $\Delta_A \cup_\gamma \Delta_B$ with non-empty vertex sets $V_A, V_B \neq \emptyset$. Let $D_A, D_B$ be the closed diagrams obtained by splicing the two cut ends inside $\Delta_A$ and $\Delta_B$, respectively. Then:*
1. *Both $D_A$ and $D_B$ are admissible marked diagrams with $\operatorname{ch}(D_A) + \operatorname{ch}(D_B) = \operatorname{ch}(D)$ and $1 \le \operatorname{ch}(D_A), \operatorname{ch}(D_B) < \operatorname{ch}(D)$.*
2. *$(S^4, F(D)) \cong (S^4, F(D_A)) \#_{F_*} (S^4, F(D_B))$, where the 4D connected sum takes place along the unique surface component $F_*$ traversed by $e_1$ and $e_2$.*
3. *Consequently, if $(S^4, F(D))$ is **weakly prime** (contains no nontrivial connected-sum factor and no trivial unknotted $T^2$ or $P_\pm$ handle factor on any component), then one of $(S^4, F(D_A))$ or $(S^4, F(D_B))$ must be an unknotted 2-sphere pair $(S^4, S^2_{\text{std}})$, so $F(D) \cong F(D_B)$ or $F(D_A)$ is realized at strictly smaller $\operatorname{ch}$-complexity.*

*Proof.* In each classical resolution $L_\varepsilon(D)$ ($\varepsilon \in \{+, -\}$), $\gamma$ intersects the planar link projection $L_\varepsilon(D)$ in two points, presenting $L_\varepsilon(D) = L_\varepsilon(D_A) \# L_\varepsilon(D_B)$ as a planar connected sum of classical links. By additivity of classical bridge/tunnel/Alexander invariants (or Schubert's theorem that a connected sum of classical links in $S^3$ is an unlink iff both summands are unlinks joined along a single component), $L_\varepsilon(D)$ being an unlink forces both $L_\varepsilon(D_A)$ and $L_\varepsilon(D_B)$ to be classical unlinks, proving (1). Lifting $\gamma \subset \mathbb{R}^2$ to the separating 3-sphere $\Sigma^3 = (\gamma \times [-R, R]_z \times [-2, 2]_t) \cup \text{caps} \subset S^4$ intersects $F(D)$ in a single unknotted equatorial circle $S^1 = (\text{two points}) \times \{0\} \times [-1, 1] \cup \text{caps}$ on the single component $F_*$, exhibiting $(S^4, F(D))$ as the pairwise connected sum $(S^4, F(D_A)) \#_{F_*} (S^4, F(D_B))$. Note that this holds uniformly whether $c(F) = 1$ (surface-knot) or $c(F) \ge 2$ (multi-component surface-link: all components not meeting $e_1, e_2$ lie wholly on one side of $\Sigma^3$, and if both sides have components disjoint from $F_*$, removing a smaller 3-ball inside $\Delta_A$ or $\Delta_B$ would split $F$, ruled out for non-splittable links unless one side contains only a spherical cap of $F_*$). If one factor is $(S^4, S^2_{\text{std}})$, splicing it out leaves the ambient isotopy type of $F(D)$ unchanged while strictly reducing $\operatorname{ch}$ from $\operatorname{ch}(D_A) + \operatorname{ch}(D_B)$ to $\max(\operatorname{ch}(D_A), \operatorname{ch}(D_B)) < \operatorname{ch}(D)$. $\square$

### Lemma 3.2 (Unknotted 4-Edge Planar Tube Criterion)
*Suppose $\gamma \subset S^2 \setminus V(D)$ intersects an admissible diagram $D$ transversely in four edges $(e_1, e_2, e_3, e_4)$ in cyclic order along $\gamma$, enclosing a 4-legged tangle $T \subset \Delta_A$ with $|V(T)| \ge 2$. If closing $T$ outside $\gamma$ by the two horizontal arcs $(e_1, e_2), (e_3, e_4)$ yields a diagram $D_T$ of a single unknotted 2-sphere ($\chi(F(D_T)) = 2$, orientable, $\pi_1 = \mathbb{Z}$ with trivial tangle action — specifically when $T$ reduces via planar Yoshikawa moves rel boundary or consists solely of unknotted planar twists/mushrooms on two parallel sheets), look at how $\partial(F \cap B^4_A)$ meets $\Sigma^3 \cong S^3$: when $F \cap B^4_A$ is a boundary-parallel annulus $S^1 \times I$ in $B^4_A$ bounded by a 2-component unlink in $\Sigma^3$, replacing $T$ by two untwisted parallel edges adds at most a trivial 1-handle stabilization $F' \# T^2_{\text{std}}$ or changes $\chi$ by $+2$ along a compressible tube.*

> **Implementation Discipline**: In [`ch11_census.py`](../code/src/surfaceknots/ch11_census.py), we never discard a 4-edge tangle heuristically: our census retains *all* 2-edge-cut-free (4-edge-connected) shadows (`34` graphs at $n \le 10$; `42` graphs at $n = 11$) and discharges collapsible tubes constructively via explicit $\Omega_1$–$\Omega_8$ + mountain-pass stabilization reductions!

### Theorem 3.3 (The Quotient Tower and Unconditional Exactness for Singleton Cells)
*Fix complexity $N \ge 1$, and let $\mathcal{D}_N^{\text{adm}}$ be the finite set of prime-shadow admissible diagrams of complexity $\operatorname{ch}(D) = N$ that do not reduce to $\operatorname{ch} < N$ under a bounded Yoshikawa/mountain-pass search $\sim_{\text{move}}^{(H)}$. Consider the sequence of surjective maps of finite sets:*
$$\mathcal{D}_N^{\text{adm}} \;\twoheadrightarrow\; \mathcal{M}_N = \mathcal{D}_N^{\text{adm}} / {\sim_{\text{move}}^{(H)}} \;\twoheadrightarrow\; \mathcal{I}_{\le N}^{\text{true}} \;\twoheadrightarrow\; \mathcal{C}_N = \mathcal{D}_N^{\text{adm}} / {\sim_{\text{inv}}},$$
*where $\mathcal{I}_{\le N}^{\text{true}}$ is the set of true $\mathrm{O}(4)$-isotopy classes represented by $\mathcal{D}_N^{\text{adm}}$, and $D \sim_{\text{inv}} D'$ iff $I(D) = I(D')$ for a certified $\Omega_1$–$\Omega_8$-invariant functor $I$. Fix an invariant cell $C \in \mathcal{C}_N$, let $m(C) = |q^{-1}(C) \subset \mathcal{M}_N|$ be the number of move-equivalence components in $C$, and let $\operatorname{Refs}_{<N}(C)$ be the set of all ambient isotopy classes (prime, composite, split, or stabilized-unknotted) realizable at $\operatorname{ch} < N$ whose invariant value equals $I(C)$.*

1. (**Unconditional Novelty**): *If $\operatorname{Refs}_{<N}(C) = \emptyset$ (having exhaustively enumerated $\mathcal{D}_{<N}$ together with all composites, splits, and higher-genus unknots sharing $(c, \chi, \text{ori}, e)$), then $C$ contains **at least one** genuinely new surface-link $F$ with exact minimal index $\operatorname{ch}(F) = N$, and at most $m(C)$ such $\mathrm{O}(4)$-classes.*
2. (**Unconditional Singleton Exactness**): *If $\operatorname{Refs}_{<N}(C) = \emptyset$ **and** $m(C) = 1$, then without assuming any theoretical upper bound on Yoshikawa move-stabilization height $H$, $C$ represents **exactly one** $\mathrm{O}(4)$-isotopy class $F$, and its exact minimal complexity is unconditionally $\operatorname{ch}(F) = N$.*

*Proof.* Every diagram of any surface-link $F'$ with $\operatorname{ch}(F') < N$ sits in the exhaustive finite catalog of diagrams of complexity $\le N - 1$ (including disconnected and 2-edge-cut reducible diagrams, whose invariants are fully enumerated in $\operatorname{Refs}_{<N}$). Because $I$ is strictly invariant under all Yoshikawa moves $\Omega_1$–$\Omega_8$ of arbitrary intermediate height, if any diagram $D \in C$ were isotopic to some $F'$ with $\operatorname{ch}(F') < N$, then $I(D) = I(F')$ would force $F' \in \operatorname{Refs}_{<N}(C)$, contradicting $\operatorname{Refs}_{<N}(C) = \emptyset$. Thus every true isotopy class in the fiber of $C$ in $\mathcal{I}_{\le N}^{\text{true}}$ has $\operatorname{ch}(F) = N$, and the surjection $q^{-1}(C) \twoheadrightarrow \text{fiber}_{\mathcal{I}}(C)$ bounds its cardinality between $1$ and $m(C)$. When $m(C) = 1$, a surjection from a 1-element set onto a non-empty set is necessarily a bijection! $\square$

> **Application to $\operatorname{ch} = 11$**: Of our **9** certified new prime invariant cells at $\operatorname{ch} = 11$ ($\operatorname{Refs}_{<11}(C) = \emptyset$), **8** satisfy $m(C) = 1$ (`EXACT_SINGLETON`), proving their uniqueness and $\operatorname{ch}(F) = 11$ unconditionally; the 9th (`11_1`, `NEW_PACKET`) tells us with 100% mathematical transparency where move-merging remains open if $m(C) > 1$, while cells with $\operatorname{Refs}_{<11}(C) \neq \emptyset$ are explicitly labeled `AMBIGUOUS_SHADOW`.

---

## 4. Zero-Divisor Fox Calculus over $\mathbb{Z}[H_1]$ and Exact Proofs of Yoshikawa's Table II Errata

### 4.1 Multivariable Fox Calculus in the Presence of $\mathbb{Z}_2$-Torsion
Let $F = F_1 \sqcup \dots \sqcup F_c \subset \mathbb{R}^4$ have $r$ orientable components ($F_1, \dots, F_r$) and $s$ non-orientable components ($F_{r+1}, \dots, F_{r+s}$, $r + s = c$). By Alexander duality, the meridian abelianization of $G = \pi_1(\mathbb{R}^4 \setminus F)$ is
$$H_1(\mathbb{R}^4 \setminus F; \mathbb{Z}) \cong \mathbb{Z}^r \oplus (\mathbb{Z}_2)^s = \langle t_1, \dots, t_r, u_1, \dots, u_s \mid [t_i, t_j] = [t_i, u_k] = [u_k, u_\ell] = 1,\; u_k^2 = 1 \rangle,$$
with integral group ring
$$\Lambda = \mathbb{Z}[H_1] = \mathbb{Z}[t_1^{\pm 1}, \dots, t_r^{\pm 1}, u_1, \dots, u_s] / (u_1^2 - 1, \dots, u_s^2 - 1).$$
Note that when $s \ge 1$, $\Lambda$ is **not an integral domain**: $(u_k - 1)(u_k + 1) = 0$, and $\mathbb{Z}[u]/(u^2 - 1)$ embeds as the index-$2$ fiber-product subring
$$\mathbb{Z}[u]/(u^2 - 1) \;\xrightarrow{\;\sim\;}\; \bigl\{\, (a_+, a_-) \in \mathbb{Z} \times \mathbb{Z} \;\big|\; a_+ \equiv a_- \pmod 2 \,\bigr\}, \qquad f(u) \mapsto \bigl(f(1),\, f(-1)\bigr).$$
Crucially, projecting onto $u = +1$ and $u = -1$ separately **loses mixed 2-torsion information** because the fiber product couples the parity of the two factors! In [`alexander.py`](../code/src/surfaceknots/invariants/alexander.py), given a Wirtinger presentation $G = \langle g_1, \dots, g_n \mid R_1, \dots, R_m \rangle$ (with Wirtinger crossing relators from $L_-(D)$ plus saddle identification relators across the marked vertices), we form the $m \times n$ Fox Jacobian $J_{ij} = \psi(\partial R_i / \partial g_j) \in \Lambda$ and extract the $k$-th **elementary ideal** $E_k(F) \subset \Lambda$ generated by all $(n - k) \times (n - k)$ minors of $J$ (with $E_k = \Lambda$ for $n - k \le 0$ and $E_k = (0)$ for $n - k > m$), simplified strictly via **integral Hermite Normal Form over $\mathbb{Z}$** on bounded monomial shifts modulo $u_k^2 - 1$ without ever dividing by $2$ or projecting to $\mathbb{Q}$.

### 4.2 Hand-Checkable Proofs of All Three Errata in Yoshikawa (1994, Table II)

Using the canonical diagrams in [`yoshikawa_table.md`](./yoshikawa_table.md) (independently confirming Ishii–Nikkuni–Oshiro 2018, §5), we give self-contained algebraic proofs of all three Table II corrections:

#### Proposition 4.1 (Erratum #1: $9_1^{1,-2}$, $T^2 \sqcup K^2$)
*For Yoshikawa's ribbon link $9_1^{1,-2}$ (`Y[1,2,3,4], Z[4,5,6,1], X[6,7,8,9], X[9,8,10,2], X[11,12,13,3], X[12,14,5,13], Z[14,15,16,7], Y[10,17,18,11], X[18,17,16,15]`), the fundamental group and first elementary ideal over $\Lambda = \mathbb{Z}[x^{\pm 1}, y]/(y^2 - 1)$ are*
$$\pi_1(\mathbb{R}^4 \setminus 9_1^{1,-2}) = \langle a, b \mid a b = b a,\; b^2 = 1 \rangle \cong \mathbb{Z} \oplus \mathbb{Z}_2, \qquad E_1(9_1^{1,-2}) = (x - 1,\; y - 1,\; 2),$$
*correcting Yoshikawa's Table II entry $\langle x, y \mid yxyx^{-1} = 1 \rangle$ and $E_1 = (x + 1, y - 1)$.*

*Proof.* Tietze-simplifying the Wirtinger + saddle presentation of $9_1^{1,-2}$ yields two generators $a$ (meridian of $T^2$, $\psi(a) = x$) and $b$ (meridian of $K^2$, $\psi(b) = y$ with $y^2 = 1$) subject to the two relators
$$r_1 = b^{-1} a^{-1} b a = 1, \qquad r_2 = b^{-1} a^{-1} b^{-1} a = 1.$$
From $r_1 = 1$ we obtain $b a = a b$, so $G$ is abelian; substituting $a^{-1} b^{-1} a = b^{-1}$ into $r_2$ gives $b^{-2} = 1$, i.e., $b^2 = 1$, whence $G \cong \mathbb{Z} \oplus \mathbb{Z}_2$ (witnessed independently by $|\operatorname{Hom}(G, S_3)| = 1 + 3 + 2 + 6 = 12$, whereas $\langle x, y \mid yxyx^{-1} = 1 \rangle$ for $7_1^{0,-2}$ has $18$ $S_3$-representations). Taking Fox derivatives of $r'_1 = a b a^{-1} b^{-1}$ and $r'_2 = b^2$ gives the $2 \times 2$ Jacobian
$$J = \begin{pmatrix} 1 - y & x - 1 \\ 0 & 1 + y \end{pmatrix}.$$
Its $1 \times 1$ minors generate $E_1 = (y - 1,\; x - 1,\; y + 1) = (x - 1,\; y - 1,\; 2)$ since $(y + 1) - (y - 1) = 2$. $\square$

#### Proposition 4.2 (Erratum #2: $10_1^{0,1}$, $S^2 \sqcup T^2$)
*For $10_1^{0,1}$, whose group presentation is $G = \langle g_1, g_2 \mid r = 1 \rangle$ with $r = g_2^{-1} g_1^{-1} g_2 g_1^{-1} g_2 g_1 g_2^{-1} g_1$ and $\psi(g_1) = x, \psi(g_2) = y \in \mathbb{Z}[x^{\pm 1}, y^{\pm 1}]$, the first elementary ideal is*
$$E_1(10_1^{0,1}) = (x - 1,\; y - 1)(x + y) = (x^2 - y^2,\; x y - x + y^2 - y),$$
*correcting Yoshikawa's Table II entry $(x - 1, y - 1)(x + 1)$ and separating $10_1^{0,1}$ from $10_2^{0,1}$ by $E_1$ alone.*

*Proof.* Reading the eight prefixes of $r = g_2^{-1} g_1^{-1} g_2 g_1^{-1} g_2 g_1 g_2^{-1} g_1$ in $\mathbb{Z}[x^{\pm 1}, y^{\pm 1}]$ gives
$$\psi\Bigl(\frac{\partial r}{\partial g_1}\Bigr) = -x^{-1} y^{-1} - x^{-2} + x^{-2} y + x^{-1}, \qquad \psi\Bigl(\frac{\partial r}{\partial g_2}\Bigr) = -y^{-1} + x^{-1} y^{-1} + x^{-2} - x^{-1}.$$
Multiplying by the units $x^2 y$ and $-x^2 y$ in $\mathbb{Z}[x^{\pm 1}, y^{\pm 1}]$ yields
$$x^2 y \cdot \psi\Bigl(\frac{\partial r}{\partial g_1}\Bigr) = x y - x + y^2 - y = (x + y)(y - 1), \qquad -x^2 y \cdot \psi\Bigl(\frac{\partial r}{\partial g_2}\Bigr) = x^2 + x y - x - y = (x + y)(x - 1).$$
Hence $E_1(10_1^{0,1}) = \bigl((x + y)(x - 1),\; (x + y)(y - 1)\bigr) = (x - 1, y - 1)(x + y)$ (or $(x - 1, y - 1)(x y + 1)$ under $y \mapsto y^{-1}$). Evaluating at $(x, y) = (-1, -1)$ annihilates $E_1(10_2^{0,1}) = (x - 1, y - 1)(x + 1)$ under every $x^{\pm 1}, y^{\pm 1}$ orientation flip while sending $(x + y)(x - 1) \mapsto (-2)(-2) = 4 \neq 0$, so $E_1(10_1^{0,1}) \neq E_1(10_2^{0,1})$. $\square$

#### Theorem 4.3 (Erratum #3 and the Mixed 2-Torsion Obstruction: $10_2^{0,-2}$ vs $10_1^{0,-2}$)
*Let $\Lambda = \mathbb{Z}[x^{\pm 1}, y]/(y^2 - 1)$ be the Alexander group ring of an $S^2 \sqcup K^2$ link ($x$ orientable meridian, $y^2 = 1$ Klein-bottle meridian). For $10_2^{0,-2}$, with group presentation*
$$G = \langle x, y_1, y_2 \mid x y_1 x^{-1} = y_2,\; y_1^2 = y_2^2,\; y_1^2 = (y_1 y_2)^2 \rangle,$$
*the exact first elementary ideal is*
$$E_1(10_2^{0,-2}) = \bigl(2x + y - 1,\; (x + 1)(y + 1),\; 4\bigr) \;\supsetneq\; (2x + y - 1,\; 4) = E_1(10_1^{0,-2}),$$
*correcting Yoshikawa's omission of $(x + 1)(y + 1)$.*

*Proof.* Eliminating $y_2 = x y_1 x^{-1}$ leaves $G = \langle x, y_1 \mid R_1 = 1,\; R_2 = 1 \rangle$ with $R_1 = y_1^2 x y_1^{-2} x^{-1}$ and $R_2 = y_1^{-1} x y_1 x^{-1} y_1 x y_1 x^{-1}$ ($\psi(x) = x$, $\psi(y_1) = y$). Computing the $2 \times 2$ Fox Jacobian in $\Lambda$ gives
$$J = \begin{pmatrix} 0 & -(x - 1)(y + 1) \\ 2(y - 1) & y(2x + y - 1) \end{pmatrix}.$$
Since $y$ is a unit ($y^2 = 1$), the $1 \times 1$ minors generate $I = \bigl(2x + y - 1,\; 2(y - 1),\; (x - 1)(y + 1)\bigr)$. Observe two exact polynomial identities in $\Lambda = \mathbb{Z}[x^{\pm 1}, y]/(y^2 - 1)$:
1. $(y + 1)(2x + y - 1) = 2x(y + 1) + (y^2 - 1) = 2x(y + 1)$; multiplying by $x^{-1}$ puts $2(y + 1) \in I$, and subtracting $2(y - 1)$ puts $4 \in I$.
2. Conversely, $(2x - y + 1)(2x + y - 1) = 4x^2 - (y - 1)^2 = 4x^2 + 2(y - 1)$, so $2(y - 1)$ and $2(y + 1) = 2(y - 1) + 4$ already lie in $(2x + y - 1,\; 4)$! Adding $2(y + 1)$ to $(x - 1)(y + 1)$ yields $(x + 1)(y + 1)$, proving $I = \bigl(2x + y - 1,\; (x + 1)(y + 1),\; 4\bigr)$.

To prove strictly that $(x + 1)(y + 1) \notin (2x + y - 1,\; 4)$ (under both orientations $x$ and $x^{-1}$), note why hand calculation easily misses it: under both character projections $y \mapsto +1$ and $y \mapsto -1$ into $\mathbb{Z}[x^{\pm 1}]$, both ideals map to the *same* images $(2)$ and $(2x - 2, 4)$. Instead, consider the finite ring $Q = \mathbb{Z}_4[x]/(x^2 + x + 1)$ (where $x^3 = 1$, so $x$ is a unit with $x^{-1} = x^2$), and define the ring homomorphism $\phi: \Lambda \to Q$ by
$$\phi(x) = x, \qquad \phi(y) = 1 + 2x \qquad \bigl(\text{well-defined since } (1 + 2x)^2 = 1 + 4x + 4x^2 \equiv 1 \pmod 4\bigr).$$
Then $\phi(2x + y - 1) = 2x + (1 + 2x) - 1 = 4x = 0$ and $\phi(4) = 0$, so $\phi\bigl(E_1(10_1^{0,-2})\bigr) = \{0\}$. Yet on $(x + 1)(y + 1)$ we obtain
$$\phi\bigl((x + 1)(y + 1)\bigr) = (x + 1)(2 + 2x) = 2(x + 1)^2 = 2(x^2 + 2x + 1) \equiv 2x \neq 0 \quad \text{in } Q!$$
Under $x \mapsto x^{-1} = x^2$, replacing $\phi(y)$ by $1 + 2x^2$ likewise kills $(2x^{-1} + y - 1, 4)$ while mapping $(x^{-1} + 1)(y + 1) \mapsto 2(x^2 + 1)^2 = 2(-x)^2 = 2x^2 \neq 0$. Thus $E_1(10_2^{0,-2}) \neq E_1(10_1^{0,-2})$ as ideals over $\Lambda$. $\square$

---

## 5. Purely 4-Dimensional Invariants Beyond $(\pi_1, E_k)$

Unlike classical knots in $S^3$ (whose complements are aspherical and determined up to mirror/orientation by peripheral $\pi_1$), surface-links in $S^4$ frequently share isomorphic fundamental groups, identical Alexander modules, and identical intrinsic topologies while remaining topologically distinct. Below are the five purely 4-dimensional invariant mechanisms implemented in `surfaceknots`.

### 5.1 Peripheral Longitudinal Holonomy on Higher-Genus Components (`linking.py`)
For a classical knot $K \subset S^3$, Artin's spun 2-knot $\tau^0(K) \subset S^4$, its 1-handle stabilization $\tau^0(K) \# T^2_{\text{std}}$, and Satoh's ribbon torus-knot $\operatorname{Tube}(K) \subset S^4$ all have the **exact same fundamental group** $\pi_1(S^3 \setminus K)$ and identical unoriented/oriented Alexander modules! What distinguishes $\operatorname{Tube}(K)$ from $\tau^0(K) \# T^2_{\text{std}}$ is the **peripheral group homomorphism**
$$\iota_*: \pi_1(\partial \nu(F_i)) \longrightarrow \pi_1(S^4 \setminus \nu(F)).$$

#### Proposition 5.1 (Peripheral Separation of $\operatorname{Tube}(3_1) = 10_1^1$ from $8_1 \# 2_1^1$)
*Let $F$ be an orientable torus-knot ($F \cong T^2$, $e(F) = 0$), so $\partial \nu(F) \cong T^3 = F \times S^1_\mu$ has $\pi_1(\partial \nu(F)) \cong \mathbb{Z}\langle \lambda_1, \lambda_2, \mu_0 \rangle$ where $\mu_0$ is a meridian and $\lambda_1, \lambda_2$ generate $\pi_1(F)$ (uniquely determined modulo $\langle \mu_0 \rangle$). Define the finite-holonomy invariant*
$$N_{\operatorname{SL}(2, \mathbb{F}_3)}^{\text{periph}}(F) = \#\bigl\{\, \rho \in \operatorname{Hom}(\pi_1(E), \operatorname{SL}(2, \mathbb{F}_3)) \;\big|\; \rho(\iota_*(\pi_1(\partial \nu(F)))) \text{ is non-cyclic, or } \operatorname{ord}(\rho(\lambda_j)) \nmid \operatorname{ord}(\rho(\mu_0)) \,\bigr\}.$$
*Then $N_{\operatorname{SL}(2, \mathbb{F}_3)}^{\text{periph}}(\tau^0(K) \# T^2_{\text{std}}) = 0$ for **every** spun 2-knot $\tau^0(K)$, whereas $N_{\operatorname{SL}(2, \mathbb{F}_3)}^{\text{periph}}(10_1^1) = 24$.*

*Proof.* On any trivial 1-handle stabilization $S \# T^2_{\text{std}}$, both generators of $\pi_1(T^2_{\text{std}})$ bound disks in the local unknotted 4-ball modulo the meridian (one inside the 1-handle core, one dual to the belt circle), so $\iota_*(\pi_1(\partial \nu(S \# T^2_{\text{std}}))) = \langle \mu_0 \rangle$ is cyclic generated by $\mu_0$, forcing $N_{\operatorname{SL}(2, \mathbb{F}_3)}^{\text{periph}} = 0$. Conversely, for $10_1^1 = \operatorname{Tube}(3_1)$, $\pi_1(E) = \langle a, b \mid a b a = b a b \rangle$, one torus generator is trivial, and the longitudinal core follows the classical 0-framed longitude $\lambda_{\text{cl}} = a^{-6}(a b)^3 \in [\pi_1, \pi_1]$. Among the $48$ non-abelian representations $\rho: \pi_1(3_1) \to \operatorname{SL}(2, \mathbb{F}_3)$ lifting the $24$ non-abelian $A_4$-representations, exactly $24$ map $\mu_0$ to an element of order $3$ (such as $A = \bigl(\begin{smallmatrix} 1 & 1 \\ 0 & 1 \end{smallmatrix}\bigr)$, $B = \bigl(\begin{smallmatrix} 1 & 0 \\ -1 & 1 \end{smallmatrix}\bigr)$) while mapping $\lambda_{\text{cl}} \mapsto A^{-6}(A B)^3 = -I$ of order $2$ (so $\rho(\lambda_{\text{cl}} \mu_0^k)$ has order $2$ or $6 \nmid 3$), whereas the other $24$ map $\mu_0 \mapsto -A$ of order $6$ with $(-A)^3 = -I = \rho(\lambda_{\text{cl}})$. Thus [`peripheral_sl2f3_profile`](../code/src/surfaceknots/invariants/linking.py) evaluates to `24` on $10_1^1$ vs `0` on $8_1 \# 2_1^1$ (and likewise separates `11_1^0,1` with `24` from `11_2^0,1` with `0` at $\operatorname{ch} = 11$). $\square$

### 5.2 Carter–Saito–Sanderson Mod-2 Double Linking $\operatorname{Dlk}(F_i, F_j)$ (`linking.py`)
In Yoshikawa's table, $8_1^{1,1}$ (0-twist torus-spun Hopf link) and $10_1^{1,1}$ (1-twist-spun Hopf link) are both orientable $T^2 \sqcup T^2$ links with **identical** $\pi_1 \cong \mathbb{Z} \oplus \mathbb{Z}$ and **identical** $E_k$ ideals for all $k \ge 0$. By Alexander duality $\widetilde{H}_q(S^4 \setminus T^2) \cong \widetilde{H}^{3-q}(T^2)$ ($H_1 = \mathbb{Z}, H_2 = \mathbb{Z}^2, H_3 = 0$), two homologically unlinked tori in $\mathbb{R}^4$ have vanishing primary linking, yet can link nontrivially via secondary Pontryagin–Thom framings in $\pi_4(S^2) \cong \mathbb{Z}_2$.

#### Theorem 5.2 (Sanderson's Double Linking via Sky-Potentials and Rokhlin's Spin Form)
*Let $F_i \sqcup F_j \subset \mathbb{R}^4$ be a 2-component sublink of disjoint orientable surfaces (or with $F_j \cong T^2$ unknotted). Let $J_{ij} \subset F_i \cap p_{xyz}^{-1}(p_{xyz}(F_j))$ be the locus of double-point curves where $F_i$ passes **under** $F_j$ in $\mathbb{R}^3$ across the movie $t \in [-1, +1]$; generically $J_{ij}$ is a collection of framed immersed circles in $F_i$. Equivalently, for unknotted tori $F_1, F_2 \cong T^2$ with $\pi_1(\mathbb{R}^4 \setminus (F_1 \sqcup F_2)) \cong \mathbb{Z} \oplus \mathbb{Z}$, let $q_{F_i}: H_1(F_i; \mathbb{Z}_2) \to \mathbb{Z}_2$ be the Rokhlin spin quadratic enhancement weighted by classical linking with $F_j$ across slices, implemented via exact sky-potential cross-difference sums along the movie in [`double_linking_mod2`](../code/src/surfaceknots/invariants/linking.py). Then $\operatorname{Dlk}(F_1, F_2) \in \mathbb{Z}_2$ is an ambient isotopy invariant satisfying:*
$$\operatorname{Dlk}(8_1^{1,1}) = (0, 0), \qquad \operatorname{Dlk}(10_1^{1,1}) = (1, 1).$$

*Proof Sketch.* Following Sanderson (1988) and Carter–Jelsovsky–Kamada–Saito (2001), the transverse intersection $M_1 \pitchfork M_2$ of Seifert 3-manifolds bounding $F_1$ and $F_2$ in $\mathbb{R}^4$ is a framed 1-manifold in $\mathbb{R}^4$, whose total Pontryagin–Thom framing in $\pi_1(\mathrm{SO}) = \mathbb{Z}_2$ modulo spin cobordism is independent of the choice of $M_1, M_2$ whenever primary triple/surface linking vanishes. During the 1-twist spin $10_1^{1,1}$, the Hopf link rotates once around its axis as $\theta \in [0, 2\pi]$, twisting the framing of the core intersection circle $S^1_\theta$ by one full revolution ($1 \in \pi_1(\mathrm{SO}(3)) = \mathbb{Z}_2$), whereas the 0-twist spin $8_1^{1,1}$ carries the product framing ($0 \in \mathbb{Z}_2$). $\square$

### 5.3 Movie-Induced Mochizuki and Kamada–Oshiro Symmetric Quandle 3-Cocycles (`quandle.py`, `tier3.py`)
Let $(X, *)$ be a finite quandle and $\theta \in Z^3_Q(X; A)$ a rack/quandle 3-cocycle satisfying $\theta(x, x, y) = \theta(x, y, y) = 0$ and the 3-cocycle identity for broken-sheet triple points:
$$\theta(x, y, z) + \theta(x * z, y * z, w) + \theta(x, z, w) = \theta(x * y, z, w) + \theta(x, y, w) + \theta(x * w, y * w, z * w).$$

#### Proposition 5.3 (Constructive Movie State-Sum Formula without Broken Surface Diagrams)
*Let $D$ be an oriented admissible `ch`-diagram, with explicit Reidemeister unlinking movies $M_-: L_-(D) \to O_{m_-}$ and $M_+: L_+(D) \to O_{m_+}$ into crossingless planar unlinks. Every $X$-coloring $C \in \operatorname{Col}_X(D)$ restricts to $X$-colorings of $L_-(D)$ and $L_+(D)$ and propagates uniquely along every R1, R2, and R3 step of $M_-$ and $M_+$. Assigning Boltzmann weight $\pm \theta(x, y, z) \in A$ to each Reidemeister III move $\tau \in M_\pm$ (where $x, y, z$ are the bottom, middle, and top sheet colors in the source region and the sign combines crossing writhes and movie time direction), the state-sum*
$$\Phi_\theta(F) = \sum_{C \in \operatorname{Col}_X(D)} \exp\Biggl(\sum_{\tau \in M_+} \varepsilon_+(\tau)\,\theta(C, \tau) \;-\; \sum_{\tau \in M_-} \varepsilon_-(\tau)\,\theta(C, \tau)\Biggr) \in \mathbb{Z}[A]$$
*equals the Carter–Jelsovsky–Kamada–Langford–Saito (2003) shadow/triple-point 3-cocycle invariant of $F(D)$, and is independent of the choice of unlinking movies $M_\pm$ and Yoshikawa moves $\Omega_1$–$\Omega_8$.*

*Proof.* Sweeping the planar projection of the movie $M_-^{-1} \circ (\text{saddles}) \circ M_+$ along $t \in [-1, +1]$ builds a generic broken surface diagram of $F(D)$ in $\mathbb{R}^3$ whose **only** triple points occur at the Reidemeister III steps of $M_-$ (traversed forward from $O_{m_-}$ to $L_-(D)$, hence negated relative to $L_- \to O_{m_-}$) and $M_+$ (traversed forward from $L_+(D)$ to $O_{m_+}$), because saddle bands and R1/R2 moves introduce zero triple points and capping crossingless unlinks introduces zero crossings or triple points. By CJKS (2003), $\Phi_\theta(F)$ depends only on the oriented ambient isotopy class of $F \subset \mathbb{R}^4$. In particular, if $\Phi_\theta(F) \notin \mathbb{Z} \cdot 1_A$, then $F$ has **positive triple-point number $t(F) > 0$** (hence cannot be a ribbon 2-knot, since Yajima proved $t(F_{\text{ribbon}}) = 0$ for ribbon 2-knots), and if $\Phi_\theta(F^*) = \overline{\Phi_\theta(F)} \neq \Phi_\theta(F)$ (as on $10_2 = \tau^2(3_1)$ with Mochizuki's 3-cocycle over $R_3$, where $\Phi_{\theta_3}(10_2) = 3 + 6u$ vs $\Phi_{\theta_3}(10_2^*) = 3 + 6u^2$), then $F$ is **4D chiral**. $\square$

#### Definition 5.4 (Kamada–Oshiro Symmetric Quandles for Non-Orientable and Mixed Surface-Links)
When $F$ has non-orientable components ($\mathbb{R}P^2, K^2$), classical oriented quandle colorings require co-orientations that fail globally around orientation-reversing loops on $F$. Kamada (2007) and Kamada–Oshiro (2010) resolved this by introducing a **symmetric quandle** $(X, *, \rho)$, where $\rho: X \to X$ is a **good involution** satisfying, for all $x, y \in X$:
$$\rho(x * y) = \rho(x) * y, \qquad x * \rho(y) = x *^{-1} y.$$
Reversing the normal co-orientation of a sheet simultaneously replaces its color $x$ by $\rho(x)$. The smallest non-dihedral symmetric quandle capturing non-abelian $\mathbb{R}P^2/K^2$ linking is the **octahedral symmetric quandle** $(Q_6, *, \rho)$, realized as the six unit-quaternion 4th roots of unity
$$Q_6 = \{\pm i,\; \pm j,\; \pm k\} \subset Q_8, \qquad x * y = y^{-1} x y, \qquad \rho(x) = x^{-1} = -x.$$
Every symmetric 3-cocycle $\theta \in Z^3_{Q,\rho}(Q_6; A)$ satisfies both the quandle 3-cocycle condition and the involution antisymmetry relations, yielding an **unoriented ambient isotopy invariant** $\Phi_\theta^{\text{sym}}(F) \in \mathbb{Z}[A]$ defined on **every** surface-link $F$ — orientable, non-orientable, or mixed — directly via Proposition 5.3. On $8_1^{-1,-1}$ ($\pi_1 \cong Q_8$), Kamada–Oshiro proved that $(Q_6, \rho)$ colorings and symmetric 3-cocycles obstruct $t = 0$ and establish $t(8_1^{-1,-1}) = 2$.

### 5.4 Second Homotopy $\pi_2(S^4 \setminus F)$ and Poincaré–Lefschetz Cover Homology (`tier3.py`)
Whereas $\pi_2(S^3 \setminus K) = 0$ for every non-split classical link (by Papakyriakopoulos's Sphere Theorem), Lomonaco (1981) emphasized that $\pi_2(S^4 \setminus F)$ is frequently nontrivial and carries rich $\mathbb{Z}[\pi_1(E)]$-module structure. By combining Fox calculus on finite quotients with Poincaré–Lefschetz duality on the 4-manifold exterior, we can compute exact second Betti numbers of regular covers — and $\operatorname{rk}_{\mathbb{Z}} \pi_2(E)$ for finite-commutator surface-knots — purely from 1-skeleton / 2-skeleton presentation data without building 3-cells.

#### Theorem 5.5 (Poincaré–Lefschetz Formula for $b_2(\widetilde{E}_\rho)$, $b_3(\widetilde{E}_\rho)$, and $\operatorname{rk}_{\mathbb{Z}}\pi_2(E)$)
*Let $F = F_1 \sqcup \dots \sqcup F_c \subset S^4$ ($F \neq \emptyset$) be any surface-link with exterior $E = S^4 \setminus \nu(F)$, and let $\rho: \pi_1(E) \twoheadrightarrow \Gamma$ be a **surjective** homomorphism onto a finite group $\Gamma$ of order $d = |\Gamma|$, defining a connected $d$-sheeted regular cover $\widetilde{E}_\rho \to E$ with boundary $\partial \widetilde{E}_\rho$. Then:*
1. *The integral first homology $H_1(\widetilde{E}_\rho; \mathbb{Z}) \cong \mathbb{Z}^{b_1(\widetilde{E}_\rho)} \oplus \operatorname{Tors} H_1(\widetilde{E}_\rho; \mathbb{Z})$ is presented by the $d m \times d n$ lifted Fox Jacobian $\partial_2^\rho: \mathbb{Z}[\Gamma]^m \to \mathbb{Z}[\Gamma]^n$ modulo $\ker(\partial_1^\rho)$.*
2. *The higher Betti numbers of $\widetilde{E}_\rho$ are determined completely by peripheral 1-homology via*
   $$b_3(\widetilde{E}_\rho) = \bigl(b_0(\partial \widetilde{E}_\rho) - 1\bigr) + b_1(\widetilde{E}_\rho) - \operatorname{rk}_{\mathbb{Q}}\!\Bigl(\operatorname{im}\bigl(H_1(\partial \widetilde{E}_\rho; \mathbb{Q}) \to H_1(\widetilde{E}_\rho; \mathbb{Q})\bigr)\Bigr),$$
   $$b_2(\widetilde{E}_\rho) = d\bigl(2 - \chi(F)\bigr) - 1 + b_1(\widetilde{E}_\rho) + b_3(\widetilde{E}_\rho), \qquad b_4(\widetilde{E}_\rho) = 0.$$
3. *In particular, if $\ker(\rho)$ is either trivial ($\pi_1(E) \cong \Gamma$ finite, so $\widetilde{E}_\rho = \widetilde{E}_{\text{univ}}$) or infinite cyclic $\ker(\rho) \cong \mathbb{Z}$ acting trivially on rational $H_2$ (as for fibered twist-spun 2-knots with finite commutator subgroup $[\pi_1(E), \pi_1(E)]$), Hurewicz's theorem gives $\operatorname{rk}_{\mathbb{Z}} \pi_2(E) = b_2(\widetilde{E}_\rho)$ when $\ker(\rho) = 1$, and $\operatorname{rk}_{\mathbb{Z}} \pi_2(E) = b_2(\widetilde{E}_\rho)$ whenever the Wang sequence of $\widetilde{E}_{\text{univ}} \to \widetilde{E}_\rho \to S^1$ has trivial monodromy on $H_2$ and $H_1(\widetilde{E}_{\text{univ}}) = 0$.*

*Proof.* Because $E \subset S^4$ is a compact connected orientable 4-manifold with non-empty boundary $\partial E = \partial \nu(F)$, its connected $d$-sheeted cover $\widetilde{E}_\rho$ is also a compact connected orientable 4-manifold with non-empty boundary $\partial \widetilde{E}_\rho$. By Poincaré–Lefschetz duality over $\mathbb{Q}$, $H_4(\widetilde{E}_\rho; \mathbb{Q}) \cong H^0(\widetilde{E}_\rho, \partial \widetilde{E}_\rho; \mathbb{Q}) = 0$ and $H_3(\widetilde{E}_\rho; \mathbb{Q}) \cong H^1(\widetilde{E}_\rho, \partial \widetilde{E}_\rho; \mathbb{Q})$. In the cohomology exact sequence
$$0 \longrightarrow H^0(\widetilde{E}_\rho; \mathbb{Q}) \xrightarrow{\;i^0\;} H^0(\partial \widetilde{E}_\rho; \mathbb{Q}) \longrightarrow H^1(\widetilde{E}_\rho, \partial \widetilde{E}_\rho; \mathbb{Q}) \longrightarrow H^1(\widetilde{E}_\rho; \mathbb{Q}) \xrightarrow{\;i^1\;} H^1(\partial \widetilde{E}_\rho; \mathbb{Q}),$$
$\operatorname{coker}(i^0)$ has dimension $b_0(\partial \widetilde{E}_\rho) - 1$ while $\ker(i^1)$ has dimension $b_1(\widetilde{E}_\rho) - \operatorname{rk}_{\mathbb{Q}}(i^1) = b_1(\widetilde{E}_\rho) - \operatorname{rk}_{\mathbb{Q}}(i_1)$, proving the formula for $b_3(\widetilde{E}_\rho)$. Multiplicativity of the Euler characteristic $\chi(\widetilde{E}_\rho) = d \cdot \chi(E) = d(2 - \chi(F)) = 1 - b_1(\widetilde{E}_\rho) + b_2(\widetilde{E}_\rho) - b_3(\widetilde{E}_\rho)$ then yields $b_2(\widetilde{E}_\rho)$.

*Canonical Examples*:
- **$2_1^{-1} = P_+$** ($\chi = 1$, $\pi_1(E) \cong \mathbb{Z}_2$, universal cover $d = 2$): $b_0(\partial \widetilde{E}) = 1$, $b_1(\widetilde{E}) = 0 \implies b_3(\widetilde{E}) = 0$ and $\operatorname{rk}_{\mathbb{Z}} \pi_2(S^4 \setminus P_+) = b_2(\widetilde{E}) = 2(2 - 1) - 1 = \mathbf{1}$.
- **$10_2 = \tau^2(3_1)$** ($\chi = 2$, $\pi_1(E) \cong \mathbb{Z} \ltimes \mathbb{Z}_3 \twoheadrightarrow S_3$, $d = 6$, $\operatorname{ord}(\rho(\mu_0)) = 2$, $\ker(\rho) = 2\mathbb{Z}$): $b_0(\partial \widetilde{E}_{S_3}) = 6/2 = 3$, $b_1 = 1$, $\operatorname{rk}(i_1) = 1 \implies b_3 = (3 - 1) + 1 - 1 = 2$ and $\operatorname{rk}_{\mathbb{Z}} \pi_2(10_2) = b_2 = 0 - 1 + 1 + 2 = \mathbf{2}$ (matching the punctured lens-space fiber $L(3, 1)^\circ$, whose universal cover $S^3 \setminus \{3\text{ pts}\}$ has $H_2 \cong \mathbb{Z}^2$).
- **$10_3 = \tau^3(3_1)$** ($\chi = 2$, $[\pi_1, \pi_1] \cong Q_8$, $\rho: \pi_1(E) \twoheadrightarrow \operatorname{SL}(2, \mathbb{F}_3)$, $d = 24$, $\operatorname{ord}(\rho(\mu_0)) = 3$, $\ker(\rho) = 3\mathbb{Z}$): $b_0(\partial \widetilde{E}_\rho) = 24/3 = 8$, $b_1 = 1$, $\operatorname{rk}(i_1) = 1 \implies b_3 = (8 - 1) + 1 - 1 = 7$ and $\operatorname{rk}_{\mathbb{Z}} \pi_2(10_3) = b_2 = 0 - 1 + 1 + 7 = \mathbf{7}$ (recovering Lomonaco's $\pi_2(\tau^3(3_1)) \cong \mathbb{Z}^7$, the punctured quaternionic space $S^3/Q_8 \setminus \{\text{pt}\}$). $\square$

### 5.5 Total Mod-2 Double Branched Cover $\Sigma_2(S^4, F)$ (`tier3.py`)
Every surface-link $F = F_1 \sqcup \dots \sqcup F_c \subset S^4$ — **whether orientable, non-orientable, or mixed** — possesses a canonical closed, oriented, smooth 4-manifold invariant: the **total mod-2 double branched cover** $\Sigma_2(S^4, F) \xrightarrow{2:1} S^4$ branched along $F$, corresponding to the unique homomorphism $\epsilon: \pi_1(S^4 \setminus F) \twoheadrightarrow \mathbb{Z}_2$ sending every Wirtinger meridian $\mu \mapsto -1 \in \{\pm 1\}$.

#### Theorem 5.6 (Fox–Plans Homology and Hirzebruch–Atiyah–Singer Signature of $\Sigma_2(S^4, F)$)
*Let $F = F_1 \sqcup \dots \sqcup F_c \subset S^4$ be any surface-link with total Euler characteristic $\chi(F) = \sum_{i=1}^c \chi(F_i)$, component normal Euler numbers $e(F_1), \dots, e(F_c)$, and universal kei coloring module invariants $(r_{\text{kei}}, T_{\text{kei}})$ from [`kei_module_invariants`](../code/src/surfaceknots/invariants/quandle.py) ($r_{\text{kei}} \ge 1$ for $F \neq \emptyset$). Then the closed oriented 4-manifold $\Sigma_2 = \Sigma_2(S^4, F)$ satisfies:*
1. **First homology**:
   $$H_1(\Sigma_2; \mathbb{Z}) \cong \mathbb{Z}^{\,r_{\text{kei}} - 1} \;\oplus\; \bigoplus_{t \in T_{\text{kei}}} \mathbb{Z}_t, \qquad b_1(\Sigma_2) = r_{\text{kei}} - 1.$$
2. **Euler characteristic and second Betti number**:
   $$\chi(\Sigma_2) = 4 - \chi(F), \qquad b_2(\Sigma_2) = 2 + 2(r_{\text{kei}} - 1) - \chi(F).$$
3. **4-Manifold signature and definite Betti numbers** (*Hirzebruch 1969; Gordon 1978*):
   $$\sigma(\Sigma_2) = -\frac{1}{2} \sum_{i=1}^c e(F_i), \qquad b_2^\pm(\Sigma_2) = \frac{b_2(\Sigma_2) \pm \sigma(\Sigma_2)}{2}.$$

*Proof.* (1) Let $\widetilde{E}_2 \to E$ be the unbranched double cover associated to $\epsilon: \pi_1(E) \twoheadrightarrow \mathbb{Z}_2 = \langle \tau \mid \tau^2 = 1 \rangle$. Filling back the $c$ lifted tubular neighborhoods completes $\widetilde{E}_2$ to $\Sigma_2$ by killing $\widetilde{\mu}_i = \mu_i^2$ for every meridian $\mu_i$. In the group presentation of $\pi_1(E)$, mapping every Wirtinger meridian $x_j$ to a reflection $(c_j, -1) \in A \rtimes_{-1} \mathbb{Z}_2$ identifies $\operatorname{Hom}_{\text{merid}}(\pi_1(E), A \rtimes_{-1} \mathbb{Z}_2) / A_{\text{trans}} \cong \operatorname{Hom}(H_1(\Sigma_2; \mathbb{Z}), A)$ (Fox 1956; Fenn–Rourke 1992): the translation action removes one free $\mathbb{Z}$ factor of trivial constant colorings, leaving $H_1(\Sigma_2; \mathbb{Z}) \cong \mathbb{Z}^{r_{\text{kei}} - 1} \oplus \bigoplus_{t \in T_{\text{kei}}} \mathbb{Z}_t$.

(2) By Riemann–Hurwitz for the 2-sheeted branched covering $\Sigma_2 \to S^4$ with branch locus $F$, $\chi(\Sigma_2) = 2\chi(S^4 \setminus F) + \chi(F) = 2\chi(S^4) - \chi(F) = 4 - \chi(F)$. Since $\Sigma_2$ is a closed connected orientable 4-manifold, Poincaré duality gives $b_0(\Sigma_2) = b_4(\Sigma_2) = 1$ and $b_3(\Sigma_2) = b_1(\Sigma_2)$, whence $b_2(\Sigma_2) = \chi(\Sigma_2) - 2 + 2b_1(\Sigma_2) = 2 + 2(r_{\text{kei}} - 1) - \chi(F)$.

(3) Let $\tau: \Sigma_2 \to \Sigma_2$ be the covering involution with fixed-point surface $\operatorname{Fix}(\tau) = \widetilde{F} \cong F$. Because $\Sigma_2 / \tau = S^4$ has $\sigma(S^4) = 0$, the eigenspace decomposition of $H^2(\Sigma_2; \mathbb{R}) = H^2_+ \oplus H^2_-$ has $H^2_+ \cong H^2(S^4; \mathbb{R}) = 0$, so $\tau^*$ acts as $-1$ on all of $H^2(\Sigma_2; \mathbb{R})$! Therefore $\operatorname{Sign}(\tau, \Sigma_2) = -\sigma(\Sigma_2)$. On the other hand, the Atiyah–Singer / Hirzebruch $G$-Signature Theorem for involutions states that $\operatorname{Sign}(\tau, \Sigma_2) = e(\nu_{\widetilde{F} \subset \Sigma_2}) = \frac{1}{2} e(\nu_{F \subset S^4}) = \frac{1}{2} \sum_{i=1}^c e(F_i)$ (since the normal bundle to $\widetilde{F}$ in $\Sigma_2$ squares under $z \mapsto z^2$ to the normal bundle of $F$ in $S^4$). Combining the two identities yields $\sigma(\Sigma_2) = -\frac{1}{2} \sum_{i=1}^c e(F_i)$. Note the sanity checks:
- For $2_1^{-1} = P_+$ ($e = +2, \chi = 1, r_{\text{kei}} = 1$): $H_1(\Sigma_2) = 0$, $b_2(\Sigma_2) = 1$, $\sigma(\Sigma_2) = -1$, so $b_2^+ = 0, b_2^- = 1$, matching $\Sigma_2(S^4, P_+) \cong \overline{\mathbb{C}P}^2$!
- For any orientable 2-sphere knot $F \cong S^2$ ($e = 0, \chi = 2, r_{\text{kei}} = 1$): $b_1(\Sigma_2) = b_2(\Sigma_2) = \sigma(\Sigma_2) = 0$ ($\Sigma_2$ is a $\mathbb{Z}_2$-homology 4-sphere with $|H_1(\Sigma_2; \mathbb{Z})| = \prod_{t \in T_{\text{kei}}} t = |\Delta_F(-1)|$, the classical determinant of $F$).
- At $\operatorname{ch} = 11$, `11_2^0,-2` ($r_{\text{kei}} = 2, T_{\text{kei}} = ()$) has $\Sigma_2$ with $b_1 = 1, b_2^+ = 1, b_2^- = 1$, whereas `11_3^0,-2` ($r_{\text{kei}} = 1, T_{\text{kei}} = (4,)$) has $\Sigma_2$ a **rational homology 4-sphere** with $H_1(\Sigma_2; \mathbb{Z}) \cong \mathbb{Z}_4$ and $b_1 = b_2 = \sigma = 0$! $\square$
