# Literature Review: Tabulating Knotted Surfaces in $\mathbb{R}^4$

This document surveys existing censuses, diagrammatic representations, invariants, and complexity cutoffs for knotted surfaces (surface-knots and surface-links) in $\mathbb{R}^4$ (equivalently $S^4$), laying the mathematical foundation for an expanded [KnotInfo](https://knotinfo.org/) / [LinkInfo](https://linkinfo.sitehost.iu.edu/)-style database.

---

## 1. Existing Tables and Censuses

Unlike classical knots—tabulated past 19–20 crossings (hundreds of millions of prime knots)—tabulation in dimension 4 is remarkably sparse. Existing censuses fall into two camps: **one general enumeration across all topologies** (stopping at $\operatorname{ch} \le 10$) and **restricted sub-tables** of special geometric families (mostly ribbon surfaces).

| Census / Authors | Scope (Topological Type) | Complexity Measure & Cutoff | Count & Status |
| :--- | :--- | :--- | :--- |
| **Yoshikawa (1994)** | All weakly prime, non-splittable surface-links (orientable & non-orientable, any # of components) | Marked graph $\operatorname{ch}$-index $\operatorname{ch}(F) \le 10$ | **23** surface-links (6 are $S^2$-knots, 2 are $T^2$-knots, 1 is $\mathbb{R}P^2$, 14 are multi-component links). **No general $\operatorname{ch} = 11$ enumeration exists in print.** |
| **T. Yasuda (2001, 2018); Kanenobu–Sumi (2019–20)** | Oriented **ribbon 2-knots** ($S^2 \hookrightarrow \mathbb{R}^4$) | Ribbon crossing number $r\text{-cr}(K) \le 4$ | **129** ribbon 2-knots ($1 + 3 + 13$ for $r\text{-cr} \le 3$; $112$ for $r\text{-cr} = 4$, after Kanenobu–Sumi found 1 missed knot $Y_{112}$). |
| **Kanenobu–Komatsu (2017); Kanenobu–Takahashi (2021)** | Oriented **ribbon 2-knots** / **1-fusion ribbon 2-knots** | Virtual arc crossings $\le 4$; 1-fusion presentation length $\le 6$ | Complete classifications up to orientation/mirror using Alexander & twisted Alexander invariants. |
| **J. Yasuda (2025)** | Oriented **2-plat 2-knots** $F(p,q)$ (plat closures of 2D 4-braids; all 1-fusion ribbon $S^2$) | Odd determinant $p = |\Delta_F(-1)| \le 19$ | Complete classification theorem for all $(p,q)$ via $(p, \tau(F))$; detects non-invertibility. |
| **Satoh (2000); Ichimori–Kanenobu (2012)** | **Ribbon $T^2$-knots** via $\operatorname{Tube}(K)$ of virtual/welded knots | Virtual crossing number $\le 4$ (and induced tables up to 6) | Surjective onto ribbon $T^2$-knots, though $\operatorname{Tube}(K)$ is not injective on virtual knots (welded moves). |
| **Kamada (2002); Nagase–Shima (2018–23)** | Oriented surface-links from **2D braids / charts** | Braid index $m \le 4$, white vertices $w = 0$ (or minimal charts with $w \le 2$) | $w = 0$ charts classify low-braid-index simple ribbon 2-knots; $w \ge 1$ is only partially classified. |
| **Meier–Zupan (2017); Allred et al. (2023); Gong et al. (2026)** | Smooth surfaces in $S^4$ via **bridge trisections / tri-plane diagrams** | Bridge number $b(F)$ and tri-plane crossing number $c_{\text{tri}}(F)$ | $b \le 3 \implies$ unknotted; $c_{\text{tri}} = 0 \implies$ unknotted orientable; $c_{\text{tri}} \le 5$ for 2-knots $\implies$ ribbon; $c_{\text{tri}}(10_2) = 6$; conjectured $\{0_1, 8_1, 10_2\}$ are all 2-knots with $c_{\text{tri}} \le 6$. |
| **Ishikawa–Koda–Naoe (2024)** | Smooth **2-knots** ($S^2 \subset S^4$) via Turaev shadows | Exterior **shadow-complexity** $\operatorname{sc}(K) \le 1$ | $\operatorname{sc} = 0$: only $0_1$. $\operatorname{sc} = 1$: already an **infinite** family $K_n$ ($n \neq 0$, $K_{-1} = 8_1$, $K_1 = 9_1$). |

---

## 2. Diagrammatic & Combinatorial Representations

For machine enumeration, storage, and invariant calculation, we need representations that are **compactly serializable**, have **algorithmically checkable admissibility conditions**, and possess **complete move calculi**.

### 2.1 Marked Graph Diagrams (`ch`-Diagrams) & Dual Graphs
- **Geometry**: By Kawauchi–Shibuya–Suzuki (1982), every surface-link $F \subset \mathbb{R}^4$ admits a *hyperbolic splitting* where all Morse minima lie in $x_4 = -1$, all saddles (hyperbolic points) lie in $x_4 = 0$, and all maxima lie in $x_4 = +1$. Projecting the $x_4 = 0$ slice to $\mathbb{R}^2$ gives a 4-valent spatial graph diagram $D$ with two vertex types: **classical crossings** ($c(D)$) and **marked vertices** ($h(D)$ saddles, each equipped with a small bar indicating the $x_4 > 0$ smoothing).
- **Admissibility**: $D$ represents an embedded surface-link iff both the positive resolution $L_+(D)$ and negative resolution $L_-(D)$ are **classical unlinks** in $\mathbb{R}^3$. Because trivial disk systems bounding an unlink in $B^4$ are unique up to smooth isotopy rel boundary, $(D, L_\pm)$ uniquely determines $F \subset S^4$.
- **Move Calculus**: **8 Yoshikawa moves** ($\Omega_1$–$\Omega_5$ preserving $\chi$, plus $\Omega_6$–$\Omega_8$ canceling/sliding saddles with caps; completeness proved by Swenton 2001 and Kearton–Kurlin 2008).
- **Machine Encodings**:
  1. **Extended Planar Code (PD code)**: Standard 4-tuple edge codes `X[a, b, c, d]` for crossings plus `M[a, b, c, d]` for marked vertices (where `(a,c)` vs `(b,d)` encodes the marker).
  2. **Yoshikawa Edge-Labeled Planar Graphs $\Gamma$**: Checkerboard-coloring $S^2 \setminus D$ gives a planar graph (reducible to a *simple* graph without loops or multi-edges without increasing $\operatorname{ch}$-index, Yoshikawa Theorem 3.5) with integer/signed edge weights. Ideal for exhaustive enumeration!
  3. **Surface Singular Braid Words** (Jabłonowski 2013, 2017): Marked graph diagrams in closed braid form on $n$ strands using generators $\sigma_i^{\pm 1}$ and marked singular generators $\tau_i, \delta_i$. Compact 1-line string representation (given for all 23 Yoshikawa entries in `arXiv:1602.08650`).
  4. **Marked Graph Mosaics** (Choi–Nelson 2022): $m \times m$ grids over 16 tile types.

### 2.2 Banded Unlink Diagrams $(L, B)$
- **Geometry**: An unlink diagram $L \subset \mathbb{R}^2$ together with framed bands $B$ attached to $L$ such that band surgery $L_B$ is also an unlink. Shrinking each band to a marked vertex gives an exact bijection with `ch`-diagrams.
- **Move Calculus**: Band slides, band swims, cup/cap cancellations, and Kirby handle slides (Hughes–Kim–Miller 2020).
- **Why keep alongside `ch`-diagrams**: Directly produces **Kirby diagrams** for the exterior $S^4 \setminus \nu(F)$ (dotted circles for $L$, $0$-framed 2-handles from bands) and for branched covers $\Sigma_k(S^4, F)$, which are essential for computing 4-manifold and skein lasagna invariants.

### 2.3 Tri-Plane Diagrams $(D_1, D_2, D_3)$ (Bridge Trisections)
- **Geometry** (Meier–Zupan 2017): Three $b$-strand trivial tangle diagrams in $D^2$ sharing $2b$ boundary endpoints such that $D_1 \cup \overline{D}_2$, $D_2 \cup \overline{D}_3$, and $D_3 \cup \overline{D}_1$ are classical unlinks with $c_1, c_2, c_3$ components.
- **Move Calculus**: Interior Reidemeister moves inside $D_i$, mutual boundary braid transpositions across all three tangles simultaneously, and stabilization/destabilization.
- **Machine Encoding**: Three tangle PD codes or three $2b$-strand half-braid/matching tuples. Every `ch`-diagram in bridge position converts algorithmically into a tri-plane diagram (Allred et al. 2023, §4).

### 2.4 Broken Surface Diagrams & Motion Pictures
- **Geometry**: Generic projection $p(F) \subset \mathbb{R}^3$ with double curves, isolated Whitney umbrellas (branch points), and triple points, with under-sheets broken along double curves; complete **7 Roseman moves** (1998), or **15 Carter–Saito movie moves** on 1D slices.
- **Tradeoff**: Indispensable geometrically (defines triple point number $t(F)$, sheet number $\operatorname{sh}(F)$, and original quandle 3-cocycle state sums), **poor as a primary enumeration format** because encoding 2D sheet combinatorics and checking liftability (Blankinship–Fox / Giller obstructions) is heavy. Instead, an admissible `ch`-diagram + explicit Reidemeister unlinking sequences $L_\pm(D) \to O_\pm$ *induces* a broken surface diagram automatically (each R-III move becomes 1 triple point; each R-I move becomes 1 branch point).

### 2.5 2D Braids and Chart Descriptions
- **Geometry** (Kamada 2002): Finite labeled graphs $\Gamma \subset D^2$ with degree-1 black vertices (branch points), degree-4 crossings, and degree-6 white vertices (triple points), encoding the monodromy of a degree-$m$ simple branched covering $F \to S^2$; complete **C-moves** (CI–CIII) + Markov moves.
- **Limitation**: **Only represents oriented surfaces**—cannot represent $\mathbb{R}P^2$, Klein bottles, or non-orientable links. Best included as an auxiliary annotation column for orientable entries.

---

## 3. Taxonomy of Invariants for a KnotInfo-Style Surface Table

To build a rich, discriminatory database, we organize invariants into three computational tiers:

### Tier 1: Directly Computable in Code from a `ch`-Diagram (Fast, Exact)
1. **Intrinsic Topology**:
   - Number of components $c$, component homeomorphism types ($S^2$, $\Sigma_g$, $\#_k \mathbb{R}P^2$).
   - Orientability (orientable, non-orientable, or mixed link).
   - Euler characteristic: $\chi(F) = |O_+| + |O_-| - h(D)$ (or $c_1 + c_2 + c_3 - b$ from a tri-plane diagram), where $|O_\pm|$ are component counts of $L_\pm(D)$.
   - **Normal Euler number** $e(F)$ (and vector $(e(F_1), \dots, e(F_c))$ per component): $0$ for orientable components; satisfies Whitney–Massey congruence $e(F_i) \equiv 2\chi(F_i) \pmod 4$ and $|e(F_i)| \le 2(2 - \chi(F_i))$ when unknotted; directly computable from signed crossing writhe on `ch`-diagrams (Kamada 1989) or pairwise writhe $w(D_1 \cup \overline{D}_2) + w(D_2 \cup \overline{D}_3) + w(D_3 \cup \overline{D}_1)$ on tri-plane diagrams (Joseph et al. 2022).
2. **Fundamental Group $\pi_1(S^4 \setminus F)$ & Derived Algebraic Invariants**:
   - **Wirtinger presentation**: Over-arcs of $D$ generate $\pi_1$, classical crossings give standard Wirtinger relations $x_k = x_i^{\pm 1} x_j x_i^{\mp 1}$, and marked vertices identify opposing arcs!
   - **Abelianization** $H_1(S^4 \setminus F; \mathbb{Z}) \cong \mathbb{Z}^{c_{\text{or}}} \oplus (\mathbb{Z}/2\mathbb{Z})^{c_{\text{non}}}$.
   - **Finite-index subgroup / quotient counts**: Homomorphism counts $|\operatorname{Hom}(\pi_1, G)|$ to small finite groups ($S_3, S_4, A_4, A_5, D_{2p}, Q_8, \operatorname{PSL}(2, \mathbb{F}_q)$)—extremely fast and powerful isomorphism filter!
   - **Alexander Module & Elementary Ideals $E_k(F)$**: Via Fox free calculus on the Wirtinger presentation over $\mathbb{Z}[H_1]$ (or $\mathbb{Z}[H_1/\text{torsion}]$). Yields the **Alexander polynomial** $\Delta_F(t)$ (when $E_1$ is principal; note $\Delta_F(t)$ need *not* be symmetric in 4D, e.g. $\Delta_{9_1}(t) = t - 2$!), **surface determinant** $|\Delta_F(-1)|$, and **Nakanishi index**.
   - **Twisted Alexander Ideals / Polynomials** (Wada / Lin; Ishii–Nikkuni–Oshiro 2018): Associated to representations $\rho: \pi_1 \to \operatorname{SL}(2, \mathbb{F}_p)$ or $\operatorname{GL}(n, R)$; separates Yoshikawa surface-links whose untwisted $E_k$ coincide.
3. **Quandle, Symmetric Quandle, Biquandle & Cocycle State-Sums**:
   - **Coloring counts**: Fox $p$-colorings $|\operatorname{Col}_{R_p}(D)|$, finite quandle colorings (for oriented $F$), involutory quandle / **kei** colorings, **symmetric quandle** colorings $(X, \rho)$ (for non-orientable $F$), and **biquandle** colorings (Ashihara 2012 directly from `ch`-diagrams).
   - **Quandle 3-Cocycle Invariants $\Phi_\theta(F) \in \mathbb{Z}[A]$** (Carter–Jelsovsky–Kamada–Langford–Saito 1999): Computable *directly from `ch`-diagrams* via the **Kamada–Kim–Lee (2015)** algorithm (or by tracking R-III moves in unlinking $L_\pm(D)$)! Crucial because $\pi_1$ and Alexander invariants are blind to mirror images and often to orientation reversal, whereas $\Phi_\theta$:
     - Detects **non-invertibility** ($F \not\cong -F$) and **chirality** ($F \not\cong F^*$), e.g., for twist-spun trefoils $10_2, 10_3$.
     - Gives sharp **lower bounds on triple point number** $t(F)$ (including for non-orientable links via symmetric quandle 3-cocycles, Kamada–Oshiro 2010).

### Tier 2: Diagrammatic & Geometric Complexity Invariants (Exact Values or `[lower, upper]` Intervals, KnotInfo-Style)
- **$\operatorname{ch}$-index** $\operatorname{ch}(F) = \min_D (c(D) + h(D))$, plus minimal hyperbolic/saddle number $h_{\min}(F)$ and minimal Morse critical point count $\operatorname{crit}(F)$.
- **Triple point number** $t(F)$: Minimum triple points in a broken surface diagram ($t(F) = 0 \iff$ *pseudo-ribbon*; for 2-knots, Teragaito/Yajima proved $t(K)=0 \iff$ ribbon, and Satoh proved $t(K) \neq 1$; Kawamura proved $t(K) \notin \{2,3\}$ for 2-knots, so $t(10_2) = 4$ is the smallest nonzero value for a 2-knot!).
- **Bridge number** $b(F)$ and minimal patch vector $(c_1, c_2, c_3)$ in a bridge trisection.
- **Tri-plane crossing number** $c_{\text{tri}}(F)$ and **concentrated crossing number** $cc(F)$ (Allred et al. 2023; Gong et al. 2026).
- **Meridional rank** $\operatorname{mr}(F)$: Minimal meridional generators of $\pi_1(E(F))$; satisfies $\operatorname{mr}(F) \le \min(c_1, c_2, c_3) < b(F)$.
- **2D braid index** $\operatorname{br}(F)$ and minimal chart white-vertex count $w(F)$ (for orientable $F$).
- **Sheet number** $\operatorname{sh}(F)$ (Saito–Satoh) and **1-handle unknotting / stabilization number** $u_1(F)$ (Hosokawa–Kawauchi: minimal internal 1-handles making an orientable $F$ unknotted).
- **Ribbon / fusion number** $\operatorname{rf}(F)$ and **ribbon crossing number** $r\text{-cr}(F)$ (when ribbon).
- **Kirby–Thompson invariants** $L(F), L^\star(F)$ (pants / dual-curve complex length of minimal bridge trisections).
- **Mosaic number** $m(F)$ (Choi–Nelson 2022) and **shadow-complexity** $\operatorname{sc}(F)$ (Ishikawa–Koda–Naoe 2024).

### Tier 3: Deeper 4D Homotopy, Branched-Cover, and Smooth Invariants
Why Tier 3 matters: In $\mathbb{R}^3$, prime knots are determined by $\pi_1(S^3 \setminus K)$ up to mirror (Gordon–Luecke). **In $\mathbb{R}^4$, $\pi_1(S^4 \setminus F)$ does NOT even determine the homotopy type of the complement, let alone $F$!**
- *Example 1*: $8_1^{1,1}$ (spun Hopf link) and $10_1^{1,1}$ (1-twist-spun Hopf link) both have $\pi_1 \cong \mathbb{Z} \oplus \mathbb{Z}$ (same as $6_1^{0,1}$!).
- *Example 2*: Suciu constructed infinitely many distinct ribbon 2-knots all sharing the trefoil group $\pi_1(S^4 \setminus K) \cong \pi_1(S^3 \setminus 3_1)$, separated by $\pi_2$ as a $\mathbb{Z}[\pi_1]$-module!
- *Example 3*: Gordon proved that $k$-twist-spun knots are generally not determined by their exterior alone (Gluck twists can preserve the exterior while changing or preserving the embedding).

Key Tier 3 invariants:
1. **Second Homotopy Module $\pi_2(E(F))$ & Lomonaco's Algebraic 2-Type** $(\pi_1, \pi_2, k \in H^3(\pi_1; \pi_2))$: Computable directly from `ch`-diagrams via Lomonaco's (1981) cellular model of the universal cover! Also the **Farber–Levine torsion pairing** on $H_1(\widetilde{E})$.
2. **Cyclic & Double Branched Covers $\Sigma_d(S^4, F)$**: Closed 4-manifolds whose homology $H_1(\Sigma_2; \mathbb{Z})$, $b_2(\Sigma_2)$, **intersection form** $Q_{\Sigma_2}$ (signature $\sigma$, definiteness, parity), and $\pi_1(\Sigma_2)$ distinguish surfaces with identical $\pi_1(E(F))$ (e.g., $\Sigma_2$ separates $8_1^{1,1}$ from $10_1^{1,1}$).
3. **Structural / Geometric Booleans & Symmetries**:
   - `is_unknotted`, `is_ribbon`, `is_pseudo_ribbon`, `is_spun`, `is_twist_spun` (with `(base_knot, twists)`), `is_fibered` (by Andrews–Sumners, the Artin $0$-spin of a classical knot $K$ is fibered iff $K$ is fibered; by Zeeman's theorem, for every $k \ge 1$ the $k$-twist-spin of *any* classical knot is fibered by punctured 3-manifolds, with $k = 1$ always unknotted).
   - **Full 4-way Symmetry Type** (for oriented surfaces, analogous to KnotInfo): Fully amphicheiral / positive amphicheiral / negative amphicheiral / reversible / chiral & non-invertible, tracking $(F, -F, F^*, -F^*)$ (where $F^*$ is 4D reflection and $-F$ reverses surface orientation).
   - **Gluck Twist Data** (for 2-knots): Whether $S^4_K$ is known diffeomorphic to $S^4$ (automatically true for all ribbon and twist-spun 2-knots) and whether the Gluck-twisted 2-knot $K^*$ is equivalent to $K$.
4. **Categorified / Smooth 4D Invariants**:
   - Classical closed-surface Khovanov cobordism maps $\operatorname{Kh}(\emptyset) \to \operatorname{Kh}(\emptyset)$ vanish or depend only on genus (Tanaka / Rasmussen / Jacobsson), *unless* enhanced via **Morrison–Walker–Wedrich skein lasagna modules** $\mathcal{S}_0^2(S^4 \setminus \nu(F))$ (computable from Kirby 2-handle diagrams of the exterior!), **equivariant / Lee / dotted Bar-Natan evaluations**, or **sliced concordance invariants** ($s, \tau, \Upsilon$ of equatorial cross-sections).

---

## 4. Is `ch`-Index the "Morally Right" Cutoff?

### 4.1 Why `ch`-Index IS Uniquely Powerful for Exhaustive Enumeration
Among *all* known complexity measures on knotted surfaces in $\mathbb{R}^4$, $\operatorname{ch}(F) = \min_D(c(D) + h(D))$ has one decisive virtue that almost every geometric competitor lacks:

> **Uniform Finiteness Theorem (Yoshikawa 1994)**: For every integer $N \ge 0$, there are **only finitely many** ambient isotopy classes of non-splittable surface-links $F \subset \mathbb{R}^4$ (across *all* genera, orientabilities, and component counts simultaneously) with $\operatorname{ch}(F) \le N$.

Compare what goes wrong if you try to tabulate using any other single complexity cutoff $\le N$:
1. **Triple point number $t(F) \le N$**: **Infinite at $N = 0$!** Every ribbon 2-knot, ribbon torus-knot, and spun knot has $t(F) = 0$.
2. **Tri-plane crossing number $c_{\text{tri}}(F) \le N$ alone**: **Infinite at $N = 0$** (all unknotted orientable $\Sigma_g$ have $c_{\text{tri}} = 0$), infinite at $N = 1$ (all unknotted $\mathcal{P}_{n,n}$ have $c_{\text{tri}} = 1$), and even restricting to 2-knots ($S^2$), **it is currently an open problem whether $\{K \cong S^2 : c_{\text{tri}}(K) \le N\}$ is finite for $N \ge 6$** because no upper bound on minimal bridge number $b(K)$ in terms of $c_{\text{tri}}(K)$ is known (Gong, Lewis-Monkman, & Osnes 2026)! Only fixing *both* $(b, c_{\text{tri}})$ gives a finite search space.
3. **Bridge number $b(F) \le N$ alone**: **Infinite at $N = 4$!** Every spun 2-bridge knot, twist-spun 2-bridge knot, and 2-plat ribbon 2-knot has $b(F) = 4$.
4. **2D braid index $\operatorname{br}(F) \le N$ alone**: **Infinite at $N = 3$** (infinitely many ribbon 2-knots have braid index 3), and misses all non-orientable surfaces.
5. **Shadow-complexity $\operatorname{sc}(F) \le N$**: **Infinite already at $N = 1$!** (Ishikawa–Koda–Naoe 2024).

Moreover, Yoshikawa's simple-planar-graph reduction (Lemma 3.9 & Theorem 3.5) bounds the number of checkerboard vertices by $\lfloor (N + 2)/2 \rfloor$:
- Up to $\operatorname{ch} \le 10$: simple planar graphs with **$\le 6$ vertices** and $\le 10$ total weight.
- At $\operatorname{ch} = 11$: **still at most $\lfloor 13/2 \rfloor = 6$ vertices!**
- At $\operatorname{ch} = 12$: at most $7$ vertices.

### 4.2 Why `ch`-Index Is "Morally Flawed" as a Standalone Conceptual Complexity
Despite its algorithmic magic, `ch`-index suffers from **four serious conceptual defects** if treated uncritically like classical crossing number:

1. **Conflation of Intrinsic Topology ($h$) with Knotting ($c$)**:
   Because $\operatorname{ch}(D) = c(D) + h(D)$ weights a Morse saddle ($h$) equally with a projection crossing ($c$), intrinsic genus eats up the complexity budget before any knotting can happen:
   - Unknotted orientable $\Sigma_g$ already costs $\operatorname{ch}(\Sigma_g) = 2g$.
   - Unknotted non-orientable $\#_k \mathbb{R}P^2$ costs $k + |e|/2$ (or $k+1$ if $e=0$).
   - Consequently, at $\operatorname{ch} \le 10$, **not a single knotted genus-2 surface or knotted connected Klein bottle appears**, and only *one* knotted torus ($10_1^1$) barely squeezes in at $\operatorname{ch} = 10$!
2. **Flooding by Multi-Component Links of Unknotted Surfaces**:
   Linking an $S^2$ with an unknotted $T^2$ or Klein bottle requires very few crossings and saddles: **14 out of Yoshikawa's 23 entries (61%) are multi-component links** whose individual components are completely unknotted! By contrast, KnotInfo and LinkInfo strictly separate knots ($c=1$) from links ($c \ge 2$).
3. **Failure of Additivity & Connected-Sum Collapse in Dimension 4**:
   In $\mathbb{R}^3$, prime decomposition is unique (Schubert), and crossing number is conjectured additive. In $\mathbb{R}^4$:
   - **Non-orientable connected sum destroys knotting**: By Viro (1973), $(3\text{-twist-spun } 3_1) \# \mathbb{R}P^2_+ \cong \mathbb{R}P^2_+$, so $\operatorname{ch}(10_3 \# 2_1^{-1}) = 2 < 10 + 2$! Even among unknots, $\operatorname{ch}(\mathbb{R}P^2_+ \# \mathbb{R}P^2_-) = 3 < 2 + 2$.
   - **Orientable 1-handle stabilization vs. trivial connected sum**: Attaching a *nontrivial* internal 1-handle can kill meridian relations and unknot any orientable surface in $\mathbb{R}^4$ (Hosokawa–Kawauchi 1979; Boyle 1988)—for instance, 1-fusion ribbon 2-knots like $8_1$ unknot after a single nontrivial internal 1-handle. By contrast, *local connected sum* $K \# T^2_{\text{std}}$ (attaching a *trivial* 1-handle inside a small 4-ball) preserves $\pi_1(E(K \# T^2_{\text{std}})) \cong \pi_1(E(K))$; however, Baykur–Sunukjian (2016) showed that for exotic smoothly knotted surfaces, even a single *trivial* stabilization $\# T^2_{\text{std}}$ can dissolve smooth knotting! Whether $\operatorname{ch}(F_1 \# F_2) = \operatorname{ch}(F_1) + \operatorname{ch}(F_2)$ holds for all orientable $F_1, F_2$ in $S^4$ remains open (Yoshikawa Question 2.9), which is why handle-reducible surfaces ($F' \# T^2_{\text{std}}$) must be separated from handle-irreducible ones.
   - **Unique prime factorization for 2-knots ($S^2 \subset S^4$) is OPEN**: Unique prime decomposition holds for $n = 1$ (Schubert 1949) and $n \ge 3$ (Milnor 1962), and is open *uniquely* in dimension $n = 2$. Even smooth primality of the unknot $0_1$ is equivalent to the Smooth Unknotting Conjecture: if $K_1 \# K_2 \cong 0_1$, then $\pi_1(E(K_1)) *_{\mathbb{Z}} \pi_1(E(K_2)) \cong \mathbb{Z}$ forces $\pi_1(E(K_i)) \cong \mathbb{Z}$, which implies $K_i$ is *topologically* unknotted (Freedman 1984), whereas *smooth* unknottedness remains open. Because $\operatorname{ch}$-additivity is also unknown, Yoshikawa had to define **"weakly prime"** ($F \neq F_1 \# F_2$ with $0 < \operatorname{ch}(F_i) < \operatorname{ch}(F)$).
4. **Excess Complexity Penalty for Twist-Spinning vs Tri-Plane Crossing Number**:
   Compare $8_1$ (spun trefoil) and $10_2$ (2-twist-spun trefoil): both have **minimal bridge number $b = 4$** and **tri-plane crossing number $\le 6$** ($c_{\text{tri}}(10_2) = 6$ exact by Gong et al. 2026; $c_{\text{tri}}(8_1) \le 6$), yet `ch`-index penalizes $10_2$ up to $\operatorname{ch} = 10$ because realizing triple points ($t(10_2) = 4$) via hyperbolic splittings forces extra crossings in $x_4 = 0$.

### 4.3 Alternate Provably Finite-to-One Complexity Measures

While *single-parameter* geometric invariants ($t(F)$, $b(F)$, $c_{\text{tri}}(F)$, $\operatorname{br}(F)$, $\operatorname{sc}(F)$) fail uniform finiteness, **`ch`-index is NOT the only provably finite-to-one complexity measure on knotted surfaces in $\mathbb{R}^4$**. At least five genuinely distinct finite-to-one families exist, each with dramatically different algorithmic tradeoffs and orderings of "small" knots:

| Complexity Measure | Uniformly Finite-to-One? | Topology vs. Knotting Separation | Primary Algorithmic Bottleneck |
| :--- | :--- | :--- | :--- |
| **1. Marked graph `ch`-index** $c(D) + h(D)$ | **Yes** (all topologies & links) | Conflates genus ($h$) with crossings ($c$) | Unlink testing on $L_\pm(D)$ + equivalence under 8 Yoshikawa moves |
| **2. Composite tri-plane** $\alpha b(F) + c_{\text{tri}}(F)$ (or $b + cc$) | **Yes** (for any $\alpha > 0$, all topologies & links) | **Clean**: $c_{\text{tri}} = 0$ for all unknotted orientable $\Sigma_g$ ($b = 2g+1$) | Larger raw branching of 3 ordered $2b$-endpoint tangles (mitigated if $c_{\text{tri}} = cc$) |
| **3. Complement ideal pentachora** $f_4(S^4 \setminus F)$ (*Regina*) | **Yes on 2-knots** ($\le 2$-to-1 per exterior by Gluck 1962) | **Intrinsic to complement**; Cappell–Shaneson knots have $f_4 = 2$, $8_1$ has $f_4 = 6$! | **Smooth $S^4$-filling recognition** & 4D Pachner-graph explosion at $f_4 \ge 7$ |
| **4. Affine / PL vertex count** $v_{\text{PL}}(F) \subset \mathbb{R}^4$ (4D stick #) | **Yes** (semi-algebraic finiteness, Thom–Milnor) | Penalizes genus via triangulation vertices | Local flatness (unknotted vertex links) + doubly-exponential CAD in $\mathbb{R}^{4V}$ |
| **5. Full 2D shadow / chart+hoops / grid mosaic $m(F)$** | **Yes** once *all* $\mathbb{R}^2$ crossings/hoops are counted | Mosaic/braid discretize `ch`; charts are orientable-only | Global $\mathbb{R}^4$ liftability (for 2D shadows) or $16^{m^2}$ grid blowup (mosaics) |

Key mathematical details on each alternative:

1. **Composite Bridge-Trisection Complexity $\alpha b(F) + c_{\text{tri}}(F)$**:
   - **Why finite**: Fixing $b \le B$ ($2b$ boundary points) and $c_1 + c_2 + c_3 \le C$ restricts each $D_i$ to finitely many planar 4-valent tangle diagrams without closed components. Checking pairwise unlink admissibility $(D_i \cup \overline{D}_j)$ has the exact same complexity as checking $L_\pm(D)$ in `ch`-diagrams.
   - **Conceptual advantage**: Unknotted orientable surfaces of *any* genus $g$ satisfy $c_{\text{tri}}(\Sigma_g) = 0$ and $b(\Sigma_g) = 2g + 1$, cleanly decoupling genus from knotting, and placing $8_1$ and $10_2$ at the same complexity $(b, c_{\text{tri}}) = (4, 6)$.
   - **Subtlety with pants distance**: Replacing $c_{\text{tri}}(F)$ by the Kirby–Thompson pants distance $L(F)$—i.e., $b(F) + L(F)$—is **not** manifestly finite-to-one, because vertices of the pants graph $\mathcal{P}(\Sigma_{0, 2b})$ have **infinite degree** (an $A$-move allows arbitrarily many Dehn twists around curves in a pants decomposition).
2. **Ideal Pentachoron Number $f_4(S^4 \setminus F)$ (3D SnapPy Analogue in *Regina*)**:
   - In 3D, SnapPy censuses hyperbolic knots by ideal tetrahedron count of $S^3 \setminus K$ (justified by Gordon–Luecke). In 4D, one enumerates $n$-pentachoron ideal triangulations of $E = S^4 \setminus \nu(F)$:
     - **Finiteness on 2-knots**: There are finitely many $n$-pentachoron gluings with $\partial E \cong S^2 \times S^1$, and by **Gluck's Theorem (1962)** ($\pi_0 \operatorname{Diff}(S^2 \times S^1) \cong (\mathbb{Z}/2)^3$), a given exterior $E$ admits **at most two** smooth fillings to $S^4$ (trivial vs. Gluck twist). *(For $g \ge 1$, e.g. $\partial E \cong T^3$, $\operatorname{MCG}(T^3) \cong \operatorname{GL}(3, \mathbb{Z})$ is infinite, though requiring the filling to have $H_*(S^4)$ restricts slopes.)*
     - **Wild reordering of small 2-knots**: **Budney–Burton–Hillman (2012)** triangulated nontrivial **Cappell–Shaneson 2-knot** complements using just **$f_4 = 2$ pentachora** (despite enormous unknown `ch`-index!), whereas **Burke, Burton, Ray, Spreer, Thompson, & Zymaris (2026)** proved the spun trefoil $8_1$ has **exact minimum $f_4(S^4 \setminus 8_1) = 6$**, with upper bounds $f_4 \le 10$ for spun $5_1$ (`5a2`), $f_4 \le 12$ for spun `7a7` / `8n3`, and $f_4 \le 16$ for spun $4_1$ ($10_1$).
     - **Why it cannot be our primary enumerator yet**: Filtering arbitrary 4-manifold triangulations requires deciding whether Dehn-filling $\partial E$ yields *smooth $S^4$*. Even closed orientable 4-manifold triangulations have only been censused up to $f_4 \le 6$ (Burke–Burton–Spreer 2026), because simplifying $S^4$ via 4D Pachner moves can require $10^5$–$10^6$ intermediate pentachora! However, $f_4(S^4 \setminus F)$ is an extraordinary **invariant column** to compute via *Regina* for entries in our table.
3. **Intrinsic PL / Affine Polyhedral Vertex Count $v_{\text{PL}}(F)$ in $\mathbb{R}^4$**:
   - Fixing $V$ vertices in $\mathbb{R}^4$ leaves finitely many abstract surface triangulations $K$; the locus of locally flat, non-self-intersecting affine realizations $\mathcal{R}(K) \subset \mathbb{R}^{4V}$ is semi-algebraic (defined by simplex disjointness inequalities once vertex links are unknotted), hence has **finitely many path-components** (ambient PL isotopy classes) by Oleinik–Petrovsky–Thom–Milnor. Enumerating path-components via cylindrical algebraic decomposition in $4V \ge 32$ real variables is computationally infeasible compared to discrete diagrams.
4. **Full 2D Shadow Complexity & Grid/Braid Discretizations**:
   - Generic projections $F \to \mathbb{R}^2$ have smooth fold curves, cusps, and transverse/tangential double-curve crossings; bounding the total number of 0D singularities in $\mathbb{R}^2$ is finite-to-one (unlike 3D triple + branch points), as is bounding $m + w + b + \text{edges/hoops}$ in Kamada charts (without edge/hoop bounds, concentric hoops between free edges encode arbitrary braids of $w = b = 0$). However, deciding global $\mathbb{R}^4$-liftability of a 2D fold/crossing shadow is far harder than unlink testing.
   - **Marked graph mosaics** $m(F)$ (Choi–Nelson 2022) and **surface singular braid** `strands + length` (Jabłonowski 2013) are provably finite-to-one discretizations of `ch`-diagrams, trading Yoshikawa's compact planar graphs ($\le 6$–$7$ vertices) for coarser grid/word search spaces.

### 4.4 Moral Verdict: How We Should Structure Our Expanded Table

Rather than abandoning `ch`-index or trusting a flat `ch`-ordered list blindly, the "morally right" architecture combines **topological stratification** with **dual `(ch, tri-plane)` indexing**:

1. **Stratify Primarily by Topological Type (KnotInfo vs LinkInfo separation)**:
   - **Sub-table A (2-Knots)**: Connected $S^2 \hookrightarrow \mathbb{R}^4$ ($0_1, 8_1, 9_1, 10_1, 10_2, 10_3, \dots$).
   - **Sub-table B (Orientable Genus-$g$ Surface-Knots)**: Connected $T^2, \Sigma_g \hookrightarrow \mathbb{R}^4$, restricted to *handle-irreducible* knots ($F \neq F' \# T^2_{\text{std}}$) ($2_1^1, 10_1^1, \dots$).
   - **Sub-table C (Non-Orientable Surface-Knots)**: Connected $\#_k \mathbb{R}P^2 \hookrightarrow \mathbb{R}^4$, restricted to *$\mathbb{R}P^2$-irreducible* knots.
   - **Sub-table D (Multi-Component Surface-Links)**: Non-splittable $c \ge 2$ surface-links, indexed by component topology vector.
2. **Record Bi-Graded `(c(D), h(D))` Complexity & Excess `ch`-Index**:
   - Define **excess `ch`-index** $\widetilde{\operatorname{ch}}(F) = \operatorname{ch}(F) - \operatorname{ch}_{\text{unknot}}(\text{top}(F), e(F))$, measuring *pure knotting/linking complexity* above the topological floor.
   - Within each topological stratum (e.g. $S^2$-knots, where $h(D) \ge 2$), enumerating by $\operatorname{ch}(F) = 11, 12, \dots$ via simple planar graphs with $\le 6$–$7$ vertices is finite, canonical, and practical.
3. **Run a Parallel `(b, c_tri)` Census & Record Complement Pentachora $f_4$**:
   - Enumerating 4-bridge tri-plane diagrams $(4; 2, 2, 2)$ by $c_{\text{tri}} = 6, 7, 8, \dots$ provides an independent finite complexity filter that puts twist-spun / non-ribbon 2-knots on equal footing with ribbon 2-knots, while recording *Regina* complement triangulation bounds $f_4(S^4 \setminus F)$ and parametric families (2-plat, twist-spun, Satoh tubes).

---

## 5. Key Open Conjectures

### 5.1 Foundational ("Dream") Conjectures — Ranked by Tractability

1. **Bridge-Boundedness by Tri-Plane Crossing Number *(Most Tractable to Prove)*** *(Gong–Lewis-Monkman–Osnes 2026)*:
   - *Statement*: Does $c_{\text{tri}}(F) \le C$ bound minimal bridge number $b(F) \le f(C, \chi(F))$ for connected (or prime) surfaces? **If YES**, pure $c_{\text{tri}}$ at fixed topology is finite-to-one *on its own*!
   - *Why a concrete proof path exists (2D/3D combinatorics, zero 4D gauge-theory obstructions)*:
     - **Base case $C = 0$ is already proved** (Meier–Zupan 2017; diagrammatic proof in Allred et al. 2023): *every* crossingless tri-plane diagram, regardless of how huge $b$ is, destabilizes down to $b = 2g + 1$.
     - **Immediate meridional-rank bound**: Gluing the $3b$ tangle intervals at their $2b$ shared bridge points $P$ yields a connected 3-regular multigraph $\Gamma = D_1 \cup_P D_2 \cup_P D_3$ (cycle rank $b + 1 \ge 2$). From a basepoint directly above the bridge sphere $\Sigma$, Wirtinger meridians are continuous across $P$ and jump *only* at the $C$ under-crossings; cutting $\Gamma$ at $C$ interior edge points yields at most $C + 1$ pieces (at most $C$ if $C \ge 1$), proving **$\operatorname{mr}(F) \le c_{\text{tri}}(F) + 1$** (and $\operatorname{mr}(F) \le c_{\text{tri}}(F)$ for $c_{\text{tri}} \ge 1$). This already rules out all high-meridional-rank families (e.g. $\#_n 8_1$).
     - **Normal-curve / Catalan band-weight reduction**: At most $2C$ strands contain a crossing, marking $|A| \le 4C$ endpoints on $\partial D^2$ and partitioning the remaining $\ge 3b - 2C$ crossingless chords into $O(C)$ parallel planar bands. If a band of width $w$ has cycle lengths $(2p, 2q, 2r)$ in $(U_{12}, U_{23}, U_{31})$, its Euler characteristic contribution is $(\frac{1}{p} + \frac{1}{q} + \frac{1}{r} - 1)w$. Hyperbolic triples ($\frac{1}{p}+\frac{1}{q}+\frac{1}{r} < 1$) force $g \to \infty$, while $(1, \infty, \infty)$-bands give parallel 1-bridge loops; the entire problem reduces to destabilizing nested $(1,\infty,\infty)$ rainbow bands and ruling out Euclidean wallpaper weaves $(p,q,r) \in \{(3,3,3), (2,4,4), (2,3,6)\}$.
2. **Orientable Complexity Additivity *(Most Tractable to Disprove — or Prove for $b$ via Thin Position)***:
   - *Statement* (**Yoshikawa Q2.9**; **Allred–Zupan Q6.1**): For orientable $F_1, F_2$, is $\operatorname{ch}(F_1 \# F_2) = \operatorname{ch}(F_1) + \operatorname{ch}(F_2)$, $c_{\text{tri}}(F_1 \# F_2) = c_{\text{tri}}(F_1) + c_{\text{tri}}(F_2)$, and $b(F_1 \# F_2) = b(F_1) + b(F_2) - 1$?
   - *Two concrete paths*: For $b$, intersecting a decomposing $S^3$ with a bridge trisection (weak reducibility / Scharlemann-style cut-and-paste) gives a plausible route toward $(-1)$-additivity—*or* toward a **counterexample**: **Joseph–Pongtanapaisan (2023)** proved that meridional rank $\operatorname{mr}$ collapses under connected sum of 2-knots, so either the 4D Meridional Rank Conjecture fails or $b$ fails $(-1)$-additivity! Searching for subadditivity of `ch` or `c_tri` on opposite twist-spins is a finite computer search.
3. **Unknotting & Gluck Twists *(Tractable in TOP for $g=1,2$ and for Smooth Counterexamples; Mysterious for General Smooth Proofs)***:
   - *Topological $g \in \{1, 2\}$*: Conway–Powell (2023) already proved $\pi_1 \cong \mathbb{Z} \implies$ topologically unknotted for $g \ge 3$ via modified surgery; closing $g = 1, 2$ is a concrete unstable $\mathbb{Z}[\mathbb{Z}]$-surgery cancellation problem.
   - *Smooth counterexample frontier*: **Miyazawa (2023)** (extended by **Kang–Park–Taniguchi 2025**, **Baraglia 2026**) used **Real Seiberg–Witten theory** to construct infinitely many **exotic smooth $\mathbb{R}P^2 \hookrightarrow S^4$** (topologically unknotted with $\pi_1 \cong \mathbb{Z}/2$, smoothly knotted, with branched covers tied to Gluck twists!). Finding an analogous exotic $S^2 \subset S^4$ with $\pi_1 \cong \mathbb{Z}$ via skein lasagna (Ren–Willis) or real Floer invariants would *disprove* Smooth Unknotting—whereas *proving* Smooth Unknotting or $S^4_K \cong S^4$ for all 2-knots sits at the Smooth 4D Poincaré wall.
4. **Unique Prime Factorization of Smooth (or TOP) 2-Knots *(Deceptively Blocked Geometrically)***:
   - Known for $n = 1$ (Schubert) and $n \ge 3$ (Milnor, via high-dimensional engulfing/$h$-cobordism), open uniquely for $n = 2$. Two decomposing 3-spheres $S^3_A, S^3_B \subset S^4$ intersect in a closed **2-manifold $\Sigma_g$** ($g > 0$!) rather than circles, which cannot always be compressed away from the slice knot $S^3_B \cap K$, and 4-manifold exteriors lack a JSJ decomposition (nor does $G_1 *_{\mathbb{Z}} G_2$ enjoy unique group-theoretic amalgamation over $\mathbb{Z}$).
5. **Computable Global Bounds on Yoshikawa Move Sequences *(Likely Andrews–Curtis / Decidability Blocked)***:
   - Unlike 3D Reidemeister bounds (Lackenby), a banded unlink diagram $(L, B)$ directly encodes 2-handle slides over 1-handles (**Andrews–Curtis moves**) plus 2-handle stabilizations ($\Omega_6, \Omega_7$). Any computable upper bound on peak `ch`-index for unknotting would immediately decide Smooth 2-Knot Unknotting and brushes up against smooth $S^4$ recognition.


### 5.2 Concrete Low-Complexity Conjectures Our Project Can Directly Attack
1. **Exact $c_{\text{tri}}(8_1) = 6$ & the $c_{\text{tri}} \le 6$ 2-Knot Census** *(Gong–Lewis-Monkman–Osnes 2026, Conjectures 3.6–3.7)*: That $c_{\text{tri}}(8_1) = 6$, and that $\{0_1, 8_1, 10_2\}$ are the **only** 2-knots with $c_{\text{tri}}(K) \le 6$ (directly testable at $b = 4, 5$ by exhaustive tri-plane enumeration!).
2. **Concentrated Crossing Equality $c_{\text{tri}}(F) = cc(F)$ & Simultaneous $(b, c_{\text{tri}})$ Minimization** *(Allred–Zupan 2023, Questions 6.2–6.3)*: Can all crossings of a minimal tri-plane diagram always be concentrated into a single tangle $D_3$ (leaving $D_1, D_2$ as crossingless Catalan matchings)? If true on low-bridge diagrams, tri-plane tabulation accelerates by orders of magnitude! Also tests whether $10_2^{0,-2}$ fails to minimize $b$ and $c_{\text{tri}}$ simultaneously.
3. **Smallest Knotted $\mathbb{R}P^2$-Irreducible Projective Plane** *(Yoshikawa 1994 §5; Allred–Zupan 2023 Q6.5)*: None exists with $\operatorname{ch} \le 10$! Finding the exact minimum $\operatorname{ch}$ (at $\operatorname{ch} = 11, 12, 13$?) and $c_{\text{tri}}$ of a knotted $\mathbb{R}P^2$ (e.g., Kinoshita's 2-twist-spun trefoil sum or twisted $\mathbb{R}P^2$'s) is a prime target for our expansion.
4. **Open Invariants Inside Yoshikawa's Table**:
   - Is $10_2^{0,-2}$ non-ribbon with $t(10_2^{0,-2}) > 0$? (Attackable via computer search over finite **symmetric quandle 3-cocycles**, Kamada–Oshiro 2010.)
   - Is the pseudo-ribbon link $10_1^{1,1}$ ($t = 0$) genuinely non-ribbon, and what is the exact value of $t(10_1^{-1,-1}) \in [2, 12]$?

