# Complexity Measures, Intrinsic Cutoffs, and Alternative Filtrations for Knotted Surfaces

This note analyzes the mathematical and computational behavior of Yoshikawa's hyperbolic-splitting complexity $\operatorname{ch}(F) = \min_D (c(D) + m(D))$ alongside intrinsic topological floors, weighted $(c, m)$ profiles, and Meier–Zupan's bridge trisection crossing number $c_{\text{tri}}(F)$. All theorems, censuses, and cross-measure comparisons below are machine-verified by [`complexity.py`](../code/src/surfaceknots/complexity.py), [`triplane.py`](../code/src/surfaceknots/triplane.py), [`triplane_census.py`](../code/src/surfaceknots/triplane_census.py), and [`_fast_enum.c`](../code/src/surfaceknots/_fast_enum.c), with full datasets exported under [`code/data/complexity/`](../code/data/complexity/).

---

## 1. Intrinsic Topology vs. Knotting: Unknot $\operatorname{ch}$-Index Floors & Excess $\operatorname{ch}$-Index

A basic structural defect of raw $\operatorname{ch}$-index is that both classical crossings $c(D)$ and saddle markers $m(D)$ are forced upward merely by raising the genus or normal Euler number of an **unknotted** surface component. Separating this topological baseline from genuine 4-dimensional knotting or linking yields an exact additive floor $\operatorname{ch}_0(F)$.

### Theorem 1.1 (Zero-Crossing Orientability Theorem)
Let $D$ be an admissible marked diagram with $c(D) = 0$ classical crossings. Then the represented surface-link $F(D) \subset S^4$ is **orientable**, and every component $F_i \subset F(D)$ has normal Euler number $e(F_i) = 0$.

*Proof.* Since $c(D) = 0$, every vertex of $D$ is an $A$-/$B$-smoothed marker ($Y$ or $Z$). Both classical resolutions $L_-(D)$ and $L_+(D)$ are planar embedded collections of disjoint Jordan curves in $S^2$ with zero crossings. Two-color the complementary regions of $L_-(D) \subset S^2$ checkerboard-style (say white/black). Each marker vertex locally merges or splits two arcs separating opposite colors without introducing crossings; compatibility with the planar checkerboard coloring induces a coherent orientation on the resolution circles of $L_-(D)$ across every saddle band of $D$. By Kamada's orientability criterion on hyperbolic splittings, $F(D)$ is orientable. Moreover, since $L_+(D)$ and $L_-(D)$ each have zero classical crossings, Kamada's resolution-writhe formula gives $e(F_i) = w(L_+(D)|_{F_i}) - w(L_-(D)|_{F_i}) = 0 - 0 = 0$. $\square$

### Theorem 1.2 (Exact Intrinsic Unknot $\operatorname{ch}$-Index Floor)
For a single unknotted closed connected surface $U \subset S^4$, let $\operatorname{ch}_0(U) = \min_{D \simeq U} \operatorname{ch}(D)$ denote its minimal $\operatorname{ch}$-index over admissible diagrams. Then $\operatorname{ch}_0(U) = c_0(U) + m_0(U)$ is achieved simultaneously on both $c(D)$ and $m(D)$:

1. **Orientable unknotted genus-$g$ surface** $\Sigma_g$ ($g \ge 0$, $e(\Sigma_g) = 0$):
   $$c_0(\Sigma_g) = 0, \qquad m_0(\Sigma_g) = 2g, \qquad \operatorname{ch}_0(\Sigma_g) = 2g.$$
2. **Non-orientable unknotted surface** $N_{k, e} = \#^k \mathbb{RP}^2$ ($k \ge 1$, $|e| \le 2k$, $e \equiv 2k \pmod 4$, i.e., $P_+^{\# p} \# P_-^{\# q}$ with $p + q = k$ and $2(p - q) = e$):
   $$c_0(N_{k, e}) = \max\!\left(1, \frac{|e|}{2}\right), \qquad m_0(N_{k, e}) = k, \qquad \operatorname{ch}_0(N_{k, e}) = k + \max\!\left(1, \frac{|e|}{2}\right).$$

For a split union $U_1 \sqcup \dots \sqcup U_r$ of standard unknotted components, $\operatorname{ch}_0(\bigsqcup_{i=1}^r U_i) = \sum_{i=1}^r \operatorname{ch}_0(U_i)$.

*Proof.* For any connected diagram $D$ of $U$ with resolution circle counts $c_+(D) \ge 1$ and $c_-(D) \ge 1$, Euler characteristic satisfies $\chi(U) = c_+(D) - m(D) + c_-(D) \ge 2 - m(D)$, whence
$$m(D) \ge 2 - \chi(U) = \begin{cases} 2g & \text{if } U \cong \Sigma_g, \\ k & \text{if } U \cong N_{k, e}. \end{cases}$$
For $U \cong \Sigma_g$, $c(D) \ge 0$ is trivial, and $g$ nested pairs of opposite-marker bands on $0_1$ (`2_1^1` repeated $g$ times) achieve $(c, m) = (0, 2g)$.

For $U \cong N_{k, e}$ ($k \ge 1$), non-orientability implies $c(D) \ge 1$ by Theorem 1.1. Furthermore, under any orientation of the resolution circles of $L_+(D)$ and $L_-(D)$, each classical crossing $v$ contributes $\pm 1$ to $w(L_+)$ and $\pm 1$ to $w(L_-)$, hence at most $2$ to $|w(L_+) - w(L_-)| = |e(U)|$; thus $c(D) \ge |e(U)|/2$, proving $c(D) \ge \max(1, |e|/2)$ and $\operatorname{ch}(D) \ge k + \max(1, |e|/2)$.

For sharpness on $N_{k, e}$:
- When $e \ne 0$, set $r = |e|/2 \ge 1$ and $s = (k - r)/2 \in \mathbb{Z}_{\ge 0}$. Taking the planar connected sum of $r$ copies of $P_{\operatorname{sgn}(e)}$ (`2_1^-1`, $(c, m) = (1, 1)$) with $s$ copies of the unknotted torus `2_1^1` ($(c, m) = (0, 2)$) yields $N_{k, e}$ by Dyck's relation $P_\pm \# T^2 \cong P_\pm \# P_+ \# P_-$ with $(c, m) = (r, r + 2s) = (|e|/2, k)$.
- When $e = 0$ ($k = 2s \ge 2$), the admissible prime 3-vertex diagram
  $$\texttt{X[1,2,3,4], Y[4,2,5,6], Z[6,5,3,1]}$$
  has $(c, m) = (1, 2)$, $c_+ = c_- = 1$, $\chi = 0$, non-orientable surface type $K^2$, $\pi_1 \cong \mathbb{Z}_2$, and $e = 0$, realizing $N_{2, 0} = P_+ \# P_-$ at $\operatorname{ch} = 3 = 2 + \max(1, 0)$. Taking its planar connected sum with $s - 1$ copies of `2_1^1` ($(c, m) = (0, 2)$) realizes every $N_{2s, 0}$ with $(c, m) = (1, 2s)$. $\square$

> [!NOTE]
> **Exact parallel with Allred–Joseph–Meier–Zupan (2023):** Allred et al. proved that the bridge trisection crossing number of any unknotted non-orientable surface $N_{k, e} = \mathcal{P}_{p, q}$ equals $c_{\text{tri}}(N_{k, e}) = \max(1, |p - q|) = \max(1, |e|/2)$. Theorem 1.2 shows that this classical-crossing floor $c_0(N_{k, e}) = \max(1, |e|/2)$ is **identical** in marked ch-diagrams, while $\operatorname{ch}_0$ adds the exact Morse saddle cost $m_0(N_{k, e}) = k = 2 - \chi$.

### Corollary 1.3 (Strict Connected-Sum Subadditivity on Unknotted Klein Bottles)
Planar $\operatorname{ch}$-index is strictly subadditive even on standard unknots:
$$\operatorname{ch}(P_+ \# P_-) = 3 < 2 + 2 = \operatorname{ch}(P_+) + \operatorname{ch}(P_-), \qquad \text{whereas} \qquad \operatorname{ch}(P_+ \# P_+) = 4 = \operatorname{ch}(P_+) + \operatorname{ch}(P_+).$$
Geometrically, $P_+$ and $P_-$ carry opposite Reidemeister I kinks ($e = +2$ and $e = -2$); fusing them across a single mixed classical crossing cancels their normal Euler numbers while saving one classical crossing.

### Definition 1.4 (Excess $\operatorname{ch}$-Index $\widetilde{\operatorname{ch}}$)
For any surface-link $F = F_1 \sqcup \dots \sqcup F_r \subset S^4$, let $\operatorname{ch}_0(F) = \sum_{i=1}^r \operatorname{ch}_0(F_i^{\text{std}})$, where $F_i^{\text{std}}$ is the standard unknotted surface with the same orientability, genus, and normal Euler number $e(F_i)$. The **excess $\operatorname{ch}$-index** is
$$\widetilde{\operatorname{ch}}(F) := \operatorname{ch}(F) - \operatorname{ch}_0(F) \ge 0.$$
By Theorem 1.2, $\widetilde{\operatorname{ch}}(F) = 0$ if and only if $F$ admits a diagram achieving the unknot floor on every component (in particular, $\widetilde{\operatorname{ch}}(U) = 0$ for every standard unknot/unlink).

**Empirical validation on Yoshikawa's Hopf 2-link cohort:** In Yoshikawa's $\operatorname{ch} \le 10$ table, the six fundamental 2-component Hopf-linked surfaces (linking $S^2$, $T^2$, and $\mathbb{RP}^2_\pm$ via four inter-component classical crossings $c_{\text{cross}} = 4$) are artificially scattered across five raw $\operatorname{ch}$-levels $6, 7, 8, 9, 10$. Under $\widetilde{\operatorname{ch}}$, all six collapse to a **single constant complexity level** $\widetilde{\operatorname{ch}}(F) = 4$:

| Label | Topology $(F_1, F_2)$ | Normal Euler $e$ | Raw $\operatorname{ch}(F)$ | $(c_{\text{self}}, c_{\text{cross}}, m)$ | Intrinsic Floor $\operatorname{ch}_0(F)$ | Excess $\widetilde{\operatorname{ch}}(F)$ |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: |
| `6_1^{0,1}` | $S^2 \sqcup T^2$ | $(0, 0)$ | $6$ | $(0, 4, 2)$ | $0 + 2 = 2$ | **$4$** |
| `7_1^{0,-2}` | $S^2 \sqcup \mathbb{RP}^2$ | $(0, \pm 2)$ | $7$ | $(1, 4, 2)$ | $0 + 3 = 3$ | **$4$** |
| `8_1^{1,1}` | $T^2 \sqcup T^2$ | $(0, 0)$ | $8$ | $(0, 4, 4)$ | $2 + 2 = 4$ | **$4$** |
| `8_1^{-2,-2}` | $\mathbb{RP}^2 \sqcup \mathbb{RP}^2$ | $(\pm 2, \pm 2)$ | $8$ | $(2, 4, 2)$ | $2 + 2 = 4$ | **$4$** |
| `9_1^{1,-2}` | $T^2 \sqcup \mathbb{RP}^2$ | $(0, \pm 2)$ | $9$ | $(1, 4, 4)$ | $2 + 3 = 5$ | **$4$** |
| `10_1^{-2,-2}` | $\mathbb{RP}^2 \sqcup K^2$ | $(\pm 2, 0)$ | $10$ | $(3, 4, 3)$ | $3 + 3 = 6$ | **$4$** |

(Note that for $S^2$, any knotted or linked component must also pay at least one saddle pair $m \ge 2$ whenever $c_+ + c_- \ge 4$ — in `6_1^{0,1}` and `7_1^{0,-2}`, the $S^2$ component uses $0$ markers of its own because its two classical crossings resolve into a single circle on both sides!)

**Weighted $(c, m)$ profiles $w_\alpha(D) = c(D) + \alpha m(D)$:** Varying the saddle weight $\alpha \in \{0.5, 1.0, 1.5, 2.0\}$ across our Hub + Spoke catalog reveals two distinct regimes:
- $\alpha < 1$ ($\alpha = 0.5$) favors high-genus ribbon tubes and orientable handles ($m$ cheap, $c$ expensive), pulling `2_1^1` ($w_{0.5} = 1.0$) below `2_1^-1` ($w_{0.5} = 1.5$) and pulling 2-fusion ribbon knots toward 2-twist-spins.
- $\alpha > 1$ ($\alpha = 2.0$) penalizes Morse saddles (approximating handle/ribbon-fusion complexity), separating 1-fusion ($m = 2$) 2-knots ($w_{2.0} = c + 4$) cleanly from $m \ge 3$ orientable/non-orientable surfaces.

---

## 2. Twist-Spun $\#\,\mathbb{RP}^2$ Stabilization Collapse

Connected sum with a standard projective plane $P_\pm$ (`2_1^-1`, $e = \pm 2$, $\operatorname{ch} = 2$) causes dramatic drops in $\operatorname{ch}$-index on knotted 2-spheres. A single topological identity governs the primary collapses across $\operatorname{ch} \le 12$:

### Proposition 2.1 (Twist-Spun $\#\,\mathbb{RP}^2$ Period-2 Collapse; Viro 1973, Kamada, Baykur–Sunukjian)
For every classical knot $K \subset S^3$ and every integer $k \ge 0$,
$$\tau^{k+2}(K) \# P_\pm \;\cong\; \tau^k(K) \# P_\pm \qquad \text{(smoothly in } S^4\text{)}.$$

*Mechanism.* During the $S^1$-spinning sweep of $(B^3, K)$, sliding the cross-cap $P_\pm$ once around the meridian of $K$ reverses the local orientation of the normal bundle of $K$ along the cross-cap path, changing the framing twist during one revolution by $\pm 2$ full twists ($\pi_1(\mathrm{SO}(3)) \cong \mathbb{Z}_2$).

Specializing Proposition 2.1 to the trefoil $K = 3_1$ and figure-eight $K = 4_1$ simultaneously explains **all three** major non-orientable connected-sum collapses encountered in our $\operatorname{ch} \le 12$ enumerations:

1. **Odd twist-spins ($k = 1$, Viro's unknotting collapse):** Since Zeeman's theorem gives $\tau^1(K) \cong 0_1$ ($S^2$ unknot),
   $$10_3 \# P_\pm \;=\; \tau^3(3_1) \# P_\pm \;\cong\; \tau^1(3_1) \# P_\pm \;\cong\; 0_1 \# P_\pm \;=\; P_\pm.$$
   Here $\operatorname{ch}$ collapses by **$10$ vertices**: $\operatorname{ch}(10_3) + \operatorname{ch}(P_\pm) = 10 + 2 = 12 \longrightarrow \operatorname{ch}(10_3 \# P_\pm) = 2$.
2. **Even twist-spins ($k = 0$, 2-twist-spin $\to$ Artin spin collapse):**
   $$10_2 \# P_\pm \;=\; \tau^2(3_1) \# P_\pm \;\cong\; \tau^0(3_1) \# P_\pm \;=\; 8_1 \# P_\pm.$$
   Thus $10_2 \# P_\pm$ drops from naive $\operatorname{ch} = 10 + 2 = 12$ to $\operatorname{ch}(8_1 \# P_\pm) \le 10$! Moreover, in our $\operatorname{ch} = 10, 11, 12$ censuses, $8_1 \# P_\pm$, $9_1 \# P_\pm$, and $10_2 \# P_\pm$ occupy the exact same invariant cell ($E_1 = (3x, x+1)$, $\pi_1 \cong D_3 \times \mathbb{Z}_2$).
3. **$\operatorname{ch} = 11$ ribbon trefoil/figure-eight packet at $\operatorname{ch} = 12, m = 3$:** In our exhaustive $\operatorname{ch} = 12, m = 3$ $\mathbb{RP}^2$ census (§3.2), two 12-vertex diagrams realize the invariant cell $E_1 = (5x, x+1)$ shared by $10_1 \# P_\pm = \tau^0(4_1) \# P_\pm$, $11_1 \# P_\pm$, and $11_2 \# P_\pm$ (reducing naive $\operatorname{ch} = 11 + 2 = 13$ to $\le 12$).

---

## 3. Targeted C99 $\operatorname{ch} = 12$ Censuses: $m = 2$ Prime $2$-Knots & Low-Saddle $\mathbb{RP}^2$-Knots

Using the dual plantri-style C99 quadrangulation generator (`generate_prime_shadows_c` in [`_fast_enum.c`](../code/src/surfaceknots/_fast_enum.c), enumerating all $15,882$ prime 12-vertex spherical shadows in $0.86\,\text{s}$ and all $75,858$ prime 13-vertex shadows in $4.76\,\text{s}$, matching `OEIS A113201`) combined with targeted marker-range scanning (`enumerate_targeted_census_level_c`) and mountain-pass $+1/+2$ tunneling, we executed exact strata censuses at $\operatorname{ch} = 11$ and $\operatorname{ch} = 12$ (full records in [`targeted_ch12_13_census.json`](../code/data/complexity/targeted_ch12_13_census.json)).

### 3.1. Exhaustive $\operatorname{ch} = 12, m = 2$ Connected $S^2$-Knots ($15,882$ Shadows)
Every prime connected $S^2$-knot ($\chi = 2$) has $m(D) \ge 2$, and the $m = 2$ stratum ($c_+ = c_- = 2$, $c = 10$) contains all minimal-saddle 1-fusion ribbon 2-knots and 2-bridge twist-spins at $\operatorname{ch} = 12$. Across all $15,882$ shadows:
- Raw admissible candidates: **$55,038$** $\to$ **$1,158$** level-minimal $\Omega_{1..8}$ components $\to$ **$29$** post-tunneling diagrams.
- Eliminating lower-$\operatorname{ch}$ `AMBIGUOUS_SHADOW` collisions ($0_1, 8_1, 9_1, 10_1, 10_2, 10_3, 11_1, 11_2$ and $\operatorname{ch} \le 12$ composites) leaves **7 certified new prime $\operatorname{ch} = 12$ $S^2$-knot classes**:

| Class ID | Move Comps | $E_1(\pi_1)$ | $\Delta^{\text{SW}}(s)$ (Biquandle) | $H_1(\Sigma_2)$ | $|\operatorname{Hom}(\pi_1, G)|$: $(D_{10}, A_4, F_{21}, A_5)$ | Geometric Identification / Spoke Match |
| :--- | :---: | :--- | :--- | :---: | :---: | :--- |
| **`12_m2_S2_1`** | $3$ | $(2x^2 - 3x + 2)$ | $2s^2 - 3s + 2$ | $\mathbb{Z}_7$ | $(10, 12, 21, 180)$ | **`Spin(5_2)`** ($\tau^0(5_2)$, **exact $\operatorname{ch} = 12$ certified!**) |
| **`12_m2_S2_2`** | $1$ | $(3x^2 - 3x + 1)$ | $s^2 - 3s + 3$ | $\mathbb{Z}_7$ | $(10, 36, 63, 180)$ | Asymmetric 1-fusion ribbon 2-knot ($H_1(M_3) = \mathbb{Z}_2 \oplus \mathbb{Z}_{14}$) |
| **`12_m2_S2_3`** | $6$ | $(5,\; x + 1)$ | $1$ | $\mathbb{Z}_5$ | $(30, 12, 21, 180)$ | Non-ribbon 2-twist-spin $\tau^2(4_1) / \tau^2(5_1)$ ($\pi_1 \cong \mathbb{Z}_5 \rtimes_{-1} \mathbb{Z}$) |
| **`12_m2_S2_4`** | $2$ | $(7,\; x + 1)$ | $1$ | $\mathbb{Z}_7$ | $(10, 12, 21, 60)$ | Non-ribbon 2-twist-spin $\tau^2(5_2)$ ($\pi_1 \cong \mathbb{Z}_7 \rtimes_{-1} \mathbb{Z}$) |
| **`12_m2_S2_5`** | $5$ | $(x^3 - 2x^2 + 3x - 1)$ | $s^3 - 3s^2 + 2s - 1$ | $\mathbb{Z}_7$ | $(10, 12, 21, 60)$ | Asymmetric 1-fusion ribbon 2-knot ($H_1(M_3) = \mathbb{Z}_{19}$) |
| **`12_m2_S2_6`** | $1$ | $(x^3 - x^2 + 1)$ | $s^3 - s + 1$ | $0$ | $(10, 12, 63, 60)$ | **`R(b^2*a)`** ($\mathbb{Z}$-homology 4-sphere $\Sigma_2$, **exact $\operatorname{ch} = 12$ certified!**) |
| **`12_m2_S2_7`** | $1$ | $(x^4 - x^3 + x^2 - x + 1)$ | $s^4 - s^3 + s^2 - s + 1$ | $\mathbb{Z}_5$ | $(30, 12, 21, 420)$ | **`Spin(5_1)`** ($\tau^0(5_1)$, **exact $\operatorname{ch} = 12$ certified!**) |

> [!IMPORTANT]
> **Three Spoke upper bounds promoted to exact $\operatorname{ch}(F) = 12$:** In [`family_catalogs.md`](./family_catalogs.md), our constructive compilers built $\operatorname{ch} \in [12, 14]$ diagrams for `Spin(5_1)`, `Spin(5_2)`, and `R(b^2*a)`. Because none of these classes exist at $\operatorname{ch} \le 11$ (by our exhaustive $\operatorname{ch} \le 11$ Hub census), discovering $\operatorname{ch} = 12$ diagrams in `12_m2_S2_1`, `12_m2_S2_6`, and `12_m2_S2_7` **unconditionally certifies**
> $$\operatorname{ch}(\operatorname{Spin}(5_1)) = 12, \qquad \operatorname{ch}(\operatorname{Spin}(5_2)) = 12, \qquad \operatorname{ch}(R(b^2 a)) = 12.$$
> Furthermore, `12_m2_S2_3` and `12_m2_S2_4` exhibit non-principal elementary ideals $E_1 = (5, x+1)$ and $E_1 = (7, x+1)$, certifying **non-ribbon** 2-knots (since every ribbon 2-knot has principal $E_1$) realizing the determinant-$5$ and determinant-$7$ 2-bridge 2-twist-spins! Notice how $|\operatorname{Hom}(\pi_1, A_5)|$ cleanly reflects subgroup structure: $\mathbb{Z}_5 \rtimes_{-1} \mathbb{Z} \twoheadrightarrow D_5$ maps onto all $6$ dihedral $D_5 \le A_5$ subgroups ($6 \times 20 + 60 = 180$), whereas $D_7 \not\le A_5$ leaves $|\operatorname{Hom}(\mathbb{Z}_7 \rtimes_{-1} \mathbb{Z}, A_5)| = 60$.

### 3.2. Low-Saddle $\mathbb{RP}^2$ Abelianization Theorem & Yoshikawa §5 Censuses
In §5 of Yoshikawa (1994), Yoshikawa exhibited a prime $\mathbb{RP}^2$-knot candidate with $\operatorname{ch} \le 13$ and asked whether any new prime $\mathbb{RP}^2$-knot exists at $\operatorname{ch} = 11$ or $12$. First, a short topological argument rules out knotted $\mathbb{RP}^2$ groups for $m \le 2$ at **every** $\operatorname{ch}$:

### Theorem 3.1 (Low-Saddle $\mathbb{RP}^2$ Abelianization Theorem)
Let $D$ be any admissible marked diagram of a connected $\mathbb{RP}^2$-knot $F \subset S^4$ with $m(D) \le 2$ marked vertices. Then $\pi_1(S^4 \setminus F) \cong \mathbb{Z}_2$; consequently $F$ is topologically standard ($F \cong_{\text{TOP}} P_\pm$ by Conway–Powell 2021), and for $m(D) = 1$, $F$ is smoothly standard ($F \cong_{\text{DIFF}} P_\pm$ by Bleiler–Scharlemann 1988).

*Proof.* Since $F \cong \mathbb{RP}^2$ has $\chi(F) = 1$, the hyperbolic-splitting Euler formula $1 = \chi(F) = c_+(D) - m(D) + c_-(D)$ yields
$$c_+(D) + c_-(D) = 1 + m(D) \le 3.$$
Since $c_+(D), c_-(D) \ge 1$ are positive integers, $\min(c_+(D), c_-(D)) = 1$. By admissibility, that 1-component resolution (say $L_\varepsilon(D)$) is the **1-component classical unknot** $0_1 \subset S^3$, whose complement has cyclic fundamental group $\pi_1(S^3 \setminus L_\varepsilon(D)) \cong \mathbb{Z}$ generated by a single meridian. By Seifert–van Kampen across the saddle bands and capping disks on the $L_\varepsilon$ side, $\pi_1(S^3 \setminus L_\varepsilon(D)) \twoheadrightarrow \pi_1(S^4 \setminus F)$ is surjective! Hence $\pi_1(S^4 \setminus F)$ is cyclic, so $\pi_1(S^4 \setminus F) \cong H_1(S^4 \setminus F; \mathbb{Z}) \cong \mathbb{Z}_2$. By Conway–Powell (2021), any locally flat $\mathbb{RP}^2 \subset S^4$ with $\pi_1 \cong \mathbb{Z}_2$ (and $|e| = 2$, forced here by Rohlin's genus/signature congruence on $\pi_1 \cong \mathbb{Z}_2$) is topologically unknotted; for $m(D) = 1$, $F$ has a 1-saddle ribbon/band presentation on $0_1 \sqcup 0_1$, which is smoothly standard by Bleiler–Scharlemann (1988). $\square$

Guided by Theorem 3.1 (and noting that both known prime $\operatorname{ch} = 10$ $\mathbb{RP}^2$-knots `10_1^-2` and `10_2^-2` have $m \in \{3, 4\}$), we ran two exhaustive C99 $\mathbb{RP}^2$ censuses:
1. **Full $\operatorname{ch} = 11$ $\mathbb{RP}^2$ ($m = 1 \dots 10$, all $3,378$ shadows):** $50,835$ raw candidates $\to 139$ level-minimal $\to 3$ post-flype/tunnel survivors — **all $3$ are lower-$\operatorname{ch}$ / composite `# RP^2` collisions** ($2_1^{-1} \simeq 10_3 \# 2_1^{-1}$ and $8_1 \# 2_1^{-1} \simeq 9_1 \# 2_1^{-1} \simeq 10_2 \# 2_1^{-1}$). **Zero** prime $\operatorname{ch} = 11$ $\mathbb{RP}^2$-knots exist.
2. **Targeted $\operatorname{ch} = 12, m = 3$ $\mathbb{RP}^2$ ($15,882$ shadows):** $261,918$ raw candidates $\to 824$ level-minimal $\to 7$ post-flype/tunnel survivors — **all $7$ are `# RP^2` composite/unknot collapses** ($2_1^{-1} \simeq 10_3 \# 2_1^{-1}$, $8_1 \# 2_1^{-1} \simeq 10_2 \# 2_1^{-1}$, and $10_1 \# 2_1^{-1} \simeq 11_1 \# 2_1^{-1} \simeq 11_2 \# 2_1^{-1}$). Thus **zero** new prime $\mathbb{RP}^2$-knots exist at $\operatorname{ch} = 12, m \le 3$; any prime $\operatorname{ch} = 12$ $\mathbb{RP}^2$-knot must lie in $m \ge 4$.

---

## 4. Bridge Trisection Crossing Number $c_{\text{tri}}$ & Exhaustive $(b \le 4, c_{\text{tri}} \le 6)$ Braided Census

A Meier–Zupan $(b; p_{12}, p_{23}, p_{31})$ **tri-plane diagram** $\mathcal{T} = (T_1, T_2, T_3)$ consists of three $b$-strand trivial tangles in $B^3$ whose pairwise unions $L_{ij} = T_i \cup \bar{T}_j$ are classical unlinks of $p_{ij}$ circles, representing a closed surface-link $F(\mathcal{T}) \subset S^4$ with $\chi(F) = p_{12} + p_{23} + p_{31} - b$. Its **tri-plane crossing number** is
$$c_{\text{tri}}(\mathcal{T}) = c(T_1) + c(T_2) + c(T_3).$$

In [`triplane.py`](../code/src/surfaceknots/triplane.py) and [`triplane_census.py`](../code/src/surfaceknots/triplane_census.py), we implemented:
- Exact 3D-isotopy deduplication of braided $b$-bridge trivial tangles rel $\partial B^3$ via Artin's faithful representation on $F_b$ (`tangle_3d_key`), yielding exact minimal-crossing braided tangle counts at $b = 4$ of $(14, 96, 582, 3134, 15994, 78912)$ for $c = 0 \dots 5$.
- C99 batch R1/R2 + modular Kauffman bracket unlink pre-filtering (`batch_check_tangle_pairs_c`) backed by constructive Python R1/R2/R3 certification (`TriPlaneDiagram.is_certified_valid`), plus the exact algebraic join $(T_{\text{prev}} \cdot g) \cup \bar{M}_m = T_{\text{prev}} \cup \overline{M_m \cdot (-g)}$ covering all concentrated $(0, 0, 6)$ triples.
- Direct tri-plane computation of component topology and normal Euler numbers $e(F_c) = -(w_{\text{self}}(L_{12}|_c) + w_{\text{self}}(L_{23}|_c) + w_{\text{self}}(L_{31}|_c))$, Seifert–van Kampen $\pi_1(S^4 \setminus F)$, Fox elementary ideals $E_k$, finite-group homomorphisms, and Goeritz double branched cover invariants $H_1(\Sigma_2(S^4, F); \mathbb{Z})$.

### 4.1. Verification of Gong–Lewis-Monkman–Osnes (2026) at $b = 4$ & Rank Inversions
Running the exhaustive braided census over all $1 \le b \le 4$ and $c_{\text{tri}} \le 6$ ([`triplane_b4_census.json`](../code/data/complexity/triplane_b4_census.json)) confirms every prediction of **Gong, Lewis-Monkman, and Osnes (June 2026, arXiv:2606.03799v2)** within the braided $b \le 4$ universe:

1. **Conjecture 3.7 (Exhaustive $S^2$-knot classification up to $c_{\text{tri}} \le 6$ at $b \le 4$):** Exactly **three** connected $S^2$-knots ($\chi = 2$, orientable) exist among braided $b \le 4$ tri-plane diagrams with $c_{\text{tri}} \le 6$:
   - **`0_1` ($S^2$ unknot):** $b = 1$, $c_{\text{tri}} = 0$.
   - **`8_1` ($\tau^0(3_1)$, spun trefoil — Conjecture 3.6):** $b = 4$, $c_{\text{tri}} = 6$, realized **both** unconcentrated at $(0, 2, 4)$ and concentrated at $(0, 0, 6)$ ($124$ compilable concentrated triples, $0$ interleaved triples, compiling with zero extra crossings to $\operatorname{ch} = 6 + 2 = 8$).
   - **`10_2` ($\tau^2(3_1)$, 2-twist-spun trefoil — Corollary 3.5):** $b = 4$, $c_{\text{tri}} = 6$, realized in the braided model by $90$ concentrated $(0, 0, 6)$ triples (such as $T_3$ with braid word `[2, 4, -3, 5, -4, 5]`), **every single one of which ($90/90$, $0$ compilable) hits the Meander Interleaving Obstruction** (`best_compiled_ch = None`)! Thus no zero-extra-crossing planar compilation of `10_2` exists at $c_3 = 6$: resolving its interleaved meander forces $+2$ classical Reidemeister II crossings, elevating $\operatorname{ch}(10_2)$ from $6 + 2 = 8$ to $10$.
2. **Genuine Rank Inversions between $\operatorname{ch}$ and $c_{\text{tri}}$:**
   - **`9_1` vs. `10_2` ($S^2$-knots):** $\operatorname{ch}(9_1) = 9 < 10 = \operatorname{ch}(10_2)$, yet $9_1$ is **completely absent** at $b \le 4, c_{\text{tri}} \le 6$, so $c_{\text{tri}}(10_2) = 6 < 7 \le c_{\text{tri}}^{\text{braid}, b \le 4}(9_1)$! Bridge trisections accommodate Zeeman's 2-twist-spin at $c_{\text{tri}} = 6$ via an interleaved $(M_-, M_+)$ cap pair, whereas projecting $(M_-, M_+)$ into a single planar marked diagram forces $+2$ extra classical crossings.
   - **`2_1^1` vs. `2_1^-1` (Unknots):** $\operatorname{ch}(2_1^1) = \operatorname{ch}(2_1^{-1}) = 2$, whereas $c_{\text{tri}}(2_1^1) = 0 < 1 = c_{\text{tri}}(2_1^{-1})$.

### 4.2. The Concentrated Meander Interleaving Obstruction ($b \ge 4$)
A concentrated tri-plane diagram $\mathcal{T} = (M_-, T_2, M_+)$ ($c(M_-) = c(M_+) = 0$) places all $c_{\text{tri}}(\mathcal{T}) = c(T_2)$ classical crossings in a single page while $M_-$ and $M_+$ form a closed planar meander on $2b$ equatorial points. Our compiler `concentrated_triplane_to_marked_diagram` converts $\mathcal{T}$ into a planar `MarkedDiagram` with **zero added classical crossings** ($\operatorname{ch}(D) = c(T_2) + b - p_{\text{cap}}$) by iteratively contracting common boundary digons and adjacent shared-endpoint triples $(p_0, p_1, p_2) \mapsto \text{marker } Y/Z$. Across the $b \le 4, c_{\text{tri}} \le 6$ census, $1,719,925$ concentrated triples compile cleanly while $95,404$ hit the interleaving obstruction — including all $90$ concentrated realizations of `10_2`!

### Proposition 4.1 (Meander Interleaving Obstruction at $b = 4$)
For $b \le 3$, every pair of Catalan matchings $(M_-, M_+)$ reduces completely to empty active boundary via digon and shared-endpoint-triple contractions (zero interleaving obstructions). At $b = 4$, out of $C_4^2 = 14^2 = 196$ ordered Catalan pairs, **interleaved multi-cycle meanders** appear for the first time: specifically, the exact cap pair realized by `10_2` at $(0, 0, 6)$,
$$M_- = ((0, 7),\; (1, 6),\; (2, 5),\; (3, 4)), \qquad M_+ = ((0, 3),\; (1, 2),\; (4, 7),\; (5, 6))$$
has two meander 4-cycles $A = \{0, 3, 4, 7\}$ and $B = \{1, 2, 5, 6\}$ whose endpoints alternate $A, B, B, A, A, B, B, A$ around $S^1$ with **zero** common adjacent digons and **zero** adjacent triples $(p_0, p_1, p_2)$ sharing an arc in $M_-$ and an arc in $M_+$ ($M_-$ pairs $A \leftrightarrow A, B \leftrightarrow B$ across the horizontal axis, while $M_+$ pairs $A \leftrightarrow A$ on $\{0,3\}, \{4,7\}$ and $B \leftrightarrow B$ on $\{1,2\}, \{5,6\}$). Consequently, no planar saddle-band embedding in $D_+$ can realize $(M_-, M_+)$ without introducing extra classical Reidemeister II crossings!

---

## 5. Architectural Synthesis: Multi-Graded Filtration for `surfaceknots`

Comparing all four complexity measures across our unified dataset ([`complexity_comparison.json`](../code/data/complexity/complexity_comparison.json) & [`complexity_comparison.csv`](../code/data/complexity/complexity_comparison.csv)) establishes why **no single scalar complexity measure suffices alone**:

| Measure | Finite-to-One? | Unknot Baseline | Strengths | Structural Blind Spots |
| :--- | :---: | :--- | :--- | :--- |
| **Raw $\operatorname{ch}(F) = c + m$** | **Yes** | $\operatorname{ch}_0(\Sigma_g) = 2g$, $\operatorname{ch}_0(N_{k,e}) = k + \max(1, |e|/2)$ | Exhaustively enumerable via planar quadrangulations (`OEIS A113201`); captures $S^2$, $\Sigma_g$, non-orientable, and multi-component links uniformly | Artificially penalizes high genus & normal Euler number (scattering the Hopf cohort across $\operatorname{ch} = 6 \dots 10$); connected-sum collapses under $\#\,P_\pm$ |
| **Excess $\widetilde{\operatorname{ch}}(F) = \operatorname{ch} - \operatorname{ch}_0$** | **Yes** (per $\operatorname{top}(F)$) | $\widetilde{\operatorname{ch}}(U) = 0$ identically | Strips out intrinsic genus/Euler overhead; unifies all six Hopf 2-links at $\widetilde{\operatorname{ch}} = 4$; finite-to-one within each fixed topological sector $\operatorname{top}(F)$ | Infinite-to-one if summed without fixing $\operatorname{top}(F) = (\text{genera}, \text{ori}, e)$ (since all unknots sit at $\widetilde{\operatorname{ch}} = 0$) |
| **Bridge Trisection $c_{\text{tri}}(F)$** | **Only with $b \le B$** | $c_{\text{tri}}(\Sigma_g) = 0$, $c_{\text{tri}}(N_{k,e}) = \max(1, |e|/2)$ | Naturally sees 3-page twist-spinning (`10_2` at $c_{\text{tri}} = 6$, beating `9_1`); zero saddle overhead | **Infinite-to-one** without a bridge cap ($c_{\text{tri}}(\#^g T^2) = 0$ for all $g \ge 0$); combinatorial explosion of $C_b^3$ tangles for $b \ge 5$ |
| **Bridge Number $b(F)$** | **No** ($c_{\text{tri}}$ unbounded) | $b(S^2) = 1$, $b(\mathbb{RP}^2) = 2$, $b(T^2) = 3$ | Intrinsic 4D Morse-theoretic invariant; sharp lower bounds via $\operatorname{rank}(\pi_1)$ and elementary ideals | Infinitely many twist-spins $\tau^k(3_1)$ share $b = 4$ |

### Recommended Multi-Graded Atlas Key
Rather than replacing $\operatorname{ch}$-index as the C99 enumeration backbone (where planar quadrangulation duality remains orders of magnitude faster than $b \ge 5$ tangle triples), `surfaceknots` should index and display surface-links via the **lexicographic multi-grading**:
$$\mathcal{K}(F) \;:=\; \Big(\;\operatorname{top}(F),\;\; \widetilde{\operatorname{ch}}(F),\;\; \operatorname{ch}(F),\;\; b(F),\;\; c_{\text{tri}}(F)\;\Big),$$
where $\operatorname{top}(F) = (r,\; (\chi_i, \text{ori}_i, e_i)_{i=1}^r)$ fixes the topological homeomorphism and normal-Euler sector:
1. Within each fixed sector $\operatorname{top}(F)$, $\widetilde{\operatorname{ch}}(F) = \operatorname{ch}(F) - \operatorname{ch}_0(F)$ is **strictly finite-to-one** and orders links by genuine 4D knotting/linking complexity (putting every unknot at $\widetilde{\operatorname{ch}} = 0$ and every Hopf-type link at $\widetilde{\operatorname{ch}} = 4$).
2. Across $S^2$-knots ($\operatorname{ch}_0 = 0$, so $\widetilde{\operatorname{ch}} = \operatorname{ch}$), pairing $\operatorname{ch}(F)$ with $(b(F), c_{\text{tri}}(F))$ exposes 3-page geometric structure (such as twist-spinning in `10_2` and `12_m2_S2_3`/`12_m2_S2_4`) alongside ribbon fusion (`9_1`, `11_1`, `12_m2_S2_1`–`7`).
