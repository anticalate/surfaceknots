# Complexity Measures, Intrinsic Cutoffs, and Alternative Filtrations for Knotted Surfaces

This note analyzes the mathematical and computational behavior of Yoshikawa's hyperbolic-splitting complexity $\operatorname{ch}(F) = \min_D (c(D) + m(D))$ alongside intrinsic topological floors, weighted $(c, m)$ profiles, and Meier–Zupan's bridge trisection crossing number $c_{\text{tri}}(F)$. All theorems, censuses, and cross-measure comparisons below are verified by [`complexity.py`](../code/src/surfaceknots/complexity.py), [`triplane.py`](../code/src/surfaceknots/triplane.py), [`triplane_census.py`](../code/src/surfaceknots/triplane_census.py), and [`_fast_enum.c`](../code/src/surfaceknots/_fast_enum.c), with exported datasets under [`code/data/complexity/`](../code/data/complexity/) and [`expanded_ch12.json`](../code/data/expanded_ch12.json).

---

## 1. Intrinsic Topology vs. Knotting: Unknot $\operatorname{ch}$-Index Floors & Excess $\operatorname{ch}$-Index

A basic structural defect of raw $\operatorname{ch}$-index is that both classical crossings $c(D)$ and saddle markers $m(D)$ are forced upward merely by raising the genus or normal Euler number of an **unknotted** surface component. Separating this topological baseline from genuine 4-dimensional knotting or linking yields an exact additive floor $\operatorname{ch}_0(F)$.

### Theorem 1.1 (Zero-Crossing Orientability Theorem)
Let $D$ be an admissible marked diagram with $c(D) = 0$ classical crossings. Then the represented surface-link $F(D) \subset S^4$ is **orientable**, and every component $F_i \subset F(D)$ has normal Euler number $e(F_i) = 0$.

*Proof.* Since $c(D) = 0$, every vertex of $D$ is an $A$-/$B$-smoothed marker ($Y$ or $Z$). Both classical resolutions $L_-(D)$ and $L_+(D)$ are planar embedded collections of disjoint Jordan curves in $S^2$ with zero crossings. Two-color the complementary regions of $L_-(D) \subset S^2$ checkerboard-style (say white/black). Each marker vertex locally merges or splits two arcs separating opposite colors without introducing crossings; compatibility with the planar checkerboard coloring induces a coherent orientation on the resolution circles of $L_-(D)$ across every saddle band of $D$. By Kamada's orientability criterion on hyperbolic splittings, $F(D)$ is orientable. Moreover, since $L_-(D)$ and $L_+(D)$ each have zero classical crossings, Kamada's resolution-writhe formula under the codebase's time-first orientation $\mathbb{R}_t \times \mathbb{R}^3_{x,y,z}$ gives $e(F_i) = w(L_-(D)|_{F_i}) - w(L_+(D)|_{F_i}) = 0 - 0 = 0$. $\square$

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

For $U \cong N_{k, e}$ ($k \ge 1$), non-orientability implies $c(D) \ge 1$ by Theorem 1.1. Furthermore, under any orientation of the resolution circles of $L_-(D)$ and $L_+(D)$, each classical crossing $v$ contributes $\pm 1$ to $w(L_-)$ and $\pm 1$ to $w(L_+)$, hence at most $2$ to $|w(L_-) - w(L_+)| = |e(U)|$; thus $c(D) \ge |e(U)|/2$, proving $c(D) \ge \max(1, |e|/2)$ and $\operatorname{ch}(D) \ge k + \max(1, |e|/2)$.

For sharpness on $N_{k, e}$:
- When $e \ne 0$, set $r = |e|/2 \ge 1$ and $s = (k - r)/2 \in \mathbb{Z}_{\ge 0}$. Taking the planar connected sum of $r$ copies of $P_{\operatorname{sgn}(e)}$ (`2_1^-1`, $(c, m) = (1, 1)$) with $s$ copies of the unknotted torus `2_1^1` ($(c, m) = (0, 2)$) yields $N_{k, e}$ by Dyck's relation $P_\pm \# T^2 \cong P_\pm \# P_+ \# P_-$ with $(c, m) = (r, r + 2s) = (|e|/2, k)$.
- When $e = 0$ ($k = 2s \ge 2$), the admissible prime 3-vertex diagram (already noted by Yoshikawa 1994, Remark 2.8; stored as `MINIMAL_KLEIN_BOTTLE_E0_EPD` in [`complexity.py`](../code/src/surfaceknots/complexity.py))
  $$\texttt{X[1,2,3,4], Y[5,6,2,1], Y[5,4,3,6]}$$
  has $(c, m) = (1, 2)$, $c_+ = c_- = 1$ (each resolution is a 1-component classical unknot with one Reidemeister I kink), $\chi = 0$, non-orientable surface type $K^2$, $\pi_1 \cong \mathbb{Z}_2$, and $e = 0$, realizing $N_{2, 0} = P_+ \# P_-$ at $\operatorname{ch} = 3 = 2 + \max(1, 0)$. Taking its planar connected sum with $s - 1$ copies of `2_1^1` ($(c, m) = (0, 2)$) realizes every $N_{2s, 0}$ with $(c, m) = (1, 2s)$. $\square$

> [!NOTE]
> **Exact parallel with Allred et al. (2023):** Allred et al. (`arXiv:2302.08431`) proved that the bridge trisection crossing number of any unknotted non-orientable surface $N_{k, e} = \mathcal{P}_{p, q}$ equals $c_{\text{tri}}(N_{k, e}) = \max(1, |p - q|) = \max(1, |e|/2)$. Theorem 1.2 shows that this classical-crossing floor $c_0(N_{k, e}) = \max(1, |e|/2)$ is **identical** in marked ch-diagrams, while $\operatorname{ch}_0$ adds the exact Morse saddle cost $m_0(N_{k, e}) = k = 2 - \chi$.

### Corollary 1.3 (Strict Connected-Sum Subadditivity on Unknotted Klein Bottles; cf. Yoshikawa 1994, Remark 2.8)
Planar $\operatorname{ch}$-index is strictly subadditive even on standard unknots:
$$\operatorname{ch}(P_+ \# P_-) = 3 < 2 + 2 = \operatorname{ch}(P_+) + \operatorname{ch}(P_-), \qquad \text{whereas} \qquad \operatorname{ch}(P_+ \# P_+) = 4 = \operatorname{ch}(P_+) + \operatorname{ch}(P_+).$$
Geometrically, $P_+$ and $P_-$ carry opposite Reidemeister I kinks ($e = +2$ and $e = -2$); fusing them across a single mixed classical crossing cancels their normal Euler numbers while saving one classical crossing (certified constructively via an explicit 4-step Yoshikawa move path $\Omega_5 \to \Omega_7 \to \Omega_3 \to \Omega_1^-$ in [`certify_klein_bottle_subadditivity`](../code/src/surfaceknots/complexity.py)).

### Definition 1.4 (Excess $\operatorname{ch}$-Index $\widetilde{\operatorname{ch}}$)
For any surface-link $F = F_1 \sqcup \dots \sqcup F_r \subset S^4$, let $\operatorname{ch}_0(F) = \sum_{i=1}^r \operatorname{ch}_0(F_i^{\text{std}})$, where $F_i^{\text{std}}$ is the standard unknotted surface with the same orientability, genus, and normal Euler number $e(F_i)$. The **excess $\operatorname{ch}$-index** is
$$\widetilde{\operatorname{ch}}(F) := \operatorname{ch}(F) - \operatorname{ch}_0(F) \ge 0.$$
By Theorem 1.2, $\widetilde{\operatorname{ch}}(F) = 0$ if and only if $F$ admits a diagram achieving the unknot floor on every component (in particular, $\widetilde{\operatorname{ch}}(U) = 0$ for every standard unknot/unlink).

**Empirical validation on Yoshikawa's Hopf 2-link cohort:** In Yoshikawa's $\operatorname{ch} \le 10$ table, the six fundamental 2-component Hopf-linked surfaces (linking $S^2$, $T^2$, $\mathbb{RP}^2_\pm$, and $K^2$ via four inter-component classical crossings $c_{\text{cross}} = 4$) are artificially scattered across five raw $\operatorname{ch}$-levels $6, 7, 8, 9, 10$. Under $\widetilde{\operatorname{ch}}$, all six collapse to a **single constant complexity level** $\widetilde{\operatorname{ch}}(F) = 4$:

| Label | Topology $(F_1, F_2)$ | Normal Euler $e$ | Raw $\operatorname{ch}(F)$ | $(c_{\text{self}}, c_{\text{cross}}, m)$ | Intrinsic Floor $\operatorname{ch}_0(F)$ | Excess $\widetilde{\operatorname{ch}}(F)$ |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: |
| `6_1^{0,1}` | $S^2 \sqcup T^2$ | $(0, 0)$ | $6$ | $(0, 4, 2)$ | $0 + 2 = 2$ | **$4$** |
| `7_1^{0,-2}` | $S^2 \sqcup K^2$ | $(0, 0)$ | $7$ | $(1, 4, 2)$ | $0 + 3 = 3$ | **$4$** |
| `8_1^{1,1}` | $T^2 \sqcup T^2$ | $(0, 0)$ | $8$ | $(0, 4, 4)$ | $2 + 2 = 4$ | **$4$** |
| `8_1^{-1,-1}` | $\mathbb{RP}^2 \sqcup \mathbb{RP}^2$ | $(+2, -2)$ | $8$ | $(2, 4, 2)$ | $2 + 2 = 4$ | **$4$** |
| `9_1^{1,-2}` | $T^2 \sqcup K^2$ | $(0, 0)$ | $9$ | $(1, 4, 4)$ | $2 + 3 = 5$ | **$4$** |
| `10_1^{-2,-2}` | $K^2 \sqcup K^2$ | $(0, 0)$ | $10$ | $(2, 4, 4)$ | $3 + 3 = 6$ | **$4$** |

(Note that for $S^2$, any knotted or linked component must also pay at least one saddle pair $m \ge 2$ whenever $c_+ + c_- \ge 4$ — in `6_1^{0,1}` and `7_1^{0,-2}`, the $S^2$ component uses $0$ markers of its own because its inter-component crossings resolve into a single circle on both sides.)

**Weighted $(c, m)$ profiles $w_\alpha(D) = c(D) + \alpha m(D)$:** Varying the saddle weight $\alpha \in \{0.5, 1.0, 1.5, 2.0\}$ across our Hub + Spoke catalog reveals two distinct regimes:
- $\alpha < 1$ ($\alpha = 0.5$) favors high-genus ribbon tubes and orientable handles ($m$ cheap, $c$ expensive), pulling `2_1^1` ($w_{0.5} = 1.0$) below `2_1^-1` ($w_{0.5} = 1.5$) and pulling 2-fusion ribbon knots toward 2-twist-spins.
- $\alpha > 1$ ($\alpha = 2.0$) penalizes Morse saddles (approximating handle/ribbon-fusion complexity), separating 1-fusion ($m = 2$) 2-knots ($w_{2.0} = c + 4$) cleanly from $m \ge 3$ orientable/non-orientable surfaces.

---

## 2. Viro's $\#\,\mathbb{RP}^2$ Unknotting Collapse vs. Jabłonowski's Open Problem

Connected sum with a standard projective plane $P_\pm$ (`2_1^-1`, $e = \pm 2$, $\operatorname{ch} = 2$) causes striking collapses of 4D knotting invariants—some proved smoothly, others remaining open problems in 4-manifold topology:

### Theorem 2.1 (Viro's Unknotting Collapse, Viro 1973; cf. Yoshikawa 1994, Remark 2.8)
For the 3-twist-spun trefoil $10_3 = \tau^3(3_1)$ and standard projective plane $P_\pm = 2_1^{-1}$,
$$10_3 \# P_\pm \;=\; \tau^3(3_1) \# P_\pm \;\cong\; P_\pm \qquad \text{(smoothly in } S^4\text{)}.$$
Here $\operatorname{ch}$ collapses by **$10$ vertices**: $\operatorname{ch}(10_3) + \operatorname{ch}(P_\pm) = 10 + 2 = 12 \longrightarrow \operatorname{ch}(10_3 \# P_\pm) = 2$. Note that in [`complexity.py`](../code/src/surfaceknots/complexity.py) (`literature_true_ch`, exported in [`complexity_comparison.json`](../code/data/complexity/complexity_comparison.json)), this smooth equivalence is recorded as a literature theorem (Viro 1973) rather than a diagrammatic move certificate: bounded bidirectional Yoshikawa move search up to $\operatorname{ch} \le 15$ ([`move_search.py`](../code/src/surfaceknots/move_search.py)) does not yet cross the 3-twist-spin stabilization barrier to connect $10_3 \# P_+$ to $P_+$. Algebraically, connected sum with $P_\pm$ kills meridional squares in $\pi_1(S^4 \setminus (F \# P_\pm)) \cong \pi_1(S^4 \setminus F) / \langle\langle \mu^2 \rangle\rangle$: since $|\Delta_{10_3}(-1)| = 1$, $\pi_1(10_3)/\langle\langle \mu^2 \rangle\rangle \cong \mathbb{Z}_2$.

### Open Problem 2.2 (Jabłonowski 2015, Remark 3.15: $10_2 \# P_\pm \stackrel{?}{\cong} 8_1 \# P_\pm$)
For even twist-spins, Jabłonowski (`arXiv:1505.07995`, Remark 3.15; cf. Kamada) posed as an **open problem** whether the 2-twist-spun trefoil and Artin-spun trefoil become ambient isotopic after connected sum with $P_\pm$:
$$10_2 \# P_\pm \;=\; \tau^2(3_1) \# P_\pm \;\stackrel{?}{\cong}\; \tau^0(3_1) \# P_\pm \;=\; 8_1 \# P_\pm.$$
Algebraically, $8_1 \# P_\pm$, $9_1 \# P_\pm$, and $10_2 \# P_\pm$ all have isomorphic fundamental group $\pi_1 \cong S_3$, identical Alexander ideal $E_1 = (3x, x + 1) \subset \mathbb{Z}[x]/(x^2 - 1)$, and identical values across all 10 layers of our invariant filtration (and similarly at $\operatorname{ch} = 12, m = 3$, $10_1 \# P_\pm$, $11_1 \# P_\pm$, and $11_2 \# P_\pm$ share $\pi_1 \cong D_{10}$ and $E_1 = (5x, x + 1)$). Bidirectional Yoshikawa move search ([`move_search.py`](../code/src/surfaceknots/move_search.py), verified by the independent port-labeled checker `verify_move_path_independent`) explores the bounded stabilization frontier from $10_2 \# P_+$ and $8_1 \# P_+$ without finding a connecting move path, leaving Jabłonowski's question open.

---

## 3. C99 $\operatorname{ch} = 12$ Censuses: The `11` Prime $S^2$-Knots & Low-Saddle $\mathbb{RP}^2$-Knots

Using the dual plantri-style C99 quadrangulation generator (`generate_prime_shadows_c` in [`_fast_enum.c`](../code/src/surfaceknots/_fast_enum.c), enumerating all `15,882` prime 12-vertex spherical shadows in $0.86\,\text{s}$ and all **`77,185`** prime 13-vertex shadows in $4.76\,\text{s}$, matching `OEIS A113201`) combined with targeted marker-range scanning (`enumerate_targeted_census_level_c`) and mountain-pass $+1/+2$ tunneling, we executed exact strata censuses at $\operatorname{ch} = 11$ and $\operatorname{ch} = 12$ ([`targeted_ch12_13_census.json`](../code/data/complexity/targeted_ch12_13_census.json) and the full $m = 1 \dots 11$ census [`expanded_ch12.json`](../code/data/expanded_ch12.json)).

### 3.1. Complete $\operatorname{ch} = 12$ Connected $S^2$-Knots (`11` Classes, All Realized at $m = 2$)
Every prime connected $S^2$-knot ($\chi = 2$) has $m(D) \ge 2$, and across the entire $\operatorname{ch} = 12$ census ($m = 1 \dots 11$), **all** surviving prime $S^2$-knot diagrams occur in the minimal-saddle $m = 2$ stratum ($c_+ = c_- = 2$, $c = 10$):
- Raw admissible $m = 2$ candidates: **`55,038`** $\to$ **`1,158`** level-minimal $\Omega_{1..8}$ components $\to$ **`29`** post-tunneling/post-flype diagrams (`22` across new prime classes + `7` lower-$\operatorname{ch}$ `AMBIGUOUS_SHADOW` collisions).
- *(Audit update vs. pre-Gröbner `targeted_ch12_13_census.json`)*: Replacing $\pm 1$-evaluation `E1` hashes with strong integral Gröbner bases, adding Mochizuki's $R_5$ 3-cocycle state sums, and adding regular $F_{21}/F_{57}$ Frobenius-cover Betti profiles expands the initial 7-bucket pre-audit grouping into **`11` certified new prime $\operatorname{ch} = 12$ $S^2$-knot invariant classes** (`6` `EXACT_SINGLETON` + `5` `NEW_PACKET`; full proofs in [`ch12_expansion.md`](./ch12_expansion.md)):

| Hub ID | Comps (`Status`) | $E_1(\pi_1)$ | $H_1(\Sigma_2)$ | $|\operatorname{Hom}(\pi_1, G)|$: $(D_{10}, A_4, F_{21}, A_5)$ | Separating Witness & Identification |
| :--- | :---: | :--- | :---: | :---: | :--- |
| **`12_1`** | `3` (`PACKET`) | $(2x^2 - 3x + 2)$ | $\mathbb{Z}_7$ | $(10, 12, 21, 180)$ | **`Spin(5_2)`** ($\tau^0(5_2)$, **`MOVE_CERTIFIED`** exact $\operatorname{ch} = 12$) |
| **`12_2`** | `2` (`PACKET`) | $(3,\; x^2 + 1)$ | $0$ | $(10, 12, 21, 60)$ | **Non-ribbon, non-solvable $\mathbb{F}_9$ $2$-knot** (not $\tau^k(K)$ for $2 \le k \le 7$) |
| **`12_3`** | `1` (`SINGLE`) | $(3x^2 - 3x + 1)$ | $\mathbb{Z}_7$ | $(10, 36, 63, 180)$ | Asymmetric 1-fusion-type $2$-knot ($H_1(M_3) = \mathbb{Z}_2 \oplus \mathbb{Z}_{14}$) |
| **`12_4`** | `5` (`PACKET`) | $(5,\; x + 1)$ | $\mathbb{Z}_5$ | $(30, 12, 21, 180)$ | Non-ribbon $\tau^2(4_1)$ cohort ($\Phi_{R_5} = 5 + 10u^2 + 10u^3$) |
| **`12_5`** | `1` (`SINGLE`) | $(5,\; x + 1)$ | $\mathbb{Z}_5$ | $(30, 12, 21, 180)$ | Non-ribbon $\tau^2(5_1)$ ($\Phi_{R_5} = 5 + 10u + 10u^4$) |
| **`12_6`** | `2` (`PACKET`) | $(7,\; x + 1)$ | $\mathbb{Z}_7$ | $(10, 12, 21, 60)$ | Non-ribbon $\tau^2(5_2)$ cohort ($\pi_1 \cong \mathbb{Z}_7 \rtimes_{-1} \mathbb{Z}$) |
| **`12_7`** | `4` (`PACKET`) | $(x^3 - 2x^2 + 3x - 1)$ | $\mathbb{Z}_7$ | $(10, 12, 21, 60)$ | **`R(b^-1*a*b*a^-1*b*a)`** (`INVARIANT_MATCH`; $F_{57}$ $\operatorname{Tors} = \mathbb{Z}_{11}^2 \oplus \mathbb{Z}_{209}$) |
| **`12_8`** | `1` (`SINGLE`) | $(x^3 - 2x^2 + 3x - 1)$ | $\mathbb{Z}_7$ | $(10, 12, 21, 60)$ | **Non-ribbon Cappell–Shaneson `CS_A1`** ($F_{57}$ $\operatorname{Tors} = \mathbb{Z}_{19}$; top. up to Gluck twist) |
| **`12_9`** | `1` (`SINGLE`) | $(x^3 - x^2 + 1)$ | $0$ | $(10, 12, 63, 60)$ | **`R(b^2*a)`** (`INVARIANT_MATCH`, $\mathbb{Z}HS^4$ $\Sigma_2$, $\operatorname{ch} = 12$) |
| **`12_10`** | `1` (`SINGLE`) | $(x^3 - x^2 + 2x - 1)$ | $\mathbb{Z}_5$ | $(30, 12, 63, 180)$ | **Non-ribbon Cappell–Shaneson `CS_A0`** ($F_{21}$ $\operatorname{Tors} = \mathbb{Z}_7$ vs `11_2`; top. up to Gluck twist) |
| **`12_11`** | `1` (`SINGLE`) | $(x^4 - x^3 + x^2 - x + 1)$ | $\mathbb{Z}_5$ | $(30, 12, 21, 420)$ | **`Spin(5_1)`** ($\tau^0(5_1)$, **`MOVE_CERTIFIED`** exact $\operatorname{ch} = 12$) |

> [!IMPORTANT]
> **Degrees of certainty across the `11` $\operatorname{ch} = 12$ $S^2$-knots:** `12_1 = Spin(5_2)` and `12_11 = Spin(5_1)` are **`MOVE_CERTIFIED`** (provable ambient isotopy to compiled diagrams, with $\operatorname{ch} \ge 12$ conditional on $\operatorname{ch} \le 11$ completeness); `12_7 = R(b^-1*a*b*a^-1*b*a)` and `12_9 = R(b^2*a)` are **`INVARIANT_MATCH`** identifications across all 10 filtration layers; `12_2, 12_4, 12_5, 12_6` are **provably non-ribbon** via non-principal $E_1$; and `12_8` (`CS_A1`) and `12_10` (`CS_A0`) are **provably non-ribbon** via Ore-domain deficiency $\operatorname{def}(G_i) \le 0$, with exteriors homeomorphic to $CS(A_1)$ and $CS(A_0)$ (topological type determined up to the two non-isotopic Gluck fillings $K_{A_i}^0 \not\cong_{\mathrm{TOP}} K_{A_i}^1$).

### 3.2. Low-Saddle $\mathbb{RP}^2$ Abelianization Theorem & Exhaustive $\operatorname{ch} \le 12$ $\mathbb{RP}^2$ Censuses
Yoshikawa's $\operatorname{ch} \le 10$ table contains only a **single** connected $\mathbb{RP}^2$-knot — the standard projective plane `2_1^-1` at $(c, m) = (1, 1)$ — while reducible connected sums $K \# P_\pm$ begin at $\operatorname{ch}(8_1 \# 2_1^{-1}) = 10$ ($m = 3$), and smoothly irreducible $\mathbb{RP}^2 \subset S^4$ were constructed by Hughes–Kim–Miller–Nahm (2026, `arXiv:2605.12921`). First, a short topological argument rules out knotted $\mathbb{RP}^2$ groups for $m \le 2$ at **every** $\operatorname{ch}$:

### Theorem 3.1 (Low-Saddle $\mathbb{RP}^2$ Abelianization Theorem)
Let $D$ be any admissible marked diagram of a connected $\mathbb{RP}^2$-knot $F \subset S^4$ with $m(D) \le 2$ marked vertices. Then $\pi_1(S^4 \setminus F) \cong \mathbb{Z}_2$; consequently $F$ is topologically standard ($F \cong_{\mathrm{TOP}} P_\pm$ by Conway–Orson–Powell 2023, `arXiv:2306.12305`), and for $m(D) = 1$, $F$ is smoothly standard ($F \cong_{\mathrm{DIFF}} P_\pm$ by Bleiler–Scharlemann 1988).

*Proof.* Since $F \cong \mathbb{RP}^2$ has $\chi(F) = 1$, the hyperbolic-splitting Euler formula $1 = \chi(F) = c_+(D) - m(D) + c_-(D)$ yields
$$c_+(D) + c_-(D) = 1 + m(D) \le 3.$$
Since $c_+(D), c_-(D) \ge 1$ are positive integers, $\min(c_+(D), c_-(D)) = 1$. By admissibility, that 1-component resolution (say $L_\varepsilon(D)$) is the **1-component classical unknot** $0_1 \subset S^3$, whose complement has cyclic fundamental group $\pi_1(S^3 \setminus L_\varepsilon(D)) \cong \mathbb{Z}$ generated by a single meridian. By Seifert–van Kampen across the saddle bands and capping disks on the $L_\varepsilon$ side, $\pi_1(S^3 \setminus L_\varepsilon(D)) \twoheadrightarrow \pi_1(S^4 \setminus F)$ is surjective, so $\pi_1(S^4 \setminus F) \cong H_1(S^4 \setminus F; \mathbb{Z}) \cong \mathbb{Z}_2$. By Conway–Orson–Powell (2023, `arXiv:2306.12305`), any locally flat $\mathbb{RP}^2 \subset S^4$ with $\pi_1 \cong \mathbb{Z}_2$ (where $|e| = 2$ follows from $e \equiv 2\chi \equiv 2 \pmod 4$ and $|e| \le 2(2-\chi) = 2$) is topologically unknotted. When $m(D) = 1$, $c_+(D) + c_-(D) = 2$ forces $c_+(D) = c_-(D) = 1$: both resolutions $L_-(D)$ and $L_+(D)$ are single-component classical unknots $0_1$, giving $F$ a Morse presentation with $c_- + m + c_+ = 1 + 1 + 1 = 3$ critical points, which is smoothly standard ($F \cong_{\mathrm{DIFF}} P_\pm$) by Bleiler–Scharlemann (1988). $\square$

Guided by Theorem 3.1, our C99 $\mathbb{RP}^2$ censuses establish:
1. **Full $\operatorname{ch} = 11$ $\mathbb{RP}^2$ ($m = 1 \dots 10$, all `3,378` shadows; [`targeted_ch12_13_census.json`](../code/data/complexity/targeted_ch12_13_census.json)):** `97,824` raw candidates $\to$ `496` level-minimal $\to$ `4` C99-tunneled $\to$ `3` post-flype survivors — **all `3` collide with lower-$\operatorname{ch}$ / composite `# RP^2` classes** ($2_1^{-1} \cong 10_3 \# 2_1^{-1}$ and $8_1 \# 2_1^{-1} / 9_1 \# 2_1^{-1} / 10_2 \# 2_1^{-1}$). **Zero** new prime $\mathbb{RP}^2$-knots exist at $\operatorname{ch} = 11$.
2. **Targeted $m = 3$ & Full $m = 1 \dots 11$ $\operatorname{ch} = 12$ $\mathbb{RP}^2$ (all `15,882` shadows; [`targeted_ch12_13_census.json`](../code/data/complexity/targeted_ch12_13_census.json) & [`expanded_ch12.json`](../code/data/expanded_ch12.json)):** At $m = 3$, `261,918` raw candidates $\to$ `824` level-minimal $\to$ `8` C99-tunneled $\to$ `7` post-flype survivors, and across **all** $m = 1 \dots 11$ in `expanded_ch12.json`, every surviving $\mathbb{RP}^2$ diagram collides with $2_1^{-1}$ or reducible `# RP^2` connected sums ($\pi_1 \cong \mathbb{Z}_2, S_3, D_{10}$). Thus **zero** new prime $\mathbb{RP}^2$-knots exist anywhere at $\operatorname{ch} \le 12$: any prime (or irreducible) $\mathbb{RP}^2$-knot must satisfy $\operatorname{ch}(F) \ge 13$.

---

## 4. Bridge Trisection Crossing Number $c_{\text{tri}}$ & $(b \le 4, c_{\text{tri}} \le 6)$ Braided Census

A Meier–Zupan $(b; p_{12}, p_{23}, p_{31})$ **tri-plane diagram** $\mathcal{T} = (T_1, T_2, T_3)$ consists of three $b$-strand trivial tangles in $B^3$ whose pairwise unions $L_{ij} = T_i \cup \bar{T}_j$ are classical unlinks of $p_{ij}$ circles, representing a closed surface-link $F(\mathcal{T}) \subset S^4$ with $\chi(F) = p_{12} + p_{23} + p_{31} - b$. Its **tri-plane crossing number** is
$$c_{\text{tri}}(\mathcal{T}) = c(T_1) + c(T_2) + c(T_3).$$

In [`triplane.py`](../code/src/surfaceknots/triplane.py) and [`triplane_census.py`](../code/src/surfaceknots/triplane_census.py), we implemented:
- Invariant-key deduplication of braided $b$-bridge trivial tangles via strand-endpoint permutations and finite-group arc-coloring multisets (`D_6, D_8, D_{10}, A_4` in `tangle_3d_key` — a sound invariant quotient finer than 3D tangle isotopy rel $\partial B^3$), yielding `(14, 96, 582, 3134, 15994, 78912)` braided `tangle_3d_key` buckets at $b = 4$ for $c = 0 \dots 5$.
- C99 batch R1/R2 + modular Kauffman bracket unlink pre-filtering (`batch_check_tangle_pairs_c`) backed by constructive Python R1/R2/R3 certification (`TriPlaneDiagram.is_certified_valid`), plus the exact algebraic join $(T_{\text{prev}} \cdot g) \cup \bar{M}_m = T_{\text{prev}} \cup \overline{M_m \cdot (-g)}$ covering concentrated $(0, 0, 6)$ triples. *(Dataset note: the stored [`triplane_b4_census.json`](../code/data/complexity/triplane_b4_census.json) reflects a pre-audit run with `36,166,398` ingested triples and `26,808,485` triples skipped by `pruned_shared_untouched_arch_b4` — which artificially raised split-unlink minima such as $0_1^{\sqcup 4}$ before that filter was removed in `triplane_census.py` to scan all `62,974,883` valid braided $b \le 4, c_{\text{tri}} \le 6$ triples.)*
- Direct tri-plane computation of component topology and normal Euler numbers, Seifert–van Kampen $\pi_1(S^4 \setminus F)$, Fox elementary ideals $E_k$, finite-group homomorphisms, and Goeritz double branched cover invariants $H_1(\Sigma_2(S^4, F); \mathbb{Z})$.

### 4.1. Consistency with Gong–Lewis-Monkman–Osnes (2026) on the Braided $b \le 4$ Locus & Rank Inversions
Within the **braided $b \le 4$ model** (noting that Gong, Lewis-Monkman, and Osnes, June 2026, `arXiv:2606.03799v2`, formulate Conjectures 3.6 & 3.7 over all bridge numbers $b$ and general unbraided trivial tangles), the census up to $c_{\text{tri}} \le 6$ finds:

1. **Connected $S^2$-knots up to $c_{\text{tri}} \le 6$ among braided $b \le 4$ tri-planes (cf. GLMO Conjecture 3.7):** Exactly **three** connected $S^2$-knots ($\chi = 2$, orientable) appear in the braided $b \le 4, c_{\text{tri}} \le 6$ universe:
   - **`0_1` ($S^2$ unknot):** $b = 1$, $c_{\text{tri}} = 0$.
   - **`8_1` ($\tau^0(3_1)$, spun trefoil — cf. Conjecture 3.6):** $b = 4$, $c_{\text{tri}} \le 6$ (realized at $c_{\text{tri}} = 6$ both unconcentrated at $(0, 2, 4)$ and concentrated at $(0, 0, 6)$, where all `124` concentrated triples are non-interleaved and compile with zero extra crossings to $\operatorname{ch} = 6 + 2 = 8$).
   - **`10_2` ($\tau^2(3_1)$, 2-twist-spun trefoil — Corollary 3.5):** $b = 4$, $c_{\text{tri}} = 6$, realized in the braided model by `90` concentrated $(0, 0, 6)$ triples (such as $T_3$ with braid word `[2, 4, -3, 5, -4, 5]`), all `90` of which hit the concentrated Meander Interleaving Obstruction (`best_compiled_ch = None`) under direct zero-crossing cap contraction.
2. **Genuine Rank Inversions between $\operatorname{ch}$ and $c_{\text{tri}}$:**
   - **`9_1` vs. `10_2` ($S^2$-knots):** $\operatorname{ch}(9_1) = 9 < 10 = \operatorname{ch}(10_2)$, yet $9_1$ is **completely absent** from braided $b \le 4, c_{\text{tri}} \le 6$ tri-planes ($c_{\text{tri}}^{\text{braid}, b \le 4}(9_1) \ge 7 > 6 = c_{\text{tri}}(10_2)$).
   - **`2_1^1` vs. `2_1^-1` (Unknots):** $\operatorname{ch}(2_1^1) = \operatorname{ch}(2_1^{-1}) = 2$, whereas $c_{\text{tri}}(2_1^1) = 0 < 1 = c_{\text{tri}}(2_1^{-1})$.

### 4.2. The Concentrated Meander Interleaving Obstruction ($b \ge 4$)
A concentrated tri-plane diagram $\mathcal{T} = (M_-, T_2, M_+)$ ($c(M_-) = c(M_+) = 0$) places all $c_{\text{tri}}(\mathcal{T}) = c(T_2)$ classical crossings in a single page while $M_-$ and $M_+$ form a closed planar meander on $2b$ equatorial points. Our compiler `concentrated_triplane_to_marked_diagram` converts $\mathcal{T}$ into a planar `MarkedDiagram` with **zero added classical crossings** ($\operatorname{ch}(D) = c(T_2) + b - p_{\text{cap}}$) whenever $(M_-, M_+)$ reduces completely via common boundary digons and adjacent shared-endpoint triples $(p_0, p_1, p_2) \mapsto \text{marker } Y/Z$.

### Proposition 4.1 (Meander Interleaving Obstruction to Direct Zero-Crossing Compilation at $b = 4$)
For $b \le 3$, every pair of Catalan matchings $(M_-, M_+)$ reduces completely to empty active boundary via digon and shared-endpoint-triple contractions (zero interleaving obstructions). At $b = 4$, out of $C_4^2 = 14^2 = 196$ ordered Catalan pairs, **interleaved multi-cycle meanders** appear for the first time: specifically, the cap pair realized on `10_2`'s concentrated $(0, 0, 6)$ diagrams,
$$M_- = ((0, 7),\; (1, 6),\; (2, 5),\; (3, 4)), \qquad M_+ = ((0, 3),\; (1, 2),\; (4, 7),\; (5, 6))$$
has two meander 4-cycles $A = \{0, 3, 4, 7\}$ and $B = \{1, 2, 5, 6\}$ whose endpoints alternate $A, B, B, A, A, B, B, A$ around $S^1$ with **zero** common adjacent digons and **zero** adjacent triples $(p_0, p_1, p_2)$ sharing an arc in $M_-$ and an arc in $M_+$. Consequently, *that fixed concentrated cap pair* cannot be contracted into planar saddle bands in $D_+$ without either introducing extra classical Reidemeister II crossings or changing the tri-plane presentation. *(Caution: meander interleaving on a particular concentrated diagram obstructs only direct compilation of that diagram, not $\operatorname{ch}(F)$ globally — e.g., $0_1$ and $8_1^{-1,-1}$ also admit interleaved concentrated $(0, 0, 6)$ diagrams while achieving $\operatorname{ch} = 0$ and $\operatorname{ch} = 8$ via non-interleaved diagrams; the global lower bound $\operatorname{ch}(10_2) = 10$ is certified independently by the exhaustive $\operatorname{ch} \le 9$ census.)*

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
Rather than replacing $\operatorname{ch}$-index as the C99 enumeration backbone (where planar quadrangulation duality remains orders of magnitude faster than $b \ge 5$ tangle triples), `surfaceknots` indexes surface-links via the **lexicographic multi-grading**:
$$\mathcal{K}(F) \;:=\; \Big(\;\operatorname{top}(F),\;\; \widetilde{\operatorname{ch}}(F),\;\; \operatorname{ch}(F),\;\; b(F),\;\; c_{\text{tri}}(F)\;\Big),$$
where $\operatorname{top}(F) = (r,\; (\chi_i, \text{ori}_i, e_i)_{i=1}^r)$ fixes the topological homeomorphism and normal-Euler sector:
1. Within each fixed sector $\operatorname{top}(F)$, $\widetilde{\operatorname{ch}}(F) = \operatorname{ch}(F) - \operatorname{ch}_0(F)$ is **strictly finite-to-one** and orders links by genuine 4D knotting/linking complexity (putting every unknot at $\widetilde{\operatorname{ch}} = 0$ and every Hopf-type link at $\widetilde{\operatorname{ch}} = 4$).
2. Across $S^2$-knots ($\operatorname{ch}_0 = 0$, so $\widetilde{\operatorname{ch}} = \operatorname{ch}$), pairing $\operatorname{ch}(F)$ with $(b(F), c_{\text{tri}}(F))$ exposes 3-page geometric structure (such as twist-spinning in `10_2` and `12_4`–`12_6`) alongside ribbon fusion, non-solvable non-ribbon modules, and Cappell–Shaneson fiberings (`9_1`, `11_1`, `11_2`, `12_1`–`12_11`).
