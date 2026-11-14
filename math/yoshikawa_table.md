# Machine-Verified & Annotated Yoshikawa Table ($\operatorname{ch}(F) \le 10$)

Yoshikawa (1994) enumerated all **23 weakly prime, non-splittable surface-links** in $\mathbb{R}^4$ (up to ambient isotopy, surface orientation, and $\mathrm{O}(4)$ mirror image) with $\operatorname{ch}(F) \le 10$. Every entry below has been independently regenerated from 4-regular planar quadrangulation shadows (`A113201`), certified admissible via explicit Reidemeister unlink movies $(M_-, M_+)$, and annotated by the `surfaceknots` computational suite (exported in [`yoshikawa_ch10.json`](../code/data/yoshikawa_ch10.json) and [`yoshikawa_ch10.csv`](../code/data/yoshikawa_ch10.csv)).

Notation $I_k^{g_1, \dots, g_c}$:
- $I = \operatorname{ch}(F)$ is the minimal `ch`-index $\min_D (c(D) + h(D))$.
- $k$ is the index within fixed $(I; g_1, \dots, g_c)$.
- $g_i \ge 0$ denotes an orientable component $\Sigma_{g_i}$ of genus $g_i$ ($0 \mapsto S^2$, $1 \mapsto T^2$; omitted for connected $S^2$-knots).
- $g_i < 0$ denotes a non-orientable component $\#_{|g_i|} \mathbb{R}P^2$ of non-orientable genus $|g_i|$ ($-1 \mapsto \mathbb{R}P^2$, $-2 \mapsto K^2 = \mathbb{R}P^2 \# \mathbb{R}P^2$ Klein bottle).

---

## 1. Complete Machine-Verified Table ($\operatorname{ch}(F) \le 10$)

Sources for geometric/trisection literature columns:
- **Topology, Constructions**: Yoshikawa (1994), *Osaka J. Math.* 31(3), 497–522; **Table II Alexander/group errata independently verified and corrected below** (matching Ishii–Nikkuni–Oshiro 2018, `arXiv:1503.07969`, §5).
- **Ribbon / Pseudo-Ribbon & Triple Point Number $t(F)$**: Cazet (2022), `arXiv:2205.11120` (compiling Satoh 2001, Satoh–Shima 2004/2005, Kamada–Oshiro 2009).
- **Minimal Bridge Number $b(F)$, Tri-Plane Crossing Number $c_{\text{tri}}(F)$, Normal Euler Number $e(F)$**: Allred et al. & Zupan (2023), `arXiv:2302.08431`, and Gong–Lewis-Monkman–Osnes (2026), `arXiv:2606.03799` (proving $c_{\text{tri}}(10_2) = 6$).
- **Kirby–Thompson Invariants $L(F), L^\star(F)$**: Blair et al. (2022), Aranda et al. (2022, 2024), and Taniguchi (2025), `arXiv:2511.14695`.

### Part A: Connected Surface-Knots ($c = 1$; 9 entries: six $S^2$, two $T^2$, one $\mathbb{R}P^2$)

| Label | Type | Ori. | $\chi$ | $e(F)$ | Construction | Ribbon? | $t(F)$ | $b(F)$ | $c_{\text{tri}}(F)$ | $L(F)$ | $L^\star(F)$ | $\pi_1(\mathbb{R}^4 \setminus F)$ | Verified $E_1(F) \subset \mathbb{Z}[H_1]$ | Machine Distinguishing Invariants |
| :--- | :--- | :---: | :---: | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- | :--- | :--- |
| $0_1$ | $S^2$ | Yes | $2$ | $0$ | Unknotted 2-sphere | Yes | $0$ | $1$ | $0^c$ | $0$ | $0$ | $\mathbb{Z}$ | $(1)$ | `S3=6, Q8=8`, $\Phi_3 = (3, 3)$ |
| $8_1$ | $S^2$ | Yes | $2$ | $0$ | Spun trefoil $\tau^0(3_1)$ | Yes | $0$ | $4$ | $\le 6^c$ | $15$ | $12$ | $\pi_1(S^3 \setminus 3_1)$ | $(x^2 - x + 1)$ | `S3=12, A5=360`, $\Phi_3 = (9, 9)$, `pSL2=0` |
| $9_1$ | $S^2$ | Yes | $2$ | $0$ | 1-fusion ribbon 2-knot (slice of $6_1$) | Yes | $0$ | $4$ | $\le 7^c$ ($\ge 6$) | $15\text{–}16$ | $12\text{–}13$ | $\operatorname{BS}(1,2) = \langle a, x \mid x a x^{-1} = a^2 \rangle$ | $(x - 2)$ | `S3=12, F21=63`, $s$-biquandle $(s - 2)$ |
| $10_1$ | $S^2$ | Yes | $2$ | $0$ | Spun figure-eight $\tau^0(4_1)$ | Yes | $0$ | $4$ | $\le 8^c$ ($\ge 6$) | det. | det. | $\pi_1(S^3 \setminus 4_1)$ | $(x^2 - 3x + 1)$ | `D10=30, A4=36, A5=300`, $s$-biquandle $(s^2 - 3s + 1)$ |
| $10_2$ | $S^2$ | Yes | $2$ | $0$ | 2-twist-spun trefoil $\tau^2(3_1)$ | **No** | $4$ | $4$ | **$6$** | det. | $12$ | $\mathbb{Z} \ltimes \mathbb{Z}_3$ | $(x + 1, 3)$ | **Chiral / Non-invertible**: $\Phi_3 = (3 + 6u, 3 + 6u^2)$ |
| $10_3$ | $S^2$ | Yes | $2$ | $0$ | 3-twist-spun trefoil $\tau^3(3_1)$ | **No** | $6$ | $4$ | $6\text{–}11$ | $15\text{–}34$ | $12\text{–}34$ | $[\pi_1, \pi_1] \cong Q_8$ | $(x^2 + x + 1, 2)$ | `S3=6, A4=36, A5=180` |
| $2_1^1$ | $T^2$ | Yes | $0$ | $0$ | Unknotted torus | Yes | $0$ | $3^\dagger$ | $0^c$ | $0$ | $0$ | $\mathbb{Z}$ | $(1)$ | `S3=6`, `pSL2=0` |
| $10_1^1$ | $T^2$ | Yes | $0$ | $0$ | Tubed trefoil $\operatorname{Tube}(3_1)$ | Yes | $0$ | $6$ | $\le 11$ | $\le 24$ | $\le 18$ | $\pi_1(S^3 \setminus 3_1)$ | $(x^2 - x + 1)$ | **`pSL2=24`** (`ord(P_0)=6 > 3=ord(mu_0)` vs `0` on $8_1 \# 2_1^1$!) |
| $2_1^{-1}$ | $\mathbb{R}P^2$ | No | $1$ | $\pm 2$ | Standard projective plane $P_\pm$ | **No** (pseudo) | $0$ | $2$ | **$1^c$** | $0$ | $0$ | $\mathbb{Z}_2$ | $(1)$ | `H_1 = Z_2`, `S3=4, Q8=2` |

> $^\dagger$ **Literature Erratum Caught (`Allred et al. 2023`, Table 1)**: Allred et al.'s Table 1 prints `Bridge # = 1` for $2_1^1$ ($T^2$), a copy-paste typo from $0_1$: as proved in their own Theorem 1.2 (p. 10), $\chi(T^2) = c_1 + c_2 + c_3 - b = 0$ with $c_i \ge 1$ forces $b(2_1^1) = 3$.

---

### Part B: Multi-Component Surface-Links ($c \ge 2$; 14 entries — all individual components certified unknotted!)

| Label | Components | Ori. | $\chi$ | $e(F_i)$ | Geometric / Sublink Notes | Ribbon? | $t(F)$ | $b(F)$ | $c_{\text{tri}}(F)$ | Verified $\pi_1(\mathbb{R}^4 \setminus F)$ | Verified $E_1(F) \subset \mathbb{Z}[H_1]$ | Machine Distinguishing Invariants |
| :--- | :--- | :---: | :---: | :---: | :--- | :---: | :---: | :---: | :---: | :--- | :--- | :--- |
| $6_1^{0,1}$ | $S^2 \sqcup T^2$ | Yes | $2$ | $(0, 0)$ | Spun Hopf link | Yes | $0$ | $4$ | $\le 4^c$ | $\mathbb{Z} \oplus \mathbb{Z}$ | $(x - 1, y - 1)$ | `S3=18, Q8=40`, `pSL2=69` |
| $7_1^{0,-2}$ | $S^2 \sqcup K^2$ | Mixed | $2$ | $(0, 0)$ | Simplest mixed ribbon link | Yes | $0$ | $4$ | $\le 5^c$ | $\langle x, y \mid yxyx^{-1} = 1 \rangle \cong \mathbb{Z} \rtimes \mathbb{Z}$ | $(x + 1, y - 1)$ | `H_1 = Z + Z_2`, `S3=18, Q8=40, Dic16=112` |
| $8_1^{1,1}$ | $T^2 \sqcup T^2$ | Yes | $0$ | $(0, 0)$ | 0-twist torus-spun Hopf link | Yes | $0$ | $6$ | $\le 6$ | $\mathbb{Z} \oplus \mathbb{Z}$ | $(x - 1, y - 1)$ | **Sanderson $\operatorname{Dlk} = (0, 0)$** |
| $8_1^{-1,-1}$ | $\mathbb{R}P^2 \sqcup \mathbb{R}P^2$ | No | $2$ | $(+2, -2)$ | Linked $P_+ \sqcup P_-$, meridian order $4$ | **No** | **$2$** | $4$ | $\le 8$ | $Q_8$ (order 8) | $(x + 1, y + 1, 2)$ | `H_1 = Z_2^2`, `S3=10, Q8=28, Dic16=52` |
| $9_1^{0,1}$ | $S^2 \sqcup T^2$ | Yes | $2$ | $(0, 0)$ | Ribbon $S^2 \sqcup T^2$ link | Yes | $0$ | $4$ | $\le 8^c$ | 2-gen, 1-relator ribbon group | $(x - 1, y - 1)(y - 1)$ | `S3=24, Q8=64, F21=189` |
| $9_1^{1,-2}$ | $T^2 \sqcup K^2$ | Mixed | $0$ | $(0, 0)$ | Ribbon $T^2 \sqcup K^2$ link | Yes | $0$ | $6$ | $\le 7$ | **$\mathbb{Z} \oplus \mathbb{Z}_2$** (**Erratum #1 fixed!**) | **$(x - 1, y - 1, 2)$** (**Erratum #1 fixed!**) | `S3=12, Q8=16, Dic16=32` (abelian $\pi_1$!) |
| $10_1^{0,1}$ | $S^2 \sqcup T^2$ | Yes | $2$ | $(0, 0)$ | Ribbon $S^2 \sqcup T^2$ link | Yes | $0$ | $4$ | $\le 8$ | $\langle x, y \mid y^{-1} x^{-1} y x^{-1} y x y^{-1} x = 1 \rangle$ | **$(x - 1, y - 1)(x + y)$** (**Erratum #2 fixed!**) | `S3=30`, $\Phi_3=(3,3,3,3)$, $s$-biq $(2s-2, s^3-s^2+s-1)$ |
| $10_2^{0,1}$ | $S^2 \sqcup T^2$ | Yes | $2$ | $(0, 0)$ | Ribbon $S^2 \sqcup T^2$ link | Yes | $0$ | $4$ | $\le 8^c$ | $\langle x, y \mid x^2 y x^{-2} y^{-1} = 1 \rangle$ | $(x - 1, y - 1)(x + 1)$ | `S3=30`, $\Phi_3=(9,9,9,9)$, $s$-biq $(s^2 - 1)$ |
| $10_1^{1,1}$ | $T^2 \sqcup T^2$ | Yes | $0$ | $(0, 0)$ | 1-twist-spun Hopf link | Pseudo | $0$ | $6$ | $\le 6$ | $\mathbb{Z} \oplus \mathbb{Z}$ | $(x - 1, y - 1)$ | **Sanderson $\operatorname{Dlk} = (1, 1)$** (separates from $8_1^{1,1}$!) |
| $10_1^{0,0,1}$ | $S^2 \sqcup S^2 \sqcup T^2$ | Yes | $4$ | $(0, 0, 0)$ | 3-component ribbon link: $(S^2_1, S^2_2)$ split; both $(S^2_i, T^2) \cong 6_1^{0,1}$ | Yes | $0$ | $5$ | $\le 8$ | $\langle x, y, z \mid y^{-1} x^{-1} z x y z^{-1} = 1 \rangle$ | $(0)$ | `S3=108, Q8=320, A5=18000` |
| $10_1^{0,-2}$ | $S^2 \sqcup K^2$ | Mixed | $2$ | $(0, 0)$ | Ribbon $S^2 \sqcup K^2$ link | Yes | $0$ | $4$ | $\le 8^c$ | $\langle x, y \mid x^{-1} y^{-1} x y x^{-1} y x y = 1 \rangle$ | $(2x + y - 1, 4)$ | `S3=12, D8=48, Q8=16, Dic16=96` |
| $10_2^{0,-2}$ | $S^2 \sqcup K^2$ | Mixed | $2$ | $(0, 0)$ | Non-ribbon candidate ($e(K^2)=0$) | ? | $0\text{–}10$ | $4$ | $\le 8^{*,c}$ | $\langle x, y_1, y_2 \mid x y_1 x^{-1} = y_2,\, y_1^2 = y_2^2 = (y_1 y_2)^2 \rangle$ | **$(2x + y - 1, (x + 1)(y + 1), 4)$** (**Erratum #3 fixed!**) | Separated from $10_1^{0,-2}$ by $E_1$ (`(x+1)(y+1)` generator)! |
| $10_1^{-1,-1}$ | $\mathbb{R}P^2 \sqcup \mathbb{R}P^2$ | No | $2$ | $(+2, -2)$ | Non-orientable link, meridian order $4$ | **No** | $2\text{–}12$ | $4$ | $\le 8^c$ | $\operatorname{Dic}_{16} = \langle 2, 2, 4 \rangle$ (order 16) | $(x + 1, y + 1, 4)$ | `S3=10, D8=36, Q8=4, Dic16=36, S4=76` |
| $10_1^{-2,-2}$ | $K^2 \sqcup K^2$ | No | $0$ | $(0, 0)$ | Ribbon link of two Klein bottles | Yes | $0$ | $6$ | $\le 8$ | $Q_8$ (order 8) | $(x + 1, y + 1, 2)$ | `S3=10, Q8=28`, $e=(0,0)$ vs $(+2,-2)$ on $8_1^{-1,-1}$ |

---

## 2. Independent Verification of All Three Yoshikawa Table II Errata

Running `check_e1_ideal_against_reference` (exact integral Hermite Normal Form over $\mathbb{Z}[t_1^{\pm 1}, \dots, t_r^{\pm 1}, u_1, \dots, u_s]/(u_j^2 - 1)$ modulo finite group-ring quotients + evaluation fingerprints, independently cross-checked on both $L_-(D)$ and $L_+(D)$) across all 23 canonical diagrams confirms exact agreement on **20 of 23** entries and isolates **all three errata in Yoshikawa's (1994) Table II** noted by Ishii–Nikkuni–Oshiro (2018, `arXiv:1503.07969`, §5):

1. **Erratum #1 — $9_1^{1,-2}$ ($T^2 \sqcup K^2$)**:
   - *Yoshikawa Table II printed*: $\pi_1 = \langle x, y \mid yxyx^{-1} = 1 \rangle$ and $E_1 = (x + 1, y - 1)$ (identical to $7_1^{0,-2}$).
   - *Machine-verified exact value*: $\pi_1(\mathbb{R}^4 \setminus 9_1^{1,-2}) = \langle x, y \mid xy = yx,\; y^2 = 1 \rangle \cong \mathbb{Z} \oplus \mathbb{Z}_2$ (**abelian**; $|\operatorname{Hom}(\pi_1, S_3)| = 12$ vs $18$ for $7_1^{0,-2}$), and
     $$E_1(9_1^{1,-2}) = (x - 1,\; y - 1,\; 2) = (x - 1,\; y + 1,\; 2) \subset \mathbb{Z}[x^{\pm 1}, y]/(y^2 - 1).$$
2. **Erratum #2 — $10_1^{0,1}$ ($S^2 \sqcup T^2$)**:
   - *Yoshikawa Table II printed*: $E_1 = (x - 1, y - 1)(x + 1)$ (making $E_1(10_1^{0,1})$ appear identical to $E_1(10_2^{0,1})$).
   - *Machine-verified exact value*: From $\pi_1 = \langle x, y \mid y^{-1} x^{-1} y x^{-1} y x y^{-1} x = 1 \rangle$, Fox calculus gives
     $$E_1(10_1^{0,1}) = (x^2 - y^2,\; xy - x + y^2 - y) = (x - 1,\; y - 1)(x + y)$$
     (up to $y \leftrightarrow y^{-1}$). In particular, $E_1$ alone already separates $10_1^{0,1}$ from $10_2^{0,1}$!
3. **Erratum #3 — $10_2^{0,-2}$ ($S^2 \sqcup K^2$)**:
   - *Yoshikawa Table II printed*: $E_1 = (2x + y - 1, 4)$ (making $E_1(10_2^{0,-2})$ appear identical to $E_1(10_1^{0,-2})$).
   - *Machine-verified exact value*: Yoshikawa omitted the generator $(x + 1)(y + 1)$; integral HNF over $\mathbb{Z}[x^{\pm 1}, y]/(y^2 - 1)$ proves
     $$E_1(10_2^{0,-2}) = \bigl(2x + y - 1,\; (x + 1)(y + 1),\; 4\bigr) \subsetneq (2x + y - 1,\; 4) = E_1(10_1^{0,-2}).$$
     Once restored, $E_1$ alone already separates $10_1^{0,-2}$ from $10_2^{0,-2}$!

*(Note: Once Errata #1–#3 are fixed, the **only** genuine $(E_1, \pi_1, \text{topology})$ collision in Yoshikawa's entire table is $8_1^{1,1}$ vs $10_1^{1,1}$, which our `double_linking_mod2` module cleanly separates via Sanderson's invariant $\operatorname{Dlk} = (0, 0)$ vs $(1, 1)$; meanwhile, the connected-sum subtlety $10_1^1 = \operatorname{Tube}(3_1)$ vs $8_1 \# 2_1^1$ is cleanly separated by our peripheral $\operatorname{SL}(2, \mathbb{F}_3)$ longitude invariant `24` vs `0`!)*

---

## 3. Tier-3 4D Invariants, Symmetry Certificates, & Nearest Separating Witnesses

Beyond classical $\pi_1$ and Fox–Alexander ideals, [`yoshikawa_ch10.json`](../code/data/yoshikawa_ch10.json) records four purely 4-dimensional invariants (proved in [`foundations_and_invariants.md`](./foundations_and_invariants.md), §5) together with each entry's **`cell_status`** (`EXACT_SINGLETON` across all 23 entries) and **`nearest_separating_witness`** under our 10-layer invariant filtration:

1. **Total Mod-2 Double Branched Cover $\Sigma_2(S^4, F)$**: Exact $(b_1,\; \operatorname{Tors} H_1,\; b_2,\; |\sigma|,\; \{b_2^+, b_2^-\})$ via Fox–Plans + Hirzebruch–Atiyah–Singer $\sigma(\Sigma_2) = -\frac{1}{2} e(F)$. Notice how $\Sigma_2$ evaluates at the canonical meridional specialization $(t_i, u_j) = (-1, -1)$, separating $10_1^{0,1}$ ($\Sigma_2$ a $\mathbb{Q}$-homology 4-sphere with $H_1 = \mathbb{Z}_4, b_2 = 0$) from $10_2^{0,1}$ ($b_1(\Sigma_2) = 1, b_2(\Sigma_2) = 2$) even before consulting meridional $E_1$ generators!
2. **Lomonaco $\pi_2(E)$ & Regular Finite Cover Betti Profiles**: Exact Poincaré–Lefschetz Betti numbers $(b_1, \operatorname{Tors} H_1, b_2, b_3)$ on regular finite covers $\widetilde{E}_\rho \to E = S^4 \setminus \nu(F)$, certifying $\operatorname{rk}_{\mathbb{Z}} \pi_2(E)$ whenever $\ker(\rho)$ is trivial or infinite cyclic with trivial Alexander determinant.
3. **Kamada–Oshiro Octahedral Symmetric Quandle $(Q_6, *, \rho)$**: Unoriented coloring count $|\operatorname{Col}_{(Q_6, \rho)}(D)|$ and 3-cocycle state-sum distribution in $H^3_{Q,\rho}(Q_6; \mathbb{Z}_2) \cong \mathbb{Z}_2^{12}$. On $8_1^{-1,-1}$ vs $10_1^{-2,-2}$ (both having $\pi_1 \cong Q_8$ and $24$ colorings), all $24$ colorings of $8_1^{-1,-1}$ carry nontrivial 3-cocycle weight (`0 + 24u`), **machine-certifying $t(8_1^{-1,-1}) > 0$** directly from unlinking movies without trisections or chart moves!
4. **4D Symmetry Certificates**: Constructive $\mathrm{SO}(4)$ Yoshikawa move paths connecting $D$ to its spatial mirror $R_z(D)$ (`X0 <-> X1`) vs algebraic chirality/non-invertibility obstructions from normal Euler numbers and Mochizuki's $R_3$ 3-cocycle. Remarkably, **17 of the 23** Yoshikawa entries are constructively $\mathrm{SO}(4)$-amphicheiral (including $\tau^0(3_1) = 8_1$ and $\tau^3(3_1) = 10_3$, despite $3_1$ being chiral in $S^3$), whereas $2_1^{-1}$ is chiral by $e = \pm 2$ and $10_2 = \tau^2(3_1)$ is both non-invertible ($+10_2 \not\cong -10_2$) and oriented-chiral ($(+10_2)^* \not\cong +10_2$).

| Label | $\Sigma_2(S^4, F)$: $(\chi, b_1, \operatorname{Tors} H_1, b_2, |\sigma|)$ | Certified $\operatorname{rk}_{\mathbb{Z}}\pi_2(E)$ & Key Cover $(b_1, b_2, b_3)$ | $(Q_6, \rho)$ State Sum | 4D Symmetry Certificate | Nearest Neighbor & Separating Layer |
| :--- | :--- | :--- | :---: | :--- | :--- |
| $0_1$ | $(2, 0, 0, 0, 0)$ ($\cong S^4$) | $0$ ($\pi_1 \cong \mathbb{Z}$) | `6 + 0u` | $\mathrm{SO}(4)$-amphicheiral | $8_1$ at Layer 3 (`E1`) |
| $2_1^1$ | $(4, 0, 0, 2, 0)$ ($\cong S^2 \times S^2$) | — ($Z_2$: $(1, 4, 0)$) | `6 + 0u` | $\mathrm{SO}(4)$-amphicheiral | $10_1^1$ at Layer 3 (`E1`) |
| $2_1^{-1}$ | $(3, 0, 0, 1, 1)$ ($\cong \overline{\mathbb{C}P}^2$ / $\mathbb{C}P^2$) | $1$ ($Z_2$: $(0, 1, 0)$) | `0` | **Chiral** ($e = \pm 2$) | $2_1^1$ at Layer 1 (`topology`) |
| $6_1^{0,1}$ | $(2, 0, \mathbb{Z}_2, 0, 0)$ ($\mathbb{Q}HS^4$) | — ($Z_2$: $(2, 2, 1), (3, 4, 2)$) | `12 + 0u` | $\mathrm{SO}(4)$-amphicheiral | $9_1^{0,1}$ at Layer 3 (`E1`) |
| $7_1^{0,-2}$ | $(2, 0, \mathbb{Z}_2, 0, 0)$ ($\mathbb{Q}HS^4$) | — ($Q_8$: $(3, 3, 1)$) | `24 + 0u` | Open (no obstruction) | $10_1^{0,-2}$ at Layer 3 (`E1`) |
| $8_1$ | $(2, 0, \mathbb{Z}_3, 0, 0)$ ($\mathbb{Q}HS^4$) | — ($S_3$: $(1, 2, 2)$) | `6 + 0u` | $\mathrm{SO}(4)$-amphicheiral | $9_1$ at Layer 3 (`E1`) |
| $8_1^{1,1}$ | $(4, 0, \mathbb{Z}_2, 2, 0)$ | — ($Z_2$: $(2, 6, 1)$) | `12 + 0u` | $\mathrm{SO}(4)$-amphicheiral | $10_1^{1,1}$ at **Layer 9 (`double_triple_linking`)** |
| $8_1^{-1,-1}$ | $(2, 0, \mathbb{Z}_2, 0, 0)$ ($\mathbb{Q}HS^4$) | **$0$** ($Q_8$: $(0, 0, 1)$) | **`0 + 24u`** ($t > 0$!) | $\mathrm{SO}(4)$-amphicheiral | $10_1^{-1,-1}$ at Layer 3 (`E1`) |
| $9_1$ | $(2, 0, \mathbb{Z}_3, 0, 0)$ ($\mathbb{Q}HS^4$) | — ($S_3$: $(1, 2, 2)$, tors $\mathbb{Z}_3$) | `6 + 0u` | Open (no obstruction) | $8_1$ at Layer 3 (`E1`) |
| $9_1^{0,1}$ | $(2, 0, \mathbb{Z}_4, 0, 0)$ ($\mathbb{Q}HS^4$) | — ($A_4$: $(7, 12, 6)$) | `36 + 0u` | Open (no obstruction) | $10_1^{0,1}$ at Layer 3 (`E1`) |
| $9_1^{1,-2}$ | $(4, 0, \mathbb{Z}_2, 2, 0)$ | — ($Z_2$: $(1, 5, 1)$) | `0` | Open (no obstruction) | $9_1^{0,1}$ at Layer 1 (`topology`) |
| $10_1$ | $(2, 0, \mathbb{Z}_5, 0, 0)$ ($\mathbb{Q}HS^4$) | — ($A_4$: $(4, 9, 6)$) | `6 + 0u` | $\mathrm{SO}(4)$-amphicheiral | $10_2$ at Layer 3 (`E1`) |
| $10_2$ | $(2, 0, \mathbb{Z}_3, 0, 0)$ ($\mathbb{Q}HS^4$) | **$2$** ($S_3$: $(1, 2, 2)$, tors $0$) | `6 + 0u` | **Oriented-chiral & non-invertible** | $10_1$ at Layer 3 (`E1`) |
| $10_3$ | $(2, 0, 0, 0, 0)$ ($\mathbb{Z}HS^4$!) | **$7$** ($\operatorname{SL}(2, \mathbb{F}_3)$: $(1, 7, 7)$) | `6 + 0u` | $\mathrm{SO}(4)$-amphicheiral | $10_1$ at Layer 3 (`E1`) |
| $10_1^1$ | $(4, 0, \mathbb{Z}_3, 2, 0)$ | — ($\operatorname{SL}(2, \mathbb{F}_3)$: $(4, 54, 3)$) | `6 + 0u` | $\mathrm{SO}(4)$-amphicheiral | $2_1^1$ at Layer 3 (`E1`) |
| $10_1^{0,1}$ | $(2, 0, \mathbb{Z}_4, 0, 0)$ ($\mathbb{Q}HS^4$) | — ($S_3$: $(5, 8, 4)$) | `36 + 0u` | $\mathrm{SO}(4)$-amphicheiral | $10_2^{0,1}$ at **Layer 4 (`Sigma2`)** |
| $10_2^{0,1}$ | $(2, 1, 0, 2, 0)$ | — ($S_3$: $(5, 10, 6)$) | `36 + 0u` | Open (no obstruction) | $10_1^{0,1}$ at **Layer 4 (`Sigma2`)** |
| $10_1^{1,1}$ | $(4, 0, \mathbb{Z}_2, 2, 0)$ | — ($Z_2$: $(2, 6, 1)$) | `12 + 0u` | $\mathrm{SO}(4)$-amphicheiral | $8_1^{1,1}$ at **Layer 9 (`double_triple_linking`)** |
| $10_1^{0,0,1}$ | $(0, 1, \mathbb{Z}_2, 0, 0)$ | — ($S_3$: 5 covers) | `120 + 0u` | $\mathrm{SO}(4)$-amphicheiral | $10_1^{-1,-1}$ at Layer 1 (`topology`) |
| $10_1^{0,-2}$ | $(2, 0, \mathbb{Z}_4, 0, 0)$ ($\mathbb{Q}HS^4$) | — ($\operatorname{Dic}_{16}$: tors $\mathbb{Z}_2^2 \oplus \mathbb{Z}_4$) | `0` | Open (no obstruction) | $10_2^{0,-2}$ at Layer 3 (`E1`) |
| $10_2^{0,-2}$ | $(2, 0, \mathbb{Z}_4, 0, 0)$ ($\mathbb{Q}HS^4$) | — ($\operatorname{Dic}_{16}$: $(3, 5, 3), (3, 7, 5)$) | `0` | $\mathrm{SO}(4)$-amphicheiral | $10_1^{0,-2}$ at Layer 3 (`E1`) |
| $10_1^{-1,-1}$ | $(2, 0, \mathbb{Z}_4, 0, 0)$ ($\mathbb{Q}HS^4$) | **$2$** ($\operatorname{Dic}_{16}$: $(0, 2, 3)$) | `0` | $\mathrm{SO}(4)$-amphicheiral | $8_1^{-1,-1}$ at Layer 3 (`E1`) |
| $10_1^{-2,-2}$ | $(4, 0, \mathbb{Z}_2, 2, 0)$ | **$16$** ($Q_8$: $(0, 16, 1)$) | `24 + 0u` | $\mathrm{SO}(4)$-amphicheiral | $10_1^{-1,-1}$ at Layer 1 (`topology`) |

---

## 4. Canonical Machine-Certified EPD Codes (`ch <= 10`)

Every EPD string below is directly parseable via `MarkedDiagram.from_epd(epd)` (`X[a,b,c,d]` has `(a,c)` as the under-strand and `(b,d)` as the over-strand in CCW order; `Y` resolves `(a,d),(b,c)` in $L_+$ and `(a,b),(c,d)` in $L_-$; `Z` is opposite):

| Label | Canonical EPD |
| :--- | :--- |
| $0_1$ | `O[1]` |
| $2_1^1$ | `Y[1,2,3,4], Z[4,3,2,1]` |
| $2_1^{-1}$ | `Y[1,2,3,4], X[4,3,2,1]` |
| $6_1^{0,1}$ | `Y[1,2,3,4], Z[4,5,6,1], X[6,7,8,9], X[9,8,10,2], X[10,11,12,3], X[11,7,5,12]` |
| $7_1^{0,-2}$ | `Z[1,2,3,4], Y[5,6,7,8], X[6,9,10,11], X[12,13,14,2], X[5,14,13,9], X[8,7,4,3], X[12,1,11,10]` |
| $8_1$ | `X[1,2,3,4], X[4,5,6,1], X[8,6,5,7], X[3,9,10,11], X[11,10,12,13], Z[8,14,15,16], Y[16,15,9,2], X[14,7,13,12]` |
| $8_1^{1,1}$ | `X[1,2,3,4], X[4,5,6,1], Z[6,7,8,9], Y[9,8,10,2], Z[3,11,12,13], Y[13,12,14,5], X[15,16,7,14], X[16,15,11,10]` |
| $8_1^{-1,-1}$ | `X[1,2,3,4], X[4,5,6,1], Y[7,8,9,2], X[10,8,7,6], Y[3,11,12,13], X[13,12,14,5], X[15,16,10,14], X[16,15,11,9]` |
| $9_1$ | `Z[1,2,3,4], Y[4,3,5,6], X[1,8,9,7], X[8,10,11,9], X[12,13,14,15], X[13,16,7,11], X[12,5,2,16], X[17,18,10,6], X[15,14,18,17]` |
| $9_1^{0,1}$ | `X[1,2,3,4], X[4,5,6,1], X[7,8,2,6], X[3,9,10,11], Y[11,12,13,5], Y[13,12,14,15], X[15,16,17,7], X[8,17,16,18], X[14,10,9,18]` |
| $9_1^{1,-2}$ | `Y[1,2,3,4], Z[4,5,6,1], X[6,7,8,9], X[9,8,10,2], X[11,12,13,3], X[12,14,5,13], Z[14,15,16,7], Y[10,17,18,11], X[18,17,16,15]` |
| $10_1$ | `X[1,2,3,4], X[5,6,7,8], Y[5,9,4,10], X[11,12,2,13], Y[3,14,15,16], X[16,17,6,10], X[9,8,13,1], X[17,18,19,7], X[18,15,20,19], X[20,14,12,11]` |
| $10_2$ | `X[1,2,3,4], X[4,5,6,1], X[6,7,8,9], Y[9,10,11,2], Y[12,13,14,15], X[16,13,12,3], X[15,17,7,5], X[18,19,10,8], X[19,20,16,11], X[17,14,20,18]` |
| $10_3$ | `X[1,2,3,4], X[4,5,6,1], X[6,7,8,9], Y[9,10,11,2], X[12,13,14,3], Y[14,15,16,5], X[17,8,7,16], X[18,19,10,17], X[11,20,13,12], X[19,18,15,20]` |
| $10_1^1$ | `Z[1,2,3,4], Y[4,5,6,1], X[6,7,8,9], X[9,10,11,2], X[12,13,14,3], X[15,16,5,14], Z[16,17,18,7], X[18,19,10,8], Y[11,19,20,12], X[13,20,17,15]` |
| $10_1^{0,1}$ | `X[1,2,3,4], X[4,5,6,1], X[7,6,5,8], X[9,3,2,10], X[9,11,12,13], X[13,12,14,15], X[15,14,16,17], Z[7,18,19,20], Y[20,19,11,10], X[18,8,17,16]` |
| $10_2^{0,1}$ | `X[1,2,3,4], X[4,5,6,1], X[6,7,8,9], X[9,10,11,2], X[12,13,10,14], X[15,16,5,17], X[18,8,7,16], Y[18,15,19,14], Z[13,20,3,11], X[12,19,17,20]` |
| $10_1^{1,1}$ | `X[1,2,3,4], X[4,5,6,7], Y[7,8,9,1], X[10,11,2,9], Z[3,11,12,5], Z[6,13,14,15], X[15,16,17,8], Y[17,18,19,10], X[19,20,13,12], X[20,18,16,14]` |
| $10_1^{0,0,1}$ | `Z[1,2,3,4], Y[5,6,2,1], X[7,8,9,10], X[11,7,10,6], X[11,5,12,13], X[14,8,13,12], X[9,15,16,17], X[17,16,18,3], X[18,19,20,4], X[19,15,14,20]` |
| $10_1^{0,-2}$ | `X[1,2,3,4], X[4,5,6,1], X[7,8,2,9], X[8,10,11,3], X[10,12,5,11], X[13,14,15,12], X[16,17,18,19], Y[18,17,6,15], Z[20,9,16,19], X[7,20,14,13]` |
| $10_2^{0,-2}$ | `X[1,2,3,4], X[4,5,6,7], X[7,8,9,1], Y[9,10,11,2], X[11,12,5,3], Y[6,13,14,15], X[15,16,17,8], X[18,19,10,17], X[19,20,13,12], X[20,18,16,14]` |
| $10_1^{-1,-1}$ | `X[1,2,3,4], X[4,5,6,7], Y[7,8,9,1], X[10,11,2,9], X[11,12,5,3], Z[6,13,14,15], X[15,16,17,8], X[18,19,10,17], X[19,20,13,12], X[20,18,16,14]` |
| $10_1^{-2,-2}$ | `Z[1,2,3,4], Y[4,5,6,1], X[6,7,8,9], X[9,8,10,11], X[12,3,2,11], X[13,14,15,12], X[14,16,5,15], Z[16,17,18,7], Y[10,19,20,13], X[20,19,18,17]` |
