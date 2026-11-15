# Parametric Surface-Link Family Compilers & Hub-and-Spoke Cross-Join

Beyond exhaustive bottom-up enumeration by $\operatorname{ch}$-index ($\operatorname{ch} \le 10$ in [`yoshikawa_table.md`](./yoshikawa_table.md) and $\operatorname{ch} = 11$ in [`ch11_expansion.md`](./ch11_expansion.md)), infinite families of 2-knots and surface-links in $\mathbb{R}^4$ are parameterized by classical braids, bridge plats, and free-group ribbon words.

The module [`families.py`](../code/src/surfaceknots/families.py) implements **deterministic, planar, admissibility-certified `MarkedDiagram` compilers** for three major family spokes — exported in [`family_catalog.json`](../code/data/families/family_catalog.json), [`family_catalog.csv`](../code/data/families/family_catalog.csv), and [`hub_realizations.json`](../code/data/families/hub_realizations.json) ($61$ compiled instances across $15$ realized $\operatorname{ch} \le 11$ Hub entries and $41$ pairwise-distinct $\operatorname{ch} \ge 12$ invariant classes):

- **Spoke 1 (Spun & Symmetric-Union 2-Knots)**: Artin's spun 2-knots $\tau^0(K)$ and Eisermann–Lamm $m$-twisted symmetric-union ribbon 2-knots $U_m(K)$ for $b$-bridge classical knots $K$.
- **Spoke 2 (1-Fusion Ribbon 2-Knots)**: General alternating-half-plane planar compiler for 1-fusion ribbon 2-knots $R(w)$ ($w \in F(a, b)$), realizing $\pi_1(S^4 \setminus R(w)) = \langle a, b \mid w^{-1} a w = b \rangle$.
- **Spoke 3 (Satoh Tubes, Mixed Sphere–Torus Links, & Twisted Klein-Bottle Tubes)**: Satoh's ribbon torus-links $\operatorname{Tube}(L)$, mixed $S^2 \sqcup T^2$ links, Brunnian 3-component surface-links, and non-orientable Klein-bottle twisted plat links.

---

## 1. Constructive Diagram Compilers & Complexity Bounds

### 1.1 The Twisted Ribbon 1-Handle Pinch Theorem (`compile_plat_surface`)

Let $L = \widehat{\beta}^{\,\text{plat}} \subset \mathbb{R}^3$ be a classical knot or link presented in $b$-bridge plat position by a braid word $\beta \in B_{2b}$ of length $|\beta|$ on vertical strands $0, 1, \dots, 2b - 1$, closed at the top and bottom by the $b$ bridge caps $(2k, 2k+1)$ ($0 \le k < b$).

Given a subset $S \subseteq \{0, 1, \dots, b - 1\}$ of bridge necks and an integer twist assignment $m_k \in \mathbb{Z}$ for each $k \in S$, `compile_plat_surface` constructs the marked diagram $D(\beta; \{m_k\}_{k \in S})$ by placing $\beta$ in the upper half, $\beta^{-1}$ in the lower half, closing $(2k, 2k+1)$ by standard top and bottom caps, and inserting across columns $(2k, 2k+1)$ at the equator ($k \in S$) a **horizontal `Z-Y` 1-handle gadget** with $|m_k|$ classical crossings (`X0` if $m_k > 0$, `X1` if $m_k < 0$) wired in series along the two horizontal rails between a left `Z`-marker on strand $2k$ and a right `Y`-marker on strand $2k+1$.

**Theorem 1.1 (Planarity, Admissibility, and Constructive $\operatorname{ch}$-Bound).**
*For every $b \ge 1$, every braid word $\beta \in B_{2b}$, and every neck-twist assignment $\{m_k\}_{k \in S}$:*
1. **Strict Planarity & Complexity**: *The 4-regular graph $D = D(\beta; \{m_k\}_{k \in S})$ is strictly planar (`_is_planar_4reg_diagram(D) == True`) with exact vertex complexity*
   $$\operatorname{ch}(D) \;=\; 2|\beta| \;+\; \sum_{k \in S} \bigl(2 + |m_k|\bigr).$$
2. **Unconditional Admissibility**: *Both resolutions $L_+(D)$ and $L_-(D)$ reduce constructively via $|m_k|$ Reidemeister I untwists followed by $|\beta|$ free Reidemeister II cancellations to the $b$-component classical unlink.*
3. **Normal Euler Number & Orientability Parity**: *Every saddle neck contributes zero net resolution-writhe difference ($w(L_+(D)) - w(L_-(D)) = 0$), so every component $F_i \subset F(D)$ has normal Euler number $e(F_i) = 0$. When $m_k$ is even, the $k$-th neck preserves coherent surface orientation (attaching an orientable ribbon 1-handle); when $m_k$ is odd, the $|m_k|$ half-twists reverse local sheet orientation across the 1-handle, yielding non-orientable Klein-bottle ($K^2$, $\chi = 0, e = 0$) components.*

*Proof.* Planarity holds because each braid generator $\sigma_i^{\pm 1}$ acts on adjacent columns $(i, i+1)$ with cyclic CCW port order $(\text{NE} = i+1_{\text{in}}, \text{NW} = i_{\text{in}}, \text{SW} = i_{\text{out}}, \text{SE} = i+1_{\text{out}})$, each neck gadget lies strictly inside the equatorial corridor between columns $2k$ and $2k+1$, and the $b$ top/bottom caps nest disjointly in the outer regions.

At the $k$-th equatorial gadget (`vz = Z` with CCW ports `(0=top_rail, 1=sL_in, 2=sL_out, 3=bot_rail)` and `vy = Y` with CCW ports `(0=sR_in, 1=top_rail, 2=bot_rail, 3=sR_out)`):
- In $L_+(D)$, `Z` pairs `(0,1), (2,3)` and `Y` pairs `(0,3), (1,2)`. Thus strand $2k+1$ passes straight vertically (`vy.0 -- vy.3`), while strand $2k$ enters `vz.1 -- vz.0`, traverses the $|m_k|$ horizontal crossings, reflects around the end cap `vy.1 -- vy.2`, traverses the $|m_k|$ crossings back, and exits `vz.3 -- vz.2`. This forms a cascade of $|m_k|$ innermost Reidemeister I kinks on strand $2k$, collapsing to two parallel vertical strands $(2k, 2k+1)$.
- Symmetrically, in $L_-(D)$, `Z` pairs `(1,2), (0,3)` (vertical pass-through on strand $2k$, end cap `vz.0 -- vz.3`) and `Y` pairs `(0,1), (2,3)`, collapsing the $|m_k|$ kinks onto strand $2k+1$.

Once every equatorial neck collapses to parallel vertical strands, $\beta$ sits directly atop $\beta^{-1}$, canceling letter-by-letter via $|\beta|$ Reidemeister II moves to the $b$ unnested bridge circles. Finally, every crossing in $\beta$ cancels the writhe of its mirror partner in $\beta^{-1}$ in both $L_\pm(D)$, while the $|m_k|$ neck crossings appear as positive/negative self-crossings of a single component in both $L_+(D)$ and $L_-(D)$ with identical signs, so $w(L_+(D)) = w(L_-(D)) = 0$, giving $e(F_i) = 0$ by Kamada's formula. $\blacksquare$

**Corollary 1.2 (Uniform $\operatorname{ch}$-Index Upper Bounds).** Let $c_{\text{plat}}(L)$ denote the minimal crossing number among $b$-bridge plat braid representatives of a classical knot or link $L$:
- **Artin Spin ($S = \{0, \dots, b - 2\},\; m_k = 0$)**: $\operatorname{ch}(\tau^0(K)) \le 2\,c_{\text{plat}}(K) + 2(b - 1)$. For every 2-bridge knot $K$ ($b = 2$, $c_{\text{plat}}(K) = c(K)$ via odd-length Conway continued fractions),
  $$\operatorname{ch}\bigl(\tau^0(K)\bigr) \;\le\; 2\,c(K) + 2.$$
- **Twisted Symmetric-Union 2-Knot ($S = \{0, \dots, b - 2\},\; m_0 = m,\; m_{k>0} = 0$)**: $\operatorname{ch}(U_m(K)) \le 2\,c_{\text{plat}}(K) + 2(b - 1) + |m|$; for 2-bridge $K$ with $m = \pm 1$, $\operatorname{ch}(U_{\pm 1}(K)) \le 2\,c(K) + 3$.
- **Satoh Ribbon Torus-Link ($S = \{0, \dots, b - 1\},\; m_k = 0$)**: $\operatorname{ch}(\operatorname{Tube}(L)) \le 2\,c_{\text{plat}}(L) + 2b$; for 2-bridge $L$, $\operatorname{ch}(\operatorname{Tube}(L)) \le 2\,c(L) + 4$.

---

### 1.2 The NORTH/SOUTH Planar 1-Fusion Ribbon Word Compiler (`compile_1fusion_ribbon_word`)

Every 1-fusion ribbon 2-knot $R(w) \subset \mathbb{R}^4$ is obtained by taking a 2-component trivial 2-link $O_a \sqcup O_b = \partial(D_a \sqcup D_b)$ with meridian generators $a, b$ and fusing $O_a$ to $O_b$ along a single ribbon 1-handle whose core follows a reduced word $w = g_1^{\varepsilon_1} \cdots g_\ell^{\varepsilon_\ell} \in F(a, b)$ ($g_j \in \{a, b\}$, $\varepsilon_j \in \{\pm 1\}$), piercing the ribbon disks $D_{g_j}$ in transverse ribbon intersections.

Naive single-half-plane routing of consecutive powers such as $b^2$ wraps a spiral around the outer arc of $O_b$, trapping the fusion attachment point inside a closed planar region and violating planarity on subsequent $a$-letters. Our compiler `compile_1fusion_ribbon_word` eliminates all topological trapping via **NORTH/SOUTH half-plane alternation**:

**Theorem 1.3 (Planar $4|w| + 2$ Compiler & Fox Polynomial Formula).**
*For every word $w \in F(a, b)$ of syllable length $|w| = \ell$, routing successive ribbon piers alternately through the upper (`NORTH`) and lower (`SOUTH`) exterior half-planes across the West tip of $O_a$ (`insert(0)`) and East tip of $O_b$ (`append`) produces a strictly planar, admissible marked diagram $D_w$ satisfying:*
1. $\operatorname{ch}(D_w) = 4|w| + 2$ ($4|w|$ classical crossings and $2$ opposite-type saddle markers `Z, Y`), representing an oriented 2-sphere $R(w) \cong S^2$.
2. $\pi_1(S^4 \setminus R(w)) \cong \langle a, b \mid w^{-1} a w = b \rangle$.
3. *By Fox free differential calculus on $r = w^{-1} a w b^{-1}$, the single-variable Alexander polynomial of $R(w)$ is given up to $\pm t^{\pm k}$ by*
   $$\Delta_{R(w)}(t) \;\doteq\; 1 \;+\; (t - 1)\left.\frac{\partial w}{\partial a}\right|_{a = b = t}.$$

*Proof.* In broken-surface movie coordinates, a 2-rail ribbon band pierces the disk $D_g$ with sign $\varepsilon \in \{\pm 1\}$ whenever both rails cross one boundary arc of $O_g$ **over** ($\varepsilon = +1$: `X1` on horizontal entry) or **under** ($\varepsilon = -1$: `X0` on horizontal entry) and cross the opposite boundary arc of $O_g$ with the opposite over/under sense (`X0` vs `X1`). Crucially, the 4D ribbon intersection sign depends *only* on the over/under sequence of the two boundary crossings — **not** on whether the band approaches $O_g$ from `NORTH` (top-to-bottom) or `SOUTH` (bottom-to-top). Alternating `NORTH` and `SOUTH` after each letter places every new crossing pair on the outermost boundary of $O_a$ (West) or $O_b$ (East), keeping both `cur_west, cur_east` and the terminal `Y`-attachment exposed on the infinite exterior face at every step.

Differentiating $r = w^{-1} a w b^{-1}$ in $\mathbb{Z}[F(a, b)]$ via $\partial(uv)/\partial a = \partial u/\partial a + u\,\partial v/\partial a$ and $\partial(w^{-1})/\partial a = -w^{-1}\,\partial w/\partial a$ yields
$$\frac{\partial r}{\partial a} \;=\; -w^{-1}\frac{\partial w}{\partial a}\,a w b^{-1} + w^{-1} + w^{-1} a \frac{\partial w}{\partial a} b^{-1}.$$
Under abelianization $\mathfrak{a}\colon \mathbb{Z}[F(a, b)] \to \mathbb{Z}[t^{\pm 1}]$ ($a, b \mapsto t$, so $\mathfrak{a}(w) = t^{\sigma(w)}$ cancels $\mathfrak{a}(w^{-1}) = t^{-\sigma(w)}$), this simplifies to
$$\mathfrak{a}\!\left(\frac{\partial r}{\partial a}\right) \;=\; t^{-1}\!\left(1 + (t - 1)\left.\frac{\partial w}{\partial a}\right|_{a=b=t}\right), \qquad \mathfrak{a}\!\left(\frac{\partial r}{\partial b}\right) \;=\; -\mathfrak{a}\!\left(\frac{\partial r}{\partial a}\right)$$
(since $(a - 1)\partial_a r + (b - 1)\partial_b r = r - 1 \mapsto 0$), so $E_1(R(w))$ is principal and generated by $1 + (t - 1)\partial_a w|_{a=b=t}$. $\blacksquare$

---

## 2. Hub Cross-Join Identifications ($\operatorname{ch} \le 11$)

Cross-joining all $61$ compiled family instances against the exhaustive $\operatorname{ch} \le 11$ Hub (`23` Yoshikawa $\operatorname{ch} \le 10$ classes + `10` new $\operatorname{ch} = 11$ classes) geometrically identifies **$15$ distinct Hub classes** across $20$ family parameterizations — including **`MOVE_CERTIFIED` constructive Yoshikawa move equivalences for $12$ Hub classes** (notably identifying the new $\operatorname{ch} = 11$ classes $11_1$, $11_2$, and $11_3^{0,-2}$!):

| Hub Label | $\operatorname{ch}$ | Topology | Realizing Spoke Instance(s) | Match Mode | Geometric Interpretation |
| :--- | :---: | :---: | :--- | :--- | :--- |
| $0_1$ | $0$ | $S^2$ | `Spin(0_1)`, `R(1)` | `MOVE_CERTIFIED` | Unknotted 2-sphere |
| $2_1^1$ | $2$ | $T^2$ | `Tube(0_1)` | `MOVE_CERTIFIED` | Unknotted Satoh torus $\operatorname{Tube}(\text{unknot})$ |
| $6_1^{0,1}$ | $6$ | $S^2 \sqcup T^2$ | `Mixed(L2a1)` | `MOVE_CERTIFIED` | Mixed sphere–torus Hopf link |
| $7_1^{0,-2}$ | $7$ | $S^2 \sqcup K^2$ | `KleinMixed(L2a1,1)` | `MOVE_CERTIFIED` | Half-twisted Klein–sphere Hopf link |
| $8_1$ | $8$ | $S^2$ | `Spin(3_1)`, `R(b*a)` | `MOVE_CERTIFIED` | Artin spun trefoil $\tau^0(3_1) \cong R(ba)$ |
| $8_1^{1,1}$ | $8$ | $T^2 \sqcup T^2$ | `Tube(L2a1)` | `MOVE_CERTIFIED` | Satoh tubed Hopf torus-link $\operatorname{Tube}(L2a1)$ |
| $9_1$ | $9$ | $S^2$ | `SymUnion(3_1,1)` (`MOVE`), `R(b*a^-1)` (`INV`) | `MOVE_CERTIFIED` | 1-twisted symmetric union $U_1(3_1) \cong R(ba^{-1})$ |
| $9_1^{1,-2}$ | $9$ | $T^2 \sqcup K^2$ | `KleinMixed(L2a1,10)` | `MOVE_CERTIFIED` | Torus–Klein-bottle Hopf link |
| $10_1$ | $10$ | $S^2$ | `Spin(4_1)` ($c=10$), `R(b^-1*a*b*a^-1)` ($c=18$) | `INVARIANT_MATCH` | Artin spun figure-eight $\tau^0(4_1) \cong R([b, a^{-1}])$ (`EXACT_SINGLETON` Hub cell) |
| $10_1^1$ | $10$ | $T^2$ | `Tube(3_1)` ($c=10$) | `INVARIANT_MATCH` | Satoh spun-torus $\operatorname{Tube}(3_1)$ (`EXACT_SINGLETON` Hub cell) |
| $10_1^{0,1}$ | $10$ | $S^2 \sqcup T^2$ | `Mixed(L4a1)` | `MOVE_CERTIFIED` | Mixed sphere–torus Solomon's link |
| $10_1^{-2,-2}$ | $10$ | $K^2 \sqcup K^2$ | `KleinTube(L2a1,11)` | `MOVE_CERTIFIED` | Double Klein-bottle Hopf link (via C mountain-pass merge) |
| **$11_1$** | **$11$** | $S^2$ | `SymUnion(4_1,1)` (`MOVE`, $c=11$), `R(b^-1*a*b^-1*a^-1)` (`INV`) | **`MOVE_CERTIFIED`** | **1-twisted symmetric union $U_1(4_1)$ of the figure-eight knot!** |
| **$11_2$** | **$11$** | $S^2$ | `R(b^-1*a*b*a)` ($c=18$) | `INVARIANT_MATCH` | **Asymmetric 4-syllable 1-fusion ribbon 2-knot $R(b^{-1}aba)$** (`EXACT_SINGLETON`) |
| **$11_3^{0,-2}$** | **$11$** | $S^2 \sqcup K^2$ | `KleinMixed(L4a1,1)` ($c=11$) | **`MOVE_CERTIFIED`** | **Half-twisted Klein–sphere Solomon surface-link!** |

### Residual Non-Elementary Census Entries ($\operatorname{ch} \le 11$)
Following the Batch Geometric Identification protocol in [`roadmap.md`](./roadmap.md), subtracting the $15$ Spoke-realized Hub entries highlights the **genuinely 4-dimensional / twist-spun / non-ribbon** prime entries up to $\operatorname{ch} \le 11$:
- **Non-zero normal Euler number ($e(F_i) \neq 0$, obstructed from all untwisted ribbon/plat families by Theorem 1.1(3))**: $2_1^{-1}$, $8_1^{-1,-1}$, $10_1^{-1,-1}$.
- **Positive Zeeman twist-spins ($\tau^m(K)$ with $m \ge 2$, non-ribbon with nontrivial 3-cocycles / $\pi_2$ ranks)**: $10_2 \cong \tau^2(3_1)$, $10_3 \cong \tau^3(3_1)$, $10_1^{1,1}$ (1-twist-spun Hopf link with $\operatorname{Dlk} \neq 0$).
- **Higher-fusion / twisted torus-sphere & non-orientable links**: $9_1^{0,1}$, $10_2^{0,1}$, $10_1^{0,0,1}$, $10_1^{0,-2}$, $10_2^{0,-2}$ (non-ribbon candidate), and at $\operatorname{ch} = 11$: the first prime knotted Klein bottle **$11_1^{-2}$**, the links $11_1^{0,1}, 11_2^{0,1}, 11_1^{0,-2}, 11_2^{0,-2}, 11_1^{1,-2}$, and the 3-component non-orientable link **$11_1^{0,0,-2}$**.

---

## 3. Certified Exact $\operatorname{ch}(F) = 12$ Determinations

Combining Theorem 1.1's constructive upper bounds with the exhaustiveness of our $\operatorname{ch} \le 11$ prime census, $\operatorname{ch} = 11$ residual shadow log, and $\operatorname{ch} \le 12$ composite/split/unknot reference catalog yields **five unconditional $\operatorname{ch}(F) = 12$ theorems** without enumerating the $\sim 1.5 \times 10^4$ prime shadows at $n = 12$:

**Theorem 3.1 (Five Exact $\operatorname{ch}(F) = 12$ Surface-Links).**
*The following five weakly prime surface-knots and surface-links have **exact** Yoshikawa ch-index $\operatorname{ch}(F) = 12$:*
1. **$\tau^0(5_1)$ (`Spin(5_1)`)**: Artin spun $(2,5)$-torus knot ($S^2$, $\Delta(t) = t^4 - t^3 + t^2 - t + 1$, $H_1(\Sigma_2) = \mathbb{Z}_5$).
2. **$\tau^0(5_2)$ (`Spin(5_2)`)**: Artin spun three-twist knot ($S^2$, $\Delta(t) = 2t^2 - 3t + 2$, $H_1(\Sigma_2) = \mathbb{Z}_7$).
3. **$\operatorname{Tube}(4_1)$ (`Tube(4_1)`)**: Satoh ribbon torus of the figure-eight knot ($T^2$, $\Delta(t) = t^2 - 3t + 1$, $H_1(\Sigma_2) = \mathbb{Z}_5, b_2(\Sigma_2) = 2$).
4. **$\operatorname{Tube}(L4a1)$ (`Tube(L4a1)`)**: Satoh 2-torus link of Solomon's $(2,4)$-torus link ($T^2 \sqcup T^2$, $E_1 = (x^2 - y^2,\; xy - x + y^2 - y)$, $H_1(\Sigma_2) = \mathbb{Z}_4, b_2(\Sigma_2) = 2$).
5. **$\operatorname{Mixed}(L5a1)$ (`Mixed(L5a1)`)**: Mixed $S^2 \sqcup T^2$ Whitehead link ($E_1 = ((x - 1)(y - 1)(xy - x - y + \dots)) = (x^2y - x^2 - 2xy + 2x + y - 1,\; xy^2 - 2xy + x - y^2 + 2y - 1)$, $H_1(\Sigma_2) = \mathbb{Z}_8$).

*Proof.* For each of the five surfaces $F$, `compile_plat_surface` constructs an explicit planar, Reidemeister-unlink-certified marked diagram $D$ with $\operatorname{ch}(D) = 12$ (`Spin(5_1)` and `Spin(5_2)` have $2 \cdot 5 + 2 = 12$; `Tube(4_1)` and `Tube(L4a1)` have $2 \cdot 4 + 4 = 12$; `Mixed(L5a1)` has $2 \cdot 5 + 2 = 12$), establishing $\operatorname{ch}(F) \le 12$.

Conversely, suppose $\operatorname{ch}(F) \le 11$, so $F$ admits some minimal diagram $D_{\min}$ with $\operatorname{ch}(D_{\min}) \le 11$.
- $D_{\min}$ cannot be disconnected/split: for $c = 1$ (`Spin(5_1)`, `Spin(5_2)`, `Tube(4_1)`) $F$ is connected; for `Tube(L4a1)` and `Mixed(L5a1)`, $E_1(F) \neq (0)$, whereas every split 2-component orientable link satisfies $E_1 = (0)$.
- $D_{\min}$ cannot be a nontrivial connected sum $F_1 \# F_2$ with $\operatorname{ch}(F_1) + \operatorname{ch}(F_2) \le 12$: our exhaustive reference catalog `build_reference_catalog(max_composite_ch=12)` enumerates every connected sum of $\operatorname{ch} \le 10$ prime classes and standard surfaces up to total $\operatorname{ch} \le 12$, and none matches the 10-layer invariant key $I(F)$!
  > **Subtlety on $\operatorname{Tube}(4_1)$ vs $10_1 \# 2_1^1$**: Both $\operatorname{Tube}(4_1)$ and the $\operatorname{ch} = 12$ composite torus $\tau^0(4_1) \# T^2_{\text{unknotted}} = 10_1 \# 2_1^1$ have homeomorphism type $T^2$, identical group $\pi_1 \cong \pi_1(S^3 \setminus 4_1)$, identical $E_1 = (x^2 - 3x + 1)$, identical $\Sigma_2$ invariants $(b_1 = 0, H_1 = \mathbb{Z}_5, b_2 = 2, \sigma = 0)$, and identical kei, bikei, and Alexander biquandle invariants! However, our Tier-3 Poincaré–Lefschetz regular $A_4$-cover Betti invariant (Layer 5) and peripheral $\operatorname{SL}(2, \mathbb{F}_3)$ holonomy (Layer 10) separate them cleanly (`test_tube_4_1_separated_from_composite_10_1_conn_2_1_1` in [`test_families.py`](../code/tests/test_families.py)):
  > - **$A_4$-cover Betti profile $(b_1, \operatorname{Tors} H_1, b_2, b_3)$**: $(4, \mathbb{Z}_2^2, 30, 3)$ on $\operatorname{Tube}(4_1)$ vs $(4, \mathbb{Z}_2^2, 33, 6)$ on $10_1 \# 2_1^1$!
  > - **Non-meridional peripheral $\operatorname{SL}(2, \mathbb{F}_3)$ representations**: $24$ on $\operatorname{Tube}(4_1)$ vs $0$ on $10_1 \# 2_1^1$!
- Therefore a minimal $\operatorname{ch} \le 11$ diagram of $F$ would have to appear in our exhaustive $\operatorname{ch} \le 10$ Yoshikawa catalog, our $\operatorname{ch} = 11$ prime catalog (`10` classes), or our $\operatorname{ch} = 11$ `AMBIGUOUS_SHADOW` residual log (`93` diagrams). Machine comparison confirms $I(F)$ is disjoint from all three sets, proving $\operatorname{ch}(F) \ge 12$ and hence $\operatorname{ch}(F) = 12$. $\blacksquare$

---

## 4. Complete Catalog Summary of Spokes 1–3 (`family_catalog.csv`)

### 4.1 Spoke 1: Artin Spun 2-Knots $\tau^0(K)$ & Symmetric-Union 2-Knots $U_1(K)$ ($22$ instances)

| `family_id` | Classical Base | $[\operatorname{ch}_{\min}, \operatorname{ch}_{\max}]$ | Exact? | Hub Match | $E_1(F)$ Generator | $\Sigma_2$: $\operatorname{Tors} H_1$ | Nearest Neighbor (Separating Layer) |
| :--- | :--- | :---: | :---: | :---: | :--- | :---: | :--- |
| `Spin(0_1)` | $0_1$ ($0/1$) | $[0, 0]$ | **Yes** | $0_1$ (`MOVE`) | $1$ | $0$ | $8_1$ (`E1`) |
| `Spin(3_1)` | $3_1$ ($3/1$) | $[8, 8]$ | **Yes** | $8_1$ (`MOVE`) | $x^2 - x + 1$ | $\mathbb{Z}_3$ | $9_1$ (`E1`) |
| `Spin(4_1)` | $4_1$ ($5/3$) | $[10, 10]$ | **Yes** | $10_1$ (`INV`) | $x^2 - 3x + 1$ | $\mathbb{Z}_5$ | $10_2$ (`E1`) |
| **`Spin(5_1)`** | $5_1$ ($5/1$) | **$[12, 12]$** | **Yes** | — | $x^4 - x^3 + x^2 - x + 1$ | $\mathbb{Z}_5$ | `Spin(5_2)` (`E1`) |
| **`Spin(5_2)`** | $5_2$ ($7/4$) | **$[12, 12]$** | **Yes** | — | $2x^2 - 3x + 2$ | $\mathbb{Z}_7$ | `Spin(5_1)` (`E1`) |
| `Spin(6_1)` | $6_1$ ($9/5$) | $[12, 14]$ | No | — | $2x^2 - 5x + 2$ | $\mathbb{Z}_9$ | `R(b^2*a)` (`E1`) |
| `Spin(6_2)` | $6_2$ ($11/7$) | $[12, 14]$ | No | — | $x^4 - 3x^3 + 3x^2 - 3x + 1$ | $\mathbb{Z}_{11}$ | `R(b^2*a)` (`E1`) |
| `Spin(6_3)` | $6_3$ ($13/8$) | $[12, 14]$ | No | — | $x^4 - 3x^3 + 5x^2 - 3x + 1$ | $\mathbb{Z}_{13}$ | `R(b^2*a)` (`E1`) |
| `Spin(7_1)`–`Spin(7_7)` | $7_1 \dots 7_7$ | $[12, 16]$ | No | — | $\Delta_{7_k}(x)$ (all 7 distinct) | $\mathbb{Z}_7 \dots \mathbb{Z}_{21}$ | Pairwise separated at `E1` |
| `SymUnion(3_1,1)` | $3_1$, $m=1$ | $[9, 9]$ | **Yes** | $9_1$ (`MOVE`) | $x - 2$ | $\mathbb{Z}_3$ | $8_1$ (`E1`) |
| **`SymUnion(4_1,1)`** | $4_1$, $m=1$ | **$[11, 11]$** | **Yes** | **$11_1$ (`MOVE`)** | $2x^2 - 2x + 1$ | $\mathbb{Z}_5$ | $11_2$ (`E1`) |
| `SymUnion(5_1,1)` | $5_1$, $m=1$ | $[12, 13]$ | No | — | $2x - 3$ | $\mathbb{Z}_5$ | `SymUnion(5_2,1)` (`E1`) |
| `SymUnion(5_2,1)` | $5_2$, $m=1$ | $[12, 13]$ | No | — | $2x^2 - 4x + 1$ | $\mathbb{Z}_7$ | `SymUnion(5_1,1)` (`E1`) |
| `SymUnion(6_1,1)` | $6_1$, $m=1$ | $[12, 15]$ | No | — | $3x^2 - 4x + 2$ | $\mathbb{Z}_9$ | `SymUnion(6_2,1)` (`E1`) |
| `SymUnion(6_2,1)` | $6_2$, $m=1$ | $[12, 15]$ | No | — | $x^3 - 4x^2 + 4x - 2$ | $\mathbb{Z}_{11}$ | `SymUnion(6_1,1)` (`E1`) |
| `SymUnion(6_3,1)` | $6_3$, $m=1$ | $[12, 15]$ | No | — | $2x^3 - 6x^2 + 4x - 1$ | $\mathbb{Z}_{13}$ | `SymUnion(6_1,1)` (`E1`) |

### 4.2 Spoke 2: 1-Fusion Ribbon 2-Knots $R(w)$ ($14$ instances)

| `family_id` | Word $w \in F(a, b)$ | $[\operatorname{ch}_{\min}, \operatorname{ch}_{\max}]$ | Hub Match | Fox $\Delta_{R(w)}(t)$ | $\Sigma_2$: $\operatorname{Tors} H_1$ | Nearest Neighbor (Separating Layer) |
| :--- | :--- | :---: | :---: | :--- | :---: | :--- |
| `R(1)` | $1$ | $[0, 0]$ | $0_1$ (`MOVE`) | $1$ | $0$ | $8_1$ (`E1`) |
| `R(b*a)` | $ba$ | $[8, 8]$ | $8_1$ (`MOVE`) | $t^2 - t + 1$ | $\mathbb{Z}_3$ | $9_1$ (`E1`) |
| `R(b*a^-1)` | $ba^{-1}$ | $[9, 10]$ | $9_1$ (`INV`) | $t - 2$ | $\mathbb{Z}_3$ | $8_1$ (`E1`) |
| `R(b^-1*a*b*a^-1)` | $b^{-1}aba^{-1}$ | $[10, 18]$ | $10_1$ (`INV`) | $t^2 - 3t + 1$ | $\mathbb{Z}_5$ | $10_2$ (`E1`) |
| `R(b^-1*a*b^-1*a^-1)` | $b^{-1}ab^{-1}a^{-1}$ | $[11, 18]$ | **$11_1$ (`INV`)** | $t^2 - 2t + 2$ | $\mathbb{Z}_5$ | $11_2$ (`E1`) |
| `R(b^-1*a*b*a)` | $b^{-1}aba$ | $[11, 18]$ | **$11_2$ (`INV`)** | $t^3 - t^2 + 2t - 1$ | $\mathbb{Z}_5$ | $11_1$ (`E1`) |
| `R(b^2*a)` | $b^2 a$ | $[12, 14]$ | — | $t^3 - t + 1$ | $0$ ($\mathbb{Z}HS^4$) | `R(b^2*a^-2)` (`kei_bikei_abq`) |
| `R(b^2*a^-1)` | $b^2 a^{-1}$ | $[12, 14]$ | — | $t^2 - t - 1$ | $0$ ($\mathbb{Z}HS^4$) | $10_3$ (`pi2_covers`) |
| `R(b^3*a)` | $b^3 a$ | $[12, 18]$ | — | $t^4 - t + 1$ | $\mathbb{Z}_3$ | $10_2$ (`pi2_covers`) |
| `R(b^3*a^-1)` | $b^3 a^{-1}$ | $[12, 18]$ | — | $t^3 - t^2 - 1$ | $\mathbb{Z}_3$ | `R(b^3*a)` (`finite_quotients`) |
| `R(b^2*a^2)` | $b^2 a^2$ | $[12, 18]$ | — | $t^4 - t^2 + 1$ | $0$ ($\mathbb{Z}HS^4$) | `R(b^2*a^-1)` (`finite_quotients`) |
| `R(b^2*a^-2)` | $b^2 a^{-2}$ | $[12, 18]$ | — | $t^2 - 2$ | $0$ ($\mathbb{Z}HS^4$) | `R(b^2*a)` (`kei_bikei_abq`) |
| `R(b^-1*a^2*b*a^-1)` | $b^{-1}a^2 ba^{-1}$ | $[12, 22]$ | — | $t^3 - 2t^2 - t + 1$ | $0$ ($\mathbb{Z}HS^4$) | `R(b^2*a*b^-1*a)` (`E1`) |
| `R(b^2*a*b^-1*a)` | $b^2 ab^{-1}a$ | $[12, 22]$ | — | $t^3 - 2t^2 + 2$ | $\mathbb{Z}_3$ | `R(b^-1*a^2*b*a^-1)` (`E1`) |

> [!NOTE]
> Notice the **integer homology 4-sphere branched covers** ($\Sigma_2(S^4, R(w))$ with $H_1 = 0, b_2 = 0$) whenever $|\Delta_{R(w)}(-1)| = 1$ (`R(b^2*a)`, `R(b^2*a^-1)`, `R(b^2*a^2)`, `R(b^2*a^-2)`, `R(b^-1*a^2*b*a^-1)`), and notice how non-reciprocal words with trivial mod-2 torsion (`R(b^2*a)` vs `R(b^2*a^-2)`) are still separated cleanly by our universal bikei / oriented Alexander biquandle invariants (`kei_bikei_abq`) and finite-quotient / cover Betti profiles!

### 4.3 Spoke 3: Satoh Tubes, Mixed Sphere–Torus Links, & Klein-Bottle Twisted Plats ($25$ instances)

| `family_id` | Topology | $[\operatorname{ch}_{\min}, \operatorname{ch}_{\max}]$ | Exact? | Hub Match | $\Sigma_2$: $(\operatorname{Tors} H_1, b_2)$ | Periph. $\operatorname{SL}(2, \mathbb{F}_3)$ | Highlights |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| `Tube(0_1)` | $T^2$ | $[2, 2]$ | **Yes** | $2_1^1$ (`MOVE`) | $(0, 2)$ | $0$ | Unknotted Satoh torus |
| `Tube(3_1)` | $T^2$ | $[10, 10]$ | **Yes** | $10_1^1$ (`INV`) | $(\mathbb{Z}_3, 2)$ | $24$ | Separated from $8_1 \# 2_1^1$ by peripheral $\operatorname{SL}(2, \mathbb{F}_3)$ |
| **`Tube(4_1)`** | $T^2$ | **$[12, 12]$** | **Yes** | — | $(\mathbb{Z}_5, 2)$ | $24$ | **Exact $\operatorname{ch} = 12$**; separated from $10_1 \# 2_1^1$ by $A_4$-cover $b_2=30$ vs $33$ & periph. $\operatorname{SL}(2, \mathbb{F}_3)$ |
| `Tube(5_1)` | $T^2$ | $[12, 14]$ | No | — | $(\mathbb{Z}_5, 2)$ | $0$ | Separated from `Spin(5_1) # 2_1^1` by $D_{10}$-cover Betti |
| `Tube(5_2)`, `Tube(6_1)` | $T^2$ | $[12, 14\text{–}16]$ | No | — | $(\mathbb{Z}_7, 2), (\mathbb{Z}_9, 2)$ | $0$ | Ribbon spun tori of twist knots |
| `Tube(L2a1)` | $T^2 \sqcup T^2$ | $[8, 8]$ | **Yes** | $8_1^{1,1}$ (`MOVE`) | $(\mathbb{Z}_2, 2)$ | $122$ | Satoh Hopf 2-torus link ($\operatorname{Dlk} = 0$) |
| **`Tube(L4a1)`** | $T^2 \sqcup T^2$ | **$[12, 12]$** | **Yes** | — | $(\mathbb{Z}_4, 2)$ | $148$ | **Exact $\operatorname{ch} = 12$** Solomon 2-torus link |
| `Tube(L5a1)` | $T^2 \sqcup T^2$ | $[12, 14]$ | No | — | $(\mathbb{Z}_8, 2)$ | $96$ | Linking-number-zero Whitehead 2-torus link |
| `Tube(L6a1)`–`Tube(L6a3)` | $T^2 \sqcup T^2$ | $[12, 16]$ | No | — | $(\mathbb{Z}_{12}, 2), (\mathbb{Z}_{10}, 2), (\mathbb{Z}_6, 2)$ | $148, 138, 162$ | All three 6-crossing 2-bridge torus links pairwise separated at `E1` |
| **`Tube(L6a4)`** | $T^2 \sqcup T^2 \sqcup T^2$ | $[12, 22]$ | No | — | $(\mathbb{Z}_4^2, 2)$ | $2664$ | **Brunnian 3-torus Borromean surface-link** (`is_brunnian = True`) |
| `Mixed(L2a1)` | $S^2 \sqcup T^2$ | $[6, 6]$ | **Yes** | $6_1^{0,1}$ (`MOVE`) | $(\mathbb{Z}_2, 0)$ | $69$ | Mixed Hopf sphere–torus link |
| `Mixed(L4a1)` | $S^2 \sqcup T^2$ | $[10, 10]$ | **Yes** | $10_1^{0,1}$ (`MOVE`) | $(\mathbb{Z}_4, 0)$ | $86$ | Mixed Solomon sphere–torus link |
| **`Mixed(L5a1)`** | $S^2 \sqcup T^2$ | **$[12, 12]$** | **Yes** | — | $(\mathbb{Z}_8, 0)$ | $48$ | **Exact $\operatorname{ch} = 12$** mixed Whitehead link |
| `Mixed(L6a1)`–`Mixed(L6a3)` | $S^2 \sqcup T^2$ | $[12, 14]$ | No | — | $(\mathbb{Z}_{12}, 0), (\mathbb{Z}_{10}, 0), (\mathbb{Z}_6, 0)$ | $86, 69, 93$ | Mixed 6-crossing sphere–torus links |
| **`Mixed(L6a4,01)`** | $S^2 \sqcup T^2 \sqcup T^2$ | $[12, 16]$ | No | — | $(\mathbb{Z}_4^2, 0)$ | $1776$ | **Brunnian $S^2 \sqcup T^2 \sqcup T^2$ Borromean link** reduced from `raw = 20` to `ch = 16`! |
| `KleinMixed(L2a1,1)` | $S^2 \sqcup K^2$ | $[7, 7]$ | **Yes** | $7_1^{0,-2}$ (`MOVE`) | $(\mathbb{Z}_2, 0)$ | $0$ | Half-twisted Hopf Klein–sphere link |
| `KleinMixed(L2a1,10)` | $T^2 \sqcup K^2$ | $[9, 9]$ | **Yes** | $9_1^{1,-2}$ (`MOVE`) | $(\mathbb{Z}_2, 2)$ | $9$ | Hopf torus–Klein-bottle link |
| `KleinTube(L2a1,11)` | $K^2 \sqcup K^2$ | $[10, 10]$ | **Yes** | $10_1^{-2,-2}$ (`MOVE`) | $(\mathbb{Z}_2, 2)$ | $0$ | Double Klein-bottle Hopf link ($\operatorname{rk}_{\mathbb{Z}}\pi_2 = 16$) |
| **`KleinMixed(L4a1,1)`** | $S^2 \sqcup K^2$ | **$[11, 11]$** | **Yes** | **$11_3^{0,-2}$ (`MOVE`)** | $(\mathbb{Z}_4, 0)$ | $0$ | Geometrically identifies $\operatorname{ch} = 11$ class $11_3^{0,-2}$! |
| `KleinMixed(L5a1,1)` | $S^2 \sqcup K^2$ | $[12, 13]$ | No | — | $(\mathbb{Z}_8, 0)$ | $0$ | Half-twisted Whitehead Klein–sphere link |
