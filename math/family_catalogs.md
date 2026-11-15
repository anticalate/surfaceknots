# Parametric Surface-Link Family Compilers & Hub-and-Spoke Cross-Join

Beyond exhaustive bottom-up enumeration by $\operatorname{ch}$-index ($\operatorname{ch} \le 10$ in [`yoshikawa_table.md`](./yoshikawa_table.md), $\operatorname{ch} = 11$ in [`ch11_expansion.md`](./ch11_expansion.md), and $\operatorname{ch} = 12$ in [`ch12_expansion.md`](./ch12_expansion.md)), infinite families of 2-knots and surface-links in $\mathbb{R}^4$ are parameterized by classical braids, bridge plats, and free-group ribbon words.

The module [`families.py`](../code/src/surfaceknots/families.py) implements **deterministic, planar, admissibility-certified `MarkedDiagram` compilers** for three major family spokes — exported in [`family_catalog.json`](../code/data/families/family_catalog.json), [`family_catalog.csv`](../code/data/families/family_catalog.csv), and [`hub_realizations.json`](../code/data/families/hub_realizations.json) (**`64` compiled instances** across `22` Spoke 1, `17` Spoke 2, and `25` Spoke 3 parameterizations: `20` instances realizing `15` $\operatorname{ch} \le 11$ Hub entries plus `44` pairwise-distinct $\operatorname{ch} \ge 12$ invariant classes, of which **`7` promote to $\operatorname{ch} = 12$ Hub classes** in [`expanded_ch12.json`](../code/data/expanded_ch12.json) — `5` `MOVE_CERTIFIED` and `2` `INVARIANT_MATCH` — leaving **`37` distinct $\operatorname{ch} \ge 13$ Spoke classes**, including **`4` certified at exact $\operatorname{ch}(F) = 13$**):

- **Spoke 1 (Spun & Symmetric-Union 2-Knots — `22` instances)**: Artin's spun 2-knots $\tau^0(K)$ and Eisermann–Lamm $m$-twisted symmetric-union ribbon 2-knots $U_m(K)$ for $b$-bridge classical knots $K$.
- **Spoke 2 (1-Fusion Ribbon 2-Knots — `17` instances)**: General alternating-half-plane planar compiler for 1-fusion ribbon 2-knots $R(w)$ ($w \in F(a, b)$), realizing $\pi_1(S^4 \setminus R(w)) = \langle a, b \mid w^{-1} a w = b \rangle$.
- **Spoke 3 (Satoh Tubes, Mixed Sphere–Torus Links, & Twisted Klein-Bottle Tubes — `25` instances)**: Satoh's ribbon torus-links $\operatorname{Tube}(L)$, mixed $S^2 \sqcup T^2$ links, Brunnian 3-component surface-links, and non-orientable Klein-bottle twisted plat links.

---

## 1. Constructive Diagram Compilers & Complexity Bounds

### 1.1 The Twisted Ribbon 1-Handle Pinch Theorem (`compile_plat_surface`)

Let $L = \widehat{\beta}^{\,\text{plat}} \subset \mathbb{R}^3$ be a classical knot or link presented in $b$-bridge plat position by a braid word $\beta \in B_{2b}$ of length $|\beta|$ on vertical strands $0, 1, \dots, 2b - 1$, closed at the top and bottom by the $b$ bridge caps $(2k, 2k+1)$ ($0 \le k < b$). Let $\pi_\beta \in S_{2b}$ map equatorial column $i$ back to its incident top-cap strand $\pi_\beta(i) \in \{0, \dots, 2b - 1\}$, and define the **cap-incidence multigraph** $\Gamma(\beta; S)$ on vertex set $V = \{0, \dots, b - 1\}$ (the $b$ unlinked 2-spheres of $D(\beta; \emptyset)$) with one edge $e_k$ joining $u_k = \lfloor \pi_\beta(2k)/2 \rfloor$ and $v_k = \lfloor \pi_\beta(2k+1)/2 \rfloor$ for each neck $k \in S$. Because the pairs $\{2k, 2k+1\}$ partition $\{0, \dots, 2b - 1\}$ and $\pi_\beta$ is a bijection, every vertex of $\Gamma(\beta; S)$ has degree at most $2$, so every connected component $\Gamma_i \subseteq \Gamma(\beta; S)$ is either a **path** ($b_1(\Gamma_i) = 0$) or a **cycle** ($b_1(\Gamma_i) = 1$, including self-loops).

Given $S \subseteq \{0, 1, \dots, b - 1\}$ and an integer twist assignment $m_k \in \mathbb{Z}$ ($k \in S$), `compile_plat_surface` constructs the marked diagram $D(\beta; \{m_k\}_{k \in S})$ by placing $\beta$ in the upper half, $\beta^{-1}$ in the lower half, closing $(2k, 2k+1)$ by standard top and bottom caps, and inserting across columns $(2k, 2k+1)$ at the equator ($k \in S$) a **horizontal `Z-Y` 1-handle gadget** with $|m_k|$ classical crossings (`X0` if $m_k > 0$, `X1` if $m_k < 0$) wired in series along the two horizontal rails between a left `Z`-marker on strand $2k$ and a right `Y`-marker on strand $2k+1$.

**Theorem 1.1 (Planarity, Admissibility, Cap-Cycle Orientability Parity, and Constructive $\operatorname{ch}$-Bound).**
*For every $b \ge 1$, every braid word $\beta \in B_{2b}$, and every neck-twist assignment $\{m_k\}_{k \in S}$:*
1. **Strict Planarity & Complexity**: *The 4-regular graph $D = D(\beta; \{m_k\}_{k \in S})$ is strictly planar (`_is_planar_4reg_diagram(D) == True`) with exact vertex complexity*
   $$\operatorname{ch}(D) \;=\; 2|\beta| \;+\; \sum_{k \in S} \bigl(2 + |m_k|\bigr).$$
2. **Unconditional Admissibility**: *Both resolutions $L_+(D)$ and $L_-(D)$ reduce constructively via $|m_k|$ Reidemeister I untwists followed by $|\beta|$ free Reidemeister II cancellations to the $b$-component classical unlink.*
3. **Vanishing Normal Euler Number & Cycle-Twist Orientability Parity**: *Every component $F_i \subset F(D)$ has vanishing normal Euler number $e(F_i) = w(L_-(D)|_{F_i}) - w(L_+(D)|_{F_i}) = 0$, and components $F_i$ correspond bijectively to connected components $\Gamma_i$ of the cap-incidence multigraph $\Gamma(\beta; S)$:*
   - *If $\Gamma_i$ is a **path** ($b_1(\Gamma_i) = 0$, e.g. Artin spins and twisted symmetric unions $U_m(K)$ where $|S| = b - 1$ forms a spanning tree on the $b$ spheres), then $F_i \cong S^2$ ($\chi = 2$) is **orientable for every integer twist assignment $\{m_k\}$** (even or odd).*
   - *If $\Gamma_i = C$ is a **cycle** ($b_1(C) = 1$, $\chi(F_i) = 0$), then $F_i$ is an **orientable torus $T^2$** if and only if $\sum_{k \in C} m_k \equiv 0 \pmod{2}$, and a **non-orientable Klein bottle $K^2$** if and only if $\sum_{k \in C} m_k \equiv 1 \pmod{2}$.*

*Proof.* Planarity holds because each braid generator $\sigma_i^{\pm 1}$ acts on adjacent columns $(i, i+1)$ with cyclic CCW port order $(\text{NE} = i+1_{\text{in}}, \text{NW} = i_{\text{in}}, \text{SW} = i_{\text{out}}, \text{SE} = i+1_{\text{out}})$, each neck gadget lies strictly inside the equatorial corridor between columns $2k$ and $2k+1$, and the $b$ top/bottom caps nest disjointly in the outer regions.

At the $k$-th equatorial gadget (`vz = Z` with CCW ports `(0=top_rail, 1=sL_in, 2=sL_out, 3=bot_rail)` and `vy = Y` with CCW ports `(0=sR_in, 1=top_rail, 2=bot_rail, 3=sR_out)`):
- In $L_+(D)$, `Z` pairs `(0,1), (2,3)` and `Y` pairs `(0,3), (1,2)`. Thus strand $2k+1$ passes straight vertically (`vy.0 -- vy.3`), while strand $2k$ enters `vz.1 -- vz.0`, traverses the $|m_k|$ horizontal crossings, reflects around the end cap `vy.1 -- vy.2`, traverses the $|m_k|$ crossings back, and exits `vz.3 -- vz.2`. This forms a cascade of $|m_k|$ innermost Reidemeister I kinks on strand $2k$, collapsing to two parallel vertical strands $(2k, 2k+1)$.
- Symmetrically, in $L_-(D)$, `Z` pairs `(1,2), (0,3)` (vertical pass-through on strand $2k$, end cap `vz.0 -- vz.3`) and `Y` pairs `(0,1), (2,3)`, collapsing the $|m_k|$ kinks onto strand $2k+1$.

Once every equatorial neck collapses to parallel vertical strands, $\beta$ sits directly atop $\beta^{-1}$, canceling letter-by-letter via $|\beta|$ Reidemeister II moves to the $b$ unnested bridge circles ($c_+(D) = c_-(D) = b$). Orienting each cap circle $u \in \{0, \dots, b - 1\}$ so top strand $2u$ runs upward and $2u + 1$ downward, neck $k \in S$ attaches a 1-handle between spheres $u_k$ and $v_k$ with $\mathbb{Z}_2$ orientation-reversal parity $s_k \equiv \pi_\beta(2k) + \pi_\beta(2k+1) + 1 + m_k \pmod{2}$. On any path component $\Gamma_i$, $H^1(\Gamma_i; \mathbb{Z}_2) = 0$, so flipping sphere orientations along the tree trivializes every $s_k$ and yields an orientable $S^2$ ($\chi = 2|V_i| - 2|E_i| = 2$). Around any cycle $C$ of length $\ell$, each vertex $u \in V(C)$ contributes both cap strands $\{2u, 2u+1\}$ ($2u + (2u+1) \equiv 1 \pmod{2}$), whence $\sum_{k \in C} (\pi_\beta(2k) + \pi_\beta(2k+1) + 1) \equiv \ell + \ell \equiv 0 \pmod{2}$ and the monodromy around $C$ is $\sum_{k \in C} s_k \equiv \sum_{k \in C} m_k \pmod{2}$, giving $T^2$ when even and $K^2$ when odd. Finally, every crossing in $\beta$ cancels the writhe of its mirror partner in $\beta^{-1}$ within each of $L_-(D)$ and $L_+(D)$, while the $|m_k|$ twist crossings in neck $k$ resolve to $|m_k|$ self-crossings of a single circle with **equal total writhe** $m_k$ in both $L_-(D)$ and $L_+(D)$ ($w(L_-(D)|_{F_i}) = w(L_+(D)|_{F_i}) = \sum_{k \in E_i} m_k$, generally non-zero when $\sum m_k \neq 0$); subtracting gives $e(F_i) = w(L_-(D)|_{F_i}) - w(L_+(D)|_{F_i}) = 0$. $\blacksquare$

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
3. *By Fox free differential calculus on $r = w^{-1} a w b^{-1}$ (with fixed oriented meridian $a, b \mapsto t$), the single-variable Alexander polynomial of $R(w)$ is given up to units $\pm t^{\pm k} \in \mathbb{Z}[t^{\pm 1}]$ by*
   $$\Delta_{R(w)}(t) \;\doteq\; 1 \;+\; (t - 1)\left.\frac{\partial w}{\partial a}\right|_{a = b = t}$$
   *(while unoriented $H_1$-generator canonicalization under $\operatorname{Aut}(\mathbb{Z}) = \{t \mapsto t,\; t \mapsto t^{-1}\}$ in `e1_generators` selects the lexicographically smaller of $\Delta_{R(w)}(t)$ and its reciprocal $\pm t^{\deg \Delta}\Delta_{R(w)}(t^{-1})$ when $\Delta_{R(w)}$ is non-palindromic).*

*Proof.* In broken-surface movie coordinates, a 2-rail ribbon band pierces the disk $D_g$ with sign $\varepsilon \in \{\pm 1\}$ whenever both rails cross one boundary arc of $O_g$ **over** ($\varepsilon = +1$: `X1` on horizontal entry) or **under** ($\varepsilon = -1$: `X0` on horizontal entry) and cross the opposite boundary arc of $O_g$ with the opposite over/under sense (`X0` vs `X1`). Crucially, the 4D ribbon intersection sign depends *only* on the over/under sequence of the two boundary crossings — **not** on whether the band approaches $O_g$ from `NORTH` (top-to-bottom) or `SOUTH` (bottom-to-top). Alternating `NORTH` and `SOUTH` after each letter places every new crossing pair on the outermost boundary of $O_a$ (West) or $O_b$ (East), keeping both `cur_west, cur_east` and the terminal `Y`-attachment exposed on the infinite exterior face at every step.

Differentiating $r = w^{-1} a w b^{-1}$ in $\mathbb{Z}[F(a, b)]$ via $\partial(uv)/\partial a = \partial u/\partial a + u\,\partial v/\partial a$ and $\partial(w^{-1})/\partial a = -w^{-1}\,\partial w/\partial a$ yields
$$\frac{\partial r}{\partial a} \;=\; -w^{-1}\frac{\partial w}{\partial a}\,a w b^{-1} + w^{-1} + w^{-1} a \frac{\partial w}{\partial a} b^{-1}.$$
Under abelianization $\mathfrak{a}\colon \mathbb{Z}[F(a, b)] \to \mathbb{Z}[t^{\pm 1}]$ ($a, b \mapsto t$, so $\mathfrak{a}(w) = t^{\sigma(w)}$ cancels $\mathfrak{a}(w^{-1}) = t^{-\sigma(w)}$), this simplifies to
$$\mathfrak{a}\!\left(\frac{\partial r}{\partial a}\right) \;=\; t^{-1}\!\left(1 + (t - 1)\left.\frac{\partial w}{\partial a}\right|_{a=b=t}\right), \qquad \mathfrak{a}\!\left(\frac{\partial r}{\partial b}\right) \;=\; -\mathfrak{a}\!\left(\frac{\partial r}{\partial a}\right)$$
(since $(a - 1)\partial_a r + (b - 1)\partial_b r = r - 1 \mapsto 0$), so $E_1(R(w))$ is principal and generated by $1 + (t - 1)\partial_a w|_{a=b=t}$. $\blacksquare$

---

## 2. Hub Cross-Join Identifications ($\operatorname{ch} \le 12$)

Cross-joining all **`64`** compiled family instances against the exhaustive $\operatorname{ch} \le 11$ Hub (`23` Yoshikawa $\operatorname{ch} \le 10$ classes + `10` $\operatorname{ch} = 11$ classes) and $\operatorname{ch} = 12$ Hub (`34` $\operatorname{ch} = 12$ classes in [`expanded_ch12.json`](../code/data/expanded_ch12.json)) identifies **`22` distinct Hub classes across `27` family parameterizations** — strictly distinguishing **`MOVE_CERTIFIED`** identifications (`17` Hub classes connected by explicit Yoshikawa move paths + $\mathrm{O}(4)$ symmetries, proving ambient isotopy) from **`INVARIANT_MATCH`** identifications (`5` Hub classes matching across all 10 filtration layers where `compiled_ch` remains above a move-search barrier):

| Hub Label | $\operatorname{ch}$ | Topology | Realizing Spoke Instance(s) (`compiled_ch`) | Match Mode | Geometric Interpretation |
| :--- | :---: | :---: | :--- | :---: | :--- |
| $0_1$ | $0$ | $S^2$ | `Spin(0_1)` ($0$), `R(1)` ($0$) | `MOVE_CERTIFIED` | Unknotted 2-sphere |
| $2_1^1$ | $2$ | $T^2$ | `Tube(0_1)` ($2$) | `MOVE_CERTIFIED` | Unknotted Satoh torus $\operatorname{Tube}(\text{unknot})$ |
| $6_1^{0,1}$ | $6$ | $S^2 \sqcup T^2$ | `Mixed(L2a1)` ($6$) | `MOVE_CERTIFIED` | Mixed sphere–torus Hopf link |
| $7_1^{0,-2}$ | $7$ | $S^2 \sqcup K^2$ | `KleinMixed(L2a1,1)` ($7$) | `MOVE_CERTIFIED` | Half-twisted Klein–sphere Hopf link |
| $8_1$ | $8$ | $S^2$ | `Spin(3_1)` ($8$), `R(b*a)` ($8$) | `MOVE_CERTIFIED` | Artin spun trefoil $\tau^0(3_1) \cong R(ba)$ |
| $8_1^{1,1}$ | $8$ | $T^2 \sqcup T^2$ | `Tube(L2a1)` ($8$) | `MOVE_CERTIFIED` | Satoh tubed Hopf torus-link $\operatorname{Tube}(L2a1)$ |
| $9_1$ | $9$ | $S^2$ | `SymUnion(3_1,1)` ($9$), `R(b*a^-1)` ($9$) | `MOVE_CERTIFIED` | 1-twisted symmetric union $U_1(3_1) \cong R(ba^{-1})$ |
| $9_1^{1,-2}$ | $9$ | $T^2 \sqcup K^2$ | `KleinMixed(L2a1,10)` ($9$) | `MOVE_CERTIFIED` | Torus–Klein-bottle Hopf link |
| $10_1$ | $10$ | $S^2$ | `Spin(4_1)` ($10$), `R(b^-1*a*b*a^-1)` ($17$) | `INVARIANT_MATCH` | Artin spun figure-eight $\tau^0(4_1) \cong R([b, a^{-1}])$ (`EXACT_SINGLETON` Hub cell) |
| $10_1^1$ | $10$ | $T^2$ | `Tube(3_1)` ($10$) | `INVARIANT_MATCH` | Satoh spun-torus $\operatorname{Tube}(3_1)$ (`EXACT_SINGLETON` Hub cell) |
| $10_1^{0,1}$ | $10$ | $S^2 \sqcup T^2$ | `Mixed(L4a1)` ($10$) | `MOVE_CERTIFIED` | Mixed sphere–torus Solomon's link |
| $10_1^{-2,-2}$ | $10$ | $K^2 \sqcup K^2$ | `KleinTube(L2a1,11)` ($10$) | `MOVE_CERTIFIED` | Double Klein-bottle Hopf link |
| **$11_1$** | **$11$** | $S^2$ | `SymUnion(4_1,1)` (`MOVE`, $11$), `R(b^-1*a*b^-1*a^-1)` (`INV`, $17$) | **`MOVE_CERTIFIED`** | 1-twisted symmetric union $U_1(4_1)$ of the figure-eight knot |
| **$11_2$** | **$11$** | $S^2$ | `R(b^-1*a*b*a)` ($17$) | `INVARIANT_MATCH` | Asymmetric 4-syllable 1-fusion ribbon 2-knot $R(b^{-1}aba)$ (`EXACT_SINGLETON`) |
| **$11_3^{0,-2}$** | **$11$** | $S^2 \sqcup K^2$ | `KleinMixed(L4a1,1)` ($11$) | **`MOVE_CERTIFIED`** | Half-twisted Klein–sphere Solomon surface-link |
| **$12_1$** | **$12$** | $S^2$ | `Spin(5_2)` ($12$) | **`MOVE_CERTIFIED`** | Artin spun three-twist knot $\tau^0(5_2)$ (`NEW_PACKET`, `3` comps) |
| **$12_7$** | **$12$** | $S^2$ | `R(b^-1*a*b*a^-1*b*a)` ($25 \to 12$) | `INVARIANT_MATCH` | 6-syllable 1-fusion ribbon twin of `CS_A1` (`NEW_PACKET`, `4` comps) |
| **$12_9$** | **$12$** | $S^2$ | `R(b^2*a)` ($14 \to 12$) | `INVARIANT_MATCH` | 3-syllable 1-fusion ribbon 2-knot with $\mathbb{Z}HS^4$ $\Sigma_2$ (`EXACT_SINGLETON`) |
| **$12_{11}$** | **$12$** | $S^2$ | `Spin(5_1)` ($12$) | **`MOVE_CERTIFIED`** | Artin spun $(2,5)$-torus knot $\tau^0(5_1)$ (`EXACT_SINGLETON`) |
| **$12_2^1$** | **$12$** | $T^2$ | `Tube(4_1)` ($12$) | **`MOVE_CERTIFIED`** | Satoh ribbon torus of the figure-eight knot (`NEW_PACKET`, `4` comps) |
| **$12_4^{0,1}$** | **$12$** | $S^2 \sqcup T^2$ | `Mixed(L5a1)` ($12$) | **`MOVE_CERTIFIED`** | Mixed sphere–torus Whitehead link (`NEW_PACKET`, `4` comps) |
| **$12_3^{1,1}$** | **$12$** | $T^2 \sqcup T^2$ | `Tube(L4a1)` ($12$) | **`MOVE_CERTIFIED`** | Satoh 2-torus Solomon link (`NEW_PACKET`, `4` comps) |

### Residual Non-Elementary Census Entries ($\operatorname{ch} \le 12$)
Subtracting the $22$ Spoke-realized Hub entries highlights the genuinely non-ribbon, twist-spun, non-zero-Euler, and non-plat prime entries up to $\operatorname{ch} \le 12$:
- **Non-zero normal Euler number ($e(F_i) \neq 0$, obstructed from all `compile_plat_surface` outputs by Theorem 1.1(3))**: $2_1^{-1}$, $8_1^{-1,-1}$, $10_1^{-1,-1}$, and at $\operatorname{ch} = 12$: $12_1^{-1,-1}, 12_2^{-1,-1}, 12_1^{1,-2}, 12_1^{-2,-2}$.
- **Positive Zeeman twist-spins & non-ribbon $2$-knots**: $10_2 \cong \tau^2(3_1)$, $10_3 \cong \tau^3(3_1)$, $10_1^{1,1}$ ($\operatorname{Dlk} \neq 0$); and at $\operatorname{ch} = 12$: the non-solvable $\mathbb{F}_9$ $2$-knot **$12_2$** ($E_1 = (3, x^2+1)$), the 2-twist-spins **$12_4 \simeq \tau^2(4_1)$**, **$12_5 \simeq \tau^2(5_1)$**, **$12_6 \simeq \tau^2(5_2)$**, and the two **Cappell–Shaneson $2$-knots $12_8$ (`CS_A1`) and $12_{10}$ (`CS_A0`)**.
- **Higher-fusion / non-plat tori, Klein bottles, & multi-component links**: $9_1^{0,1}$, $10_2^{0,1}$, $10_1^{0,0,1}$, $10_1^{0,-2}$, $10_2^{0,-2}$; at $\operatorname{ch} = 11$: **$11_1^{-2}$** ($\pi_1 \cong 2O$), $11_1^{0,1}, 11_2^{0,1}, 11_1^{0,-2}, 11_2^{0,-2}, 11_1^{1,-2}, 11_1^{0,0,-2}$; and $16$ further prime $\operatorname{ch} = 12$ classes (including `12_3` ($S^2$, $\Delta = 3x^2 - 3x + 1$), turned torus `12_1^1`, and $14$ zero-Euler multi-component links).

---

## 3. Certified Exact $\operatorname{ch}(F) = 12$ & $\operatorname{ch}(F) = 13$ Determinations

Combining Theorem 1.1's constructive upper bounds with the $\operatorname{ch} \le 11$ and $\operatorname{ch} = 12$ censuses yields **five move-certified $\operatorname{ch}(F) = 12$ theorems** (plus two 10-layer invariant-matched $\operatorname{ch} = 12$ ribbon words) and **four exact $\operatorname{ch}(F) = 13$ determinations** (all lower bounds conditional on completeness of the lower-$\operatorname{ch}$ censuses):

**Theorem 3.1 (Exact $\operatorname{ch}(F) = 12$ Surface-Links).**
*The following five weakly prime surface-knots and surface-links have constructive `compiled_ch = 12` diagrams (`MOVE_CERTIFIED` against [`expanded_ch12.json`](../code/data/expanded_ch12.json)) and **exact** Yoshikawa ch-index $\operatorname{ch}(F) = 12$ (conditional on $\operatorname{ch} \le 11$ census completeness):*
1. **$\tau^0(5_1)$ (`Spin(5_1)` = `12_11`)**: Artin spun $(2,5)$-torus knot ($S^2$, $\Delta(t) = t^4 - t^3 + t^2 - t + 1$, $H_1(\Sigma_2) = \mathbb{Z}_5$).
2. **$\tau^0(5_2)$ (`Spin(5_2)` = `12_1`)**: Artin spun three-twist knot ($S^2$, $\Delta(t) = 2t^2 - 3t + 2$, $H_1(\Sigma_2) = \mathbb{Z}_7$).
3. **$\operatorname{Tube}(4_1)$ (`Tube(4_1)` = `12_2^1`)**: Satoh ribbon torus of the figure-eight knot ($T^2$, $\Delta(t) = t^2 - 3t + 1$, $H_1(\Sigma_2) = \mathbb{Z}_5, b_2(\Sigma_2) = 2$).
4. **$\operatorname{Tube}(L4a1)$ (`Tube(L4a1)` = `12_3^{1,1}` )**: Satoh 2-torus link of Solomon's $(2,4)$-torus link ($T^2 \sqcup T^2$, $E_1 = (x^2 - y^2,\; xy - x + y^2 - y)$, $H_1(\Sigma_2) = \mathbb{Z}_4, b_2(\Sigma_2) = 2$).
5. **$\operatorname{Mixed}(L5a1)$ (`Mixed(L5a1)` = `12_4^{0,1}` )**: Mixed $S^2 \sqcup T^2$ Whitehead link ($E_1 = (x^2y - x^2 - 2xy + 2x + y - 1,\; xy^2 - 2xy + x - y^2 + 2y - 1)$, $H_1(\Sigma_2) = \mathbb{Z}_8$).

*(In addition, `R(b^2*a)` and `R(b^-1*a*b*a^-1*b*a)` match `12_9` and `12_7` across all 10 filtration layers via `INVARIANT_MATCH`, tightening their constructive bounds $[12, 14]$ and $[12, 25]$ to $\operatorname{ch} = 12$ conditional on invariant faithfulness.)*

*Proof.* For each of the five surfaces $F$, `compile_plat_surface` constructs an explicit planar, Reidemeister-unlink-certified marked diagram $D$ with $\operatorname{ch}(D) = 12$ (`Spin(5_1)` and `Spin(5_2)` have $2 \cdot 5 + 2 = 12$; `Tube(4_1)` and `Tube(L4a1)` have $2 \cdot 4 + 4 = 12$; `Mixed(L5a1)` has $2 \cdot 5 + 2 = 12$), establishing $\operatorname{ch}(F) \le 12$ unconditionally.

Conversely, suppose $\operatorname{ch}(F) \le 11$, so $F$ admits some minimal diagram $D_{\min}$ with $\operatorname{ch}(D_{\min}) \le 11$.
- $D_{\min}$ cannot be disconnected/split: for $c = 1$ (`Spin(5_1)`, `Spin(5_2)`, `Tube(4_1)`) $F$ is connected; for `Tube(L4a1)` and `Mixed(L5a1)`, $E_1(F) \neq (0)$, whereas every split 2-component orientable link satisfies $E_1 = (0)$.
- $D_{\min}$ cannot be a nontrivial connected sum $F_1 \# F_2$ with $\operatorname{ch}(F_1) + \operatorname{ch}(F_2) \le 12$: our reference catalog `build_reference_catalog(max_composite_ch=12)` enumerates every connected sum of $\operatorname{ch} \le 10$ prime classes and standard surfaces up to total $\operatorname{ch} \le 12$, and none matches the 10-layer invariant key $I(F)$.
  > **Subtlety on $\operatorname{Tube}(4_1)$ vs $10_1 \# 2_1^1$**: Both $\operatorname{Tube}(4_1)$ and the $\operatorname{ch} = 12$ composite torus $\tau^0(4_1) \# T^2_{\text{unknotted}} = 10_1 \# 2_1^1$ have homeomorphism type $T^2$, identical group $\pi_1 \cong \pi_1(S^3 \setminus 4_1)$, identical $E_1 = (x^2 - 3x + 1)$, identical $\Sigma_2$ invariants $(b_1 = 0, H_1 = \mathbb{Z}_5, b_2 = 2, \sigma = 0)$, and identical kei, bikei, and Alexander biquandle invariants. However, Layer 6 regular $A_4$-cover Betti profiles and Layer 8 peripheral $\operatorname{SL}(2, \mathbb{F}_3)$ holonomy separate them (`test_tube_4_1_separated_from_composite_10_1_conn_2_1_1` in [`test_families.py`](../code/tests/test_families.py)):
  > - **$A_4$-cover Betti profile $(b_1, \operatorname{Tors} H_1, b_2, b_3)$**: $(4, \mathbb{Z}_2^2, 30, 3)$ on $\operatorname{Tube}(4_1)$ vs $(4, \mathbb{Z}_2^2, 33, 6)$ on $10_1 \# 2_1^1$.
  > - **Non-meridional peripheral $\operatorname{SL}(2, \mathbb{F}_3)$ representations**: $24$ on $\operatorname{Tube}(4_1)$ vs $0$ on $10_1 \# 2_1^1$.
- Machine comparison confirms $I(F)$ is disjoint from all `23` Yoshikawa $\operatorname{ch} \le 10$ classes, all `10` $\operatorname{ch} = 11$ prime classes, and all `93` $\operatorname{ch} = 11$ `AMBIGUOUS_SHADOW` residuals, proving $\operatorname{ch}(F) \ge 12$ and hence $\operatorname{ch}(F) = 12$. $\blacksquare$

**Corollary 3.2 (Four Certified Exact $\operatorname{ch}(F) = 13$ Surface-Links).**
*The four Spoke families **`SymUnion(5_1, 1)`** ($U_1(5_1)$, $S^2$, $\Delta \doteq 2x - 3$), **`SymUnion(5_2, 1)`** ($U_1(5_2)$, $S^2$, $\Delta \doteq 2x^2 - 4x + 1$), **`R(b^2*a^-1)`** ($S^2$, $\Delta \doteq x^2 + x - 1$), and **`KleinMixed(L5a1, 1)`** ($S^2 \sqcup K^2$, $H_1(\Sigma_2) = \mathbb{Z}_8$) each admit an explicit planar, Reidemeister-certified diagram of complexity `compiled_ch = 13` (`families.py`), have irreducible / non-composite 10-layer invariants disjoint from all $\operatorname{ch} \le 13$ composites, and have **zero** invariant matches across the entire $\operatorname{ch} \le 12$ census (`expanded_ch12.json` entries and `742` residuals). Conditional on completeness of the $\operatorname{ch} \le 12$ census, all four satisfy **exact** $\operatorname{ch}(F) = 13$.*

---

## 4. Complete Catalog Summary of Spokes 1–3 (`family_catalog.json` + $\operatorname{ch} = 12$ Hub Cross-Join)

*(In the tables below, $[\operatorname{ch}_{\min}, \operatorname{ch}_{\max}]$ reflects the sharpened post-$\operatorname{ch}=12$ bounds, while the "Separating Witness" column matches the `family_catalog.json` strong-Gröbner witness assignment.)*

### 4.1 Spoke 1: Artin Spun 2-Knots $\tau^0(K)$ & Symmetric-Union 2-Knots $U_1(K)$ (`22` instances)

| `family_id` | Classical Base | $[\operatorname{ch}_{\min}, \operatorname{ch}_{\max}]$ | Exact? | Hub Match | $E_1(F)$ Generator (`family_catalog.json`) | $\Sigma_2$: $\operatorname{Tors} H_1$ | Nearest Neighbor (Separating Layer) |
| :--- | :--- | :---: | :---: | :---: | :--- | :---: | :--- |
| `Spin(0_1)` | $0_1$ ($0/1$) | $[0, 0]$ | **Yes** | $0_1$ (`MOVE`) | $1$ | $0$ | $8_1$ (`E1`) |
| `Spin(3_1)` | $3_1$ ($3/1$) | $[8, 8]$ | **Yes** | $8_1$ (`MOVE`) | $x^2 - x + 1$ | $\mathbb{Z}_3$ | $9_1$ (`E1`) |
| `Spin(4_1)` | $4_1$ ($5/3$) | $[10, 10]$ | **Yes** | $10_1$ (`INV`) | $x^2 - 3x + 1$ | $\mathbb{Z}_5$ | $10_2$ (`E1`) |
| **`Spin(5_1)`** | $5_1$ ($5/1$) | **$[12, 12]$** | **Yes** | **$12_{11}$ (`MOVE`)** | $x^4 - x^3 + x^2 - x + 1$ | $\mathbb{Z}_5$ | `Spin(5_2)` (`E1`) |
| **`Spin(5_2)`** | $5_2$ ($7/4$) | **$[12, 12]$** | **Yes** | **$12_1$ (`MOVE`)** | $2x^2 - 3x + 2$ | $\mathbb{Z}_7$ | `Spin(5_1)` (`E1`) |
| `Spin(6_1)` | $6_1$ ($9/5$) | $[13, 14]$ | No | — | $2x^2 - 5x + 2$ | $\mathbb{Z}_9$ | `R(b^2*a)` (`E1`) |
| `Spin(6_2)` | $6_2$ ($11/7$) | $[13, 14]$ | No | — | $x^4 - 3x^3 + 3x^2 - 3x + 1$ | $\mathbb{Z}_{11}$ | `R(b^2*a)` (`E1`) |
| `Spin(6_3)` | $6_3$ ($13/8$) | $[13, 14]$ | No | — | $x^4 - 3x^3 + 5x^2 - 3x + 1$ | $\mathbb{Z}_{13}$ | `R(b^2*a)` (`E1`) |
| `Spin(7_1)`–`Spin(7_7)` | $7_1 \dots 7_7$ | $[13, 16]$ | No | — | $\Delta_{7_k}(x)$ (all 7 distinct) | $\mathbb{Z}_7 \dots \mathbb{Z}_{21}$ | Pairwise separated at `E1` |
| `SymUnion(3_1,1)` | $3_1$, $m=1$ | $[9, 9]$ | **Yes** | $9_1$ (`MOVE`) | $x - 2$ | $\mathbb{Z}_3$ | $8_1$ (`E1`) |
| **`SymUnion(4_1,1)`** | $4_1$, $m=1$ | **$[11, 11]$** | **Yes** | **$11_1$ (`MOVE`)** | $2x^2 - 2x + 1$ | $\mathbb{Z}_5$ | $11_2$ (`E1`) |
| **`SymUnion(5_1,1)`** | $5_1$, $m=1$ | **$[13, 13]$** | **Yes** | — ($\operatorname{ch}=13$) | $2x - 3$ | $\mathbb{Z}_5$ | `R(b^2*a^-1)` (`E1`) |
| **`SymUnion(5_2,1)`** | $5_2$, $m=1$ | **$[13, 13]$** | **Yes** | — ($\operatorname{ch}=13$) | $2x^2 - 4x + 1$ | $\mathbb{Z}_7$ | `R(b^2*a^-1)` (`E1`) |
| `SymUnion(6_1,1)` | $6_1$, $m=1$ | $[13, 15]$ | No | — | $3x^2 - 4x + 2$ | $\mathbb{Z}_9$ | `SymUnion(6_2,1)` (`E1`) |
| `SymUnion(6_2,1)` | $6_2$, $m=1$ | $[13, 15]$ | No | — | $x^3 - 4x^2 + 4x - 2$ | $\mathbb{Z}_{11}$ | `SymUnion(6_1,1)` (`E1`) |
| `SymUnion(6_3,1)` | $6_3$, $m=1$ | $[13, 15]$ | No | — | $2x^3 - 6x^2 + 4x - 1$ | $\mathbb{Z}_{13}$ | `SymUnion(6_1,1)` (`E1`) |

### 4.2 Spoke 2: 1-Fusion Ribbon 2-Knots $R(w)$ (`17` instances)

| `family_id` | Word $w \in F(a, b)$ | $[\operatorname{ch}_{\min}, \operatorname{ch}_{\max}]$ | `compiled_ch` | Hub Match | Fox $\Delta_{R(w)}(t)$ (and `E1` if $t \leftrightarrow t^{-1}$ reversed) | $\Sigma_2$: $\operatorname{Tors} H_1$ | Nearest Neighbor (Separating Layer) |
| :--- | :--- | :---: | :---: | :---: | :--- | :---: | :--- |
| `R(1)` | $1$ | $[0, 0]$ | $0$ | $0_1$ (`MOVE`) | $1$ | $0$ | $8_1$ (`E1`) |
| `R(b*a)` | $ba$ | $[8, 8]$ | $8$ | $8_1$ (`MOVE`) | $t^2 - t + 1$ | $\mathbb{Z}_3$ | $9_1$ (`E1`) |
| `R(b*a^-1)` | $ba^{-1}$ | $[9, 9]$ | $9$ | **$9_1$ (`MOVE`)** | $t - 2$ (`2x - 1`) | $\mathbb{Z}_3$ | $8_1$ (`E1`) |
| `R(b^-1*a*b*a^-1)` | $b^{-1}aba^{-1}$ | $[10, 10]$ | $17$ | $10_1$ (`INV`) | $t^2 - 3t + 1$ | $\mathbb{Z}_5$ | $10_2$ (`E1`) |
| `R(b^-1*a*b^-1*a^-1)` | $b^{-1}ab^{-1}a^{-1}$ | $[11, 11]$ | $17$ | **$11_1$ (`INV`)** | $2t^2 - 2t + 1$ (`x^2 - 2x + 2`) | $\mathbb{Z}_5$ | $11_2$ (`E1`) |
| `R(b^-1*a*b*a)` | $b^{-1}aba$ | $[11, 11]$ | $17$ | **$11_2$ (`INV`)** | $t^3 - t^2 + 2t - 1$ (`x^3 - 2x^2 + x - 1`) | $\mathbb{Z}_5$ | $11_1$ (`E1`) |
| **`R(b^2*a)`** | $b^2 a$ | **$[12, 12]$** | $14$ | **$12_9$ (`INV`)** | $t^3 - t^2 + 1$ (`x^3 - x + 1`) | $0$ ($\mathbb{Z}HS^4$) | `Spin(6_1)` (`E1`) |
| **`R(b^2*a^-1)`** | $b^2 a^{-1}$ | **$[13, 13]$** | $13$ | — ($\operatorname{ch}=13$) | $t^2 - t - 1$ (`x^2 + x - 1`) | $0$ ($\mathbb{Z}HS^4$) | `SymUnion(5_1,1)` (`E1`) |
| `R(b^3*a)` | $b^3 a$ | $[13, 18]$ | $18$ | — | $t^4 - t^3 + 1$ (`x^4 - x + 1`) | $\mathbb{Z}_3$ | `R(b^2*a^-2)` (`E1`) |
| `R(b^3*a^-1)` | $b^3 a^{-1}$ | $[13, 17]$ | $17$ | — | $t^3 - t^2 - 1$ (`x^3 + x - 1`) | $\mathbb{Z}_3$ | `Spin(7_1)` (`E1`) |
| `R(b^2*a^2)` | $b^2 a^2$ | $[13, 18]$ | $18$ | — | $t^4 - t^2 + 1$ | $0$ ($\mathbb{Z}HS^4$) | `R(b^2*a^-2)` (`E1`) |
| `R(b^2*a^-2)` | $b^2 a^{-2}$ | $[13, 18]$ | $18$ | — | $t^2 - 2$ (`2x^2 - 1`) | $0$ ($\mathbb{Z}HS^4$) | `R(b^2*a^2)` (`E1`) |
| `R(b^-1*a^2*b*a^-1)` | $b^{-1}a^2 ba^{-1}$ | $[13, 21]$ | $21$ | — | $t^3 + t^2 - 2t - 1$ (`x^3 - x^2 - 2x + 1`) | $0$ ($\mathbb{Z}HS^4$) | `R(b^-1*a^2*b^-1*a^-1)` (`E1`) |
| `R(b^2*a*b^-1*a)` | $b^2 ab^{-1}a$ | $[13, 22]$ | $22$ | — | $2t^3 - 2t^2 + 1$ (`2x^3 - 2x^2 + 1`, recip. $t^3 - 2t + 2$) | $\mathbb{Z}_3$ | `R(b^-2*a*b*a^-1)` (`E1`) |
| **`R(b^-1*a*b*a^-1*b*a)`** | $b^{-1}aba^{-1}ba$ | **$[12, 12]$** | $25$ | **$12_7$ (`INV`)** | $t^3 - 2t^2 + 3t - 1$ (`x^3 - 3x^2 + 2x - 1`) | $\mathbb{Z}_7$ | `R(b^-2*a*b*a^-1)` (`E1`) |
| `R(b^-2*a*b*a^-1)` | $b^{-2}aba^{-1}$ | $[13, 22]$ | $22$ | — | $2t - 1$ | $\mathbb{Z}_3$ | $9_1$ (`finite_quotients`) |
| `R(b^-1*a^2*b^-1*a^-1)` | $b^{-1}a^2b^{-1}a^{-1}$ | $[13, 21]$ | $21$ | — | $1$ (trivial Alexander module!) | $0$ ($\mathbb{Z}HS^4$) | $0_1$ (`finite_quotients`) |

> [!NOTE]
> **Integer homology 4-sphere covers & Alexander-trivial / Alexander-twin 1-fusion knots:** Six Spoke 2 words (`R(b^2*a)`, `R(b^2*a^-1)`, `R(b^2*a^2)`, `R(b^2*a^-2)`, `R(b^-1*a^2*b*a^-1)`, and `R(b^-1*a^2*b^-1*a^-1)`) satisfy $|\Delta_{R(w)}(-1)| = 1$, so their double branched covers $\Sigma_2(S^4, R(w))$ are integer homology 4-spheres ($H_1 = 0, b_2 = 0$). Strong integral Gröbner canonicalization of $E_1$ separates the first five directly at Layer 3 (`E1`), while `R(b^-1*a^2*b^-1*a^-1)` has **trivial** Alexander polynomial $\Delta(t) \doteq 1$ (a 1-fusion analogue of Kinoshita–Terasaka, separated from `0_1` at Layer 5 `finite_quotients`) and `R(b^-2*a*b*a^-1)` shares $E_1 = (2x - 1)$ with `9_1` (separated at Layer 5 `finite_quotients`).

### 4.3 Spoke 3: Satoh Tubes, Mixed Sphere–Torus Links, & Klein-Bottle Twisted Plats (`25` instances)

| `family_id` | Topology | $[\operatorname{ch}_{\min}, \operatorname{ch}_{\max}]$ | Exact? | Hub Match | $\Sigma_2$: $(\operatorname{Tors} H_1, b_2)$ | Periph. $\operatorname{SL}(2, \mathbb{F}_3)$ | Highlights & Separation Status |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| `Tube(0_1)` | $T^2$ | $[2, 2]$ | **Yes** | $2_1^1$ (`MOVE`) | $(0, 2)$ | $0$ | Unknotted Satoh torus |
| `Tube(3_1)` | $T^2$ | $[10, 10]$ | **Yes** | $10_1^1$ (`INV`) | $(\mathbb{Z}_3, 2)$ | $24$ | Separated from $8_1 \# 2_1^1$ by peripheral $\operatorname{SL}(2, \mathbb{F}_3)$ |
| **`Tube(4_1)`** | $T^2$ | **$[12, 12]$** | **Yes** | **$12_2^1$ (`MOVE`)** | $(\mathbb{Z}_5, 2)$ | $24$ | **Exact $\operatorname{ch} = 12$**; separated from $10_1 \# 2_1^1$ by $A_4$-cover $b_2=30$ vs $33$ & periph. $\operatorname{SL}(2, \mathbb{F}_3)$ |
| `Tube(5_1)` | $T^2$ | $[13, 14]$ | No | — | $(\mathbb{Z}_5, 2)$ | $0$ | Pairwise distinct from `Tube(5_2)` at `E1` *(collides across all 10 current layers with composite `Spin(5_1) # 2_1^1`; see warning below)* |
| `Tube(5_2)`, `Tube(6_1)` | $T^2$ | $[13, 14\text{–}16]$ | No | — | $(\mathbb{Z}_7, 2), (\mathbb{Z}_9, 2)$ | $0$ | Ribbon spun tori of twist knots *(likewise collide on current 10 layers with `Spin(K) # 2_1^1`)* |
| `Tube(L2a1)` | $T^2 \sqcup T^2$ | $[8, 8]$ | **Yes** | $8_1^{1,1}$ (`MOVE`) | $(\mathbb{Z}_2, 2)$ | $122$ | Satoh Hopf 2-torus link ($\operatorname{Dlk} = 0$) |
| **`Tube(L4a1)`** | $T^2 \sqcup T^2$ | **$[12, 12]$** | **Yes** | **$12_3^{1,1}$ (`MOVE`)** | $(\mathbb{Z}_4, 2)$ | $148$ | **Exact $\operatorname{ch} = 12$** Solomon 2-torus link |
| `Tube(L5a1)` | $T^2 \sqcup T^2$ | $[13, 14]$ | No | — | $(\mathbb{Z}_8, 2)$ | $96$ | Linking-number-zero Whitehead 2-torus link |
| `Tube(L6a1)`–`Tube(L6a3)` | $T^2 \sqcup T^2$ | $[13, 16]$ | No | — | $(\mathbb{Z}_{12}, 2), (\mathbb{Z}_{10}, 2), (\mathbb{Z}_6, 2)$ | $148, 138, 162$ | All three 6-crossing 2-bridge torus links pairwise separated at `E1` |
| **`Tube(L6a4)`** | $T^2 \sqcup T^2 \sqcup T^2$ | $[13, 22]$ | No | — | $(\mathbb{Z}_4^2, 2)$ | $2664$ | **Brunnian 3-torus Borromean surface-link** (`is_brunnian = True`) |
| `Mixed(L2a1)` | $S^2 \sqcup T^2$ | $[6, 6]$ | **Yes** | $6_1^{0,1}$ (`MOVE`) | $(\mathbb{Z}_2, 0)$ | $69$ | Mixed Hopf sphere–torus link |
| `Mixed(L4a1)` | $S^2 \sqcup T^2$ | $[10, 10]$ | **Yes** | $10_1^{0,1}$ (`MOVE`) | $(\mathbb{Z}_4, 0)$ | $86$ | Mixed Solomon sphere–torus link |
| **`Mixed(L5a1)`** | $S^2 \sqcup T^2$ | **$[12, 12]$** | **Yes** | **$12_4^{0,1}$ (`MOVE`)** | $(\mathbb{Z}_8, 0)$ | $48$ | **Exact $\operatorname{ch} = 12$** mixed Whitehead link |
| `Mixed(L6a1)`–`Mixed(L6a3)` | $S^2 \sqcup T^2$ | $[13, 14]$ | No | — | $(\mathbb{Z}_{12}, 0), (\mathbb{Z}_{10}, 0), (\mathbb{Z}_6, 0)$ | $86, 69, 93$ | Mixed 6-crossing sphere–torus links |
| **`Mixed(L6a4,01)`** | $S^2 \sqcup T^2 \sqcup T^2$ | $[13, 16]$ | No | — | $(\mathbb{Z}_4^2, 0)$ | $1776$ | **Brunnian $S^2 \sqcup T^2 \sqcup T^2$ Borromean link** reduced from `raw = 20` to `ch = 16` |
| `KleinMixed(L2a1,1)` | $S^2 \sqcup K^2$ | $[7, 7]$ | **Yes** | $7_1^{0,-2}$ (`MOVE`) | $(\mathbb{Z}_2, 0)$ | $0$ | Half-twisted Hopf Klein–sphere link |
| `KleinMixed(L2a1,10)` | $T^2 \sqcup K^2$ | $[9, 9]$ | **Yes** | $9_1^{1,-2}$ (`MOVE`) | $(\mathbb{Z}_2, 2)$ | $9$ | Hopf torus–Klein-bottle link |
| `KleinTube(L2a1,11)` | $K^2 \sqcup K^2$ | $[10, 10]$ | **Yes** | $10_1^{-2,-2}$ (`MOVE`) | $(\mathbb{Z}_2, 2)$ | $0$ | Double Klein-bottle Hopf link ($\operatorname{rk}_{\mathbb{Z}}\pi_2 = 16$) |
| **`KleinMixed(L4a1,1)`** | $S^2 \sqcup K^2$ | **$[11, 11]$** | **Yes** | **$11_3^{0,-2}$ (`MOVE`)** | $(\mathbb{Z}_4, 0)$ | $0$ | Geometrically identifies $\operatorname{ch} = 11$ class $11_3^{0,-2}$ |
| **`KleinMixed(L5a1,1)`** | $S^2 \sqcup K^2$ | **$[13, 13]$** | **Yes** | — ($\operatorname{ch}=13$) | $(\mathbb{Z}_8, 0)$ | $0$ | **Exact $\operatorname{ch} = 13$** half-twisted Whitehead Klein–sphere link |

> [!WARNING]
> **Current Invariant Blindness on `Tube(K)` vs `Spin(K) # 2_1^1` for $K \in \{5_1, 5_2, 6_1\}$:** While `Tube(3_1)` vs `8_1 # 2_1^1` and `Tube(4_1)` vs `10_1 # 2_1^1` are separated cleanly by regular $A_4$-cover Betti numbers and peripheral $\operatorname{SL}(2, \mathbb{F}_3)$ holonomy, the 2-bridge knots $K \in \{5_1, 5_2, 6_1\}$ have determinants $5, 7, 9$ and admit **zero** non-abelian representations into $A_4$ or $\operatorname{SL}(2, \mathbb{F}_3)$ (and $D_{10}, D_{14}, D_{18}$ are not currently in `_COVER_TARGET_GROUPS`). Consequently, `Tube(K)` and the composite torus `Spin(K) # 2_1^1` (`ch = 14, 14, 16`) currently **collide across all 10 layers** of `compute_full_invariant_record` — motivating dihedral $D_{2p}$ regular-cover Betti profiles and peripheral $\operatorname{SL}(2, \mathbb{F}_p)$ ($p \in \{5, 7\}$) holonomy for $\operatorname{ch} \ge 14$ torus censuses.
