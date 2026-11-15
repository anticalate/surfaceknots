# Expanded Census of Knotted Surfaces & Surface-Links at $\operatorname{ch}(F) = 11$

Extending Yoshikawa's (1994) $\operatorname{ch}(F) \le 10$ census by one vertex requires crossing a steep combinatorial jump: from $803$ prime 4-regular spherical shadows at $n = 10$ to **$3{,}378$** shadows at $n = 11$ (`OEIS A113201`, dual to simple planar quadrangulations with $13$ vertices), corresponding to $3{,}378 \times 4^{11} \approx 1.417 \times 10^{10}$ raw vertex assignments across all non-extreme marker counts $m \in \{1, \dots, 10\}$ ($m = 0$ and $m = 11$ being unknotted spheres and planar $\mathbb{R}^3$ surfaces, respectively).

Using the C99 + Python enumeration, Reidemeister unlink certification, multi-pass + post-flype mountain-pass stabilization pipeline, and Tier 1–3 4D invariant suite in `surfaceknots` (exported in [`expanded_ch11.json`](../code/data/expanded_ch11.json) and [`expanded_ch11.csv`](../code/data/expanded_ch11.csv)), we obtain an exhaustive survivor pool of **$105$** post-flype-tunnel `ch = 11` diagrams, which stratifies cleanly under our two-tier `Move-Component / Invariant Cell` schema (relying on Yoshikawa's $\operatorname{ch} \le 10$ classification, verified computationally by our blind C99 census through $\operatorname{ch} \le 10$, plus full Yoshikawa Definition 5.1 composite checks across all $\operatorname{ch}(A), \operatorname{ch}(B) \le 10$ connected sums) into:
- **10 certified new, weakly prime, non-splittable $\operatorname{ch} = 11$ invariant classes** ($12$ diagrams: **$9$ `EXACT_SINGLETON`** classes with $\operatorname{ch}(F) = 11$ exactness and **$1$ `NEW_PACKET`** $11_1.\{\mathrm{a},\mathrm{b},\mathrm{c}\}$ with $m = 3$ move-components whose fundamental groups are directly isomorphic — comprising three connected surface-knots $11_1, 11_2, 11_1^{-2}$, six 2-component surface-links $11_1^{0,1}, 11_2^{0,1}, 11_1^{0,-2}, 11_2^{0,-2}, 11_3^{0,-2}, 11_1^{1,-2}$, and one 3-component surface-link $11_1^{0,0,-2}$), all pairwise separated from one another and from every $\operatorname{ch} \le 10$, unknot, split, and $\operatorname{ch}(A), \operatorname{ch}(B) \le 10$ connected-sum class; and
- **93 documented `AMBIGUOUS_SHADOW` residual diagram collisions** (across $21$ collision cells) whose full 10-layer invariant profiles match known $\operatorname{ch} \le 10$ Yoshikawa entries ($19$), higher-genus standard unknots ($65$), connected sums ($8$), or split unions ($1$) lacking a bounded $\le 13$-vertex mountain-pass simplification path.

---

## 1. Exhaustive $\operatorname{ch} = 11$ Enumeration & Simplification Funnel

Every stage below preserves exact topological completeness (discarding a diagram **only** when an explicit planar 2-edge cut, disconnected projection, Kauffman bracket non-unlink obstruction on $L_\pm(D)$, or constructive Yoshikawa $\Omega_1$–$\Omega_8$ / 4D Tait flype / certified unknotted-tube simplification exhibits a split decomposition, nontrivial connected-sum decomposition, or $\operatorname{ch} \le 10$ representative; note that pre-tunnel counts are Kauffman-bracket-filtered candidates while **all** tunneled survivors carry explicit Reidemeister $\mathrm{R1}/\mathrm{R2}/\mathrm{R3}$ unlink certificates on $(L_-(D), L_+(D))$):

| Pipeline Stage | Count | Pruning / Certification Mechanism |
| :--- | ---: | :--- |
| Prime 4-regular spherical shadows ($n = 11$) | $3{,}378$ | Recursive quadrangulation generator (`OEIS A113201`) |
| Raw vertex assignments ($4^{11}$ per shadow) | $14{,}168{,}358{,}912$ | Marked bigon ($\Omega_6'$), linking-parity, & 2-stage $O(c \cdot 2^c)$ Fast Kauffman Bracket Butterfly Transform $\bmod (10^9 + 7)$ across $m = 1 \dots 10$ |
| Bracket + immediate 1-step filtered candidates (`m = 1..10`) | $78{,}378{,}425$ | ($66{,}017{,}754$ in $m = 2 \dots 7$ in stored `expanded_ch11.json`; $75{,}125{,}028$ in the live C99 pre-filter counter with identical survivor classes) |
| Level $\Omega_4/\Omega_5/\Omega_7/\Omega_8$ components explored | $26{,}164{,}808$ | ($21{,}483{,}971$ in $m = 2 \dots 7$; $26{,}036{,}915$ components immediately reduced) |
| Globally merged level-irreducible candidates (`m = 1..10`) | $39{,}400$ | ($39{,}179$ in $m = 2 \dots 7$; canonical spherical BFS deduplication modulo $\mathrm{O}(4)$ + Eulerian 2-edge-cut filter) |
| After progressive $+1/+2$ C99 mountain-pass cascade (`256/256`) | $205$ | 12–13-vertex $\Omega_1^+ / \Omega_6^+ / \Omega_2^+$ mountain-pass BFS + **constructive Reidemeister $\mathrm{R1}/\mathrm{R2}/\mathrm{R3}$ unlink certification** |
| After certified 4-edge cyclic-tube filter (`84` pruned) | $121$ | 4-edge planar annulus cuts bounding a certified unknotted $S^2$ punctured twice (`116` pre-stratification raw survivors in `ch11_raw_survivors.json` after enhanced tube pruning) |
| After 4D Tait flype canonicalization & post-flype C99 tunnel (`384/256`) | **$105$** | 4D Tait flypes unlock additional $\le 13$-vertex mountain-pass collapses (`12` certified new prime diagrams + `93` residual collisions) |
| **Certified new weakly prime $\operatorname{ch} = 11$ invariant classes** | **$10$** ($12$ diags) | **$9$ `EXACT_SINGLETON`** ($m = 1$) + **$1$ `NEW_PACKET`** ($11_1$, $m = 3$); zero overlap with $\operatorname{ch} \le 10$, unknot, split, and all $\operatorname{ch}(A), \operatorname{ch}(B) \le 10$ connected sums |

*(Note on `EXACT_SINGLETON` $\operatorname{ch}$-exactness: because prime-shadow enumeration excludes diagrams with a planar 2-edge cut, certifying $\operatorname{ch}(F) = N$ for a new singleton class requires both zero overlap with $\operatorname{ch} < N$ prime/unknot/split classes and zero overlap with all connected sums $A \# B$ of smaller-$\operatorname{ch}$ surfaces—supplied here by the exhaustive Yoshikawa Definition 5.1 composite check across $\operatorname{ch}(A), \operatorname{ch}(B) \le 10$, including exact proper-sublink profiles separating $11_2^{0,1}$ from $6_1^{0,1} \# 8_1$ and $11_1^{0,-2}$ from $7_1^{0,-2} \# 8_1$.)*

---

## 2. Certified New Weakly Prime Surface-Knots & Surface-Links ($\operatorname{ch}(F) = 11$)

Notation follows Yoshikawa's convention $11_k^{g_1, \dots, g_c}$ ($g_i \ge 0$ for orientable $\Sigma_{g_i}$, omitted for connected $S^2$-knots; $g_i < 0$ for non-orientable $\#_{|g_i|} \mathbb{R}P^2$, so $-2 \mapsto K^2 = \mathbb{R}P^2 \# \mathbb{R}P^2$).

### Part A: Connected Surface-Knots ($c = 1$; 3 certified new weakly prime classes)

| Label | Cell Status | Type | Ori. | $\chi$ | $e(F)$ | $(c, h)$ | $b(F)$ Bounds | Dehn-Simplified $\pi_1(\mathbb{R}^4 \setminus F)$ | Verified $E_1(F) \subset \mathbb{Z}[H_1]$ | Distinguishing Machine Invariants |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- | :--- | :--- |
| $11_1$ | `NEW_PACKET` ($m=3$) | $S^2$ | Yes | $2$ | $0$ | $(9, 2)$ | $4\text{–}9$ | $\langle x_1, x_2 \mid x_2^{-1} x_1^{-1} x_2^{-1} x_1 x_2 x_1 x_2^{-1} x_1^{-1} x_2 x_1 = 1 \rangle$ | $(2x^2 - 2x + 1)$ | `D10=30, F21=21, A5=180`; kei $\mathbb{Z} \oplus \mathbb{Z}_5$; $s$-biq $(s^2 - 2s + 2)$; ABQ torsion $(5, 13, 65)$; $\Phi_3 = (3, 3)$ |
| $11_2$ | `EXACT_SINGLETON` | $S^2$ | Yes | $2$ | $0$ | $(9, 2)$ | $4\text{–}9$ | $\langle x_1, x_2 \mid x_2^{-1} x_1^{-1} x_2 x_1^{-1} x_2^{-1} x_1^{-1} x_2 x_1 x_2^{-1} x_1 = 1 \rangle$ | $(x^3 - x^2 + 2x - 1)$ | `D10=30, F21=63, A5=180`; kei $\mathbb{Z} \oplus \mathbb{Z}_5$; $s$-biq $(s^3 - 2s^2 + s - 1)$; ABQ torsion $(5, 7, 35)$; $\Phi_3 = (3, 3)$ |
| $11_1^{-2}$ | `EXACT_SINGLETON` | $K^2$ | No | $0$ | $0$ | $(7, 4)$ | $4\text{–}10$ | $\pi_1 \cong 2O$ (binary octahedral group of order $48$, meridian order $8$) | $(3x,\; x + 1)$ | **First prime knotted $K^2$ in $\operatorname{ch}$-order**: `S3=10, S4=58, A5=76`; $\Sigma_2$: $H_1 = \mathbb{Z}_3, b_2 = 2$; kei $\mathbb{Z} \oplus \mathbb{Z}_3$ |

#### Mathematical Highlights (Part A)
1. **Non-Classical Asymmetric Alexander 2-Knots ($11_1$ and $11_2$) & Explicit Packet Group Isomorphism on `11_1.{a,b,c}`**: Every spun or twist-spun classical knot $\tau^k(K) \subset \mathbb{R}^4$ inherits a reciprocal Alexander ideal ($\Delta(x^{-1}) \doteq \Delta(x)$) or torsion quotient of $\Delta_K(x)$. Both $11_1$ ($\Delta(x) = 2x^2 - 2x + 1$, non-monic, $|\Delta(-1)| = 5$) and $11_2$ ($\Delta(x) = x^3 - x^2 + 2x - 1$, monic cubic and asymmetric, $|\Delta(-1)| = 5$) have single length-10 2-generator Wirtinger relators and are separated by $E_1$, $|\operatorname{Hom}(\pi_1, F_{21})| = 21$ vs $63$, and Alexander biquandle torsion $13$ vs $7$. Moreover, for the three unmerged $\le 13$-vertex move components `11_1.{a,b,c}`:
   - `11_1.a` and `11_1.b` ($G_b = \langle a, b \mid (b a b^{-1} a^{-1})^{-1} a (b a b^{-1} a^{-1}) = b \rangle$, i.e., $R(b a b^{-1} a^{-1})$) differ by a single Nielsen generator substitution; and
   - `11_1.b` and `11_1.c` ($G_c = \langle a, b \mid (b a b^{-1} a)^{-1} a (b a b^{-1} a) = b \rangle$, i.e., $R(b a b^{-1} a)$) have **isomorphic knot groups via mutually inverse meridional-reversing homomorphisms**:
     $$\phi\colon G_b \xrightarrow{\;\sim\;} G_c,\quad \phi(a) = a^{-1},\quad \phi(b) = (b^{-1} a)\,b^{-1}\,(a^{-1} b),$$
     $$\psi\colon G_c \xrightarrow{\;\sim\;} G_b,\quad \psi(a) = a^{-1},\quad \psi(b) = (a^{-1} b a)\,b^{-1}\,(a^{-1} b^{-1} a),$$
     sending every meridian to the inverse of a conjugate meridian (`test_soundness_and_ch12.py`). Whether `11_1.b` and `11_1.c` are ambient isotopic (or orientation-reversed ribbon cousins sharing the same group) remains a sharp open move-merging question!
2. **First Weakly Prime Knotted Klein Bottle in $\operatorname{ch}$-Order ($11_1^{-2}$, $\pi_1 \cong 2O$)**: While knotted $e = 0$ Klein bottles in $\mathbb{R}^4$ with arbitrary positive even meridian order were constructed geometrically by Yoshikawa (1998) (with smooth irreducibility of the order-4 example proved by Meier 2025), Yoshikawa's $\operatorname{ch} \le 10$ table contains **zero** prime connected Klein bottles (the unknotted $K^2(e=0) = P_+ \# P_-$ has $\operatorname{ch} = 3$ and is composite). At $\operatorname{ch} = 11$, $11_1^{-2}$ appears as an `EXACT_SINGLETON` with $e(F) = 0$ and $E_1 = (3, x + 1) \subset \mathbb{Z}[x]/(x^2 - 1)$. Exact coset enumeration shows that $\pi_1(S^4 \setminus 11_1^{-2})$ is finite of **order $48$**, isomorphic to the **binary octahedral group $2O \cong \operatorname{SL}(2, \mathbb{F}_3) \cdot 2$**, in which the meridian has **order $8$**! Because any connected-sum decomposition $F_0 \# P_\pm$ or $F_0 \# K^2_{\text{std}}$ forces a meridian of order $2$ (from $\pi_1(S^4 \setminus P_\pm) \cong \mathbb{Z}_2$), $11_1^{-2}$ admits **no** unknotted $\mathbb{R}P^2$ or $K^2$ summand, and full Yoshikawa Definition 5.1 verification across all $\operatorname{ch}(A), \operatorname{ch}(B) \le 10$ connected sums confirms that $11_1^{-2}$ is **weakly prime** (`S4 = 58, A5 = 76` vs `S4 = 34, A5 = 76` on $K \# K^2_{\text{std}}$ for $K \in \{8_1, 9_1, 10_2\}$).

---

### Part B: Multi-Component Surface-Links ($c \in \{2, 3\}$; 7 certified new weakly prime classes)

By extracting proper component subdiagrams and reducing under Yoshikawa moves, every $S^2$ component reduces to $\operatorname{ch} = 0$ ($0_1$), every $T^2$ component reduces to $\operatorname{ch} \in \{2, 4, 7\} < 10$ (smoothly unknotted $2_1^1$), and every $K^2(e=0)$ component reduces below the $\operatorname{ch} = 8$ threshold of knotted surfaces (smoothly unknotted $P_+ \# P_-$). In particular, exact sublink invariant profiles (`exact_sublink_invariant_profile`) separate $11_2^{0,1}$ from the $\operatorname{ch} = 14$ composite $6_1^{0,1} \#_{S^2} 8_1$ (whose $S^2$ sublink is $8_1 \neq 0_1$) and separate $11_1^{0,-2}$ from $7_1^{0,-2} \#_{S^2} 8_1$. Thus **every individual component across all seven new $\operatorname{ch} = 11$ links is smoothly unknotted**, and every new multi-component link is an **`EXACT_SINGLETON`** ($m = 1$).

| Label | Components | Ori. | $\chi$ | $(c, h)$ | $b(F)$ | Dehn-Simplified $\pi_1(\mathbb{R}^4 \setminus F)$ | Verified $E_1(F) \subset \mathbb{Z}[H_1]$ | Distinguishing Machine Invariants |
| :--- | :--- | :---: | :---: | :---: | :---: | :--- | :--- | :--- |
| $11_1^{0,1}$ | $S^2 \sqcup T^2$ | Yes | $2$ | $(9, 2)$ | $4\text{–}9$ | $\langle x_1, x_2 \mid x_2^{-1} x_1^{-1} x_2 x_1^{-1} x_2^{-1} x_1 x_2^{-1} x_1 x_2 x_1^{-1} x_2 x_1 = 1 \rangle$ | $(x - 1,\; y - 1)(xy - x - y)$ | `S3=24, A4=72, F21=147, S4=192, A5=600`; `sl2=264/93`; $\Sigma_2$: $H_1 = \mathbb{Z}_6$ |
| $11_2^{0,1}$ | $S^2 \sqcup T^2$ | Yes | $2$ | $(9, 2)$ | $4\text{–}9$ | $\langle x_1, x_2 \mid x_2^{-1} x_1^{-1} x_2 x_1^{-1} x_2^{-1} x_1^{-1} x_2 x_1 x_2^{-1} x_1 x_2 x_1 = 1 \rangle$ | $(y^3 - 2y^2 + 2y - 1,\; xy^2 - xy + x - y^2 + y - 1)$ | `S3=30, A4=120, F21=105, S4=384, A5=1500`; `sl2=456/141`; $\Sigma_2$: $H_1 = \mathbb{Z}_6$ |
| $11_1^{0,-2}$ | $S^2 \sqcup K^2$ | Mixed | $2$ | $(9, 2)$ | $4\text{–}9$ | $\langle x_1, x_2 \mid x_2^{-1} x_1^{-1} x_2^{-1} x_1 x_2 x_1^{-1} x_2^{-1} x_1 x_2^{-1} x_1^{-1} x_2 x_1 = 1 \rangle$ | $(3y - 3,\; 3xy + 3y,\; xy + x + y + 1)$ | `S3=30, Q8=40, Dic16=112, S4=384, A5=1500`; $\Sigma_2$: $H_1 = \mathbb{Z}_6$; $|\operatorname{Col}_{Q_6}| = 24$ |
| $11_2^{0,-2}$ | $S^2 \sqcup K^2$ | Mixed | $2$ | $(9, 2)$ | $4\text{–}9$ | $\langle x_1, x_2 \mid x_2^{-2} x_1^{-1} x_2^2 x_1^{-1} = 1 \rangle$ | $(x^2 + y,\; (x + 1)(y - 1))$ | `S3=18, Q8=16, Dic16=64, S4=168, A5=300`; **$\Sigma_2$: $(b_1=1, b_2=2)$**; kei $\mathbb{Z}^2$ |
| $11_3^{0,-2}$ | $S^2 \sqcup K^2$ | Mixed | $2$ | $(9, 2)$ | $4\text{–}9$ | $\langle x_1, x_2 \mid x_2^{-1} x_1^{-1} x_2^{-1} x_1^{-1} x_2 x_1 x_2 x_1^{-1} = 1 \rangle$ | $(x^2 + y,\; (x - 1)(y - 1))$ | `S3=18, Q8=16, Dic16=64, S4=168, A5=300`; **$\Sigma_2$: $\mathbb{Q}HS^4$ ($H_1 = \mathbb{Z}_4$)**; kei $\mathbb{Z} \oplus \mathbb{Z}_4$ |
| $11_1^{1,-2}$ | $T^2 \sqcup K^2$ | Mixed | $0$ | $(7, 4)$ | $6\text{–}10$ | $\langle x_1, x_2 \mid x_2^{-1} x_1^{-1} x_2 x_1^{-1} = 1,\; x_1^{-4} = 1 \rangle \cong \mathbb{Z}_4 \rtimes \mathbb{Z}$ | $(4y,\; 3y + 1,\; xy + 1)$ | `S3=12, Q8=40, Dic16=80, S4=96, A5=120`; $\Sigma_2$: $H_1 = \mathbb{Z}_2, b_2 = 2$; $|\operatorname{Col}_{Q_6}| = 24$ |
| $11_1^{0,0,-2}$ | $S^2 \sqcup S^2 \sqcup K^2$ | Mixed | $4$ | $(9, 2)$ | $5\text{–}10$ | $\langle x_1, x_2, x_3 \mid (x_2^{-1} x_3)^{-1} x_1^{-1} (x_2^{-1} x_3) = x_1 \rangle \cong K * \mathbb{Z}$ | $(0)$ | **Non-orientable twin of $10_1^{0,0,1}$**: $|\operatorname{Col}_{Q_6}| = 96 \neq 144$; $\mathbb{Z}_3$-cover $b_2 = 1 \neq 3$; both $(K^2, S_i^2)$ sublinks are $7_1^{0,-2}$ |

#### Mathematical Highlights (Part B)
1. **The Algebraic Twin Pair $11_2^{0,-2}$ vs $11_3^{0,-2}$**: Both links share identical homomorphism counts into every test group (`S3, D8, Q8, D10, A4, Dic12, Dic16, F21, S4, SL2F3, A5`) and identical single-variable Alexander polynomial $\Delta(t) = t^2 + 1$, yet are separated independently at Layer 3 (meridian-pinned `E1` over $\mathbb{Z}[x^{\pm 1}, y]/(y^2 - 1)$ with $x = [\mu_{S^2}]$ and $y = [\mu_{K^2}]$), Layer 4 (**double branched cover $\Sigma_2$**: $(b_1 = 1, b_2 = 2)$ vs a rational homology 4-sphere with $H_1 = \mathbb{Z}_4, b_2 = 0$), and Layer 7 (universal kei $\mathbb{Z}^2$ vs $\mathbb{Z} \oplus \mathbb{Z}_4$).
2. **Tier-3 Discovery of $11_1^{0,0,-2}$ (Masquerading as $0_1 \sqcup 7_1^{0,-2}$)**: In $11_1^{0,0,-2}$, two unlinked 2-spheres $S_1^2, S_2^2$ both thread the orientation-reversing loop of a single Klein bottle $K^2$ (making $11_1^{0,0,-2}$ the exact non-orientable analogue of Yoshikawa's lone 3-component prime link $10_1^{0,0,1}$). Because the orientation-preserving loop of $K^2$ is unthreaded, $\pi_1(\mathbb{R}^4 \setminus 11_1^{0,0,-2}) = \langle x_1, x_2, x_3 \mid (x_2^{-1} x_3)^{-1} x_1^{-1} (x_2^{-1} x_3) = x_1 \rangle$ is **isomorphic as an abstract group** (via the non-meridional Nielsen change of basis $u = x_2^{-1} x_3$) to $K * \mathbb{Z} \cong \pi_1(\mathbb{R}^4 \setminus (0_1 \sqcup 7_1^{0,-2}))$, rendering all unbased group, Alexander, kei/biquandle, and $\Sigma_2$ invariants identical to the split union $0_1 \sqcup 7_1^{0,-2}$. Our Tier-3 peripheral and sheet invariants separate them decisively:
   - **Regular Cover Peripheral Betti Numbers (Theorem 5.5)**: Across $\mathbb{Z}_2, \mathbb{Z}_3, S_3, Q_8, A_4, \operatorname{Dic}_{16}$ covers, $\partial \widetilde{E}_\rho \to \widetilde{E}_\rho$ has strictly higher peripheral cycle rank on $11_1^{0,0,-2}$ than on $0_1 \sqcup 7_1^{0,-2}$ (e.g., $b_2 = 1$ vs $b_2 = 3$ in $\mathbb{Z}_3$-covers).
   - **Kamada–Oshiro Octahedral Symmetric Quandle $(Q_6, *, \rho)$**: $|\operatorname{Col}_{(Q_6, \rho)}(11_1^{0,0,-2})| = 96$, whereas multiplicativity under split union forces $|\operatorname{Col}_{(Q_6, \rho)}(0_1 \sqcup 7_1^{0,-2})| = 6 \times 24 = 144$.
   - **Proper Sublinks**: Both $(K^2, S_1^2)$ and $(K^2, S_2^2)$ reduce to non-split $7_1^{0,-2}$, while $(S_1^2, S_2^2)$ reduces to $0_1 \sqcup 0_1$, ruling out both splittability and connected-sum decomposability.

---

## 3. Tier-3 4D Invariants, Two-Tier Cell Status, & Separating Witnesses ($\operatorname{ch} = 11$)

| Label | `cell_status` ($m$) | $\Sigma_2(S^4, F)$ $(\chi, b_1, \operatorname{Tors} H_1, b_2, |\sigma|)$ | Cover Betti Highlight $(b_1, \operatorname{Tors} H_1, b_2, b_3)$ | $\Phi_{(Q_6, \rho)}$ | Nearest Catalog Neighbor & Separating Layer |
| :--- | :--- | :--- | :--- | :---: | :--- |
| $11_1$ | `NEW_PACKET` ($3$) | $(2, 0, \mathbb{Z}_5, 0, 0)$ ($\mathbb{Q}HS^4$) | $\mathbb{Z}_2$: $(1, \mathbb{Z}_5, 0, 0)$; $\mathbb{Z}_3$: $(1, \mathbb{Z}_{13}, 0, 0)$ | `6 + 0u` | $11_2$ at Layer 3 (`E1`) |
| $11_2$ | `EXACT_SINGLETON` ($1$) | $(2, 0, \mathbb{Z}_5, 0, 0)$ ($\mathbb{Q}HS^4$) | $F_{21}$: $(1, \mathbb{Z}_2^5 \oplus \mathbb{Z}_{14}, 6, 6)$; $\mathbb{Z}_3$: $(1, \mathbb{Z}_7, 0, 0)$ | `6 + 0u` | $11_1$ at Layer 3 (`E1`) |
| $11_1^{-2}$ | `EXACT_SINGLETON` ($1$) | $(4, 0, \mathbb{Z}_3, 2, 0)$ | $S_3$: $(0, \mathbb{Z}_2^2, 13, 2)$; $\mathbb{Z}_2$: $(0, \mathbb{Z}_3, 3, 0)$ | `0` | $11_1$ at Layer 1 (`topology`) |
| $11_1^{0,1}$ | `EXACT_SINGLETON` ($1$) | $(2, 0, \mathbb{Z}_6, 0, 0)$ ($\mathbb{Q}HS^4$) | $S_3$: $(4, \mathbb{Z}_3, 8, 5)$; $F_{21}$: $(8, \mathbb{Z}_7, 20, 13)$ | `12 + 0u` | $11_2^{0,1}$ at Layer 3 (`E1`) |
| $11_2^{0,1}$ | `EXACT_SINGLETON` ($1$) | $(2, 0, \mathbb{Z}_6, 0, 0)$ ($\mathbb{Q}HS^4$) | $S_3$: $(6, 0, 9, 4),\; (6, 0, 10, 5)$ | `12 + 0u` | $11_1^{0,1}$ at Layer 3 (`E1`) |
| $11_1^{0,-2}$ | `EXACT_SINGLETON` ($1$) | $(2, 0, \mathbb{Z}_6, 0, 0)$ ($\mathbb{Q}HS^4$) | $Q_8$: $(2, \mathbb{Z}_3^2, 3, 2)$; $\operatorname{Dic}_{16}$: $(2, \mathbb{Z}_3^2, 5, 4)$ | `24 + 0u` | $11_2^{0,-2}$ at Layer 3 (`E1`) |
| $11_2^{0,-2}$ | `EXACT_SINGLETON` ($1$) | $(2, 1, 0, 2, 0)$ | $\operatorname{Dic}_{16}$: $(4, 0, 7, 4)$ | `0` | $11_1^{0,-2}$ at Layer 3 (`E1`) |
| $11_3^{0,-2}$ | `EXACT_SINGLETON` ($1$) | $(2, 0, \mathbb{Z}_4, 0, 0)$ ($\mathbb{Q}HS^4$) | $\operatorname{Dic}_{16}$: $(4, 0, 8, 5)$ | `0` | $11_1^{0,-2}$ at Layer 3 (`E1`) |
| $11_1^{1,-2}$ | `EXACT_SINGLETON` ($1$) | $(4, 0, \mathbb{Z}_2, 2, 0)$ | $Q_8$: $(1, 0, 18, 2)$ | `24 + 0u` | $9_1^{1,-2}$ at Layer 3 (`E1`) |
| $11_1^{0,0,-2}$ | `EXACT_SINGLETON` ($1$) | $(0, 1, \mathbb{Z}_2, 0, 0)$ | $\mathbb{Z}_3$: $(4, \mathbb{Z}_2^3, 1, 4)$; $Q_8$: $(10, 0, 2, 9)$ | `96 + 0u` | $10_1^{0,0,1}$ at Layer 1 (`topology`); $0_1 \sqcup 7_1^{0,-2}$ at **Layer 6 (`pi2_covers`)** |

---

## 4. Canonical EPD Codes for Certified New $\operatorname{ch} = 11$ Classes

| Label | Canonical EPD |
| :--- | :--- |
| $11_1$ (`11_1.a`) | `X[1,2,3,4], X[4,5,6,1], X[6,7,8,9], Y[9,10,11,2], X[12,13,14,3], X[15,7,5,14], X[16,17,10,8], X[11,18,19,12], Z[13,20,21,15], X[21,22,17,16], X[20,19,18,22]` |
| $11_2$ | `X[1,2,3,4], X[4,5,6,1], X[7,8,9,6], X[10,3,2,9], Y[10,11,12,13], X[13,14,15,5], X[16,15,14,17], X[8,18,19,20], Y[20,21,12,11], X[22,17,21,19], X[16,22,18,7]` |
| $11_1^{-2}$ | `X[1,2,3,4], X[4,5,6,1], Y[6,7,8,9], X[9,10,11,2], Y[3,12,13,14], X[15,16,5,14], Y[16,17,8,7], X[17,18,19,20], X[21,11,10,20], Y[21,22,13,12], X[22,19,18,15]` |
| $11_1^{0,1}$ | `X[1,2,3,4], X[4,5,6,1], X[7,8,2,6], Y[3,9,10,11], Z[11,10,12,13], X[14,15,16,13], X[17,7,5,16], X[18,19,9,8], X[19,18,20,12], X[20,17,21,22], X[22,21,15,14]` |
| $11_2^{0,1}$ | `X[1,2,3,4], X[4,5,6,1], X[7,3,2,8], Z[9,8,10,11], X[13,14,15,12], X[16,17,14,13], X[18,19,7,9], X[12,20,19,18], Y[11,10,6,16], X[20,15,22,21], X[17,5,21,22]` |
| $11_1^{0,-2}$ | `X[1,2,3,4], X[4,5,6,1], X[7,8,2,6], X[9,10,11,3], X[10,12,5,11], X[13,14,8,7], Z[15,13,16,17], X[18,19,20,21], X[19,22,9,20], X[15,18,21,14], Y[17,16,12,22]` |
| $11_2^{0,-2}$ | `X[1,2,3,4], X[4,5,6,1], X[6,7,8,9], X[10,11,9,12], X[13,3,2,11], X[14,15,5,16], X[17,8,7,15], Y[17,14,18,19], Z[19,18,20,12], X[21,20,16,22], X[10,21,22,13]` |
| $11_3^{0,-2}$ | `Y[1,2,3,4], Z[4,5,6,1], X[6,7,8,2], X[9,10,11,3], X[12,13,14,11], X[13,15,16,14], X[15,17,5,16], X[17,12,18,7], X[8,18,19,20], X[20,19,21,22], X[22,21,10,9]` |
| $11_1^{1,-2}$ | `Z[1,2,3,4], Y[4,5,6,1], X[6,7,8,2], X[3,9,10,11], Z[11,12,13,14], Y[14,13,15,16], X[16,17,18,5], X[17,19,7,18], X[20,10,9,8], X[20,21,22,12], X[15,22,21,19]` |
| $11_1^{0,0,-2}$ | `Y[1,2,3,4], Y[5,6,7,8], X[1,4,6,5], X[9,10,11,3], X[7,11,12,13], X[8,13,14,15], X[15,16,9,2], X[17,18,19,10], X[12,19,18,20], X[14,20,21,22], X[22,21,17,16]` |

*(Move-component sub-IDs `11_1.a, 11_1.b, 11_1.c` and their EPDs are recorded in `move_component_ids` and `move_component_epds` inside [`expanded_ch11.json`](../code/data/expanded_ch11.json).)*

---

## 5. Stratification of the $93$ `AMBIGUOUS_SHADOW` Residual Diagram Collisions

Because 4D diagram simplification lacks a monotonic Reidemeister theorem (hard surface-unknot `ch`-diagrams can require $> 2$ stabilizations to unknot; cf. Jabłonowski 2013/2017), we report every surviving 11-vertex diagram that resists $\le 13$-vertex post-flype mountain-pass simplification as `cell_status = "AMBIGUOUS_SHADOW"` in `residual_collisions` inside [`expanded_ch11.json`](../code/data/expanded_ch11.json). Notably, Dehn relator shortening (`Tietze I`) certifies algebraically that the lone $S^2$ unknot-profile residual has $\pi_1(\mathbb{R}^4 \setminus F) \cong \mathbb{Z}$ (`<x1 | >`), proving it **topologically unknotted** in $S^4$ by Freedman's Theorem (1982), while both $\mathbb{R}P^2$ residuals have $\pi_1 \cong \mathbb{Z}_2$ (`<x1 | x1^-2 = 1>`):

| Collision Stratum | Diagram Count | Breakdown by Matched Reference Class |
| :--- | ---: | :--- |
| **Higher-genus unknot profiles** (`unknot:*`) | $65$ | `#4RP^2(e=0)` ($22$), `K^2(e=0)` ($17$), `#3RP^2(e=±2)` ($8$), `Sigma_2(e=0)` ($7$), `#5RP^2(e=±2)` ($4$), `#4RP^2(e=±4)` ($3$), `K^2(e=±4)` ($2$), `#3RP^2(e=±6)` ($1$), `Sigma_3(e=0)` ($1$) |
| **Lower-`ch` Yoshikawa classes** (`ch<=10:*`) | $19$ | `2_1^1` ($12$), `2_1^{-1}` ($2$, $\pi_1 \cong \mathbb{Z}_2$), `0_1` ($1$, $\pi_1 \cong \mathbb{Z}$, topologically unknotted by Freedman), `10_2` ($1$), `10_1^1` ($1$), `9_1^{0,1}` ($1$), `10_2^{0,-2}` ($1$) |
| **Composite classes** (`composite:*`) | $8$ | `8_1 # K^2(e=0)` ($4$), `2_1^1 # 9_1^{0,1}` ($2$), `2_1^1 # 6_1^{0,1}` ($1$), `2_1^{-1} # 8_1` ($\pi_1 \cong S_3$, $1$) |
| **Split unions** (`split:*`) | $1$ | `0_1 u 2_1^1` ($1$) |
| **Total `AMBIGUOUS_SHADOW` residuals** | **$93$** | Full EPDs, `collision_matches`, and Tier 1–3 records preserved in `expanded_ch11.json["residual_collisions"]` |

---

## 6. Concentrated Bridge Trisection Cross-Checks ($b \le 4$)

Alongside the `ch = 11` marked-diagram enumeration, `triplane.py` independently verifies low-bridge Meier–Zupan $(b; p_{12}, p_{23}, p_{31})$ tri-plane diagrams via constructive Reidemeister unlink movies on all three pairwise unions $(T_1 \cup \overline{T}_2, T_2 \cup \overline{T}_3, T_3 \cup \overline{T}_1)$ (recorded in `census_stats.triplane_b_le_4_crosscheck` inside [`expanded_ch11.json`](../code/data/expanded_ch11.json)):
- **Crossingless $c_{\text{tri}} = 0$ ($b = 1, 2, 3$)**: Exhaustively recovers the unique $(1; 1, 1, 1)$ unknotted sphere ($\chi = 2$), $(2; 1, 1, 2)$ ($\chi = 2$) and $(2; 2, 2, 2)$ ($S^2 \sqcup S^2$, $\chi = 4$), and for $b = 3$ the standard $(3; 1, 1, 1)$ unknotted torus $2_1^1$ ($\chi = 0$) together with $\chi \in \{2, 4, 6\}$ unlinks — confirming Meier–Zupan's classification of crossingless $b \le 3$ trisections.
- **Concentrated $c_{\text{tri}} \le 2$ ($b = 2, 3$)**: Exhaustive braid insertions in $T_1$ recover $(2; 1, 1, 1)$ at $c_{\text{tri}} = 1$ ($\chi = 1$, standard projective plane $2_1^{-1}$ realizing $c_{\text{tri}}(P_\pm) = 1$) and $(3; 1, 1, 2)$ at $c_{\text{tri}} = 1$ ($\chi = 1$, stabilized $P_\pm$) plus $(3; 2, 2, 2)$ ($\chi = 3$, $S^2 \sqcup P_\pm$), with no knotted surfaces below $b = 4$.
- **Bridge $b = 4$ Spun Trefoil $(4; 2, 2, 2)$**: Certifies $6$ explicit 6-crossing concentrated $(4; 2, 2, 2)$ tri-plane models ($\chi = 6 - 4 = 2$, $S^2$) matching $c_{\text{tri}}(8_1) \le 6$ alongside Allred et al. (2023) and Gong–Lewis-Monkman–Osnes (2026).
