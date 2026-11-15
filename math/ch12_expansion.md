# Complete $\operatorname{ch}(F) = 12$ Census, Spoke Promotions, and Two Cappell–Shaneson $2$-Knots

Extending the bottom-up marked-diagram census from $\operatorname{ch} = 11$ ($3{,}378$ prime shadows) to $\operatorname{ch} = 12$ crosses an $18.8\times$ jump in raw configuration space: **$15{,}882$** prime 4-regular spherical shadows (`OEIS A113201(12)`, dual to simple 14-vertex planar quadrangulations) across $15{,}882 \times 4^{12} \approx 2.665 \times 10^{11}$ raw vertex assignments.

Using the 64-core C99 + OpenMP enumeration and mountain-pass tunneling engine ([`_fast_enum.c`](../code/src/surfaceknots/_fast_enum.c), [`ch12_census.py`](../code/src/surfaceknots/ch12_census.py)) together with our audit-hardened 10-layer 4D invariant filtration (including meridian-pinned strong Gröbner Fox ideals over $\mathbb{Z}[H_1]$, Mochizuki $R_5$ dihedral 3-cocycle state sums, and regular $F_{21}$ & $F_{57}$ cover Betti profiles), we scanned **every non-extreme saddle sector $m \in \{1, \dots, 11\}$** ($m = 0$ having zero saddles, hence an $S^2$ unlink, and $m = 12$ having zero classical crossings, hence an unknotted planar surface in $\mathbb{R}^3$ by Theorem 1.1 of [`complexity_measures.md`](./complexity_measures.md)). Full machine-verified outputs are exported in [`ch12_raw_survivors.json`](../code/data/ch12_raw_survivors.json), [`expanded_ch12.json`](../code/data/expanded_ch12.json), and [`expanded_ch12.csv`](../code/data/expanded_ch12.csv):
- **34 certified new, weakly prime, non-splittable $\operatorname{ch} = 12$ invariant classes** across **$63$ diagrams**: **$20$ `EXACT_SINGLETON`** classes ($m = 1$) and **$14$ `NEW_PACKET`** classes ($m \in \{2, 3, 4, 5\}$), comprising $11$ connected $S^2$-knots (`12_1`–`12_11`), $2$ connected $T^2$-knots (`12_1^1`, `12_2^1`), $0$ connected $\mathbb{R}P^2$ or $K^2$ prime knots, $19$ two-component surface-links, and $2$ three-component surface-links (`12_1^{0,0,1}`, `12_1^{0,1,1}`), conditional on completeness of the $\operatorname{ch} \le 11$ reference tables; and
- **55 documented `AMBIGUOUS_SHADOW` collision cells** across **$742$ residual diagrams** ($607$ matching higher-genus standard unknots or unknot-composites, $90$ matching $\operatorname{ch} \le 10$ Yoshikawa classes, $40$ matching pure composites, $4$ matching $\operatorname{ch} = 11$ classes, and $1$ split union), after **$104$** of the $909$ raw C99-tunneled survivors are pruned/merged during 4D Tait flypes, 4-edge cyclic-tube checks, and post-flype C99 tunneling (`post_flype_unique_diagrams = 805`).

---

## 1. Complete $\operatorname{ch} = 12$ ($m = 1 \dots 11$) Enumeration & Tunneling Funnel

> [!IMPORTANT]
> **Pre-tunnel vs. post-tunnel certification semantics:** In `per_m_stats` below, `raw_cands`, `components`, and `level_kept` count prime-shadow diagrams surviving marked-bigon ($\Omega_6'$) pruning, linking-parity checks, and the two-stage $O(c \cdot 2^c)$ Fast Kauffman Bracket Butterfly Transform $\bmod (10^9 + 7)$ on $(L_-(D), L_+(D))$ (a strictly necessary condition for classical unlinkedness). Every survivor after C99 mountain-pass tunneling (`intra_m_tunneled_kept = 1,261` and global cross-$m$ `cyclic_tube_filtered_kept = 909`) carries **explicit constructive Reidemeister $\mathrm{R1}/\mathrm{R2}/\mathrm{R3}$ unlink certificates** on both $L_-(D)$ and $L_+(D)$ (`all_tunneled_survivors_reidemeister_unlink_certified = True`).

### 1.1 Per-Saddle Sector Breakdown across All $15{,}882$ Prime Shadows (`ch12_raw_survivors.json`)

| Saddle Count $m$ | Classical $c = 12 - m$ | Bracket-Filtered `raw_cands` | Level $\Omega_{4,5,7,8}$ `components` | Level-Minimal `kept` | Intra-$m$ Tunneled (`+1/+2`) | Wall Time ($64$ cores) |
| :---: | :---: | ---: | ---: | ---: | ---: | ---: |
| $1$ | $11$ | $37{,}576$ | $17{,}726$ | $730$ | **$0$** | $3.81\,\text{s}$ |
| $2$ | $10$ | $730{,}771$ | $279{,}699$ | $13{,}653$ | **$211$** | $41.07\,\text{s}$ |
| $3$ | $9$ | $6{,}576{,}533$ | $1{,}973{,}739$ | $32{,}160$ | **$203$** | $104.56\,\text{s}$ |
| $4$ | $8$ | $39{,}481{,}522$ | $12{,}042{,}142$ | $95{,}683$ | **$612$** | $580.80\,\text{s}$ |
| $5$ | $7$ | $137{,}856{,}031$ | $44{,}501{,}987$ | $101{,}641$ | **$120$** | $800.64\,\text{s}$ |
| $6$ | $6$ | $285{,}896{,}789$ | $98{,}423{,}366$ | $84{,}222$ | **$104$** | $1{,}218.76\,\text{s}$ |
| $7$ | $5$ | $340{,}907{,}322$ | $125{,}352{,}072$ | $26{,}994$ | **$5$** | $1{,}110.98\,\text{s}$ |
| $8$ | $4$ | $232{,}313{,}665$ | $90{,}888{,}537$ | $5{,}077$ | **$6$** | $721.81\,\text{s}$ |
| $9$ | $3$ | $91{,}547{,}119$ | $37{,}588{,}745$ | $58$ | **$0$** | $199.98\,\text{s}$ |
| $10$ | $2$ | $20{,}672{,}204$ | $8{,}682{,}310$ | $1$ | **$0$** | $54.03\,\text{s}$ |
| $11$ | $1$ | $2{,}397{,}072$ | $1{,}006{,}026$ | $0$ | **$0$** | $9.76\,\text{s}$ |
| **Total ($m = 1 \dots 11$)** | **$1 \dots 11$** | **$1{,}158{,}416{,}604$** | **$420{,}756{,}349$** | **$360{,}219$** | **$1{,}261$** | **$4{,}846.20\,\text{s}$** |

Notice empirically that **every** diagram in the extreme sectors $m \in \{1, 9, 10, 11\}$ collapses completely to `0` survivors under intra-$m$ $+1/+2$ tunneling, confirming that all non-collapsible $\operatorname{ch} = 12$ prime-shadow survivors reside in $m \in \{2, \dots, 8\}$. Combining all $1{,}261$ intra-$m$ survivors across $m = 1 \dots 11$, applying cross-$m$ C99 mountain-pass tunneling (`194.19s`) and parallel 4-edge cyclic-tube pruning yields **$909$** Reidemeister-certified post-tube survivors (`ch12_raw_survivors.json`). Subsequent 4D Tait flype simplification (`10` directly reduced, `3` split/cut), 4-edge cyclic-tube checks (`12` pruned), canonical deduplication (`855` pre-post-flype-tunnel diagrams), and post-flype C99 $+1/+2$ tunneling (`depth 384/256`, `50` pruned) eliminate **$104$** diagrams in total, leaving **$805$** post-flype unique diagrams (`63` certified new prime diagrams across `34` invariant classes + `742` residual diagrams across `55` `AMBIGUOUS_SHADOW` cells) after full Yoshikawa Definition 5.1 composite screening ($\operatorname{ch}(F_1) \le 11, \operatorname{ch}(F_2) \le 11$, summand diagram sum up to $22$) in [`expanded_ch12.json`](../code/data/expanded_ch12.json).

### 1.2 Summary Table of All $34$ Certified New Weakly Prime $\operatorname{ch} = 12$ Classes

#### Part A: Connected Surface-Knots ($c = 1$; $11$ $S^2$-knots + $2$ $T^2$-knots, $0$ prime $\mathbb{R}P^2$ or $K^2$ knots)

| Label | `cell_status` ($m$) | Type | $(c, h)$ | Verified $E_1(F) \subset \mathbb{Z}[H_1]$ | $\Sigma_2$: $(\operatorname{Tors} H_1, b_2)$ | Distinguishing 4D Invariants & Realizations |
| :--- | :--- | :---: | :---: | :--- | :---: | :--- |
| $12_1$ | `NEW_PACKET` ($3$) | $S^2$ | $(10, 2)$ | $(2x^2 - 3x + 2)$ | $(\mathbb{Z}_7, 0)$ | **`Spin(5_2)`** ($\tau^0(5_2)$, `MOVE_CERTIFIED`); `A5=180` |
| $12_2$ | `NEW_PACKET` ($2$) | $S^2$ | $(10, 2)$ | $(3,\; x^2 + 1)$ | $(0, 0)$ | **Non-ribbon $\mathbb{F}_9$ 2-knot** ($\mathbb{Z}HS^4$ $\Sigma_2$; excluded from $\tau^k(K)$ for $2 \le k \le 7$) |
| $12_3$ | `EXACT_SINGLETON` ($1$) | $S^2$ | $(10, 2)$ | $(3x^2 - 3x + 1)$ | $(\mathbb{Z}_7, 0)$ | `A4=36, F21=63, S4=48, A5=180`; $F_{21}$-cover $(1, \mathbb{Z}_2^2 \oplus \mathbb{Z}_4^2 \oplus \mathbb{Z}_{28}, 6, 6)$ |
| $12_4$ | `NEW_PACKET` ($5$) | $S^2$ | $(10, 2)$ | $(5,\; x + 1)$ | $(\mathbb{Z}_5, 0)$ | **$\tau^2(4_1)$ cohort** (`D10=30, A5=180`); Mochizuki $R_5$: `5 + 10*u^2 + 10*u^3` |
| $12_5$ | `EXACT_SINGLETON` ($1$) | $S^2$ | $(10, 2)$ | $(5,\; x + 1)$ | $(\mathbb{Z}_5, 0)$ | **$\tau^2(5_1)$ cohort** (`D10=30, A5=180`); Mochizuki $R_5$: `5 + 10*u + 10*u^4` |
| $12_6$ | `NEW_PACKET` ($2$) | $S^2$ | $(10, 2)$ | $(7,\; x + 1)$ | $(\mathbb{Z}_7, 0)$ | Non-principal determinant-$7$ ($\tau^2(5_2)$ cohort; $\pi_1 \cong \mathbb{Z}_7 \rtimes_{-1} \mathbb{Z}$) |
| $12_7$ | `NEW_PACKET` ($4$) | $S^2$ | $(10, 2)$ | $(x^3 - 2x^2 + 3x - 1)$ | $(\mathbb{Z}_7, 0)$ | **`R(b^-1*a*b*a^-1*b*a)`** (`INVARIANT_MATCH`); $F_{57}$-cover $(1, \mathbb{Z}_{11}^2 \oplus \mathbb{Z}_{209}, 18, 18)$ |
| $12_8$ | `EXACT_SINGLETON` ($1$) | $S^2$ | $(10, 2)$ | $(x^3 - 2x^2 + 3x - 1)$ | $(\mathbb{Z}_7, 0)$ | **Cappell–Shaneson `CS_A1`** ($\mathbb{Z}^3 \rtimes_{A_1} \mathbb{Z}$, non-ribbon); $F_{57}$-cover $(1, \mathbb{Z}_{19}, 18, 18)$ |
| $12_9$ | `EXACT_SINGLETON` ($1$) | $S^2$ | $(10, 2)$ | $(x^3 - x^2 + 1)$ | $(0, 0)$ | **`R(b^2*a)`** (`INVARIANT_MATCH`, $\mathbb{Z}HS^4$ $\Sigma_2$); $F_{21}$-cover $(1, \mathbb{Z}_2^2 \oplus \mathbb{Z}_{14}, 6, 6)$ |
| $12_{10}$ | `EXACT_SINGLETON` ($1$) | $S^2$ | $(10, 2)$ | $(x^3 - x^2 + 2x - 1)$ | $(\mathbb{Z}_5, 0)$ | **Cappell–Shaneson `CS_A0`** ($\mathbb{Z}^3 \rtimes_{A_0} \mathbb{Z}$, non-ribbon); $F_{21}$-cover $(1, \mathbb{Z}_7, 6, 6)$ |
| $12_{11}$ | `EXACT_SINGLETON` ($1$) | $S^2$ | $(10, 2)$ | $(x^4 - x^3 + x^2 - x + 1)$ | $(\mathbb{Z}_5, 0)$ | **`Spin(5_1)`** ($\tau^0(5_1)$, `MOVE_CERTIFIED`); `D10=30, A5=420` |
| $12_1^1$ | `NEW_PACKET` ($4$) | $T^2$ | $(8, 4)$ | $(3,\; x + 1)$ | $(\mathbb{Z}_3, 2)$ | Non-principal $T^2$-knot (`S3=12, Dic12=24, S4=48, A5=120`; `pSL2=0`, $S_3$-cover $b_2 = 14$) |
| $12_2^1$ | `NEW_PACKET` ($4$) | $T^2$ | $(8, 4)$ | $(x^2 - 3x + 1)$ | $(\mathbb{Z}_5, 2)$ | **`Tube(4_1)`** (`MOVE_CERTIFIED`); `pSL2=24`, $A_4$-cover $b_2 = 30 \neq 33$ |

#### Part B: Multi-Component Surface-Links ($c \in \{2, 3\}$; $21$ classes across $34$ diagrams)

| Label | `cell_status` ($m$) | Components ($e$) | $(c, h)$ | Verified $E_1(F) \subset \mathbb{Z}[H_1]$ | Distinguishing Notes & Spoke Identifications |
| :--- | :--- | :--- | :---: | :--- | :--- |
| $12_1^{0,0}$ | `EXACT_SINGLETON` ($1$) | $S^2 \sqcup S^2$ $(0,0)$ | $(10, 2)$ | $(0)$ | First prime $S^2 \sqcup S^2$ 2-link in $\operatorname{ch}$-order! ($\Sigma_2$: $b_1 = 1, b_2 = 0$; $|\operatorname{Col}_{Q_6}| = 36$) |
| $12_1^{0,1}$–$12_6^{0,1}$ | $3$ `SING` + $3$ `PKT` | $S^2 \sqcup T^2$ $(0,0)$ | $(10, 2)$ | $6$ pairwise-distinct ideals | **`12_4^{0,1} = Mixed(L5a1)`** (`MOVE_CERTIFIED`, $m=4$); **`12_5^{0,1} = T(bab^{-1}a^{-1}b)`** ($m=2$) |
| $12_1^{1,1}$–$12_3^{1,1}$ | $3$ `NEW_PACKET` ($2,2,4$) | $T^2 \sqcup T^2$ $(0,0)$ | $(8, 4)$ | $(x-1,y-1)$; $(x-1,y-1)^2$; $(x-1,y-1)(x+y)$ | **`12_3^{1,1} = Tube(L4a1)`** (`MOVE_CERTIFIED`, $m=4$); $12_1^{1,1}$ separated from $10_1^{1,1}$ at Layer 6 (`pi2_covers`) |
| $12_1^{0,-2}$–$12_3^{0,-2}$ | $2$ `SING` + $1$ `PKT` | $S^2 \sqcup K^2$ $(0,0)$ | $(10, 2)$ | $3$ ideals ($\Sigma_2$: $\mathbb{Z}_6$) | $12_1^{0,-2} = K(baba^{-1}b^{-1})$ ($m=1$), $12_2^{0,-2} = K(bab^{-1}ab)$ ($m=2$, sep. from $11_1^{0,-2}$ by `S4=216` vs `384`), $12_3^{0,-2}$ ($m=1$) |
| $12_1^{1,-2}$ | `EXACT_SINGLETON` ($1$) | $T^2 \sqcup K^2$ $(0, \pm 4)$ | $(8, 4)$ | $(2y,\; y + 1,\; xy + 1)$ | Extremal normal Euler $|e(K^2)| = 4$ on the Klein-bottle component! |
| $12_1^{-2,-2}$–$12_3^{-2,-2}$ | $3$ `EXACT_SINGLETON` ($1$) | $K^2 \sqcup K^2$ | $(8, 4)$ | $(2xy, xy+1, xy+x)$; $(4xy, \dots)$ | $12_1^{-2,-2}$ has $e = (-4, +4)$; $12_2^{-2,-2}, 12_3^{-2,-2}$ have $e = (0, 0)$ (separated at `E1`) |
| $12_1^{-1,-1}, 12_2^{-1,-1}$ | $2$ `EXACT_SINGLETON` ($1$) | $\mathbb{R}P^2 \sqcup \mathbb{R}P^2$ $(-2, +2)$ | $(10, 2)$ | $(6xy,\; xy + x,\; xy + y)$ | Both certify $t(F) > 0$ via $(Q_6, \rho) = \texttt{0 + 24u}$; separated at Layer 6 ($Q_8$-cover torsion $\mathbb{Z}_3$ vs $\mathbb{Z}_3^3$)! |
| $12_1^{0,0,1}$ | `EXACT_SINGLETON` ($1$) | $T^2 \sqcup S^2 \sqcup S^2$ $(0,0,0)$ | $(10, 2)$ | $(0)$ | $\Sigma_2$: $(b_1 = 1, \operatorname{Tors} = \mathbb{Z}_4, b_2 = 0)$; $|\operatorname{Col}_{Q_6}| = 216$ |
| $12_1^{0,1,1}$ | `EXACT_SINGLETON` ($1$) | $T^2 \sqcup S^2 \sqcup T^2$ $(0,0,0)$ | $(8, 4)$ | $(y^2-2y+1, xy-x-y+1, yz-y-z+1)$ | $\Sigma_2$: $(b_1 = 0, \operatorname{Tors} = \mathbb{Z}_2^2, b_2 = 0)$; $|\operatorname{Col}_{Q_6}| = 24$ |

---

## 2. Two Cappell–Shaneson $2$-Knots at $\operatorname{ch} = 12, m = 2$ (`CS_A0` & `CS_A1`)

The single most striking mathematical discovery of the $\operatorname{ch} = 12$ census occurs in the minimal-saddle $(c, m) = (10, 2)$ connected $S^2$ stratum. Prior to adding regular $F_{21}$ and $F_{57}$ finite-cover homology to Layer 6 (`pi2_covers` in `tier3.py`), two non-ribbon **Cappell–Shaneson $2$-knots** (`12_10` and `12_8` in [`expanded_ch12.json`](../code/data/expanded_ch12.json)) were hidden inside invariant cells occupied by 1-fusion ribbon $2$-knots sharing the *exact same* Alexander module, double branched cover $\Sigma_2$, 3-fold cyclic cover $M_3$, kei/biquandle invariants, and small finite-group homomorphism counts.

### Theorem 2.1 (Discovery, Group Structure, Non-Ribbonness, and Frobenius-Cover Separation of `CS_A0` & `CS_A1`)
*Consider the two 12-vertex $(c, m) = (10, 2)$ admissible $S^2$-knot diagrams (`12_10` and `12_8`, both `EXACT_SINGLETON` cells):*
- **`12_10` (`CS_A0`)**:
  $$\texttt{X[1,2,3,4], X[5,6,4,7], X[6,8,9,1], Y[9,10,11,2], X[12,13,14,3], Y[14,15,16,7],}$$
  $$\texttt{X[5,17,18,8], X[17,19,10,18], X[11,20,21,12], X[13,21,22,15], X[23,24,19,16], X[24,23,22,20]}$$
- **`12_8` (`CS_A1`)**:
  $$\texttt{X[1,2,3,4], X[4,5,6,1], X[7,8,9,6], Y[9,10,11,2], X[3,12,13,14], Y[14,15,16,17],}$$
  $$\texttt{X[18,7,5,17], X[8,18,19,10], X[20,21,12,11], X[21,22,15,13], X[23,24,19,16], X[24,23,22,20]}$$

*Both diagrams satisfy the following properties (machine-verified in [`test_soundness_and_ch12.py`](../code/tests/test_soundness_and_ch12.py)):*
1. **Exact Cappell–Shaneson Semidirect Product Groups**: *Let $G_i = \pi_1(S^4 \setminus \mathrm{CS}_{A_i})$. Then $G_i$ is metabelian with free abelian commutator subgroup $G_i' = [G_i, G_i] \cong \mathbb{Z}^3$ and splits as*
   $$G_i \;\cong\; \mathbb{Z}^3 \rtimes_{A_i} \mathbb{Z}, \qquad i \in \{0, 1\},$$
   *where Gompf's matrices $A_0, A_1 \in \operatorname{SL}(3, \mathbb{Z})$ (Cappell–Shaneson 1976; Gompf 1991) satisfy $\det(A_i) = +1$ and $\det(I - A_i) = -1$ with characteristic polynomials*
   $$\det(tI - A_0) \;=\; t^3 - 2t^2 + t - 1 \quad \bigl(\det(tI - A_0^{-1}) = t^3 - t^2 + 2t - 1,\;\; \det(I - A_0^{-1}) = +1\bigr),$$
   $$\det(tI - A_1) \;=\; t^3 - 3t^2 + 2t - 1 \quad \bigl(\det(tI - A_1^{-1}) = t^3 - 2t^2 + 3t - 1,\;\; \det(I - A_1^{-1}) = +1\bigr).$$
2. **Strict Non-Ribbonness ($\operatorname{def}(G_i) \le 0$)**: *Neither $\mathrm{CS}_{A_0}$ nor $\mathrm{CS}_{A_1}$ is a ribbon $2$-knot (nor $\pi_1$-equivalent to any ribbon $2$-knot), and both have exact Yoshikawa index $\operatorname{ch}(\mathrm{CS}_{A_0}) = \operatorname{ch}(\mathrm{CS}_{A_1}) = 12$ (certified conditional on completeness of the $\operatorname{ch} \le 11$ census).*
3. **Topological Exterior Identification & Gluck Non-Reflexivity**: *By Hillman (1989, Ch. 6, Theorem 6.11), both $\mathrm{CS}_{A_0}$ and $\mathrm{CS}_{A_1}$ have exterior homeomorphic to the punctured-$T^3$-fibered Cappell–Shaneson exterior $CS(A_i)$, hence are topologically isotopic to one of the two Gluck fillings $K_{A_i}^0, K_{A_i}^1$ on $S^2 \times D^2$. Because $\operatorname{tr}(A_0) = 2 \ge 0$ and $\operatorname{tr}(A_1) = 3 \ge 0$, $CS(A_i)$ is **Gluck non-reflexive** (Cappell–Shaneson 1976; Gordon 1976): $K_{A_i}^0 \not\cong_{\mathrm{TOP}} K_{A_i}^1$ are two topologically non-isotopic $2$-knots sharing homeomorphic exteriors (complement-based invariants identify the exterior $CS(A_i)$ up to Gluck twist and do not select between the two fillings or certify smooth equivalence).*
4. **Regular Non-Abelian Frobenius-Cover Separation from Ribbon Alexander-Twins**:
   - *`12_10` ($\mathrm{CS}_{A_0}$) and the $\operatorname{ch} = 11$ 1-fusion ribbon $2$-knot `11_2` $= R(b^{-1}aba)$ share identical $E_1 = (x^3 - x^2 + 2x - 1)$, identical $\Sigma_2$ ($H_1 = \mathbb{Z}_5, b_2 = 0$), identical 3-fold cyclic cover $M_3$ ($H_1 = \mathbb{Z}_7$), and identical hom counts into all 10 Layer-5 groups, but are separated by their regular $21$-sheeted Frobenius $F_{21}$-cover homology $(b_1, \operatorname{Tors} H_1, b_2, b_3)$:*
     $$F_{21}\text{-cover}(\mathrm{CS}_{A_0}) = (1,\; \mathbb{Z}_7,\; 6,\; 6) \quad \neq \quad (1,\; \mathbb{Z}_2^5 \oplus \mathbb{Z}_{14},\; 6,\; 6) = F_{21}\text{-cover}(11_2).$$
   - *`12_8` ($\mathrm{CS}_{A_1}$) and its $\operatorname{ch} = 12$ 1-fusion ribbon twin `12_7` (matching $R(b^{-1}aba^{-1}ba)$ across all 10 layers) share identical $E_1 = (x^3 - 2x^2 + 3x - 1)$, identical $\Sigma_2$ ($H_1 = \mathbb{Z}_7, b_2 = 0$), and identical 3-fold cyclic cover $M_3$ ($H_1 = \mathbb{Z}_{19}$), but are separated by their regular $57$-sheeted Frobenius $F_{57}$-cover homology:*
     $$F_{57}\text{-cover}(\mathrm{CS}_{A_1}) = (1,\; \mathbb{Z}_{19},\; 18,\; 18) \quad \neq \quad (1,\; \mathbb{Z}_{11}^2 \oplus \mathbb{Z}_{209},\; 18,\; 18) = F_{57}\text{-cover}\bigl(R(b^{-1}aba^{-1}ba)\bigr).$$

*Proof.*
**(1) Group presentations and $G_i' \cong \mathbb{Z}^3$.** Both resolutions $L_\pm(\mathrm{CS}_{A_i})$ reduce to the 2-component classical unlink ($\chi = 2 - 2 + 2 = 2$, $S^2$). From the 10-sheet, 9-relator broken-surface Wirtinger presentation (`tietze_simplify` produces a 2-generator, 4-relator presentation, which GAP's Tietze routines further reduce to 2 generators and 2 relators), introducing the meridian $x = x_1$ and commutator-subgroup generators
$$u \;=\; x_2 x_1^{-1}, \qquad v \;=\; x u x^{-1}, \qquad w \;=\; x v x^{-1}$$
eliminates every Wirtinger generator into $x^{\pm 1}$ times a word in $\{u, v, w\}$. Knuth–Bendix completion under wreath-product order ($x$ largest; independently verified in GAP on both simplified and raw Wirtinger presentations) terminates in a **finite confluent rewriting system** (`32` rules from simplified, `80` rules from raw Wirtinger):
- All twelve commutator rewrites among $\{u^{\pm 1}, v^{\pm 1}, w^{\pm 1}\}$ reduce to equality ($[u, v] = [u, w] = [v, w] = 1$).
- Conjugation by $x$ acts on the normal-form basis $(u, v, w)$ by the companion matrix of $A_i^{-1}$ (while $x^{-1}$-conjugation acts by $A_i$):
  - On `CS_A0`: $x u x^{-1} = v$, $x v x^{-1} = w$, and $x w x^{-1} = u v^{-2} w$ ($\det(tI - A_0^{-1}) = t^3 - t^2 + 2t - 1$; $\det(tI - A_0) = t^3 - 2t^2 + t - 1$), with unique normal forms $x^k u^a v^b w^c$ ($k, a, b, c \in \mathbb{Z}$).
  - On `CS_A1`: $x u x^{-1} = v$, $x v x^{-1} = w$, and $x w x^{-1} = u v^{-2} w^3$ ($\det(tI - A_1^{-1}) = t^3 - 2t^2 + 3t - 1$; $\det(tI - A_1) = t^3 - 3t^2 + 2t - 1$), with unique normal forms $x^k u^a v^b w^c$.
Hence $G_i' = \langle u, v, w \rangle \cong \mathbb{Z}^3$ is free abelian of rank $3$ and $G_i \cong \mathbb{Z}^3 \rtimes_{A_i} \mathbb{Z}$ with $\det(A_i) = +1$ and $\det(I - A_i) = -1$.

**(2) Deficiency obstruction to ribbonness and exact $\operatorname{ch} = 12$.** Every ribbon $2$-knot $R \subset S^4$ admits a Wirtinger ribbon presentation with $g$ meridional generators and $g - 1$ band relators, hence has group deficiency $\operatorname{def}(\pi_1(S^4 \setminus R)) \ge 1$ (indeed $= 1$). By contrast, $G_i = \mathbb{Z}^3 \rtimes_{A_i} \mathbb{Z}$ is the fundamental group of the closed orientable aspherical mapping torus $M_i^4 = (T^3 \times [0, 1]) / ((z, 1) \sim (A_i z, 0))$, so $M_i^4 \simeq BG_i$ and
$$H_4(G_i; \mathbb{Z}) \;\cong\; H_4(M_i^4; \mathbb{Z}) \;\cong\; \mathbb{Z} \;\neq\; 0.$$
Suppose for contradiction that $G_i$ admitted a presentation $\langle g_1, \dots, g_n \mid r_1, \dots, r_m \rangle$ with deficiency $n - m \ge 1$, and let $K$ be its 2-dimensional presentation CW-complex. Because $G_i$ is torsion-free poly-$\mathbb{Z}$, the group ring $\mathbb{Z}[G_i]$ is a Noetherian Ore domain admitting a skew-field of fractions $\mathcal{D}$. In the cellular chain complex of the universal cover $\widetilde{K}$, $\partial_1\colon \mathcal{D}^n \to \mathcal{D}$ has image of right $\mathcal{D}$-dimension $1$ (since $G_i \neq 1$), so $\ker(\partial_1) \cong \mathcal{D}^{n-1}$. Because $H_1(\widetilde{K}; \mathbb{Z}) = 0$ and $\mathcal{D}$ is flat over $\mathbb{Z}[G_i]$, $\partial_2 \otimes_{\mathbb{Z}[G_i]} \mathcal{D}\colon \mathcal{D}^m \to \ker(\partial_1) \cong \mathcal{D}^{n-1}$ is surjective. Under $m \le n - 1$, surjectivity forces $m = n - 1$ and $\partial_2 \otimes \mathcal{D}$ to be an isomorphism; since $\mathbb{Z}[G_i]^m \hookrightarrow \mathcal{D}^m$ is injective, $\ker(\partial_2) = 0$, so $\pi_2(K) \cong H_2(\widetilde{K}; \mathbb{Z}) = 0$. This would make the 2-complex $K$ aspherical ($K \simeq BG_i$), contradicting $H_4(BG_i; \mathbb{Z}) \cong \mathbb{Z} \neq 0 = H_4(K; \mathbb{Z})$. Thus $\operatorname{def}(G_i) \le 0$ (achieved at $0$ by the 2-generator, 2-relator presentations), proving $\mathrm{CS}_{A_0}$ and $\mathrm{CS}_{A_1}$ are strictly **non-ribbon**. Since $\Delta_i(t)$ is irreducible over $\mathbb{Z}[t^{\pm 1}]$ and neither group appears anywhere in our $\operatorname{ch} \le 11$ census (`expanded_ch11.json`), both are prime with exact $\operatorname{ch} = 12$.

**(3) & (4) Hillman exterior theorem, Gluck non-reflexivity, and regular $F_{3p}$-cover torsion.** Hillman (1989, *2-Knots and Their Groups*, Theorem 6.11) proved that any $2$-knot in $S^4$ with group $\mathbb{Z}^3 \rtimes_A \mathbb{Z}$ ($A \in \operatorname{SL}(3, \mathbb{Z})$, $\det(I - A) = \pm 1$) has exterior homeomorphic to $CS(A)$, and Cappell–Shaneson (1976) / Gordon (1976) showed that for $\operatorname{tr}(A) \ge 0$ the two Gluck fillings on $CS(A)$ yield topologically non-isotopic $2$-knots. Note that abelian cyclic covers cannot separate $\mathrm{CS}_{A_i}$ from its ribbon twin because $\Sigma_k$ and $M_k$ depend only on the Alexander module $\mathbb{Z}[t^{\pm 1}]/(\Delta_i(t))$. Instead, $\det(tI - A_0^{-1}) \bmod 7$ has root $t \equiv 2$ of order $3$ in $\mathbb{F}_7^\times$, yielding a unique surjective quotient onto $F_{21} = \mathbb{Z}_7 \rtimes \mathbb{Z}_3$, while $\det(tI - A_1^{-1}) \bmod 19$ has root $t \equiv 11$ of order $3$ in $\mathbb{F}_{19}^\times$, yielding a unique surjective quotient onto $F_{57} = \mathbb{Z}_{19} \rtimes \mathbb{Z}_3$. Because $G_i'' = 1$ for metabelian $\mathbb{Z}^3 \rtimes_{A_i} \mathbb{Z}$, the regular $F_{3p}$-cover of $S^4 \setminus \mathrm{CS}_{A_i}$ factors cleanly with $\operatorname{Tors} H_1 = \mathbb{Z}_7$ ($i = 0$) and $\mathbb{Z}_{19}$ ($i = 1$), whereas on the 1-fusion ribbon twins $R(w)$ the non-trivial second commutator subgroup $G(R(w))'' \neq 1$ contributes extra torsion ($\mathbb{Z}_2^5 \oplus \mathbb{Z}_{14}$ on `11_2` in $F_{21}$, and $\mathbb{Z}_{11}^2 \oplus \mathbb{Z}_{209}$ on `12_7` in $F_{57}$), separating both pairs unconditionally. $\blacksquare$

> [!NOTE]
> **Comparison with Budney–Burton–Hillman (2012, `arXiv:1109.3899`):** Budney, Burton, and Hillman constructed an ideal triangulation with **$f_4 = 2$ pentachora** of the Cappell–Shaneson $2$-knot exterior $CS(A_{-2})$, whose characteristic polynomial is $\Delta(t) = t^3 + t^2 - 1$ (unit-associated reciprocal $-t^3 \Delta(t^{-1}) = t^3 - t - 1$, $\det(I - A_{-2}) = +1$). That elementary ideal $(x^3 + x^2 - 1) = (x^3 - x - 1)$ is **completely absent** from our exhaustive $\operatorname{ch} \le 12$ census already at Layer 3 (`E1`), certifying $\operatorname{ch}(\mathrm{CS}_{\mathrm{BBH}}) \ge 13$ (conditional on lower-$\operatorname{ch}$ completeness). By contrast, `12_9` carries $E_1 = (x^3 - x^2 + 1) = (x^3 - x + 1)$ (matching the 1-fusion ribbon $2$-knot `R(b^2*a)`, with $F_{21}$-cover $(1, \mathbb{Z}_2^2 \oplus \mathbb{Z}_{14}, 6, 6)$ and distinct from `R(b^3*a^-1)` whose polynomial $t^3 + t - 1$ has $|\Delta(-1)| = 3 \neq 1$): indeed, **no** Cappell–Shaneson matrix $B \in \operatorname{SL}(3, \mathbb{Z})$ can realize $E_1 = (x^3 - x^2 + 1)$, since both monic unit-associated cubics $t^3 - t^2 + 1$ and $t^3 - t + 1$ have constant term $-\det(B) = +1 \implies \det(B) = -1 \notin \operatorname{SL}(3, \mathbb{Z})$ (and $t \mapsto -t$ is not a unit substitution in $\mathbb{Z}[t^{\pm 1}]$).

---

## 3. Further Non-Ribbon $2$-Knots & Mochizuki $R_5$ Splitters at $\operatorname{ch} = 12$

### 3.1 The Non-Ribbon $\mathbb{F}_9$ Alexander-Module $2$-Knot (`12_2`, $E_1 = (3,\; x^2 + 1)$)
Two $(c, m) = (10, 2)$ $S^2$-survivors (merged into the `NEW_PACKET` cell `12_2` of size $2$ once strong Gröbner $E_1$ replaced $\pm 1$ evaluations, which previously misfiled them against $0_1$ because $|3|$ and $|(\pm 1)^2 + 1| = 2$ are coprime: $\gcd(3, 2) = 1$) carry first elementary ideal
$$E_1(G) \;=\; (3,\; x^2 + 1) \;\subset\; \mathbb{Z}[x^{\pm 1}], \qquad \mathbb{Z}[x^{\pm 1}] / E_1(G) \;\cong\; \mathbb{F}_9.$$
- **Strict Non-Ribbonness**: Because $\mathbb{Z}[x^{\pm 1}]$ is a UFD in which $3$ and $x^2 + 1$ are coprime non-units, $E_1(G) = (3, x^2 + 1)$ is **not principal**, so $\operatorname{def}(G) \le 0$ and the $2$-knot is **non-ribbon**.
- **Isomorphic Knot Group Across Both Packet Diagrams**: Simplifying both stored Wirtinger presentations yields the same meridian-preserving 2-generator, 2-relator group
  $$G \;\cong\; \bigl\langle a, y \;\big|\; a y a = y^3,\;\; (y a^{-1})^3 = 1 \bigr\rangle \qquad (a \text{ a meridian}).$$
- **Exclusion of $k$-Twist-Spins for All $2 \le k \le 7$**: If $F \cong \tau^k(K)$ were a Zeeman $k$-twist-spin, $\mu^k$ would be central in $G$, so $Q_k := G / \langle\!\langle \mu^k \rangle\!\rangle \cong G' \rtimes \mathbb{Z}_k$ would have derived subgroup $[Q_k, Q_k] \cong G'$ mapping onto the abelianized Alexander module $G'/G'' \cong \mathbb{Z}_3^2$. Exact coset enumeration gives $|Q_k| = 2, 3, 108, 5, 720, 2352$ for $k = 2, 3, 4, 5, 6, 7$:
  - For $k \in \{2, 3, 5\}$, $[Q_k, Q_k] = 1$; for $k = 6$, $[Q_6, Q_6] \cong \operatorname{SL}(2, 5)$ (perfect, proving $G$ is **non-solvable**); for $k = 7$, $[Q_7, Q_7] \cong \operatorname{SL}(2, 7)$ (perfect). None maps onto $\mathbb{Z}_3^2$.
  - For $k = 4$, centrality of $\mu^4$ in $G$ would imply $\langle\!\langle \mu^8 \rangle\!\rangle$ has index $\le 2$ in $\langle\!\langle \mu^4 \rangle\!\rangle$, forcing $|Q_8| \le 2 |Q_4| = 216$; instead $|Q_8| = 13{,}824 \gg 216$, excluding $k = 4$ as well.

### 3.2 Mochizuki $R_5$ 3-Cocycle Separation of $\tau^2(4_1)$ (`12_4`) vs $\tau^2(5_1)$ (`12_5`)
In the non-principal ideal cell $E_1 = (5, x + 1)$ ($\pi_1 \cong \mathbb{Z}_5 \rtimes_{-1} \mathbb{Z}$), both the 2-twist-spun figure-eight $\tau^2(4_1)$ and the 2-twist-spun $(2,5)$-torus knot $\tau^2(5_1)$ have **isomorphic** fundamental groups ($\mathbb{Z}_5 \rtimes_{-1} \mathbb{Z}$), identical Alexander ideals, identical finite-cover homology, identical $R_3$ colorings ($|\operatorname{Col}_{R_3}| = 5^0 \cdot 3^1$ trivial $\implies 3$), and identical $R_5$ coloring counts ($|\operatorname{Col}_{R_5}(D)| = 5^2 = 25$ on all $6$ diagrams, as $\operatorname{Col}_{R_5}(D)$ is an $\mathbb{F}_5$-affine space). Evaluating Mochizuki's $R_5$ dihedral 3-cocycle state sum (`mochizuki_3cocycle_invariant`, stored as `mochizuki_c5_mod5`) splits the $6$ diagrams cleanly into **two distinct invariant classes** purely via their Boltzmann phases in $\mathbb{Z}[u]/(u^5 - 1)$:
- **`12_4` ($\tau^2(4_1)$ cohort, `NEW_PACKET` of `5` move components)**: $\Phi_{R_5} = 5 + 10 u^2 + 10 u^3$ ($5$ monochromatic $+ 20$ non-monochromatic colorings, all weighted by $u^{\pm 2}$); and
- **`12_5` ($\tau^2(5_1)$, `EXACT_SINGLETON` of `1` diagram)**: $\Phi_{R_5} = 5 + 10 u + 10 u^4$ ($5$ monochromatic $+ 20$ non-monochromatic colorings, weighted by $u^{\pm 1}$), matching Asami–Satoh's (2005) dihedral 3-cocycle separation of 2-twist-spun 2-bridge knots with isomorphic knot groups.

---

## 4. Exact $\operatorname{ch}(F) = 12$ Spoke Promotions, Exact $\operatorname{ch}(F) = 13$ Certifications, & Refutation of the Universal Klein–Torus Cord Law

### 4.1 Exact $\operatorname{ch} = 12$ Promotions vs. Certified Exact $\operatorname{ch} = 13$ Non-Matches
Cross-matching [`expanded_ch12.json`](../code/data/expanded_ch12.json) (`promoted_exact_ch12_spokes`) against the `64` compiled Spoke instances in [`family_catalog.json`](../code/data/families/family_catalog.json) promotes **exactly `7` Spoke families** to $\operatorname{ch} = 12$ Hub realizations (`5` by explicit Yoshikawa move paths and `2` by 10-layer invariant matching, with lower bounds $\operatorname{ch}(F) \ge 12$ conditional on $\operatorname{ch} \le 11$ census completeness):
- **`MOVE_CERTIFIED` (`5` families — provable ambient isotopy via explicit Yoshikawa move paths from `compiled_ch = 12`):**
  1. **`12_1` = `Spin(5_2)`** ($\tau^0(5_2)$, Spoke 1, $S^2$, `NEW_PACKET` of `3` move components).
  2. **`12_11` = `Spin(5_1)`** ($\tau^0(5_1)$, Spoke 1, $S^2$, `EXACT_SINGLETON`).
  3. **`12_2^1` = `Tube(4_1)`** (Spoke 3, $T^2$, `NEW_PACKET` of `4` move components; separated from composite $10_1 \# 2_1^1$ by $A_4$-cover $b_2 = 30 \neq 33$ and peripheral $\operatorname{SL}(2, \mathbb{F}_3) = 24 \neq 0$).
  4. **`12_4^{0,1}` = `Mixed(L5a1)`** (Spoke 3, $S^2 \sqcup T^2$, `NEW_PACKET` of `4` move components).
  5. **`12_3^{1,1}` = `Tube(L4a1)`** (Spoke 3, $T^2 \sqcup T^2$, `NEW_PACKET` of `4` move components).
- **`INVARIANT_MATCH` (`2` families — heuristic 10-layer invariant matches tightening upper bounds down to $\operatorname{ch} = 12$):**
  6. **`12_7` = `R(b^-1*a*b*a^-1*b*a)`** (Spoke 2, $S^2$, `compiled_ch = 25` $\to$ matched to `NEW_PACKET` `12_7` of `4` move components at $\operatorname{ch} = 12$).
  7. **`12_9` = `R(b^2*a)`** (Spoke 2, $S^2$, $\mathbb{Z}HS^4$ cover, `compiled_ch = 14` $\to$ matched to `EXACT_SINGLETON` `12_9` at $\operatorname{ch} = 12$).
- **Four Newly Certified Exact $\operatorname{ch}(F) = 13$ Families (Zero $\operatorname{ch} = 12$ Matches):** Conversely, the remaining `37` of `44` previously unmatched Spoke classes (`44 - 7 = 37`) have **zero** 10-layer invariant matches anywhere in `expanded_ch12.json` (neither among the `34` new prime classes nor among the `742` residual diagrams), sharpening their lower bounds to $\operatorname{ch}(F) \ge 13$ (conditional on $\operatorname{ch} \le 12$ completeness). In particular, the **four** Spoke families compiled with upper bound `compiled_ch = 13` — **`SymUnion(5_1, 1)`** ($U_1(5_1)$, $\Delta \doteq 2x - 3$), **`SymUnion(5_2, 1)`** ($U_1(5_2)$, $\Delta \doteq 2x^2 - 4x + 1$), **`R(b^2*a^-1)`** ($\Delta \doteq x^2 + x - 1$), and **`KleinMixed(L5a1, 1)`** ($S^2 \sqcup K^2$, $H_1(\Sigma_2) = \mathbb{Z}_8$) — are completely absent at $\operatorname{ch} \le 12$, certifying their **exact Yoshikawa index $\operatorname{ch}(F) = 13$**!

### 4.2 Counterexample Refuting the Universal Klein–Torus Cord Law $\operatorname{ch}(K(w)) = \operatorname{ch}(T(w)) + 1$
For two unknotted spheres linked by a single orientable tube $T(w)$ ($S^2 \sqcup T^2$) versus a single non-orientable tube $K(w)$ ($S^2 \sqcup K^2$) along a cord word $w \in F(a, b)$ (`analyze_marked_bigon_1handles` in [`move_search.py`](../code/src/surfaceknots/move_search.py)), inserting a half-twist across a ribbon 1-handle satisfies $\operatorname{ch}(K(w)) = \operatorname{ch}(T(w)) + 1$ on six small words:
- $w \in \{b,\; bab^{-1},\; b^2,\; bab\}$ (`6 -> 7`, `9 -> 10`, `10 -> 11`, `10 -> 11`), and
- $w \in \{baba^{-1}b^{-1},\; bab^{-1}ab\}$, whose orientable $S^2 \sqcup T^2$ links `11_1^{0,1}` and `11_2^{0,1}` have $\operatorname{ch} = 11$ and whose non-orientable $S^2 \sqcup K^2$ partners appear at $\operatorname{ch} = 12$ (`12_1^{0,-2}` and `12_2^{0,-2}`, illustrating the tube-visible upper bound $\operatorname{ch}(K(w)) \le \operatorname{ch}(T(w)) + 1$ whenever a minimal diagram of $T(w)$ displays its marked-bigon tube).

However, the exhaustive $\operatorname{ch} = 11$ and $\operatorname{ch} = 12$ censuses **refute** $\operatorname{ch}(K(w)) = \operatorname{ch}(T(w)) + 1$ as a general identity! For the 5-syllable word
$$w \;=\; b a b^{-1} a^{-1} b \;\in\; F(a, b):$$
- Its non-orientable Klein-bottle link $K(b a b^{-1} a^{-1} b)$ is realized already at $\operatorname{ch} = 11$ by **`11_1^{0,-2}`** ($\Sigma_2$: $H_1 = \mathbb{Z}_6$; exact $\operatorname{ch}(K(w)) = 11$).
- Its orientable torus partner $T(b a b^{-1} a^{-1} b)$ ($S^2 \sqcup T^2$, $\Sigma_2$: $H_1 = \mathbb{Z}_6$) is **completely absent from $\operatorname{ch} \le 11$** (disjoint from all Hub entries and all `93` $\operatorname{ch} = 11$ residuals) and appears for the first time at $\operatorname{ch} = 12$ as **`12_5^{0,1}`** (`NEW_PACKET` of `2` move components), giving
  $$\operatorname{ch}\bigl(K(b a b^{-1} a^{-1} b)\bigr) \;=\; 11 \;=\; \operatorname{ch}\bigl(T(b a b^{-1} a^{-1} b)\bigr) - 1 \;<\; 12 \;=\; \operatorname{ch}\bigl(T(b a b^{-1} a^{-1} b)\bigr).$$
Thus replacing an orientable 1-handle by a non-orientable 1-handle along the same cord word $w$ can **strictly decrease** $\operatorname{ch}$ by $1$!
