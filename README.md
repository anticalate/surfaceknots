# Knotted Surfaces in $\mathbb{R}^4$ — An Atlas & Computational Census

A mathematical and computational research repository constructing a machine-verified atlas of **knotted surfaces and surface-links in $\mathbb{R}^4$**—orientable and non-orientable, connected and multi-component, up to $\mathrm{O}(4)$ equivalence—uniting exhaustive planar `ch`-diagram and braided bridge-trisection censuses with constructive parametric families, a 10-layer 4D topological invariant engine, and automated vector rendering.

> **Live Interactive WebAtlas:** **[https://anticalate.github.io/surfaceknots/](https://anticalate.github.io/surfaceknots/)** &nbsp;·&nbsp; **[Slides (PDF)](https://anticalate.github.io/surfaceknots/slides.pdf)** &nbsp;·&nbsp; **[Report (PDF)](https://anticalate.github.io/surfaceknots/report.pdf)**

| ![8_1](./pics/svg/8_1.svg) | ![10_2](./pics/svg/10_2.svg) | ![11_1^0,0,-2](./pics/svg/11_1_sup_0_0_m2.svg) |
| :---: | :---: | :---: |
| **`8_1` — Artin Spun Trefoil** ($S^2$, $\operatorname{ch}=8$)<br>$\pi_1$-ribbon · $c_{\text{tri}}\le 6$ · amphicheiral | **`10_2` — 2-Twist-Spun Trefoil** ($S^2$, $\operatorname{ch}=10$)<br>non-ribbon · $c_{\text{tri}}=6$ · oriented-chiral & non-invertible ($\Phi_{\theta_3}$) | **`11_1^{0,0,-2}` — Klein–Sphere–Sphere Link** ($K^2 \sqcup S^2 \sqcup S^2$, $\operatorname{ch}=11$)<br>3-component non-orientable · $e = (0, 0, 0)$ |

![10_2 hyperbolic splitting movie](./pics/movies/10_2_movie.svg)
*Pixel-registered hyperbolic splitting movie $L_-(D) \to D \to L_+(D)$ for the $(-)$-amphicheiral (oriented-chiral, unoriented $\mathrm{SO}(4)$-amphicheiral) 2-twist-spun trefoil `10_2`, certifying classical unlinks at $t = \pm \varepsilon$ and $\chi(S^2) = c_- - m + c_+ = 2 - 2 + 2 = 2$.*

---

## Highlights & Current Census Scope

- **Yoshikawa's Table ($\operatorname{ch} \le 10$) — Machine Verification, Blind Recovery & Errata Proofs**:
  - Complete machine verification of all **23** weakly prime surface-links with $\operatorname{ch}(F) \le 10$ with zero invariant collisions across our 10-layer filtration, backed by blind $\operatorname{ch} \le 9$ (`294` shadows) and $\operatorname{ch} \le 10$ (`1,097` cumulative shadows) C99 recovery censuses ([`math/yoshikawa_table.md`](./math/yoshikawa_table.md), [`code/data/yoshikawa_ch10.json`](./code/data/yoshikawa_ch10.json)).
  - Meridian-pinned strong Gröbner & zero-divisor Fox calculus proofs over $\mathbb{Z}[H_1]$ reproducing both Ishii–Nikkuni–Oshiro (2018) Table II corrections (`9_1^{1,-2}` and `10_2^{0,-2}`, while confirming `10_1^{0,1}` agrees via $y \mapsto y^{-1}$).
- **Complete $\operatorname{ch} = 11$ & $\operatorname{ch} = 12$ Censuses — Including Two Cappell–Shaneson $2$-Knots**:
  - Exact prime 4-regular spherical shadow generation (`OEIS A113201`: `3,378` shadows at $\operatorname{ch}=11$; `15,882` at $\operatorname{ch}=12$; `77,185` at $\operatorname{ch}=13$; `393,075` at $\operatorname{ch}=14$) paired with a C99 Fast Kauffman Bracket Butterfly pre-filter, constructive Reidemeister $\mathrm{R1}/\mathrm{R2}/\mathrm{R3}$ unlink certification, and multi-pass $\Omega_1^+/\Omega_6^+/\Omega_2^+$ mountain-pass stabilization tunneling ([`math/ch11_expansion.md`](./math/ch11_expansion.md), [`math/ch12_expansion.md`](./math/ch12_expansion.md)).
  - **10** certified new weakly prime invariant classes at $\operatorname{ch} = 11$ (`12` diagrams across `9` `EXACT_SINGLETON` + `1` `NEW_PACKET` cells, including `11_1.{a,b,c}` with certified group isomorphism, `11_2`, the first weakly prime knotted $e=0$ Klein bottle in $\operatorname{ch}$-order `11_1^{-2}` with $\pi_1 \cong 2O$ of order $48$, and `11_1^{0,0,-2}`, plus `93` `AMBIGUOUS_SHADOW` residuals).
  - **34** certified new weakly prime invariant classes at $\operatorname{ch} = 12$ (`63` diagrams across `20` `EXACT_SINGLETON` + `14` `NEW_PACKET` cells, plus `55` `AMBIGUOUS_SHADOW` cells across `742` residuals; [`code/data/expanded_ch12.json`](./code/data/expanded_ch12.json)), uncovering **two non-ribbon Cappell–Shaneson $2$-knots** (`12_10 = CS_A0` and `12_8 = CS_A1` with $\pi_1 \cong \mathbb{Z}^3 \rtimes_{A_i} \mathbb{Z}$, separated from their 1-fusion ribbon twins by regular $F_{21}$ and $F_{57}$ Frobenius-cover homology; topological exterior determined up to the two non-isotopic Gluck fillings), non-solvable non-ribbon $\mathbb{F}_9$ $2$-knots (`12_2`, $E_1 = (3, x^2+1)$), Mochizuki $R_5$ 3-cocycle separation of $\tau^2(4_1)$ (`12_4`) vs $\tau^2(5_1)$ (`12_5`), and the cord counterexample $\operatorname{ch}(K(bab^{-1}a^{-1}b)) = 11 < 12 = \operatorname{ch}(T(bab^{-1}a^{-1}b))$ refuting $\operatorname{ch}(K(w)) = \operatorname{ch}(T(w)) + 1$.
- **Hub-and-Spoke Parametric Families (`104` Distinct Classes Across Hub + Spokes)**:
  - Constructive planar `ch`-diagram compilers (`64` instances across Spokes 1–3) for Artin spins $\tau^0(K)$, $m$-twisted symmetric-union 2-knots $U_m(K)$, general NORTH/SOUTH 1-fusion ribbon words $R(w)$, Satoh tubes $\operatorname{Tube}(L)$, mixed sphere–torus links, and non-orientable Klein-bottle links ([`math/family_catalogs.md`](./math/family_catalogs.md)).
  - Cross-joining against the **67** $\operatorname{ch} \le 12$ Hub classes (`23 + 10 + 34`) promotes **7** Spoke families to $\operatorname{ch} = 12$ (`5` `MOVE_CERTIFIED`: `Spin(5_1)`, `Spin(5_2)`, `Tube(4_1)`, `Tube(L4a1)`, `Mixed(L5a1)`; `2` `INVARIANT_MATCH`: `R(b^2*a)`, `R(b^-1*a*b*a^-1*b*a)`), certifies **exact $\operatorname{ch}(F) = 13$** on **4** unmatched `compiled_ch = 13` families (`SymUnion(5_1,1)`, `SymUnion(5_2,1)`, `R(b^2*a^-1)`, `KleinMixed(L5a1,1)`, conditional on lower-$\operatorname{ch}$ completeness), and yields **37** distinct $\operatorname{ch} \ge 13$ Spoke classes (**104** distinct surface-link classes in total).
- **Complexity Measures & $(b \le 4, c_{\text{tri}} \le 6)$ Braided Tri-Plane Census**:
  - Exact intrinsic unknot floors $\operatorname{ch}_0(F)$ and excess complexity $\widetilde{\operatorname{ch}}(F) = \operatorname{ch}(F) - \operatorname{ch}_0(F)$, strict connected-sum subadditivity ($\operatorname{ch}(P_+ \# P_-) = 3 < 4$, Yoshikawa Remark 2.8), Viro's smooth collapse $10_3 \# P_\pm \cong P_\pm$, and Jabłonowski's open problem $10_2 \# P_\pm \stackrel{?}{\cong} 8_1 \# P_\pm$ ([`math/complexity_measures.md`](./math/complexity_measures.md)).
  - C99 enumeration of **62,974,883** certified valid braided $b \le 4$ tri-plane triples up to $c_{\text{tri}} \le 6$, consistent with Gong–Lewis-Monkman–Osnes (2026) Conjectures 3.6 & 3.7 on the braided $b \le 4$ locus and characterizing the concentrated Meander Interleaving Obstruction to direct zero-crossing cap contraction.
- **Audit-Hardened 10-Layer 4D Invariant Engine & Move Verifier**:
  - Topology $(\chi, \text{orientability}, e(F_i))$ under time-first orientation $\mathbb{R}_t \times \mathbb{R}^3_{x,y,z}$ ($e = w(L_-) - w(L_+)$), Wirtinger $\pi_1(S^4 \setminus F)$ & finite-group hom counts, meridian-pinned strong Gröbner multivariable Fox ideals $E_k \subset \mathbb{Z}[H_1]$, twisted $S_3$ Alexander ranks, kei / Alexander quandle / biquandle colorings, Mochizuki $R_3, R_5$ & Kamada–Oshiro $(Q_6, \rho)$ 3-cocycle state sums, peripheral $\operatorname{SL}(2, \mathbb{F}_3)$ holonomy, Carter–Saito–Satoh triple-linking $t_{ijk}$, Rohlin–Matsumoto horizontal Seifert membrane double-linking $\operatorname{Dlk} \in \mathbb{Z}_2$, regular finite-cover Betti profiles ($\mathbb{Z}_2, \mathbb{Z}_3, S_3, Q_8, A_4, \operatorname{Dic}_{16}, F_{21}, \operatorname{SL}(2, \mathbb{F}_3), F_{57}$), double branched covers $\Sigma_2(S^4, F)$, and an independent port-labeled Yoshikawa move-path verifier (`move_search.py`).

---

## Repository Map

- **[`math/`](./math/)** — Self-contained theory, proofs, census tables, and roadmap ([`math/README.md`](./math/README.md)):
  - [`literature_review.md`](./math/literature_review.md) · [`foundations_and_invariants.md`](./math/foundations_and_invariants.md) · [`yoshikawa_table.md`](./math/yoshikawa_table.md) · [`ch11_expansion.md`](./math/ch11_expansion.md) · [`ch12_expansion.md`](./math/ch12_expansion.md) · [`family_catalogs.md`](./math/family_catalogs.md) · [`complexity_measures.md`](./math/complexity_measures.md) · [`roadmap.md`](./math/roadmap.md)
- **[`code/`](./code/)** — Python 3.11 + C99 `surfaceknots` package, test suite, and JSON/CSV datasets ([`code/README.md`](./code/README.md)):
  - Core & C99 kernels: [`diagrams.py`](./code/src/surfaceknots/diagrams.py), [`moves.py`](./code/src/surfaceknots/moves.py), [`move_search.py`](./code/src/surfaceknots/move_search.py), [`unlink.py`](./code/src/surfaceknots/unlink.py), [`_fast_enum.c`](./code/src/surfaceknots/_fast_enum.c), [`fast_enum.py`](./code/src/surfaceknots/fast_enum.py)
  - Censuses, families, complexity & tri-planes: [`ch11_census.py`](./code/src/surfaceknots/ch11_census.py), [`ch12_census.py`](./code/src/surfaceknots/ch12_census.py), [`targeted_census.py`](./code/src/surfaceknots/targeted_census.py), [`families.py`](./code/src/surfaceknots/families.py), [`complexity.py`](./code/src/surfaceknots/complexity.py), [`triplane.py`](./code/src/surfaceknots/triplane.py), [`triplane_census.py`](./code/src/surfaceknots/triplane_census.py)
  - 4D invariants: [`invariants/`](./code/src/surfaceknots/invariants/) (`fundamental_group.py`, `alexander.py`, `quandle.py`, `linking.py`, `tier3.py`)
  - Verified datasets: [`code/data/`](./code/data/) (`yoshikawa_ch10.*`, `expanded_ch11.*`, `expanded_ch12.*`, `families/*`, `complexity/*`)
- **[`pics/`](./pics/)** — Automated IRLS $p$-Laplacian Tutte SVG & standalone TikZ diagrams and hyperbolic splitting movies ([`pics/README.md`](./pics/README.md)).
- **[`expo/`](./expo/)** — Self-contained interactive KnotInfo-style Web Atlas ([`expo/surfaceknots_atlas.html`](./expo/surfaceknots_atlas.html)), research manuscript ([`expo/report.tex`](./expo/report.tex), [`expo/report.pdf`](./expo/report.pdf)), and Beamer presentation ([`expo/slides.tex`](./expo/slides.tex), [`expo/slides.pdf`](./expo/slides.pdf), [`expo/speaker_notes.md`](./expo/speaker_notes.md)).
- **[`REFS.md`](./REFS.md)** — Annotated bibliography of censuses, move calculi, 4D invariants, bridge trisections, and open problems.
- **[`HANDOFF.md`](./HANDOFF.md)** — Engineering notes, C99 scaling guardrails, and prioritized compute playbook for multi-core high-performance enumeration.

---

## Quickstart

Requires Python $\ge 3.11$ and `gcc` (`_fast_enum.c` auto-compiles to `_fast_enum.so` via `gcc -O3 -march=native` on first import):

```bash
# Bootstrap virtual environment with uv and run the full 39-test verification suite
uv venv code/.venv
uv pip install --python code/.venv/bin/python -e "code[dev]"
PYTHONPATH=code/src code/.venv/bin/pytest -v
```

Inspect a surface-link and compute its full 10-layer invariant profile in Python:

```python
from surfaceknots import MarkedDiagram, YOSHIKAWA_BY_LABEL, compute_full_invariant_record

d = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_2"].canonical_epd)
rec = compute_full_invariant_record(d, label="10_2")
print(rec["e1_generators"], rec["mochizuki_c3_mod3"], rec["symmetry"]["summary"])
# ['3', 'x + 1'] ['3 + 6*u', '3 + 6*u^2'] Oriented-chiral & non-invertible (Mochizuki 3-cocycle)
```
