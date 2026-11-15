# `surfaceknots` — Computational Enumeration & Invariants of Knotted Surfaces in $\mathbb{R}^4$

Python 3.11 + C99 package for generating marked ch-diagrams (Yoshikawa diagrams) and bridge trisection tri-plane diagrams, certifying classical unlink resolutions, simplifying under Yoshikawa moves $\Omega_1$–$\Omega_8$ and 4D Tait flypes, and computing algebraic, homological, peripheral, quandle/biquandle, and 4D linking invariants for knotted surfaces and surface-links in $\mathbb{R}^4$.

## Environment & Test Suite

Using the local `uv` virtual environment in `code/.venv`:

```bash
# Run the full 39-test pytest verification suite
PYTHONPATH=src .venv/bin/pytest -q
```

On first import of `surfaceknots.fast_enum`, the C99 enumeration and tunneling kernel [`src/surfaceknots/_fast_enum.c`](./src/surfaceknots/_fast_enum.c) (`4,853` lines) is automatically compiled via `gcc -O3 -march=native -shared -fPIC` into `_fast_enum.so`.

## Package Architecture

- [`quadrangulations.py`](./src/surfaceknots/quadrangulations.py): Exact recursive generation of simple planar quadrangulations and dual prime 4-regular spherical link shadows (OEIS `A113201`: `1, 1, 2, 3, 9, 18, 62, 198, 803, 3378, 15882, 77185` for $n = 2 \dots 13$).
- [`diagrams.py`](./src/surfaceknots/diagrams.py): Combinatorial `ClassicalLinkDiagram` and `MarkedDiagram` representations, Extended Planar Diagram (EPD) parser/serializer, negative/positive Lomonaco–Yoshikawa resolutions $L_-(D)$ and $L_+(D)$, exact surface topology ($\chi$, orientability, component genera, and normal Euler numbers $e(F_i) = w(L_-(D_i)) - w(L_+(D_i))$ under time-first orientation $\mathbb{R}_t \times \mathbb{R}^3_{x,y,z}$), and $\mathrm{O}(4)$ mirror / orientation involutions.
- [`unlink.py`](./src/surfaceknots/unlink.py): Constructive Reidemeister $\mathrm{R1}/\mathrm{R2}/\mathrm{R3}$ unlink certifier producing explicit `UnlinkCertificate` move movies (zero uncertified polynomial heuristics).
- [`moves.py`](./src/surfaceknots/moves.py): Exact planar port-rewiring implementations of Yoshikawa moves $\Omega_1$–$\Omega_8$, planar 2-edge cut ($S^3 \subset S^4$ connected-sum) detection, cyclic unknotted-tube pruning, and 4D Tait flypes.
- [`_fast_enum.c`](./src/surfaceknots/_fast_enum.c) & [`fast_enum.py`](./src/surfaceknots/fast_enum.py): Thread-safe, GIL-free C99 kernel combining marked-bigon pruning, linking-parity filtering, two-stage $O(c \cdot 2^c)$ Fast Kauffman Bracket Butterfly Transforms modulo $10^9 + 7$, canonical spherical BFS hashing across all $8$ $\mathrm{O}(4)$/reflection symmetries, bitmask Eulerian 2-edge-cut filtering, and multi-pass $+1/+2$ ($\Omega_1^+, \Omega_6^+, \Omega_2^+$) mountain-pass stabilization tunneling.
- [`triplane.py`](./src/surfaceknots/triplane.py): Meier–Zupan `(b; p_12, p_23, p_31)` `TriPlaneDiagram` and `TangleDiagram` constructors, Catalan matching generator, constructive pairwise unlink certification, tri-plane surface topology ($\chi$, orientability, normal Euler numbers), Wirtinger $\pi_1(S^4 \setminus F)$ presentation, Goeritz double branched cover $\Sigma_2(S^4, F)$ homology, concentrated $(0, 0, c_3) \to$ `MarkedDiagram` compiler with meander-interleaving detection, and rigorous bridge bounds.
- [`triplane_census.py`](./src/surfaceknots/triplane_census.py): Braided $(b \le 4, c_{\text{tri}} \le 6)$ tri-plane C99 enumeration (`62,974,883` valid braided triples) with C99 tangle-pair unlink prefiltering, concentrated $(0, 0, c_3)$ algebraic join, evaluation of Gong–Lewis-Monkman–Osnes (2026) Conjectures 3.6 & 3.7 on the braided $b \le 4$ locus, and Hub/Spoke $\mathrm{O}(4)$ cross-matching.
- [`complexity.py`](./src/surfaceknots/complexity.py): Exact intrinsic unknot floor $\operatorname{ch}_0(F) = c_0(F) + m_0(F)$, excess $\operatorname{ch}$-index $\widetilde{\operatorname{ch}}(F) = \operatorname{ch}(F) - \operatorname{ch}_0(F)$, weighted $(c + \alpha m)$ complexity profiles, constructive Klein-bottle subadditivity ($\operatorname{ch}(P_+ \# P_-) = 3 < 4$), Viro's $10_3 \# P_\pm \cong P_\pm$ collapse vs Jabłonowski's open problem diagnostics, and multi-measure comparison dataset export.
- [`move_search.py`](./src/surfaceknots/move_search.py): Independent port-labeled Yoshikawa move-path verifier (`verify_move_path_independent`), 1-handle cord double-coset classifier (`analyze_marked_bigon_1handles` over $F_2$, $B_3$, $\operatorname{BS}(1,2)$, and $\pi_1(4_1)$), and parallel bidirectional mountain-pass move search.
- [`invariants/`](./src/surfaceknots/invariants/):
  - [`fundamental_group.py`](./src/surfaceknots/invariants/fundamental_group.py): Broken-surface Wirtinger presentations of $\pi_1(S^4 \setminus F)$, Tietze simplification, $H_1$ invariant factors, and $|\operatorname{Hom}(\pi_1, G)|$ counts (with cut-vertex block decomposition and greedy relator-completing search) for $G \in \{S_3, D_8, Q_8, D_{10}, A_4, \operatorname{Dic}_{12}, \operatorname{Dic}_{16}, F_{21}, S_4, A_5\}$.
  - [`alexander.py`](./src/surfaceknots/invariants/alexander.py): Fox free calculus over meridian-pinned integral group rings $\mathbb{Z}[H_1]$ with $\mathbb{Z}_2$-torsion variables ($u_j^2 = 1$) for non-orientable components, computing canonical strong Gröbner bases and Smith normal forms over $\mathbb{Z}$ without rational division, plus twisted $S_3$ Alexander ranks.
  - [`quandle.py`](./src/surfaceknots/invariants/quandle.py): Kei ($\mathbb{Z}_2[t]/(t+1)$), Alexander quandle ($\mathbb{F}_p[t^{\pm 1}]$), bikei, and Alexander biquandle invariants, plus Mochizuki's dihedral $R_3, R_5$ 3-cocycle state-sum invariants certifying oriented 4D chirality (`10_2` is oriented-chiral & non-invertible, while unoriented $\mathrm{SO}(4)$-amphicheiral) and separating $\tau^2(4_1)$ from $\tau^2(5_1)$ via Boltzmann phases.
  - [`linking.py`](./src/surfaceknots/invariants/linking.py): Subsurface diagram extraction & proper sublink analysis, peripheral longitudinal holonomy in $\operatorname{SL}(2, \mathbb{F}_3)$ (separating ribbon tubes $\operatorname{Tube}(K)$ from twist-spun spheres $\tau^0(K)$ for $K \in \{3_1, 4_1\}$), Carter–Saito–Satoh mod-2 triple-linking $t_{ijk}$, and mod-2 double-linking $\operatorname{Dlk}(F_i, F_j) \in \mathbb{Z}_2$ via Rohlin–Matsumoto's horizontal Seifert membrane quadratic pairing (`compute_torus_dlk_pair`).
  - [`tier3.py`](./src/surfaceknots/invariants/tier3.py): Tier-3 4D invariants — Lomonaco second homotopy $\pi_2(S^4 \setminus F)$ $\mathbb{Z}$-ranks and regular finite-cover Betti profiles ($\mathbb{Z}_2, \mathbb{Z}_3, S_3, Q_8, A_4, \operatorname{Dic}_{16}, F_{21}, \operatorname{SL}(2, \mathbb{F}_3), F_{57}$, separating both $\operatorname{ch}=12$ Cappell–Shaneson 2-knots from their ribbon Alexander-twins) via Poincaré–Lefschetz duality, double branched cover $\Sigma_2(S^4, F)$ homology $H_1(\Sigma_2; \mathbb{Z})$, $b_2(\Sigma_2)$, and $|\sigma(\Sigma_2)|$, Kamada–Oshiro octahedral symmetric quandle $(Q_6, \rho)$ colorings and non-orientable 3-cocycle state sums, and $\mathrm{O}(4)$ chirality / orientation symmetry profiles.
- [`yoshikawa_data.py`](./src/surfaceknots/yoshikawa_data.py): Canonical diagrams, 10-layer invariant filtration, two-tier equivalence-cell classification (`EXACT_SINGLETON`, `NEW_PACKET`, `AMBIGUOUS_SHADOW`), nearest separating witnesses, and automated proof reproducing the two Ishii–Nikkuni–Oshiro (2018) Table II errata (`9_1^{1,-2}` and `10_2^{0,-2}`).
- [`ch11_census.py`](./src/surfaceknots/ch11_census.py), [`ch12_census.py`](./src/surfaceknots/ch12_census.py), & [`targeted_census.py`](./src/surfaceknots/targeted_census.py): Connected-sum / split-union / unknot reference catalogs, Yoshikawa Definition 5.1 weak-primality verification, post-flype C mountain-pass tunneling, and complete stratification/export of the $\operatorname{ch} = 11$, $\operatorname{ch} = 12$, and targeted $\operatorname{ch} \in \{11, 12\}$ censuses.
- [`drawing.py`](./src/surfaceknots/drawing.py): Automated planar layout engine (face-area maximizing Tutte harmonic embedding with $p$-Laplacian IRLS vertex repulsion and Bézier parallel-edge routing), pixel-registered 3-panel hyperbolic splitting movie triptychs ($L_-(D) \to D \to L_+(D)$), and vector SVG + standalone TikZ/PGF exporters for `/pics`.
- [`families.py`](./src/surfaceknots/families.py): Constructive planar diagram compilers (`64` instances) for Spoke 1 ($\tau^0(K)$ Artin spins and $U_m(K)$ twisted symmetric-union 2-knots), Spoke 2 (NORTH/SOUTH alternating planar 1-fusion ribbon words $R(w)$), and Spoke 3 ($\operatorname{Tube}(L)$ Satoh tubes, mixed sphere–torus links, and non-orientable Klein-bottle twisted links), with automated Hub cross-join and certified $[\operatorname{ch}_{\min}, \operatorname{ch}_{\max}]$ bounds.
- [`atlas.py`](./src/surfaceknots/atlas.py): Self-contained interactive KnotInfo-style Web Atlas generator (`../expo/surfaceknots_atlas.html`) uniting all **104** distinct surface-link classes across the $\operatorname{ch} \le 12$ Hub (`67` classes: `23 + 10 + 34`) and $\operatorname{ch} \ge 13$ Spokes (`37` classes) with client-side IRLS $p$-Laplacian Tutte SVG/movie rendering, LaTeX/TikZ export, a 10-layer pairwise invariant separator, and a live custom EPD sandbox.

## Machine-Verified Datasets (`code/data/`)

- [`data/yoshikawa_ch10.json`](./data/yoshikawa_ch10.json) & [`data/yoshikawa_ch10.csv`](./data/yoshikawa_ch10.csv): Complete machine-verified annotation of all **23** weakly prime surface-links with $\operatorname{ch} \le 10$ (zero invariant collisions across all 23 entries; all 23 `EXACT_SINGLETON`).
- [`data/expanded_ch11.json`](./data/expanded_ch11.json) & [`data/expanded_ch11.csv`](./data/expanded_ch11.csv): Expanded $\operatorname{ch} = 11$ census stratified into **10** certified new weakly prime invariant classes (`12` diagrams across `9` `EXACT_SINGLETON` + `1` `NEW_PACKET` cells) and `93` documented `AMBIGUOUS_SHADOW` diagram collisions.
- [`data/expanded_ch12.json`](./data/expanded_ch12.json) & [`data/expanded_ch12.csv`](./data/expanded_ch12.csv): Complete $\operatorname{ch} = 12$ ($m = 1 \dots 11$, `15,882` shadows) census stratified into **34** certified new weakly prime invariant classes (`63` diagrams across `20` `EXACT_SINGLETON` + `14` `NEW_PACKET` cells, including `CS_A0`, `CS_A1`, `(3, x^2+1)`, and $\tau^2(4_1)$ vs $\tau^2(5_1)$), `55` `AMBIGUOUS_SHADOW` cells (`742` diagrams), and `7` promoted exact $\operatorname{ch} = 12$ Spoke realizations.
- [`data/families/family_catalog.json`](./data/families/family_catalog.json), [`data/families/family_catalog.csv`](./data/families/family_catalog.csv), & [`data/families/hub_realizations.json`](./data/families/hub_realizations.json): Parametric Spoke 1–3 surface-link catalog (`64` compiled instances; cross-joined with `expanded_ch12.json` to yield `27` Hub-matched instances across `22` Hub entries, `37` $\operatorname{ch} \ge 13$ Spoke classes, and `4` certified exact $\operatorname{ch}(F) = 13$ classes).
- [`data/complexity/triplane_b4_census.json`](./data/complexity/triplane_b4_census.json), [`data/complexity/targeted_ch12_13_census.json`](./data/complexity/targeted_ch12_13_census.json), [`data/complexity/complexity_comparison.json`](./data/complexity/complexity_comparison.json), & [`data/complexity/complexity_comparison.csv`](./data/complexity/complexity_comparison.csv): Braided $(b \le 4, c_{\text{tri}} \le 6)$ tri-plane census, targeted $\operatorname{ch} \in \{11, 12\}$ C99 enumerations, and multi-measure complexity comparison tables.

To regenerate datasets, galleries, and the interactive Web Atlas from `code/` via [`scripts/reproduce_datasets.py`](./scripts/reproduce_datasets.py):

```bash
PYTHONPATH=src .venv/bin/python scripts/reproduce_datasets.py --yoshikawa
PYTHONPATH=src .venv/bin/python scripts/reproduce_datasets.py --ch11-enumerate --ch11-export
PYTHONPATH=src .venv/bin/python scripts/reproduce_datasets.py --ch12-export
PYTHONPATH=src .venv/bin/python scripts/reproduce_datasets.py --gallery --extended-gallery
PYTHONPATH=src .venv/bin/python scripts/reproduce_datasets.py --families --atlas
```

## Quickstart Example

```python
from surfaceknots import MarkedDiagram, YOSHIKAWA_BY_LABEL, compute_full_invariant_record

d = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_2"].canonical_epd)
rec = compute_full_invariant_record(d, label="10_2")
print(rec["e1_generators"], rec["mochizuki_c3_mod3"], rec["symmetry"]["summary"])
# ['3', 'x + 1'] ['3 + 6*u', '3 + 6*u^2'] Oriented-chiral & non-invertible (Mochizuki 3-cocycle)
```
