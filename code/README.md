# `surfaceknots` — Computational Enumeration & Invariants of Knotted Surfaces in $\mathbb{R}^4$

Python 3.11 + C99 package for generating marked ch-diagrams (Yoshikawa diagrams) and bridge trisection tri-plane diagrams, certifying classical unlink resolutions, simplifying under Yoshikawa moves $\Omega_1$–$\Omega_8$ and 4D Tait flypes, and computing algebraic, homological, peripheral, quandle/biquandle, and 4D linking invariants for knotted surfaces and surface-links in $\mathbb{R}^4$.

## Environment & Test Suite

Using the local `uv` virtual environment in `code/.venv`:

```bash
# Run the full pytest verification suite
PYTHONPATH=src .venv/bin/pytest -q
```

On first import of `surfaceknots.fast_enum`, the C99 enumeration and tunneling kernel [`src/surfaceknots/_fast_enum.c`](./src/surfaceknots/_fast_enum.c) is automatically compiled via `gcc -O3 -march=native -shared -fPIC` into `_fast_enum.so`.

## Package Architecture

- [`quadrangulations.py`](./src/surfaceknots/quadrangulations.py): Exact recursive generation of simple planar quadrangulations and dual prime 4-regular spherical link shadows (OEIS `A113201`: `1, 1, 2, 3, 9, 18, 62, 198, 803, 3378` for $n = 2 \dots 11$).
- [`diagrams.py`](./src/surfaceknots/diagrams.py): Combinatorial `ClassicalLinkDiagram` and `MarkedDiagram` representations, Extended Planar Diagram (EPD) parser/serializer, negative/positive Lomonaco–Yoshikawa resolutions $L_-(D)$ and $L_+(D)$, exact surface topology ($\chi$, orientability, component genera, and normal Euler numbers $e(F_i)$ via resolution-writhe difference), and $\mathrm{O}(4)$ mirror / orientation involutions.
- [`unlink.py`](./src/surfaceknots/unlink.py): Constructive Reidemeister $\mathrm{R1}/\mathrm{R2}/\mathrm{R3}$ unlink certifier producing explicit `UnlinkCertificate` move movies (zero uncertified polynomial heuristics).
- [`moves.py`](./src/surfaceknots/moves.py): Exact planar port-rewiring implementations of Yoshikawa moves $\Omega_1$–$\Omega_8$, planar 2-edge cut ($S^3 \subset S^4$ connected-sum) detection, cyclic unknotted-tube pruning, and 4D Tait flypes.
- [`_fast_enum.c`](./src/surfaceknots/_fast_enum.c) & [`fast_enum.py`](./src/surfaceknots/fast_enum.py): Fast C99 kernel combining marked-bigon pruning, linking-parity filtering, two-stage $O(c \cdot 2^c)$ Fast Kauffman Bracket Butterfly Transforms modulo $10^9 + 7$, canonical spherical BFS hashing across all $8$ $\mathrm{O}(4)$/reflection symmetries, bitmask Eulerian 2-edge-cut filtering, and multi-pass $+1/+2$ ($\Omega_1^+, \Omega_6^+, \Omega_2^+$) mountain-pass stabilization tunneling.
- [`triplane.py`](./src/surfaceknots/triplane.py): Meier–Zupan `(b; p_12, p_23, p_31)` `TriPlaneDiagram` and `TangleDiagram` constructors, Catalan matching generator, constructive pairwise unlink certification, and rigorous bridge bounds from marked diagrams.
- [`invariants/`](./src/surfaceknots/invariants/):
  - [`fundamental_group.py`](./src/surfaceknots/invariants/fundamental_group.py): Broken-surface Wirtinger presentations of $\pi_1(S^4 \setminus F)$, Tietze simplification, $H_1$ invariant factors, and $|\operatorname{Hom}(\pi_1, G)|$ counts for $G \in \{S_3, A_4, S_4, A_5, Q_8, \operatorname{Dic}_{16}\}$.
  - [`alexander.py`](./src/surfaceknots/invariants/alexander.py): Fox free calculus over integral group rings $\mathbb{Z}[H_1]$ with $\mathbb{Z}_2$-torsion variables ($u_j^2 = 1$) for non-orientable components, computing multivariable elementary ideals $E_k$ via integral Hermite/Smith normal forms without rational division, plus twisted $S_3$ Alexander ranks.
  - [`quandle.py`](./src/surfaceknots/invariants/quandle.py): Kei ($\mathbb{Z}_2[t]/(t+1)$), Alexander quandle ($\mathbb{F}_p[t^{\pm 1}]$), bikei, and Alexander biquandle invariants, plus Mochizuki's dihedral 3-cocycle state-sum invariant $\Phi_{\theta_3}(D) \in \mathbb{Z}[\mathbb{Z}_3]$ certifying 4D chirality.
  - [`linking.py`](./src/surfaceknots/invariants/linking.py): Subsurface diagram extraction & proper sublink analysis, peripheral longitudinal holonomy in $\operatorname{SL}(2, \mathbb{F}_3)$ (separating ribbon tubes $\operatorname{Tube}(K)$ from twist-spun spheres $\tau^0(K)$), Carter–Saito–Satoh mod-2 triple-linking $t_{ijk}$, and Sanderson's mod-2 double-linking invariant $\operatorname{Dlk}(F_i, F_j) \in \mathbb{Z}_2$ via cellular sky-potential transport and Rokhlin's spin quadratic form.
  - [`tier3.py`](./src/surfaceknots/invariants/tier3.py): Tier-3 4D invariants — Lomonaco second homotopy $\pi_2(S^4 \setminus F)$ $\mathbb{Z}$-ranks and regular finite-cover Betti profiles via Poincaré–Lefschetz duality, double branched cover $\Sigma_2(S^4, F)$ homology $H_1(\Sigma_2; \mathbb{Z})$, $b_2(\Sigma_2)$, and $|\sigma(\Sigma_2)|$, Kamada–Oshiro octahedral symmetric quandle $(X_6, \rho)$ colorings and non-orientable 3-cocycle state sums, and full $\mathrm{O}(4)$ chirality / orientation symmetry profiles.
- [`yoshikawa_data.py`](./src/surfaceknots/yoshikawa_data.py): Canonical diagrams, 10-layer invariant filtration, two-tier equivalence-cell classification (`EXACT_SINGLETON`, `NEW_PACKET`, `AMBIGUOUS_SHADOW`), nearest separating witnesses, and automated proof of all three Yoshikawa (1994) Table II errata (`9_1^{1,-2}`, `10_1^{0,1}`, `10_2^{0,-2}`).
- [`ch11_census.py`](./src/surfaceknots/ch11_census.py): Connected-sum / split-union / unknot reference catalog up to $\operatorname{ch} \le 12$, post-flype C mountain-pass tunneling, and exact stratification/export of the expanded $\operatorname{ch} = 11$ census.
- [`drawing.py`](./src/surfaceknots/drawing.py): Automated planar layout engine (face-area maximizing Tutte harmonic embedding with $p$-Laplacian IRLS vertex repulsion and Bézier parallel-edge routing), pixel-registered 3-panel hyperbolic splitting movie triptychs ($L_-(D) \to D \to L_+(D)$), and vector SVG + standalone TikZ/PGF exporters for `/pics`.
- [`families.py`](./src/surfaceknots/families.py): Constructive planar diagram compilers for Spoke 1 ($\tau^0(K)$ Artin spins and $U_m(K)$ twisted symmetric-union 2-knots), Spoke 2 (NORTH/SOUTH alternating planar 1-fusion ribbon words $R(w)$), and Spoke 3 ($\operatorname{Tube}(L)$ Satoh tubes, mixed sphere–torus links, and non-orientable Klein-bottle twisted links), with automated Hub cross-join and certified $[\operatorname{ch}_{\min}, \operatorname{ch}_{\max}]$ bounds.
- [`atlas.py`](./src/surfaceknots/atlas.py): Self-contained interactive KnotInfo-style Web Atlas generator (`../expo/surfaceknots_atlas.html`) uniting all **74** distinct surface-link classes ($33$ Hub + $41$ Spoke $\operatorname{ch} \ge 12$) with client-side IRLS $p$-Laplacian Tutte SVG/movie rendering, LaTeX/TikZ export, a 10-layer pairwise invariant separator, and a live custom EPD sandbox.

## Machine-Verified Datasets (`code/data/`)

- [`data/yoshikawa_ch10.json`](./data/yoshikawa_ch10.json) & [`data/yoshikawa_ch10.csv`](./data/yoshikawa_ch10.csv): Complete machine-verified annotation of all **23** weakly prime surface-links with $\operatorname{ch} \le 10$ (zero invariant collisions across all 23 entries; all 23 certified `EXACT_SINGLETON`).
- [`data/expanded_ch11.json`](./data/expanded_ch11.json) & [`data/expanded_ch11.csv`](./data/expanded_ch11.csv): Expanded $\operatorname{ch} = 11$ census stratified into **10** certified new weakly prime invariant classes ($12$ diagrams: $11_1, 11_2, 11_1^{-2}$, six 2-component links, and the 3-component non-orientable link $11_1^{0,0,-2}$) and $93$ documented `AMBIGUOUS_SHADOW` diagram collisions.
- [`data/families/family_catalog.json`](./data/families/family_catalog.json), [`data/families/family_catalog.csv`](./data/families/family_catalog.csv), & [`data/families/hub_realizations.json`](./data/families/hub_realizations.json): Parametric Spoke 1–3 surface-link catalog cross-joined against the $\operatorname{ch} \le 11$ Hub, including certified exact $\operatorname{ch}(F) = 12$ proofs.

To regenerate datasets, galleries, and the interactive Web Atlas from `code/` via [`scripts/reproduce_datasets.py`](./scripts/reproduce_datasets.py):

```bash
PYTHONPATH=src .venv/bin/python scripts/reproduce_datasets.py --yoshikawa
PYTHONPATH=src .venv/bin/python scripts/reproduce_datasets.py --ch11-enumerate --ch11-export
PYTHONPATH=src .venv/bin/python scripts/reproduce_datasets.py --gallery --extended-gallery
PYTHONPATH=src .venv/bin/python scripts/reproduce_datasets.py --families --atlas
```

## Quickstart Example

```python
from surfaceknots import MarkedDiagram, YOSHIKAWA_BY_LABEL, compute_full_invariant_record

d = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL["10_2"].canonical_epd)
rec = compute_full_invariant_record(d, label="10_2")
print(rec["e1_generators"], rec["mochizuki_c3_mod3"], rec["symmetry"]["summary"])
# ['3', 'x + 1'] ['3 + 6*u', '3 + 6*u^2'] chiral (no mirror symmetry certified)
```
