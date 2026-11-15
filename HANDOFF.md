# Handoff Guide & 64-Core EPYC Compute Playbook

This guide summarizes the current state of `surfaceknots`, architectural invariants, low-level C99 scaling limits, deferred project decisions, and prioritized compute campaigns for high-core-count machines (**64 × AMD EPYC 7502 / Zen 2**).

---

## 1. Environment Bootstrap & Sanity Check (< 60s)

From the repository root:

```bash
# 1. Check CPU topology & compiler
nproc && lscpu | grep -E "Model name|Thread|Core|Socket"
gcc --version

# 2. Create venv & install package + dev dependencies via uv
uv venv code/.venv
uv pip install --python code/.venv/bin/python -e "code[dev]"

# 3. Compile _fast_enum.so (-O3 -march=native enables Zen 2 AVX2/POPCNT/BMI2) & run all 39 tests
rm -f code/src/surfaceknots/_fast_enum.so
PYTHONPATH=code/src code/.venv/bin/pytest -v
```

All **39** tests in `code/tests/` should pass cleanly.

---

## 2. Current Verified Repository State, Editorial Standards & Deferred Decisions

### Current Verified State
1. **Yoshikawa ($\operatorname{ch} \le 10$) Hub** ([`math/yoshikawa_table.md`](./math/yoshikawa_table.md), [`code/data/yoshikawa_ch10.json`](./code/data/yoshikawa_ch10.json)): All **23** weakly prime surface-links pairwise separated with zero invariant collisions across our 10-layer filtration (backed by blind $\operatorname{ch} \le 9$ and $\operatorname{ch} \le 10$ C99 recovery checks), plus live zero-divisor Fox calculus & strong integral Gröbner proofs reproducing the **2** Ishii–Nikkuni–Oshiro (2018) Table II errata (`9_1^{1,-2}` and `10_2^{0,-2}`, while confirming `10_1^{0,1}` agrees up to $y \mapsto y^{-1}$).
2. **Expanded $\operatorname{ch} = 11$ Hub** ([`math/ch11_expansion.md`](./math/ch11_expansion.md), [`code/data/expanded_ch11.json`](./code/data/expanded_ch11.json)): Exhaustive enumeration over all `3,378` prime shadows (`OEIS A113201(11)`) across $m = 1 \dots 10$, yielding **10** certified new weakly prime invariant classes (`12` diagrams across `9` `EXACT_SINGLETON` + `1` `NEW_PACKET` cells, including packet `11_1.{a,b,c}` with certified group isomorphism, the first prime $e=0$ Klein bottle in $\operatorname{ch}$-order `11_1^{-2}` with $\pi_1 \cong 2O$ of order $48$, and 3-component `11_1^{0,0,-2}`) plus `93` documented `AMBIGUOUS_SHADOW` diagram collisions.
3. **Complete $\operatorname{ch} = 12$ Hub & Parametric Spoke Families** ([`math/ch12_expansion.md`](./math/ch12_expansion.md), [`math/family_catalogs.md`](./math/family_catalogs.md), [`code/data/expanded_ch12.json`](./code/data/expanded_ch12.json), [`code/data/families/family_catalog.json`](./code/data/families/family_catalog.json)):
   - Exhaustive `15,882`-shadow $\operatorname{ch} = 12$ ($m = 1 \dots 11$) census (`909` C-tunneled raw survivors $\to$ `805` post-flype/post-filter diagrams) stratified into **34** certified new weakly prime invariant classes (`63` diagrams across `20` `EXACT_SINGLETON` + `14` `NEW_PACKET` cells) and `55` `AMBIGUOUS_SHADOW` cells (`742` diagrams).
   - Uncovers **two distinct Cappell–Shaneson 2-knots** (`12_10 = CS_A0` and `12_8 = CS_A1`, separated from their ribbon Alexander-twins by regular $F_{21}$ and $F_{57}$ Frobenius-cover homology; topological exterior determined up to the two non-isotopic Gluck fillings), non-solvable non-ribbon $\mathbb{F}_9$ 2-knots (`12_2`, $E_1 = (3, x^2+1)$), Mochizuki $R_5$ 3-cocycle phase separation of $\tau^2(4_1)$ (`12_4`) vs $\tau^2(5_1)$ (`12_5`), and the cord counterexample $\operatorname{ch}(K(bab^{-1}a^{-1}b)) = 11 < 12 = \operatorname{ch}(T(bab^{-1}a^{-1}b))$ refuting $\operatorname{ch}(K(w)) = \operatorname{ch}(T(w)) + 1$.
   - Cross-matching against the `64` compiled Spoke instances promotes **7** Spokes to $\operatorname{ch} = 12$ (`5` `MOVE_CERTIFIED`: `Spin(5_1)`, `Spin(5_2)`, `Tube(4_1)`, `Tube(L4a1)`, `Mixed(L5a1)`; `2` `INVARIANT_MATCH`: `R(b^2*a)`, `R(b^-1*a*b*a^-1*b*a)`), certifies **exact $\operatorname{ch}(F) = 13$** on **4** `compiled_ch = 13` families (`SymUnion(5_1,1)`, `SymUnion(5_2,1)`, `R(b^2*a^-1)`, `KleinMixed(L5a1,1)`), and yields **104** distinct surface-link classes (`67` Hub + `37` Spoke $\operatorname{ch} \ge 13$).
4. **Complexity Measures, $\mathbb{RP}^2$ Censuses & Braided Tri-Plane Census** ([`math/complexity_measures.md`](./math/complexity_measures.md), [`code/data/complexity/`](./code/data/complexity/)):
   - Intrinsic unknot floors $\operatorname{ch}_0(F)$, excess $\operatorname{ch}$-index $\widetilde{\operatorname{ch}}(F) = \operatorname{ch}(F) - \operatorname{ch}_0(F)$, Klein-bottle subadditivity ($\operatorname{ch}(P_+ \# P_-) = 3 < 4$, Yoshikawa Remark 2.8), Viro's smooth collapse ($10_3 \# P_\pm \cong P_\pm$), and Jabłonowski's open problem ($10_2 \# P_\pm \stackrel{?}{\cong} 8_1 \# P_\pm$, `arXiv:1505.07995`, Remark 3.15).
   - Exhaustive $\operatorname{ch} \le 12$ $\mathbb{RP}^2$ censuses confirming zero new prime $\mathbb{RP}^2$-knots up to $\operatorname{ch} \le 12$ (contextualized against Conway–Orson–Powell 2023 and Hughes–Kim–Miller–Nahm 2026, `arXiv:2605.12921`).
   - Braided $b \le 4, c_{\text{tri}} \le 6$ C99 tri-plane enumeration (`62,974,883` valid braided triples after removing the pre-audit `pruned_shared_untouched_arch_b4` shortcut), consistent with Gong–Lewis-Monkman–Osnes (2026) Conjectures 3.6 & 3.7 on the braided $b \le 4$ locus and characterizing the concentrated Meander Interleaving Obstruction.
5. **Visual Gallery, Interactive Web Atlas, Manuscript & Beamer Slides** ([`pics/`](./pics/), [`expo/`](./expo/)): IRLS $p$-Laplacian Tutte SVG + TikZ figures, pixel-registered 3-panel splitting movies $L_-(D) \to D \to L_+(D)$, research manuscript (`report.tex` / `report.pdf`), Beamer presentation (`slides.tex` / `slides.pdf` / `speaker_notes.md`), and self-contained Web Atlas (`atlas.py`).

### Rigor Rules & Dataset Synchronization Notes
- **Zero uncertified unlink heuristics**: Every candidate diagram returned by `_fast_enum.c` (or `triplane_census.py`) MUST be certified by explicit Reidemeister $\mathrm{R1}/\mathrm{R2}/\mathrm{R3}$ unlinking sequences (`certify_classical_unlink` in `unlink.py`). Never rely on polynomial necessary conditions alone to claim admissibility.
- **Explicit degrees of certainty**: Always distinguish:
  1. **`MOVE_CERTIFIED`** (provable ambient isotopy via an explicit independently verified Yoshikawa move path) vs **`INVARIANT_MATCH`** (heuristic match across all 10 filtration layers);
  2. **Conditional lower bounds** ($\operatorname{ch}(F) \ge N$ conditional on completeness of lower-$\operatorname{ch}$ tables $\operatorname{Refs}_{<N}$) vs **unconditional move-height $H$ bounds**;
  3. **Topological exterior equivalence up to Gluck twist** (e.g. `CS_A0`, `CS_A1` where $\operatorname{tr}(A_m) \ge 0$ produces two topologically non-isotopic Gluck fillings $K_{A_m}^0 \not\cong_{\mathrm{TOP}} K_{A_m}^1$) vs **smooth ambient isotopy**.
- **Phase 2 Downstream Dataset Regeneration**: Whenever `yoshikawa_data.py`, `families.py`, `triplane_census.py`, or `expanded_ch12.json` is updated, regenerate downstream artifacts in order (`reproduce_datasets.py --families --atlas`, `complexity_comparison.*`, and unpruned `triplane_b4_census.json`) so `family_catalog.json` (`hub_max_ch = 12`) and `expo/surfaceknots_atlas.html` (`104` classes) stay byte-for-byte synchronized with `expanded_ch12.json`.

---

## 3. C99 Engine Internals & 64-Core EPYC Scaling Guardrails

### Threading vs Multiprocessing Profile
- **C99 kernels (`_fast_enum.c`, `4,853` lines) are 100% thread-safe & GIL-free**:
  - There is **zero mutable global/static state** in [`code/src/surfaceknots/_fast_enum.c`](./code/src/surfaceknots/_fast_enum.c), and `ctypes.CDLL` automatically releases the Python GIL around every C call.
  - `ThreadPoolExecutor(max_workers=64)` drives all 64 Zen 2 cores at 100% native utilization during shadow generation, candidate scanning, C mountain-pass tunneling, and C tri-plane triple scanning—with zero shared-memory copy or pickle overhead (`sh_buf` is shared read-only via pointer offsets).
- **Load-balance shadow slices (`num_slices >> 64`)**:
  - Survivor counts across `OEIS A113201` prime shadows (`803` at $n=10$; `3,378` at $n=11$; `15,882` at $n=12$; `77,185` at $n=13$; `393,075` at $n=14$) are heavy-tailed. Partition `sh_buf` into **`512`–`2048` small slices** dispatched over `ThreadPoolExecutor(max_workers=64)` rather than `64` equal static chunks so straggler shadows never leave 63 cores idle.
- **Parallel chunked `tunnel_diagrams_c` + binary-tree reduction**:
  - At $\operatorname{ch} \ge 12$, run `tunnel_diagrams_c` in parallel chunks across `ThreadPoolExecutor(max_workers=64)` followed by a binary-tree reduction (`merge_and_filter_diagrams_c` and `tunnel_and_merge_diagrams_c` are completely re-entrant).
- **Pure-Python 10-layer invariant stratification IS GIL-bound — use `ProcessPoolExecutor(64)`**:
  - `compute_full_invariant_record` (integral multivariable Fox Gröbner/SNF over $\mathbb{Z}[H_1]$, Mochizuki/Kamada–Oshiro 3-cocycle state sums, $\Sigma_2$ Goeritz covers, regular cover Betti profiles, finite-group hom counts) runs in pure Python. Always map survivor invariant evaluation across `ProcessPoolExecutor(max_workers=64)` and deduplicate simplified $\pi_1$ presentations first.

### Compile-Time Constants in `code/src/surfaceknots/_fast_enum.c`
1. **`MAX_N 20` (line 25) & `uint32_t` masks in `CandidateDiagram`**:
   - Supports diagrams up to $20$ vertices out of the box (so $+2$ stabilization tunneling $\operatorname{ch} \to \operatorname{ch} + 2$ works cleanly through $\operatorname{ch} = 18$).
2. **`MAX_STATES (1 << 15)` (line 26) — Classical Crossing Limit $c \le 15$**:
   - Sizes the Fast Kauffman Bracket Butterfly arrays (`c <= 15`, safe out of the box for all $\operatorname{ch} \le 16, m \ge 1$ and $\operatorname{ch} = 17, m \ge 2$).
   - **Guardrail**: Before scanning any sector with $c = \operatorname{ch} - m \ge 16$, heap-allocate (`malloc`/`free`) the butterfly arrays inside `_fast_enum.c` (or call `threading.stack_size(64 * 1024 * 1024)` before creating threads) so worker stacks never overflow.
3. **`MAX_QV 22` (line 3308)**:
   - Quadrangulation dualizer supports up to $\operatorname{ch} = 20$ ($F = \operatorname{ch} + 2 \le 22$), and `_A113201_HINTS` in `fast_enum.py` already includes $n \le 14$ (`393,075`).
4. **`TP_HT_SIZE 1048576` (`1 << 20`, line 4481)**:
   - Open-addressed hash table for 128-bit `(e_vec, H1_SNF, Tietze_rel_hash)` tri-plane buckets, already heap-allocated via `malloc` (`max_buckets <= TP_HT_SIZE / 2 = 524,288`). Always assert `returned_buckets < max_buckets` in Python when scanning $c_{\text{tri}} \ge 7, 8$ or $b = 5$.

---

## 4. Prioritized 64-Core Compute Playbook (What to Run Next)

### Campaign 1 (~1–2 hours on 64 cores): Targeted & Complete $\operatorname{ch} = 13$ (`77,185` Shadows) & $\operatorname{ch} = 14$ (`393,075` Shadows) Sweeps
Use `enumerate_targeted_census_level_c` (which filters `m_min..m_max`, `target_components`, `target_orientability`, and $\chi(F)$ in $O(n)$ *before* Kauffman butterfly transforms) on three high-value topological sectors:
1. **Irreducible vs Reducible $\mathbb{RP}^2$ in $\mathbb{R}^4$ (`target_components=1, target_orientability=0, target_chi_min=1, target_chi_max=1`)**:
   - *Mathematical hook*: Hughes–Kim–Miller–Nahm (`arXiv:2605.12921`, May 2026) disproved Kinoshita's conjecture by constructing an irreducible $\mathbb{RP}^2 \subset S^4$ whose peripheral map has kernel $\mathbb{Z}_2$ (meridian order $> 2$). Through $\operatorname{ch} \le 12$ (all $m$), zero new prime $\mathbb{RP}^2$-knots exist (every survivor has $\pi_1 \cong \mathbb{Z}_2, S_3, D_{10}$, matching $2_1^{-1}$ or $K \# P_\pm$). Pinning down the minimum $\operatorname{ch}$-index of a prime or irreducible $\mathbb{RP}^2$ at $\operatorname{ch} \ge 13$ is a sharp open target.
2. **$m = 2$ Prime $2$-Knots at $\operatorname{ch} = 13, 14$ (`m_min=2, m_max=2, target_components=1, target_chi_min=2, target_chi_max=2`)**:
   - Directly tests `MOVE_CERTIFIED` realizations of the `4` certified exact $\operatorname{ch} = 13$ Spoke $2$-knots/links (`SymUnion(5_1,1)`, `SymUnion(5_2,1)`, `R(b^2*a^-1)`), searches for Budney–Burton–Hillman's Cappell–Shaneson $2$-knot $\mathrm{CS}_{\mathrm{BBH}}$ ($\Delta(t) = t^3 + t^2 - 1$, proven $\operatorname{ch} \ge 13$), and resolves exact $\operatorname{ch} \in \{13, 14\}$ on `Spin(6_1)`, `Spin(6_2)`, `Spin(6_3)`.
3. **Orientable Tori & Composite-Torus Splitters at $\operatorname{ch} \ge 13$**:
   - Add dihedral $D_{2p}$ ($D_{10}, D_{14}, D_{18}$) regular-cover Betti profiles and peripheral $\operatorname{SL}(2, \mathbb{F}_5)$ / $\operatorname{SL}(2, \mathbb{F}_7)$ holonomy to `tier3.py` / `linking.py` to break the 10-layer invariant collision between `Tube(K)` and `Spin(K) # 2_1^1` for $K \in \{5_1, 5_2, 6_1\}$.

### Campaign 2 (~30–90 min on 64 cores): Crushing Equivalence Packets (`11_1.{a,b,c}` & the `14` $\operatorname{ch} = 12$ `NEW_PACKET` Cells)
Attack unresolved `NEW_PACKET` and `AMBIGUOUS_SHADOW` cells from both sides simultaneously:
- **From above (move equivalence)**: Run parallel 64-thread C bidirectional BFS / Tabu-guided mountain-pass search through $\operatorname{ch} + 2$ and $\operatorname{ch} + 3$ ($\Omega_1^+, \Omega_6^+, \Omega_2^+$ + level $\Omega_3/\Omega_4/\Omega_5/\Omega_8$ + 4D Tait flypes, `max_bfs_states = 65536..262144`) on `11_1.{a,b,c}`, `12_1` (`Spin(5_2)`), `12_2` ($(3, x^2+1)$), `12_4` ($\tau^2(4_1)$), `12_6` ($\tau^2(5_2)$), and `12_7` ($R(b^{-1}aba^{-1}ba)$), verifying every discovered path via `verify_move_path_independent`.
- **From below (stronger algebraic splitters)**: Add hom-count and twisted-invariant layers over `ProcessPoolExecutor(64)` for non-solvable / linear targets ($\operatorname{PSL}(2, 7)$ of order $168$, $\operatorname{PSL}(2, 8)$, $\operatorname{SL}(2, \mathbb{F}_5)$, $S_5$) and $R_7$ dihedral 3-cocycle state sums (e.g. testing whether the two diagrams in `12_6` split into $\tau^2(5_2)$ vs another determinant-$7$ 2-twist-spin).

### Campaign 3 (~1–4 hours on 64 cores): Braided Tri-Plane Census at $(b = 4, c_{\text{tri}} = 7, 8)$ & Unconditional Non-Braided $b = 4, c_{\text{tri}} \le 6$
1. **Regenerate unpruned `triplane_b4_census.json` & extend braided $b = 4$ to $c_{\text{tri}} \in \{7, 8\}$**: Dispatch independent `(c1, c2, c3)` C99 `scan_triplane_level_c` calls concurrently across `ThreadPoolExecutor(max_workers=64)` (`TP_HT_SIZE = 1,048,576` is already heap-allocated!) to determine exact $c_{\text{tri}}^{\text{braid}, b \le 4}(9_1) \ge 7$.
2. **General (Non-Braided) Planar 8-Boundary 4-Tangle Enumeration at $c_{\text{tri}} \le 6$**: Augment braided tangle generation in `triplane_census.py` with general planar 8-boundary tangle shadow enumeration up to $6$ crossings to test GLMO (2026) Conjectures 3.6 & 3.7 unconditionally at $b = 4$ without the braided hypothesis.

### Campaign 4 (~30–60 min on 64 cores): Parallel Spoke Compiler $\operatorname{ch}_{\max}$ Tightening & `MOVE_CERTIFIED` Bridges
- Dispatch 64-core randomized Reidemeister/ribbon-word pre-simplification + simulated-annealing Yoshikawa reduction (`simplify_compiled_diagram` in `families.py`) to close the move gap on `INVARIANT_MATCH` entries (`10_1 = R(b^-1*a*b*a^-1)`, `11_1 = R(b^-1*a*b^-1*a^-1)`, `11_2 = R(b^-1*a*b*a)`, `12_7 = R(b^-1*a*b*a^-1*b*a)`, `12_9 = R(b^2*a)`), upgrading heuristic invariant matches into provable `MOVE_CERTIFIED` isotopies.

---

## 5. Ready-to-Paste Starter Prompt for a 64-Core Compute Session

```text
We are continuing research on `surfaceknots` on a 64-core AMD EPYC 7502 machine.
First read `HANDOFF.md`, `README.md`, `math/README.md`, and `code/README.md`, verify `code/.venv` compiles `_fast_enum.c` with `-O3 -march=native` cleanly, and confirm all 39 pytest tests pass.
Then inspect the 64-core scaling notes and compute campaigns in `HANDOFF.md` (Sections 3 & 4), ensure Phase 2 downstream exports (`family_catalog.json`, `triplane_b4_census.json`, `surfaceknots_atlas.html`) are synchronized with `expanded_ch12.json`, and launch Campaign 1 (targeted ch = 13 `77,185`-shadow censuses for m = 2 2-knots and RP^2-knots)!
```
