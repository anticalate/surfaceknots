# Handoff Guide & 64-Core EPYC Compute Playbook

This guide summarizes the current state of `surfaceknots`, architectural invariants, low-level C99 scaling limits, deferred project decisions, and prioritized compute campaigns for picking up work on a high-core-count machine (**64 × AMD EPYC 7502 / Zen 2**).

---

## 1. Environment Bootstrap & Sanity Check (< 60s)

All prior work—including the hyperbolic splitting movie / Web Atlas branch (`exp/atlas-and-movies`) and the complexity / C99 targeted $\operatorname{ch}=12$ / $b \le 4$ tri-plane census branch (`exp/complexity-and-cutoffs`)—is **already merged into `main`**.

From the repository root (`/home/antigravity/anticalate/surfaceknots`):

```bash
# 1. Check CPU topology & compiler
nproc && lscpu | grep -E "Model name|Thread|Core|Socket"
gcc --version

# 2. Create venv & install package + dev dependencies via uv
uv venv code/.venv
uv pip install --python code/.venv/bin/python -e "code[dev]"

# 3. Compile _fast_enum.so (-O3 -march=native enables Zen 2 AVX2/POPCNT/BMI2) & run all 29 tests
rm -f code/src/surfaceknots/_fast_enum.so
PYTHONPATH=code/src code/.venv/bin/pytest -v
```

All **29** tests in `code/tests/` should pass in ~10 seconds.

---

## 2. Current State, Editorial Standards & Deferred Decisions

### What Is Already Verified & Committed on `main`
1. **Yoshikawa ($\operatorname{ch} \le 10$) Hub** ([`math/yoshikawa_table.md`](./math/yoshikawa_table.md), [`code/data/yoshikawa_ch10.json`](./code/data/yoshikawa_ch10.json)): All **23** weakly prime surface-links certified `EXACT_SINGLETON` with zero invariant collisions across our 10-layer filtration, plus live zero-divisor Fox calculus proofs of all **3** Yoshikawa Table II errata (`9_1^{1,-2}`, `10_1^{0,1}`, `10_2^{0,-2}`).
2. **Expanded $\operatorname{ch} = 11$ Hub** ([`math/ch11_expansion.md`](./math/ch11_expansion.md), [`code/data/expanded_ch11.json`](./code/data/expanded_ch11.json)): Exhaustive enumeration over all `3,378` prime shadows (`OEIS A113201(11)`), yielding **10** certified new prime invariant classes (`12` diagrams, including packet `11_1.{a,b,c}`, the first prime $e=0$ Klein bottle `11_1^{-2}`, and 3-component `11_1^{0,0,-2}`) plus `93` documented `AMBIGUOUS_SHADOW` diagram collisions.
3. **Parametric Spoke Families** ([`math/family_catalogs.md`](./math/family_catalogs.md), [`code/data/families/family_catalog.json`](./code/data/families/family_catalog.json)): Constructive planar diagram compilers (`families.py`) for Spoke 1 ($\tau^0(K)$ Artin spins, $U_m(K)$ twisted symmetric unions), Spoke 2 ($R(w)$ NORTH/SOUTH 1-fusion ribbon words), and Spoke 3 ($\operatorname{Tube}(L)$ Satoh tubes, mixed sphere–torus links, Klein-bottle twisted links), cross-joined against the Hub with certified $[\operatorname{ch}_{\min}, \operatorname{ch}_{\max}]$ bounds.
4. **Complexity Measures, Targeted $\operatorname{ch} = 12$ & $(b \le 4, c_{\text{tri}} \le 6)$ Tri-Plane Census** ([`math/complexity_measures.md`](./math/complexity_measures.md), [`code/data/complexity/`](./code/data/complexity/)):
   - Intrinsic unknot floors $\operatorname{ch}_0(F)$, excess $\operatorname{ch}$-index $\widetilde{\operatorname{ch}}(F) = \operatorname{ch}(F) - \operatorname{ch}_0(F)$, Klein-bottle subadditivity ($\operatorname{ch}(P_+ \# P_-) = 3 < 4$), and Viro collapse ($\tau^{k+2}(K) \# P_\pm \cong \tau^k(K) \# P_\pm$).
   - Targeted C99 $\operatorname{ch} = 12$ scans (`targeted_ch12_13_census.json`): **7** certified new $m=2$ prime $2$-knots (`12_m2_1`–`12_m2_7`), promoting `R(b^2*a)`, `R(b^2*a^-1)`, and `SymUnion(5_2,1)` to exact $\operatorname{ch}(F) = 12$, plus low-saddle ($m \le 3$) non-orientable $\mathbb{RP}^2$ censuses.
   - Exhaustive C99 scan of **36,166,398** certified valid braided $b = 4$ tri-plane triples up to $c_{\text{tri}} \le 6$ (`triplane_b4_census.json`), confirming Gong–Lewis-Monkman–Osnes (2026) Conjectures 3.6 & 3.7 on the braided $b \le 4$ locus and uncovering the **Meander Interleaving Obstruction** on `10_2` (`90/90` concentrated $(0,0,6)$ realizations interleaved).
5. **Visual Gallery & Interactive Web Atlas** ([`pics/`](./pics/), [`expo/surfaceknots_atlas.html`](./expo/surfaceknots_atlas.html)): IRLS $p$-Laplacian Tutte SVG + TikZ figures, pixel-registered 3-panel splitting movies $L_-(D) \to D \to L_+(D)$, and a self-contained Web Atlas (`atlas.py`).

### Rigor Rules & Deferred Design Decisions (Important Context)
- **Zero uncertified unlink heuristics**: Every candidate diagram returned by `_fast_enum.c` (or `triplane_census.py`) MUST be certified by explicit Reidemeister $\mathrm{R1}/\mathrm{R2}/\mathrm{R3}$ unlinking sequences (`certify_classical_unlink` in `unlink.py`). Never rely on polynomial necessary conditions alone to claim admissibility.
- **Honest equivalence tiers (`EXACT_SINGLETON` vs `NEW_PACKET` vs `AMBIGUOUS_SHADOW`)**: Never silently collapse diagrams that share invariants unless connected by certified Yoshikawa moves + 4D Tait flypes + $\mathrm{O}(4)$ symmetries, and never claim two diagrams represent distinct surfaces without a separating invariant witness.
- **Deferred systematic relabeling (do in ONE unified pass later)**: Provisional labels (`11_1`, `11_1.{a,b,c}`, `12_m2_1..7`, `Spin(...)`, `R(...)`, `Tube(...)`) are intentionally left as provisional/construction tags while enumeration continues. **Do not piecemeal-rename labels mid-census**—wait until full $\operatorname{ch} = 12$ (and targeted $\operatorname{ch} = 13, 14$) tables stabilize, then apply a single systematic multi-graded naming scheme across the entire atlas.
- **Web Atlas & Gallery synchronization**: `expo/surfaceknots_atlas.html` (`atlas.py`) currently renders the `74` Hub ($\operatorname{ch} \le 11$) + Spoke ($\operatorname{ch} \ge 12$) classes, whereas `complexity_comparison.json` already tracks `78` distinct classes (adding the non-Spoke `12_m2_*` 2-knots and promoting `3` Spokes to exact $\operatorname{ch}=12$). Once full $\operatorname{ch} = 12$ finishes, fold the complete $\operatorname{ch} = 12$ Hub + updated Spoke bounds into `atlas.py` and `pics/` in one clean update.

---

## 3. C99 Engine Internals & 64-Core EPYC Scaling Guardrails

### Threading vs Multiprocessing Profile
- **C99 kernels (`_fast_enum.c`, `4,391` lines) are 100% thread-safe & GIL-free**:
  - There is **zero mutable global/static state** anywhere in [`code/src/surfaceknots/_fast_enum.c`](./code/src/surfaceknots/_fast_enum.c) (the sole `static` array is `static const int slot_refl[4]`), and `ctypes.CDLL` automatically releases the Python GIL around every C call.
  - Consequently, `ThreadPoolExecutor(max_workers=64)` drives all 64 Zen 2 cores at 100% native utilization during shadow generation, candidate scanning, C mountain-pass tunneling, and C tri-plane triple scanning—with zero shared-memory copy or pickle overhead (`sh_buf` is shared read-only via pointer offsets).
- **Load-balance shadow slices (`num_slices >> 64`)**:
  - Survivor counts across `OEIS A113201` prime shadows (`803` at $n=10$; `3,378` at $n=11$; `15,882` at $n=12$; `77,185` at $n=13$; `393,075` at $n=14$) are heavy-tailed. In `enumerate_census_level_candidates` and `enumerate_targeted_census_level_candidates` ([`fast_enum.py`](./code/src/surfaceknots/fast_enum.py)), partition `sh_buf` into **`512`–`2048` small slices** dispatched over `ThreadPoolExecutor(max_workers=64)` rather than `64` equal static chunks so straggler shadows never leave 63 cores idle.
- **Upgrade the post-scan `merge_diagrams_c` / `cascade_tunnel_diagrams_c` step from 1 thread to a 64-core tree reduction**:
  - Currently in `fast_enum.py`, worker shards scan in parallel, but the final `merge_diagrams_c` and `cascade_tunnel_diagrams_c` calls run on a **single thread** over `shard_diags`.
  - At $\operatorname{ch} \ge 12$, run `tunnel_diagrams_c` in parallel chunks across `ThreadPoolExecutor(max_workers=64)` followed by a binary-tree reduction (see `_run_parallel_chunks` in [`code/scripts/reproduce_datasets.py`](./code/scripts/reproduce_datasets.py)). Both `merge_and_filter_diagrams_c` and `tunnel_and_merge_diagrams_c` are completely re-entrant!
- **Pure-Python 10-layer invariant stratification IS GIL-bound — use `ProcessPoolExecutor(64)`**:
  - `compute_full_invariant_record` (integral multivariable Fox SNF over $\mathbb{Z}[H_1]$, Mochizuki/Kamada–Oshiro 3-cocycle state sums, $\Sigma_2$ Goeritz covers, finite-group hom counts) runs in pure Python. Always map survivor invariant evaluation across `ProcessPoolExecutor(max_workers=64)` and deduplicate simplified $\pi_1$ presentations first (`_TIETZE_CACHE` pattern in `triplane_census.py`).

### Compile-Time Constants in `code/src/surfaceknots/_fast_enum.c` (Check Before Scaling Up!)
1. **`MAX_N 20` (line 22) & `uint32_t` masks in `CandidateDiagram` (lines 25–29)**:
   - Supports diagrams up to $20$ vertices out of the box (so $+2$ stabilization tunneling $\operatorname{ch} \to \operatorname{ch} + 2$ works cleanly through $\operatorname{ch} = 18$).
2. **`MAX_STATES (1 << 14)` (line 23) — Classical Crossing Limit $c \le 14$**:
   - Sizes the four stack-allocated Fast Kauffman Bracket Butterfly arrays (`int64_t f_p3[MAX_STATES]`, etc., ~`524 KB` per thread vs Linux's `8 MB` default `pthread` stack).
   - Safe out of the box for all $\operatorname{ch} \le 15, m \ge 1$ and $\operatorname{ch} = 16, m \ge 2$.
   - **Guardrail**: Before scanning any sector with $c = \operatorname{ch} - m \ge 15$, heap-allocate (`malloc`/`free`) the butterfly arrays inside `_fast_enum.c` (or call `threading.stack_size(64 * 1024 * 1024)` before creating threads) so worker stacks never overflow!
3. **`MAX_QV 22` & `code[112]` in `QGraph` (lines 2881–2890)**:
   - Quadrangulation dualizer supports up to $\operatorname{ch} = 20$ ($F = \operatorname{ch} + 2 \le 22$, canonical code length $5\,\operatorname{ch} + 2 \le 102 < 112$). Add `14: 393075` (and $n = 15$ once counted) to `_A113201_HINTS` in `fast_enum.py` for single-pass buffer allocation.
4. **`TP_HT_SIZE 8192` (line 4041) & `max_buckets = 4096` (`scan_triplane_level_c`)**:
   - Open-addressed hash table aggregating lossless 128-bit `(e_vec, H1_SNF, Tietze_rel_hash)` buckets per `(c1, c2, c3)` tri-plane level.
   - Sized for $(b \le 4, c_{\text{tri}} \le 6)$ where every level had $< 200$ buckets.
   - **Guardrail**: Before scanning $c_{\text{tri}} \ge 7, 8$ or $b = 5$, **heap-allocate `TpLevelHashTable` with `TP_HT_SIZE 262144` and `max_buckets = 131072`**, and assert `returned_buckets < max_buckets` in Python so `tp_record_resolved` never silently drops new buckets!

---

## 4. Prioritized 64-Core Compute Playbook (What to Run Next!)

Ordered from immediate high-certainty wins (~15 min) to ambitious multi-hour campaigns that produce genuinely new 4-dimensional topology results:

### Campaign 1 (~10–20 min on 64 cores): Complete Unconditional $\operatorname{ch} = 12$ Prime Census Across All Topologies
- **Goal**: Produce the **first complete unconditional $\operatorname{ch} = 12$ prime surface-link table** across all saddle counts $m \in \{2, \dots, 8\}$, all orientabilities, and all component counts ($S^2$, $T^2$, $\mathbb{RP}^2$, $K^2$, $\#_k \mathbb{RP}^2$, and 2/3/4-component surface-links).
- **Why**: On 2 vCPUs, the $m=2$ ($S^2$) and low-saddle $\mathbb{RP}^2$ sectors finished in ~3 minutes. With 64 EPYC cores + parallelized chunked C mountain-pass tunneling (`depth_plus1=256, depth_plus2=256, max_bfs_states=4096`) + `ProcessPoolExecutor(64)` 10-layer invariant stratification against composites/unlinks up to $\operatorname{ch} \le 12$ (note `build_reference_catalog` in `ch11_census.py` already builds composites up to $\operatorname{ch} \le 12$!), the entire unconditional $\operatorname{ch}=12$ census is completely within reach.
- **Payoff**: Certifies every prime surface-link at $\operatorname{ch} = 12$, promotes all matching Spoke families (`Tube(4_1)`, `Spin(5_1)`, `Spin(5_2)`, `Tube(L4a1)`, `Mixed(L5a1)`, etc.) directly from the bottom-up Hub, and automatically tightens lower bounds to $\operatorname{ch}_{\min}(F) \ge 13$ for every unmatched Spoke!

### Campaign 2 (~1–3 hours on 64 cores): Targeted $\operatorname{ch} = 13$ (`77,185` shadows) & $\operatorname{ch} = 14$ (`393,075` shadows) Sector Sweeps
Use `enumerate_targeted_census_level_c` (which filters `m_min..m_max`, `target_components`, `target_orientability`, and $\chi(F)$ in $O(n)$ *before* Kauffman butterfly transforms) on three high-value topological sectors:
1. **First Prime Knotted $\mathbb{RP}^2$ in $\mathbb{R}^4$ (`target_components=1, target_orientability=0, target_chi_min=1, target_chi_max=1`)**:
   - *Mathematical hook*: Through $\operatorname{ch} \le 11$ (all $m$) and $\operatorname{ch} = 12$ ($m \le 3$), **zero** prime knotted $\mathbb{RP}^2$-knots exist (`2_1^{-1}` is unknotted, and even-twist $\tau^{2k}(K) \# P_\pm$ collapse to $K \# P_\pm$ Composite!). Pinning down the exact minimum $\operatorname{ch}$-index where a prime knotted $\mathbb{RP}^2$ (such as $\tau^1(3_1) \# P_\pm$ or a roll-spun $\mathbb{RP}^2$) first appears is a crisp, publishable discovery.
2. **$m = 2$ Prime $2$-Knots at $\operatorname{ch} = 13, 14$ (`m_min=2, m_max=2, target_components=1, target_chi_min=2, target_chi_max=2`)**:
   - Certifies exact $\operatorname{ch} \in \{13, 14\}$ for 1-fusion ribbon words $R(w)$ of length $\le 5$, twisted symmetric unions, and twist-spun knots ($\tau^2(4_1)$, $\tau^3(3_1)$, etc.).
3. **Orientable Tori & 2-Component Links at $\operatorname{ch} = 13$ (`target_orientability=1, target_chi_min=0, target_chi_max=2`)**:
   - Captures spun/turned tori, Satoh tubes of 5- and 6-crossing knots/links, and mixed sphere–torus links.

### Campaign 3 (~30–90 min on 64 cores): Crushing Equivalence Packets (`11_1.{a,b,c}` & New $\operatorname{ch} = 12$ Packets)
Attack unresolved `NEW_PACKET` and `AMBIGUOUS_SHADOW` cells from both sides simultaneously:
- **From above (move equivalence)**: Write a parallel 64-thread C bidirectional BFS / Tabu-guided mountain-pass search through $\operatorname{ch} + 2$ and $\operatorname{ch} + 3$ ($\Omega_1^+, \Omega_6^+, \Omega_2^+$ + level $\Omega_3/\Omega_4/\Omega_5/\Omega_8$ + 4D Tait flypes, `max_bfs_states = 65536..262144`) testing whether `11_1.a`, `11_1.b`, `11_1.c` connect over a $\operatorname{ch}=13$ or $14$ saddle pass, and sweeping the `93` $\operatorname{ch}=11$ `AMBIGUOUS_SHADOW` diagrams against canonical Hub representatives.
- **From below (stronger algebraic splitters)**: Add fast hom-count and twisted-invariant layers over `ProcessPoolExecutor(64)` for non-solvable / linear targets ($\operatorname{PSL}(2, 7)$ of order $168$, $\operatorname{PSL}(2, 8)$, $\operatorname{SL}(2, \mathbb{F}_5)$, $S_5$), twisted Alexander polynomials over $\operatorname{SL}(2, \mathbb{F}_p)$ representations, and higher dihedral/tetrahedral quandle 3-cocycle state sums ($R_5, R_7$) to split genuine 4D mutant or chirality pairs.

### Campaign 4 (~1–4 hours on 64 cores): Braided Tri-Plane Census at $(b = 4, c_{\text{tri}} = 7, 8)$ & Exact $c_{\text{tri}}^{\text{braid}}(9_1)$
- **Goal**: Extend [`triplane_census.py`](./code/src/surfaceknots/triplane_census.py) beyond Gong–Lewis-Monkman–Osnes's $c_{\text{tri}} = 6$ frontier to $c_{\text{tri}} \in \{7, 8\}$ at $b = 4$ (and concentrated $(0, 0, c_3)$ at $b = 5$).
- **How**: First widen `TP_HT_SIZE` in `_fast_enum.c` (see Guardrail 4 above), then dispatch the independent `(c1, c2, c3)` C99 `scan_triplane_level_c` calls concurrently across `ThreadPoolExecutor(max_workers=64)`.
- **Payoff**: At $c_{\text{tri}} \le 6$, `8_1` ($c_{\text{tri}}=4$) and `10_2` ($c_{\text{tri}}=6$) appeared, while the 1-fusion ribbon trefoil `9_1` ($\operatorname{ch}=9$) was completely absent! Determining exact $c_{\text{tri}}^{\text{braid}}(9_1) \in \{7, 8\}$—alongside the full $c_{\text{tri}} = 7, 8$ surface-link spectrum and meander-interleaving obstruction rates—directly advances the 2026 bridge-trisection literature.

### Campaign 5 (~30–60 min on 64 cores): Parallel Spoke Compiler $\operatorname{ch}_{\max}$ Tightening
- Many Spoke entries in `family_catalog.json` carry conservative upper bounds (e.g. $[12, 16]$ or $[12, 22]$) from single-presentation compilation + shallow greedy simplification.
- Dispatching 64-core randomized Reidemeister/braid-word pre-simplification + deep simulated-annealing Yoshikawa reduction (`simplify_compiled_diagram` in `families.py`) across equivalent classical braid/plat/ribbon words will pull $\operatorname{ch}_{\max}$ down sharply—meeting $\operatorname{ch}_{\min} \in \{13, 14\}$ from Campaigns 1–2 to certify dozens of new **exact** $\operatorname{ch}(F)$ values!

### Campaign 6 (Stretch Breakthroughs)
1. **Unconditional (Non-Braided) $b = 4, c_{\text{tri}} \le 6$ Tangle Enumeration**: Replace pure-braid closure generation in `enumerate_minimal_braided_tangles` with general planar 8-boundary tangle shadow generation up to $6$ crossings. Combined with our existing C99 triple scanner, this removes the sole remaining hypothesis ("braided") separating our $b \le 4, c_{\text{tri}} \le 6$ census from an **unconditional computer-assisted proof of GLMO (2026) Conjectures 3.6 & 3.7 at $b = 4$**!
2. **Ideal 4D Triangulation ($f_4$ Pentachora) Compiler**: Build a combinatorial `MarkedDiagram` / `TriPlaneDiagram` $\to$ 4-manifold complement triangulation builder to compare $\operatorname{ch}(F)$, $\widetilde{\operatorname{ch}}(F)$, $c_{\text{tri}}(F)$, and ideal pentachoron complexity $f_4(S^4 \setminus F)$ (Burke–Burton–Ray–Spreer–Thompson–Zymaris 2026).

---

## 5. Ready-to-Paste Starter Prompt for the New 64-Core Session

```text
We are continuing research on `surfaceknots` on a 64-core AMD EPYC 7502 machine.
First read `HANDOFF.md`, `README.md`, `math/README.md`, and `code/README.md`, bootstrap `code/.venv` with `uv` so `_fast_enum.c` compiles with `-O3 -march=native` for this Zen 2 CPU, and verify that all 29 pytest tests pass.
Then inspect the 64-core scaling notes and compute campaigns in `HANDOFF.md` (Sections 3 & 4), upgrade the parallel chunking / tree-merge harness in `fast_enum.py` for 64 cores, and let's start with Campaign 1 (the complete unconditional ch = 12 prime surface-link census across all topologies)!
```
