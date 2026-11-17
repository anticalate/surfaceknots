# Post-Audit Roundup: Deep Move Bridges, Geometric Identification of $12_2$, & Profiling (`ch <= 12 -> 13`)

This directory preserves all scripts, C99 harnesses, verified JSON move trajectories, and mathematical proof scripts from the three post-presentation deep-dive tracks, so everything can be resumed seamlessly on any machine after `git clone`.

---

## 1. `move_bridges/` — Bidirectional Mountain-Pass Search & Residual Audit

- **Engine & Scripts:** [`wide_search.c`](./move_bridges/wide_search.c), [`test_wider_move_bridges.py`](./move_bridges/test_wider_move_bridges.py), [`test_wider_residual_collapse.py`](./move_bridges/test_wider_residual_collapse.py), [`audit_soundness_checks.py`](./move_bridges/audit_soundness_checks.py)
- **Verified JSON Certificates:** [`wide_packet_results.json`](./move_bridges/wide_packet_results.json), [`residual_stage1_ch11.json`](./move_bridges/residual_stage1_ch11.json), [`residual_stage1_ch12.json`](./move_bridges/residual_stage1_ch12.json), [`wide_residual_results_ch11.json`](./move_bridges/wide_residual_results_ch11.json), [`wide_residual_results_ch12.json`](./move_bridges/wide_residual_results_ch12.json)
- **Key Verified Results (`100%` checked step-by-step via `verify_move_path_independent`):**
  1. **Zero False Primes (`75 / 75` Negative Control):** None of the `12` ($\operatorname{ch}=11$) or `63` ($\operatorname{ch}=12$) prime diagrams reduced below $\operatorname{ch}$, split, or hit a 2-cut across $>7.8\times 10^6$ mountain-pass states up to $\operatorname{ch}=15..16$.
  2. **`5 / 15` Multi-Diagram `NEW_PACKET`s Merged (`46 -> 40` components), Promoting Two to `EXACT_SINGLETON`:**
     - **`12_2` ($S^2$, $\mathbb{F}_9$ non-ribbon $2$-knot): `2 -> 1` (`EXACT_SINGLETON`)** via a 16-step $+3$ pass (`(12,10,2) -> (14,12,2) -> (13,11,2) -> (15,13,2) -> (13,11,2) -> (12,10,2)`).
     - **`12_2^{1,1}` ($T^2 \sqcup T^2$): `2 -> 1` (`EXACT_SINGLETON`)** via a 19-step $+3$ pass at $\max\operatorname{ch}=15$.
     - **`12_4^{0,1}`: `4 -> 2`** and **`12_6^{0,1}`: `3 -> 2`** at $+0$ ($\max\operatorname{ch}=12$) via interleaved 4D Tait flypes; **`12_3^{1,1}`: `4 -> 3`** at $+3$.
     - Updates $\operatorname{ch}=12$ prime breakdown from `20` singletons / `14` packets (`63` diagrams) to **`22` singletons / `12` packets (`57` diagrams)**.
  3. **Residual Collapse (`742 -> 503` diagrams / `488` move components at $\operatorname{ch}=12$; `93 -> 74` / `71` at $\operatorname{ch}=11$):**
     - Eliminates **`12 / 55`** ambiguous cells at $\operatorname{ch}=12$ (including **all** $\operatorname{ch}=11$ $S^2$ residuals `11_1` & `11_2`, `6_1^{0,1}`, `7_1^{0,-2}`, `10_2^{0,-2}`, `9_1^{0,1}`, `split:0_1_u_K^2`, `Sigma_3`, `#7RP^2`, and 3 composite cells) and **`5 / 21`** cells at $\operatorname{ch}=11$.
     - Notes that `676 / 742` ($91.1\%$) of $\operatorname{ch}=12$ residuals match unknots (`285` pure + `391` non-orientable unknots algebraically colliding with $10_3 \# (\#_k\mathbb{RP}^2)$, which smoothly equals $\#_k\mathbb{RP}^2$ by Viro's 1973 theorem).
  4. **Bottlenecks Identified for Integration into `moves.py` / `_fast_enum.c`:**
     - `has_certified_cyclic_unknotted_tube` (`moves.py:1471`) hardcoded `red.num_vertices in (8, 9)` for $\operatorname{ch}=11$; extending to `<= 10` (with recursive verification of both summands on any planar 2-cut!) safely collapses `153` $\operatorname{ch}=12$ unknotted 1-handle residuals.
     - Interleaving `list_4d_tait_flypes` inside level-BFS and adding two-step $+1$ stabilizations (`(c+1, m+1)` and `(c, m+2)`) + shared bidirectional frontier hash tables in `_fast_enum.c`.

---

## 2. `mystery_12_2_and_cs/` — Complete Geometric Identification of $12_2$, Cappell–Shaneson Proofs, & SvK Speedup

- **Scripts & Patch:** [`audit_12_2.py`](./mystery_12_2_and_cs/audit_12_2.py), [`certify_12_2_iso.py`](./mystery_12_2_and_cs/certify_12_2_iso.py), [`check_inverse_12_2_and_disks.py`](./mystery_12_2_and_cs/check_inverse_12_2_and_disks.py), [`test_constructions_12_2.py`](./mystery_12_2_and_cs/test_constructions_12_2.py), [`audit_cs.py`](./mystery_12_2_and_cs/audit_cs.py), [`svk_meridian_speedup.patch`](./mystery_12_2_and_cs/svk_meridian_speedup.patch)
- **Key Verified Results:**
  1. **Geometric & Algebraic Resolution of Open Question 8 (`12_2`, $E_1 = (3, t^2+1)$):**
     - Certified over $\mathbb{Z}[t^{\pm 1}]$ via explicit $7\times 7$ Fox minors $m_1 = t^2 - 3t + 1$, $m_2 = 2t^2 - 3t + 2 \in E_1$ ($m_2 - m_1 = t^2 + 1$, $t^{-1}(m_2 - 2m_1) = 3$) that $E_1(12_2) = (3, t^2+1)$ and $\operatorname{def}(G) = 0$.
     - Certified bidirectional 24-edge isomorphism:
       $$G(12_2) \cong \langle a, q \mid q^3 = 1,\; aqa = qaqaq \rangle \cong (\mathbb{Z}_3 * \mathbb{Z}_3) \rtimes_\varphi \mathbb{Z} \cong \pi_1(S^3 \setminus 4_1)/\langle\!\langle (zx^{-1})^3 \rangle\!\rangle \cong \pi_1(S^3 \setminus 5_2)/\langle\!\langle (zx^{-1})^3 \rangle\!\rangle.$$
     - **Two explicit 4D geometric constructions of $12_2$:**
       - *Asymmetric Ribbon-Disk Doubling:* Smoothing the two markers of $12_2$ in opposite orders splits $(S^4, 12_2) = (B^4_-, D_-) \cup_{(S^3, K_{\text{eq}})} (B^4_+, D_+)$ along **10-crossing composite ribbon knots** $K_{\text{eq}}^{(1)} = 4_1 \# 4_1$ and $K_{\text{eq}}^{(2)} = 5_2 \# (-5_2^*)$, gluing two 1-fusion ribbon disks whose pushout shifts the second meridian by $(zx^{-1})^3$ (explaining $c = 5+5 = 10,\; m = 1+1 = 2 \implies \operatorname{ch} = 12$!).
       - *Fibered Punctured-$(L(3,1) \# L(3,1))$ Bundle (Companion of Kanenobu's 2-Cable $C_2(\tau^2(3_1))$):* Shares the exact abstract group $\langle y, q \mid q^3 = 1, y^2 q y^{-2} = q^{-1} \rangle$ with Kanenobu's 2-cable $C_2(\tau^2(3_1))$, but with non-periodic meridian $a = q^{-1} y$ (`|G/<<a^4>>| = 108` vs `|G/<<y^4>>| = \infty`), obtained by composing the periodic 2-cable monodromy with a **Hendriks–Laudenbach puncture-slide** around an order-3 $L(3,1)$ core loop.
     - **Complete Exclusion Theorem:** Proved $12_2$ is **not** any twist-spin $\tau^k(K)$ ($k \ge 1$), **not** any Litherland twist-roll spin $\tau^{k,\ell}(K)$ ($k \ge 1, \ell \in \mathbb{Z}$, by Equivariant Sphere Theorem + normal-bundle rotation), and **not** any Pao/Plotnick branched twist-spin $\tau^{k,s}(K)$ (infinite-order fiber monodromy).
  2. **Pure-Python Confluent Wreath Knuth–Bendix Proofs on $12_{10}$ (`CS_A0`) and $12_8$ (`CS_A1`):**
     - 36-rule confluent systems certifying $\pi_1(12_{10}) \cong \mathbb{Z}^3 \rtimes_{A_0} \mathbb{Z}$ ($\operatorname{tr} A_0 = 2$) and $\pi_1(12_8) \cong \mathbb{Z}^3 \rtimes_{A_1} \mathbb{Z}$ ($\operatorname{tr} A_1 = 3$), analytic + diagrammatic $F_{21}$ ($\mathbb{Z}_7$) and $F_{57}$ ($\mathbb{Z}_{19}$, $\det(I - A_1^3) = -19$) cover homology, plus discovery that equatorial slices of $12_{10}$ and $12_8$ glue pairs of 1-fusion ribbon disks with the exact groups of $11_2$ and $12_7$!
  3. **Seifert–van Kampen Meridian-Vector Speedup (`svk_meridian_speedup.patch`, applied to `ch12_census.py` & `families.py`):**
     - Replaces brute-force Wirtinger + backtracking across $493{,}521$ connected-sum pairs (`28.9 hours` single-core!) with exact Seifert–van Kampen conjugacy-class meridian dot products (`1.59 seconds` across all $993 \times 10$ finite-group pairs, verified `100.000%` exact on `960` random pairs).
     - Also identified why `atlas.py` listed `105` vs `104` classes before regenerating `family_catalog.json`: `Tube(L4a1)` in the pre-audit `family_catalog.json` carried a stale `e1_fingerprint_hash`, double-counting `Tube(L4a1) = 12_3^{1,1}` (`67 + 37 = 104` once regenerated; note `build_reference_catalog` in `ch11_census.py` can also apply the same SvK pre-filter for sub-minute `export_family_datasets`).

---

## 3. `c99_profile/` — Cycle-Accurate TSC Stage Profiler (`ch <= 12` & Roadmap to `ch = 13`)

- **Profiler & Data:** [`_prof_enum.c`](./c99_profile/_prof_enum.c), [`profile_stages.py`](./c99_profile/profile_stages.py), [`baseline_profile.json`](./c99_profile/baseline_profile.json)
- **Key Findings from `baseline_profile.json`:**
  1. **Zero C Kauffman-Bracket False Positives Through $\operatorname{ch} \le 11$ & $\operatorname{ch}=12, m\in\{2,6\}$:** `false_positives_rejected_by_py == 0` across every profiled slice (`certified_kept == merged_unverified`).
  2. **Two Distinct Bottleneck Regimes Across Saddle Count $m$:**
     - **Low-saddle regime ($m = 1, 2$, e.g. $S^2$ $2$-knots!):** At $\operatorname{ch}=12, m=2$ (`15,882` shadows, `16` threads), C enumeration takes **only `43.0 s` wall** (`91.0%` in Kauffman bracket butterfly), while single-process Python `certify_classical_unlink` on the `13,653` unmerged survivors takes **`244.7 s` (`84%` of wall time!)**. Parallelizing Python `certify_classical_unlink` (or porting R1/R2/R3 unlink certification to C) + in-C R1/R2 + pairwise linking-number pre-reduction makes the **entire $\operatorname{ch}=13, m=2$ $S^2$ sector (`77,185` shadows) runnable in ~5–8 minutes on 64 cores**!
     - **Middle-saddle regime ($m = 5, 6, 7$, which dominates full-$\operatorname{ch}$ runtime):** At $\operatorname{ch}=11, m=6$, Kauffman bracket butterfly is **only `10.4%`** of C CPU time (and `5.3%` at $m=7$), whereas **`stage_bfs_canon_dedup` (`56.5%`) + `stage_b_yosh1step` (`31.0%`) consume `87.5%` of C CPU time**! Thus full-$\operatorname{ch}=13$ speedup requires bit-packed canonical-key / immediate-reducing-move pre-checks *before* full BFS component expansion, rather than bracket micro-optimizations alone.
