# Audit of `surfaceknots` (state of `main` at `3089bb0`)

Auditor: the agent writing in `fable.md`, Sun 15 Nov 2026. Untracked working document; not for commit as-is.
Status of this file: complete as of 10:25 EST. An independent read at 10:15 corrected several counts
(§5, §6.5) and one statement (§3.1); the corrections are marked where they occur. Sections 1–9 describe `main` at `3089bb0`; §10 and §12
record what changed on Carbon's branch afterwards and what I re-tested (§12 is the latest); §11 is my
cross-check of a third-party audit, which corrected one of my own conclusions (§7.3). Everything was run from a private clone,
so the shared tree was not touched (apart from this file, `fable.md`, `audit_scripts/` and `expo/report.*`).
A plain-language version of the conclusions is in `expo/report.pdf`.

How to read the evidence tags:

- **[ran]** I executed something and saw the result. The scripts are in `audit_scripts/` (§8); logs and
  larger outputs are under
  `/tmp/claude-29349/-afs-csail-mit-edu-u-e-ejli-surfaceknots/cbe399a3-0fb9-5aba-8e99-ce46cd6219aa/scratchpad/`,
  which is local disk and will not survive a reboot.
- **[read]** I read the code or the source paper.
- **[lit]** Checked against the published literature by a bibliography pass; the two items I rely on most
  (Yoshikawa's Table II and the Kinoshita counterexample) I re-read myself.

---

## 1. Verdict in one page

The engineering is real, and the central result held up: the ch=11 census reproduces from the committed
code, no unsound step was found on the way to it, and its ten classes are distinct and weakly prime. The
write-up around it is a mixture of correct computation, restated literature presented as new, and a
handful of claims that are wrong. It overstates what has been shown in most headline sentences, one
invariant layer is broken, and the tests do not protect against any of this because they compare outputs
to stored copies of the same outputs.

What holds up:

| Claim | Status |
| :--- | :--- |
| Every diagram in every dataset is admissible (both resolutions are unlinks) | **Confirmed independently** with SnapPy on 249 distinct diagrams, 498 resolutions [ran] |
| The stored counts of homomorphisms `Hom(pi_1, G)` for the 23 + 10 table entries and 93 residuals | **Confirmed independently** with a from-scratch Wirtinger presentation and counter, G in {S3, Q8, A4, S4}: 126 of 126 agree [ran] |
| The pipeline, run blind, recovers Yoshikawa's ch=9 and ch=10 rows | **Confirmed** (§5): no unexpected non-cyclic group profile at ch 9 or 10 [ran] |
| The C kernel is memory-safe on the paths the tests and a full ch=9 scan exercise | **Confirmed** under AddressSanitizer + UBSan [ran] |
| Two of Yoshikawa's Table II entries (`9_1^{1,-2}`, `10_2^{0,-2}`) are misprinted and the repo's values are right | **Confirmed**, but this is Ishii–Nikkuni–Oshiro (2018), not new [lit] |
| No ch=11 or ch=12 enumeration exists in print | **No counterexample found** [lit] |
| The ch=11 result (105 survivors, ten classes, 93 unresolved) | **Reproduces** in a clean 75-minute re-run from commit 3089bb0 (§7.8) [ran] |
| The ten ch=11 classes are weakly prime | **Confirmed** by a finite check the repo had not done (§6.5) [ran] |
| No unlink certification is silently lost at ch=11 | **Confirmed**: 156,353,366 candidates, zero failures [ran] |
| The implemented moves are isotopies | No violation in 851,040 applications along random walks, incl. 583 of Ω8 [ran] |

What does not hold up:

| Claim in the repo | Problem | Severity |
| :--- | :--- | :--- |
| "All 3 Table II errata resolved" | There are 2. "Erratum #2" (`10_1^{0,1}`) compares against a value Yoshikawa did not print (§2.1) | Wrong |
| Prop 2.1 of `complexity_measures.md`: `tau^{k+2}(K) # P = tau^k(K) # P` "Viro 1973" | Open question in the literature. Only `10_3 # P = P` is a theorem (§2.2) | Wrong |
| "10 certified new **weakly prime** classes at ch=11"; `11_1^{-2}` "strictly prime" | The repo's check covers only sums of total ch <= 12; weak primeness (Yoshikawa Def 5.1) needs all sums of two surfaces of ch <= 10. **I completed that check (§6.5): all ten are weakly prime**, but two of them have the same group as a composite the repo never compared with. "Strictly prime" is not shown (§3.1) | Incomplete in the repo; now done |
| "First prime knotted RP^2" campaign; description of Yoshikawa §5 | Premise wrong (`tau^1(K)` is the unknot); §5 text is invented; Kinoshita's conjecture was disproved in May 2026 and the repo does not know (§2.3) | Wrong / stale |
| "7 certified new prime 2-knots at ch=12" | The dataset has no generator in the repo (§4.1), though calling the census function directly does reproduce its 29 survivors (55,038 candidates, 1,158 level-irreducible); the classes are known ribbon knots, spins and twist-spins (§2.4) | Not reproducible |
| "Confirming GLMO Conjectures 3.6 & 3.7 on the braided b<=4 locus" | The conjectures quantify over all bridge numbers. The census is consistent with them, nothing more (§2.5) | Overclaim |
| "Meander Interleaving Obstruction forcing `10_2` up to ch=10" (Prop 4.1) | The combinatorial fact is true; the conclusion does not follow. The same obstruction occurs for diagrams of the unknotted sphere (§7.4) | Overclaim |
| "10-layer invariant engine", all layers isotopy invariants | Layer 9 (`double_linking_mod2`) changes under one Reidemeister II move (§7.1). It is the only layer separating `8_1^{1,1}` from `10_1^{1,1}` | **Bug** |
| "E_1 over Z[H_1]" as a key layer | The key uses an evaluation fingerprint, which filed two knotted ch=12 spheres as "matches the unknot"; and for links with both orientable and non-orientable components it depends on the presentation (§7.3) | **Bug** |
| "36,166,398 certified valid b=4 triples" | Produced with an unsound pruning rule; the unpruned run has 62,974,883 and the same classes (§7.4) | Number wrong, conclusion survives |
| ch=11 funnel table (75,125,028 raw survivors, ...) | The first two rows do not reproduce (66,017,754). Everything from the global merge on does, including 105 survivors, the ten classes and the 93 residuals, in a clean 75-minute re-run (§7.8) | Two stale numbers; **result reproduces** |
| "Zero uncertified heuristics" | Failed certifications and failed invariant computations are turned into verdicts silently (§7.2). Not observed to fire | Latent |
| "All 23 certified `EXACT_SINGLETON`" | True by construction: 23 hand-entered diagrams. Run blind, the pipeline leaves `10_1`, `10_2`, `10_3` as 2 unmerged diagrams each and `10_1^1` as 3 (§5) | Misleading |
| "Self-contained" web atlas | Loads Tailwind from `www.gstatic.com/antigravity/web/dev/tailwindcss.min.js` | Minor |

Things the audit turned up that are worth having:

- **`pi_1(S^4 \ 11_1^{-2})` is the binary octahedral group 2O (order 48), with meridians of order 8**
  (§6.1). The repo describes it only as "a 2-generator 3-relator group". This proves, with no enumeration,
  that `11_1^{-2}` has no standard projective-plane or unknotted Klein-bottle summand.
- **`11_1^{-2}` is the spun trefoil with one non-orientable 1-handle attached** (§6.2), and several of the
  93 unresolved ch=11 diagrams are also "knotted sphere + 1-handle". An exact cord-class computation
  settles three of them and identifies three more with each other (§6.3, §9 item 2).
- **Seventeen of the 33 entries through ch=11 are two or three unknotted spheres with one tube**, and can
  be named by a word read exactly from the diagram: R(w), T(w), K(w) (§6.4). The names predicted which
  Klein-bottle links would appear at ch=12, and both predicted ones are in Carbon's partial ch=12 data.
- **Two of the ten new classes have the same group and Alexander ideal as composites of table entries**
  (`11_2^{0,1}` and `6_1^{0,1} # 8_1`; `11_1^{0,-2}` and `7_1^{0,-2} # 8_1`), and differ from them only in
  their components (§6.5).
- **A ch=12 diagram whose group is a Cappell–Shaneson group** (Carbon's observation; I proved the group
  by Knuth–Bendix completion in GAP; §10). A twelve-vertex diagram, in the standard 4-sphere, of a 2-knot
  with group Z ⋉ Z^3.

---

## 2. Literature: what is wrong, stale, or already known

### 2.1 The errata [lit, re-read]

Yoshikawa (1994), Table II, p. 521, prints:

| Entry | Printed group | Printed `E_1` |
| :--- | :--- | :--- |
| `9_1^{1,-2}` | `<x,y : yxyx^{-1}>` | `(x+1, y-1)` |
| `10_1^{0,1}` | `<x,y : x^{-1}y^{-1}x^{-1}yxyxy^{-1}>` | `(x-1, y-1)(xy+1)` |
| `10_2^{0,1}` | `<x,y : x^2 y x^{-2} y^{-1}>` | `(x-1, y-1)(x+1)` |
| `10_2^{0,-2}` | `<x,y_1,y_2 : x y_1 x^{-1} y_2^{-1}, y_1^2 = y_2^2 = (y_1y_2)^2>` | `(2x+y-1, 4)` |

- The repo claims Yoshikawa printed `(x-1,y-1)(x+1)` for `10_1^{0,1}`. He printed that for `10_2^{0,1}`.
- The repo's "corrected" `(x-1,y-1)(x+y)` equals the printed `(x-1,y-1)(xy+1)` after `y -> y^{-1}`.
- Ishii–Nikkuni–Oshiro (Osaka J. Math. 55 (2018), arXiv:1503.07969): "We correct three mistakes: the knot
  group of `9_1^{1,-2}`, its first Alexander ideal and the first Alexander ideal of `10_2^{0,-2}`."
- Where it appears: `README.md` highlights; `HANDOFF.md` §2 item 1; `math/yoshikawa_table.md` header, Part B
  row and §2; `math/foundations_and_invariants.md` §4.2 Prop 4.2; `code/README.md`;
  `code/src/surfaceknots/yoshikawa_data.py` (`yoshikawa_printed_e1` for `10_1^0,1`);
  `code/tests/test_yoshikawa_table.py::test_yoshikawa_three_table_ii_errata_live`, which passes only because
  the "printed" value in the data file is wrong.

### 2.2 The "Viro collapse" [lit]

- Theorem (Viro 1973; Yoshikawa Remark 2.8): `10_3 # P = P`. The repo's use of this is fine, except that
  `complexity.py:288` sets `certified_true_ch = 2` by matching the label `"10_3"`, not by a move path.
- Not a theorem: `tau^{k+2}(K) # P = tau^k(K) # P` for all K, k. Jabłonowski (arXiv:1505.07995, Remark 3.15)
  states that equivalence of `tau^0(3_1) # P` and `tau^2(3_1) # P` is open; Bae–Carter–Choi–Kim
  (arXiv:1407.6028) pose the general question. Baykur–Sunukjian prove stable statements only.
- So `10_2 # P = 8_1 # P` is unproven, and §2 of `complexity_measures.md` must separate one theorem from
  one open problem. This is also an opportunity: see §9, item 3.

### 2.3 Projective planes and Kinoshita's conjecture [lit, re-read]

- Hughes–Kim–Miller–Nahm, "An irreducible real projective plane in the 4-sphere", arXiv:2605.12921
  (13 May 2026): a counterexample to Kinoshita's conjecture. Not cited anywhere in the repo.
- `complexity_measures.md` §3.2 says Yoshikawa §5 "exhibited a prime RP^2-knot candidate with ch <= 13 and
  asked whether any new prime RP^2-knot exists at ch = 11 or 12". No such passage exists. §5 says only that
  the existence of an irreducible projective plane is open.
- The same section refers to "both known prime ch=10 RP^2-knots `10_1^-2` and `10_2^-2`". There are no
  such entries. The only connected non-orientable entry with ch <= 10 is `2_1^{-1}`.
- `HANDOFF.md` Campaign 2 item 1 proposes `tau^1(3_1) # P` as a prime knotted RP^2. `tau^1(K)` is unknotted
  (Zeeman), and any `K # P` is reducible by definition.
- "Conway–Powell (2021)" for RP^2 with group Z_2 should be Conway–Orson–Powell, arXiv:2306.12305.
- Yoshikawa constructed irreducible e=0 Klein bottles (Proc. AMS 126 (1998) 3727–3731), so `11_1^{-2}` is
  at most "first in ch-order".

### 2.4 Results presented as new that are in the literature [lit]

- `complexity_measures.md` Thm 1.2 and Cor 1.3 (unknot ch floors; `ch(P+ # P-) = 3`): Yoshikawa 1994,
  Remark 2.8.
- The new 2-knots `11_1`, `11_2`, `12_m2_S2_1..7` are, by the repo's own identification, 1-fusion ribbon
  2-knots `R(w)`, spun knots and 2-twist-spins. These are tabulated by T. Yasuda, Kanenobu–Komatsu,
  Kanenobu–Sumi, and Kanenobu–Takahashi (Topology Appl. 301 (2021) 107521). What may be new is their
  ch-index. No cross-match against those tables has been done.
- The two real errata: Ishii–Nikkuni–Oshiro 2018.

### 2.5 Gong–Lewis-Monkman–Osnes (arXiv:2606.03799) [lit]

Exact statements: Thm 3.3 (every 2-knot with a bridge trisection of at most five crossings is ribbon, no
bridge restriction); Cor 3.5 (`c_tri(tau^2(3_1)) = 6`); Conj 3.6 (`c_tri(spun trefoil) = 6`); Conj 3.7
(six or fewer crossings implies unknot, spun trefoil or 2-twist-spun trefoil).

- The repo quotes them accurately but then says its braided b <= 4 census "confirms" 3.6 and 3.7. It
  cannot: both quantify over all bridge numbers and all (not only braided) tangle diagrams.
- `yoshikawa_table.md` lists `c_tri >= 6` for `9_1` and `10_1`. That is Conjecture 3.7, not a result.
- `README.md` (image caption) and `HANDOFF.md` Campaign 4 say `c_tri(8_1) = 4`. That would disprove
  Conjecture 3.6. `complexity_measures.md` and `yoshikawa_table.md` say 6.
- `literature_review.md` §4.1 attributes a finiteness question to GLMO that the paper does not pose.

### 2.6 `REFS.md` repairs [lit]

Wrong arXiv ids (the id given resolves to an unrelated paper):

| Citation | Given | Correct |
| :--- | :--- | :--- |
| Kearton–Kurlin | 0802.2775 | 0801.3647 |
| Budney–Burton–Hillman | 1108.5961 | 1109.3899 |
| Meier–Zupan | 1507.08625 | 1507.08370 |
| Joseph–Meier–Miller–Zupan | 2105.03523 | 2112.11557 |
| Conway–Powell | 2009.13426 | 2009.13461 |
| Baykur–Sunukjian | 1309.2543 | 1504.04086 |

Invented or wrong attribution:

- "Conway–Powell, 'Topological unknotting of surfaces in the 4-sphere', JEMS 25(9)": no such paper. The
  real one is "Embedded surfaces with infinite cyclic knot group", Geom. Topol. 27 (2023) 739–821.
- "Ishikawa–Koda–Naoe": H. Naoe alone, "Shadows of 2-knots and complexity", AGT 24 (2024), arXiv:2208.09173.
- "Kang–Park–Taniguchi (2025)": J. Kuhrman, arXiv:2507.03798.
- Authors: Sherry Gong, Samuel Lewis-Monkman, Jesse Osnes; Finn Thompson, Orion Zymaris; Seonmi Choi;
  Taizo Kanenobu and Kota Takahashi (article 107521, not 107778).
- "Allred–Joseph–Meier–Zupan (2023)" is the wrong author list for Allred et al., arXiv:2302.08431.
- Kim 2020 (immersed surface-links) is Topology Appl. 279 (2020) 107241, not JKTR.
- "Kawamura proved t(K) not in {2,3}": Satoh, Osaka J. Math. 42 (2005).
- Double linking: Sanderson, J. London Math. Soc. 35 (1987); Carter–Kamada–Saito–Satoh, AGT 1 (2001).
- "Unique prime decomposition holds for n >= 3 (Milnor 1962)": false (Bayer, Comment. Math. Helv. 55 (1980)).
- "Uniform Finiteness Theorem (Yoshikawa 1994)": a remark in §4, not a stated theorem.

Confirmed as cited (no action): Yoshikawa 1994 (23 entries; Q2.4, Q2.9, Thm 3.5, Lemma 3.9, Def 5.1);
OEIS A113201 values 803, 3378, 15882, 77185, 393075; Burke–Burton–Ray–Spreer–Thompson–Zymaris
(arXiv:2607.05923) Table 1; Allred et al. Table 1 including the misprinted bridge number of `2_1^1`;
Cazet's triple-point bounds; Mochizuki 2003; Kamada–Oshiro 2010; Satoh 2000; Lomonaco 1981.

---

## 3. Logic of the census claims

### 3.1 What `EXACT_SINGLETON` and "new class" actually certify [read]

`stratify_and_export_ch11` (`ch11_census.py:307`) groups the surviving ch=11 diagrams by invariant key. A
key absent from `build_reference_catalog` is "new"; a new key with one surviving diagram is
`EXACT_SINGLETON`. Granting that the moves and invariants are sound (§7), this proves:

> Among surfaces with ch <= 11, exactly one has this invariant key, and its ch-index is 11.

(Correction after an independent read at 10:15: the uniqueness half also needs the key to match no
connected sum of two surfaces of smaller ch-index. A composite of ch-index exactly 11 all of whose
11-vertex diagrams have a two-edge cut is never enumerated, since only prime shadows are. For the ten
classes this is what §6.5 checks; for `11_2^{0,1}` and `11_1^{0,-2}` the key has to be read together with
the knot types of the components.)

By itself it does not prove the class is weakly prime, and it says nothing about primeness. Yoshikawa's
definitions (1994, §5, read from the scan):

- *irreducible*: not `F1 # F2` with `F2` a standard `P+`, `P-` or standard torus;
- *prime*: not a sum of two surfaces neither of which is a trivial 2-knot ("in general it is not easy to
  determine whether a given surface is prime");
- *weakly prime* (Def 5.1): not `F1 # F2` with `ch(F1) < ch(F)` and `ch(F2) < ch(F)`.

So for a class with ch = 11, weak primeness means: not a sum of two surfaces each of ch <= 10. That is a
**finite** condition, because there are finitely many surfaces of ch <= 10. The repo checks a smaller set:
`build_reference_catalog` contains sums `A # B` of table entries and small unknots with
`ch(A) + ch(B) <= 12`. Sums such as `8_1 # 8_1` (16), `6_1^{0,1} # 8_1` (14), `8_1 # 9_1` (17) are not
in it. They cannot be skipped on the grounds that the sum "obviously" has large ch: ch is strictly
subadditive (the repo's own `ch(P+ # P-) = 3`, and `ch(10_3 # P) = 2`).

I carried out the full check; it is in §6.5. All ten classes are weakly prime, but two of the links
share their group and `E_1` with composites that the repo's catalogue did not contain.

"Strictly prime" (the wording used for `11_1^{-2}`) is a different and much stronger statement, since a
summand may have any ch-index. Where a real argument exists, give it:

- 2-knots with irreducible Alexander polynomial and a one-relator group (`11_1`, `11_2`): a short argument
  is available but has to be written.
- `11_1^{-2}`: no standard P or K^2 summand, by meridian order (§6.1).
- Links whose components are all unknotted and whose proper sublinks are identified (`11_1^{0,0,-2}`):
  the sublink argument in `ch11_expansion.md` is fine for non-splitness.

### 3.2 What the ch=11 count means [read]

"10 new classes" is a lower bound on new classes only in the sense of new invariant keys. The 93
`AMBIGUOUS_SHADOW` diagrams are ch=11 diagrams that the search could not reduce and whose invariants match
something known. Each is either a hard diagram of the matched surface or a new surface that these
invariants cannot see. The write-up says this in `ch11_expansion.md` §5 but not in the README, which
states the count flatly.

Breakdown of the 93, from the committed `expanded_ch11.json` (first recorded match of each): 80 match an
unknotted surface (65 of higher genus, labelled `unknot:`; 12 the unknotted torus `2_1^1`, 2 the
projective plane `2_1^{-1}`, 1 the sphere `0_1`), 4 a knotted table entry (`9_1^{0,1}`, `10_2`,
`10_2^{0,-2}`, `10_1^1`), 8 a connected sum, 1 a split union.

At ch=10, where the answer is known, the same pipeline leaves 17 unresolved diagrams with cyclic group and
one `S^2 u T^2` diagram with the group of `6_1^{0,1}` (§5). So "unresolved, matches something known" is the
normal state of affairs for this search depth, not evidence of anything exotic.

### 3.3 Dependence on Yoshikawa's table [read]

Lower bounds ("not realised at ch <= 10") rest on Yoshikawa's hand enumeration being complete. The blind
runs in §5 support it at ch 9 and 10 at the level of group hom counts. This dependence should be stated
once, plainly, in `foundations_and_invariants.md` §3.

The tube-pruning rule `has_certified_cyclic_unknotted_tube` (`moves.py:1336`) also leans on the table: it
treats any connected surface diagram that reduces to <= 7 vertices as unknotted, which is true only because
the table has no knotted connected surface below ch 8.

---

## 4. Reproducibility

### 4.1 Datasets without a generator [ran]

- `code/data/complexity/targeted_ch12_13_census.json`: nothing in `src/`, `scripts/` or `tests/` writes it.
  It carries the 7 ch=12 2-knot classes, the three "promoted to exact ch=12" spokes and both RP^2 censuses.
  The file name says `ch12_13`; the content is ch 11 and 12.
- `math/complexity_measures.md` §3.2 quotes "50,835 raw candidates -> 139 level-minimal" for the ch=11 RP^2
  run. The JSON has `raw_cands: 97824`, `kept: 496`.
- `reproduce_datasets.py` has no flag for the tri-plane census or the complexity comparison; they are
  reachable only by calling `run_triplane_b4_census` / `export_complexity_comparison_datasets` by hand.

### 4.2 Environment

- `uv venv code/.venv` as written in README and HANDOFF picks Python 3.10 on this machine and the install
  fails (`requires-python >= 3.11`). It needs `--python 3.11`.
- No dependency pins or lock file. No provenance (git hash, parameters, runtime) inside the JSON outputs.
- `HANDOFF.md` §1 gives the repository root as `/home/antigravity/anticalate/surfaceknots`.
- Two EPD strings in `triplane_b4_census.json` (`O[1]O[1]`, `O[1]O[1]O[1]`) are not parseable by
  `MarkedDiagram.from_epd`: the writer drops the separator between free loops.

### 4.3 Tests [read]

29 tests pass in about 10 s. Their weakness is structural:

- Assertions on headline numbers read the committed JSON (`len(entries) == 10`,
  `exact_singletons_count == 9`, all 23 `EXACT_SINGLETON`). They detect accidental edits, not errors.
- `test_yoshikawa_three_table_ii_errata_live` encodes the wrong printed value (§2.1) and passes.
- `test_viro_collapse_and_rp2_stabilization` asserts a value that the code sets from a label (§2.2).
- There is no blind enumeration test at any ch. One at ch <= 9 costs 12 s (§5).
- There is no test that a move preserves invariants on a large random sample, and no differential test of
  the C moves against the Python ones.

---

## 5. Blind recovery of Yoshikawa's table [ran]

`enumerate_targeted_prime_ch_census(ch, m_min=1, m_max=ch-1, filter_cyclic_tubes=False)`, all topologies,
default tunnelling depths; then `|Hom(pi_1, G)|` for G in {S3, Q8, A4, S4, A5} on every survivor.

**ch = 9** (198 shadows, 489,273 raw candidates, 12.6 s on 8 workers): 12 survivors.

| Profile | Diagrams | Meaning |
| :--- | ---: | :--- |
| S^2, non-cyclic | 1 | `9_1` |
| T^2 u S^2, non-cyclic | 1 | `9_1^{0,1}` |
| T^2 u K^2, non-cyclic | 1 | `9_1^{1,-2}` |
| one-generator group: T^2 (2), genus 2 (3), K^2 e=0 (4) | 9 | unknot-profile leftovers; 6 removed by the tube rule, 3 Klein bottles remain |

**ch = 10** (803 shadows, 6,256,765 raw candidates, 179 s): 54 survivors.

| Profile | Diagrams | Meaning |
| :--- | ---: | :--- |
| S^2, three distinct non-cyclic profiles | 2 + 2 + 2 | `10_1`, `10_2`, `10_3`, each as two unmerged diagrams |
| T^2, trefoil-group profile | 3 | `10_1^1`, three unmerged diagrams |
| T^2 u S^2, S3-count 30 | 2 | `10_1^{0,1}` and `10_2^{0,1}` (same hom counts; split by `E_1`) |
| S^2 u K^2 | 2 | `10_1^{0,-2}` and `10_2^{0,-2}` (same hom counts; split by `E_1`) |
| T^2 u T^2, group Z^2 | 1 | `10_1^{1,1}` |
| T^2 u S^2 u S^2 | 1 | `10_1^{0,0,1}` |
| RP^2 u RP^2 | 1 | `10_1^{-1,-1}` |
| K^2 u K^2 | 1 | `10_1^{-2,-2}` |
| T^2 u S^2, group Z^2 | 1 | not in the ch=10 row: profile of `6_1^{0,1}` (unreduced diagram, presumably) |
| one-generator group, 6 topological types | 36 | unknot-profile leftovers; 19 removed by the tube rule, 17 remain |

(Recounted at 10:16 from `audit_scripts/tube_ch10.json`: 17 diagrams for the 12 table classes, 1 with the
`6_1^{0,1}` profile, 36 with a one-generator group; an earlier version of this table said 35.)

Conclusions:

- Every ch=9 and ch=10 class of Yoshikawa's table appears, and no non-cyclic profile appears that is not in
  the table. This is good evidence for both the pipeline and the table.
- The search does **not** merge all diagrams of a class even at ch=10. So `NEW_PACKET (m=3)` for `11_1` is
  the expected behaviour for a single class, and `EXACT_SINGLETON` is partly luck.
- In both runs the tube rule fired only on diagrams with a one-generator group, as it should.

Suggested permanent tests: ch <= 9 on every run; ch = 10 marked slow.

---

## 6. New observations from the audit

### 6.1 The group of `11_1^{-2}` [ran]

From the presentation stored in `expanded_ch11.json`,
`<x1, x2 | x1x2x1 = x2x1x2, x2^-1 x1^-1 x2 x1^-1 x2^-1 x1, x2^-2 x1^-1 x2^2 x1^-1>`:

- Order 48 by coset enumeration. Element orders 1:1, 2:1, 3:8, 4:18, 6:8, 8:12. One involution, so the
  group is the binary octahedral group 2O and not GL(2,3). Derived series 48 > 24 > 8 > 2 > 1.
- Both meridians have order 8 and each is conjugate to its inverse.
- Any two of the three relators already present the group (Carbon's observation; I re-checked it).
- Centre Z_2, quotient S_4. This is why every hom count in the table matched S_4.
- Confirmed at 09:55 with GAP from my own Wirtinger presentation of the diagram
  (`audit_scripts/finite_groups.py`): `IdGroup = [48, 28]`, the binary octahedral group; one conjugacy
  class of meridians, of order 8.
- The same script gives the group of every other entry whose components are all non-orientable (the
  only entries whose group can be finite): `2_1^{-1}`: Z_2; `8_1^{-1,-1}`: Q_8; `10_1^{-1,-1}`: the
  generalised quaternion group Q_16; `10_1^{-2,-2}`: Q_8. In the three links the meridians have order 4
  and the two components give two conjugacy classes. This agrees with Yoshikawa's Prop. 5.3 for the
  first two (meridians of order 4, hence no standard projective-plane summand).

Consequence: a standard P summand or an unknotted Klein-bottle summand forces a meridian of order 2.
So `11_1^{-2}` is not `K # P+ # P-` for any 2-knot K of any ch-index. This replaces the "strictly prime"
sentence in `ch11_expansion.md` §2, which compared hom counts with three composites.

### 6.2 `11_1^{-2}` is the spun trefoil plus a non-orientable 1-handle [ran]

The stored diagram has two marked bigons (vertices (2,6) and (4,9)). Deleting either and reconnecting the
outer darts gives an admissible 9-vertex S^2 diagram, which reduces to 8 vertices and which
`tunnel_diagrams_c` merges with the canonical `8_1`. In group terms the surface is
`<trefoil group | b^2 a b^-2 = a^-1>`: one relation of the form "meridian = conjugate of inverse meridian",
which is what a non-orientable 1-handle contributes. The handle is non-trivial, since the trivial one
gives `8_1 # K^2` with group S_3.

Caveat: this uses the repo's C tunnelling for the identification of the base with `8_1`, so it is as
sound as §7.2 says that is. (By Yoshikawa's table any knotted 2-sphere with ch <= 8 is `8_1` anyway.)

Carbon has proposed a description as a "twist-turn-spun trefoil". Until a diagram compiled from that
construction is matched to `11_1^{-2}` by moves, or a reference is found, that should be recorded as a
conjecture.

Literature [lit]: Yoshikawa, "The order of a meridian of a knotted Klein bottle" (Proc. AMS 126, 1998)
shows every even meridian order occurs. His order-4 example is a twisted tube on the 2-twist-spin of the
pretzel knot P(-2,3,3); by the presentation in Meier (arXiv:2510.06358) its group also has order 48
(checked by coset enumeration by the lookup agent, not by me), with meridian of order 4. So `11_1^{-2}`
(meridian order 8) is not that surface, but it is the same kind of object: a knotted sphere with one
twisted tube. Whether it is one of Yoshikawa's order-8 examples is open here; his paper was not
accessible tonight.

### 6.3 Residuals that are visibly "base + 1-handle" [ran]

Scanning all 188 catalogue diagrams for marked bigons:

- 12 residuals match the unknotted torus; at least two are `8_1` or `9_1` plus an orientable 1-handle
  whose result has group Z.
- 4 residuals match `8_1 # K^2`; they are `8_1` or `9_1` plus a non-orientable 1-handle with group S_3.
- One residual matching `10_1^1` is `8_1` plus an orientable handle with the trefoil group.

See §9 item 2 for why this matters.

---

### 6.4 Names: fifteen entries are two unknotted spheres with one tube [ran]

Take two unknotted, unlinked spheres A and B, with free group F(a,b) of the complement. One tube along a
cord word w gives `R(w)` (tube from A to B: a ribbon 2-knot, group `<a,b | a = w b w^-1>`), `T(w)` (both
ends on A, orientable: a torus linked with B, `a = w a w^-1`) or `K(w)` (both ends on A, non-orientable: a
Klein bottle linked with B, `a = w a^-1 w^-1`). `audit_scripts/handle_free.py` finds a marked bigon whose
base is the two-sphere unlink (searching level-move neighbours if the stored diagram shows none), realises
F(a,b) by Sanov's matrices in SL(2,Z), assigns a matrix to every arc, and reads the cord word from the
meridians at the two feet. `handle_free_groupcheck.py` confirms for all 16 diagrams that the presented
group has the same hom counts into S3, Q8, A4, S4 as the diagram.

| Entry | Name | Entry | Name |
| :--- | :--- | :--- | :--- |
| `8_1` | R(b a) | `6_1^{0,1}` | T(b) |
| `9_1` | R(b a^-1) | `9_1^{0,1}` | T(b a b^-1) |
| `11_1.b` | R(b a b^-1 a^-1) | `10_1^{0,1}` | T(b a b) |
| `11_1.c` | R(b a b^-1 a) | `10_2^{0,1}` | T(b^2) |
| `11_2` | R(b a b a^-1) | `11_1^{0,1}` | T(b a b a^-1 b^-1) |
| `7_1^{0,-2}` | K(b) | `11_2^{0,1}` | T(b a b^-1 a b) |
| `10_1^{0,-2}` | K(b a b^-1) | `11_2^{0,-2}` | K(b^2) |
| `11_3^{0,-2}` | K(b a b) | `11_1^{0,-2}` | K(b a b^-1 a^-1 b) |

- Agreement with the repo's family compiler where both exist: `8_1` = R(ba), `9_1` = R(ba^-1),
  Mixed(L4a1) = `10_1^{0,1}` = T(bab), KleinMixed(L4a1,1) = `11_3^{0,-2}` = K(bab). The family name of
  `11_1`, R(b^-1 a b^-1 a^-1), is the word of `11_1.c` up to symmetry.
- `ch(K(w)) = ch(T(w)) + 1` in the four cases where both occur (w = b: 6, 7; bab^-1: 9, 10; b^2: 10, 11;
  bab: 10, 11).
- The prediction that K(b a b a^-1 b^-1) and K(b a b^-1 a b) first appear at ch=12 was tested on Carbon's
  partial `ch12_raw_survivors.json` (two-saddle diagrams, 211 of them; `name_ch12.log`): both are there,
  with hom counts not seen for any S^2 u K^2 class at ch <= 11. Six of six. The same file shows a new
  S^2 u T^2 word, T(b a b^-1 a^-1 b), and Mixed(L5a1) reading as T(b a b^-1 a^-1 b^-1 a b).
- The two three-component entries fit the same pattern, read from their diagrams by
  `handle_free3.py` (F_3 realised as an index-2 subgroup of Sanov's group): `10_1^{0,0,1}` is three
  unknotted spheres with an orientable tube on one of them, cord b c (once around each of the other
  two); `11_1^{0,0,-2}` is the same with a non-orientable tube. ch 10 and 11: seven of seven for
  "the Klein-bottle version costs one more vertex".
- Links with a handle on each component (`handle_ab.py`): `8_1^{1,1}` = T(b) + orientable tube on B with
  cord a; `9_1^{1,-2}` = T(b) + non-orientable tube on B with cord a (the base group is Z^2, so the cord
  class is an integer; it is 1 in both); `10_1^{-2,-2}` = K(b) + non-orientable tube on B;
  `11_1^{1,-2}` = T(b a b^-1) + non-orientable tube on B. The canonical diagram of `10_1^{1,1}` shows no
  tube, so this does not yet give an independent separation of `8_1^{1,1}` from `10_1^{1,1}`.
  (The group quoted for `11_1^{1,-2}` in the report, Z_4 ⋊ Z, is read off the repository's stored
  presentation `<x1, x2 | x2^-1 x1^-1 x2 x1^-1, x1^-4>`; the stored presentations agree with my
  independent hom counts on every diagram, but I did not derive this group separately.)
- `10_2^{0,-2}` shows no tube in any diagram examined (the table calls it a non-ribbon candidate).
  `10_1` is ribbon but its diagrams show no bigon. `10_2`, `10_3` are not ribbon.
- **The `11_1` packet contains two different cord words**, R(b a b^-1 a^-1) and R(b a b^-1 a).
  Alexander polynomials t^2-2t+2 and 2t^2-2t+1 (mirror images); equal hom counts into S5, A6, PSL(2,7),
  AGL(1,p) for p <= 41 (`packet11.py`). Literature [lit]: in Kanenobu's notation the second is
  R(1,1,-1,1) = Yasuda's Y3_4 and the first is R(1,1,-1,-1), whose Alexander polynomial occurs only for
  the mirror image Y3_4* among ribbon 2-knots with at most four ribbon crossings (a complete list). A
  ribbon 2-knot is its own mirror image with orientation reversed, so as unoriented surfaces the two
  words agree. So the packet is one class, by inference from the published classification; the move
  search has not confirmed it.
- Scope of the theorem I am using [lit]: Boyle (1988) classifies *orientable* 1-handles on an *oriented
  surface-knot* (any genus) by double cosets of the peripheral subgroup, with g and g^-1 giving the
  same handle; non-orientable handles are Kamada (arXiv:1403.0690, Case 2); surface-*links* are Kamada,
  "Tensor products of quandles and 1-handles attached to surface-links" (arXiv:2003.12977), in terms of
  quandles, with the remark that the knot-group arguments do not apply directly. My double-coset
  reading for a two-sphere base is the natural translation (cords up to homotopy) but I have not found
  it stated in print. The names are checked against groups and against independent constructions; the
  claim "same word, same surface" for links rests on this.

### 6.5 Weak primeness of the ten ch=11 classes: done [ran]

Yoshikawa's Def 5.1: not a sum of two surfaces each of ch <= 10.

- `11_1`, `11_2`: `E_1` is multiplicative under connected sum; `(2x^2-2x+1)` and `(x^3-x^2+2x-1)` are
  irreducible, differ from the principal ideals of `8_1`, `9_1`, `10_1`, and are not contained in
  `(x+1,3)` or `(x^2+x+1,2)`.
- `11_1^{-2}`: meridian of order 8 (§6.1); every connected non-orientable surface with ch <= 10 has a
  meridian of order 2.
- The seven link classes (`audit_scripts/weak_prime.py`): summands = all table entries, the unknotted
  Klein bottle, `8_1 # T^2`, `8_1 # P`; every pair, every pair of components, both splice orientations,
  both mirrors; keep sums with the topology of a target. 458 sums; hom counts into S3, Q8, A4, S4, D8,
  D10, F21 compared.
  - `11_1^{0,0,-2}`: no sum has its topology.
  - `11_1^{0,1}`, `11_2^{0,-2}`, `11_3^{0,-2}`, `11_1^{1,-2}`: no sum has their hom counts.
  - **`11_2^{0,1}` has the hom counts of `6_1^{0,1} # 8_1`, and `11_1^{0,-2}` those of
    `7_1^{0,-2} # 8_1`.** These composites have summand total 14 and 15, so the repo's catalogue
    (total <= 12) never saw them. The groups are in fact isomorphic: from the name T(b a b^-1 a b), with
    c = b a b^-1 the relation becomes a c a = c a c, and b' = (ac)^-1 b commutes with a, so the group is
    B_3 amalgamated with Z^2 along `<a>`, the group of `6_1^{0,1} # 8_1`. `E_1` agrees too.
  - They are nevertheless different surfaces: in the composite one component is knotted (the torus is
    `8_1 # T^2`, or the sphere is `8_1`; the Klein bottle is `8_1 # K^2`), while every component of the
    census class is unknotted (the repo reduces them to <= 7 vertices; in the names above, setting b = 1
    gives a cord that is a power of a). In the repo's key the only sound layer that separates them is
    `pi2_covers` (`weak_prime2.py`).
- So all ten are weakly prime, on the same two assumptions as everything else (Yoshikawa's table is
  complete; the blind runs show no other non-cyclic profile at ch <= 10). The lesson for ch=12: build the
  reference catalogue from all pairs of surfaces of ch <= 11, not from a total-ch cutoff, and compare
  sublinks.

## 7. Code soundness

Four deep checks were run, one each on: unlink certification and the Yoshikawa moves (A); the enumeration
filters and the C kernel (B); the invariant layers (C); the tri-plane census and complexity module (D).
They were stopped before completion to save budget; what was not reached is listed in §7.7. Reproducers
are in `audit_A/` … `audit_D/` under
`/tmp/claude-29349/-afs-csail-mit-edu-u-e-ejli-surfaceknots/cbe399a3-0fb9-5aba-8e99-ce46cd6219aa/scratchpad/`
(local disk; copy what you want to keep). Items marked **[re-ran]** I reproduced myself.

**Summary.** No unsound move and no unsound discard was found in the ch-diagram census path. One
invariant layer is not an invariant, and a second (E_1) depends on the presentation for mixed links,
which I missed until a third-party audit pointed at it (§7.3, §11). Several components turn "could not decide" into a verdict without
recording it. One layer of the key is much weaker than documented and has already hidden a class.

### 7.1 Critical: `double_linking_mod2` is not an isotopy invariant **[re-ran]**

`invariants/linking.py:902-1224`; enters the key at `yoshikawa_data.py:625-633` (layer 9).

Each of these two diagrams reduces to the canonical table diagram by one `find_omega2_reduction` (the
canonical keys are equal), and the value changes:

| Surface | Canonical value | Value on the variant |
| :--- | :---: | :---: |
| `8_1^{1,1}` | (0,0) | (0,1) |
| `10_1^{1,1}` | (1,1) | (0,1) |

```
8_1^1,1:  Z[1,2,3,4], X[6,4,7,5], X[8,9,10,11], X[12,11,10,13], Y[9,14,6,15], X[5,7,17,16], Y[13,15,16,12], X[18,8,19,20], Z[19,17,3,2], X[20,1,14,18]
10_1^1,1: Z[1,2,3,4], X[5,6,7,8], X[10,11,2,9], X[12,13,14,15], X[16,17,6,18], X[8,20,3,19], Z[14,21,16,15], X[21,23,4,22], X[9,1,23,24], X[24,13,12,10], Y[17,22,20,7], Y[19,11,18,5]
```

- Rate over 120 random move-equivalent variants per surface: 21 (`8_1^{1,1}`), 28 (`10_1^{1,1}`),
  12 (`6_1^{0,1}`), 6 (`9_1^{1,-2}`), 7 (`10_1^{0,0,1}`) changed.
- Cause: line 1208 uses `(cross_ab // 2) & 1`, half a crossing count, where a linking number mod 2 is
  needed; lines 1181-1202 and 1210 add an Arf term; lines 1134-1214 use the first slice found.
- Consequence 1: layer 9 is the only layer separating `8_1^{1,1}` from `10_1^{1,1}`. The repo's code does
  not establish that separation, so "all 23 entries pairwise separated" is false as a statement about
  the code. (The two links are distinct; that is Yoshikawa's remark with references [15], [20].)
- Consequence 2: a non-invariant layer can split one surface into two cells. Any novelty claim that
  rests on layer 9 is void. The three ch=11 classes with a torus component differ from every
  same-topology table entry in hom counts, so they do not rest on it; the composite and split
  references were not checked.
- A corrected form (`audit_C/mydlk.py`: sum over the other slice circles of lk mod 2, no Arf term, on the
  non-zero class in the kernel of the linking functional) gives 0 and 1 for the two links and was
  unchanged on 600 of 600 variants. Not proved equal to Sanderson's invariant.
- `foundations_and_invariants.md` Thm 5.2 offers a "proof sketch" that is not a proof of anything the
  code computes.

### 7.2 High: undecided cases become verdicts

1. **Unlink certification.** `certify_classical_unlink` returns `is_unlink=False`, reason
   `unlinkedness_not_certified`, after a bounded search (`unlink.py:717-721`: 512 R3 states, one
   stabilisation), and rejects at <= 7 crossings on a cited theorem that could not be verified
   (`unlink.py:665`). Every caller reads only `.is_unlink`: `fast_enum.py:297, 365, 413, 465`;
   `ch11_census.py:78`; `moves.py:1403`; `complexity.py:104, 216`; `families.py:552`. In
   `merge_diagrams_c` and `tunnel_diagrams_c` the diagram being certified represents a whole merged
   class, so one false negative silently removes a class.
   Measured: 0 failures in the pipeline for ch 6–10 (11,486 calls at ch=10); 0 in 1,034,918 calls over
   all ch=9 candidates; brute force over all marker counts for n <= 9 certifies every bracket survivor
   (6,905,465 at n=9); 0 in the first 17.2M of about 75M ch=11 candidates. So it has not been observed
   to fire. The module docstring's "zero false negatives up to 12 crossings" is unsupported.
2. **Movie-based invariants.** `mochizuki_3cocycle_invariant` returns None above 729 colourings or when
   `find_unlink_movie` fails (`quandle.py:888-889, 898-899`); the triple-linking profile returns None
   (`linking.py:659-660`); Dlk omits pairs (`linking.py:1085, 1105, 1218`); the Q6 routine raises
   (`tier3.py:1342`). A diagram where the search fails gets a different key from another diagram of the
   same surface. Not triggered in our runs. The docstring at `yoshikawa_data.py:562-566`
   ("100% presentation-independent") is false for layers 9 and 10.

### 7.3 High: the E1 layer is a fingerprint, and has hidden a class **[re-ran]**

- The key uses evaluations at t = ±1 and reductions mod t²−1, t³−1 (`alexander.py:313-397`), not the
  ideal. It is weak.
- **Correction (09:30 EST): it is also not presentation-independent for links with both orientable and
  non-orientable components.** I had written here that the layer "cannot separate a wrong pair"; that
  was wrong. The third-party audit found the defect in Carbon's new Gröbner key, and I then reproduced
  it with my own script on the *original* code at `3089bb0` (`audit_scripts/groebner_presentation_dep.py`):
  after the change of generators g_o = y_1, g_n = y_1^-1 y_2 on the simplified presentation (component
  meridians carried along), `e1_generators` and the `evaluation_fingerprint` of `11_2^{0,-2}` become
  those of `11_3^{0,-2}` and vice versa. The cause is that the code can assign the coordinate x y, not
  x, to the orientable meridian; the two recorded ideals differ exactly by x -> x y. Seven other mixed
  links tested did not change. So for such links layer 3 can split one surface into two cells, and the
  statement in `ch11_expansion.md` that E_1 separates `11_2^{0,-2}` from `11_3^{0,-2}` cannot be relied
  on as computed. **The ten ch=11 classes do not depend on it**: each of the five with mixed
  orientability differs from every table entry and every other new class of its topology in hom counts,
  the double branched cover or the kei module (the twin pair in the last two), and the comparison with
  composites in §6.5 used hom counts and components only.
- **How far the defect reaches** (`audit_scripts/presentation_invariance.py`, original code): I wrapped
  `wirtinger_presentation` so that every simplified presentation is changed by three random Nielsen
  moves, with relators and component meridians rewritten, and recomputed the full key six times for
  each of the 32 non-trivial table and ch=11 entries. In 186 completed variants the only layer that
  changed was E_1, and only on `11_2^{0,-2}` and `11_3^{0,-2}` (twice silently; six more times the
  code's own cross-check between the two resolutions raised "Resolution mismatch in E_1 fingerprint").
  Every other layer, including the cover-homology and peripheral layers, was unchanged. So the defect is
  confined to E_1 on that pair in this data, and in a census run it would more often stop with an
  assertion than split a class.
- In `targeted_ch12_13_census.json`, entries 3 and 10 of `ch12_m2_s2_knots.ambiguous_shadow_collisions`
  are filed as matching `0_1`. Their stored `e1_generators` are `['3', 'x**2 + 1']` and their
  presentations have two generators. They are knotted spheres (144 homomorphisms with non-abelian image
  into F_9 ⋊ Z_4, `audit_D/audit_amb2.py`). **So the ch=12, m=2 slice has at least 8 keys not seen
  below, not 7.** Grouping its 29 survivors by (E1 generators, hom counts, Σ2) gives 13 cells against
  the repo's 10.
- At ch=11 nothing is hidden this way: all 80 residuals matched to an unknotted surface have a
  one-generator presentation, and the other 13 have the same `e1_generators` as their reference.
- `10_1^{0,1}` and `10_2^{0,1}` have the same E1 fingerprint (`3f51b7ba8ceeef1f`); the code separates
  them at Σ2. "E_1 alone separates them" is true of the ideals and false of the code.
- `exact_zh1_ideal_equals` (`yoshikawa_data.py:311-361`) returns True on a necessary condition (Hermite
  forms in Z[t]/(t^m−1), m = 3, 4, 5). "Exact agreement on 20 of 23" is not a proof of equality.
- `12_m2_S2_3` (six unmerged diagrams, labelled "`tau^2(4_1)` / `tau^2(5_1)`"): no computed field differs
  between the six. Only the mod-3 Mochizuki cocycle is implemented, which cannot tell these apart.

### 7.4 High: the tri-plane census

Re-run from the repo reproduces the JSON exactly (316 s on one core): 36,166,398 triples, 22 classes.

- **The arch pruning is unsound** (`triplane_census.py:869-880, 1248-1259`; `_fast_enum.c` ~4250-4285).
  The cap test `tau[u]==u+1 and len(top_w[u])==1` is wrong: at b=4, 21,076 of 69,118 flagged bits are
  not boundary-parallel caps and 8,415 real caps are missed. Even with a correct cap test, nothing
  checks that the destabilised tangle is trivial. It removes 26,808,485 valid triples; the trefoil-group
  class keeps 344 triples and has 942 without pruning, and `b(8_1) = 4` means none of those can be
  destabilisations.
- **The conclusion survives.** With arch pruning off (630 s): 62,974,883 valid triples, the same 22
  classes, the same three connected-sphere classes (group Z: 8,094,033 triples; trefoil group: 942;
  group of `10_2`: 90). So within the braided b <= 4 universe the statement is true. The figure
  36,166,398 and `min_c_tri` for split classes (the 4-component unlink is listed at 4, it is 0) are
  artefacts of the pruning.
- `tangle_3d_key` (`triplane_census.py:42-69`) is described as exact 3D-isotopy deduplication. It is
  strictly finer. True isotopy-class counts at b=4 for c = 0..5 are 14, 56, 208, 764, 2810, 10490; the
  repo reports 14, 96, 582, 3134, 15994, 78912. Harmless for exhaustiveness.
- No pair is accepted on the bracket alone: of 32,152,352 pair decisions at b=4, 11,386,807 were reduced
  to zero crossings by Reidemeister moves in C, 20,736,830 rejected, 28,715 certified in Python. No
  certificates are stored.
- Identification of classes is by invariant key only. `8_1` has one move-certified witness; `10_2` is
  matched by key, and the same key matches the spoke `R(b^3*a)`. Spheres with group Z at b=4 are
  topologically unknotted; smooth unknottedness is not shown.
- `c_tri(8_1)`: the data say 6. "4" in the README and HANDOFF is the bridge number.
- **Meander obstruction.** The combinatorial statement is true: an exhaustive search over crossing-free
  band systems agrees with the greedy compiler on all 4, 25 and 196 ordered Catalan pairs for b = 2, 3, 4,
  and exactly four pairs fail at b=4. The topological conclusions are not established: the same failure
  hits 12,822 triples of the unknotted sphere and 37,442 of the standard projective plane. It is a
  property of a diagram. It does not force `ch(10_2) = 10`.
- The compiler itself is fine: 950 random compilable triples and 37 stored witnesses give admissible
  marked diagrams with matching invariants.

### 7.5 Medium and low

- **Marker range.** The scan covers m in [2, ch−4] only (`fast_enum.py:255-258, 486-489, 692-695`;
  `reproduce_datasets.py:83`). The docs state all 4^11 assignments. The excluded ranges were scanned
  through the pipeline for ch 3–11 and nothing survives (ch=11: m in {0,1}: 7,237 raw, 0 after
  tunnelling; m in {8..11}: 12,385,093 raw, 0), so ch <= 11 is unaffected. Unproven for ch >= 12.
- **The ch=11 scan counters are stale; the result reproduces.** With the script's partition and
  arguments: raw survivors 66,017,754 (documented 75,125,028); level components 22,155,075 (documented
  29,916,769). The intermediate files are not in the repo, and `--ch11-export` falls back to the stored
  `expanded_ch11.json` when they are missing (`reproduce_datasets.py:327-334`). A clean re-run of the
  whole pipeline matches the documentation from the global merge onward (§7.8).
- **Tube pruning** (`moves.py:1336-1449`). The rule is a correct use of a real theorem when the base has
  cyclic group, but (a) it tries two reconnections with no check of marker types, and only the one
  joining the two outer darts at each marker is the handle base; (b) "cyclic" is inferred from hom
  counts, vertex counts and an empty tunnelling result rather than certified. Observed behaviour is
  correct: every prune at ch 9 and 10 (5 and 18 in the default pipeline) went through the right
  reconnection and hit a diagram with the hom counts of Z or Z_2; on 188 catalogue diagrams the wrong
  reconnection never passed the numeric checks **[re-ran]**. The 79 prunes at ch=11 were audited afterwards and all hit cyclic-profile diagrams (§7.8).
  `foundations_and_invariants.md` Lemma 3.2 is not a theorem as stated.
- **Printed Klein-bottle diagram.** `complexity_measures.md` Thm 1.2 prints
  `X[1,2,3,4], Y[4,2,5,6], Z[6,5,3,1]`, which is inadmissible (both resolutions are 2-component links
  with non-zero linking number). The code uses `X[1,2,3,4], Y[5,6,2,1], Y[5,4,3,6]`, which is admissible,
  and the 4-step move path to `P+ # P-` preserves invariants.
- `find_planar_connected_sum_decomposition` (`moves.py:1010-1023`) calls a factor non-trivial when it has
  8 or more vertices after greedy reduction; `is_certified_composite` can be a false label. The discard
  is still justified.
- `certified_chiral_by_c3` (`yoshikawa_data.py:459-461`) detects non-invertibility, not chirality.
- EPD: free loops are lost on round trip when vertices are present; `O[1]O[1]` is unparseable;
  `from_epd` skips malformed tokens and checks neither planarity nor that `opp` is an involution.
- `_fast_enum.c:758` skips sectors with more than 14 crossings without error; tunnel-budget comments
  (2522-2523) disagree with the code (2340-2343); allocations are unchecked; the `.so` rebuild can race.
- README names the promoted spokes as `R(b^2a)`, `R(b^2a^-1)`, `SymUnion(5_2,1)`; the JSON and
  `complexity_measures.md` say `Spin(5_2)`, `R(b^2*a)`, `Spin(5_1)`.

### 7.6 Checked and found correct

- **Shadows.** 1, 1, 2, 3, 9, 18, 62, 198, 803, 3378, 15882, 77185, 393075 for n = 2..14. C and Python
  generators agree up to isomorphism for n <= 11. Every C shadow with n <= 11 is planar, connected,
  loop-free, without 2-edge cut, pairwise non-isomorphic. Excluding loops is justified (a loop is an
  Ω1 or Ω6 kink). Completeness rests on the generation theorem for quadrangulations and the count match.
- **Filters.** Kauffman bracket: an exact necessary condition, no overflow, zero admissible diagrams
  rejected by brute force for n <= 9. Pairwise linking number zero: necessary. "Linking parity":
  vacuous. Fixing the first marker and first crossing type: quotient by the two mirrors. Bigon and
  one-step prunes: every removed candidate has an immediate reduction or a validated neighbour with
  one (n <= 8). 2-edge-cut test agrees with brute force.
- **Dedup and buffers.** No hash-only comparison; canonical keys are lossless and partition about
  49,000 diagrams exactly as the Python key does. Output counts are checked and calls re-run with larger
  buffers. Every cap causes over-keeping. The only static is a const table.
- **Classical moves.** R1, R2, R3 and their inverses in `unlink.py`: 9,814 random-walk applications and
  1,366 one-step neighbours preserve planarity, component count, Jones polynomial and Fox colourings.
  No false positive on fifteen non-trivial knots and links and their mirrors.
- **Yoshikawa moves.** Ω1–Ω7 verified by slot-level derivation. C level moves give the same neighbour
  sets as `moves.py` on about 17,700 admissible diagrams. About 16,000 applications (826 flypes, 47 tube
  flips, 14,834 stabilisations, 124 level moves) preserve topology, H1 and hom counts into S3, Q8, A4.
  The "4D Tait flype" is the classical flype applied to a ball containing whole bands, and the tube
  flip exchanges two pictures of the same 1-handle; both are judged genuine isotopies. Merged classes
  have constant invariants at n = 6..9.
- **Topology formulas.** Components, Euler characteristic, orientability and
  e = w(L−) − w(L+) check out on the table.
- **Group.** Three independent Wirtinger implementations (mine and two agents') agree with the repo's
  hom counts on every diagram tried.
- **Other layers.** Layers 1–8 and 10 unchanged on 96 move-equivalent variants (three per surface, 32
  surfaces), and again on 384 variants (twelve per surface, four moves each including level moves)
  **[re-ran]**. In the second run the only layer that changed was `double_triple_linking`, 13 times, on
  `8_1^{1,1}`, `10_1^{1,1}`, `9_1^{1,-2}` and the new class `11_2^{0,1}`. These tests vary the diagram
  by moves. They do not vary the group presentation independently, and so did not catch the
  presentation-dependence of the E_1 layer on mixed-orientability links (§7.3, correction). Mochizuki R_3 values match the literature for `10_2`, `8_1`, `9_1`, `10_3`. The cover
  Betti formula (Thm 5.5) and its peripheral subgroup were first checked by reading only; they have since
  been run through the presentation test of §7.3 (186 changes of generators, unchanged) **[re-ran]**.
- **Mirrors.** The enumeration keeps one diagram from each orbit of "switch all crossings" and "switch
  all markers", so the key has to be unchanged by both. It is: 96 comparisons over the 32 surfaces, no
  layer changed, on the original code and on Carbon's tree of 09:46
  (`audit_scripts/mirror_invariance.py`) **[re-ran]**.
- **Admissibility of stored data.** All 249 distinct diagrams, by SnapPy **[re-ran]**.

### 7.7 Not reached

- Ω8 case B (`moves.py:771-846`) was not verified by hand. (Random walks were run afterwards
  **[re-ran]**: 851,040 move applications in walks from 31 surfaces, each checked for planarity,
  admissibility, topology, H1 and hom counts, with zero violations; this includes 42,103 applications
  of Ω7, 583 of Ω8, 254,437 flypes and 95,109 tube flips.)
- (Done afterwards **[re-ran]**: the 79 tube prunes at ch=11, §7.8; certification over all ch=11 candidates,
  `audit_A/t3_full_candidates.py` in twelve shards: 156,353,366 candidates, zero failures. So the silent
  drop of §7.2 item 1 did not occur at ch=11. A clean re-run of the whole ch=11 pipeline from commit
  3089bb0 was also started; see §7.8.)
- False-negative rate of the unlink certifier on hard unknot diagrams with 8–12 crossings.
- A full invariance test at 20 or more variants per surface with level moves
  (`audit_C/invtest.py`; twelve per surface were run on the original code and sixteen on Carbon's tree, and mirrors
  separately, §7.6); the Q6 cocycle basis and the "0 + 24u" value; triple linking on any non-zero
  example; `alexander_biquandle_invariants` and `twisted_s3_invariants`.
- (Done afterwards with GAP **[re-ran]**: the groups of the entries that can have a finite group, §6.1.)
- Whether braided form loses any tangle diagram with at most six crossings at b=4; the half-twist and
  T_i = T_j prunings in the tri-plane census (read only).
- What the two ch=12 spheres with E1 = (3, x²+1) are, and whether they are one knot. Carbon's note of
  08:01 argues that they are not ribbon (the ideal is not principal, so the group has no presentation
  of deficiency one) and not k-twist-spins for k <= 7 (from the orders of G/<<mu^k>>; see §10). I had not
  checked either; at 10:10 I did **[re-ran]**, and both hold:
  - the ideal, recomputed by a Fox calculus written from scratch on the stored presentations
    (`audit_scripts/f9_fox.py`), is (3, x²+1) for both diagrams, with 3 exhibited in it. Its quotient is
    the field of nine elements, finite and non-zero, so it is not principal;
  - GAP on my own Wirtinger presentation (`audit_scripts/f9_knot.py`) gives H_1 = Z + Z_3 + Z_3 for the
    4- and 8-fold cyclic covers of the complement and Z for the 2-, 3-, 5-, 6-fold ones, and
    |G/<<mu^k>>| = 2, 3, 108, 720, 13824 for k = 2, 3, 4, 6, 8, the values Carbon quotes. Identical for
    the two diagrams;
  - the quotient of order 108 has derived subgroup of order 27 and second derived subgroup of order 3,
    and the quotient of order 720 contains SL(2,5). So the commutator subgroup of G is not abelian and
    G is not soluble.
  - **the two diagrams have the same group, by a map taking a meridian to a meridian** (10:24,
    `audit_scripts/f9_same_group.py`): GAP's simplified presentations of both rewrite, with p the product
    of the two generators, to `< a, p | p^3, a p a = p a p a p >`, where a is a meridian; each relator
    of each presentation reduces to the empty word in the other in at most two rewrites. This does not
    prove that the diagrams are one knot, but nothing separates them.
  Which named 2-knot this is remains open.

---

### 7.8 Clean re-run of ch=11 from the committed code [re-ran]

`scripts/reproduce_datasets.py --ch11-enumerate --ch11-export --data-dir <fresh dir>` at commit 3089bb0,
unmodified, two threads.

| Stage | This run | Documented |
| :--- | ---: | ---: |
| Shadows | 3,378 | 3,378 |
| Raw survivors of the scan | 66,017,754 | 75,125,028 |
| Level components explored / discarded | 22,155,075 / 22,061,181 | 29,916,769 / 29,835,370 |
| Level-irreducible after global merge | 39,185 | 39,184 |
| After pass A | 6,899 | 6,899 |
| After pass B | 712 | 712 |
| After passes C1, C2 | 335, 195 | 335, 195 |
| After the tube filter | 116 | 114 in the table, 116 in the JSON |
| After flype merge and final tunnelling | 105 | 105 |
| New classes (singletons + packets) | 10 (9 + 1 packet of 3) | 10 (9 + 1 packet of 3) |
| Unresolved diagrams, cells | 93, 21 | 93, 21 |

About 75 minutes of wall time. **The result reproduces**: the same ten labels with the same hom counts
and `E_1`, identical representative diagrams for nine of the ten (a different diagram of
`11_2^{0,1}`), the same packet of three for `11_1`, and the same stratification of the 93 unresolved
diagrams (their representatives differ in some cells). Only the scan counters differ from the
documented ones (and agree exactly with the agent's independent re-run), so the first two rows of the
funnel table in `ch11_expansion.md` came from other code or parameters and should be replaced. The
output is in `audit_scripts/repro11/`.

Tube prunes at ch=11 (`audit_scripts/tube11.py` on the 195 pre-filter diagrams): 79 pruned, 116 kept, and
every pruned diagram has the hom counts of Z or Z_2 into S3, Q8, A4, S4, A5.

## 8. Scripts worth keeping

Mine are copied to `audit_scripts/` in the repo root (untracked), with a README giving the command for
each. The agents' scripts are copied to `audit_scripts/agents/audit_A/` … `audit_D/` (their logs and
larger outputs remain under the scratch path given in §7). The ones I would turn into permanent tests:

| Script | Use |
| :--- | :--- |
| `audit_scripts/snappy_unlink.py` | External admissibility check of every stored diagram |
| `audit_scripts/indep_pi1.py` | Independent Wirtinger presentation and hom counter |
| `audit_scripts/tube_check.py` | Blind census at one ch with hom-count profile of each survivor |
| `audit_scripts/handle_cord2.py`, `handle_cord_bs12.py` | 1-handle cord classes for bases `8_1`, `9_1` |
| `audit_scripts/handle_free.py`, `handle_free3.py`, `handle_ab.py`, `handle_free_groupcheck.py` | Names R(w), T(w), K(w): cord words over two or three unknotted spheres; second tube on T(b); group check of each name |
| `audit_scripts/name_ch12.py` | The same naming run over Carbon's partial ch=12 two-saddle survivors |
| `audit_scripts/weak_prime.py`, `weak_prime2.py` | Weak-primeness check of the ch=11 link classes |
| `audit_scripts/tube11.py`, `audit_scripts/repro11/` | Tube-prune audit at ch=11; output of the clean ch=11 re-run |
| `audit_scripts/cs_check2.py`, `cs_check4.py` | Quotient orders and hom counts for the candidate Cappell–Shaneson diagram |
| `audit_scripts/presentation_invariance.py` | Invariance of each key layer under random changes of generators of the group presentation (this is the test that catches the E_1 defect of §7.3) |
| `audit_scripts/groebner_presentation_dep.py` | The specific change of generators that swaps the E_1 data of `11_2^{0,-2}` and `11_3^{0,-2}` |
| `audit_scripts/verifier_adversarial.py` | False path, near misses and genuine moves through Carbon's move-path verifier (§12) |
| `audit_scripts/f9_same_group.py` | Proof by rewriting that the two (3, x²+1) spheres have the same group |
| `audit_scripts/f9_fox.py`, `f9_knot.py` | Independent Alexander ideal and cover homology of the two ch=12 spheres with ideal (3, x²+1) (§7.7) |
| `audit_scripts/finite_groups.py` | Order and structure of the groups of the all-non-orientable entries, by GAP (§6.1) |
| `audit_scripts/mirror_invariance.py` | The key under the two mirror operations the enumeration quotients by |
| `audit_scripts/dlk_many_markers.py` | Double linking under marked-kink stabilisation, to reach the many-marker branch (§12) |
| `audit_scripts/packet11_iso_verify.py`, `packet11_iso.py` | Proof by rewriting that the two cord words of the `11_1` packet have isomorphic groups |
| `audit_scripts/gap/` | Knuth–Bendix proofs that two ch=12 groups are Z ⋉ Z^3 (§10) |
| `audit_scripts/cs_hand_proof_check.py` | Checks the starting relations of the third party's hand proof of the same statement (§11) |
| `audit_C/t4.py`, `audit_C/invtest.py` | Invariance of each key layer under random move-equivalent variants (this is the test that catches Dlk) |
| `audit_C/mydlk.py` | Candidate corrected double-linking form |
| `audit_A/t1_classical_moves.py`, `t4_yoshikawa_moves.py` | Property tests of Reidemeister and Yoshikawa moves |
| `audit_A/t2_census_rejections.py`, `audit_B/certcount.py` | Count certification failures inside the pipeline |
| `audit_B/t_mrange.py` | Scan the marker ranges the default census skips |
| `audit_B/ch11_scan.py` | Scan stage of ch=11 with the script's own partition (gives 66,017,754) |
| `audit_D/run_variant.py` | Tri-plane census with arch pruning switched off |
| `audit_D/audit_meander.py` | Exhaustive check of the meander statement |
| `audit_D/audit_amb2.py` | Shows the two "unknot-matched" ch=12 spheres are knotted |

---

## 9. How to make this better

Ordered by value for effort. Items 1–4 are about making what exists true and checkable; 5–9 are
mathematical extensions; 10–14 are engineering.

### Make the existing claims true

**1. Weak primeness: done for ch=11 (§6.5); build it into the pipeline for ch=12.** The reference
catalogue must contain all sums of two surfaces of ch <= N-1 with the topology of a candidate, built
from minimal diagrams (the Klein bottle with e = 0 has a 3-vertex diagram that is not a planar sum of
table diagrams), and the key must include sublink data: `11_2^{0,1}` and `11_1^{0,-2}` have the same
group and `E_1` as composites of table entries and differ from them only in their components.

**2. Resolve residuals with 1-handle calculus instead of deeper search (days).** A marked bigon in a
diagram exhibits the surface as (base) + (1-handle). By Boyle (Trans. AMS 306, 1988) the handle on a 2-knot
is determined by a double coset of the meridian subgroup, and that coset can be computed exactly from the
diagram when the base group has a faithful matrix representation. The prototype in `audit_scripts/` does
this for bases `8_1` (B_3 in SL(2,Z) x Z) and `9_1` (BS(1,2) in the affine group of Z[1/2]). It already:

- proves ch=11 residual #89 is `10_1^1` and residuals #31 and #56 are unknotted tori;
- shows residuals #7, #55, #79 are one and the same surface, `9_1` + non-orientable h(a^3) =
  `8_1` + non-orientable h(Delta^2), the Klein-bottle analogue of `Tube(3_1)`;
- merges the three ch=10 diagrams of `10_1^1` that the move search leaves apart;
- gives `11_1^{-2}` = `8_1` + non-orientable h(b^2).

The free-group case (base = two unknotted spheres) is done too and gives names to sixteen diagrams
(§6.4), a conjecture `ch(K(w)) = ch(T(w)) + 1` confirmed in six cases, and a way to label ch=12 survivors
as they come out (`handle_free.name_surface`).

To extend: base `10_1` (figure-eight group in SL(2, Z[omega])); bases of positive genus (larger peripheral
subgroup; Kamada, "Cords and 1-handles attached to surface-knots", 2014); multi-component bases. Before
this goes into a paper, re-read the hypotheses of Boyle's theorem and Kamada's treatment of non-orientable
handles; I have relied on them from memory.

**3. Turn the open question of §2.2 into a target.** `tau^0(3_1) # P` versus `tau^2(3_1) # P` is open in
print. Diagrams of both exist in the repo (`8_1 # P` at 10 vertices, `10_2 # P` at 12) and they share an
invariant cell. A bidirectional move search between them, with the path stored, would be a publishable
result if it succeeds. If it fails, nothing is lost. The same goes for "`8_1` + non-orientable
h(longitude) = `8_1 # K^2`?" from item 2.

**4. Reframe the RP^2 campaign (§2.3).** Two well-posed halves: (a) for every chi = 1 survivor up to ch N,
an explicit move path to a 2-edge-cut diagram with a P summand, which proves "every projective plane with
ch <= N is reducible"; (b) flag any chi = 1 survivor whose meridian has order greater than 2 in some finite
quotient, which certifies irreducibility by the argument of Yoshikawa's Prop 5.3 and of
Hughes–Kim–Miller–Nahm. Also well-posed and open: the ch-index of the HKMN projective plane.

### Mathematical extensions

**5. Cross-match with the ribbon 2-knot tables.** `11_1`, `11_2` and the ch=12 2-knots are `R(w)`, spins
and twist-spins. Kanenobu–Takahashi (2021) classify 1-fusion ribbon 2-knots of length up to six. Matching
gives each census 2-knot a standard name and gives that literature a new column (ch-index).

**6. Split `tau^2(4_1)` from `tau^2(5_1)`.** The ch=12 cell `12_m2_S2_3` (group Z_5 ⋊ Z, six unmerged
diagrams) is described as "`tau^2(4_1)` / `tau^2(5_1)`". These are different 2-knots (fibres L(5,2) and
L(5,1)) with the same group. The Mochizuki 3-cocycle of the dihedral quandle R_5 distinguishes twist-spins
of 2-bridge knots (Asami–Satoh); the repo has it only for R_3. Implementing R_p for p = 5, 7 is small and
also gives chirality and triple-point lower bounds.

**7. Exact ideal comparison.** Replace "Hermite normal form on bounded monomial shifts plus evaluation
fingerprints" for `E_k` over `Z[x^{±1}, y]/(y^2-1)` by Gröbner bases over Z (Singular, Macaulay2 or Sage).
That makes equality of ideals a decision procedure rather than a heuristic.

**8. Identify `11_1^{-2}` geometrically.** Group 2O, meridian order 8, spun trefoil plus one non-orientable
handle. Compare with Yoshikawa's irreducible Klein bottles (Proc. AMS 126, 1998) and test Carbon's
deform-spinning conjecture by compiling a diagram and matching it by moves.

**9. Use the tri-plane side for what it is good at.** Bridge-number upper bounds for the new classes (the
tables currently say "4–9"), and a statement of exactly which diagram universe the b <= 4 census covers.
An unconditional b = 4, c_tri <= 6 census (all trivial tangles, not only braided words) would be a real
contribution to Conjectures 3.6 and 3.7 of GLMO, though still only for b <= 4.

### Engineering

**10. Certificates, not assertions.** Store, for every discard after the level-irreducible stage (39,184
diagrams at ch=11) and every `MOVE_CERTIFIED` identification, the move path; write a checker of a few
hundred lines that validates a path move by move and shares no code with the search. Then the census is
checkable without trusting the 4,400-line C kernel.

**11. Tests that can fail.** A blind ch <= 9 recovery on every run (12 s), ch = 10 marked slow; a
property test that random sequences of every move type preserve every key layer; **a change-of-generators
test for every layer computed from a group presentation** (the test that exposes the E_1 defect of §7.3;
my move-based tests did not; `audit_scripts/presentation_invariance.py` now runs it on every layer); a
test that the key is unchanged by the two mirror operations (`mirror_invariance.py`); a
differential test of C moves against Python moves; an external admissibility check (SnapPy) on every
stored diagram; for any move-path verifier, known non-equivalent pairs that it must reject
(`Tube(3_1)` vs `8_1 # T^2` is the first). Delete or rewrite the tests that assert stored numbers.

**12. Reproducibility.** One command per dataset, including the tri-plane and ch=12 ones; provenance
(commit, parameters, wall time, core count) inside each JSON; `--python 3.11` in the bootstrap line; pinned
dependencies; bundle Tailwind into the atlas.

**13. A computer algebra system for groups.** GAP (conda-forge `gap-defaults`) did in seconds what the
sympy routines could not do at all: orders and structure of the finite groups (§6.1), Knuth–Bendix
proofs of the two Z ⋊ Z^3 groups (§10), homology of cyclic covers and quotients by powers of a meridian
(§7.7). Add it to the environment and use it for the group layer of any certificate.

**14. No sentence ahead of its artefact.** At 10:13 `complexity_measures.md` said a collapse was "verified
diagrammatically ... in `move_search_results.json`", a file that did not exist. Generate the documents'
numbers and status words from the data files, or check in a script that fails when a cited file is
missing.

### Exposition

Separate three kinds of statement everywhere: theorems from the literature (with a correct citation),
computations (with the assumptions they rest on), and conjectures. Remove exclamation marks and the words
"certified", "unconditional" and "proves" wherever the thing certified is an invariant match. State once,
in one place, the two assumptions the whole census rests on: Yoshikawa's table is complete, and the
implemented moves are isotopies. `expo/report.tex` is a draft of that plainer account.

---

## 10. After the audit: what changed on Carbon's branch, and what I re-tested (as of 08:45 EST)

Everything above describes `main` at `3089bb0`. During the night Carbon worked on the branch
`exp/epyc64-campaigns` (uncommitted working tree when I looked) and reports having fixed the items in §7.
I copied that working tree to a scratch snapshot at about 08:25 and ran the audit's tests against it.
I did not review the new code line by line.

| Item | Test | Result on the snapshot |
| :--- | :--- | :--- |
| Double linking (§7.1) | my two counterexamples; 600 move-equivalent variants (the test the old code failed 74 times) | fixed; 600/600 unchanged |
| All ten key layers | 512 variants (16 per surface, 32 surfaces, four moves including level moves) | no layer changed |
| Errata data (§2.1) | read the diff | `10_1^{0,1}` no longer flagged; printed value corrected; test renamed to two errata |
| Tube rule (§7.5) | condition evaluated on 204 marked bigons in 405 diagrams | rejects no genuine handle |
| Test suite | `pytest` | 28 pass, 1 fails (`test_atlas`: stored HTML is one character shorter than the regenerated one) |

Not re-tested: the new Gröbner-basis ideal comparison, the tri-state unlink certificate, the sublink
profile in the key, the full marker range, and the ch=12 census as a whole. (Later the same morning I
did test the ideal comparison under changes of presentation, §12, and checked the specific ch=12 claims
below.)

Carbon's ch=12 claims, and what I could check (`audit_scripts/cs_check2.py`, `cs_check4.py`):

- **A ch=12 two-saddle diagram with the Alexander polynomial of `11_2` but a different group.**
  Confirmed. `X[1,2,3,4], X[5,6,4,7], X[6,8,9,1], Y[9,10,11,2], X[12,13,14,3], Y[14,15,16,7],
  X[5,17,18,8], X[17,19,10,18], X[11,20,21,12], X[13,21,22,15], X[23,24,19,16], X[24,23,22,20]`
  has |G/<<mu^3>>| = 21 where `11_2` has 168. (It was already in `targeted_ch12_13_census.json`, filed as
  matching `11_2`.)
- **Its group is Z ⋉ Z^3 with characteristic polynomial t^3-t^2+2t-1: proved** (computer-assisted,
  09:00 EST). Evidence first: hom counts into S5, A6, PSL(2,7), AGammaL(1,8) and AGL(1,p) for p = 7, 11,
  13, 29, 43 equal those of an explicit presentation of that group and differ from `11_2` (672 against
  1176 into AGammaL(1,8)). Then the proof: Knuth–Bendix completion in GAP 4.15.1 (installed in a private
  conda environment under my scratch directory). With x a meridian, u = (second meridian)·x^-1,
  v = x u x^-1, w = x v x^-1 and a wreath-product ordering, the relations complete to a finite confluent
  rewriting system (32 rules from the repo's stored 4-relator presentation; 80 rules from my own
  24-generator Wirtinger presentation of the diagram, in which every edge generator reduces to
  x^{±1} times a word in u, v, w). In it u, v, w commute, x w x^-1 = w v^-2 u (resp. w^2 v^-1 u for the
  other meridian orientation), and the only other rules are inverses and the conjugation rules, so
  normal forms are x^k u^a v^b w^c. Scripts and logs: `audit_scripts/gap/`. The diagram is admissible by
  SnapPy and is a sphere with two minima, two saddles and two maxima.
  Since det(A - I) = ±1 this is a Cappell–Shaneson group; t^3-2t^2+t-1 (the reciprocal of the polynomial above, i.e. the
  other meridian orientation) is x^3-nx^2+(n-1)x-1 with n = 2, the characteristic polynomial of a
  Cappell–Shaneson matrix of trace 2 in the usual normalisation (checked at 10:30 against the statement
  quoted from arXiv:2404.05096 in a search result; I did not read the paper). The second diagram's
  t^3-3t^2+2t-1 is the trace-3 case. (Added 10:42: a search result summarises Theorem 16.11 of Hillman's *Four-manifolds,
  geometries and knots* as: the group Z^3 ⋊ Z of a Cappell–Shaneson knot has deficiency 0 exactly when it
  is 2-generated. These two groups are 2-generated, so by that statement they have deficiency 0 and the
  knots are not ribbon. I did not read the theorem itself.) Consequences that follow from the group (deficiency zero, hence not ribbon; Hillman's
  classification up to Gluck reconstruction) should be cited from the sources before they are written
  as theorems. The discovery is Carbon's; the proof of the group is mine.
- **The same recipe over all 29 stored ch=12 two-saddle spheres** (`audit_scripts/gap/cs_all.py`, log
  alongside; independent Wirtinger presentation, GAP Knuth–Bendix, 150 s limit each; 28 of 29 complete; the
  one that does not, `12_m2_S2_0` with Alexander polynomial x^4-x^3+x^2-x+1, cannot have a group Z ⋉ Z^3,
  whose Alexander polynomial is a cubic):
  - exactly two diagrams have Z ⋉ Z^3: the one above, and the fifth stored diagram of the class with
    `E_1 = (x^3-2x^2+3x-1)`,
    `X[1,2,3,4], X[4,5,6,1], X[7,8,9,6], Y[9,10,11,2], X[3,12,13,14], Y[14,15,16,17], X[18,7,5,17],
    X[8,18,19,10], X[20,21,12,11], X[21,22,15,13], X[23,24,19,16], X[24,23,22,20]`,
    with x w x^-1 = w^3 v^-2 u, i.e. characteristic polynomial t^3-3t^2+2t-1 (reciprocal of
    x^3-2x^2+3x-1): the next polynomial of the Cappell–Shaneson family;
  - the other four diagrams stored in that class have non-commuting u, v, w, so **the class
    "`12_m2_S2_5`" of `complexity_measures.md` (five move components) is two different knots**, a ribbon
    knot and a Cappell–Shaneson-type knot with the same Alexander polynomial. The same happens with
    `11_2`: one of the three ch=12 diagrams filed as "matches `11_2`" is the first Cappell–Shaneson-type
    knot, the other two are `11_2`;
  - all six diagrams of the class with `E_1 = (5, x+1)` have group Z_5 ⋊ Z and both with `(7, x+1)` have
    Z_7 ⋊ Z, as expected for 2-twist-spins; the group cannot separate `tau^2(4_1)` from `tau^2(5_1)`.
- **The (3, x^2+1) sphere is not a k-twist-spin for k <= 7.** Carbon's claim, which I first judged
  only partly supported and have since recomputed (§7.7; `audit_scripts/f9_knot.py`). If the knot were
  a k-twist-spin, mu^k would be central, G/<<mu^k>> would be G' ⋊ Z_k, and its derived subgroup would
  be G', which must map onto the Alexander module Z_3 + Z_3. GAP gives |G/<<mu^k>>| = 2, 3, 108, 5,
  720, 2352, 13824 for k = 2..8 (Carbon's values), with derived subgroup trivial for k = 2, 3, 5,
  SL(2,5) for k = 6 and SL(2,7) for k = 7: none maps onto Z_3 + Z_3. For k = 4 the derived subgroup has
  order 27 and passes that test; there Carbon's second argument applies (a central mu^4 would force
  |G/<<mu^8>>| <= 2·108, and it is 13824). For k = 8 the derived subgroup of the quotient has order 1728 and
  abelianisation Z_3 + Z_3 (GAP, 10:54), so the first test does not exclude an 8-twist-spin, and the
  doubling test would need the quotient by mu^16, which I did not compute. k >= 8 and the
  twist-roll-spins claimed in the note are therefore not covered.

Carbon's notes are in `carbon.md`; my replies, with the detail behind each line above, are in `fable.md`.

---

## 11. The third-party audit of 09:05 EST, and what I checked of it

A separate read-only audit (`third_party_audit.md` under
`~/.gemini/antigravity-cli/brain/1dadd8f0-.../`, with scripts in its `scratch/`) looked at Carbon's live
working tree and at my files. I read it and its scripts; I did not run its scripts. Status of each point:

**On Carbon's working tree**

| Finding | My check |
| :--- | :--- |
| `verify_move_path_independent` (`move_search.py`) accepts a one-step "path" between `Tube(3_1)` and `8_1 # T^2`, because cut darts are all replaced by the same token | **Reproduced** on my 08:25 snapshot of Carbon's tree with the two diagrams given: `is_valid=True, step_kinds=('OMEGA_3_4_7',)`. The two diagrams are different surfaces: peripheral SL(2,3) counts 24 and 0, and my cord-class computation gives the longitude class for the first and the trivial handle for the second. **No move-path certificate validated only by that version of the function can be trusted.** Carbon has since rewritten it; the rewritten version rejects this path (§12). |
| `canonical_zh1_groebner_key` depends on the presentation (orientable meridians are not constrained mod 2), swapping the keys of `11_2^{0,-2}` and `11_3^{0,-2}` under a Tietze change; and the Z-Buchberger can hang | **Reproduced** with my own script (`audit_scripts/groebner_presentation_dep.py`) on the 08:25 snapshot, **and on the original code at `3089bb0`**: the defect predates Carbon's change and was missed by my audit (§7.3, corrected). My invariance tests vary the diagram by moves and do not exercise it; "no layer changed" is necessary, not sufficient. The hang was not re-run. Carbon's tree of 09:46 passes my presentation test (§12). |
| Funnel counters at ch=12 are recorded before unlink certification | Not re-run; consistent with §7.2 and with what I saw at ch=11 (where certification of all candidates found no failure). Final survivors are certified. |
| New `double_linking_mod2` drops a pair silently when a component has more than 12 markers | Not re-run on that version. Carbon's tree of 09:46 has a separate branch for this case, which passes a direct test (§12). |

**On my files**

| Finding | Action |
| :--- | :--- |
| `expo/report.tex` states the handle classification without the connected-base condition at the point where it is first used, and does not say that the trivial twisted tube gives the Klein bottle with e = 0 | Fixed: the sentence in §5.2 now says "when the base is connected", points to the caveat, and states the e = 0 fact. (The caveat for link bases, citing Kamada arXiv:2003.12977, was already in §5.3.) |
| `\bibitem{K14}` without pages; `\bibitem{Kan}` without volume, year, pages | K14: pages 595–609 restored (two sources now agree). Kan: I rendered the first page of the repository PDF; it shows title, author and 14 pages but no volume or year, so the entry says so instead of guessing. |
| `audit_scripts/handle_free.py::bigon_bases` does not check marker-type compatibility | Added. All sixteen names and the three-sphere and two-tube readings are unchanged. (For the record, the other parity with chi + 2 is also a tube on the same base along the same cord, attached below the saddle level; such diagrams reduce by Ω6 and none occurs in any data set here.) |

**New mathematics in it, and how it relates to mine**

- A hand proof that the commutator subgroup of the Cappell–Shaneson candidate is Z^3, plus an argument
  that such a group has no deficiency-one presentation. **Checked line by line at 10:40**
  (`audit_scripts/cs_hand_proof_check.py` and by hand). Their `idx=191` is the same diagram as the first
  one of §10 (I had taken it for a different one: the diagram string printed in their report is another
  rendering, which parses as a Klein bottle in the original code and is rejected as non-planar by
  Carbon's current parser; the index in Carbon's survivor file points to my diagram). Rewriting the
  repository's four relators with g1 = t, g2 = t v_{-1} gives exactly their starting relations: R0 is
  (E1), R1 is v3 = v1^-1 v0 v2 v1^-1, and R2 shifted down, with v_{-1} = v0 v2 v0 v1^-1 from R1, is (E2).
  Their steps 3 and 4 then follow as written, so v0, v1, v2 commute, and conjugation by t acts as in my
  rewriting system (x w x^-1 = w v^-2 u). So the two proofs are of the same statement on the same
  diagram by different methods (theirs by hand from the repository's presentation, mine by Knuth–Bendix
  from that presentation and from an independent one). I did not check their deficiency argument.
- Explicit mutually inverse homomorphisms between the groups of `11_1.b` = R(b a b^-1 a^-1) and
  `11_1.c` = R(b a b^-1 a), sending the meridian to its inverse. This is a group isomorphism, consistent
  with my reading of the ribbon-knot classification (§6.4: the two are mirror images, hence one
  unoriented surface). **Verified at 09:42** (`audit_scripts/packet11_iso_verify.py`): each of the six
  identities (the two relators map to 1, and both composites fix a and b) reduces to the empty word in
  at most three replacements of a subword by the complementary part of the relator, so the maps are
  mutually inverse isomorphisms. My own search (`packet11_iso.py`) confirms that no automorphism of the
  free group of the simple kind (a to a^±1, b to a conjugate of b^±1 by a word of length at most 6)
  does this, as the third party says. Isomorphic groups alone do not make the surfaces equal; the
  classification argument does, for these knots.
- It identifies the polynomial family as Delta_n(t) = t^3 - (n+1)t^2 + (n+2)t - 1, with n = 0 for this knot
  and n = 1 for Carbon's second candidate, matching what I recalled (Gompf's A_0) up to t -> 1/t.

---

## 12. Second re-test of Carbon's working tree (snapshot of 09:46 EST)

Carbon kept editing after §10 and §11 were written (`alexander.py`, `linking.py`, `move_search.py`
and others changed between 08:59 and 09:32). I took a new snapshot and ran the tests again.

| Item | Test | Result on the 09:46 snapshot |
| :--- | :--- | :--- |
| Move invariance of the key | 512 variants, as in §10 | no layer changed |
| Presentation invariance of the key (§7.3, §11) | `audit_scripts/presentation_invariance.py`, 640 variants (20 random changes of generators for each of 32 surfaces) | **no layer changed, no assertion raised.** The E_1 defect that the third-party audit found is fixed in this tree, as far as this test can see. |
| The false move path of §11 | `audit_scripts/verifier_adversarial.py` | **now rejected**, by the invariant check and by the step classifier on its own |
| Near misses for the step classifier | same script: one crossing or one marker flipped after a genuine move; 2,488 flips, 164 admissible | none accepted with a changed key (each accepted one is literally an output of a move generator, or a symmetry) |
| Genuine moves through the step classifier | same script and `verifier_o78.py`, moves from the repo's own generators, on the table diagrams and along random walks | all accepted: Ω1/Ω6 64/64, Ω2 64/64, Ω3/Ω4 27/27 and 60/60, Ω5 15/15 and 60/60, Ω7 150/150, Ω8 60/60, flype 52/52, tube flip 23/23 |
| Double linking with more than 12 markers on a torus (§11) | `audit_scripts/dlk_many_markers.py`: marked kinks added to each of the 11 links with a torus component until the diagram has 40 markers; 110 diagrams | value unchanged in all 110, no pair lost; the new branch for more than 12 markers ran 131 times. (For four links, `9_1^{0,1}`, `10_1^{0,1}`, `10_2^{0,1}`, `11_1^{1,-2}`, the function returns no value at all, before and after.) |
| Test suite | `pytest` | 28 pass, 1 fails (`test_atlas`: `spoke_total` is 41, the test wants at least 44) |

Notes.

1. The rewritten `verify_move_path_independent` accepts a step only if the second diagram is, up to
   symmetry, an output of one of the repo's move generators applied to the first, and a port-matching
   check passes. So it is as sound as the generators, which §7.7 records as tested on the original code (851,040
   applications, no violation) and which the 512-variant test above exercises on this tree. It is not an
   independent implementation of the moves.
2. A snapshot I took fifteen minutes earlier rejected every Ω3/Ω4 step on the table diagrams (the port
   check counted the triangle's neighbours as modified). Carbon had fixed that by the time I reported
   it; the table above is for the fixed version.
3. With `allow_o4=True` (the default) a mirror image counts as a step. A path verified that way shows
   equality up to mirror image, which should be said wherever such a certificate is quoted.

I did not review the new Alexander-ideal code or the new verifier by reading them in full; the
statements above are test results.

Tests (snapshot of 10:43). Carbon added `code/tests/test_soundness_and_ch12.py` with seven tests that
follow §9 item 11: blind recovery at ch <= 9 and at ch = 10, invariance of the ideal key under changes
of presentation, rejection of the false move path, the double-linking counterexamples and the
many-marker case, the `11_1` packet isomorphism and the group 2O, and the two Cappell–Shaneson
diagrams. On my copy the suite gives 35 passed, 1 failed in 112 s; the failure is still `test_atlas`,
now on `total_distinct_classes >= 95`. I read the names and docstrings of the new tests, and the bodies
of two. One remark: the verifier test rejects `10_1^1` against a planar sum `8_1 # T^2`, two diagrams
that are far apart, and its near-miss block sits behind two `if`s that may not execute. The pair that
actually fooled the old verifier (two 11-vertex diagrams differing on three vertices,
`audit_scripts/verifier_adversarial.py`, first lines) is not in the test and should be.

`HANDOFF.md` (edited by Carbon at about 10:45) is headed "What Is Verified & Committed". At that time
nothing on the branch was committed: `git log` ended at `3089bb0` and the changes were in the working
tree.

Last regression (snapshot of 10:46, all five scripts): presentation changes 256/256 unchanged; moves
512/512 unchanged; mirrors 96/96 unchanged; false path rejected, all genuine moves accepted, no near
miss accepted with a changed key; double linking with many markers 66/66 unchanged.

Still stale in the documents at 10:51 (a text search, not a full read): the root `README.md` is unchanged
from the commit and still quotes 36,166,398 tri-plane triples and the meander obstruction "forcing"
`ch(10_2) = 10`; `HANDOFF.md` has "zero uncertified unlink heuristics" as a rule, which is fine as a rule; `math/complexity_measures.md` still has "7
certified new prime" ch=12 classes and the citation of the missing `move_search_results.json`.

Diagram strings. Carbon changed `MarkedDiagram.from_epd` / `to_epd`. I checked that 60 strings from the
committed data (the 23 table entries, the ten ch=11 entries, 25 unresolved diagrams and the two
Cappell–Shaneson diagrams quoted in §10) are read as the same diagram by the original parser and by
Carbon's (same vertices, component types, normal Euler numbers and hom counts into S3 and A4;
`audit_scripts/epd_compat.py`). So the diagrams quoted in this file and in the report mean the same
thing in both trees. The string printed in the third-party report for its `idx=191` is the exception
(§11).

A move path for Viro's theorem. `complexity_measures.md` cites a path from `10_3 # P` to `P` as verified.
No such result exists in the tree. I ran Carbon's own bidirectional search for it on my copy for fifteen
minutes (24 workers, beam 4000, at most 15 vertices, eight starting diagrams): it did not connect in
that time. That says nothing about whether a path exists; it says the sentence is not yet earned.

Documents (checked at 10:20, uncommitted): `REFS.md` now has the six corrected arXiv identifiers and none
of the old ones; `math/yoshikawa_table.md` and `foundations_and_invariants.md` speak of two errata;
`complexity_measures.md` restricts Viro's theorem to `10_3 # P` and states the general question as open;
`REFS.md` and `literature_review.md` cite the Kinoshita counterexample; `ch11_expansion.md` gives the
group of `11_1^{-2}` as 2O. Not yet right: `complexity_measures.md` cites `move_search_results.json`,
which does not exist, as verification of a move path (§9 item 14), and §3.1 of the same file still
announces "7 certified new prime" ch=12 classes. I did not re-read these documents in full.
