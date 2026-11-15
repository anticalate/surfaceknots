# audit_scripts (untracked; written by Fable during the 15 Nov 2026 audit)

Independent checks used in `AUDIT.md`. None of them modifies the repo. Run from `code/` with the repo venv
unless noted. `SRC=src`, `CODE=.` below.

| Script | What it checks | How to run |
| :--- | :--- | :--- |
| `snappy_unlink.py` | Both resolutions of every EPD in every JSON are unlinks, using SnapPy/spherogram only (needs `pip install snappy`) | `python snappy_unlink.py SRC CODE` |
| `indep_pi1.py` | Stored `hom_counts` (S3, Q8, A4, S4) against a from-scratch Wirtinger-with-saddles presentation and counter | `python indep_pi1.py SRC CODE` |
| `tube_check.py` | Blind census at one ch, all topologies, no tube pruning; hom-count profile of every survivor | `PYTHONPATH=src python tube_check.py <ch> <workers> <outdir>` |
| `tube_opts.py` | For every marked bigon, both reconnections tried by `has_certified_cyclic_unknotted_tube`; flags a base whose group is not a quotient-source of the diagram's | `PYTHONPATH=src python tube_opts.py` |
| `handle_cord.py`, `handle_cord2.py` | Exact cord class of a 1-handle on the spun trefoil (faithful rep of B_3 in SL(2,Z) x Z) | `python handle_cord2.py SRC CODE <dir with tube_ch*.json>` |
| `handle_cord_bs12.py` | Same for base `9_1` (BS(1,2) in the affine group of Z[1/2]) | `python handle_cord_bs12.py SRC CODE <dir with tube_ch*.json>` |

`tube_ch9.json`, `tube_ch10.json`: survivors of the blind ch=9 and ch=10 runs with their hom counts.
`indep_pi1_lib.py`, `tube_opts_lib.py` are the importable halves of the scripts of the same name.

`agents/audit_{A,B,C,D}/`: scripts written by the four code-audit agents (moves and unlink; C kernel and
enumeration; invariants; tri-plane and complexity). They were written to run from a clone as
`cd code && PYTHONPATH=src .venv/bin/python <script> <args>`, some with `PYTHONPATH=src:<their dir>` for
shared helpers (`lib.py`, `common.py`, `variants.py`). See AUDIT.md §8 for which is which.

Added later the same morning (see AUDIT.md §8 for the one-line purpose of each): `handle_free.py`,
`handle_free3.py`, `handle_ab.py`, `handle_free_groupcheck.py`, `packet11.py`, `name_ch12.py` (+ log),
`weak_prime.py`, `weak_prime2.py`, `tube11.py`, `cs_parse.py`, `cs_check2.py`, `cs_check4.py`, and
`repro11/` (outputs of the clean ch=11 re-run). Most take the arguments `SRC CODE DIR` =
`code/src`, `code`, and this directory.

Added after the third-party audit: `presentation_invariance.py`, `verifier_adversarial.py`, `verifier_o78.py`, `dlk_many_markers.py` and `mirror_invariance.py` (run from
`code/` with `PYTHONPATH=src`; the first takes the number of variants per surface),
`groebner_presentation_dep.py`, `packet11_iso_verify.py` and `packet11_iso.py` (plain Python, no
arguments), and `gap/` (GAP 4.15 scripts and logs for the two Z ⋉ Z^3 groups).

`finite_groups.py SCRATCH` and `f9_knot.py SCRATCH` need GAP and the scratch layout (clone under
`SCRATCH/sk`, GAP under `SCRATCH/gapenv`, `indep_pi1_lib.py` under `SCRATCH/mine`); edit the paths at
their top to run them elsewhere. `f9_fox.py` needs only sympy.
`cs_hand_proof_check.py`, `f9_same_group.py`: run from `code/` with `PYTHONPATH=src` (the first) or with plain Python (the second).
