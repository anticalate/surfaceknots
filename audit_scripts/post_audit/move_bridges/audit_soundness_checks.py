import json
from pathlib import Path
from collections import Counter
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.yoshikawa_data import YOSHIKAWA_TABLE, YOSHIKAWA_BY_LABEL, compute_full_invariant_record, invariant_discrimination_key
from surfaceknots.ch11_census import planar_connected_sum, _oriented_diagram_variants, _is_admissible

def _component_rep_darts(d: MarkedDiagram) -> list[int]:
    dc, _ = d._dart_surface_components()
    rep: dict[int, int] = {}
    for hd in range(4 * d.num_vertices):
        rep.setdefault(dc[hd], hd)
    return [rep[k] for k in sorted(rep)]
from surfaceknots.moves import (
    canonical_marked_key,
    list_level_yoshikawa_neighbors,
    list_4d_tait_flypes,
    _list_omega3_and_4_neighbors,
    _list_omega5_neighbors,
    _list_omega7_neighbors,
    _list_omega8_neighbors,
    _list_marked_tube_flip_neighbors,
)

DATA_DIR = Path("code/data")
ch11_data = json.load((DATA_DIR / "expanded_ch11.json").open())
ch12_data = json.load((DATA_DIR / "expanded_ch12.json").open())

print("1. Checking indecomposability status of all 10 ch=11 and 34 ch=12 prime entries:")
for ch, entries in [(11, ch11_data["entries"]), (12, ch12_data["entries"])]:
    ind_counts = Counter(str(e.get("indecomposability")) for e in entries)
    print(f"  ch={ch} ({len(entries)} entries):")
    for k, v in ind_counts.most_common():
        print(f"    {v:2d}x {k}")

print("\n2. Checking if planar_connected_sum EVER fails _is_admissible on prime atoms up to naive_ch=16:")
atoms = []
for e in YOSHIKAWA_TABLE:
    if e.ch_index > 0:
        atoms.append((e.label, MarkedDiagram.from_epd(e.canonical_epd)))
for e in ch11_data["entries"]:
    atoms.append((e["label"], MarkedDiagram.from_epd(e["epd"])))

inadmissible_cnt = 0
total_splices = 0
for i, (la, da) in enumerate(atoms):
    rep_a = _component_rep_darts(da)
    for j in range(i, len(atoms)):
        lb, db = atoms[j]
        if da.num_vertices + db.num_vertices > 16:
            continue
        for lb_v, db_v in _oriented_diagram_variants(lb, db):
            rep_b = _component_rep_darts(db_v)
            for hda in rep_a:
                for hdb in rep_b:
                    for fl in (False, True):
                        total_splices += 1
                        cs = planar_connected_sum(da, hda, db_v, hdb, flip_splice=fl)
                        if not _is_admissible(cs):
                            inadmissible_cnt += 1
                            print(f"  WARNING: inadmissible connected sum: {la}(hd={hda}) # {lb_v}(hd={hdb}) fl={fl} (v={cs.num_vertices})")
print(f"  Tested {total_splices} connected-sum splices up to naive_ch=16: inadmissible={inadmissible_cnt}")

print("\n3. Checking if varying dart choice WITHIN the same connected component ever changes invariant_discrimination_key:")
# Test on non-invertible / chiral / multi-component examples: 10_2, 11_2, 9_1^0,1, 10_1^0,-2
test_labels = ["8_1", "10_2", "9_1^0,1", "10_1^0,-2", "2_1^1", "2_1^-1"]
test_diags = {lbl: MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL[lbl].canonical_epd) for lbl in test_labels}
test_diags["11_2"] = MarkedDiagram.from_epd(next(e["epd"] for e in ch11_data["entries"] if e["label"] == "11_2"))

dart_diff_found = 0
for la in ["10_2", "11_2", "9_1^0,1"]:
    da = test_diags[la]
    dc_a, _ = da._dart_surface_components()
    for lb in ["2_1^1", "2_1^-1", "8_1"]:
        db = test_diags[lb]
        if da.num_vertices + db.num_vertices > 14:
            continue
        # Collect keys using ALL darts on component 0 vs just _component_rep_darts[0]
        keys_rep = set()
        keys_all = set()
        rep0_a = _component_rep_darts(da)[0]
        for lb_v, db_v in _oriented_diagram_variants(lb, db):
            rep0_b = _component_rep_darts(db_v)[0]
            for fl in (False, True):
                cs = planar_connected_sum(da, rep0_a, db_v, rep0_b, flip_splice=fl)
                keys_rep.add(invariant_discrimination_key(compute_full_invariant_record(cs, label="r", include_sublinks=False)))
            # Sample 6 darts on comp 0 of da and 4 darts on comp 0 of db_v
            darts_a0 = [d for d in range(4*da.num_vertices) if dc_a[d] == 0][:6]
            dc_b, _ = db_v._dart_surface_components()
            darts_b0 = [d for d in range(4*db_v.num_vertices) if dc_b[d] == 0][:4]
            for hda in darts_a0:
                for hdb in darts_b0:
                    for fl in (False, True):
                        cs = planar_connected_sum(da, hda, db_v, hdb, flip_splice=fl)
                        keys_all.add(invariant_discrimination_key(compute_full_invariant_record(cs, label="a", include_sublinks=False)))
        if keys_all - keys_rep:
            dart_diff_found += 1
            print(f"  DIFF FOUND for {la} # {lb}: extra keys = {len(keys_all - keys_rep)}")
print(f"  Intra-component dart choice invariance check: diffs={dart_diff_found}")

print("\n4. Checking move symmetry (B in move(A) => A in move(B) up to isomorphism) on sample diagrams:")
sample_epds = [e["epd"] for e in ch11_data["entries"]] + [e["epd"] for e in ch12_data["entries"][:15]] + [r["epd"] for r in ch12_data["residual_collisions"][:25]]
move_fns = [
    ("omega3_4", _list_omega3_and_4_neighbors),
    ("omega5", _list_omega5_neighbors),
    ("omega7", _list_omega7_neighbors),
    ("omega8", _list_omega8_neighbors),
    ("tube_flip", _list_marked_tube_flip_neighbors),
    ("tait_flype", list_4d_tait_flypes),
]
asym_cnt = 0
for epd in sample_epds:
    d = MarkedDiagram.from_epd(epd)
    dk = canonical_marked_key(d, up_to_mirror_and_reverse=False)
    for mname, mfn in move_fns:
        for nb in mfn(d):
            back_keys = {canonical_marked_key(b, up_to_mirror_and_reverse=False) for b in mfn(nb)}
            if dk not in back_keys:
                asym_cnt += 1
                print(f"  ASYMMETRY in {mname} on {epd}!")
print(f"  Move symmetry check across {len(sample_epds)} diagrams: asymmetries={asym_cnt}")
