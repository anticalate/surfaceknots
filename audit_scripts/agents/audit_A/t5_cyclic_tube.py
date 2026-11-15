"""Q2: audit has_certified_cyclic_unknotted_tube on the pre-filter census lists.
For every diagram it prunes, record WHICH 2-gon / which smoothing fired and whether the 2-gon is a 'tube' configuration,
and check whether pi_1 of the pruned diagram is consistent with 'F(cand) + trivial handle' (must be cyclic-like)."""
import json
import pickle
import sys
import time

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from lib import *  # noqa
import surfaceknots.fast_enum as fe
from surfaceknots import MarkedDiagram
from surfaceknots import moves as M
from surfaceknots.invariants.fundamental_group import finite_quotient_profile, wirtinger_presentation
from surfaceknots.unlink import certify_classical_unlink

ch = int(sys.argv[1])
HERE = __file__.rsplit("/", 1)[0]
cache = f"{HERE}/tunneled_ch{ch}.pkl"
try:
    tunneled = pickle.load(open(cache, "rb"))
except FileNotFoundError:
    t0 = time.time()
    tunneled, stats = fe.enumerate_prime_ch_census(ch, num_workers=2, filter_cyclic_tubes=False)
    print("census (no tube filter)", stats, f"{time.time()-t0:.0f}s", flush=True)
    pickle.dump(tunneled, open(cache, "wb"))


def instrumented(d, _depth=0):
    """Verbatim logic of moves.has_certified_cyclic_unknotted_tube, returning (verdict, info)."""
    n = d.num_vertices
    if n < 5 or _depth > 3:
        return False, None
    st = d.surface_topology()
    if st.num_components != 1:
        return False, None
    for u in range(n):
        if d.vtype[u] not in ("Y", "Z"):
            continue
        for su in range(4):
            d_a = d.opp[4 * u + su]
            d_b = d.opp[4 * u + ((su + 1) & 3)]
            w = d_a >> 2
            if w <= u or (d_b >> 2) != w or d.vtype[w] not in ("Y", "Z"):
                continue
            sw_a = d_a & 3
            sw_b = d_b & 3
            if ((sw_b + 1) & 3) != sw_a:
                continue
            tube_cfg = (d.vtype[u] == d.vtype[w]) == ((su & 1) == (sw_a & 1))
            p0, p1 = 4 * u + ((su + 2) & 3), 4 * u + ((su + 3) & 3)
            p2, p3 = 4 * w + ((sw_b + 2) & 3), 4 * w + ((sw_b + 3) & 3)
            remap = [v - (1 if v > u else 0) - (1 if v > w else 0) for v in range(n)]
            m = n - 2
            nvt = tuple(d.vtype[v] for v in range(n) if v not in (u, w))
            for pi, pairs in enumerate(([(p0, p1), (p2, p3)], [(p0, p3), (p1, p2)])):
                partner = {x: y for x, y in pairs} | {y: x for x, y in pairs}
                nopp = [0] * (4 * m)
                valid = True
                for v in range(n):
                    if v in (u, w):
                        continue
                    nv = remap[v]
                    for s in range(4):
                        od = d.opp[4 * v + s]
                        steps = 0
                        while (od >> 2) in (u, w):
                            od = d.opp[partner[od]]
                            steps += 1
                            if steps > 8:
                                valid = False
                                break
                        if not valid:
                            break
                        nopp[4 * nv + s] = (remap[od >> 2] << 2) | (od & 3)
                    if not valid:
                        break
                if not valid:
                    continue
                cand = MarkedDiagram(num_vertices=m, opp=tuple(nopp), vtype=nvt, free_loops=0)
                if not M._is_planar_4reg_diagram(cand):
                    continue
                if not (certify_classical_unlink(cand.resolve(1)).is_unlink and certify_classical_unlink(cand.resolve(-1)).is_unlink):
                    continue
                cst = cand.surface_topology()
                if cst.num_components != 1 or cst.euler_characteristic != st.euler_characteristic + 2:
                    continue
                wp_c = wirtinger_presentation(cand)
                if wp_c.num_generators > 1:
                    fqp_c = finite_quotient_profile(wp_c, ("S3", "A4", "F21", "A5"))
                    if fqp_c["S3"] > 6 or fqp_c["A4"] > 12 or fqp_c["F21"] > 21 or fqp_c["A5"] > 60:
                        continue
                info = {"u": u, "w": w, "pairing": "loop(p0p1)(p2p3)" if pi == 0 else "thru(p0p3)(p1p2)", "tube_cfg": tube_cfg,
                        "cand_epd": cand.to_epd(), "cand_surf": surf_inv(cand)}
                red, sp, cp = M.simplify_and_explore_yoshikawa(cand, max_states=24)
                if sp:
                    continue
                if cp or red.num_vertices <= 7:
                    info["via"] = "composite" if cp else f"reduced_to_{red.num_vertices}"
                    return True, info
                if red.num_vertices in (8, 9):
                    t_kept = fe.tunnel_diagrams_c(red.num_vertices, [red], depth_plus1=96, depth_plus2=64, max_bfs_states=512)
                    if len(t_kept) == 0:
                        info["via"] = "tunnel_empty"
                        return True, info
                    v2, i2 = instrumented(t_kept[0], _depth + 1)
                    if v2:
                        info["via"] = "recursive"
                        info["sub"] = i2
                        return True, info
    return False, None


CYC = {"Z": {"S3": 6, "Q8": 8, "A4": 12, "S4": 24}, "Z2": {"S3": 4, "Q8": 2, "A4": 4, "S4": 10}}
out = []
npruned = 0
for d in tunneled:
    orig = M.has_certified_cyclic_unknotted_tube(d)
    v, info = instrumented(d)
    assert v == orig, "instrumented copy disagrees with original"
    if not v:
        continue
    npruned += 1
    prof = {g: indep_hom_count(d, g) for g in ("S3", "Q8", "A4", "S4")}
    cyc = "Z" if prof == CYC["Z"] else ("Z2" if prof == CYC["Z2"] else "NONCYCLIC")
    rec = {"epd": d.to_epd(), "surf": surf_inv(d), "prof": prof, "pi1_class": cyc, **info}
    out.append(rec)
    print(json.dumps(rec), flush=True)
print(f"ch={ch}: tunneled {len(tunneled)} pruned {npruned};",
      "NONCYCLIC pruned:", sum(1 for r in out if r["pi1_class"] == "NONCYCLIC"),
      "tube_cfg+thru pruned:", sum(1 for r in out if r["tube_cfg"] and r["pairing"].startswith("thru")))
json.dump(out, open(f"{HERE}/t5_pruned_ch{ch}.json", "w"))
