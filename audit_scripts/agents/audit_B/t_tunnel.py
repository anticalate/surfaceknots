"""Tunnel soundness: stabilized states valid; merged classes have constant invariants; Yoshikawa reps never killed."""
import sys, time, ctypes, random, json
import aud
from aud import L, arrs, mk
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.fast_enum import enumerate_census_level_candidates
from surfaceknots.yoshikawa_data import YOSHIKAWA_TABLE

HOOK = ctypes.CFUNCTYPE(None, ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int))
captured = []
def _hook(n, m, opp, vt):
    captured.append(mk(m, [opp[i] for i in range(4*m)], [vt[i] for i in range(m)]))
hook = HOOK(_hook)

def stab_states(d, d1=1, d2=1):
    captured.clear()
    L.x_set_stab_hook(hook)
    n = d.num_vertices; o, v = arrs(d)
    nb = ctypes.c_int(0)
    bo = (ctypes.c_int * (1024*4*n))(); bv = (ctypes.c_int * (1024*n))()
    r = L.x_run_probe(n, o, v, d1, d2, ctypes.byref(nb), bo, bv, 1024)
    L.x_set_stab_hook(ctypes.cast(None, HOOK))
    bridges = [mk(n, bo[i*4*n:(i+1)*4*n], bv[i*n:(i+1)*n]) for i in range(nb.value)]
    return r, list(captured), bridges

def tunnel_classes(n, diags, d1, d2, mbfs):
    num = len(diags)
    io = (ctypes.c_int * (num*4*n))(*[x for d in diags for x in d.opp])
    iv = (ctypes.c_int * (num*n))(*[aud.VTI[t] for d in diags for t in d.vtype])
    oo = (ctypes.c_int * ((num+16)*4*n))(); ov = (ctypes.c_int * ((num+16)*n))()
    root = (ctypes.c_int * num)(); dead = (ctypes.c_int * num)()
    L.x_set_root_out(root, dead)
    k = L.tunnel_and_merge_diagrams_c(n, num, io, iv, d1, d2, mbfs, oo, ov, num+16)
    L.x_set_root_out(None, None)
    outs = [mk(n, oo[i*4*n:(i+1)*4*n], ov[i*n:(i+1)*n]) for i in range(k)]
    return list(root), list(dead), outs

if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "stab":
        n = int(sys.argv[2]); cnt = int(sys.argv[3])
        level, st = enumerate_census_level_candidates(n, num_workers=2)
        random.seed(7); samp = level if len(level) <= cnt else random.sample(level, cnt)
        tot = dict(diag=0, states=0, nonplanar=0, nonadm=0, inv=0, killed=0, bridges=0, bridge_inv=0, greedy_inv=0, greedy_checked=0)
        for d in samp:
            i0 = aud.inv(d); tot["diag"] += 1
            r, states, bridges = stab_states(d, 1, 1)
            tot["killed"] += r
            seen = set()
            for s in states:
                k = aud.canon_key(s)
                if k in seen: continue
                seen.add(k); tot["states"] += 1
                if not aud.is_planar(s): tot["nonplanar"] += 1; print("NONPLANAR", d.to_epd(), s.to_epd()); continue
                if not aud.admissible(s): tot["nonadm"] += 1; print("NONADM", s.to_epd()); continue
                if aud.inv(s) != i0: tot["inv"] += 1; print("INV", d.to_epd(), s.to_epd())
                # greedy on level neighbours of stabilized states
                for x in aud.level_neighbors(s)[:6]:
                    g = aud.greedy_reduce(x); tot["greedy_checked"] += 1
                    if g.num_vertices == 0:
                        continue
                    if aud.inv(g) != i0: tot["greedy_inv"] += 1; print("GREEDY_INV", x.to_epd(), g.to_epd(), g.free_loops)
            for b in bridges:
                tot["bridges"] += 1
                if aud.inv(b) != i0: tot["bridge_inv"] += 1; print("BRIDGE_INV", d.to_epd(), b.to_epd())
        print(f"stab n={n}", tot)
    elif mode == "classes":
        n = int(sys.argv[2]); d1 = int(sys.argv[3]); d2 = int(sys.argv[4]); mbfs = int(sys.argv[5])
        level, st = enumerate_census_level_candidates(n, num_workers=2)
        t = time.time()
        root, dead, outs = tunnel_classes(n, level, d1, d2, mbfs)
        invs = [aud.inv(d) for d in level]
        cls = {}
        for i, r in enumerate(root): cls.setdefault(r, []).append(i)
        bad = 0; multi = 0
        for r, mem in cls.items():
            if len(mem) > 1: multi += 1
            if len(set(invs[i] for i in mem)) > 1:
                bad += 1; print("CLASS INV MISMATCH root", r, "dead", dead[mem[0]], [level[i].to_epd() for i in mem[:3]])
        # dead classes whose invariants equal a ch==n Yoshikawa entry
        yk = {}
        for e in YOSHIKAWA_TABLE:
            if e.ch_index == n: yk[aud.inv(MarkedDiagram.from_epd(e.canonical_epd))] = e.label
        dead_y = {}; live_y = {}
        for r, mem in cls.items():
            lab = yk.get(invs[mem[0]])
            if lab:
                (dead_y if dead[mem[0]] else live_y).setdefault(lab, 0)
                if dead[mem[0]]: dead_y[lab] += 1
                else: live_y[lab] += 1
        print(f"classes n={n} d=({d1},{d2},{mbfs}) inputs={len(level)} classes={len(cls)} multi-member={multi} live_out={len(outs)} class_inv_mismatch={bad} dead_classes_with_chN_yoshikawa_inv={dead_y} live={live_y} time={time.time()-t:.1f}s")
    elif mode == "yosh":
        for e in YOSHIKAWA_TABLE:
            d = MarkedDiagram.from_epd(e.canonical_epd); n = d.num_vertices
            if n < 3: continue
            root, dead, outs = tunnel_classes(n, [d], 256, 192, 2048)
            print(e.label, "survives" if len(outs) == 1 else "KILLED!!", flush=True)
