import sys, json, random, time, traceback, os
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variants import *
from surfaceknots import YOSHIKAWA_BY_LABEL, MarkedDiagram
from surfaceknots.yoshikawa_data import compute_full_invariant_record, invariant_filtration_layers
from surfaceknots.unlink import certify_classical_unlink

def layers(d):
    rec = compute_full_invariant_record(d, "x", include_sublinks=False, compute_so4_symmetry_bfs=False)
    return dict(invariant_filtration_layers(rec)), rec

def work(args):
    lbl, epd, seed, nsteps, use_level, allowed = args
    rng = random.Random(seed)
    d0 = MarkedDiagram.from_epd(epd)
    try:
        d, path = random_variant(d0, rng, nsteps=nsteps, use_level=use_level, allowed=allowed)
    except Exception as ex:
        return (lbl, seed, "GENFAIL", repr(ex), None, None)
    if not is_planar(d):
        return (lbl, seed, "NONPLANAR", path, d.to_epd(), None)
    try:
        L, rec = layers(d)
    except Exception as ex:
        return (lbl, seed, "EXC", path, d.to_epd(), traceback.format_exc()[-600:])
    return (lbl, seed, "OK", path, d.to_epd(), {k: repr(v) for k, v in L.items()})

if __name__ == "__main__":
    N = int(sys.argv[1]); nsteps = int(sys.argv[2]); use_level = sys.argv[3] == "level"
    out = sys.argv[4]
    allowed = sys.argv[5].split(",") if len(sys.argv) > 5 and sys.argv[5] != "-" else None
    only = sys.argv[6].split(";") if len(sys.argv) > 6 else None
    ds = {lbl: e.canonical_epd for lbl, e in YOSHIKAWA_BY_LABEL.items() if lbl != "0_1"}
    j = json.load(open('data/expanded_ch11.json'))
    for e in j['entries']: ds[e['label']] = e['epd']
    if only: ds = {k: v for k, v in ds.items() if k in only}
    base = {}
    for lbl, epd in ds.items():
        L, _ = layers(MarkedDiagram.from_epd(epd))
        base[lbl] = {k: repr(v) for k, v in L.items()}
    jobs = [(lbl, epd, 1000 * i + hash(lbl) % 997, nsteps, use_level, allowed) for lbl, epd in ds.items() for i in range(N)]
    res = []
    with Pool(7) as p:
        for r in p.imap_unordered(work, jobs, chunksize=2):
            res.append(r)
    stats = {}; fails = []
    nok = 0
    for lbl, seed, st, path, epd, L in res:
        if st != "OK":
            fails.append((lbl, seed, st, path, epd, L)); continue
        nok += 1
        for k, v in L.items():
            if v != base[lbl][k]:
                stats.setdefault(k, []).append((lbl, seed, path, epd, base[lbl][k][:300], v[:300]))
    json.dump({"stats": stats, "fails": fails, "nok": nok, "njobs": len(jobs)}, open(out, "w"), indent=1)
    print("jobs", len(jobs), "ok", nok, "fails", len(fails))
    from collections import Counter
    print("fail kinds", Counter((f[2]) for f in fails))
    for k, v in stats.items():
        print("LAYER", k, "mismatches", len(v), "surfaces", sorted({x[0] for x in v}))
