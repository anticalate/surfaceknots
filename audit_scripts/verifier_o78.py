"""Do genuine Omega_7 / Omega_8 (and Omega_3/4, Omega_5) steps pass Carbon's step classifier?
Random walks from the table and ch=11 diagrams (stabilise, poke, slide), collecting Omega_7 and Omega_8
neighbours from the repo's generators along the way.  Run from code/ with PYTHONPATH=src."""
import json, random, sys, collections, time
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots import move_search as ms
from surfaceknots import moves as mv
from surfaceknots.yoshikawa_data import YOSHIKAWA_TABLE
rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
c11 = json.load(open("data/expanded_ch11.json"))
pool = [MarkedDiagram.from_epd(e.canonical_epd) for e in YOSHIKAWA_TABLE if e.canonical_epd]
pool += [MarkedDiagram.from_epd(r["epd"]) for r in c11["entries"]]
pool = [d for d in pool if d.num_vertices >= 6]
gens = {"O34": mv._list_omega3_and_4_neighbors, "O5": mv._list_omega5_neighbors,
        "O7": mv._list_omega7_neighbors, "O8": mv._list_omega8_neighbors}
walkers = [mv.list_omega1_6_stabilizations, mv.list_omega2_pokes, mv._list_omega3_and_4_neighbors,
           mv._list_omega5_neighbors, mv._list_omega7_neighbors]
res = collections.Counter(); t0 = time.time(); cap = {"O34": 60, "O5": 60, "O7": 150, "O8": 60}
seen = collections.Counter()
while time.time() - t0 < float(sys.argv[2]) and (seen["O7"] < cap["O7"] or seen["O8"] < cap["O8"]):
    d = rng.choice(pool); n0 = d.num_vertices
    for step in range(12):
        for tag, g in gens.items():
            if seen[tag] >= cap[tag]: continue
            L = g(d)
            if not L: continue
            nb = rng.choice(L); seen[tag] += 1
            k = ms._indep_is_single_local_move(d, nb, True)
            res[(tag, k)] += 1
        w = rng.choice(walkers)
        if d.num_vertices >= n0 + 3 and w in walkers[:2]: w = walkers[2]
        L = w(d)
        if L: d = rng.choice(L)
print("genuine steps by generator and classifier verdict:")
for k, v in sorted(res.items(), key=str): print("  ", k, v)
