import certcount, json, time, sys, aud
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.fast_enum import enumerate_prime_ch_census
from surfaceknots.moves import has_certified_cyclic_unknotted_tube
from surfaceknots.yoshikawa_data import YOSHIKAWA_TABLE
N = int(sys.argv[1])
t = time.time()
surv, st = enumerate_prime_ch_census(N, num_workers=2, filter_cyclic_tubes=False)
print("N", N, st, "certcounts", dict(certcount.COUNTS), "time %.0f" % (time.time() - t))
json.dump([d.to_epd() for d in surv], open(f"tunneled_nofilter_{N}.json", "w"))
ytab = {}
for e in YOSHIKAWA_TABLE:
    d = MarkedDiagram.from_epd(e.canonical_epd)
    if d.num_vertices: ytab.setdefault(aud.inv(d), []).append(e.label)
for d in surv:
    if has_certified_cyclic_unknotted_tube(d):
        i = aud.inv(d)
        print("TUBE-PRUNED topo", i[0], "S3", dict(i[1])["S3"], "A5", dict(i[1])["A5"], "matches", ytab.get(i), d.to_epd())
