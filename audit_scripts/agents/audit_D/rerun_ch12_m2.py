import time, pickle, os, json
from collections import defaultdict
from surfaceknots.fast_enum import enumerate_targeted_prime_ch_census
from surfaceknots.yoshikawa_data import compute_full_invariant_record, invariant_discrimination_key
t = time.time()
diags, stats = enumerate_targeted_prime_ch_census(ch=12, m_min=2, m_max=2, target_components=1, target_orientability=1,
                                                  target_chi_min=2, target_chi_max=2, num_workers=8)
print("elapsed", time.time() - t, stats, flush=True)
recs = []
for d in diags:
    r = compute_full_invariant_record(d, label="x", include_sublinks=False)
    recs.append((d.to_epd(), r))
pickle.dump(recs, open(os.path.dirname(os.path.abspath(__file__)) + "/ch12_m2_survivors.pkl", "wb"))
weak = defaultdict(list); strong = defaultdict(list)
for epd, r in recs:
    weak[str(invariant_discrimination_key(r))].append(epd)
    strong[(str(r["e1_generators"]), str(sorted(r["hom_counts"].items())), str(r["sigma2_h1_torsion"]))].append(epd)
print("survivor diagrams:", len(recs), "| classes by repo 10-layer key:", len(weak), "| classes by (E1 generators, hom counts, Sigma2):", len(strong))
for k, v in sorted(strong.items()): print(len(v), k[0], k[2])
