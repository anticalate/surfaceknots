"""Blind end-to-end run of the default pipeline at ch=N; compare survivors with Yoshikawa table by cheap invariants."""
import sys, time, json
import aud
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.fast_enum import enumerate_prime_ch_census
from surfaceknots.yoshikawa_data import YOSHIKAWA_TABLE

N = int(sys.argv[1]); workers = int(sys.argv[2])
kw = {}
if len(sys.argv) > 3: kw = json.loads(sys.argv[3])
t = time.time()
surv, stats = enumerate_prime_ch_census(N, num_workers=workers, **kw)
el = time.time() - t
ytab = {}
for e in YOSHIKAWA_TABLE:
    d = MarkedDiagram.from_epd(e.canonical_epd)
    if d.num_vertices == 0: continue
    ytab.setdefault(aud.inv(d), []).append((e.label, d.num_vertices))
groups = {}
for d in surv:
    groups.setdefault(aud.inv(d), []).append(d)
print(f"N={N} time={el:.1f}s stats={stats} survivors={len(surv)} inv-classes={len(groups)}")
exp = sorted(e.label for e in YOSHIKAWA_TABLE if e.ch_index == N)
found = set()
for k, ds in sorted(groups.items(), key=lambda kv: str(kv[0])):
    lab = ytab.get(k)
    tag = "UNMATCHED" if lab is None else ",".join(f"{l}(ch{c})" for l, c in lab)
    if lab:
        for l, c in lab:
            if c == N: found.add(l)
    print(f"   {len(ds):3d} diagrams  topo={k[0]}  -> {tag}")
print("   expected at ch=N:", exp)
print("   missing expected:", sorted(set(exp) - found))
json.dump([d.to_epd() for d in surv], open(f"blind_surv_{N}.json", "w"))
