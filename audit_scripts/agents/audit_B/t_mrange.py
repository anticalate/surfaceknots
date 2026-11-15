"""Marker-count ranges excluded by the default pipeline (m<2 and m>ch-4): do they contain non-reducible classes?"""
import sys, time, json
import aud
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.fast_enum import enumerate_census_level_candidates, cascade_tunnel_diagrams_c
from surfaceknots.moves import has_certified_cyclic_unknotted_tube
N = int(sys.argv[1])
for (a, b) in [(0, 1), (max(2, N - 4) + 1, N)]:
    t = time.time()
    lvl, st = enumerate_census_level_candidates(N, m_min=a, m_max=b, num_workers=2)
    tun = cascade_tunnel_diagrams_c(N, lvl) if lvl else []
    filt = [d for d in tun if not (N >= 9 and has_certified_cyclic_unknotted_tube(d))]
    groups = {}
    for d in filt: groups.setdefault(aud.inv(d), []).append(d)
    print(f"N={N} m in [{a},{b}]: raw={st['raw_cands']} level_kept={len(lvl)} tunneled={len(tun)} after_tube={len(filt)} time={time.time()-t:.1f}s")
    for k, ds in groups.items():
        print("    ", len(ds), "diagrams topo", k[0], "homS3", dict(k[1]).get('S3'), "m=", sorted(set(sum(1 for t in d.vtype if t in 'YZ') for d in ds)), ds[0].to_epd())
    json.dump([d.to_epd() for d in filt], open(f"mrange_{N}_{a}_{b}.json", "w"))
