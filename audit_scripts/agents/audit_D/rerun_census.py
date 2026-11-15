import sys, json, time
from pathlib import Path
from surfaceknots.triplane_census import run_triplane_b4_census
dd = Path(sys.argv[1])
t=time.time()
res = run_triplane_b4_census(dd, max_c_tri=6)
print("elapsed", time.time()-t)
print(json.dumps(res["enumeration_stats"], indent=1))
print(json.dumps(res["gong_lewis_monkman_osnes_2026_verification"], indent=1))
