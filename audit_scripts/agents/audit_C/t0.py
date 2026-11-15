import time, sys
from surfaceknots import YOSHIKAWA_BY_LABEL, MarkedDiagram
from surfaceknots.yoshikawa_data import compute_full_invariant_record, invariant_filtration_layers
for lbl, e in YOSHIKAWA_BY_LABEL.items():
    d = MarkedDiagram.from_epd(e.canonical_epd)
    t=time.time()
    rec = compute_full_invariant_record(d, lbl, include_sublinks=False, compute_so4_symmetry_bfs=False)
    print(lbl, round(time.time()-t,2), flush=True)
    for k,v in invariant_filtration_layers(rec):
        print("   ",k, str(v)[:300])
