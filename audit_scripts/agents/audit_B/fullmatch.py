import sys, json, time
from concurrent.futures import ProcessPoolExecutor
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.yoshikawa_data import YOSHIKAWA_TABLE, compute_full_invariant_record, invariant_discrimination_key
def rec(a):
    epd, tag = a
    try:
        r = compute_full_invariant_record(MarkedDiagram.from_epd(epd), label=tag, include_sublinks=False)
        return tag, invariant_discrimination_key(r), None
    except Exception as e:
        return tag, None, repr(e)[:200]
if __name__ == "__main__":
    N = int(sys.argv[1]); path = sys.argv[2]
    surv = json.load(open(path))
    tasks = [(e.canonical_epd, "Y:" + e.label) for e in YOSHIKAWA_TABLE if e.ch_index > 0] + [(s, f"S{i}") for i, s in enumerate(surv)]
    t = time.time()
    with ProcessPoolExecutor(6) as p: out = list(p.map(rec, tasks, chunksize=1))
    ykeys = {}
    for tag, k, err in out:
        if tag.startswith("Y:"): ykeys.setdefault(k, []).append(tag[2:])
    cnt = {}
    for tag, k, err in out:
        if tag.startswith("Y:"): continue
        lab = "ERR " + err if err else ",".join(ykeys.get(k, ["UNMATCHED"]))
        cnt[lab] = cnt.get(lab, 0) + 1
    print(f"N={N} full-invariant match of {len(surv)} survivors (time {time.time()-t:.0f}s):")
    for k in sorted(cnt): print("   ", cnt[k], k)
    exp = set(e.label for e in YOSHIKAWA_TABLE if e.ch_index == N)
    print("   missing expected ch=N classes:", sorted(exp - set(cnt)))
