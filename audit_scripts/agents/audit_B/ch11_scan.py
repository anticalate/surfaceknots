"""Replicates Step 1-2 of scripts/reproduce_datasets.py::run_ch11_enumeration with identical partitioning
(12 batches x 6 micro-parts, args (11, part, 2, 7, 384, False), per-batch merge_diagrams_c(.., 512)),
but runs the micro-parts on more threads. Writes the same checkpoint files into a scratch data dir."""
import sys, json, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import certcount
from surfaceknots.fast_enum import _get_lib, _run_c_census_shard, merge_diagrams_c
from surfaceknots.quadrangulations import generate_prime_shadows

data_dir = Path(sys.argv[1]); threads = int(sys.argv[2])
data_dir.mkdir(parents=True, exist_ok=True)
_get_lib()
t0 = time.time()
opp_list = [list(s.opp) for s in generate_prime_shadows(11)]
(data_dir / "_shadows_ch11.json").write_text(json.dumps(opp_list))
num_batches = 12
batch_size = (len(opp_list) + num_batches - 1) // num_batches
tasks = []
for b_idx in range(num_batches):
    sub = [tuple(x) for x in opp_list[b_idx * batch_size:(b_idx + 1) * batch_size]]
    micro_sz = max(1, (len(sub) + 5) // 6)
    parts = [sub[i:i + micro_sz] for i in range(0, len(sub), micro_sz)]
    for p_idx, p in enumerate(parts):
        if p: tasks.append((b_idx, p_idx, p))
def run(t):
    b, pi, p = t
    ts = time.time()
    d_list, st = _run_c_census_shard((11, p, 2, 7, 384, False))
    print(f"  part b{b}.{pi}: shadows={len(p)} raw={st[0]} comps={st[1]} disc={st[2]} kept={len(d_list)} {time.time()-ts:.0f}s", flush=True)
    return b, pi, d_list, st
with ThreadPoolExecutor(max_workers=threads) as ex:
    results = list(ex.map(run, tasks))
print(f"scan done in {time.time()-t0:.0f}s", flush=True)
shard_state = {"completed": {}, "stats": {}}
for b_idx in range(num_batches):
    rs = [r for r in results if r[0] == b_idx]
    batch_diags = []; raw_c = comp_c = disc_c = 0
    for _, _, d_list, st in rs:
        batch_diags.extend(d_list); raw_c += st[0]; comp_c += st[1]; disc_c += st[2]
    before = dict(certcount.COUNTS)
    merged_b = merge_diagrams_c(11, batch_diags, max_bfs_states=512)
    shard_state["completed"][str(b_idx)] = [d.to_epd() for d in merged_b]
    shard_state["stats"][str(b_idx)] = [raw_c, comp_c, disc_c, len(batch_diags), len(merged_b)]
    print(f"[Step 2] Batch {b_idx+1}/12: raw={raw_c}, comps={comp_c}, kept_in={len(batch_diags)}, merged={len(merged_b)}", flush=True)
(data_dir / "_ch11_shards.json").write_text(json.dumps(shard_state))
tot = [sum(v[i] for v in shard_state["stats"].values()) for i in range(5)]
print(f"TOTAL raw={tot[0]} comps={tot[1]} disc={tot[2]} kept_in={tot[3]} merged_sum={tot[4]} time={time.time()-t0:.0f}s certcounts={dict(certcount.COUNTS)}", flush=True)
