import sys, json, time, os
from pathlib import Path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tc_variant
dd = Path(sys.argv[1])
t=time.time()
res = tc_variant.run_triplane_b4_census(dd, max_c_tri=int(sys.argv[2]) if len(sys.argv)>2 else 6)
print("ARCH_MODE", os.environ.get("ARCH_MODE"), "elapsed", time.time()-t)
print(json.dumps(res["enumeration_stats"], indent=1))
print(json.dumps(res["gong_lewis_monkman_osnes_2026_verification"], indent=1))
for e in res["entries"]:
    inv = e["invariants"]
    print(f"  c_tri={e['min_c_tri']} (conc={e['min_concentrated_c_tri']}, unconc={e['min_unconcentrated_c_tri']}) b={e['min_bridge_in_census']} comp_ch={e['best_compiled_ch']} n={e['valid_triple_count']} cc={e['conc_compilable_count']} ci={e['conc_interleaved_count']} topo={inv['component_types']} e={inv['normal_euler_numbers']} E1={inv['e1_generators']} pi1={inv['pi1_presentation'][:70]} matches={e['catalog_matches']} triples={e['crossing_triples_seen']}")
