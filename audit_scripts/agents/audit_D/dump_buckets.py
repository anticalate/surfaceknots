"""Run the census enumeration (tc_variant, ARCH_MODE env) and dump every prekey bucket with its representative,
full-strength group checks, and which o4 class it merges into."""
import sys, os, pickle, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tc_variant
from surfaceknots.triplane import compute_triplane_invariant_record, triplane_o4_key, triplane_fundamental_group
t = time.time()
buckets, stats = tc_variant.enumerate_triplane_up_to_b4(max_c_tri=6)
print("mode", os.environ.get("ARCH_MODE"), "elapsed", time.time() - t, {k: v for k, v in stats.items() if not k.startswith("tangle")}, flush=True)
out = []
for pk, bkt in buckets.items():
    rec = compute_triplane_invariant_record(bkt.rep_tp)
    ok = triplane_o4_key(rec)
    out.append(dict(prekey=pk, rep=tc_variant._serialize_triplane(bkt.rep_tp), n=bkt.valid_triple_count, min_b=bkt.min_b, min_c=bkt.min_c_tri,
                    triples=sorted(bkt.crossing_triples_seen), o4=ok, comp=rec["component_types"], e=rec["normal_euler_numbers"],
                    pi1=rec["pi1_presentation"], e1=rec["e1_generators"], hom=rec["hom_counts"],
                    min_unconc=tc_variant._serialize_triplane(bkt.min_unconc_tp) if bkt.min_unconc_tp else None,
                    min_conc=tc_variant._serialize_triplane(bkt.min_conc_tp) if bkt.min_conc_tp else None,
                    best_epd=bkt.best_compiled_md.to_epd() if bkt.best_compiled_md else None,
                    best_src=tc_variant._serialize_triplane(bkt.best_compiled_source_tp) if bkt.best_compiled_source_tp else None,
                    cc=bkt.conc_compilable_count, ci=bkt.conc_interleaved_count))
pickle.dump(out, open(os.path.dirname(os.path.abspath(__file__)) + f"/buckets_{os.environ.get('ARCH_MODE','orig')}.pkl", "wb"))
for o in sorted(out, key=lambda o: (str(o["comp"]), str(o["e"]), o["pi1"])):
    print(o["comp"], o["e"], "n=%d" % o["n"], "minb=%d minc=%d" % (o["min_b"], o["min_c"]), "E1=", o["e1"], "pi1=", o["pi1"][:100], "| prekey:", str(o["prekey"])[:90])
