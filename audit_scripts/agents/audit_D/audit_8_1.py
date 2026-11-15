import json
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.moves import canonical_marked_key, list_level_yoshikawa_neighbors
from surfaceknots.yoshikawa_data import YOSHIKAWA_BY_LABEL
d = json.load(open("data/complexity/triplane_b4_census.json"))
for e in d["entries"]:
    for m in e["catalog_matches"]:
        if m.startswith("Hub:") and e["best_compiled_epd"] and m[4:] in YOSHIKAWA_BY_LABEL and e["best_compiled_ch"] >= 6:
            md = MarkedDiagram.from_epd(e["best_compiled_epd"])
            ref = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL[m[4:]].canonical_epd)
            k1 = canonical_marked_key(md, up_to_mirror_and_reverse=True); k2 = canonical_marked_key(ref, up_to_mirror_and_reverse=True)
            # BFS over level moves from compiled diagram
            seen = {k1}; q = [md]; hit = (k1 == k2); i = 0
            while i < len(q) and not hit and len(seen) < 20000:
                cur = q[i]; i += 1
                for nb in list_level_yoshikawa_neighbors(cur):
                    nk = canonical_marked_key(nb, up_to_mirror_and_reverse=True)
                    if nk == k2: hit = True; break
                    if nk not in seen: seen.add(nk); q.append(nb)
            print(m, "compiled ch", md.ch_index, "(c,m)=", md.num_classical, md.num_marked, "| identical diagram:", k1 == k2, "| connected by level Yoshikawa moves:", hit, "| states explored:", len(seen))
