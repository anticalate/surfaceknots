# For every named diagram: |Hom(F(a,b)/<<x = w y^s w^-1>>, G)| must equal |Hom(pi_1(diagram), G)|, G in S3,Q8,A4,S4.
import sys, json, importlib.util
spec=importlib.util.spec_from_file_location('hf', sys.argv[3]+'/handle_free.py'); hf=importlib.util.module_from_spec(spec); spec.loader.exec_module(hf)
sys.path.insert(0, sys.argv[3])
from indep_pi1_lib import GROUPS, count_homs
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.moves import list_level_yoshikawa_neighbors, canonical_marked_key
from surfaceknots.yoshikawa_data import YOSHIKAWA_TABLE
def pres_homs(X,Y,W,S,G):
    n,mul,inv=G; tot=0
    for a in range(n):
        for b in range(n):
            img={'a':a,'A':inv[a],'b':b,'B':inv[b]}
            w=None
            for ch in W: w=img[ch] if w is None else mul[w][img[ch]]
            x=img[X]; y=img[Y] if S==1 else inv[img[Y]]
            rhs=y if w is None else mul[mul[w][y]][inv[w]]
            tot+= (x==rhs)
    return tot
code=sys.argv[2]
pool=[(e.label,e.canonical_epd) for e in YOSHIKAWA_TABLE]
c11=json.load(open(code+'/data/expanded_ch11.json'))
for r in c11['entries']:
    for i,epd in enumerate(r['move_component_epds']): pool.append((r['label']+(f".{'abc'[i]}" if len(r['move_component_epds'])>1 else ''),epd))
rows=[]
for lab,epd in pool:
    d0=MarkedDiagram.from_epd(epd)
    if d0.num_vertices<4: continue
    seen={canonical_marked_key(d0):d0}; q=[d0]; i=0; found=None
    while i<len(q) and i<400 and found is None:
        cur=q[i]; i+=1
        for u,w,cl,nm in hf.analyse(cur):
            if cl is not None: found=(cl,nm,cur); break
        if found: break
        for nb in list_level_yoshikawa_neighbors(cur):
            k=canonical_marked_key(nb)
            if k not in seen: seen[k]=nb; q.append(nb)
    if not found: continue
    (X,Y,W,S),nm,cur=found
    ok=all(pres_homs(X,Y,W,S,GROUPS[g])==count_homs(d0,GROUPS[g]) for g in ('S3','Q8','A4','S4'))
    rows.append((lab,nm,ok,d0.num_vertices))
    print(f"{lab:14s} ch={d0.num_vertices:2d}  {nm:34s} group check: {'OK' if ok else 'FAIL'}",flush=True)
json.dump(rows,open(sys.argv[3]+'/handle_names.json','w'),indent=1)
