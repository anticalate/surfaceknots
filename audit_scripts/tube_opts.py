import json, sys
from surfaceknots import MarkedDiagram, YOSHIKAWA_TABLE
from surfaceknots.moves import _is_planar_4reg_diagram, has_certified_cyclic_unknotted_tube
from surfaceknots.unlink import certify_classical_unlink
from surfaceknots.invariants.fundamental_group import wirtinger_presentation, finite_quotient_profile
GR=("S3","Q8","A4","S4","A5")
def homs(d):
    return dict(finite_quotient_profile(wirtinger_presentation(d),GR))
def options(d):
    n=d.num_vertices; st=d.surface_topology(); res=[]
    for u in range(n):
        if d.vtype[u] not in ("Y","Z"): continue
        for su in range(4):
            d_a=d.opp[4*u+su]; d_b=d.opp[4*u+((su+1)&3)]; w=d_a>>2
            if w<=u or (d_b>>2)!=w or d.vtype[w] not in ("Y","Z"): continue
            sw_a=d_a&3; sw_b=d_b&3
            if ((sw_b+1)&3)!=sw_a: continue
            p0,p1=4*u+((su+2)&3),4*u+((su+3)&3); p2,p3=4*w+((sw_b+2)&3),4*w+((sw_b+3)&3)
            remap=[v-(1 if v>u else 0)-(1 if v>w else 0) for v in range(n)]; m=n-2
            nvt=tuple(d.vtype[v] for v in range(n) if v not in (u,w))
            for oi,pairs in enumerate(([(p0,p1),(p2,p3)],[(p0,p3),(p1,p2)])):
                partner={x:y for x,y in pairs}|{y:x for x,y in pairs}
                nopp=[0]*(4*m); valid=True
                for v in range(n):
                    if v in (u,w): continue
                    nv=remap[v]
                    for s in range(4):
                        od=d.opp[4*v+s]; steps=0
                        while (od>>2) in (u,w):
                            od=d.opp[partner[od]]; steps+=1
                            if steps>8: valid=False; break
                        if not valid: break
                        nopp[4*nv+s]=(remap[od>>2]<<2)|(od&3)
                    if not valid: break
                if not valid: continue
                cand=MarkedDiagram(num_vertices=m,opp=tuple(nopp),vtype=nvt,free_loops=0)
                if not _is_planar_4reg_diagram(cand): continue
                if not (certify_classical_unlink(cand.resolve(1)).is_unlink and certify_classical_unlink(cand.resolve(-1)).is_unlink): continue
                cst=cand.surface_topology()
                if cst.num_components!=1 or cst.euler_characteristic!=st.euler_characteristic+2: continue
                res.append((u,w,d.vtype[u],d.vtype[w],oi,cand))
    return res
pool=[]
for e in YOSHIKAWA_TABLE: pool.append((e.label,MarkedDiagram.from_epd(e.canonical_epd)))
ch11=json.load(open('data/expanded_ch11.json'))
for r in ch11['entries']:
    for i,epd in enumerate(r['move_component_epds']): pool.append((r['label']+f'.{i}',MarkedDiagram.from_epd(epd)))
for i,r in enumerate(ch11['residual_collisions']): pool.append((f'resid{i}:'+r['collision_matches'][0],MarkedDiagram.from_epd(r['epd'])))
fam=json.load(open('data/families/family_catalog.json'))
ents=fam['entries'] if isinstance(fam,dict) and 'entries' in fam else fam
def walk(o):
    if isinstance(o,dict):
        if 'epd' in o and isinstance(o['epd'],str): yield (o.get('family_id') or o.get('label') or '?', o['epd'])
        for k,v in o.items():
            if k.endswith('epd') and k!='epd' and isinstance(v,str) and v.startswith(('X','Y','Z')): yield ((o.get('family_id') or '?')+':'+k, v)
            yield from walk(v)
    elif isinstance(o,list):
        for v in o: yield from walk(v)
seen=set()
for lab,epd in walk(fam):
    if epd in seen: continue
    seen.add(epd)
    try: pool.append(('fam:'+lab,MarkedDiagram.from_epd(epd)))
    except Exception as ex: pass
print('pool',len(pool),flush=True)
bad=0; nopt=0
for lab,d in pool:
    if d.num_vertices<5 or d.surface_topology().num_components!=1: continue
    opts=options(d)
    if not opts: continue
    hd=homs(d)
    for (u,w,tu,tw,oi,cand) in opts:
        nopt+=1
        hc=homs(cand)
        viol=[g for g in GR if hc[g]<hd[g]]
        tag='VIOLATION' if viol else 'ok'
        if viol: bad+=1
        print(f"{tag:9s} {lab:28s} n={d.num_vertices} bigon({u}{tu},{w}{tw}) opt{oi} base_n={cand.num_vertices} hom_d={hd} hom_base={hc} tubeflag={has_certified_cyclic_unknotted_tube(d) if viol else ''}",flush=True)
print('options checked',nopt,'violations',bad)
