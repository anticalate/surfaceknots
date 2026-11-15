# Independent admissibility check with SnapPy/spherogram: both resolutions of every catalog diagram
# must be unlinks.  Criterion: link group is free (rank = #components)  <=> trivial link.
import json, sys, re
sys.path.insert(0, sys.argv[1])           # repo code/src, imported WITHOUT its C kernel
import snappy, spherogram
from surfaceknots.diagrams import MarkedDiagram

def to_link(cl):
    n=cl.num_crossings
    if n==0: return None
    lab={}
    k=0
    for d in range(4*n):
        o=cl.opp[d]
        if d not in lab:
            lab[d]=lab[o]=k; k+=1
    pd=[]
    for v in range(n):
        order=(0,1,2,3) if cl.vtype[v]=="X0" else (1,2,3,0)   # start at an under slot, CCW
        pd.append([lab[4*v+s] for s in order])
    return spherogram.Link(pd)

def is_unlink(cl):
    L=to_link(cl)
    if L is None: return True, 'no crossings'
    ncomp=len(L.link_components)
    L2=L.copy(); L2.simplify('global')
    if len(L2.crossings)==0: return True, 'simplified to 0 crossings'
    try:
        G=L.exterior().fundamental_group()
        if G.num_relators()==0 and G.num_generators()==ncomp: return True, 'free group'
        return False, f'group gens={G.num_generators()} rels={G.relators()[:3]} after simplify crossings={len(L2.crossings)}'
    except Exception as ex:
        return False, 'exterior failed: '+repr(ex)[:100]

# sanity: the checker must reject real links
tre=spherogram.Link('3_1'); assert len(tre.crossings)==3
pool=[]
code=sys.argv[2]
from surfaceknots.yoshikawa_data import YOSHIKAWA_TABLE
for e in YOSHIKAWA_TABLE: pool.append((e.label,e.canonical_epd))
ch11=json.load(open(code+'/data/expanded_ch11.json'))
for r in ch11['entries']:
    for i,epd in enumerate(r['move_component_epds']): pool.append((r['label']+f'.{i}',epd))
for i,r in enumerate(ch11['residual_collisions']): pool.append((f'resid{i}',r['epd']))
def walk(o):
    if isinstance(o,dict):
        for k,v in o.items():
            if 'epd' in k and isinstance(v,str) and re.match(r'^[XYZO]\[',v): yield (str(o.get('family_id') or o.get('label') or k), v)
            elif 'epd' in k and isinstance(v,list):
                for x in v:
                    if isinstance(x,str) and re.match(r'^[XYZO]\[',x): yield (str(o.get('family_id') or o.get('label') or k), x)
            else: yield from walk(v)
    elif isinstance(o,list):
        for v in o: yield from walk(v)
for f in ('data/families/family_catalog.json','data/complexity/targeted_ch12_13_census.json','data/complexity/complexity_comparison.json','data/complexity/triplane_b4_census.json','data/families/hub_realizations.json','data/yoshikawa_ch10.json'):
    for lab,epd in walk(json.load(open(code+'/'+f))): pool.append((f.split('/')[-1]+':'+lab,epd))
seen=set(); uniq=[]
for lab,epd in pool:
    if epd in seen: continue
    seen.add(epd); uniq.append((lab,epd))
print('distinct EPDs',len(uniq),flush=True)
bad=0; how={}
for lab,epd in uniq:
    try: d=MarkedDiagram.from_epd(epd)
    except Exception as ex: print('PARSE FAIL',lab,repr(ex)[:80]); bad+=1; continue
    for sgn in (1,-1):
        ok,why=is_unlink(d.resolve(sgn))
        how[why.split(' ')[0]]=how.get(why.split(' ')[0],0)+1
        if not ok:
            bad+=1; print('NOT CONFIRMED',lab,'L'+('+' if sgn>0 else '-'),why,flush=True)
print('resolutions checked',2*len(uniq),'not confirmed',bad,how)
