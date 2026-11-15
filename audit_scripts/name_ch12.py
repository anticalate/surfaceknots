import sys, json, importlib.util, collections
spec=importlib.util.spec_from_file_location('hf', sys.argv[3]+'/handle_free.py'); hf=importlib.util.module_from_spec(spec); spec.loader.exec_module(hf)
sys.path.insert(0, sys.argv[3])
from indep_pi1_lib import GROUPS, count_homs
from surfaceknots.diagrams import MarkedDiagram
d=json.load(open(sys.argv[3]+'/ch12_raw_snapshot.json'))
print({k:(v if not isinstance(v,(dict,list)) else type(v).__name__) for k,v in d.items()})
for m,st in d['per_m_stats'].items(): print('m',m,st)
epds=d['per_m_level_epds'].get('2',[])
print('m=2 diagrams',len(epds),flush=True)
rows=[]
for i,epd in enumerate(epds):
    dg=MarkedDiagram.from_epd(epd); st=dg.surface_topology()
    typ='+'.join(sorted(c.homeomorphism_type for c in st.components))
    e=tuple(sorted(c.normal_euler for c in st.components))
    hom=tuple(count_homs(dg,GROUPS[g]) for g in ('S3','Q8','A4','S4'))
    names={}
    if typ in ('S^2','S^2+T^2','K^2+S^2'):
        names,_=hf.name_surface(dg, 40)
    rows.append((typ,e,hom,tuple(sorted(names)),epd))
    print(i,typ,e,hom,sorted(names),flush=True)
json.dump(rows,open(sys.argv[3]+'/ch12_m2_named.json','w'))
agg=collections.defaultdict(lambda: [0,set()])
for typ,e,hom,names,epd in rows:
    agg[(typ,e,hom)][0]+=1; agg[(typ,e,hom)][1].update(names)
print('\n=== cells (type, e, homs S3/Q8/A4/S4): diagrams, names ===')
for k,(n,nm) in sorted(agg.items()): print(k[0],k[1],k[2],'| diagrams',n,'|',sorted(nm) if nm else '-')
