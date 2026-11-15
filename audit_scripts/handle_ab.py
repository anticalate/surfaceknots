"""Second tube on a link whose group is Z^2 (base = 6_1^{0,1} = T(b), or a diagram of it).
For a tube with both feet on one component X of the base, the cord class in <mu_X>\Z^2/<mu_X> = Z is the
exponent of the OTHER component's meridian picked up by carrying the meridian along X from foot to foot."""
import sys, json, importlib.util
spec=importlib.util.spec_from_file_location('hf', sys.argv[3]+'/handle_free.py'); hf=importlib.util.module_from_spec(spec); spec.loader.exec_module(hf)
sys.path.insert(0, sys.argv[3])
from indep_pi1_lib import relations, GROUPS, count_homs
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.moves import list_level_yoshikawa_neighbors, canonical_marked_key
def neg(v): return (-v[0],-v[1])
def add(u,v): return (u[0]+v[0],u[1]+v[1])
def ab_rep(cand):
    k,eid,rels=relations(cand); comp,nc=cand._dart_surface_components()
    if nc!=2: return None
    ecomp={e:comp[x] for x,(e,s) in eid.items()}
    asg=[None]*k
    def val(x):
        e,s=eid[x]; v=asg[e]; return None if v is None else (v if s==1 else neg(v))
    def setv(x,g):
        e,s=eid[x]; g2=g if s==1 else neg(g)
        if asg[e] is None: asg[e]=g2; return True
        return asg[e]==g2
    for c in (0,1):
        e=next(e for e in range(k) if ecomp[e]==c); asg[e]=(1,0) if c==0 else (0,1)
    ch=True
    while ch:
        ch=False
        for r in rels:
            a,b=(r[1],r[2]) if r[0]=='eq' else (r[1],r[2])     # abelian: conj relation is also g(out) = -g(inn)
            x,y=val(a),val(b)
            if x is not None and y is not None:
                if x!=neg(y): return None
            elif x is not None:
                setv(b,neg(x)); ch=True
            elif y is not None:
                setv(a,neg(y)); ch=True
    if any(v is None for v in asg): return None
    return eid,rels,asg,comp
def transport(d,eid,rels,asg,start):
    def val(x):
        e,s=eid[x]; v=asg[e]; return v if s==1 else neg(v)
    conj={}
    for r in rels:
        if r[0]=='conj':
            _,out,inn,over,eps=r; O=val(over); O=O if eps>0 else neg(O)
            conj[inn]=(out,O); conj[out]=(inn,neg(O))
    W={start:(0,0)}; q=[start]
    while q:
        x=q.pop(0); w=W[x]; v=x>>2; s=x&3; nb=[(d.opp[x],w)]
        if d.vtype[v] in ('Y','Z'): nb+=[(4*v+t,w) for t in range(4) if t!=s]
        elif x in conj:
            y,O=conj[x]; nb.append((y,add(O,w)))
        else: nb.append((x^2,w))
        for y,wy in nb:
            if y not in W: W[y]=wy; q.append(y)
    return W
def analyse(d):
    out=[]
    for u,w,cand,dA,dB in hf.bigon_bases(d):
        cst=cand.surface_topology()
        if cst.num_components!=2: continue
        types=[c.homeomorphism_type for c in cst.components]
        hom=tuple(count_homs(cand,GROUPS[g]) for g in ('S3','Q8','A4','S4'))
        if hom!=(18,40,48,120): out.append(f"base {types} has hom counts {hom}: not Z^2"); continue
        rep=ab_rep(cand)
        if rep is None: out.append('no abelian rep'); continue
        eid,rels,asg,comp=rep
        cA,cB=comp[dA],comp[dB]
        if cA!=cB: out.append('tube joins the two components'); continue
        W=transport(cand,eid,rels,asg,dA)
        k=W[dB][1-cA]          # exponent of the other component's meridian
        res=d.surface_topology(); 
        out.append(f"base {'+'.join(sorted(types))} (n={cand.num_vertices}, group Z^2); tube on the {types[cA]} component; cord class k = {abs(k)}")
    return out
code=sys.argv[2]
from surfaceknots.yoshikawa_data import YOSHIKAWA_BY_LABEL
pool=[(l,YOSHIKAWA_BY_LABEL[l].canonical_epd) for l in ('8_1^1,1','10_1^1,1','9_1^1,-2','10_1^-2,-2')]
c11={r['label']:r['epd'] for r in json.load(open(code+'/data/expanded_ch11.json'))['entries']}
pool.append(('11_1^1,-2',c11['11_1^1,-2']))
for lab,epd in pool:
    d0=MarkedDiagram.from_epd(epd); seen={canonical_marked_key(d0):d0}; q=[d0]; i=0; res=[]
    while i<len(q) and i<200:
        cur=q[i]; i+=1; res+=analyse(cur)
        for nb in list_level_yoshikawa_neighbors(cur):
            kk=canonical_marked_key(nb)
            if kk not in seen: seen[kk]=nb; q.append(nb)
    import collections
    print(lab,'| diagrams examined',i); 
    for r,n in collections.Counter(res).most_common(): print('     ',n,'x',r)
