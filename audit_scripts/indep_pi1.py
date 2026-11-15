# Independent pi_1 check.  Uses ONLY the repo's EPD parser (MarkedDiagram.from_epd -> opp, vtype);
# the Wirtinger-with-saddles presentation and the Hom counting are written from scratch here.
import json, sys, itertools
sys.path.insert(0, sys.argv[1])
from surfaceknots.diagrams import MarkedDiagram
from sympy.combinatorics import Permutation, PermutationGroup
from sympy.combinatorics.named_groups import SymmetricGroup, AlternatingGroup

def group_table(G):
    els=list(G.elements); idx={e:i for i,e in enumerate(els)}
    n=len(els)
    mul=[[idx[els[i]*els[j]] for j in range(n)] for i in range(n)]
    inv=[idx[els[i]**-1] for i in range(n)]
    return n,mul,inv
Q8=PermutationGroup([Permutation(0,1,2,3)(4,5,6,7), Permutation(0,4,2,6)(1,7,3,5)])
GROUPS={'S3':group_table(SymmetricGroup(3)),'Q8':group_table(Q8),'A4':group_table(AlternatingGroup(4)),'S4':group_table(SymmetricGroup(4))}
assert GROUPS['Q8'][0]==8

def relations(d):
    """Signed generators: dart x -> (edge_id, sign); g(x) = meridian for edge oriented INTO its vertex.
    Returns (num_edges, rels) with rels as ('eq', a, b)  meaning g(a) = g(b)^-1  (a,b darts)
    or ('conj', out, inn, over, eps) meaning g(out)^-1 = g(over)^eps g(inn) g(over)^-eps."""
    n=d.num_vertices; eid={}; k=0
    for x in range(4*n):
        if x not in eid:
            eid[x]=(k,1); eid[d.opp[x]]=(k,-1); k+=1
    rels=[]
    for v in range(n):
        t=d.vtype[v]; b=4*v
        if t in ('Y','Z'):
            # saddle: all four arms carry the same meridian, signs alternate around the vertex
            for s in range(3): rels.append(('eq',b+s,b+s+1))
        else:
            u=0 if t=='X0' else 1; o=1-u if False else (1 if t=='X0' else 0)
            rels.append(('eq',b+o,b+o+2))                      # over strand passes straight through
            eps=1 if ((o-u)%4)==1 else -1
            rels.append(('conj',b+u+2,b+u,b+o,eps))
    return k,eid,rels

def count_homs(d,G):
    n,mul,inv=G
    k,eid,rels=relations(d)
    if k==0: return n**d.free_loops
    def val(asg,x):
        e,s=eid[x]; a=asg[e]
        if a<0: return -1
        return a if s==1 else inv[a]
    def setv(asg,x,g,trail):
        e,s=eid[x]; g2=g if s==1 else inv[g]
        if asg[e]<0: asg[e]=g2; trail.append(e); return True
        return asg[e]==g2
    def propagate(asg,trail):
        changed=True
        while changed:
            changed=False
            for r in rels:
                if r[0]=='eq':
                    a,bb=val(asg,r[1]),val(asg,r[2])
                    if a>=0 and bb>=0:
                        if a!=inv[bb]: return False
                    elif a>=0:
                        if not setv(asg,r[2],inv[a],trail): return False
                        changed=True
                    elif bb>=0:
                        if not setv(asg,r[1],inv[bb],trail): return False
                        changed=True
                else:
                    _,out,inn,over,eps=r
                    O=val(asg,over)
                    if O<0: continue
                    if eps<0: O=inv[O]
                    a,bb=val(asg,out),val(asg,inn)
                    if bb>=0:
                        tgt=inv[mul[mul[O][bb]][inv[O]]]      # g(out) = (O g(inn) O^-1)^-1
                        if a>=0:
                            if a!=tgt: return False
                        else:
                            if not setv(asg,out,tgt,trail): return False
                            changed=True
                    elif a>=0:
                        tgt=mul[mul[inv[O]][inv[a]]][O]       # g(inn) = O^-1 g(out)^-1 O
                        if not setv(asg,inn,tgt,trail): return False
                        changed=True
        return True
    asg=[-1]*k
    def rec():
        try: e=asg.index(-1)
        except ValueError: return 1
        tot=0
        for g in range(n):
            trail=[]; asg[e]=g; trail.append(e)
            if propagate(asg,trail): tot+=rec()
            for t in trail: asg[t]=-1
        return tot
    return rec()*(n**d.free_loops)

code=sys.argv[2]
recs=[(r['label'],r['epd'],r['hom_counts']) for r in json.load(open(code+'/data/yoshikawa_ch10.json'))]
c11=json.load(open(code+'/data/expanded_ch11.json'))
recs+=[(r['label'],r['epd'],r['hom_counts']) for r in c11['entries']]
recs+=[(f"resid{i}",r['epd'],r['hom_counts']) for i,r in enumerate(c11['residual_collisions'])]
bad=0
for lab,epd,hc in recs:
    d=MarkedDiagram.from_epd(epd)
    mine={g:count_homs(d,GROUPS[g]) for g in ('S3','Q8','A4','S4')}
    theirs={g:hc[g] for g in mine}
    ok=mine==theirs
    if not ok: bad+=1
    if not ok or not lab.startswith('resid'): print(('OK  ' if ok else 'DIFF'),lab,mine,'' if ok else theirs,flush=True)
print('checked',len(recs),'mismatches',bad)
