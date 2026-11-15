"""Cord classes for 'base 9_1 + 1-handle'.  pi_1(9_1) = BS(1,2) = <a,x | x a x^-1 = a^2>, faithful in the affine
group of Z[1/2]:  element (k, c) acts t -> 2^k t + c.   Meridian normalised to x = (1, 0).
x^i (k,c) x^j = (i+k+j, 2^i c), so the double coset <x> W <x> is classified by the odd part of c (0 <=> trivial)."""
import json, sys
from fractions import Fraction as Fr
sys.path.insert(0, sys.argv[1]); sys.path.insert(0, sys.argv[3])
from surfaceknots.diagrams import MarkedDiagram
from indep_pi1_lib import relations
from tube_opts_lib import options

def gmul(g,h): return (g[0]+h[0], g[1]+Fr(2)**g[0]*h[1])        # (g*h)(t) = g(h(t))
def ginv(g): return (-g[0], -g[1]/Fr(2)**g[0])
ONE=(0,Fr(0))
def oddpart(c):
    if c==0: return 0
    n,dn=c.numerator,c.denominator
    while n%2==0: n//=2
    return n          # denominator is a power of 2

def find_rep(d):
    k,eid,rels=relations(d)
    def val(asg,x):
        e,s=eid[x]; v=asg[e]
        return None if v is None else (v if s==1 else ginv(v))
    def setv(asg,x,g):
        e,s=eid[x]; g2=g if s==1 else ginv(g)
        if asg[e] is None: asg[e]=g2; return True
        return asg[e]==g2
    def propagate(asg):
        ch=True
        while ch:
            ch=False
            for r in rels:
                if r[0]=='eq':
                    x,y=val(asg,r[1]),val(asg,r[2])
                    if x is not None and y is not None:
                        if x!=ginv(y): return False
                    elif x is not None:
                        if not setv(asg,r[2],ginv(x)): return False
                        ch=True
                    elif y is not None:
                        if not setv(asg,r[1],ginv(y)): return False
                        ch=True
                else:
                    _,out,inn,over,eps=r
                    O=val(asg,over)
                    if O is None: continue
                    if eps<0: O=ginv(O)
                    x,y=val(asg,out),val(asg,inn)
                    if y is not None:
                        t=ginv(gmul(gmul(O,y),ginv(O)))
                        if x is not None:
                            if x!=t: return False
                        else:
                            if not setv(asg,out,t): return False
                            ch=True
                    elif x is not None:
                        t=gmul(gmul(ginv(O),ginv(x)),O)
                        if not setv(asg,inn,t): return False
                        ch=True
        return True
    X=(1,Fr(0))
    for e1 in range(k):
        for e2 in range(k):
            if e1==e2: continue
            for s2 in ((1,Fr(1)),(1,Fr(-1)),ginv((1,Fr(1))),ginv((1,Fr(-1)))):
                asg=[None]*k; asg[e1]=X; asg[e2]=s2
                if propagate(asg) and all(v is not None for v in asg): return eid,rels,asg
    return None

def transport(d,eid,rels,asg,start):
    def val(x):
        e,s=eid[x]; v=asg[e]; return v if s==1 else ginv(v)
    conj={}
    for r in rels:
        if r[0]=='conj':
            _,out,inn,over,eps=r
            O=val(over); O=O if eps>0 else ginv(O)
            conj[inn]=(out,O); conj[out]=(inn,ginv(O))
    W={start:ONE}; q=[start]
    while q:
        x=q.pop(0); w=W[x]; v=x>>2; s=x&3
        nb=[(d.opp[x],w)]
        if d.vtype[v] in ('Y','Z'): nb+=[(4*v+t,w) for t in range(4) if t!=s]
        elif x in conj:
            y,O=conj[x]; nb.append((y,gmul(O,w)))
        else: nb.append((x^2,w))
        for y,wy in nb:
            if y not in W: W[y]=wy; q.append(y)
    return W,val

def analyse(d):
    out=[]; n=d.num_vertices
    for (u,w,tu,tw,oi,cand) in options(d):
        if oi!=0: continue
        st=cand.surface_topology()
        if st.num_components!=1 or st.components[0].homeomorphism_type!='S^2': continue
        rep=find_rep(cand)
        if rep is None: continue
        eid,rels,asg=rep
        remap=[v-(1 if v>u else 0)-(1 if v>w else 0) for v in range(n)]
        def far(x):
            y=d.opp[x]
            return None if (y>>2) in (u,w) else (remap[y>>2]<<2)|(y&3)
        A=[far(x) for x in range(4*u,4*u+4) if (d.opp[x]>>2)!=w]; B=[far(x) for x in range(4*w,4*w+4) if (d.opp[x]>>2)!=u]
        A=[x for x in A if x is not None]; B=[x for x in B if x is not None]
        if not A or not B: continue
        W,val=transport(cand,eid,rels,asg,A[0])
        mu=val(A[0]); Wb=W[B[0]]
        # normalise mu to x=(1,0): mu=(e,c), e=+-1.  If e=-1 use mu^-1.  Conjugate by translation T=(0,t): T (1,c) T^-1 = (1, c - t) -> t=c
        m=mu if mu[0]==1 else ginv(mu)
        T=(0,m[1]); Wn=gmul(gmul(T,Wb),ginv(T))
        assert gmul(gmul(T,m),ginv(T))==(1,Fr(0))
        cls=abs(oddpart(Wn[1]))                 # up to inversion / mirror: sign dropped
        out.append((u,w,cls,d.surface_topology().components[0].homeomorphism_type))
    return out

code=sys.argv[2]
c11=json.load(open(code+'/data/expanded_ch11.json'))
pool=[(f"ch11 resid{i} [{r['collision_matches'][0]}]",r['epd']) for i,r in enumerate(c11['residual_collisions'])]
for f in ('tube_ch10.json','tube_ch9.json'):
    for i,r in enumerate(json.load(open(sys.argv[3]+'/'+f))):
        if len(r['types'])==1: pool.append((f"blind {f[5:9]} #{i} {r['types'][0]} S3={r['hom']['S3']}",r['epd']))
for lab,epd in pool:
    d=MarkedDiagram.from_epd(epd)
    if d.surface_topology().num_components!=1: continue
    for (u,w,cls,typ) in analyse(d):
        print(f"{lab:55s} bigon({u},{w}) result={typ} base=9_1 cord class |odd part|={cls}  {'TRIVIAL handle' if cls==0 else ''}",flush=True)
