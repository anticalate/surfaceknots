"""Prototype: classify 1-handles visible as marked bigons in a ch-diagram, when the base is the spun trefoil.

d = base + 1-handle (marked bigon u,w).  G = pi_1(S^4 \ base) = B_3 (trefoil group).  The handle's cord class
is W in <mu>\G/<mu>, where g(B) = W g(A)^{+-1} W^{-1} is the meridian at the far attaching arc expressed by
transport ALONG THE SURFACE from the near attaching arc.  Handle trivial  <=>  W in <mu>   (Boyle 1988).
Faithful representation of B_3:  w -> (image in SL(2,Z), exponent sum); kernel of B_3 -> SL(2,Z) is <Delta^4>,
which has exponent sum 12, so the pair is injective.
"""
import json, sys, itertools
sys.path.insert(0, sys.argv[1]); sys.path.insert(0, sys.argv[3])
from surfaceknots.diagrams import MarkedDiagram
from indep_pi1_lib import relations

def mm(A,B): return ((A[0][0]*B[0][0]+A[0][1]*B[1][0], A[0][0]*B[0][1]+A[0][1]*B[1][1]),(A[1][0]*B[0][0]+A[1][1]*B[1][0], A[1][0]*B[0][1]+A[1][1]*B[1][1]))
def mi(A): return ((A[1][1],-A[0][1]),(-A[1][0],A[0][0]))
def gmul(x,y): return (mm(x[0],y[0]), x[1]+y[1])
def ginv(x): return (mi(x[0]), -x[1])
I2=((1,0),(0,1)); ONE=(I2,0)
a=(((1,1),(0,1)),1); b=(((1,0),(-1,1)),1)
assert gmul(gmul(a,b),a)==gmul(gmul(b,a),b)
def gpow(x,k):
    r=ONE; y=x if k>=0 else ginv(x)
    for _ in range(abs(k)): r=gmul(r,y)
    return r

def find_rep(d):
    """Assign (matrix, exp) to every edge generator consistently with all relations, seeds a^{+-1}, b^{+-1}."""
    k,eid,rels=relations(d)
    def val(asg,x):
        e,s=eid[x]; v=asg[e]
        if v is None: return None
        return v if s==1 else ginv(v)
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
    for e1 in range(k):
        for e2 in range(k):
            if e1==e2: continue
            for s2 in (b,ginv(b)):
                asg=[None]*k; asg[e1]=a; asg[e2]=s2
                if propagate(asg) and all(v is not None for v in asg):
                    return eid,rels,asg
    return None

def transport(d,eid,rels,asg,start):
    """BFS along the surface from dart `start`.  Returns dict dart -> W with g(dart) = W g(start)^{+-1} W^{-1}."""
    def val(x):
        e,s=eid[x]; v=asg[e]; return v if s==1 else ginv(v)
    n=d.num_vertices
    conj={}
    for r in rels:
        if r[0]=='conj':
            _,out,inn,over,eps=r
            O=val(over); O=O if eps>0 else ginv(O)
            conj[inn]=(out,O); conj[out]=(inn,ginv(O))
    W={start:ONE}; q=[start]
    while q:
        x=q.pop(0); w=W[x]; v=x>>2; s=x&3
        nbrs=[(d.opp[x],w)]
        if d.vtype[v] in ('Y','Z'):
            nbrs+=[(4*v+t,w) for t in range(4) if t!=s]
        elif x in conj:
            y,O=conj[x]; nbrs.append((y,gmul(O,w)))
        else:
            nbrs.append((x^2,w))
        for y,wy in nbrs:
            if y not in W: W[y]=wy; q.append(y)
    return W,val

def analyse(d, bigon_filter=None):
    from tube_opts_lib import options
    out=[]
    n=d.num_vertices
    for (u,w,tu,tw,oi,cand) in options(d):
        if oi!=0: continue
        st=cand.surface_topology()
        if st.num_components!=1 or st.components[0].homeomorphism_type!='S^2': continue
        rep=find_rep(cand)
        if rep is None: out.append((u,w,'base group is not B_3 (no rep found)')); continue
        eid,rels,asg=rep
        # locate attaching edges A (through u) and B (through w) in cand
        remap=[v-(1 if v>u else 0)-(1 if v>w else 0) for v in range(n)]
        def outer(v):
            return [x for x in range(4*v,4*v+4) if (d.opp[x]>>2) not in (u,w) or (d.opp[x]>>2)==v]
        def far(x):
            y=d.opp[x]
            return None if (y>>2) in (u,w) else (remap[y>>2]<<2)|(y&3)
        A=[far(x) for x in range(4*u,4*u+4) if (d.opp[x]>>2)!=w]; B=[far(x) for x in range(4*w,4*w+4) if (d.opp[x]>>2)!=u]
        A=[x for x in A if x is not None]; B=[x for x in B if x is not None]
        if not A or not B: out.append((u,w,'degenerate attaching arcs')); continue
        W,val=transport(cand,eid,rels,asg,A[0])
        mu=val(A[0]); Wb=W[B[0]]
        gB=val(B[0])
        sign=None
        for sgn,m in ((1,mu),(-1,ginv(mu))):
            if gmul(gmul(Wb,m),ginv(Wb))==gB: sign=sgn
        assert sign is not None, 'transport inconsistent'
        k=Wb[1]*mu[1]            # mu has exponent +-1
        trivial = (gpow(mu,k)==Wb)
        ori = d.surface_topology().components[0].homeomorphism_type
        out.append((u,w,f"result={ori} cord exponent={k} cord trivial={trivial} W={Wb[0]}"))
    return out

if __name__=='__main__':
    code=sys.argv[2]
    from surfaceknots.yoshikawa_data import YOSHIKAWA_BY_LABEL
    pool=[('10_1^1 (Tube(3_1), expect NON-trivial)',YOSHIKAWA_BY_LABEL['10_1^1'].canonical_epd)]
    c11=json.load(open(code+'/data/expanded_ch11.json'))
    for r in c11['entries']:
        if r['label']=='11_1^-2': pool.append(('11_1^-2 (expect NON-trivial)',r['epd']))
    for i,r in enumerate(c11['residual_collisions']): pool.append((f"resid{i}: {r['collision_matches'][0]}",r['epd']))
    for lab,epd in pool:
        d=MarkedDiagram.from_epd(epd)
        if d.surface_topology().num_components!=1: continue
        res=analyse(d)
        for r in res: print(lab,'| bigon',r[0],r[1],'|',r[2],flush=True)
