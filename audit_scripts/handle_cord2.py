"""Normalised cord classes for 'spun trefoil + 1-handle' diagrams.  See handle_cord.py for the set-up.
Normal form: conjugate so that the meridian at the near attaching arc is a = [[1,1],[0,1]] (exponent +1).
For W = ([[p,q],[r,s]], e):  a^i W a^j = ([[p+ir, *],[r, s+jr]], e+i+j).  Complete double-coset invariants:
  r != 0 : (r, p mod |r|, s mod |r|, K) with K = e*r - p - s          r == 0 : W = (+-a^k, e)  -> ('central', sign, e-k)
Classes are reported up to inversion (Boyle: h(g) = h(g^-1)) and up to the mirror automorphism a->a^-1, b->b^-1."""
import json, sys
sys.path.insert(0, sys.argv[1]); sys.path.insert(0, sys.argv[3])
from surfaceknots.diagrams import MarkedDiagram
from handle_cord import *

def conj_to_a(mu):
    """Return C in SL(2,Z) and flag psi such that C M C^-1 = a after optionally applying psi (a->a^-1,b->b^-1)."""
    M,e=mu
    psi = (e==-1)
    if psi: M=((M[0][0],-M[0][1]),(-M[1][0],M[1][1]))
    # M parabolic, trace 2: fixed primitive vector v with M v = v
    x,y = M[0][1], 1-M[0][0]
    if x==0 and y==0: x,y = 1-M[1][1], M[1][0]
    from math import gcd
    g=gcd(x,y); x//=g; y//=g
    # C with C v = e1:  C = [[s,t],[-y,x]] with s*x+t*y = 1
    def egcd(a_,b_):
        if b_==0: return (1,0) if a_>0 else (-1,0)
        u,v=egcd(b_,a_%b_); return (v,u-(a_//b_)*v)
    s,t=egcd(x,y)
    if s*x+t*y!=1: s,t=-s,-t
    C=((s,t),(-y,x))
    R=mm(mm(C,M),mi(C))
    if R==((1,-1),(0,1)):                      # got a^-1: flip with diag(1,-1)... use v -> -v
        C=((-s,-t),(y,-x)); R=mm(mm(C,M),mi(C))
    if R!=((1,1),(0,1)):
        # conjugate by diag(-1,1)-type is not in SL2; fall back: reflect orientation
        C2=((s,t),(-y,x)); R2=mm(mm(C2,M),mi(C2))
        raise ValueError(f'normalisation failed {R}')
    return C,psi

def key_of(Wm,e):
    (p,q),(r,s)=Wm
    if r==0:
        sign=1 if p==1 else -1; k=q*sign      # W = sign * a^k
        return ('r=0', sign, e-k)
    if r<0: pass
    m=abs(r)
    return (r, p%m, s%m, e*r-p-s)

def norm_class(W,mu):
    C,psi=conj_to_a(mu)
    Wm,e=W
    if psi: Wm=((Wm[0][0],-Wm[0][1]),(-Wm[1][0],Wm[1][1])); e=-e
    Wn=mm(mm(C,Wm),mi(C))
    ks=[]
    for (X,ex) in ((Wn,e),(mi(Wn),-e)):
        ks.append(key_of(X,ex))
        Xp=((X[0][0],-X[0][1]),(-X[1][0],X[1][1]))      # mirror automorphism psi, then re-normalise mu=a^-1 -> a via psi again is identity; use D-conjugation
        # psi(W) relative to meridian psi(a)=a^-1; conjugating everything by the orientation flip gives key of (D X D, -ex) w.r.t. a
        ks.append(key_of(Xp,-ex))
    trivial = (Wn[1][0]==0 and Wn==gpow(a,e)[0])
    return min(ks, key=repr), trivial

def analyse2(d):
    from tube_opts_lib import options
    out=[]; n=d.num_vertices
    for (u,w,tu,tw,oi,cand) in options(d):
        if oi!=0: continue
        st=cand.surface_topology()
        if st.num_components!=1 or st.components[0].homeomorphism_type!='S^2': continue
        rep=find_rep(cand)
        if rep is None: out.append((u,w,None,None,'base not B_3')); continue
        eid,rels,asg=rep
        remap=[v-(1 if v>u else 0)-(1 if v>w else 0) for v in range(n)]
        def far(x):
            y=d.opp[x]
            return None if (y>>2) in (u,w) else (remap[y>>2]<<2)|(y&3)
        A=[far(x) for x in range(4*u,4*u+4) if (d.opp[x]>>2)!=w]; B=[far(x) for x in range(4*w,4*w+4) if (d.opp[x]>>2)!=u]
        A=[x for x in A if x is not None]; B=[x for x in B if x is not None]
        if not A or not B: continue
        W,val=transport(cand,eid,rels,asg,A[0])
        mu=val(A[0]); Wb=W[B[0]]; gB=val(B[0])
        orient = +1 if gmul(gmul(Wb,mu),ginv(Wb))==gB else -1
        assert orient==1 or gmul(gmul(Wb,ginv(mu)),ginv(Wb))==gB
        try: key,triv=norm_class(Wb,mu)
        except ValueError as ex: out.append((u,w,None,None,str(ex))); continue
        out.append((u,w,key,triv,'handle sign(g(B) vs g(A))=%+d'%orient))
    return out

if __name__=='__main__':
    code=sys.argv[2]
    from surfaceknots.yoshikawa_data import YOSHIKAWA_BY_LABEL
    pool=[('TABLE 10_1^1 = Tube(3_1)',YOSHIKAWA_BY_LABEL['10_1^1'].canonical_epd)]
    c11=json.load(open(code+'/data/expanded_ch11.json'))
    for r in c11['entries']:
        if r['label']=='11_1^-2': pool.append(('TABLE 11_1^-2',r['epd']))
    for i,r in enumerate(c11['residual_collisions']): pool.append((f"ch11 resid{i} [{r['collision_matches'][0]}]",r['epd']))
    try:
        for i,r in enumerate(json.load(open(sys.argv[3]+'/tube_ch10.json'))):
            if len(r['types'])==1: pool.append((f"blind ch10 #{i} {r['types'][0]} S3={r['hom']['S3']}",r['epd']))
        for i,r in enumerate(json.load(open(sys.argv[3]+'/tube_ch9.json'))):
            if len(r['types'])==1: pool.append((f"blind ch9 #{i} {r['types'][0]} S3={r['hom']['S3']}",r['epd']))
    except Exception as ex: print('no blind files',ex)
    fam=json.load(open(code+'/data/families/family_catalog.json'))
    def walk(o):
        if isinstance(o,dict):
            if isinstance(o.get('epd'),str): yield (str(o.get('family_id')),o['epd'])
            for v in o.values(): yield from walk(v)
        elif isinstance(o,list):
            for v in o: yield from walk(v)
    for lab,epd in walk(fam):
        if lab in ('Tube(3_1)',): pool.append(('FAMILY '+lab,epd))
    for lab,epd in pool:
        d=MarkedDiagram.from_epd(epd)
        if d.surface_topology().num_components!=1: continue
        for (u,w,key,triv,note) in analyse2(d):
            print(f"{lab:55s} bigon({u},{w}) class={key} trivial={triv} {note}",flush=True)
