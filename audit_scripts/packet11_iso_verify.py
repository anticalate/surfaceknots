"""Check the isomorphism G(R(b a b^-1 a^-1)) -> G(R(b a b^-1 a)) stated in the third-party audit:
    phi(a) = a^-1, phi(b) = (b^-1 a) b^-1 (a^-1 b);   psi(a) = a^-1, psi(b) = (a^-1 b a) b^-1 (a^-1 b^-1 a).
Each of the six identities (phi(r1) = 1 in G2, psi(r2) = 1 in G1, psi phi = id on a, b, phi psi = id on
a, b) is proved by an explicit chain of rewrites: a subword u of the current word is replaced by v^-1
where u v is a cyclic rotation of the relator or its inverse and |u| >= |v|.  Each rewrite is valid in
the group, so reaching the empty word is a proof.  Breadth-first search with a visited set.
"""
import sys
from collections import deque
def red(w):
    o=[]
    for x in w:
        if o and o[-1]==-x: o.pop()
        else: o.append(x)
    return tuple(o)
def inv(w): return tuple(-x for x in reversed(w))
def sub(ia, ib, w):
    o=[]
    for l in w:
        p = ia if abs(l)==1 else ib
        o.extend(p if l>0 else inv(p))
    return red(o)
def rel(w): return red((-1,)+w+(2,)+inv(w))
r1=rel((2,1,-2,-1)); r2=rel((2,1,-2,1))
phi=((-1,), red((-2,1,-2,-1,2)))
psi=((-1,), red((-1,2,1,-2,-1,-2,1)))
def rules(r):
    out={}
    for rr in (r, inv(r)):
        n=len(rr)
        for i in range(n):
            rot=rr[i:]+rr[:i]
            for k in range((n+1)//2, n+1):
                out.setdefault(rot[:k], set()).add(inv(rot[k:]))
    return out
def trivial(w, r, limit=400000, slack=2):
    R=rules(r); maxk=len(r); w=red(w)
    if not w: return 0
    seen={w}; q=deque([(w,0)]); bound=len(w)+slack
    while q and len(seen)<limit:
        cur,dist=q.popleft()
        n=len(cur)
        for i in range(n):
            for k in range((maxk+1)//2, min(maxk,n-i)+1):
                u=cur[i:i+k]
                for v in R.get(u,()):
                    nw=red(cur[:i]+v+cur[i+k:])
                    if not nw: return dist+1
                    if len(nw)<=bound and nw not in seen:
                        seen.add(nw); q.append((nw,dist+1))
    return None
name=lambda w: ' '.join({1:'a',-1:'a^-1',2:'b',-2:'b^-1'}[x] for x in w) or '1'
comp=lambda f,g,w: sub(*g, sub(*f, w))     # g after f
checks=[('phi(r1) = 1 in G2', sub(*phi,r1), r2),
        ('psi(r2) = 1 in G1', sub(*psi,r2), r1),
        ('psi(phi(a)) = a in G1', red(comp(phi,psi,(1,))+(-1,)), r1),
        ('psi(phi(b)) = b in G1', red(comp(phi,psi,(2,))+(-2,)), r1),
        ('phi(psi(a)) = a in G2', red(comp(psi,phi,(1,))+(-1,)), r2),
        ('phi(psi(b)) = b in G2', red(comp(psi,phi,(2,))+(-2,)), r2)]
ok=True
for nm,w,r in checks:
    res=trivial(w,r)
    print(f'{nm:24s} word length {len(w):3d}:', 'PROVED in %d rewrites'%res if res is not None else 'not settled by this search', flush=True)
    ok &= res is not None
print('ALL SIX PROVED' if ok else 'INCOMPLETE')
