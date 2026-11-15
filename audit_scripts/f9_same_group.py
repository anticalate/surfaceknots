"""The two ch=12 spheres with Alexander ideal (3, x^2+1) have the same group, by a map of meridians.

GAP's simplified presentations of the two groups (independent Wirtinger presentations, f9_knot.py):
  G1 = < a, c | (c^-1 a^-1)^3,  c a c a^-1 c^-1 a c a c^-1 a^-1 >            (a = f1, c = f3)
  G2 = < a, b | (a b)^3,        b a^2 b a^-1 b^-1 a^-1 b a^2 >               (a = f1, b = f6)
Claim: with p = a c (resp. p = a b) both are  H = < a, p | p^3, a p a = p a p a p >.
Proof by rewriting: each relator of one presentation, rewritten in the generators of the other, is
reduced to the empty word by replacing a subword u by v^-1 where u v is a cyclic rotation of a relator
or its inverse and |u| >= |v|.  Every such step is valid in the group."""
from collections import deque
def red(w):
    o=[]
    for x in w:
        if o and o[-1]==-x: o.pop()
        else: o.append(x)
    return tuple(o)
def inv(w): return tuple(-x for x in reversed(w))
def sub(m, w):
    o=[]
    for l in w:
        p=m[abs(l)]; o.extend(p if l>0 else inv(p))
    return red(o)
def rules(rels):
    out={}
    for r in rels:
        for rr in (r, inv(r)):
            n=len(rr)
            for i in range(n):
                rot=rr[i:]+rr[:i]
                for k in range((n+1)//2, n+1): out.setdefault(rot[:k], set()).add(inv(rot[k:]))
    return out
def trivial(w, rels, limit=300000, slack=4):
    R=rules(rels); maxk=max(len(r) for r in rels); mink=min((len(r)+1)//2 for r in rels); w=red(w)
    if not w: return 0
    seen={w}; q=deque([(w,0)]); bound=len(w)+slack
    while q and len(seen)<limit:
        cur,dist=q.popleft(); n=len(cur)
        for i in range(n):
            for k in range(mink, min(maxk,n-i)+1):
                for v in R.get(cur[i:i+k],()):
                    nw=red(cur[:i]+v+cur[i+k:])
                    if not nw: return dist+1
                    if len(nw)<=bound and nw not in seen: seen.add(nw); q.append((nw,dist+1))
    return None
A,X=1,2   # generator 1 = a ; generator 2 = c, b or p depending on the presentation
G1=[red((-2,-1)*3), (2,1,2,-1,-2,1,2,1,-2,-1)]
G2=[red((1,2)*3),   (2,1,1,2,-1,-2,-1,2,1,1)]
H =[(2,2,2), red((1,2,1)+inv((2,1,2,1,2)))]
to_H   ={1:(1,), 2:(-1,2)}      # c or b  ->  a^-1 p
from_H ={1:(1,), 2:(1,2)}       # p       ->  a c  or  a b
ok=True
for name,G in (("G1",G1),("G2",G2)):
    for i,r in enumerate(G):
        res=trivial(sub(to_H,r),H); ok&=res is not None
        print(f"relator {i+1} of {name} is trivial in H:", "PROVED in %d rewrites"%res if res is not None else "not settled")
    for i,r in enumerate(H):
        res=trivial(sub(from_H,r),G); ok&=res is not None
        print(f"relator {i+1} of H is trivial in {name}:", "PROVED in %d rewrites"%res if res is not None else "not settled")
# the two substitutions are mutually inverse in the free group: p -> a c -> a a^-1 p = p ; c -> a^-1 p -> a^-1 a c = c
assert sub(to_H, sub(from_H,(2,)))==(2,) and sub(from_H, sub(to_H,(2,)))==(2,)
print("substitutions are mutually inverse in the free group: True")
print("ALL PROVED: G1 = H = G2, with the meridian a going to a" if ok else "INCOMPLETE")
