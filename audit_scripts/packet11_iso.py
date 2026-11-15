"""Is the group of R(b a b^-1 a^-1) isomorphic to that of R(b a b^-1 a) by a map of meridians?

G(w) = < a, b | a = w b w^-1 >.  Search for a free-group automorphism phi of F(a,b) sending each of
a, b to a conjugate of a^e or b^e (e = +1 or -1, the same for both) such that phi(r1) is conjugate in
the free group to r2 or r2^-1.  Such a phi induces an isomorphism G1 -> G2 carrying meridians to
meridians (e = +1) or to inverse meridians (e = -1).  Everything is checked in the free group, so
no word problem is needed.  Basis test: (x, y) is a basis of F2 iff [x, y] is conjugate to [a, b]^+-1
(Nielsen).
"""
from itertools import product
def red(w):
    o=[]
    for x in w:
        if o and o[-1]==-x: o.pop()
        else: o.append(x)
    return tuple(o)
def inv(w): return tuple(-x for x in reversed(w))
def cyc(w):
    w=list(red(w))
    while len(w)>1 and w[0]==-w[-1]: w=w[1:-1]
    return tuple(w)
def conj_class(w):
    w=cyc(w); return min(w[i:]+w[:i] for i in range(len(w))) if w else ()
def sub(ia, ib, w):
    o=[]
    for l in w:
        p = ia if abs(l)==1 else ib
        o.extend(p if l>0 else inv(p))
    return red(o)
def rel(w): return red((-1,)+w+(2,)+inv(w))
w1=(2,1,-2,-1); w2=(2,1,-2,1)
r1, r2 = rel(w1), rel(w2)
targets={conj_class(r2):'r2', conj_class(inv(r2)):'r2^-1'}
comm=(1,2,-1,-2); basis_ok={conj_class(comm), conj_class(inv(comm))}
words=[()]
for L in range(1,7):
    words += [p for p in product((1,-1,2,-2),repeat=L) if red(p)==p]
name=lambda w: ' '.join({1:'a',-1:'a^-1',2:'b',-2:'b^-1'}[x] for x in w) or '1'
hits=[]
for e in (1,-1):
    for swap in (False, True):
        for u in words:
            ia = (2*e,) if swap else (1*e,)
            ib = red(u + ((1*e,) if swap else (2*e,)) + inv(u))
            if conj_class(sub(ia,ib,comm)) not in basis_ok: continue
            c = conj_class(sub(ia,ib,r1))
            if c in targets: hits.append((e,swap,u,targets[c]))
print('r1 =',name(r1)); print('r2 =',name(r2))
print(len(hits),'automorphisms of F2 found carrying r1 to a conjugate of r2^+-1')
for e,swap,u,t in hits[:8]:
    print('  a ->', ('b' if swap else 'a')+('' if e==1 else '^-1'), ';  b -> u', ('a' if swap else 'b')+('' if e==1 else '^-1'), 'u^-1 with u =', name(u), ';  image of r1 ~', t)
# control: the same search from r1 to r1 with e = -1 (is R(w1) invertible by such a map?) and r1 -> r2 with e=+1
print('meridian-preserving (e=+1):', sum(1 for h in hits if h[0]==1), ' meridian-inverting (e=-1):', sum(1 for h in hits if h[0]==-1))
