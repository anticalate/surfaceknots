"""Independent pi_1 of surface-link complement from a marked graph diagram.
Generators: one per EDGE of D (oriented from the smaller dart to the larger dart).
m(d) := meridian of the strand oriented OUTWARD from the vertex at dart d
      = x_e if d is the smaller dart of e (edge leaves this vertex), else x_e^{-1}.
Relations:
 classical crossing, over slots (o, o+2): strand through => m(o) * m(o+2) = 1
   under slots (u, u+2): orient under u -> u+2, over o -> o+2 (o in), incoming meridian X = m(u)^{-1},
   outgoing Y = m(u+2), over Z = m(o+2).   sign s = cross(over_dir, under_dir) ; Y = Z^{s} X Z^{-s}
 marked vertex: adjacent darts a, b = a+1:  m(a) * m(b) = 1  (alternating in/out <-> same meridian)
"""
import itertools, sys
from sympy.combinatorics.free_groups import free_group
from sympy.combinatorics.fp_groups import FpGroup

DIRV = [(1,0),(0,1),(-1,0),(0,-1)]  # slot 0=E,1=N,2=W,3=S (CCW)

def presentation(d):
    n = d.num_vertices
    eid = {}
    ne = 0
    for hd in range(4*n):
        od = d.opp[hd]
        if hd not in eid:
            eid[hd] = ne; eid[od] = ne; ne += 1
    def m(hd):  # literal (1-based signed)
        od = d.opp[hd]
        e = eid[hd] + 1
        if hd == od:
            raise ValueError
        return e if hd < od else -e
    rels = []
    for v in range(n):
        vt = d.vtype[v]
        if vt in ("X0","X1"):
            us = (0,2) if vt == "X0" else (1,3)
            os_ = (1,3) if vt == "X0" else (0,2)
            u_in, u_out = us; o_in, o_out = os_
            rels.append((m(4*v+o_in), m(4*v+o_out)))
            # directions of travel: under goes from u_in to u_out => direction = DIRV[u_out]; over direction = DIRV[o_out]
            ud = DIRV[u_out]; od_ = DIRV[o_out]
            s = od_[0]*ud[1] - od_[1]*ud[0]
            X = -m(4*v+u_in); Y = m(4*v+u_out); Z = m(4*v+o_out)
            if s == 1:
                rels.append((Z, X, -Z, -Y))
            else:
                rels.append((-Z, X, Z, -Y))
        else:
            for a in range(3):
                rels.append((m(4*v+a), m(4*v+a+1)))
    ngen = ne + d.free_loops
    return ngen, rels

def perm_group(gens):
    deg = len(gens[0]); ident = tuple(range(deg))
    elems=[ident]; seen={ident:0}; i=0
    while i < len(elems):
        p = elems[i]; i+=1
        for g in gens:
            c = tuple(p[g[k]] for k in range(deg))
            if c not in seen:
                seen[c]=len(elems); elems.append(c)
    return elems, seen

def mulp(a,b): return tuple(a[b[k]] for k in range(len(a)))
def invp(a):
    r=[0]*len(a)
    for i,x in enumerate(a): r[x]=i
    return tuple(r)

GROUPS = {
 "S3": [(1,2,0),(1,0,2)],
 "S4": [(1,2,3,0),(1,0,2,3)],
 "A4": [(1,2,0,3),(0,2,3,1)],
 "A5": [(1,2,3,4,0),(1,2,0,3,4)],
 "D8": [(1,2,3,0),(0,3,2,1)],
 "D10": [(1,2,3,4,0),(0,4,3,2,1)],
 "F21": [(1,2,3,4,5,6,0), tuple((2*i)%7 for i in range(7))],
 "Q8": [(1,2,3,0,5,6,7,4),(4,7,6,5,2,1,0,3)],  # i,j acting regularly (checked below)
}

def hom_count(ngen, rels, gname):
    """Brute force with union-find style propagation: eliminate generators using length-2 relators first."""
    elems, seen = perm_group(GROUPS[gname])
    # union-find with sign on generators through length-2 relators (a,b): x_a^{sa} x_b^{sb} = 1
    parent = list(range(ngen+1)); sgn=[1]*(ngen+1)
    def find(x):
        s=1
        while parent[x]!=x:
            s*=sgn[x]; x=parent[x]
        return x,s
    extra=[]
    for r in rels:
        if len(r)==2:
            a,b = r
            ra,sa = find(abs(a)); rb,sb = find(abs(b))
            # x_a^{sign a} = x_b^{-sign b};  x_a = root_a^{sa}
            ea = (1 if a>0 else -1)*sa; eb = -(1 if b>0 else -1)*sb
            if ra!=rb:
                parent[ra]=rb; sgn[ra]= ea*eb   # root_a^{ea} = root_b^{eb} -> root_a = root_b^{ea*eb}
            else:
                if ea!=eb: extra.append((ra,ra))  # root^2 = 1
        else:
            extra.append(r)
    roots = sorted({find(g)[0] for g in range(1,ngen+1)})
    ridx = {r:i for i,r in enumerate(roots)}
    def tr(l):
        r,s = find(abs(l)); s *= (1 if l>0 else -1)
        return (ridx[r]+1)*s
    red = [tuple(tr(l) for l in r) for r in extra]
    k = len(roots)
    # order vars by first appearance; backtrack
    by_max = [[] for _ in range(k)]
    for r in red:
        by_max[max(abs(l)-1 for l in r)].append(r)
    N=len(elems); inv=[seen[invp(e)] for e in elems]
    mul=[[seen[mulp(a,b)] for b in elems] for a in elems]
    assign=[0]*k
    def ok(i):
        for r in by_max[i]:
            acc=0
            for l in r:
                g = assign[abs(l)-1]
                acc = mul[acc][g if l>0 else inv[g]]
            if acc!=0: return False
        return True
    def bt(i):
        if i==k: return 1
        t=0
        for val in range(N):
            assign[i]=val
            if ok(i): t+=bt(i+1)
        return t
    return bt(0)

def reduced(ngen, rels):
    """eliminate generators via length-2 relators and one-occurrence relators; return (k, rels)"""
    from surfaceknots.invariants.fundamental_group import GroupPresentation, tietze_simplify
    # NOTE: uses repo tietze only as a convenience for sympy coset enumeration speed; hom counts above don't.
    p = tietze_simplify(GroupPresentation(ngen, tuple(rels), ()))
    return p.num_generators, p.relators

def fp_group(k, rels):
    F, *gens = free_group(",".join(f"x{i}" for i in range(k)))
    rr=[]
    for r in rels:
        w = F.identity
        for l in r:
            w = w * (gens[abs(l)-1] if l>0 else gens[abs(l)-1]**-1)
        rr.append(w)
    return FpGroup(F, rr), gens

def hom_count2(ngen, rels, gname, _cache={}):
    """Independent counter with unit propagation (deduce a generator when it is the only unknown, occurring once)."""
    if gname not in _cache:
        elems, seen = perm_group(GROUPS[gname])
        inv=[seen[invp(e)] for e in elems]
        mul=[[seen[mulp(a,b)] for b in elems] for a in elems]
        _cache[gname]=(len(elems),mul,inv)
    N,mul,inv=_cache[gname]
    rels=[tuple(r) for r in rels]
    relg=[[abs(l)-1 for l in r] for r in rels]
    assign=[-1]*ngen
    def ev(l):
        g=assign[abs(l)-1]
        return g if l>0 else inv[g]
    def propagate(changed):
        while True:
            prog=False
            for r,rg in zip(rels,relg):
                un=[k for k,g in enumerate(rg) if assign[g]==-1]
                if not un:
                    acc=0
                    for l in r: acc=mul[acc][ev(l)]
                    if acc!=0: return False
                elif len(un)==1:
                    pos=un[0]
                    pre=0
                    for l in r[:pos]: pre=mul[pre][ev(l)]
                    suf=0
                    for l in r[pos+1:]: suf=mul[suf][ev(l)]
                    val=mul[inv[pre]][inv[suf]]
                    assign[rg[pos]]= val if r[pos]>0 else inv[val]
                    changed.append(rg[pos]); prog=True
            if not prog: return True
    def bt():
        changed=[]
        if not propagate(changed):
            for g in changed: assign[g]=-1
            return 0
        rem=[g for g in range(ngen) if assign[g]==-1]
        if not rem:
            for g in changed: assign[g]=-1
            return 1
        best=None;bs=None
        for g in rem:
            sc=0
            for rg in relg:
                if g in rg:
                    u=len({x for x in rg if assign[x]==-1})
                    sc+= 10 if u==2 else 1
            if bs is None or sc>bs: bs=sc;best=g
        t=0
        for val in range(N):
            assign[best]=val
            t+=bt()
        assign[best]=-1
        for g in changed: assign[g]=-1
        return t
    return bt()
