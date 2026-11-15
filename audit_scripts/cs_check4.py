"""Finite-quotient comparison: the 4-relator ch=12 2-knot group vs the metabelian group
M = Z |x Z[t]/(t^3 - t^2 + 2t - 1) (a Cappell-Shaneson group) vs the ribbon knot 11_2."""
import sys, json, time
sys.path.insert(0,sys.argv[1]); sys.path.insert(0,sys.argv[2])
from sympy.combinatorics import Permutation, PermutationGroup
from sympy.combinatorics.named_groups import SymmetricGroup, AlternatingGroup
from surfaceknots.diagrams import MarkedDiagram
from indep_pi1_lib import count_homs
def table(G):
    els=list(G.elements); idx={e:i for i,e in enumerate(els)}; n=len(els)
    return n,[[idx[els[i]*els[j]] for j in range(n)] for i in range(n)],[idx[els[i]**-1] for i in range(n)]
def agl(q):
    g=next(g for g in range(2,q) if len({pow(g,k,q) for k in range(q-1)})==q-1)
    return PermutationGroup([Permutation([(x+1)%q for x in range(q)]), Permutation([(g*x)%q for x in range(q)])])
def psl27():
    def mob(a,b,c,d):
        img=[]
        for x in list(range(7))+[None]:
            num,den=(a,c) if x is None else ((a*x+b)%7,(c*x+d)%7)
            img.append(7 if den==0 else (num*pow(den,-1,7))%7)
        return Permutation(img)
    return PermutationGroup([mob(1,1,0,1),mob(0,-1,1,0)])
def agaml18():
    # AGammaL(1,8) = F_8 x| (F_8^* x| Gal): order 168, on the 8 points of F_8 = F_2[z]/(z^3+z+1)
    def mulz(v):  # multiply by z
        b0,b1,b2=v&1,(v>>1)&1,(v>>2)&1
        return (b2) | ((b0^b2)<<1) | (b1<<2)
    def frob(v):  # squaring
        r=0
        for a in range(8):
            pass
        # square via multiplication table
        def mul(a,b):
            r=0
            for i in range(3):
                if (b>>i)&1:
                    t=a
                    for _ in range(i): t=mulz(t)
                    r^=t
            return r
        return mul(v,v)
    return PermutationGroup([Permutation([v^1 for v in range(8)]),Permutation([mulz(v) for v in range(8)]),Permutation([frob(v) for v in range(8)])])
GROUPS={'S5':SymmetricGroup(5),'A6':AlternatingGroup(6),'PSL(2,7)':psl27(),'AGammaL(1,8)':agaml18(),'AGL(1,7)':agl(7),'AGL(1,11)':agl(11),'AGL(1,13)':agl(13),'AGL(1,29)':agl(29),'AGL(1,43)':agl(43)}
def homs_M(G,rev):
    n,mul,inv=G; tot=0
    def conj(a,b): return mul[mul[a][b]][inv[a]]
    for x in range(n):
        for u in range(n):
            u1=conj(x,u)
            if mul[u][u1]!=mul[u1][u]: continue
            u2=conj(x,u1)
            if mul[u][u2]!=mul[u2][u]: continue
            u3=conj(x,u2)
            if not rev: rhs=mul[mul[u2][mul[inv[u1]][inv[u1]]]][u]          # t^3 = t^2 - 2t + 1
            else:       rhs=mul[mul[mul[u2][u2]][inv[u1]]][u]               # t^3 = 2t^2 - t + 1
            tot+= (u3==rhs)
    return tot
code=sys.argv[3]
t=json.load(open(code+'/data/complexity/targeted_ch12_13_census.json'))
c11={r['label']:r for r in json.load(open(code+'/data/expanded_ch11.json'))['entries']}
cands=[r for r in t['ch12_m2_s2_knots']['ambiguous_shadow_collisions'] if r['e1_generators']==['x**3 - x**2 + 2*x - 1']]
new=[r for r in cands if r['pi1_presentation'].count('=')>=3][0]
d_new=MarkedDiagram.from_epd(new['epd']); d_rib=MarkedDiagram.from_epd(c11['11_2']['epd'])
print('ch=12 diagram:',new['epd'])
print(f"{'group':14s} {'order':>6s} {'ch12 4-rel':>11s} {'M (f)':>8s} {'M (f rev)':>10s} {'11_2':>8s}")
for name,G in GROUPS.items():
    T=table(G); assert name!='AGammaL(1,8)' or T[0]==168, T[0]
    print(f"{name:14s} {T[0]:6d} {count_homs(d_new,T):11d} {homs_M(T,False):8d} {homs_M(T,True):10d} {count_homs(d_rib,T):8d}",flush=True)
