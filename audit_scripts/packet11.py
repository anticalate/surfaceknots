import sys, itertools, time
from sympy.combinatorics import Permutation, PermutationGroup
from sympy.combinatorics.named_groups import SymmetricGroup, AlternatingGroup
def table(G):
    els=list(G.elements); idx={e:i for i,e in enumerate(els)}; n=len(els)
    mul=[[idx[els[i]*els[j]] for j in range(n)] for i in range(n)]
    inv=[idx[els[i]**-1] for i in range(n)]
    return n,mul,inv
def psl27():
    # PSL(2,7) acting on the projective line over F_7 (8 points: 0..6, inf=7)
    def mob(a,b,c,d):
        img=[]
        for x in list(range(7))+[None]:
            if x is None: num,den=a,c
            else: num,den=(a*x+b)%7,(c*x+d)%7
            img.append(7 if den==0 else (num*pow(den,-1,7))%7)
        return Permutation(img)
    return PermutationGroup([mob(1,1,0,1),mob(0,-1,1,0)])
def agl(q):   # x -> u x + v over F_q, q prime
    g=next(g for g in range(2,q) if len({pow(g,k,q) for k in range(q-1)})==q-1)
    return PermutationGroup([Permutation([(x+1)%q for x in range(q)]), Permutation([(g*x)%q for x in range(q)])])
GROUPS={'S5':SymmetricGroup(5),'A6':AlternatingGroup(6),'PSL(2,7)':psl27(),'AGL(1,5)':agl(5),'AGL(1,13)':agl(13),'AGL(1,17)':agl(17),'AGL(1,29)':agl(29),'AGL(1,37)':agl(37),'AGL(1,41)':agl(41)}
WORDS={'11_1.b R(b a b^-1 a^-1)':'baBA','11_1.c R(b a b^-1 a)':'baBa','11_2 R(b a b a^-1)':'babA','10_1 R(b^-1 a b a^-1)':'BabA'}
for gname,G in GROUPS.items():
    t=time.time(); n,mul,inv=table(G)
    out={}
    for wn,W in WORDS.items():
        tot=0
        for a in range(n):
            ia=inv[a]
            for b in range(n):
                img={'a':a,'A':ia,'b':b,'B':inv[b]}
                w=img[W[0]]
                for ch in W[1:]: w=mul[w][img[ch]]
                if mul[a][w]==mul[w][b]: tot+=1
        out[wn]=tot
    print(f"{gname:10s} order {n:5d}  "+"  ".join(f"{k.split()[0]}:{v}" for k,v in out.items()), f"({time.time()-t:.0f}s)",flush=True)
