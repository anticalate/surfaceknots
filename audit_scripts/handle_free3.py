"""Three unknotted spheres + one tube: read the cord word in F(a,b,c) from a diagram.
F_3 is realised inside Sanov's free group <A,B> < SL(2,Z) as the index-2 subgroup with free basis
a = A^2, b = B, c = A B A^-1 (Schreier generators for the transversal {1, A})."""
import sys, json, itertools, importlib.util
spec=importlib.util.spec_from_file_location('hf', sys.argv[3]+'/handle_free.py'); hf=importlib.util.module_from_spec(spec); spec.loader.exec_module(hf)
sys.path.insert(0, sys.argv[3])
from indep_pi1_lib import relations
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.moves import list_level_yoshikawa_neighbors, canonical_marked_key
mm,mi,I2=hf.mm,hf.mi,hf.I2
SA,SB=hf.GEN['a'],hf.GEN['b']
G3={'a':mm(SA,SA),'b':SB,'c':mm(mm(SA,SB),mi(SA))}
for k in list(G3): G3[k.upper()]=mi(G3[k])
INV3={k:k.swapcase() for k in G3}
def red(w):
    out=[]
    for ch in w:
        if out and out[-1]==INV3[ch]: out.pop()
        else: out.append(ch)
    return ''.join(out)
def winv(w): return ''.join(INV3[c] for c in reversed(w))
def wmat(w):
    M=I2
    for ch in w: M=mm(M,G3[ch])
    return M
def decode3(M):
    w=hf.decode(M)                     # word in Sanov letters: a,A = A^{+-1}; b,B = B^{+-1}
    if w is None: return None
    out=[]; p=0
    for ch in w:
        if ch=='b': out.append('b' if p==0 else 'c')
        elif ch=='B': out.append('B' if p==0 else 'C')
        elif ch=='a':
            if p==1: out.append('a')
            p^=1
        else:
            if p==0: out.append('A')
            p^=1
    if p!=0: return None
    return red(''.join(out))
def conj_form(w):
    i=0; j=len(w)-1
    while i<j and w[i]==INV3[w[j]]: i+=1; j-=1
    return (w[:i],w[i]) if i==j else None
def strip(w,l,r):
    w=red(w)
    while w and w[0].lower()==l: w=w[1:]
    while w and w[-1].lower()==r: w=w[:-1]
    return w
def free_rep3(cand):
    k,eid,rels=relations(cand)
    comp,nc=cand._dart_surface_components()
    if nc!=3: return None
    ecomp={e:comp[x] for x,(e,s) in eid.items()}
    def val(asg,x):
        e,s=eid[x]; v=asg[e]; return None if v is None else (v if s==1 else mi(v))
    def setv(asg,x,g):
        e,s=eid[x]; g2=g if s==1 else mi(g)
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
                        if x!=mi(y): return False
                    elif x is not None:
                        if not setv(asg,r[2],mi(x)): return False
                        ch=True
                    elif y is not None:
                        if not setv(asg,r[1],mi(y)): return False
                        ch=True
                else:
                    _,out,inn,over,eps=r
                    O=val(asg,over)
                    if O is None: continue
                    if eps<0: O=mi(O)
                    x,y=val(asg,out),val(asg,inn)
                    if y is not None:
                        t=mi(mm(mm(O,y),mi(O)))
                        if x is not None:
                            if x!=t: return False
                        elif not setv(asg,out,t): return False
                        else: ch=True
                    elif x is not None:
                        t=mm(mm(mi(O),mi(x)),O)
                        if not setv(asg,inn,t): return False
                        ch=True
        return True
    conjs=['']+[''.join(p) for L in (1,2) for p in itertools.product('aAbBcC',repeat=L) if red(''.join(p))==''.join(p)]
    def cands(e):
        g='abc'[ecomp[e]]
        for c in conjs:
            for gg in (g,g.upper()): yield wmat(red(c+gg+winv(c)))
    def rec(asg,depth):
        if all(v is not None for v in asg): return asg
        if depth==0: return None
        e=asg.index(None)
        for M in cands(e):
            a2=list(asg); a2[e]=M
            if propagate(a2):
                r=rec(a2,depth-1)
                if r is not None: return r
        return None
    seeds=[next(e for e in range(k) if ecomp[e]==i) for i in range(3)]
    for s1 in ('b','B'):
        for s2 in ('c','C'):
            asg=[None]*k; asg[seeds[0]]=G3['a']; asg[seeds[1]]=G3[s1]; asg[seeds[2]]=G3[s2]
            if not propagate(asg): continue
            r=rec(asg,2)
            if r is not None: return eid,r,comp
    return None
def analyse3(d):
    out=[]
    for u,w,cand,dA,dB in hf.bigon_bases(d):
        cst=cand.surface_topology()
        if cst.num_components!=3 or any(c.homeomorphism_type!='S^2' for c in cst.components): continue
        rep=free_rep3(cand)
        if rep is None: out.append('no free rep'); continue
        eid,asg,comp=rep
        def val(x):
            e,s=eid[x]; v=asg[e]; return v if s==1 else mi(v)
        wa,wb=decode3(val(dA)),decode3(val(dB))
        if wa is None or wb is None: out.append('decode failed'); continue
        ca,cb=conj_form(wa),conj_form(wb)
        if ca is None or cb is None: out.append('not conjugates'); continue
        (c,gx),(dd,gy)=ca,cb
        x,y=gx.lower(),gy.lower(); s=(1 if gx.islower() else -1)*(1 if gy.islower() else -1)
        W=strip(winv(c)+dd,x,y)
        out.append(f"tube from sphere {x} to sphere {y}, sign {s:+d}, cord {hf.pretty(W) if W else '1'}   (relation {x} = w {y}^{s:+d} w^-1)")
    return out
code=sys.argv[2]
from surfaceknots.yoshikawa_data import YOSHIKAWA_BY_LABEL
c11={r['label']:r['epd'] for r in json.load(open(code+'/data/expanded_ch11.json'))['entries']}
for lab,epd in (('10_1^0,0,1',YOSHIKAWA_BY_LABEL['10_1^0,0,1'].canonical_epd),('11_1^0,0,-2',c11['11_1^0,0,-2'])):
    d0=MarkedDiagram.from_epd(epd); seen={canonical_marked_key(d0):d0}; q=[d0]; i=0; res=[]
    while i<len(q) and i<300 and not [r for r in res if r.startswith('tube')]:
        cur=q[i]; i+=1; res+=analyse3(cur)
        for nb in list_level_yoshikawa_neighbors(cur):
            kk=canonical_marked_key(nb)
            if kk not in seen: seen[kk]=nb; q.append(nb)
    print(lab,'| diagrams examined',i,'|',sorted(set(res)) or 'no bigon over three spheres found',flush=True)
