"""Name a surface-link as (trivial 2-component sphere link) + (one 1-handle with cord word w).

A marked bigon (u,w) in a ch-diagram d exhibits d = base + 1-handle.  If the base is the unlink of two
spheres, pi_1(base complement) = F(a,b) and each Wirtinger generator is a conjugate c g^{+-1} c^{-1} of a
basis element.  Reading the meridians at the two feet of the handle, m_A = c x^e c^-1 and m_B = d y^f d^-1,
the handle relation is m_A = m_B^{+-1} and the cord class is  c^-1 d  in  <x>\F/<y>.
  x == y : the handle is attached to one sphere  -> torus or Klein bottle linked with the other sphere
  x != y : the handle fuses the two spheres       -> 1-fusion ribbon 2-knot R(w) (or a projective plane...)
F(a,b) is realised faithfully by Sanov's matrices A=[[1,2],[0,1]], B=[[1,0],[2,1]] in SL(2,Z).
"""
import json, sys, itertools
sys.path.insert(0, sys.argv[1]); sys.path.insert(0, sys.argv[3])
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.moves import _is_planar_4reg_diagram
from surfaceknots.unlink import certify_classical_unlink
from indep_pi1_lib import relations

def mm(A,B): return ((A[0][0]*B[0][0]+A[0][1]*B[1][0], A[0][0]*B[0][1]+A[0][1]*B[1][1]),(A[1][0]*B[0][0]+A[1][1]*B[1][0], A[1][0]*B[0][1]+A[1][1]*B[1][1]))
def mi(A): return ((A[1][1],-A[0][1]),(-A[1][0],A[0][0]))
I2=((1,0),(0,1))
GEN={'a':((1,2),(0,1)),'A':((1,-2),(0,1)),'b':((1,0),(2,1)),'B':((1,0),(-2,1))}
INV={'a':'A','A':'a','b':'B','B':'b'}
def word_mat(w):
    M=I2
    for ch in w: M=mm(M,GEN[ch])
    return M
def red(w):
    out=[]
    for ch in w:
        if out and out[-1]==INV[ch]: out.pop()
        else: out.append(ch)
    return ''.join(out)
def winv(w): return ''.join(INV[c] for c in reversed(w))
def decode(M):
    """Reduced word of M in <A,B> (ping-pong: the correct first letter is the one whose inverse shrinks M)."""
    w=[]; n=lambda X: sum(x*x for r in X for x in r)
    while M!=I2:
        best=None
        for ch in 'aAbB':
            X=mm(GEN[INV[ch]],M)
            if n(X)<n(M) and (best is None or n(X)<best[1]): best=(ch,n(X),X)
        if best is None: return None
        w.append(best[0]); M=best[2]
        if len(w)>400: return None
    return ''.join(w)
def conj_form(w):
    """w reduced, conjugate of a single letter:  w = c g c^-1  -> (c, g)."""
    i=0; j=len(w)-1
    while i<j and w[i]==INV[w[j]]: i+=1; j-=1
    if i!=j: return None
    return w[:i], w[i]

def bigon_bases(d, *, require_compatible_parity: bool = True, return_parity_flag: bool = False):
    """Yield marked-bigon 1-handle surgeries `(u, w, cand, dart_a, dart_b)` (plus `is_compatible_tube`
    when `return_parity_flag=True`).

    Marker-parity compatibility `is_compatible_tube = ((d.vtype[u] == d.vtype[w]) == ((su & 1) == (sw_a & 1)))`
    distinguishes genuine standard 1-handle tubes (where `L_+` and `L_-` of `cand` match those of `d`)
    from simultaneous saddle-pair birth/death tubes (`is_compatible_tube == False`, which admit an
    `Omega_7` slide / `Omega_6` stabilization reduction).
    """
    n=d.num_vertices; st=d.surface_topology()
    for u in range(n):
        if d.vtype[u] not in ("Y","Z"): continue
        for su in range(4):
            d_a=d.opp[4*u+su]; d_b=d.opp[4*u+((su+1)&3)]; w=d_a>>2
            if w<=u or (d_b>>2)!=w or d.vtype[w] not in ("Y","Z"): continue
            sw_a=d_a&3; sw_b=d_b&3
            if ((sw_b+1)&3)!=sw_a: continue
            is_compatible_tube = ((d.vtype[u]==d.vtype[w]) == ((su&1)==(sw_a&1)))
            if require_compatible_parity and not is_compatible_tube: continue
            p0,p1=4*u+((su+2)&3),4*u+((su+3)&3); p2,p3=4*w+((sw_b+2)&3),4*w+((sw_b+3)&3)
            if any((d.opp[p]>>2) in (u,w) for p in (p0,p1,p2,p3)): continue
            remap=[v-(1 if v>u else 0)-(1 if v>w else 0) for v in range(n)]; m=n-2
            nvt=tuple(d.vtype[v] for v in range(n) if v not in (u,w))
            partner={p0:p1,p1:p0,p2:p3,p3:p2}
            nopp=[0]*(4*m)
            for v in range(n):
                if v in (u,w): continue
                for s in range(4):
                    od=d.opp[4*v+s]
                    if (od>>2) in (u,w): od=d.opp[partner[od]]
                    nopp[4*remap[v]+s]=(remap[od>>2]<<2)|(od&3)
            cand=MarkedDiagram(num_vertices=m,opp=tuple(nopp),vtype=nvt,free_loops=0)
            if not _is_planar_4reg_diagram(cand): continue
            if not (certify_classical_unlink(cand.resolve(1)).is_unlink and certify_classical_unlink(cand.resolve(-1)).is_unlink): continue
            cst=cand.surface_topology()
            if cst.euler_characteristic!=st.euler_characteristic+2: continue
            far=lambda x: (remap[d.opp[x]>>2]<<2)|(d.opp[x]&3)
            if return_parity_flag:
                yield u,w,cand,far(p0),far(p2),is_compatible_tube
            else:
                yield u,w,cand,far(p0),far(p2)

def free_rep(cand):
    """Assign reduced words (via Sanov matrices) to all edge generators of a 2-sphere unlink diagram."""
    k,eid,rels=relations(cand)
    comp,nc=cand._dart_surface_components()
    if nc!=2: return None
    ecomp={}
    for x,(e,s) in eid.items(): ecomp[e]=comp[x]
    def val(asg,x):
        e,s=eid[x]; v=asg[e]
        return None if v is None else (v if s==1 else mi(v))
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
                        else:
                            if not setv(asg,out,t): return False
                            ch=True
                    elif x is not None:
                        t=mm(mm(mi(O),mi(x)),O)
                        if not setv(asg,inn,t): return False
                        ch=True
        return True
    conjs=['']+[''.join(p) for L in (1,2) for p in itertools.product('aAbB',repeat=L) if red(''.join(p))==''.join(p)]
    def cands(e):
        g='a' if ecomp[e]==0 else 'b'
        for c in conjs:
            for gg in (g,INV[g]):
                yield word_mat(red(c+gg+winv(c)))
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
    # seed: one edge on each component gets the basis letter itself
    e0=next(e for e in range(k) if ecomp[e]==0); e1=next(e for e in range(k) if ecomp[e]==1)
    for s1 in ('b','B'):
        asg=[None]*k; asg[e0]=GEN['a']; asg[e1]=GEN[s1]
        if not propagate(asg): continue
        r=rec(asg,3)
        if r is not None:
            return eid, r
    return None

def strip(w,left,right):
    w=red(w)
    while w and w[0].lower()==left: w=w[1:]
    while w and w[-1].lower()==right: w=w[:-1]
    return w
def flip(t,g): return ''.join(c.swapcase() if c.lower()==g else c for c in t)
def variants(w,x,y,s):
    """All equivalent descriptions (x, y, word, sign) of the handle  x = w y^s w^-1  on the trivial 2-link:
    reverse the handle, flip the orientation of either sphere, swap the names of the spheres."""
    out=set()
    for rev in (False,True):
        ww,xx,yy=(winv(w),y,x) if rev else (w,x,y)
        for fa in (False,True):
            for fb in (False,True):
                for sw in (False,True):
                    t=ww; ss=s
                    if fa:
                        t=flip(t,'a'); ss*= (-1 if (xx=='a')!=(yy=='a') else 1)
                    if fb:
                        t=flip(t,'b'); ss*= (-1 if (xx=='b')!=(yy=='b') else 1)
                    X,Y=xx,yy
                    if sw:
                        t=t.translate(str.maketrans('aAbB','bBaA')); X={'a':'b','b':'a'}[X]; Y={'a':'b','b':'a'}[Y]
                    out.add((X,Y,strip(t,X,Y),ss))
    return out
def canon(w,x,y,s):
    vs=variants(w,x,y,s)
    if x!=y: vs={v for v in vs if v[3]==1 and v[0]=='a'} or vs      # fusion: normal form  a = w b w^-1
    else:    vs={v for v in vs if v[0]=='a'} or vs                 # self-handle on sphere a
    return min(vs, key=lambda t:(len(t[2]), sum(c.isupper() for c in t[2]), t[2].swapcase()))
def pretty(w):
    if not w: return '1'
    out=[]
    for c,grp in itertools.groupby(w):
        n=len(list(grp)); base=c.lower(); e=n if c.islower() else -n
        out.append(base if e==1 else f"{base}^{{{e}}}")
    return ' '.join(out)

def analyse(d):
    res=[]
    for u,w,cand,dA,dB in bigon_bases(d):
        cst=cand.surface_topology()
        if cst.num_components!=2 or any(c.homeomorphism_type!='S^2' for c in cst.components): continue
        rep=free_rep(cand)
        if rep is None: res.append((u,w,None,'no free rep found')); continue
        eid,asg=rep
        def val(x):
            e,s=eid[x]; v=asg[e]; return v if s==1 else mi(v)
        wa,wb=decode(val(dA)),decode(val(dB))
        if wa is None or wb is None: res.append((u,w,None,'decode failed')); continue
        ca,cb=conj_form(wa),conj_form(wb)
        if ca is None or cb is None: res.append((u,w,None,'not conjugates of letters')); continue
        (c,gx),(dd,gy)=ca,cb
        x,y=gx.lower(),gy.lower()
        s=(1 if gx.islower() else -1)*(1 if gy.islower() else -1)      # handle relation: m_A = m_B
        X,Y,W,S=canon(strip(winv(c)+dd,x,y),x,y,s)
        if X==Y: name=('T' if S==1 else 'K')+f"({pretty(W)})"; rel=f"a = w a^{{{S}}} w^-1"
        else:    name=f"R({pretty(W)})"; rel="a = w b w^-1" if S==1 else "a = w b^-1 w^-1"
        res.append((u,w,(X,Y,W,S),name))
    return res

def name_surface(d, max_states=400):
    """Names from the diagram itself, or from level-equivalent diagrams if it shows no handle."""
    from surfaceknots.moves import list_level_yoshikawa_neighbors, canonical_marked_key
    seen={canonical_marked_key(d):d}; q=[d]; i=0; names={}
    while i<len(q) and i<max_states:
        cur=q[i]; i+=1
        for u,w,cl,nm in analyse(cur):
            if cl is not None: names[nm]=names.get(nm,0)+1
        if names and i>=1 and cur is d: break
        if names: break
        for nb in list_level_yoshikawa_neighbors(cur):
            k=canonical_marked_key(nb)
            if k not in seen: seen[k]=nb; q.append(nb)
    return names, i

if __name__=='__main__':
    code=sys.argv[2]
    from surfaceknots.yoshikawa_data import YOSHIKAWA_TABLE
    pool=[('TABLE '+e.label,e.canonical_epd) for e in YOSHIKAWA_TABLE]
    c11=json.load(open(code+'/data/expanded_ch11.json'))
    for r in c11['entries']:
        for i,epd in enumerate(r['move_component_epds']): pool.append((f"CH11 {r['label']}"+(f".{'abc'[i]}" if len(r['move_component_epds'])>1 else ''),epd))
    for i,r in enumerate(c11['residual_collisions']): pool.append((f"resid{i} [{r['collision_matches'][0]}]",r['epd']))
    for lab,epd in pool:
        d=MarkedDiagram.from_epd(epd)
        st=d.surface_topology()
        typ='+'.join(sorted(c.homeomorphism_type for c in st.components))
        if d.num_vertices<4 or typ not in ('S^2','S^2+T^2','K^2+S^2','RP^2'): continue
        names,explored=name_surface(d, 60 if lab.startswith('resid') else 400)
        print(f"{lab:44s} {typ:9s} {names if names else '-'}   (diagrams examined: {explored})",flush=True)
