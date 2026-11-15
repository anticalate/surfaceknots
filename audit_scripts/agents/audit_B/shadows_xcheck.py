import time, sys, itertools
from surfaceknots.fast_enum import generate_prime_shadows_c_list
from surfaceknots.quadrangulations import generate_quadrangulations, quadrangulation_to_shadow

def canon_shadow(n, opp):
    best = None
    refl = [0,3,2,1]
    opp_r = [0]*(4*n)
    for d in range(4*n):
        v,s = divmod(d,4); od = opp[d]; ov,os = divmod(od,4)
        opp_r[4*v+refl[s]] = 4*ov+refl[os]
    for o in (opp, opp_r):
        for d0 in range(4*n):
            v0,s0 = divmod(d0,4)
            num = {v0:0}; shift={v0:s0}; q=[v0]; code=[]
            i=0
            while i < len(q):
                v=q[i]; i+=1
                for k in range(4):
                    od = o[4*v+((shift[v]+k)&3)]
                    ov,os = divmod(od,4)
                    if ov not in num:
                        num[ov]=len(q); shift[ov]=os; q.append(ov)
                    code.append((num[ov], (os-shift[ov])&3))
            assert len(q)==n, "disconnected"
            t=tuple(code)
            if best is None or t<best: best=t
    return best

def check_props(n, opp):
    # involution, no loops, planar, no 2-edge cut (any sides)
    for d in range(4*n):
        assert opp[opp[d]]==d and opp[d]!=d
        assert opp[d]//4 != d//4, "loop"
    # faces
    seen=[0]*(4*n); F=0
    for d in range(4*n):
        if seen[d]: continue
        F+=1; cur=d
        while not seen[cur]:
            seen[cur]=1; od=opp[cur]; cur=(od//4)*4+((od%4)-1)%4
    assert n - 2*n + F == 2, "nonplanar"
    edges=[(d,opp[d]) for d in range(4*n) if d<opp[d]]
    for (a,b),(c,e) in itertools.combinations(edges,2):
        banned={a,b,c,e}
        vis={0}; st=[0]
        while st:
            u=st.pop()
            for s in range(4):
                dd=4*u+s
                if dd in banned: continue
                w=opp[dd]//4
                if w not in vis: vis.add(w); st.append(w)
        assert len(vis)==n, "2-edge cut"

maxn=int(sys.argv[1])
for n in range(2,maxn+1):
    t=time.time()
    cs = generate_prime_shadows_c_list(n)
    cset=set()
    for sh in cs:
        check_props(n, sh.opp)
        cset.add(canon_shadow(n, sh.opp))
    assert len(cset)==len(cs), "C dup"
    t1=time.time()
    ps = [quadrangulation_to_shadow(q) for q in generate_quadrangulations(n)]
    pset=set(canon_shadow(n, sh.opp) for sh in ps)
    print(n, "C", len(cs), "Py", len(ps), "py distinct", len(pset), "sets equal", cset==pset, "props ok", "tC %.1f tPy %.1f"%(t1-t, time.time()-t1), flush=True)
