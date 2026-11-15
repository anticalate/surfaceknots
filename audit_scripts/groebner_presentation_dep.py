# Does the E_1 key computed by Carbon's alexander.py depend on the presentation?  Apply a Nielsen change of
# generators (old g_n := y_1^-1 y_2, old g_o := y_1) to the simplified presentation, carrying the component
# meridians along, and recompute.  A correct invariant must not change.  (Idea from the third-party audit;
# this is my own implementation, run on my snapshot of Carbon's tree.)
import json, signal, sys, dataclasses
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.invariants.alexander import compute_alexander_invariants
from surfaceknots.invariants.fundamental_group import GroupPresentation, wirtinger_presentation
from surfaceknots.yoshikawa_data import YOSHIKAWA_TABLE
class TO(Exception): pass
def h(*a): raise TO()
signal.signal(signal.SIGALRM,h)
c11=json.load(open('data/expanded_ch11.json'))
pool=[(e.label,e.canonical_epd) for e in YOSHIKAWA_TABLE]+[(r['label'],r['epd']) for r in c11['entries']]
def inv(w): return tuple(-x for x in reversed(w))
def summary(al):
    out={}
    for k in ('e1_generators','evaluation_fingerprint','groebner_key','e1_groebner_key','canonical_groebner_key'):
        if hasattr(al,k): out[k]=getattr(al,k)
    return out
for lbl,epd in pool:
    d=MarkedDiagram.from_epd(epd); tp=d.surface_topology()
    if not (tp.num_components>=2 and any(c.orientable for c in tp.components) and any(not c.orientable for c in tp.components)): continue
    pa=wirtinger_presentation(d,simplify=True); n=pa.num_generators
    osc=next(i for i,c in enumerate(tp.components) if c.orientable); nsc=next(i for i,c in enumerate(tp.components) if not c.orientable)
    cnt=lambda w,g: sum(1 for x in w if abs(x)==g)%2
    odd_o=[g for g in range(1,n+1) if cnt(pa.component_meridians[osc],g)==1 and cnt(pa.component_meridians[nsc],g)==0]
    odd_n=[g for g in range(1,n+1) if cnt(pa.component_meridians[nsc],g)==1 and cnt(pa.component_meridians[osc],g)==0]
    if not (odd_o and odd_n): print(f"{lbl:14s} no suitable generator pair"); continue
    go,gn=odd_o[0],odd_n[0]; mp={}; nxt=3
    for j in range(1,n+1):
        if j==go: mp[j]=(1,)
        elif j==gn: mp[j]=(-1,2)
        else: mp[j]=(nxt,); nxt+=1
    def sub(w):
        out=[]
        for lit in w:
            wd=mp[abs(lit)]; out.extend(wd if lit>0 else inv(wd))
        return tuple(out)
    pb=GroupPresentation(num_generators=n,relators=tuple(sub(r) for r in pa.relators),component_meridians=tuple(sub(m) for m in pa.component_meridians))
    try:
        signal.alarm(20); a=compute_alexander_invariants(pres=pa,topo=tp); b=compute_alexander_invariants(pres=pb,topo=tp); signal.alarm(0)
    except TO:
        print(f"{lbl:14s} timeout"); continue
    except Exception as ex:
        signal.alarm(0); print(f"{lbl:14s} error {repr(ex)[:80]}"); continue
    sa,sb=summary(a),summary(b)
    diffs=[k for k in sa if sa[k]!=sb.get(k)]
    print(f"{lbl:14s} fields that changed under the change of generators: {diffs or 'none'}")
    for k in diffs: print(f"      {k}: {str(sa[k])[:110]}  ->  {str(sb[k])[:110]}")
