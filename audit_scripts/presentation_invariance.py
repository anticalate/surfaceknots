"""Presentation-independence test of the whole invariant key (original code).
Every call of wirtinger_presentation(..., simplify=True) is wrapped so that the returned presentation is
changed by a random sequence of Nielsen moves (g_i := g_i g_j^{+-1}, g_i := g_i^-1, swaps), with relators and
component meridians rewritten accordingly.  The group and its meridians are unchanged, so every layer of the
key must be unchanged.  Run from code/ with PYTHONPATH=src."""
import sys, json, random, dataclasses, importlib, pkgutil
import surfaceknots
from surfaceknots.invariants import fundamental_group as fg
from surfaceknots import MarkedDiagram, YOSHIKAWA_TABLE
orig = fg.wirtinger_presentation
RNG = random.Random(0); ENABLED = [False]
def inv(w): return tuple(-x for x in reversed(w))
def red(w):
    out=[]
    for x in w:
        if out and out[-1]==-x: out.pop()
        else: out.append(x)
    return tuple(out)
def nielsen(p, k):
    n=p.num_generators
    if n<2: return p
    rel=list(p.relators); mer=list(p.component_meridians)
    for _ in range(k):
        t=RNG.choice(('mul','inv','swap'))
        i,j=RNG.sample(range(1,n+1),2)
        if t=='mul':
            e=RNG.choice((1,-1)); mp=lambda g: (i, -e*j) if g==i else (g,)       # old g_i = y_i y_j^{-e}
        elif t=='inv': mp=lambda g: (-i,) if g==i else (g,)
        else: mp=lambda g: (j,) if g==i else ((i,) if g==j else (g,))
        def sub(w):
            out=[]
            for lit in w:
                wd=mp(abs(lit)); out.extend(wd if lit>0 else inv(wd))
            return red(out)
        rel=[sub(r) for r in rel]; mer=[sub(m) for m in mer]
    return dataclasses.replace(p, relators=tuple(rel), component_meridians=tuple(mer))
def wrapped(d, *a, **kw):
    p=orig(d,*a,**kw)
    simp = kw.get('simplify', True)
    if ENABLED[0] and simp: return nielsen(p, 3)
    return p
for m in list(sys.modules.values()):
    if m and getattr(m,'__name__','').startswith('surfaceknots') and getattr(m,'wirtinger_presentation',None) is orig:
        setattr(m,'wirtinger_presentation',wrapped)
from surfaceknots.yoshikawa_data import compute_full_invariant_record, invariant_filtration_layers
c11=json.load(open('data/expanded_ch11.json'))
pool=[(e.label,e.canonical_epd) for e in YOSHIKAWA_TABLE if e.label!='0_1']+[(r['label'],r['epd']) for r in c11['entries']]
def layers(d):
    rec=compute_full_invariant_record(d,"x",include_sublinks=False,compute_so4_symmetry_bfs=False)
    return {k:repr(v) for k,v in invariant_filtration_layers(rec)}
tot=0; bad={}
for lbl,epd in pool:
    d=MarkedDiagram.from_epd(epd)
    ENABLED[0]=False; base=layers(d)
    for trial in range(int(sys.argv[1])):
        ENABLED[0]=True
        try: L=layers(d)
        except Exception as ex:
            bad.setdefault('EXCEPTION',[]).append((lbl,repr(ex)[:80])); continue
        finally: ENABLED[0]=False
        tot+=1
        for k in base:
            if L[k]!=base[k]: bad.setdefault(k,[]).append(lbl)
    print(lbl,'done',flush=True)
print('variants',tot)
for k,v in bad.items():
    from collections import Counter
    print('LAYER',k,'changed:',dict(Counter(v)) if k!='EXCEPTION' else v[:6])
if not bad: print('no layer changed')
