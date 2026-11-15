import json,re
from sympy.combinatorics.free_groups import free_group
def parse(pres):
    m=re.match(r"<(.*)\|(.*)>",pres); gens=[g.strip() for g in m.group(1).split(',')]
    F,*xs=free_group(",".join(gens)); d=dict(zip(gens,xs)); rels=[]
    for r in m.group(2).split(','):
        r=r.strip()
        if not r: continue
        lhs,rhs=[s.strip() for s in r.split('=')]
        def word(s):
            w=F.identity
            for tok in s.split():
                if tok=='1': continue
                mm=re.match(r"(x\d+)(?:\^(-?\d+))?$",tok); w=w*d[mm.group(1)]**int(mm.group(2) or 1)
            return w
        rels.append(word(lhs)*word(rhs)**-1)
    return F,xs,rels
def items():
    t=json.load(open('data/complexity/targeted_ch12_13_census.json'))
    c11={r['label']:r for r in json.load(open('data/expanded_ch11.json'))['entries']}
    cands=[r for r in t['ch12_m2_s2_knots']['ambiguous_shadow_collisions'] if r['e1_generators']==['x**3 - x**2 + 2*x - 1']]
    return [('11_2 ribbon R(baba^-1)',c11['11_2']['pi1_presentation'],c11['11_2']['epd'])]+[(f"ch12 #{i} ({len(r['pi1_presentation'].split(','))-1} relators)",r['pi1_presentation'],r['epd']) for i,r in enumerate(cands)]
