"""Independent check: do the two ch=12,m=2 'AMBIGUOUS ~ 0_1' diagrams have a non-abelian quotient F_9 x| Z_4 (order 36)?"""
import json, re, itertools
d = json.load(open('../sk/code/data/complexity/targeted_ch12_13_census.json'))
def parse(pres):
    body = pres[pres.index('|')+1:pres.rindex('>')]
    rels = []
    for r in body.split(','):
        r = r.strip().replace('= 1','').strip()
        if not r: continue
        w = []
        for tok in r.split():
            m = re.fullmatch(r'x(\d+)(?:\^(-?\d+))?', tok)
            g = int(m.group(1)); e = int(m.group(2) or 1)
            w += [g if e > 0 else -g] * abs(e)
        rels.append(w)
    return rels
# F_9 = F_3[i]; elements (p,q) = p + q i.  Affine maps z -> u z + a with u in <i>.
def mul(a, b): return ((a[0]*b[0] - a[1]*b[1]) % 3, (a[0]*b[1] + a[1]*b[0]) % 3)
def add(a, b): return ((a[0]+b[0]) % 3, (a[1]+b[1]) % 3)
I = (0, 1); ONE = (1, 0); ZERO = (0, 0)
def comp(f, g):  # f o g ; maps (u,a): z->u z + a
    return (mul(f[0], g[0]), add(mul(f[0], g[1]), f[1]))
def invm(f):
    u = f[0]; ui = u
    for _ in range(2): ui = mul(ui, u)  # u^3 = u^-1 since u^4=1
    return (ui, mul(ui, ((-f[1][0]) % 3, (-f[1][1]) % 3)))
F9 = [(p, q) for p in range(3) for q in range(3)]
def count(rels, ngen=2):
    n = 0; nonab = 0
    for trans in itertools.product(F9, repeat=ngen):
        imgs = [(I, a) for a in trans]
        ok = True
        for w in rels:
            acc = (ONE, ZERO)
            for lit in w:
                m = imgs[abs(lit)-1]
                acc = comp(acc, m if lit > 0 else invm(m))
            if acc != (ONE, ZERO): ok = False; break
        if ok:
            n += 1
            if len(set(trans)) > 1: nonab += 1
    return n, nonab
for a in d['ch12_m2_s2_knots']['ambiguous_shadow_collisions']:
    rels = parse(a['pi1_presentation'])
    print(a['e1_generators'], a['collision_matches'], "homs to F9:Z4 with meridians -> order-4 affine maps (total, non-abelian image):", count(rels) if rels else 'free')
for e in d['ch12_m2_s2_knots']['new_prime_classes']:
    print(e['class_id'], e['e1_generators'], count(parse(e['pi1_presentation'])))
print("--- all homs to G = F_9 x| Z_4 (order 36): (total, nonabelian-image)")
U = [ONE, I, mul(I, I), mul(I, mul(I, I))]
G = [(u, a) for u in U for a in F9]
def count_all(rels, ngen=2):
    n = 0; nonab = 0
    for imgs in itertools.product(G, repeat=ngen):
        ok = True
        for w in rels:
            acc = (ONE, ZERO)
            for lit in w:
                m = imgs[abs(lit)-1]
                acc = comp(acc, m if lit > 0 else invm(m))
            if acc != (ONE, ZERO): ok = False; break
        if ok:
            n += 1
            if comp(imgs[0], imgs[1]) != comp(imgs[1], imgs[0]): nonab += 1
    return n, nonab
for a in d['ch12_m2_s2_knots']['ambiguous_shadow_collisions']:
    rels = parse(a['pi1_presentation'])
    if rels: print(a['e1_generators'], a['collision_matches'], count_all(rels), a['epd'])
print("unknot Z: (36, 0)")
y = json.load(open('../sk/code/data/yoshikawa_ch10.json'))
for e in y:
    if e['num_components']==1 and e['component_types']==['S^2'] and e['pi1_generators']==2:
        print(e['label'], e['e1_generators'], count_all(parse(e['pi1_presentation'])))
c11 = json.load(open('../sk/code/data/expanded_ch11.json'))
for e in c11['entries']:
    if e['num_components']==1 and e['component_types']==['S^2'] and e['pi1_generators']==2:
        print(e['label'], e['e1_generators'], count_all(parse(e['pi1_presentation'])))
