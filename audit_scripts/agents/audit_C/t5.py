import sys, json
from mypi1 import *
from surfaceknots import YOSHIKAWA_BY_LABEL, MarkedDiagram
from surfaceknots.invariants.fundamental_group import wirtinger_presentation
from sympy.combinatorics.perm_groups import PermutationGroup
from sympy.combinatorics import Permutation
from sympy.combinatorics.named_groups import SymmetricGroup
ds = {lbl: MarkedDiagram.from_epd(e.canonical_epd) for lbl,e in YOSHIKAWA_BY_LABEL.items() if lbl!="0_1"}
j = json.load(open('data/expanded_ch11.json'))
for e in j['entries']: ds[e['label']] = MarkedDiagram.from_epd(e['epd'])

def order_with(d, extra_power=None):
    ng, rels = presentation(d)
    rels = list(rels)
    if extra_power: rels.append((1,)*extra_power)   # generator x_1 is a meridian
    G, gens = fp_group(ng, rels)
    return G, gens

for lbl in ["2_1^-1","8_1^-1,-1","10_1^-1,-1","10_1^-2,-2","11_1^-2"]:
    G, gens = order_with(ds[lbl])
    o = G.order()
    print(lbl, "order", o, " repo pres:", wirtinger_presentation(ds[lbl]).to_human_string())
    P, T = G._to_perm_group()
    print("    perm group order", P.order(), "abelian", P.is_abelian, "center", P.center().order(), "derived", P.derived_subgroup().order(),
          "elt orders", sorted(__import__('collections').Counter(g.order() for g in P.elements).items()))
    if lbl == "11_1^-2":
        S4 = SymmetricGroup(4)
        from sympy.combinatorics.homomorphisms import group_isomorphism, is_isomorphic
        print("    isomorphic to S4:", is_isomorphic(P, S4))
        # meridian image
        print("    meridian (x0) order:", T(gens[0]).order(), " cycle structure in S4-action? conj class size:", len(P.conjugacy_class(T(gens[0]))))
for lbl, pw, exp in [("10_2",2,6),("10_2",4,12),("10_3",3,24),("10_3",6,48),("8_1",2,6),("8_1",3,24),("8_1",4,96),("8_1",5,600),("9_1",2,6),("9_1",3,None),("10_1",2,10),("10_1",3,None)]:
    try:
        G, gens = order_with(ds[lbl], pw)
        o = G.order()
    except Exception as ex:
        o = "fail " + repr(ex)[:60]
    print(lbl, "pi/<<mu^%d>>" % pw, "order", o, "expected", exp)
for lbl in ["8_1","9_1","10_1","10_2","10_3","6_1^0,1","8_1^1,1","10_1^1,1","10_1^1"]:
    for res in (-1,1):
        print(lbl, res, wirtinger_presentation(ds[lbl], from_resolution=res).to_human_string())
