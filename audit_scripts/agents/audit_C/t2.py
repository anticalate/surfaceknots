import sys, json
sys.path.insert(0, '.')
from mypi1 import *
from surfaceknots import YOSHIKAWA_BY_LABEL, MarkedDiagram
from surfaceknots.invariants.fundamental_group import wirtinger_presentation, finite_quotient_profile
names = ["S3","S4","A4","A5","D8","D10","F21","Q8"]
bad = 0
ds = {lbl: MarkedDiagram.from_epd(e.canonical_epd) for lbl,e in YOSHIKAWA_BY_LABEL.items() if lbl!="0_1"}
j = json.load(open('data/expanded_ch11.json'))
for e in j['entries']: ds[e['label']] = MarkedDiagram.from_epd(e['epd'])
for lbl,d in ds.items():
    ng, rels = presentation(d)
    mine = {g: hom_count2(ng, rels, g) for g in names}
    rp = {}
    for res in (-1,1):
        wp = wirtinger_presentation(d, from_resolution=res)
        fq = finite_quotient_profile(wp)
        rp[res] = {g: fq[g] for g in names}
    flag = "" if (mine==rp[-1]==rp[1]) else "  <<<<< MISMATCH"
    if flag: bad += 1
    print(lbl, mine, flag)
    if flag: print("   repo", rp)
print("mismatches", bad)
