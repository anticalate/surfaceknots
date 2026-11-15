from surfaceknots import YOSHIKAWA_BY_LABEL, MarkedDiagram
from surfaceknots.invariants.alexander import compute_alexander_invariants
from surfaceknots.invariants.fundamental_group import wirtinger_presentation
for lbl in ["10_1^0,1","10_2^0,1","10_1^0,-2","10_2^0,-2","9_1^1,-2","8_1^1,1","10_1^1,1","6_1^0,1"]:
    d = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL[lbl].canonical_epd)
    for res in (-1,1):
        al = compute_alexander_invariants(d, from_resolution=res)
        p = wirtinger_presentation(d, from_resolution=res)
        print(lbl,res, p.to_human_string(), p.component_meridians)
        print("   E1",al.e1_generators,"E2",al.e2_generators)
        print("   fp",al.evaluation_fingerprint)
