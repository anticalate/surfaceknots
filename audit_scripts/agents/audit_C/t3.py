import sys, random
from mydlk import dlk_true
from variants import *
from surfaceknots import YOSHIKAWA_BY_LABEL, MarkedDiagram
from surfaceknots.invariants.linking import double_linking_mod2
for lbl in ["6_1^0,1","8_1^1,1","10_1^1,1","9_1^1,-2","10_1^0,0,1","9_1^0,1","10_1^0,1","10_2^0,1"]:
    d = MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL[lbl].canonical_epd)
    print(lbl, "repo:", double_linking_mod2(d))
    for x in dlk_true(d): print("     mine:", x)
