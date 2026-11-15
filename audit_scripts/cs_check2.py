import json,re,sys,time
sys.path.insert(0,sys.argv[3])
from cs_parse import parse, items
from sympy.combinatorics.fp_groups import FpGroup
idx=int(sys.argv[1]); k=int(sys.argv[2])
name,pres,epd=items()[idx]
F,xs,rels=parse(pres)
G=FpGroup(F,rels+[xs[0]**k]); o=G.order()
extra=''
if o!=float('inf') and 1<o<3000:
    P,T=G._to_perm_group(); extra=f" derived series {[H.order() for H in P.derived_series()]}"
print(f"{name[:40]:40s} |G/<<x1^{k}>>| = {o}{extra}",flush=True)
