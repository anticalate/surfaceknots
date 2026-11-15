# Audit the cyclic-tube prunes at ch=11: every diagram the rule removes must have a cyclic group.
import json, sys
from surfaceknots import MarkedDiagram
from surfaceknots.moves import has_certified_cyclic_unknotted_tube, find_planar_connected_sum_decomposition
from surfaceknots.invariants.fundamental_group import wirtinger_presentation, finite_quotient_profile
epds=json.load(open(sys.argv[1])); print('input diagrams',len(epds),flush=True)
GR=("S3","Q8","A4","S4","A5")
Zp={'S3':6,'Q8':8,'A4':12,'S4':24,'A5':60}; Z2={'S3':4,'Q8':2,'A4':4,'S4':10,'A5':16}
pruned=0; cut=0; bad=[]
for epd in epds:
    d=MarkedDiagram.from_epd(epd)
    if find_planar_connected_sum_decomposition(d) is not None: cut+=1; continue
    if has_certified_cyclic_unknotted_tube(d):
        pruned+=1
        fq=dict(finite_quotient_profile(wirtinger_presentation(d),GR))
        st=d.surface_topology(); c=st.components[0]
        if fq not in (Zp,Z2): bad.append((epd,fq))
        print('pruned',c.homeomorphism_type,c.normal_euler,fq,'cyclic profile' if fq in (Zp,Z2) else 'NOT CYCLIC',flush=True)
print('pruned by tube rule',pruned,'with 2-edge cut',cut,'kept',len(epds)-pruned-cut,'non-cyclic prunes',len(bad))
