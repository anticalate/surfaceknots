import sys, json
from surfaceknots import MarkedDiagram, YOSHIKAWA_BY_LABEL
from surfaceknots.ch11_census import planar_connected_sum
from surfaceknots.yoshikawa_data import compute_full_invariant_record, invariant_discrimination_key, invariant_filtration_layers
from surfaceknots.invariants.linking import sublink_summary
def D(l): return MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL[l].canonical_epd)
def comp_dart(d,typ):
    comp,nc=d._dart_surface_components(); st=d.surface_topology()
    for x in range(4*d.num_vertices):
        if st.components[comp[x]].homeomorphism_type==typ: return x
ch11={r['label']:r for r in json.load(open('data/expanded_ch11.json'))['entries']}
cases=[('11_2^0,1','6_1^0,1','T^2'),('11_2^0,1','6_1^0,1','S^2'),('11_1^0,-2','7_1^0,-2','K^2'),('11_1^0,-2','7_1^0,-2','S^2')]
for tgt,link,ctyp in cases:
    d=D(link); s=planar_connected_sum(d,comp_dart(d,ctyp),D('8_1'),0)
    rec=compute_full_invariant_record(s,label='sum',include_sublinks=True)
    t=ch11[tgt]
    same_key = invariant_discrimination_key(rec)==invariant_discrimination_key(t)
    La=dict(invariant_filtration_layers(rec)); Lb=dict(invariant_filtration_layers(t))
    diff=[k for k in La if repr(La[k])!=repr(Lb[k])]
    print(f"{tgt} vs {link} #_{ctyp} 8_1: same repo key = {same_key}; layers that differ: {diff}")
    print("    sum sublinks:", [(p['types'],p.get('min_ch_found'),p.get('certified_split')) for p in rec.get('sublinks',{}).get('proper_sublinks',[])][:6] if isinstance(rec.get('sublinks'),dict) else rec.get('sublinks'))
    print("    tgt sublinks:", [(p['types'],p.get('min_ch_found'),p.get('certified_split')) for p in t.get('sublinks',{}).get('proper_sublinks',[])][:6] if isinstance(t.get('sublinks'),dict) else t.get('sublinks'))
    sys.stdout.flush()
