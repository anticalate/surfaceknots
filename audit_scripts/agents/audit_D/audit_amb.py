import json
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.yoshikawa_data import compute_full_invariant_record, invariant_discrimination_key, YOSHIKAWA_BY_LABEL
d=json.load(open('data/complexity/targeted_ch12_13_census.json'))
s=d['ch12_m2_s2_knots']
yosh = json.load(open('data/yoshikawa_ch10.json'))
ch11 = json.load(open('data/expanded_ch11.json'))
print(ch11.keys())
cat = {}
for e in yosh: cat.setdefault(str(invariant_discrimination_key(e)), []).append(e['label'])
for e in ch11['entries']: cat.setdefault(str(invariant_discrimination_key(e)), []).append(e['label'])
for a in s['ambiguous_shadow_collisions']:
    if a['e1_generators']==['3','x**2 + 1']:
        md = MarkedDiagram.from_epd(a['epd'])
        rec = compute_full_invariant_record(md, label='amb', include_sublinks=False)
        print(md.num_classical, md.num_marked, rec['component_types'], rec['e1_generators'], rec['e2_generators'], rec['pi1_presentation'], rec['hom_counts'], 'S2', rec['sigma2_h1_torsion'], rec.get('cover_betti_profile'))
        print('  in ch<=11 catalog:', cat.get(str(invariant_discrimination_key(rec))))
u = compute_full_invariant_record(MarkedDiagram.unknot_sphere(), label='0_1', include_sublinks=False)
print('unknot', u['e1_generators'], u['hom_counts'])
for lbl in ('10_3','10_2'):
    r=[e for e in yosh if e['label']==lbl][0]
    print(lbl, r['e1_generators'], r['hom_counts'], r['pi1_presentation'])
from surfaceknots.yoshikawa_data import invariant_filtration_layers
a=[x for x in s['ambiguous_shadow_collisions'] if x['e1_generators']==['3','x**2 + 1']][0]
rec = compute_full_invariant_record(MarkedDiagram.from_epd(a['epd']), label='amb', include_sublinks=False)
for (n1,v1),(n2,v2) in zip(invariant_filtration_layers(rec), invariant_filtration_layers(u)):
    print(n1, v1==v2, str(v1)[:150], '|', str(v2)[:150])
print(rec['e1_fingerprint_hash'], u['e1_fingerprint_hash'], rec['single_var_alex'], rec['alexander_biquandle_s'], u['alexander_biquandle_s'], rec['kei_invariant_factors'], rec['bikei_invariant_factors'])
