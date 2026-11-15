"""Weak primeness (Yoshikawa Def 5.1) of the ch=11 classes with >= 2 components.
F is weakly prime iff it is not F1 # F2 with ch(F1), ch(F2) <= 10.  For each class we list every way its
topological type can arise as a sum of two surfaces of ch <= 10 (using: the only non-split links with
ch <= 10 are table entries and their sums with knots of total ch <= 10; blind runs at ch 9, 10 show no
other non-cyclic profile), build the sums, and compare homomorphism counts into S3, Q8, A4, S4, A5-free."""
import sys, json, itertools
sys.path.insert(0, sys.argv[1])
from surfaceknots import MarkedDiagram, YOSHIKAWA_BY_LABEL
from surfaceknots.ch11_census import planar_connected_sum
from surfaceknots.unlink import certify_classical_unlink
from surfaceknots.invariants.fundamental_group import wirtinger_presentation, finite_quotient_profile
GR=("S3","Q8","A4","S4","D8","D10","F21")
def homs(d): return tuple(finite_quotient_profile(wirtinger_presentation(d),GR)[g] for g in GR)
def D(lbl): return MarkedDiagram.from_epd(YOSHIKAWA_BY_LABEL[lbl].canonical_epd)
def comp_darts(d):
    comp,nc=d._dart_surface_components(); st=d.surface_topology()
    out={}
    for x in range(4*d.num_vertices): out.setdefault(comp[x],x)
    return [(st.components[c].homeomorphism_type,st.components[c].normal_euler,x) for c,x in sorted(out.items())]
def typ(d):
    st=d.surface_topology(); return tuple(sorted((c.homeomorphism_type,abs(c.normal_euler)) for c in st.components))
def sums(d1,d2):
    for t1,e1,x1 in comp_darts(d1):
        for t2,e2,x2 in comp_darts(d2):
            for fl in (False,True):
                for mir in (False,True):
                    dd2=d2.mirror_3d() if mir else d2
                    s=planar_connected_sum(d1,x1,dd2,x2,flip_splice=fl)
                    yield (t1,t2,fl,mir),s
P=D('2_1^-1'); T=D('2_1^1')
K0=planar_connected_sum(P,0,P.mirror_3d(),0)            # unknotted Klein bottle, e = 0
knots2=['8_1','9_1','10_1','10_2','10_3']
tori=['2_1^1','10_1^1']
table_links=[l for l in YOSHIKAWA_BY_LABEL if ',' in l]
# all surfaces of ch <= 10 that can serve as a summand: table entries, unknotted T^2, P+-, K^2(e=0), 8_1#T^2, 8_1#P
summands=[(l,D(l)) for l in YOSHIKAWA_BY_LABEL if l!='0_1']+[('K2std',K0)]
summands+=[('8_1#T2',planar_connected_sum(D('8_1'),0,T,0)),('8_1#P',planar_connected_sum(D('8_1'),0,P,0))]
ch11=json.load(open('data/expanded_ch11.json'))
targets={r['label']:MarkedDiagram.from_epd(r['epd']) for r in ch11['entries']}
tinfo={l:(typ(d),homs(d)) for l,d in targets.items()}
by_type={}
for l,(t,h) in tinfo.items(): by_type.setdefault(t,[]).append(l)
print('targets by type:',{str(k):v for k,v in by_type.items()},flush=True)
checked=0; hits=[]
cnt={l:0 for l in targets}
for (l1,d1),(l2,d2) in itertools.combinations_with_replacement(summands,2):
    if d1.num_vertices+d2.num_vertices>22: continue
    seen=set()
    for key,s in sums(d1,d2):
        t=typ(s)
        if t not in by_type: continue
        if not (certify_classical_unlink(s.resolve(1)).is_unlink and certify_classical_unlink(s.resolve(-1)).is_unlink): 
            print('non-admissible sum?',l1,l2,key); continue
        h=homs(s); checked+=1
        for l in by_type[t]:
            cnt[l]+=1
            if h==tinfo[l][1]: hits.append((l,l1,l2,key)); print('HOM MATCH',l,'=?',l1,'#',l2,key,flush=True)
print('sums with a target topology checked:',checked)
print('comparisons per target:',cnt)
print('hom-count matches:',len(hits))
