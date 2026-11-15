"""For each ch=12 two-saddle sphere in targeted_ch12_13_census.json: independent Wirtinger presentation ->
GAP Knuth-Bendix with generators x, u, v=xux^-1, w=xvx^-1.  Reports the ones whose completion is confluent
with u, v, w commuting (group Z |x Z^3, or a quotient thereof) and the action of x on w."""
import sys, json, subprocess, os, re
from concurrent.futures import ThreadPoolExecutor
S=sys.argv[1]; sys.path.insert(0,S+'/sk/code/src'); sys.path.insert(0,S+'/mine')
from indep_pi1_lib import relations
from surfaceknots.diagrams import MarkedDiagram
TAIL=open(S+'/mine/cs_indep_tail.g').read()
def gap_input(epd,j):
    d=MarkedDiagram.from_epd(epd); k,eid,rels=relations(d)
    exp=[None]*k; exp[0]=1; ch=True
    def val(x):
        e,s=eid[x]; return None if exp[e] is None else exp[e]*s
    while ch:
        ch=False
        for r in rels:
            a,b=r[1],r[2]
            if val(a) is not None and val(b) is None: e,s=eid[b]; exp[e]=-val(a)*s; ch=True
            elif val(b) is not None and val(a) is None: e,s=eid[a]; exp[e]=-val(b)*s; ch=True
    def g(x):
        e,s=eid[x]; return f"e{e+1}" if s==1 else f"e{e+1}^-1"
    R=[]
    for r in rels:
        if r[0]=='eq': R.append(f"{g(r[1])}*{g(r[2])}")
        else:
            _,out,inn,over,eps=r; O=g(over) if eps>0 else f"({g(over)})^-1"
            R.append(f"{g(out)}*{O}*{g(inn)}*({O})^-1")
    names=",".join(f'"e{i+1}"' for i in range(k))
    out=[f'F := FreeGroup("u","v","w","x",{names});;','u:=F.1;; v:=F.2;; w:=F.3;; x:=F.4;;']
    out+=[f"e{i+1}:=F.{i+5};;" for i in range(k)]
    out.append("R := [\n "+",\n ".join(R)+"\n];;")
    out.append(f"m := function(i) local es; es := {exp}; return GeneratorsOfGroup(F)[i+4]^es[i]; end;;")
    out.append("Add(R, x^-1*m(1));;")
    return "\n".join(out)+"\n"+TAIL.replace("JJ",str(j)), k
def run(job):
    tag,epd,j=job
    src,k=gap_input(epd,j)
    f=f"{S}/mine/csall/{tag}_j{j}.g"; open(f,'w').write(src)
    try:
        p=subprocess.run([S+'/gapenv/bin/gap','-q','--nointeract',f],capture_output=True,text=True,timeout=150)
        out=p.stdout
    except subprocess.TimeoutExpired: return tag,j,None
    m=re.search(r"confluent with (\d+) rules\s+\[u,v\]=1: (\w+)\s+\[u,w\]=1: (\w+)\s+\[v,w\]=1: (\w+)\s+x w x\^-1 = (\S+)",out)
    return tag,j,(m.groups() if m else out[-200:])
t=json.load(open(S+'/sk/code/data/complexity/targeted_ch12_13_census.json'))
items=[]
for r in t['ch12_m2_s2_knots']['new_prime_classes']:
    eps=r.get('all_epds') or [r['epd']]
    for i,e in enumerate(eps): items.append((f"{r['label']}.{i}", e, r['e1_generators']))
for i,r in enumerate(t['ch12_m2_s2_knots']['ambiguous_shadow_collisions']):
    items.append((f"amb{i}_{r['collision_matches'][0].split(':')[-1]}", r['epd'], r['e1_generators']))
print('diagrams',len(items),flush=True)
jobs=[(tag,epd,j) for tag,epd,e1 in items for j in (2,3,5,9)]
res={}
with ThreadPoolExecutor(max_workers=24) as ex:
    for tag,j,r in ex.map(run,jobs):
        if r is not None and isinstance(r,tuple): res.setdefault(tag,[]).append((j,r))
e1={tag:e for tag,_,e in items}
for tag,_,_ in items:
    if tag in res:
        j,(n,a,b,c,act)=res[tag][0]
        print(f"{tag:22s} E1={e1[tag]}  KB confluent ({n} rules, second meridian e{j}); u,v,w commute: {a},{b},{c}; x w x^-1 = {act}")
    else:
        print(f"{tag:22s} E1={e1[tag]}  no confluent system within 150 s")
