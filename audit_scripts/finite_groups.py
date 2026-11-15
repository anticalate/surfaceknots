"""Orders of the groups of the table / ch=11 entries all of whose components are non-orientable
(the only entries whose group can be finite: H_1 = Z^(orientable comps) + Z_2^(non-orientable comps)).
Independent Wirtinger presentation (indep_pi1_lib) -> GAP: simplified presentation, Size, structure,
order of a meridian.   usage: python finite_groups.py SCRATCH"""
import sys, json, subprocess
S = sys.argv[1]; sys.path.insert(0, S + '/sk/code/src'); sys.path.insert(0, S + '/mine')
from indep_pi1_lib import relations
from surfaceknots import MarkedDiagram, YOSHIKAWA_TABLE
c11 = json.load(open(S + '/sk/code/data/expanded_ch11.json'))
pool = [(e.label, e.canonical_epd) for e in YOSHIKAWA_TABLE if e.canonical_epd] + [(r['label'], r['epd']) for r in c11['entries']]
for lbl, epd in pool:
    d = MarkedDiagram.from_epd(epd)
    if d.num_vertices == 0: continue
    topo = d.surface_topology()
    if any(c.orientable for c in topo.components): continue
    k, eid, rels = relations(d)
    def g(x):
        e, s = eid[x]; return f"F.{e+1}" if s == 1 else f"F.{e+1}^-1"
    Rl = []
    for r in rels:
        if r[0] == 'eq': Rl.append(f"{g(r[1])}*{g(r[2])}")
        else:
            _, out, inn, over, eps = r; O = g(over) if eps > 0 else f"({g(over)})^-1"
            Rl.append(f"{g(out)}*{O}*{g(inn)}*({O})^-1")
    src = f"""F := FreeGroup({k});;
G := F / [ {', '.join(Rl)} ];;
P := PresentationFpGroup(G);; TzOptions(P).printLevel := 0;; TzGoGo(P);;
H := FpGroupPresentation(P);;
Print("gens ", Length(GeneratorsOfGroup(H)), " rels ", Length(RelatorsOfFpGroup(H)), "\\n");
Print("abelian invariants ", AbelianInvariants(G), "\\n");
sz := Size(H);; Print("order ", sz, "\\n");
if sz < 100000 then
  Q := Image(IsomorphismPermGroup(H));; Print("structure ", StructureDescription(Q), "\\n");
  iso := IsomorphismPermGroup(G);;
  Print("meridian orders ", Collected(List(GeneratorsOfGroup(G), x -> Order(Image(iso, x)))), "\\n");
  Print("centre ", Size(Centre(Q)), " derived ", Size(DerivedSubgroup(Q)), "\\n");
fi;
QUIT;
"""
    f = f"{S}/mine/finite/{lbl.replace('^','_').replace(',','_')}.g"; open(f, 'w').write(src)
    try:
        p = subprocess.run([S + '/gapenv/bin/gap', '-q', '--nointeract', f], capture_output=True, text=True, timeout=240)
        out = p.stdout.strip() + (('\n  stderr: ' + p.stderr.strip()[-300:]) if p.stderr.strip() else '')
    except subprocess.TimeoutExpired:
        out = 'no answer in 240 s (coset enumeration did not finish: the group may be infinite)'
    print(lbl, [ (c.homeomorphism_type, c.normal_euler) for c in topo.components]); print('  ' + out.replace('\n', '\n  '), flush=True)
