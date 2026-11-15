"""The ch=12 two-saddle spheres with Alexander ideal (3, x^2+1): what can be said about the group?
Independent Wirtinger presentation -> GAP.   usage: python f9_knot.py SCRATCH"""
import sys, json, subprocess
S = sys.argv[1]; sys.path.insert(0, S + '/sk/code/src'); sys.path.insert(0, S + '/mine')
from indep_pi1_lib import relations
from surfaceknots import MarkedDiagram
EPDS = {
 "amb2": "X[1,2,3,4], X[4,5,6,1], X[6,7,8,9], Y[9,10,11,2], X[12,13,14,3], Y[14,15,16,5], X[17,18,19,16], X[18,8,7,19], X[20,21,10,17], X[11,22,23,24], X[24,23,13,12], X[21,20,15,22]",
 "amb9": "X[1,2,3,4], X[4,5,6,1], X[7,8,9,6], X[10,11,2,9], X[11,12,13,3], X[13,14,15,16], Y[16,17,7,5], X[18,19,20,8], Z[20,21,12,10], X[19,22,14,21], X[23,24,17,15], X[24,23,22,18]",
}
for tag, epd in EPDS.items():
    d = MarkedDiagram.from_epd(epd); k, eid, rels = relations(d)
    exp = [None] * k; exp[0] = 1; ch = True
    def val(x):
        e, s = eid[x]; return None if exp[e] is None else exp[e] * s
    while ch:
        ch = False
        for r in rels:
            a, b = r[1], r[2]
            if val(a) is not None and val(b) is None: e, s = eid[b]; exp[e] = -val(a) * s; ch = True
            elif val(b) is not None and val(a) is None: e, s = eid[a]; exp[e] = -val(b) * s; ch = True
    def g(x):
        e, s = eid[x]; return f"F.{e+1}" if s == 1 else f"F.{e+1}^-1"
    Rl = []
    for r in rels:
        if r[0] == 'eq': Rl.append(f"{g(r[1])}*{g(r[2])}")
        else:
            _, out, inn, over, eps = r; O = g(over) if eps > 0 else f"({g(over)})^-1"
            Rl.append(f"{g(out)}*{O}*{g(inn)}*({O})^-1")
    src = f"""F := FreeGroup({k});; R := [ {', '.join(Rl)} ];; G := F / R;; ex := {exp};;
Print("abelianisation ", AbelianInvariants(G), "\\n");
for n in [2,3,4,5,6,8] do
  C := CyclicGroup(IsPermGroup, n);; z := C.1;;
  h := GroupHomomorphismByImages(G, C, GeneratorsOfGroup(G), List(ex, e -> z^e));;
  Print("H_1 of the ", n, "-fold cyclic cover of the complement: ", AbelianInvariants(Kernel(h)), "\\n");
od;
for n in [2,3,4,6,8,12] do
  Q := F / Concatenation(R, [F.1^n]);;
  sz := Size(Q);; Print("|G/<<m^", n, ">>| = ", sz);
  if sz < 5000 then P := Image(IsomorphismPermGroup(Q));; Print("   ", StructureDescription(P), "  derived ", Size(DerivedSubgroup(P)), " second derived ", Size(DerivedSubgroup(DerivedSubgroup(P)))); fi;
  Print("\\n");
od;
QUIT;
"""
    f = f"{S}/mine/f9/{tag}.g"; open(f, 'w').write(src)
    try:
        p = subprocess.run([S + '/gapenv/bin/gap', '-q', '--nointeract', f], capture_output=True, text=True, timeout=400)
        out = p.stdout
    except subprocess.TimeoutExpired as ex:
        out = (ex.stdout.decode() if ex.stdout else '') + '\n[stopped after 400 s]'
    print('==', tag); print('\n'.join(l for l in out.splitlines() if not l.startswith('#I') and l.strip()), flush=True)
