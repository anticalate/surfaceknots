"""Check the inputs of the third-party hand proof for the ch=12 diagram with group Z |x Z^3
(index 191 of the m=2 survivors = the first diagram of AUDIT section 10).
The repository's simplified presentation <g1, g2 | R0..R3> is rewritten by Reidemeister-Schreier with
g1 = t, g2 = t v_{-1}, v_k = t^k v_0 t^-k.  The hand proof uses three of the resulting relations:
  (E1)  v1^-1 v0^-1 v1 v2^-1 v0 v2 = 1          (a conjugate of R0)
  (R1') v3 = v1^-1 v0 v2 v1^-1                  (a conjugate of R1)
  (E2)  v1^-1 (v2 v0) v1 (v2 v0)^-1 = 1         (from a conjugate of R2 after eliminating v3 or v_-1)
Run from code/ with PYTHONPATH=src (original code)."""
import json, sys
from surfaceknots import MarkedDiagram
from surfaceknots.invariants.fundamental_group import wirtinger_presentation
EPD = "X[1,2,3,4], X[5,6,4,7], X[6,8,9,1], Y[9,10,11,2], X[12,13,14,3], Y[14,15,16,7], X[5,17,18,8], X[17,19,10,18], X[11,20,21,12], X[13,21,22,15], X[23,24,19,16], X[24,23,22,20]"
p = wirtinger_presentation(MarkedDiagram.from_epd(EPD), simplify=True)
print("generators", p.num_generators, "relators", [tuple(r) for r in p.relators])
def rs(w):
    k = 0; out = []
    for l in w:
        if l == 1: k += 1
        elif l == -1: k -= 1
        elif l == 2: out.append((k, 1)); k += 1
        elif l == -2: out.append((k - 1, -1)); k -= 1
    assert k == 0
    m = min(i for i, _ in out)
    return [(i - m, e) for i, e in out], m
def show(w): return " ".join(f"v{i}" + ("" if e == 1 else "^-1") for i, e in w)
def cyc_class(w):
    w = list(w)
    ch = True
    while ch and len(w) > 1:
        ch = False
        if w[0][0] == w[-1][0] and w[0][1] == -w[-1][1]: w = w[1:-1]; ch = True
    rots = [tuple(w[i:] + w[:i]) for i in range(len(w))]
    invs = [tuple((i, -e) for i, e in reversed(r)) for r in rots]
    return min(rots + invs)
E1 = [(1,-1),(0,-1),(1,1),(2,-1),(0,1),(2,1)]
R1 = [(3,-1),(1,-1),(0,1),(2,1),(1,-1)]
found = {}
for n, r in enumerate(p.relators):
    w, m = rs(tuple(r))
    print(f"R{n} -> (shift {m})  {show(w)}")
    c = cyc_class(w)
    if c == cyc_class(E1): found["E1"] = n
    if c == cyc_class(R1): found["R1'"] = n
print("matches, up to cyclic rotation, inversion and shift of indices:", found)
