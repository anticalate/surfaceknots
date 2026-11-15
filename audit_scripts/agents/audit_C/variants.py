"""Generate move-equivalent diagrams.  'safe' moves: relabel/rotate, Omega1/6 stabilisation, Omega2 poke
(repo generators, verified by reading), plus optionally repo 'level' moves (Omega3,4,5,7,8, tube flip)."""
import random
from surfaceknots import MarkedDiagram
from surfaceknots.moves import (list_omega1_6_stabilizations, list_omega2_pokes,
    list_level_yoshikawa_neighbors, _list_omega3_and_4_neighbors, _list_omega5_neighbors,
    _list_omega7_neighbors, _list_omega8_neighbors, _list_marked_tube_flip_neighbors)

ROT_CL = {"X0": "X1", "X1": "X0", "Y": "Z", "Z": "Y"}

def relabel(d, rng):
    n = d.num_vertices
    perm = list(range(n)); rng.shuffle(perm)       # old v -> new perm[v]
    rot = [rng.randrange(4) for _ in range(n)]
    def nd(hd):
        v, s = hd >> 2, hd & 3
        return 4 * perm[v] + ((s + rot[v]) & 3)
    opp = [0] * (4 * n); vt = [None] * n
    for hd in range(4 * n):
        opp[nd(hd)] = nd(d.opp[hd])
    for v in range(n):
        vt[perm[v]] = d.vtype[v] if rot[v] % 2 == 0 else ROT_CL[d.vtype[v]]
    return MarkedDiagram(n, tuple(opp), tuple(vt), d.free_loops)

def is_planar(d):
    n = d.num_vertices
    if n == 0: return True
    # count faces: face walk  cur -> opp -> slot-1
    seen = [False] * (4 * n); F = 0
    for s in range(4 * n):
        if seen[s]: continue
        F += 1; cur = s
        while not seen[cur]:
            seen[cur] = True
            od = d.opp[cur]
            cur = ((od >> 2) << 2) | (((od & 3) - 1) & 3)
    # components
    par = list(range(n))
    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for hd in range(4 * n):
        a, b = f(hd >> 2), f(d.opp[hd] >> 2)
        if a != b: par[a] = b
    C = len({f(v) for v in range(n)})
    return n - 2 * n + F == 2 * C

MOVES = {
 "O1": lambda d: [x for x in list_omega1_6_stabilizations(d) if x.vtype[-1] in ("X0", "X1")],
 "O6": lambda d: [x for x in list_omega1_6_stabilizations(d) if x.vtype[-1] in ("Y", "Z")],
 "O2": list_omega2_pokes,
 "O34": _list_omega3_and_4_neighbors,
 "O5": _list_omega5_neighbors,
 "O7": _list_omega7_neighbors,
 "O8": _list_omega8_neighbors,
 "TF": _list_marked_tube_flip_neighbors,
}
SAFE = ["O1", "O6", "O2"]
LEVEL = ["O34", "O5", "O7", "O8", "TF"]

def random_variant(d, rng, nsteps=4, max_v=15, use_level=True, allowed=None):
    cur = d; path = []
    for _ in range(nsteps):
        names = list(allowed) if allowed else (SAFE + (LEVEL if use_level else []))
        rng.shuffle(names)
        if cur.num_vertices >= max_v:
            names = [x for x in names if x in LEVEL]
        done = False
        for nm in names:
            if nm == "O2" and cur.num_vertices + 2 > max_v: continue
            cands = MOVES[nm](cur)
            if cands:
                cur = rng.choice(cands); path.append(nm); done = True
                break
        if not done: break
    cur = relabel(cur, rng)
    return cur, path
