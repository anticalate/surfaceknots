"""Audit helpers: instrumented C lib wrappers."""
import ctypes, os
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.unlink import certify_classical_unlink
HERE = os.path.dirname(os.path.abspath(__file__))
L = ctypes.CDLL(os.path.join(HERE, "fe_inst.so"))
MAXN = 20
VT = ("X0", "X1", "Y", "Z")
VTI = {"X0": 0, "X1": 1, "Y": 2, "Z": 3}
IP = ctypes.POINTER(ctypes.c_int)

def arrs(d):
    n = d.num_vertices
    return (ctypes.c_int * (4 * n))(*d.opp), (ctypes.c_int * n)(*[VTI[v] for v in d.vtype])

def mk(n, opp, vt):
    return MarkedDiagram(num_vertices=n, opp=tuple(int(x) for x in opp), vtype=tuple(VT[int(v)] for v in vt), free_loops=0)

def level_neighbors(d, maxn=512):
    n = d.num_vertices
    o, v = arrs(d)
    oo = (ctypes.c_int * (maxn * 4 * n))()
    ov = (ctypes.c_int * (maxn * n))()
    c = L.x_level_neighbors(n, o, v, oo, ov, maxn)
    return [mk(n, oo[i*4*n:(i+1)*4*n], ov[i*n:(i+1)*n]) for i in range(c)]

def canon_key(d):
    n = d.num_vertices
    o, v = arrs(d)
    w = (ctypes.c_uint32 * n)()
    L.x_canon_key(n, o, v, w)
    return tuple(w)

def greedy_reduce(d):
    n = d.num_vertices
    o, v = arrs(d)
    oo = (ctypes.c_int * (4 * n))(); ov = (ctypes.c_int * n)(); lp = ctypes.c_int(0)
    m = L.x_greedy_reduce(n, o, v, oo, ov, ctypes.byref(lp))
    return MarkedDiagram(num_vertices=m, opp=tuple(oo[:4*m]), vtype=tuple(VT[x] for x in ov[:m]), free_loops=lp.value)

def has_cut(d):
    o, v = arrs(d); return L.x_has_cut(d.num_vertices, o)
def has_imm(d):
    o, v = arrs(d); return L.x_has_imm(d.num_vertices, o, v)
def has_1step(d, chk):
    o, v = arrs(d); return L.x_has_1step(d.num_vertices, o, v, chk)

def admissible(d):
    return certify_classical_unlink(d.resolve(1)).is_unlink and certify_classical_unlink(d.resolve(-1)).is_unlink

def is_planar(d):
    n = d.num_vertices; opp = d.opp
    seen = [0]*(4*n); F = 0
    for s in range(4*n):
        if seen[s]: continue
        F += 1; cur = s
        while not seen[cur]:
            seen[cur] = 1; od = opp[cur]; cur = (od//4)*4 + ((od % 4) - 1) % 4
    # connected components
    comp = 0; vis = [0]*n
    for v0 in range(n):
        if vis[v0]: continue
        comp += 1; st = [v0]; vis[v0] = 1
        while st:
            u = st.pop()
            for s in range(4):
                w = opp[4*u+s]//4
                if not vis[w]: vis[w] = 1; st.append(w)
    return n - 2*n + F == 2*comp

_FQ = None
def inv(d):
    """mirror-insensitive invariant: topology + finite quotient hom counts"""
    from surfaceknots.invariants.fundamental_group import wirtinger_presentation, finite_quotient_profile
    topo = d.surface_topology()
    comps = tuple(sorted((c.genus_symbol, abs(c.normal_euler)) for c in topo.components))
    wp = wirtinger_presentation(d, from_resolution=-1)
    fq = finite_quotient_profile(wp)
    if isinstance(fq, dict): fq = tuple(sorted(fq.items()))
    return (comps, fq)
