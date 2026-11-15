import json, sys
from surfaceknots import MarkedDiagram, YOSHIKAWA_TABLE
from surfaceknots.moves import _is_planar_4reg_diagram, has_certified_cyclic_unknotted_tube
from surfaceknots.unlink import certify_classical_unlink
from surfaceknots.invariants.fundamental_group import wirtinger_presentation, finite_quotient_profile
GR=("S3","Q8","A4","S4","A5")
def homs(d):
    return dict(finite_quotient_profile(wirtinger_presentation(d),GR))
def options(d):
    n=d.num_vertices; st=d.surface_topology(); res=[]
    for u in range(n):
        if d.vtype[u] not in ("Y","Z"): continue
        for su in range(4):
            d_a=d.opp[4*u+su]; d_b=d.opp[4*u+((su+1)&3)]; w=d_a>>2
            if w<=u or (d_b>>2)!=w or d.vtype[w] not in ("Y","Z"): continue
            sw_a=d_a&3; sw_b=d_b&3
            if ((sw_b+1)&3)!=sw_a: continue
            p0,p1=4*u+((su+2)&3),4*u+((su+3)&3); p2,p3=4*w+((sw_b+2)&3),4*w+((sw_b+3)&3)
            remap=[v-(1 if v>u else 0)-(1 if v>w else 0) for v in range(n)]; m=n-2
            nvt=tuple(d.vtype[v] for v in range(n) if v not in (u,w))
            for oi,pairs in enumerate(([(p0,p1),(p2,p3)],[(p0,p3),(p1,p2)])):
                partner={x:y for x,y in pairs}|{y:x for x,y in pairs}
                nopp=[0]*(4*m); valid=True
                for v in range(n):
                    if v in (u,w): continue
                    nv=remap[v]
                    for s in range(4):
                        od=d.opp[4*v+s]; steps=0
                        while (od>>2) in (u,w):
                            od=d.opp[partner[od]]; steps+=1
                            if steps>8: valid=False; break
                        if not valid: break
                        nopp[4*nv+s]=(remap[od>>2]<<2)|(od&3)
                    if not valid: break
                if not valid: continue
                cand=MarkedDiagram(num_vertices=m,opp=tuple(nopp),vtype=nvt,free_loops=0)
                if not _is_planar_4reg_diagram(cand): continue
                if not (certify_classical_unlink(cand.resolve(1)).is_unlink and certify_classical_unlink(cand.resolve(-1)).is_unlink): continue
                cst=cand.surface_topology()
                if cst.num_components!=1 or cst.euler_characteristic!=st.euler_characteristic+2: continue
                res.append((u,w,d.vtype[u],d.vtype[w],oi,cand))
    return res
