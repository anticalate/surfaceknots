import json
from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.moves import canonical_marked_key
from surfaceknots.unlink import certify_classical_unlink
from surfaceknots.yoshikawa_data import compute_full_invariant_record, YOSHIKAWA_BY_LABEL
from surfaceknots.complexity import certify_klein_bottle_subadditivity, MINIMAL_KLEIN_BOTTLE_E0_EPD
doc = "X[1,2,3,4], Y[4,2,5,6], Z[6,5,3,1]"
for name, epd in (("writeup", doc), ("code", MINIMAL_KLEIN_BOTTLE_E0_EPD)):
    d = MarkedDiagram.from_epd(epd)
    cp, cm = certify_classical_unlink(d.resolve(1)), certify_classical_unlink(d.resolve(-1))
    print(name, epd, "| L+:", cp.is_unlink, cp.reason, d.resolve(1).num_crossings, d.resolve(1).num_circles(), "| L-:", cm.is_unlink, cm.reason, d.resolve(-1).num_crossings, d.resolve(-1).num_circles(), "| c,m:", d.num_classical, d.num_marked)
    try:
        st = d.surface_topology()
        print("   comps:", [(c.homeomorphism_type, c.normal_euler) for c in st.components], "| c+/c-:", st.num_circles_plus, st.num_circles_minus)
        rec = compute_full_invariant_record(d, label=name, include_sublinks=False)
        print("   pi1:", rec["pi1_presentation"], rec["hom_counts"])
    except Exception as ex:
        print("   EXC", repr(ex)[:200])
cert = certify_klein_bottle_subadditivity()
print(json.dumps(cert, indent=1))
for epd in cert["move_path_epds"]:
    d = MarkedDiagram.from_epd(epd)
    st = d.surface_topology()
    rec = compute_full_invariant_record(d, label="x", include_sublinks=False)
    print(d.num_vertices, [(c.homeomorphism_type, c.normal_euler) for c in st.components], rec["pi1_presentation"])
