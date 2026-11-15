"""Do diagram strings from the committed data mean the same diagram to the parser in this tree?
Prints, for each string, (vertices, component types and normal Euler numbers, |Hom(pi_1, S3)|, |Hom(pi_1, A4)|)
computed with indep_pi1_lib, so that the output of two trees can be diffed.  usage: python epd_compat.py SCRATCH"""
import sys, json
S = sys.argv[1]; sys.path.insert(0, S + '/mine')
from indep_pi1_lib import count_homs, GROUPS
from surfaceknots.diagrams import MarkedDiagram
a = json.load(open(S + '/mine/ch11_head.json'))
strings = [(e['label'], e['epd']) for e in a['entries']] + [(r['label'], r['epd']) for r in a['residual_collisions'][:25]]
t = json.load(open(S + '/sk/code/data/yoshikawa_ch10.json'))
ent = t if isinstance(t, list) else (t.get('entries') or t.get('table') or [])
strings += [(e['label'], e.get('epd') or e.get('canonical_epd')) for e in ent if (e.get('epd') or e.get('canonical_epd'))]
strings += [("CS1", "X[1,2,3,4], X[5,6,4,7], X[6,8,9,1], Y[9,10,11,2], X[12,13,14,3], Y[14,15,16,7], X[5,17,18,8], X[17,19,10,18], X[11,20,21,12], X[13,21,22,15], X[23,24,19,16], X[24,23,22,20]"),
            ("CS2", "X[1,2,3,4], X[4,5,6,1], X[7,8,9,6], Y[9,10,11,2], X[3,12,13,14], Y[14,15,16,17], X[18,7,5,17], X[8,18,19,10], X[20,21,12,11], X[21,22,15,13], X[23,24,19,16], X[24,23,22,20]")]
G = GROUPS if isinstance(GROUPS, dict) else dict(GROUPS)
for lbl, epd in strings:
    try:
        d = MarkedDiagram.from_epd(epd)
        topo = d.surface_topology()
        print(lbl, d.num_vertices, sorted((c.homeomorphism_type, c.normal_euler) for c in topo.components), count_homs(d, G['S3']), count_homs(d, G['A4']))
    except Exception as ex:
        print(lbl, 'ERROR', repr(ex)[:80])
