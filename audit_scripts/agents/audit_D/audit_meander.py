"""Is the greedy compiler's failure a genuine obstruction?  Exhaustive search for crossing-free band systems
(k = b - p bands, no closed loops) in the disk taking planar matching M- to planar matching M+."""
import itertools, sys
from collections import defaultdict
from surfaceknots.triplane import generate_catalan_matchings
from surfaceknots.triplane_census import _meander_marker_count

def meander_cycles(b, mm, mp):
    a = {}; c = {}
    for u, v in mm: a[u] = v; a[v] = u
    for u, v in mp: c[u] = v; c[v] = u
    seen = set(); n = 0
    for s in range(2 * b):
        if s in seen: continue
        n += 1; cur = s
        while cur not in seen:
            seen.add(cur); nx = a[cur]; seen.add(nx); cur = c[nx]
    return n

def regions(b, mm):
    """Region boundary cycles of the disk cut along M-: list of lists of (arc_index, direction) with region on the left.
    direction +1: traversed lo->hi, -1: hi->lo."""
    n = 2 * b
    part = {}; arc_of = {}
    for idx, (u, v) in enumerate(mm):
        part[u] = v; part[v] = u; arc_of[u] = idx; arc_of[v] = idx
    seen = set(); regs = []
    for x0 in range(n):
        if x0 in seen: continue
        cyc = []; x = x0
        while x not in seen:
            seen.add(x)
            s = (x + 1) % n          # walk boundary circle x -> x+1 (ccw), then arc from s to part[s]
            y = part[s]
            cyc.append((arc_of[s], 1 if s < y else -1))
            x = y
        regs.append(cyc)
    return regs

def band_systems(b, mm, k):
    """Yield (pairing, loops) for every crossing-free system of k bands attached to M- (up to isotopy)."""
    regs = regions(b, mm)
    sides = [(ri, pos) for ri, cyc in enumerate(regs) for pos in range(len(cyc))]  # arc-sides
    # a band: region ri and two side positions pos1 <= pos2 within region
    band_choices = [(ri, p1, p2) for ri, cyc in enumerate(regs) for p1 in range(len(cyc)) for p2 in range(p1, len(cyc))]
    lo_hi = [tuple(sorted(p)) for p in mm]
    for bands in itertools.combinations_with_replacement(band_choices, k):
        # attachments: (band j, end e) on arc A at side (ri,pos)
        att = []
        for j, (ri, p1, p2) in enumerate(bands):
            att.append((j, 0, ri, p1)); att.append((j, 1, ri, p2))
        by_arc = defaultdict(list)
        for a in att:
            by_arc[regs[a[2]][a[3]][0]].append(a)
        arcs = sorted(by_arc)
        for perms in itertools.product(*(itertools.permutations(by_arc[A]) for A in arcs)):
            # rank of each attachment along its arc (lo -> hi)
            rank = {}
            for A, pm in zip(arcs, perms):
                for r, a in enumerate(pm): rank[(a[0], a[1])] = r
            # non-crossing check within each region: order attachments along region traversal
            ok = True
            for ri, cyc in enumerate(regs):
                seq = []
                for pos, (A, d) in enumerate(cyc):
                    here = [a for a in att if a[2] == ri and a[3] == pos]
                    here.sort(key=lambda a: rank[(a[0], a[1])] * d)
                    seq += [a[0] for a in here]
                # each band appears twice; check planar (non-crossing) & discard symmetric duplicates lazily
                st = []
                for j in seq:
                    if st and st[-1] == j: st.pop()
                    else: st.append(j)
                # non-crossing iff fully reducible by adjacent cancellation
                if st: ok = False; break
            if not ok: continue
            # build pieces: arc A cut into len+1 pieces 0..m ; piece q lies between attachment rank q-1 and q
            # nodes: (A, q).  endpoints: lo -> (A,0), hi -> (A,m)
            nxt = defaultdict(list)
            def before_after(a):
                A, d = regs[a[2]][a[3]]
                r = rank[(a[0], a[1])]
                # traversal direction d: before = piece r (if d=+1) else piece r+1
                return ((A, r), (A, r + 1)) if d == 1 else ((A, r + 1), (A, r))
            for j in range(k):
                P = att[2 * j]; Q = att[2 * j + 1]
                # ensure P precedes Q in region traversal (for self-side bands need order by rank*d)
                A_, d_ = regs[P[2]][P[3]]
                if P[3] == Q[3] and rank[(P[0], P[1])] * d_ > rank[(Q[0], Q[1])] * d_:
                    P, Q = Q, P
                bP, aP = before_after(P); bQ, aQ = before_after(Q)
                nxt[aP].append(bQ); nxt[bQ].append(aP)
                nxt[bP].append(aQ); nxt[aQ].append(bP)
            # each interior piece has 2 band connections, end pieces have 1 (+ endpoint)
            m_of = {A: len(by_arc.get(A, [])) for A in range(b)}
            start = {}
            for A, (lo, hi) in enumerate(lo_hi):
                start[lo] = (A, 0); start[hi] = (A, m_of[A])
            end_of = {v: kk for kk, v in start.items()}
            visited = set(); pairing = {}
            valid = True
            for p0 in range(2 * b):
                if p0 in pairing: continue
                node = start[p0]; prev = None
                if m_of[node[0]] == 0:
                    visited.add(node)
                    pairing[p0] = [q for q in lo_hi[node[0]] if q != p0][0]; pairing[pairing[p0]] = p0; continue
                while True:
                    visited.add(node)
                    cand = [x for x in nxt[node]]
                    if prev is not None and prev in cand: cand.remove(prev)
                    if not cand:
                        break
                    prev, node = node, cand[0]
                    # arriving at piece `node` via a band; if it is an end piece (and not same as start w/ single), finish
                    if (node[1] == 0 or node[1] == m_of[node[0]]) and len(nxt[node]) == 1:
                        visited.add(node); break
                A, q = node
                endp = lo_hi[A][0] if q == 0 else lo_hi[A][1]
                pairing[p0] = endp; pairing[endp] = p0
            allnodes = {(A, q) for A in range(b) for q in range(m_of[A] + 1)}
            rest = allnodes - visited
            loops = 0
            while rest:
                loops += 1
                node = rest.pop(); stack = [node]
                while stack:
                    x = stack.pop()
                    for y in nxt[x]:
                        if y in rest: rest.discard(y); stack.append(y)
            yield tuple(sorted(tuple(sorted((u, pairing[u]))) for u in range(2 * b) if u < pairing[u])), loops

if __name__ == "__main__":
    for b in (2, 3, 4):
        cats = generate_catalan_matchings(b)
        reach = {}
        for mm in cats:
            for k in range(0, b):
                res = defaultdict(set)
                for pairing, loops in band_systems(b, mm, k):
                    res[pairing].add(loops)
                reach[(mm, k)] = res
        agree = disagree = fail = 0; examples = []; extra = []
        for mm in cats:
            for mp in cats:
                p = meander_cycles(b, mm, mp)
                k = b - p
                greedy = _meander_marker_count(b, mm, mp)
                exists = 0 in reach[(mm, k)].get(tuple(sorted(tuple(sorted(x)) for x in mp)), set())
                if greedy is None: fail += 1
                if (greedy is not None) == exists:
                    agree += 1
                    if greedy is not None: assert greedy == k, (greedy, k)
                else:
                    disagree += 1; examples.append((mm, mp, greedy, exists))
                if greedy is None and b == 4:
                    # any crossing-free system with more bands (k+1 <= b-1) and closed loops?
                    more = {kk: sorted(reach[(mm, kk)].get(tuple(sorted(tuple(sorted(x)) for x in mp)), set())) for kk in range(k, b)}
                    extra.append((mm, mp, p, more))
        print(f"b={b}: {len(cats)**2} ordered pairs; greedy failures={fail}; greedy==exhaustive(no-loop, k=b-p bands): agree={agree} disagree={disagree}")
        for e in examples[:10]: print("   DISAGREE", e)
        for e in extra: print("   obstructed pair", e[0], e[1], "meander cycles p=%d" % e[2], "crossing-free band systems by #bands -> loop counts:", e[3])
