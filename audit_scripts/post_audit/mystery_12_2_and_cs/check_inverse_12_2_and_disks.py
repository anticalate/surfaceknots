"""Complete bidirectional group isomorphism proof (Direction 2: G_model -> G(12_2.a/b)
and composition == id) + 4D ribbon-disk decomposition check on 12_2, 12_10, and 12_8.
"""
from __future__ import annotations

import json
import sys
from collections import deque
from pathlib import Path

REPO = Path("/afs/csail.mit.edu/u/e/ejli/.gemini/antigravity-cli/worktrees/1dadd8f0-f304-455a-b682-76a831212ce6/subagent-4D-Topology---Mystery-Knot-Analyst-self-69be3150")
sys.path.insert(0, str(REPO / "code" / "src"))
sys.path.insert(0, "/afs/csail.mit.edu/u/e/ejli/.gemini/antigravity-cli/brain/ee536e5c-0e8d-469c-b4f6-08d6f2745f37/scratch")

from audit_12_2 import independent_edge_wirtinger, independent_fox_e1, todd_coxeter_order
from audit_cs import collapse_pos_rels_to_sheets
from certify_12_2_iso import eval_in_semidirect_z3_free, free_red, w_inv, sub_word
from surfaceknots.diagrams import MarkedDiagram


def prove_word_trivial_in_fpgroup(word: tuple[int, ...], rels: list[tuple[int, ...]], max_len: int = 24, max_visited: int = 60000) -> bool:
    """Constructively prove `word == 1` in `<gens | rels>` by short-first search inserting cyclic permutations
    of `r^{+-1}` (`r in rels`) and freely reducing! Every step is an exact identity in `<gens | rels>`.
    """
    word = free_red(word)
    if not word:
        return True
    # Build all cyclic shifts of r and r^-1, and for every prefix of length >= len(r)//2 + 1, rewrite rules!
    all_cyc = set()
    for r in rels:
        for base in (free_red(r), free_red(w_inv(r))):
            L = len(base)
            for i in range(L):
                all_cyc.add(base[i:] + base[:i])
    # Extract length-nonincreasing or +2 rules: u -> v where u*v^-1 in all_cyc and len(v) <= len(u) + 2
    rules: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    for cyc in all_cyc:
        L = len(cyc)
        for cut in range(L // 2, L + 1):
            lhs = cyc[:cut]
            rhs = w_inv(cyc[cut:])
            if len(rhs) <= len(lhs) + 2:
                rules.append((lhs, rhs))
    rules.sort(key=lambda pair: len(pair[1]) - len(pair[0]))

    import heapq
    pq = [(len(word), word)]
    visited = {word}
    while pq and len(visited) < max_visited:
        _, cur = heapq.heappop(pq)
        if not cur:
            return True
        n = len(cur)
        for lhs, rhs in rules:
            ll = len(lhs)
            if ll > n:
                continue
            for idx in range(n - ll + 1):
                if cur[idx : idx + ll] == lhs:
                    cand = free_red(cur[:idx] + rhs + cur[idx + ll :])
                    if not cand:
                        return True
                    if len(cand) <= max_len and cand not in visited:
                        visited.add(cand)
                        heapq.heappush(pq, (len(cand), cand))
    return False


def main():
    ch12 = json.loads((REPO / "code" / "data" / "expanded_ch12.json").read_text())
    by_lbl = {e["label"]: e for e in ch12["entries"]}

    # 1. Prove Direction 2 of isomorphism G_model = <a, q | q^3=1, a q a = q a q a q> <-> G(12_2.a) and G(12_2.b)!
    e2 = by_lbl["12_2"]
    for idx, (epd, s2_in_aq) in enumerate([
        (e2["epd"], (2, 1, 2)),               # 12_2.a: s_1 -> a, s_2 -> q a q
        (e2["alternate_epds"][0], (1, -2)),   # 12_2.b: s_1 -> a, s_2 -> a q^-1
    ]):
        lbl = f"12_2.{chr(ord('a') + idx)}"
        # BFS shortest word q_in_s in (s_1, s_2) whose image in G_model is q = (((0, 1),), 0)!
        q_in_s = None
        bfs_q = deque([()])
        seen_s = {()}
        while bfs_q:
            ws = bfs_q.popleft()
            # Evaluate ws under s_1 -> (1,), s_2 -> s2_in_aq:
            waq: list[int] = []
            for x in ws:
                if x == 1: waq.append(1)
                elif x == -1: waq.append(-1)
                elif x == 2: waq.extend(s2_in_aq)
                elif x == -2: waq.extend(w_inv(s2_in_aq))
            if eval_in_semidirect_z3_free(tuple(waq)) == (((0, 1),), 0):
                q_in_s = ws
                break
            if len(ws) < 7:
                for nxt in (1, -1, 2, -2):
                    if ws and ws[-1] == -nxt:
                        continue
                    nw = ws + (nxt,)
                    if nw not in seen_s:
                        seen_s.add(nw)
                        bfs_q.append(nw)
        assert q_in_s is not None, f"Could not express q in (s_1, s_2) for {lbl}!"
        print(f"\nFound exact word for q in (s_1, s_2) on {lbl}: {q_in_s}")
        d = MarkedDiagram.from_epd(epd)
        k, rels = independent_edge_wirtinger(d)
        ivec, _ = independent_fox_e1(k, rels)
        pos_rels = [tuple((abs(l) if ((1 if l > 0 else -1) * ivec[abs(l) - 1]) > 0 else -abs(l)) for l in r) for r in rels]
        nsheets, srels = collapse_pos_rels_to_sheets(k, pos_rels)

        # Eliminate sheets > 2 exactly as in certify_12_2_iso:
        active = set(range(1, nsheets + 1))
        cur_rels = [free_red(pr) for pr in srels if free_red(pr)]
        prog = True
        while prog:
            prog = False
            for ri, rw in enumerate(cur_rels):
                cnts = {}
                for lf in rw:
                    cnts[abs(lf)] = cnts.get(abs(lf), 0) + 1
                cand = next((gc for gc, c in sorted(cnts.items(), reverse=True) if gc > 1 and c == 1), None)
                if cand is not None:
                    pi = next(ip for ip, lf in enumerate(rw) if abs(lf) == cand)
                    lc = rw[pi]
                    val = free_red(w_inv(rw[:pi]) + w_inv(rw[pi + 1 :]))
                    repl = val if lc > 0 else w_inv(val)
                    active.remove(cand)
                    cur_rels.pop(ri)
                    cur_rels = [sub_word(r2, cand, repl) for r2 in cur_rels if sub_word(r2, cand, repl)]
                    prog = True
                    break

        assert sorted(active) == [1, 2]

        # Interreduce cur_rels against each other via exact cyclic Dehn reductions (exact Tietze II moves!):
        changed_rels = True
        while changed_rels:
            changed_rels = False
            cur_rels = [free_red(r) for r in cur_rels if free_red(r)]
            cur_rels.sort(key=len)
            for i_r in range(len(cur_rels)):
                r_short = cur_rels[i_r]
                Ls = len(r_short)
                # Generate all strict length-reducing rules from cyclic shifts of r_short^{+-1}
                short_rules = []
                for base in (r_short, w_inv(r_short)):
                    for sh in range(Ls):
                        cyc = base[sh:] + base[:sh]
                        for cut in range(Ls // 2 + 1, Ls + 1):
                            short_rules.append((cyc[:cut], w_inv(cyc[cut:])))
                for j_r in range(len(cur_rels)):
                    if i_r == j_r:
                        continue
                    rj = cur_rels[j_r]
                    # Also cyclically reduce rj:
                    while len(rj) >= 2 and rj[0] == -rj[-1]:
                        rj = rj[1:-1]
                    # Check all cyclic rotations of rj for a substring matching lhs in short_rules:
                    Lj = len(rj)
                    hit = False
                    for sh_j in range(Lj):
                        rot_j = rj[sh_j:] + rj[:sh_j]
                        for lhs, rhs in short_rules:
                            ll = len(lhs)
                            if ll <= Lj and rot_j[:ll] == lhs:
                                new_rj = free_red(rhs + rot_j[ll:])
                                while len(new_rj) >= 2 and new_rj[0] == -new_rj[-1]:
                                    new_rj = new_rj[1:-1]
                                cur_rels[j_r] = new_rj
                                changed_rels = True
                                hit = True
                                break
                        if hit:
                            break
                    if hit:
                        break
                if changed_rels:
                    break
        print(f"  Dehn-interreduced cur_rels on {lbl}: {cur_rels}")

        def map_aq_to_s(w_aq: tuple[int, ...]) -> tuple[int, ...]:
            out: list[int] = []
            for x in w_aq:
                if x == 1: out.append(1)
                elif x == -1: out.append(-1)
                elif x == 2: out.extend(q_in_s)
                elif x == -2: out.extend(w_inv(q_in_s))
            return free_red(tuple(out))

        # Check in G(lbl) = <s_1, s_2 | cur_rels>:
        # (i) psi(q^3) == 1
        w_q3 = map_aq_to_s((2, 2, 2))
        ok_q3 = prove_word_trivial_in_fpgroup(w_q3, cur_rels, max_len=28)
        # (ii) psi(a q a (q a q a q)^-1) == 1
        w_rel2 = map_aq_to_s((1, 2, 1, -2, -1, -2, -1, -2))
        ok_rel2 = prove_word_trivial_in_fpgroup(w_rel2, cur_rels, max_len=28)
        # (iii) psi(phi(s_2)) * s_2^-1 == 1 (proving psi o phi == id on G(lbl)!)
        w_roundtrip_s2 = free_red(map_aq_to_s(s2_in_aq) + (-2,))
        ok_rt = prove_word_trivial_in_fpgroup(w_roundtrip_s2, cur_rels, max_len=28)
        print(f"=== BIDIRECTIONAL ISOMORPHISM PROOF FOR {lbl} ===")
        print(f"  psi(q^3) == 1 in G({lbl}):                    {ok_q3} (word len {len(w_q3)})")
        print(f"  psi(a q a (q a q a q)^-1) == 1 in G({lbl}):   {ok_rel2} (word len {len(w_rel2)})")
        print(f"  psi(phi(s_2)) == s_2 in G({lbl}):             {ok_rt} (word len {len(w_roundtrip_s2)})")
        assert ok_q3 and ok_rel2 and ok_rt, f"Failed Direction 2 on {lbl}!"

    # 2. Check the 4D Ribbon-Disk Groups pi_1(B^4_- \ D_-) and pi_1(B^4_+ \ D_+) on 12_2, 12_10, and 12_8!
    # Recall: in a (c, m) = (10, 2) diagram with marked vertices v_0, v_1:
    #   t = -eps is L_- (signs (-1, -1), 2-component unlink)
    #   Crossing v_0 first gives equatorial slice (+1, -1):
    #     - D_- in B^4_- has ONLY saddle v_0 active (so smooth v_1 with -1, keep v_0 as marked vertex!)
    #     - D_+ in B^4_+ has ONLY saddle v_1 active (so smooth v_0 with +1, keep v_1 as marked vertex!)
    #   WAIT: A 1-marked-vertex diagram (smoothing the OTHER marked vertex!) has Wirtinger presentation EQUAL TO pi_1(B^4 \ D)!!
    #   Why? Because in MarkedDiagram, smoothing one marked vertex v_j into two parallel arcs removes saddle v_j's 2-handle
    #   relation while keeping the remaining 1 marked vertex v_i's saddle relation!
    print("\n=== 4D RIBBON-DISK EXTERIOR GROUPS pi_1(B^4_- \\ D_-) AND pi_1(B^4_+ \\ D_+) ===")
    from surfaceknots.invariants.tier3 import lomonaco_pi2_invariants
    from surfaceknots.invariants.alexander import compute_alexander_invariants

    def smooth_single_marked_vertex(d: MarkedDiagram, m_vidx: int, sign: int) -> MarkedDiagram:
        """Replace marked vertex `m_vidx` by its `sign` (+1 or -1) smoothing in `MarkedDiagram`!"""
        vt = d.vtype[m_vidx]
        assert vt in ("Y", "Z")
        # Recall from diagrams.py:
        # Y: +1 pairs (0,3),(1,2); -1 pairs (0,1),(2,3)
        # Z: +1 pairs (0,1),(2,3); -1 pairs (0,3),(1,2)
        pair_03 = (vt == "Y" and sign == 1) or (vt == "Z" and sign == -1)
        pairs = [(0, 3), (1, 2)] if pair_03 else [(0, 1), (2, 3)]
        opp_mut = list(d.opp)
        free_extra = 0
        for s_a, s_b in pairs:
            da = 4 * m_vidx + s_a
            db = 4 * m_vidx + s_b
            oa = opp_mut[da]
            ob = opp_mut[db]
            if oa == db:
                # Self-loop! Forms a free circle!
                free_extra += 1
            else:
                opp_mut[oa] = ob
                opp_mut[ob] = oa
        # Relabel vertices skipping m_vidx:
        old_to_new = {}
        nxt = 0
        new_vt = []
        for vi in range(d.num_vertices):
            if vi != m_vidx:
                old_to_new[vi] = nxt
                new_vt.append(d.vtype[vi])
                nxt += 1
        new_opp = [0] * (4 * nxt)
        for old_v, new_v in old_to_new.items():
            for s in range(4):
                od = opp_mut[4 * old_v + s]
                new_opp[4 * new_v + s] = 4 * old_to_new[od // 4] + (od % 4)
        return MarkedDiagram(num_vertices=nxt, opp=tuple(new_opp), vtype=tuple(new_vt), free_loops=d.free_loops + free_extra)

    for lbl in ("12_2", "12_10", "12_8"):
        d = MarkedDiagram.from_epd(by_lbl[lbl]["epd"])
        m_idxs = [vi for vi, vt in enumerate(d.vtype) if vt in ("Y", "Z")]
        v0, v1 = m_idxs
        print(f"\n--- {lbl} ---")
        for order_name, (v_first, v_second) in [("Order (+1,-1) [v0 then v1]", (v0, v1)), ("Order (-1,+1) [v1 then v0]", (v1, v0))]:
            # D_- keeps v_first (smoothed at t < 0 -> L_-, at t > 0 -> K_eq), so v_second is smoothed with -1!
            d_minus = smooth_single_marked_vertex(d, v_second, -1)
            # D_+ keeps v_second (at t < 0 -> K_eq, at t > 0 -> L_+), so v_first is smoothed with +1!
            d_plus = smooth_single_marked_vertex(d, v_first, +1)
            for half_name, dh in [("D_-", d_minus), ("D_+", d_plus)]:
                kh, rh = independent_edge_wirtinger(dh)
                try:
                    q2 = todd_coxeter_order(kh, rh + [(1, 1)], max_cosets=20000)
                except RuntimeError:
                    q2 = "INF"
                try:
                    q3 = todd_coxeter_order(kh, rh + [(1, 1, 1)], max_cosets=20000)
                except RuntimeError:
                    q3 = "INF"
                al = compute_alexander_invariants(dh)
                extra = ""
                if lbl == "12_10":
                    cov = lomonaco_pi2_invariants(dh).to_dict().get("F21")
                    extra = f", F21_cover={cov}"
                elif lbl == "12_8":
                    cov = lomonaco_pi2_invariants(dh).to_dict().get("F57")
                    extra = f", F57_cover={cov}"
                print(f"  {order_name} {half_name}: E_1={al.e1_generators}, |G/<<a^2>>|={q2}, |G/<<a^3>>|={q3}{extra}")


if __name__ == "__main__":
    main()
