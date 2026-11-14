"""Fundamental group pi_1(S^4 \\ F), Tietze simplification, H_1, and finite-group quotient counts.

Derives exact Wirtinger + saddle presentations directly on `MarkedDiagram` from either resolution
(`L_-(D)` or `L_+(D)`), handling orientable, non-orientable, and mixed multi-component surface-links.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import sympy as sp

from surfaceknots.diagrams import MarkedDiagram


# Words in free group on generators 0..g-1 are tuples of nonzero signed 1-based integers:
# `+k` represents `x_{k-1}`, `-k` represents `x_{k-1}^{-1}`.
Word = tuple[int, ...]


def reduce_word(w: Word) -> Word:
    """Freely reduce a word `w` by cancelling adjacent inverse pairs."""
    st: list[int] = []
    for lit in w:
        if st and st[-1] == -lit:
            st.pop()
        else:
            st.append(lit)
    return tuple(st)


def cyclically_reduce_word(w: Word) -> Word:
    """Freely and cyclically reduce a relator word `w`."""
    rw = list(reduce_word(w))
    while len(rw) >= 2 and rw[0] == -rw[-1]:
        rw.pop()
        rw.pop(0)
    return tuple(rw)


def invert_word(w: Word) -> Word:
    return tuple(-x for x in reversed(w))


@dataclass(frozen=True)
class GroupPresentation:
    """Finite presentation `<x_0, ..., x_{g-1} | r_0, ..., r_{m-1}>` of `pi_1(S^4 \\ F)`.

    `component_meridians[i]` gives a word representing a meridian of connected surface component `F_i`.
    """

    num_generators: int
    relators: tuple[Word, ...]
    component_meridians: tuple[Word, ...]

    def to_human_string(self) -> str:
        if self.num_generators == 0:
            return "< | >"
        names = [f"x{i+1}" for i in range(self.num_generators)]

        def fmt_word(w: Word) -> str:
            if not w:
                return "1"
            out = []
            i = 0
            while i < len(w):
                lit = w[i]
                exp = 1 if lit > 0 else -1
                g_idx = abs(lit) - 1
                while i + 1 < len(w) and abs(w[i + 1]) - 1 == g_idx:
                    exp += 1 if w[i + 1] > 0 else -1
                    i += 1
                if exp == 1:
                    out.append(names[g_idx])
                elif exp != 0:
                    out.append(f"{names[g_idx]}^{exp}")
                i += 1
            return " ".join(out) if out else "1"

        gens_str = ", ".join(names)
        if not self.relators:
            return f"<{gens_str} | >"
        rels_str = ", ".join(f"{fmt_word(r)} = 1" for r in self.relators)
        return f"<{gens_str} | {rels_str}>"

    def abelianization_invariant_factors(self) -> tuple[int, tuple[int, ...]]:
        """Return `(free_rank, torsion_invariants)` of `H_1(S^4 \\ F; Z)` via Smith Normal Form."""
        g = self.num_generators
        if g == 0:
            return (0, ())
        if not self.relators:
            return (g, ())
        rows = []
        for r in self.relators:
            row = [0] * g
            for lit in r:
                idx = abs(lit) - 1
                row[idx] += 1 if lit > 0 else -1
            if any(row):
                rows.append(row)
        if not rows:
            return (g, ())
        mat = sp.Matrix(rows)
        # Compute Smith normal form diagonal invariant factors over ZZ
        from sympy.matrices.normalforms import smith_normal_form

        snf = smith_normal_form(mat, domain=sp.ZZ)
        diag = [abs(int(snf[i, i])) for i in range(min(snf.rows, snf.cols))]
        nonzero = [d for d in diag if d > 0]
        free_rank = g - len(nonzero)
        torsion = tuple(sorted(d for d in nonzero if d > 1))
        return (free_rank, torsion)


def wirtinger_presentation(
    d: MarkedDiagram,
    from_resolution: Literal[-1, 1] = -1,
    simplify: bool = True,
) -> GroupPresentation:
    """Compute exact presentation of `pi_1(S^4 \\ F(D))` from `L_-(D)` (`-1`) or `L_+(D)` (`+1`)."""
    n = d.num_vertices
    dart_scomp, total_scomps = d._dart_surface_components()

    # Build internal smoothings on `0..4n-1` for chosen resolution:
    # At classical `v`: straight across `s ^ 2`, breaking over-arc iff under-strand!
    # At marked `v`: smooth along resolution without breaking over-arc!
    is_classical = [vt in ("X0", "X1") for vt in d.vtype]

    def is_under_slot(v: int, s: int) -> bool:
        vt = d.vtype[v]
        if vt == "X0":
            return (s & 1) == 0
        return (s & 1) == 1

    def marked_partner(v: int, s: int) -> int:
        vt = d.vtype[v]
        pair_03_12 = (vt == "Y" and from_resolution == 1) or (vt == "Z" and from_resolution == -1)
        if pair_03_12:
            return 3 - s
        return s ^ 1

    # 1. Orient every closed resolution cycle on `0..4n-1`:
    # `dart_dir[hd] = +1` if cycle ENTERS `hd>>2` at `hd`, `-1` if cycle EXITS `hd>>2` at `hd`.
    dart_dir = [0] * (4 * n)
    cycles: list[list[int]] = []  # list of entry darts in traversal order
    for start_d in range(4 * n):
        if dart_dir[start_d] != 0:
            continue
        cyc_entries: list[int] = []
        cur_in = start_d
        while dart_dir[cur_in] == 0:
            dart_dir[cur_in] = 1
            cyc_entries.append(cur_in)
            v, s = cur_in >> 2, cur_in & 3
            exit_s = (s ^ 2) if is_classical[v] else marked_partner(v, s)
            cur_out = (v << 2) | exit_s
            dart_dir[cur_out] = -1
            cur_in = d.opp[cur_out]
        cycles.append(cyc_entries)

    # 2. Partition each oriented cycle into unbroken OVER-ARCS (`0..g-1`):
    dart_arc = [-1] * (4 * n)
    next_gen = 0
    for cyc in cycles:
        # Find all indices `k` in `cyc` where `cyc[k]` enters an UNDER-crossing!
        under_positions = [
            k for k, hd in enumerate(cyc) if is_classical[hd >> 2] and is_under_slot(hd >> 2, hd & 3)
        ]
        if not under_positions:
            # Entire cycle never passes under any crossing -> 1 single generator!
            gid = next_gen
            next_gen += 1
            for hd in cyc:
                v, s = hd >> 2, hd & 3
                exit_s = (s ^ 2) if is_classical[v] else marked_partner(v, s)
                dart_arc[hd] = gid
                dart_arc[(v << 2) | exit_s] = gid
        else:
            m_u = len(under_positions)
            for idx_u in range(m_u):
                gid = next_gen
                next_gen += 1
                pos_start = under_positions[idx_u]
                pos_end = under_positions[(idx_u + 1) % m_u]
                # Over-arc starts at EXIT of `cyc[pos_start]` and ends at ENTRY of `cyc[pos_end]`!
                hd_start = cyc[pos_start]
                vs, ss = hd_start >> 2, hd_start & 3
                dart_arc[(vs << 2) | (ss ^ 2)] = gid

                k = (pos_start + 1) % len(cyc)
                while k != pos_end:
                    hd = cyc[k]
                    vv, sv = hd >> 2, hd & 3
                    ex_s = (sv ^ 2) if is_classical[vv] else marked_partner(vv, sv)
                    dart_arc[hd] = gid
                    dart_arc[(vv << 2) | ex_s] = gid
                    k = (k + 1) % len(cyc)

                hd_end = cyc[pos_end]
                dart_arc[hd_end] = gid

    # Free loops (`0` vertices on them): each contributes 1 generator
    free_loop_gens: list[int] = []
    for _ in range(d.free_loops):
        free_loop_gens.append(next_gen)
        next_gen += 1

    # Record one meridian generator per connected surface component `0..total_scomps-1`
    comp_meridian: list[Word | None] = [None] * total_scomps
    for hd in range(4 * n):
        sc = dart_scomp[hd]
        if comp_meridian[sc] is None:
            comp_meridian[sc] = (dart_arc[hd] + 1,)
    graph_scomps = total_scomps - d.free_loops
    for i, fg in enumerate(free_loop_gens):
        comp_meridian[graph_scomps + i] = (fg + 1,)

    # 3. Collect Wirtinger classical crossing relators and marked saddle relators:
    relators: list[Word] = []
    for v in range(n):
        vt = d.vtype[v]
        if vt in ("X0", "X1"):
            under_slots = (0, 2) if vt == "X0" else (1, 3)
            over_slots = (1, 3) if vt == "X0" else (0, 2)
            u_in = under_slots[0] if dart_dir[(v << 2) | under_slots[0]] == 1 else under_slots[1]
            u_out = u_in ^ 2
            o_in = over_slots[0] if dart_dir[(v << 2) | over_slots[0]] == 1 else over_slots[1]

            x = dart_arc[(v << 2) | u_in] + 1
            y = dart_arc[(v << 2) | u_out] + 1
            z = dart_arc[(v << 2) | o_in] + 1
            diff = (o_in - u_in) & 3
            if diff == 3:
                # Positive crossing: z x z^{-1} y^{-1} = 1
                w = cyclically_reduce_word((z, x, -z, -y))
            else:
                # Negative crossing: z^{-1} x z y^{-1} = 1
                w = cyclically_reduce_word((-z, x, z, -y))
            if w:
                relators.append(w)
        else:
            # Marked saddle vertex!
            # Check which pair `(s0, s0+1)` & `(s0+2, s0+3)` is spliced in `from_resolution`:
            # Pair (0,1)&(2,3) -> s0 = 0; Pair (1,2)&(3,0) -> s0 = 1!
            pair_03_12 = (vt == "Y" and from_resolution == 1) or (vt == "Z" and from_resolution == -1)
            s0 = 1 if pair_03_12 else 0
            d0 = (v << 2) | s0
            d2 = (v << 2) | ((s0 + 2) & 3)
            a1 = dart_arc[d0] + 1
            a2 = dart_arc[d2] + 1
            lit1 = a1 if dart_dir[d0] == 1 else -a1
            lit2 = -a2 if dart_dir[d2] == 1 else a2
            w = cyclically_reduce_word((lit1, lit2))
            if w:
                relators.append(w)

    pres = GroupPresentation(
        num_generators=next_gen,
        relators=tuple(relators),
        component_meridians=tuple(m if m is not None else () for m in comp_meridian),
    )
    return tietze_simplify(pres) if simplify else pres


def tietze_simplify(pres: GroupPresentation) -> GroupPresentation:
    """Simplify a group presentation via Tietze generator eliminations and cyclic reductions."""
    num_g = pres.num_generators
    rels = [cyclically_reduce_word(r) for r in pres.relators if cyclically_reduce_word(r)]
    meridians = [reduce_word(m) for m in pres.component_meridians]

    def sub_word(w: Word, target_gen_1based: int, repl: Word) -> Word:
        inv_repl = invert_word(repl)
        expanded: list[int] = []
        for lit in w:
            if lit == target_gen_1based:
                expanded.extend(repl)
            elif lit == -target_gen_1based:
                expanded.extend(inv_repl)
            else:
                expanded.append(lit)
        return reduce_word(tuple(expanded))

    changed = True
    while changed and num_g > 0:
        changed = False
        # Deduplicate relators up to inversion and cyclic shift
        uniq_rels: list[Word] = []
        seen_rel_keys = set()
        for r in rels:
            cr = cyclically_reduce_word(r)
            if not cr:
                continue
            inv_cr = invert_word(cr)
            # Canonical cyclic shift
            rots = [cr[i:] + cr[:i] for i in range(len(cr))] + [
                inv_cr[i:] + inv_cr[:i] for i in range(len(inv_cr))
            ]
            key = min(rots)
            if key not in seen_rel_keys:
                seen_rel_keys.add(key)
                uniq_rels.append(key)
        rels = uniq_rels

        # Look for any relator where some generator `g` (`1..num_g`) appears EXACTLY ONCE!
        # Sort candidate relators by length so short `x = y` or `x = z y z^{-1}` eliminate first!
        rels.sort(key=len)
        elim_target = None
        elim_repl: Word = ()
        elim_rel_idx = -1

        for r_idx, r in enumerate(rels):
            counts: dict[int, int] = {}
            for lit in r:
                g_abs = abs(lit)
                counts[g_abs] = counts.get(g_abs, 0) + 1
            singletons = [g_abs for g_abs, cnt in counts.items() if cnt == 1]
            if singletons:
                # Pick highest-numbered singleton generator to eliminate
                g_elim = max(singletons)
                pos = next(i for i, lit in enumerate(r) if abs(lit) == g_elim)
                lit = r[pos]
                # r[pos] * rest == 1, where rest = r[pos+1:] + r[:pos]
                rest = r[pos + 1 :] + r[:pos]
                if lit > 0:
                    # g * rest = 1 => g = rest^{-1}
                    elim_repl = invert_word(rest)
                else:
                    # g^{-1} * rest = 1 => g = rest
                    elim_repl = reduce_word(rest)
                elim_target = g_elim
                elim_rel_idx = r_idx
                break

        if elim_target is not None:
            changed = True
            new_rels: list[Word] = []
            for idx, r in enumerate(rels):
                if idx == elim_rel_idx:
                    continue
                sr = cyclically_reduce_word(sub_word(r, elim_target, elim_repl))
                if sr:
                    new_rels.append(sr)
            new_meridians = [sub_word(m, elim_target, elim_repl) for m in meridians]

            # Renumber generators > elim_target down by 1
            def shift_down(w: Word) -> Word:
                out = []
                for lit in w:
                    g_abs = abs(lit)
                    assert g_abs != elim_target
                    ng = g_abs - 1 if g_abs > elim_target else g_abs
                    out.append(ng if lit > 0 else -ng)
                return tuple(out)

            rels = [shift_down(r) for r in new_rels]
            meridians = [shift_down(m) for m in new_meridians]
            num_g -= 1

    return GroupPresentation(
        num_generators=num_g,
        relators=tuple(rels),
        component_meridians=tuple(meridians),
    )


# ---------------------------------------------------------------------------
# Fast Finite-Group Homomorphism Counter |Hom(pi_1, G)|
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class FiniteGroup:
    """Multiplication & inverse table representation of a finite group of order `n` (`0` = identity)."""

    name: str
    order: int
    mul: tuple[tuple[int, ...], ...]
    inv: tuple[int, ...]


def _group_from_permutations(name: str, gens: list[tuple[int, ...]]) -> FiniteGroup:
    """Generate a finite permutation group from generating permutations on `0..d-1`."""
    if not gens:
        return FiniteGroup(name=name, order=1, mul=((0,),), inv=(0,))
    deg = len(gens[0])
    ident = tuple(range(deg))
    elems = [ident]
    seen = {ident: 0}
    idx = 0
    while idx < len(elems):
        p = elems[idx]
        idx += 1
        for g in gens:
            comp = tuple(p[g[i]] for i in range(deg))
            if comp not in seen:
                seen[comp] = len(elems)
                elems.append(comp)
    m = len(elems)
    mul = [[0] * m for _ in range(m)]
    inv = [0] * m
    for i, pi in enumerate(elems):
        inv_p = [0] * deg
        for a, b in enumerate(pi):
            inv_p[b] = a
        inv[i] = seen[tuple(inv_p)]
        for j, pj in enumerate(elems):
            comp = tuple(pi[pj[k]] for k in range(deg))
            mul[i][j] = seen[comp]
    return FiniteGroup(name=name, order=m, mul=tuple(tuple(r) for r in mul), inv=tuple(inv))


def _dicyclic_group(name: str, n: int) -> FiniteGroup:
    """Dicyclic group `Dic_{4n} = <a, x | a^{2n} = 1, x^2 = a^n, x a x^{-1} = a^{-1}>` of order `4n`."""
    # Represent via regular left-multiplication permutations on elements `(e in 0..1, k in 0..2n-1)`:
    # `(0, k)` = `a^k`, `(1, k)` = `x * a^k`.
    order = 4 * n

    def enc(e: int, k: int) -> int:
        return e * (2 * n) + (k % (2 * n))

    # Left multiplication by `a`: `a * a^k = a^{k+1}`; `a * (x a^k) = x a^{k-1}`!
    perm_a = [0] * order
    # Left multiplication by `x`: `x * a^k = x a^k`; `x * (x a^k) = a^{n + k}`!
    perm_x = [0] * order
    for k in range(2 * n):
        perm_a[enc(0, k)] = enc(0, k + 1)
        perm_a[enc(1, k)] = enc(1, k - 1)
        perm_x[enc(0, k)] = enc(1, k)
        perm_x[enc(1, k)] = enc(0, n + k)
    return _group_from_permutations(name, [tuple(perm_a), tuple(perm_x)])


STANDARD_TEST_GROUPS: tuple[FiniteGroup, ...] = (
    # S_3 (order 6)
    _group_from_permutations("S3", [(1, 2, 0), (1, 0, 2)]),
    # D_8 (order 8)
    _group_from_permutations("D8", [(1, 2, 3, 0), (0, 3, 2, 1)]),
    # Q_8 (quaternion group of order 8 = Dic_8)
    _dicyclic_group("Q8", 2),
    # D_10 (order 10)
    _group_from_permutations("D10", [(1, 2, 3, 4, 0), (0, 4, 3, 2, 1)]),
    # A_4 (order 12)
    _group_from_permutations("A4", [(1, 2, 0, 3), (0, 2, 3, 1)]),
    # Dic_12 (order 12)
    _dicyclic_group("Dic12", 3),
    # Dic_16 (order 16 = <2,2,4>, detects 10_1^{-1,-1}!)
    _dicyclic_group("Dic16", 4),
    # F_21 = Z_7 rtimes_2 Z_3 (order 21, detects Baumslag-Solitar BS(1,2) = pi_1(9_1)!)
    _group_from_permutations(
        "F21",
        [
            (1, 2, 3, 4, 5, 6, 0),
            tuple((2 * i) % 7 for i in range(7)),
        ],
    ),
    # S_4 (order 24)
    _group_from_permutations("S4", [(1, 2, 3, 0), (1, 0, 2, 3)]),
    # A_5 (order 60)
    _group_from_permutations("A5", [(1, 2, 3, 4, 0), (1, 2, 0, 3, 4)]),
)


_CONJ_CLASSES_CACHE: dict[str, list[tuple[int, int]]] = {}


def _get_conjugacy_classes(group: FiniteGroup) -> list[tuple[int, int]]:
    cached = _CONJ_CLASSES_CACHE.get(group.name)
    if cached is not None:
        return cached
    order = group.order
    mul = group.mul
    inv = group.inv
    seen = [False] * order
    classes: list[tuple[int, int]] = []
    for x in range(order):
        if seen[x]:
            continue
        cls = set()
        for h in range(order):
            c = mul[mul[h][x]][inv[h]]
            cls.add(c)
            seen[c] = True
        classes.append((x, len(cls)))
    _CONJ_CLASSES_CACHE[group.name] = classes
    return classes


def count_homomorphisms(pres: GroupPresentation, group: FiniteGroup) -> int:
    """Exact count `|Hom(pi_1, G)|` via connected-component factoring, conjugacy symmetry, and unit solving."""
    g = pres.num_generators
    if g == 0:
        return 1
    order = group.order
    mul = group.mul
    inv = group.inv

    non_empty_rels = [r for r in pres.relators if r]
    if not non_empty_rels:
        return order**g

    # Partition generators `0..g-1` into connected components of the relator hypergraph:
    adj: list[set[int]] = [set() for _ in range(g)]
    gen_rels: list[list[int]] = [[] for _ in range(g)]
    rel_var_sets: list[set[int]] = []
    for ri, r in enumerate(non_empty_rels):
        vs = {abs(lit) - 1 for lit in r}
        rel_var_sets.append(vs)
        for v in vs:
            gen_rels[v].append(ri)
            adj[v].update(vs)

    comp_vis = [False] * g
    conj_classes = _get_conjugacy_classes(group)
    grand_total = 1

    for start_v in range(g):
        if comp_vis[start_v]:
            continue
        if not gen_rels[start_v]:
            comp_vis[start_v] = True
            grand_total *= order
            continue

        # BFS connected component of generators and collect its relators:
        comp_vars: list[int] = []
        comp_rel_ids: set[int] = set()
        q = [start_v]
        comp_vis[start_v] = True
        while q:
            u = q.pop()
            comp_vars.append(u)
            comp_rel_ids.update(gen_rels[u])
            for w in adj[u]:
                if not comp_vis[w]:
                    comp_vis[w] = True
                    q.append(w)

        # Greedy variable ordering: pick next variable completing the most relators (then highest degree):
        k = len(comp_vars)
        placed: set[int] = set()
        ordered_vars: list[int] = []
        rem_vars = set(comp_vars)
        c_rels = [non_empty_rels[ri] for ri in sorted(comp_rel_ids)]
        c_vsets = [rel_var_sets[ri] for ri in sorted(comp_rel_ids)]
        while rem_vars:
            best_v = max(
                rem_vars,
                key=lambda cand_v: (
                    sum(1 for vs in c_vsets if cand_v in vs and vs.issubset(placed | {cand_v})),
                    sum(1 for vs in c_vsets if cand_v in vs and bool(vs & placed)),
                    len(gen_rels[cand_v]),
                    -cand_v,
                ),
            )
            ordered_vars.append(best_v)
            placed.add(best_v)
            rem_vars.remove(best_v)

        old_to_new = {ov: ni for ni, ov in enumerate(ordered_vars)}
        mapped_rels: list[Word] = [
            tuple((old_to_new[lit - 1] + 1) if lit > 0 else -(old_to_new[-lit - 1] + 1) for lit in r)
            for r in c_rels
        ]

        rels_by_max_gen: list[list[Word]] = [[] for _ in range(k)]
        unique_solve: list[tuple[Word, int, bool] | None] = [None] * k
        for r in mapped_rels:
            max_idx = max(abs(lit) - 1 for lit in r)
            rels_by_max_gen[max_idx].append(r)

        for idx in range(k):
            for r in rels_by_max_gen[idx]:
                positions = [p for p, lit in enumerate(r) if abs(lit) - 1 == idx]
                if len(positions) == 1:
                    pos = positions[0]
                    unique_solve[idx] = (r, pos, r[pos] > 0)
                    break

        assign = [0] * k

        def _check_rels(idx: int) -> bool:
            for r in rels_by_max_gen[idx]:
                acc = 0
                for lit in r:
                    gv = assign[lit - 1] if lit > 0 else inv[assign[-lit - 1]]
                    acc = mul[acc][gv]
                if acc != 0:
                    return False
            return True

        def backtrack(idx: int) -> int:
            if idx == k:
                return 1
            u_info = unique_solve[idx]
            if u_info is not None:
                r_u, pos, is_pos = u_info
                pref = 0
                for lit in r_u[:pos]:
                    pref = mul[pref][assign[lit - 1] if lit > 0 else inv[assign[-lit - 1]]]
                suff = 0
                for lit in r_u[pos + 1 :]:
                    suff = mul[suff][assign[lit - 1] if lit > 0 else inv[assign[-lit - 1]]]
                val_lit = mul[inv[pref]][inv[suff]]
                assign[idx] = val_lit if is_pos else inv[val_lit]
                if _check_rels(idx):
                    return backtrack(idx + 1)
                return 0

            if idx == 0:
                tot0 = 0
                for rep_val, cls_sz in conj_classes:
                    assign[0] = rep_val
                    if _check_rels(0):
                        sub = backtrack(1)
                        if sub:
                            tot0 += cls_sz * sub
                return tot0

            tot = 0
            for val in range(order):
                assign[idx] = val
                if _check_rels(idx):
                    tot += backtrack(idx + 1)
            return tot

        comp_cnt = backtrack(0)
        if comp_cnt == 0:
            return 0
        grand_total *= comp_cnt

    return grand_total



_STANDARD_GROUPS_BY_NAME: dict[str, FiniteGroup] = {grp.name: grp for grp in STANDARD_TEST_GROUPS}


def finite_quotient_profile(
    pres: GroupPresentation,
    groups: tuple[FiniteGroup | str, ...] = STANDARD_TEST_GROUPS,
) -> dict[str, int]:
    """Compute `{group.name: |Hom(pi_1, G)|}` across standard finite target groups."""
    resolved = [_STANDARD_GROUPS_BY_NAME[g] if isinstance(g, str) else g for g in groups]
    return {grp.name: count_homomorphisms(pres, grp) for grp in resolved}

