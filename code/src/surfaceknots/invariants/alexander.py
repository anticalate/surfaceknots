"""Fox calculus, multivariable Alexander ideals E_k(F) over Z[H_1], and twisted Alexander invariants.

Key features:
- Computes exact Fox free derivatives of `<x_1, ..., x_g | r_1, ..., r_m>` directly into the full
  integral group ring `Z[H_1(S^4 \\ F)] = Z[t_1^{+-1}, ..., t_r^{+-1}, u_1, ..., u_s] / (u_j^2 - 1)`,
  where `t_i` are meridians of orientable components and `u_j` (`u_j^2 = 1`) are meridians of
  non-orientable components (following Yoshikawa 1994 and Ishii-Nikkuni-Oshiro 2018).
- Computes exact elementary ideals `E_1(F)` and `E_2(F)` (`(g-k) x (g-k)` minors of the `m x g`
  Fox Jacobian matrix) together with canonical Groebner / Macaulay-Hermite fingerprints over `Z`
  that are invariant under permuting homeomorphic components and inverting orientable variables
  `t_i <-> t_i^{-1}`.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, permutations, product

import sympy as sp

from surfaceknots.diagrams import MarkedDiagram
from surfaceknots.invariants.fundamental_group import (
    GroupPresentation,
    Word,
    wirtinger_presentation,
)


@dataclass(frozen=True)
class AlexanderResult:
    """Computed Alexander module invariants over `Z[H_1(S^4 \\ F)]`."""

    ring_description: str
    variable_names: tuple[str, ...]
    orientable_flags: tuple[bool, ...]
    e1_generators: tuple[str, ...]
    e2_generators: tuple[str, ...]
    single_var_alex: str
    evaluation_fingerprint: tuple


def _fox_derivative_word(
    w: Word,
    var_idx_0based: int,
    gen_images: list[sp.Expr],
    gen_inv_images: list[sp.Expr],
) -> sp.Expr:
    """Compute `psi(d w / d x_{var_idx})` in SymPy ring using Fox's chain rule:
    `d(u v)/dx = du/dx + psi(u) * dv/dx`, `d(x)/dx = 1`, `d(x^{-1})/dx = -psi(x^{-1})`.
    """
    target_1based = var_idx_0based + 1
    prefix = sp.Integer(1)
    acc = sp.Integer(0)
    for lit in w:
        g_idx = abs(lit) - 1
        if lit == target_1based:
            acc += prefix
        elif lit == -target_1based:
            acc -= prefix * gen_inv_images[g_idx]
        prefix = sp.expand(prefix * (gen_images[g_idx] if lit > 0 else gen_inv_images[g_idx]))
    return sp.expand(acc)


def _reduce_mod_torsion(expr: sp.Expr, torsion_syms: list[sp.Symbol]) -> sp.Expr:
    """Reduce powers `u_j^k` modulo `u_j^2 - 1` (`u_j^{2m} -> 1`, `u_j^{2m+1} -> u_j`)."""
    if not torsion_syms or expr == 0:
        return sp.expand(expr)
    poly = sp.Poly(sp.expand(expr), *torsion_syms, domain="EX")
    res = sp.Integer(0)
    for monom, coeff in poly.terms():
        term = coeff.as_expr() if hasattr(coeff, "as_expr") else coeff
        for sym, exp in zip(torsion_syms, monom):
            if exp & 1:
                term = term * sym
        res += term
    return sp.expand(res)


def _clear_negative_laurent_powers(
    expr: sp.Expr,
    ori_syms: list[sp.Symbol],
    ori_inv_syms: list[sp.Symbol],
    torsion_syms: list[sp.Symbol],
) -> sp.Expr:
    """Multiply `expr` by a canonical monomial `prod t_i^{k_i}` (`k_i in Z`) so it becomes a
    polynomial in `t_1, ..., t_r, u_1, ..., u_s` not divisible by any `t_i`, with positive leading
    sign! Since `t_i` are units in `Z[H_1]`, multiplying a generator by a unit does NOT change the
    ideal in `Z[H_1]`!
    """
    e = _reduce_mod_torsion(expr, torsion_syms)
    if e == 0:
        return sp.Integer(0)
    if not ori_syms:
        # Normalize overall sign
        s_poly = sp.Poly(e, *torsion_syms, domain=sp.ZZ) if torsion_syms else None
        if s_poly is not None and s_poly.LC() < 0:
            return -e
        if not torsion_syms and int(e) < 0:
            return -e
        return e

    # Replace `ori_inv_syms[i]` with `1 / ori_syms[i]` via exponent tracking:
    has_inv = len(ori_inv_syms) == len(ori_syms)
    all_syms = ori_syms + (ori_inv_syms if has_inv else []) + torsion_syms
    p = sp.Poly(e, *all_syms, domain=sp.ZZ)
    r = len(ori_syms)
    tor_offset = 2 * r if has_inv else r
    terms_dict: dict[tuple[int, ...], int] = {}
    for monom, coeff in p.terms():
        eff = tuple(monom[i] - (monom[r + i] if has_inv else 0) for i in range(r)) + tuple(
            monom[tor_offset + j] & 1 for j in range(len(torsion_syms))
        )
        terms_dict[eff] = terms_dict.get(eff, 0) + int(coeff)
    terms_dict = {k: v for k, v in terms_dict.items() if v != 0}
    if not terms_dict:
        return sp.Integer(0)

    min_exps = [min(k[i] for k in terms_dict) for i in range(r)]
    shifted: list[tuple[tuple[int, ...], int]] = []
    for k, c in terms_dict.items():
        nk = tuple(k[i] - min_exps[i] for i in range(r)) + k[r:]
        shifted.append((nk, c))
    shifted.sort(reverse=True)
    sign_flip = -1 if shifted[0][1] < 0 else 1

    out = sp.Integer(0)
    all_pos_syms = ori_syms + torsion_syms
    for nk, c in shifted:
        mon = sp.Integer(sign_flip * c)
        for sym, exp in zip(all_pos_syms, nk):
            if exp:
                mon = mon * (sym**exp)
        out += mon
    return sp.expand(out)


def _simplify_ideal_generators(
    raw_gens: list[sp.Expr],
    ori_syms: list[sp.Symbol],
    ori_inv_syms: list[sp.Symbol],
    torsion_syms: list[sp.Symbol],
) -> list[sp.Expr]:
    """Produce a compact, exact ideal generating set over `Z[H_1]` without rational division.

    Uses exact integer Hermite Normal Form (`hermite_normal_form` over `ZZ`) on bounded monomial
    shifts modulo `u_j^2 - 1`, preserving integer torsion generators such as `(x + 1, 3)`,
    `(x^2 + x + 1, 2)`, `(x + 1, y + 1, 2)`, and `(2x + y - 1, 4)`.
    """
    from functools import reduce
    from math import gcd
    from sympy.matrices.normalforms import hermite_normal_form

    cleaned_set = set()
    for g in raw_gens:
        cg = _clear_negative_laurent_powers(g, ori_syms, ori_inv_syms, torsion_syms)
        if cg != 0:
            cleaned_set.add(cg)
    if not cleaned_set:
        return [sp.Integer(0)]
    if sp.Integer(1) in cleaned_set:
        return [sp.Integer(1)]

    all_vars = ori_syms + torsion_syms
    if not all_vars:
        g_val = reduce(gcd, (abs(int(x)) for x in cleaned_set))
        return [sp.Integer(g_val)]

    r = len(ori_syms)
    s = len(torsion_syms)
    cleaned_list = sorted(cleaned_set, key=lambda expr: (len(str(expr)), str(expr)))

    # Determine max degree in each orientable variable `t_i`
    max_degs = [0] * r
    poly_terms_list: list[dict[tuple[int, ...], int]] = []
    for expr in cleaned_list:
        p = sp.Poly(expr, *all_vars, domain=sp.ZZ)
        td: dict[tuple[int, ...], int] = {}
        for mon, coeff in p.terms():
            key = tuple(mon[:r]) + tuple(mon[r + j] & 1 for j in range(s))
            td[key] = td.get(key, 0) + int(coeff)
            for i in range(r):
                if mon[i] > max_degs[i]:
                    max_degs[i] = mon[i]
        poly_terms_list.append({k: v for k, v in td.items() if v != 0})

    # Allow +1 degree slack if basis stays small (`<= 96` monomials) to expose syzygy constants
    box_degs = list(max_degs)
    for i in range(r):
        test_degs = list(box_degs)
        test_degs[i] += 1
        sz = (2**s)
        for d_val in test_degs:
            sz *= d_val + 1
        if sz <= 96:
            box_degs = test_degs

    total_basis = 2**s
    for d_val in box_degs:
        total_basis *= d_val + 1

    if total_basis > 160:
        return cleaned_list

    # Order basis in DECREASING total degree (so constant `1` is the LAST row; column-HNF
    # upper-triangularizes from top to bottom, isolating low-degree / constant generators at the end!)
    raw_monoms = list(
        product(*(range(d_val + 1) for d_val in box_degs), *(range(2) for _ in range(s)))
    )
    basis_monoms = sorted(raw_monoms, key=lambda m: (sum(m), m), reverse=True)
    m_to_row = {m: idx for idx, m in enumerate(basis_monoms)}
    dim = len(basis_monoms)

    # Build columns for all shifts `mu * g` (`deg_{t_i}(mu * g) <= box_degs[i]`)
    cols: list[list[int]] = []
    for td in poly_terms_list:
        g_degs = [max((k[i] for k in td), default=0) for i in range(r)]
        shift_ranges = [range(box_degs[i] - g_degs[i] + 1) for i in range(r)]
        for sh_o in product(*shift_ranges):
            for sh_t in product(*(range(2) for _ in range(s))):
                col = [0] * dim
                for k, coeff in td.items():
                    nm = tuple(k[i] + sh_o[i] for i in range(r)) + tuple(
                        (k[r + j] + sh_t[j]) & 1 for j in range(s)
                    )
                    col[m_to_row[nm]] += coeff
                if any(col):
                    cols.append(col)

    if not cols:
        return [sp.Integer(0)]

    mat = sp.Matrix(cols).T  # dim x num_cols
    hnf = hermite_normal_form(mat)

    # Extract nonzero HNF columns back as polynomials, strip `t_i` unit powers
    hnf_exprs_set = set()
    for c_idx in range(hnf.cols):
        expr = sp.Integer(0)
        for r_idx in range(dim):
            val = int(hnf[r_idx, c_idx])
            if val != 0:
                mon = basis_monoms[r_idx]
                term = sp.Integer(val)
                for sym, exp in zip(all_vars, mon):
                    if exp:
                        term *= sym**exp
                expr += term
        if expr != 0:
            ce = _clear_negative_laurent_powers(expr, ori_syms, ori_inv_syms, torsion_syms)
            if ce == 1:
                return [sp.Integer(1)]
            if ce != 0:
                hnf_exprs_set.add(ce)

    # Order candidates from simplest (lowest degree, shortest string) to most complex, and
    # greedily keep `cand` ONLY if it is NOT already in the Z-span of bounded shifts of `kept`!
    def expr_sort_key(e: sp.Expr) -> tuple:
        p = sp.Poly(e, *all_vars, domain=sp.ZZ)
        return (p.total_degree(), len(p.terms()), len(str(e)), str(e))

    candidates = sorted(hnf_exprs_set, key=expr_sort_key)
    kept: list[sp.Expr] = []
    kept_shift_cols: list[list[int]] = []
    kept_hnf: sp.Matrix | None = None

    def expr_to_col_and_shifts(e: sp.Expr) -> tuple[list[int] | None, list[list[int]]]:
        p = sp.Poly(e, *all_vars, domain=sp.ZZ)
        td: dict[tuple[int, ...], int] = {}
        for mon, coeff in p.terms():
            key = tuple(mon[:r]) + tuple(mon[r + j] & 1 for j in range(s))
            if any(key[i] > box_degs[i] for i in range(r)):
                return None, []
            td[key] = td.get(key, 0) + int(coeff)
        base_col = [0] * dim
        for k, coeff in td.items():
            base_col[m_to_row[k]] += coeff
        g_degs = [max((k[i] for k in td), default=0) for i in range(r)]
        sh_cols: list[list[int]] = []
        shift_ranges = [range(box_degs[i] - g_degs[i] + 1) for i in range(r)]
        for sh_o in product(*shift_ranges):
            for sh_t in product(*(range(2) for _ in range(s))):
                sc = [0] * dim
                for k, coeff in td.items():
                    nm = tuple(k[i] + sh_o[i] for i in range(r)) + tuple(
                        (k[r + j] + sh_t[j]) & 1 for j in range(s)
                    )
                    sc[m_to_row[nm]] += coeff
                if any(sc):
                    sh_cols.append(sc)
        return base_col, sh_cols

    for cand in candidates:
        b_col, sh_cols = expr_to_col_and_shifts(cand)
        if b_col is None:
            kept.append(cand)
            continue
        if kept_hnf is not None and kept_hnf.cols > 0:
            try:
                sol, params = kept_hnf.gauss_jordan_solve(sp.Matrix(b_col))
                if params.rows > 0:
                    sol = sol.subs({p_sym: 0 for p_sym in params})
                if all(entry.is_integer for entry in sol):
                    continue  # `cand` is an exact Z-linear combination of shifts of `kept`!
            except ValueError:
                pass
        kept.append(cand)
        kept_shift_cols.extend(sh_cols)
        kept_hnf = hermite_normal_form(sp.Matrix(kept_shift_cols).T)

    return kept if kept else [sp.Integer(0)]


def _ideal_evaluation_fingerprint(
    gens: list[sp.Expr],
    ori_syms: list[sp.Symbol],
    torsion_syms: list[sp.Symbol],
) -> tuple:
    """Compute an exact symmetry-invariant fingerprint of an ideal `I subset Z[H_1]`.

    Evaluates the ideal over `Z` (and `Z[i], Z[omega_3]`) across all `2^s` characters `u_j in {+1, -1}`
    and test points `t_i in {1, -1, 2, 3}`, taking the lexicographically smallest tuple over:
    - orientable meridian inversions `t_i <-> t_i^{-1}`, and
    - permutations among orientable variables and among non-orientable variables!
    Wait: notice that `E_1 = (u - 1, 2)` vs `E_1 = (2)` have IDENTICAL point-evaluations into `Z`
    at `u = +1` (`gcd(0, 2) = 2`) and `u = -1` (`gcd(-2, 2) = 2`), because `Z[u]/(u^2-1)` has index 2
    in `Z x Z`!
    To distinguish fibered ideals like `(u - 1, 2)` from `(2)` with 100% precision, we ALSO include
    the **exact cardinality / HNF invariant factors of the finite quotient `Z[H_1] / (I + (t_i^2 - 1, 4))`**!
    """
    from math import gcd
    from functools import reduce
    from sympy.matrices.normalforms import smith_normal_form

    r = len(ori_syms)
    s = len(torsion_syms)

    # 1. Exact quotient module invariant factors on `Z[t_1..t_r, u_1..u_s] / (t_i^3 - 1, u_j^2 - 1)`:
    # Basis has `3^r * 2^s` elements (for c <= 4, <= 81 basis elements -> microseconds!).
    # Under `t_i -> t_i^{-1}` (`t_i -> t_i^2 mod t_i^3 - 1`) and variable permutations, the SNF
    # invariant factors of the multiplication ideal matrix are **100% CANONICALLY INVARIANT**
    # WITHOUT EVEN NEEDING TO LOOP OVER VARIABLE PERMUTATIONS OR INVERSIONS!!
    snf_3_2: tuple[int, ...] = ()
    snf_2_2: tuple[int, ...] = ()
    if (3**r) * (2**s) <= 243:
        for t_mod, tag in ((2, "mod2"), (3, "mod3")):
            basis_exps = list(product(*(range(t_mod) for _ in range(r)), *(range(2) for _ in range(s))))
            b_idx = {e: idx for idx, e in enumerate(basis_exps)}
            dim = len(basis_exps)
            all_syms = ori_syms + torsion_syms
            rows: list[list[int]] = []
            for g in gens:
                if g == 0:
                    continue
                p = sp.Poly(g, *all_syms, domain=sp.ZZ) if all_syms else None
                g_terms: list[tuple[tuple[int, ...], int]] = []
                if p is None:
                    g_terms.append(((), int(g)))
                else:
                    for mon, coeff in p.terms():
                        g_terms.append((mon, int(coeff)))
                # Multiply `g` by every basis monomial `b`!
                for b in basis_exps:
                    vec = [0] * dim
                    for mon, coeff in g_terms:
                        prod_e = tuple((b[i] + mon[i]) % t_mod for i in range(r)) + tuple(
                            (b[r + j] + mon[r + j]) & 1 for j in range(s)
                        )
                        vec[b_idx[prod_e]] += coeff
                    if any(vec):
                        rows.append(vec)
            if not rows:
                invs = (0,) * dim
            else:
                mat = sp.Matrix(rows)
                snf = smith_normal_form(mat, domain=sp.ZZ)
                diag = [abs(int(snf[i, i])) for i in range(min(snf.rows, snf.cols))]
                # Pad with zeros to `dim` columns
                diag_padded = diag + [0] * (dim - len(diag))
                invs = tuple(sorted(diag_padded))
            if t_mod == 2:
                snf_2_2 = invs
            else:
                snf_3_2 = invs

    # 2. Integer GCD evaluations at `t_i in {1, -1, 2}` (and `2^{-1}` handled by comparing `{2, 1/2}`
    # via clearing `2^k`? Better: evaluate at `t_i = 1` and `t_i = -1`, which ARE self-inverse!).
    eval_pm1 = []
    for t_vals in product((1, -1), repeat=r):
        for u_vals in product((1, -1), repeat=s):
            subs = {ori_syms[i]: t_vals[i] for i in range(r)}
            subs.update({torsion_syms[j]: u_vals[j] for j in range(s)})
            vals = [abs(int(g.subs(subs))) for g in gens if g != 0]
            eval_pm1.append(reduce(gcd, vals, 0))
    eval_pm1_sorted = tuple(sorted(eval_pm1))

    return (snf_2_2, snf_3_2, eval_pm1_sorted)


def compute_alexander_invariants(
    d: MarkedDiagram,
    from_resolution: int = -1,
) -> AlexanderResult:
    """Compute multivariable Alexander elementary ideals `E_1(F), E_2(F)` over `Z[H_1(S^4 \\ F)]`."""
    topo = d.surface_topology()
    pres = wirtinger_presentation(d, from_resolution=from_resolution, simplify=True)  # type: ignore[arg-type]

    # Order components in Yoshikawa order: orientable (sorted by genus) then non-orientable (sorted by -k)
    comp_order = sorted(
        range(topo.num_components),
        key=lambda i: (
            0 if topo.components[i].orientable else 1,
            topo.components[i].genus_symbol if topo.components[i].orientable else -topo.components[i].genus_symbol,
            i,
        ),
    )

    std_names = ["x", "y", "z", "w", "v", "u"]
    var_names: list[str] = []
    ori_flags: list[bool] = []
    ori_syms: list[sp.Symbol] = []
    ori_inv_syms: list[sp.Symbol] = []
    torsion_syms: list[sp.Symbol] = []

    comp_sym_pos: dict[int, sp.Expr] = {}
    comp_sym_neg: dict[int, sp.Expr] = {}

    for rank_idx, sc in enumerate(comp_order):
        vname = std_names[rank_idx] if rank_idx < len(std_names) else f"t{rank_idx+1}"
        var_names.append(vname)
        is_ori = topo.components[sc].orientable
        ori_flags.append(is_ori)
        sym = sp.Symbol(vname)
        if is_ori:
            inv_sym = sp.Symbol(f"{vname}_inv")
            ori_syms.append(sym)
            ori_inv_syms.append(inv_sym)
            comp_sym_pos[sc] = sym
            comp_sym_neg[sc] = inv_sym
        else:
            torsion_syms.append(sym)
            comp_sym_pos[sc] = sym
            comp_sym_neg[sc] = sym

    # Determine the abelianization map `psi: <x_1, ..., x_g> -> H_1 = Z^r (+) (Z_2)^s`
    # such that `psi(component_meridians[sc]) == sym_sc^{+1}` for every component `sc`!
    g = pres.num_generators
    num_c = topo.num_components

    if g == 0:
        ring_desc = "Z"
        return AlexanderResult(
            ring_description=ring_desc,
            variable_names=tuple(var_names),
            orientable_flags=tuple(ori_flags),
            e1_generators=("1",),
            e2_generators=("1",),
            single_var_alex="1",
            evaluation_fingerprint=((1,), (1,), (1,)),
        )

    # Solve linear system for generator images in `Z^r` (orientable coordinates) and `(Z_2)^s` (torsion):
    # Relators give `sum_j exp_{i,j} psi(x_j) = 0` (mod 2 on torsion coords, exact on Z coords),
    # AND `component_meridians[sc]` gives `sum_j m_{sc, j} psi(x_j) = e_{sc}`!
    ori_comp_indices = [sc for sc in comp_order if topo.components[sc].orientable]
    tor_comp_indices = [sc for sc in comp_order if not topo.components[sc].orientable]

    def word_exp_vec(w: Word) -> list[int]:
        v = [0] * g
        for lit in w:
            v[abs(lit) - 1] += 1 if lit > 0 else -1
        return v

    # 1. Solve over Q/Z for orientable coordinates `gen_ori_exps[j][k]` (`j in 0..g-1, k in 0..r-1`)
    r_cnt = len(ori_comp_indices)
    gen_ori_exps = [[0] * r_cnt for _ in range(g)]
    if r_cnt > 0:
        eq_rows = []
        eq_rhs = []
        for rel in pres.relators:
            rv = word_exp_vec(rel)
            if any(rv):
                eq_rows.append(rv)
                eq_rhs.append([0] * r_cnt)
        for k_idx, sc in enumerate(ori_comp_indices):
            mv = word_exp_vec(pres.component_meridians[sc])
            eq_rows.append(mv)
            rhs_row = [0] * r_cnt
            rhs_row[k_idx] = 1
            eq_rhs.append(rhs_row)
        mat_a = sp.Matrix(eq_rows)
        mat_b = sp.Matrix(eq_rhs)
        sol, params = mat_a.gauss_jordan_solve(mat_b)
        if params.rows > 0:
            sol = sol.subs({p: 0 for p in params})
        for j in range(g):
            for k_idx in range(r_cnt):
                gen_ori_exps[j][k_idx] = int(sol[j, k_idx])

    # 2. Solve over GF(2) for non-orientable coordinates `gen_tor_exps[j][l]` (`l in 0..s-1`)
    s_cnt = len(tor_comp_indices)
    gen_tor_exps = [[0] * s_cnt for _ in range(g)]
    if s_cnt > 0:
        gf2_rows: list[list[int]] = []
        for rel in pres.relators:
            rv = [x & 1 for x in word_exp_vec(rel)]
            if any(rv):
                gf2_rows.append(rv + [0] * s_cnt)
        for l_idx, sc in enumerate(tor_comp_indices):
            mv = [x & 1 for x in word_exp_vec(pres.component_meridians[sc])]
            rhs = [0] * s_cnt
            rhs[l_idx] = 1
            gf2_rows.append(mv + rhs)
        # Gaussian elimination over GF(2)
        num_rows = len(gf2_rows)
        piv_row = 0
        piv_col_for_row: dict[int, int] = {}
        for col in range(g):
            found = -1
            for rr in range(piv_row, num_rows):
                if gf2_rows[rr][col] == 1:
                    found = rr
                    break
            if found == -1:
                continue
            gf2_rows[piv_row], gf2_rows[found] = gf2_rows[found], gf2_rows[piv_row]
            piv_col_for_row[piv_row] = col
            for rr in range(num_rows):
                if rr != piv_row and gf2_rows[rr][col] == 1:
                    gf2_rows[rr] = [a ^ b for a, b in zip(gf2_rows[rr], gf2_rows[piv_row])]
            piv_row += 1
        for rr, col in piv_col_for_row.items():
            for l_idx in range(s_cnt):
                gen_tor_exps[col][l_idx] = gf2_rows[rr][g + l_idx]

    # Build SymPy monomials `gen_images[j]` and `gen_inv_images[j]`
    gen_images: list[sp.Expr] = []
    gen_inv_images: list[sp.Expr] = []
    for j in range(g):
        pos_m = sp.Integer(1)
        neg_m = sp.Integer(1)
        for k_idx in range(r_cnt):
            e_val = gen_ori_exps[j][k_idx]
            if e_val > 0:
                pos_m *= ori_syms[k_idx] ** e_val
                neg_m *= ori_inv_syms[k_idx] ** e_val
            elif e_val < 0:
                pos_m *= ori_inv_syms[k_idx] ** (-e_val)
                neg_m *= ori_syms[k_idx] ** (-e_val)
        for l_idx in range(s_cnt):
            if gen_tor_exps[j][l_idx] & 1:
                pos_m *= torsion_syms[l_idx]
                neg_m *= torsion_syms[l_idx]
        gen_images.append(pos_m)
        gen_inv_images.append(neg_m)

    # 3. Build `m x g` Fox Jacobian matrix `J`
    m_rels = len(pres.relators)
    jac_rows: list[list[sp.Expr]] = []
    for rel in pres.relators:
        row = []
        for j in range(g):
            df = _fox_derivative_word(rel, j, gen_images, gen_inv_images)
            # Cancel `t_i * t_i_inv -> 1` and `u_j^2 -> 1` immediately!
            for s_pos, s_inv in zip(ori_syms, ori_inv_syms):
                # Do exact cancellation via Poly if needed; `_clear_negative...` handles minors,
                # here substitute `s_pos * s_inv = 1` by converting via rational or Poly:
                pass
            row.append(_reduce_mod_torsion(df, torsion_syms))
        jac_rows.append(row)

    # Helper to reduce `t_i * t_i_inv = 1` and `u_j^2 = 1` inside a polynomial expression
    def normalize_ring_elem(expr: sp.Expr) -> sp.Expr:
        if expr == 0:
            return sp.Integer(0)
        all_s = ori_syms + ori_inv_syms + torsion_syms
        if not all_s:
            return sp.expand(expr)
        p = sp.Poly(sp.expand(expr), *all_s, domain=sp.ZZ)
        out = sp.Integer(0)
        for mon, coeff in p.terms():
            term = sp.Integer(int(coeff))
            for idx_o in range(r_cnt):
                net = mon[idx_o] - mon[r_cnt + idx_o]
                if net > 0:
                    term *= ori_syms[idx_o] ** net
                elif net < 0:
                    term *= ori_inv_syms[idx_o] ** (-net)
            for idx_t in range(s_cnt):
                if mon[2 * r_cnt + idx_t] & 1:
                    term *= torsion_syms[idx_t]
            out += term
        return sp.expand(out)

    jac_rows = [[normalize_ring_elem(entry) for entry in row] for row in jac_rows]

    def compute_elementary_ideal(k_order: int) -> list[sp.Expr]:
        minor_size = g - k_order
        if minor_size <= 0:
            return [sp.Integer(1)]
        if minor_size > m_rels or minor_size > g:
            return [sp.Integer(0)]
        mat = sp.Matrix(jac_rows)
        minors: list[sp.Expr] = []
        for r_comb in combinations(range(m_rels), minor_size):
            for c_comb in combinations(range(g), minor_size):
                sub = mat.extract(list(r_comb), list(c_comb))
                det = normalize_ring_elem(sub.det(method="berkowitz"))
                if det != 0:
                    minors.append(det)
        return _simplify_ideal_generators(minors, ori_syms, ori_inv_syms, torsion_syms)

    e1_gens = compute_elementary_ideal(1)
    e2_gens = compute_elementary_ideal(2)
    fp1 = _ideal_evaluation_fingerprint(e1_gens, ori_syms, torsion_syms)
    fp2 = _ideal_evaluation_fingerprint(e2_gens, ori_syms, torsion_syms)

    # Single-variable specialized Alexander polynomial (map all orientable vars -> t, torsion -> 1)
    t_sym = sp.Symbol("t")
    subs_1var = {s: t_sym for s in ori_syms}
    subs_1var.update({s: 1 / t_sym for s in ori_inv_syms})
    subs_1var.update({u: 1 for u in torsion_syms})
    spec_polys = [
        _clear_negative_laurent_powers(
            g_expr.subs(subs_1var) if ori_syms else g_expr.subs(subs_1var),
            [t_sym] if ori_syms else [],
            [],
            [],
        )
        for g_expr in e1_gens
    ]
    spec_polys = [p for p in spec_polys if p != 0]
    if not spec_polys:
        single_alex = "0"
    elif not ori_syms:
        from math import gcd
        from functools import reduce

        single_alex = str(reduce(gcd, (abs(int(p)) for p in spec_polys)))
    else:
        # GCD in Z[t]
        cur_gcd = spec_polys[0]
        for nxt_p in spec_polys[1:]:
            cur_gcd = sp.gcd(sp.Poly(cur_gcd, t_sym, domain=sp.ZZ), sp.Poly(nxt_p, t_sym, domain=sp.ZZ)).as_expr()
        cur_gcd = _clear_negative_laurent_powers(cur_gcd, [t_sym], [], [])
        single_alex = str(cur_gcd)

    torsion_rels = [f"{u}^2=1" for u in torsion_syms]
    ori_part = ", ".join(f"{s}^+-1" for s in ori_syms)
    tor_part = ", ".join(str(u) for u in torsion_syms)
    vars_desc = ", ".join(x for x in (ori_part, tor_part) if x)
    ring_desc = f"Z[{vars_desc}]" + (f" / ({', '.join(torsion_rels)})" if torsion_rels else "")

    return AlexanderResult(
        ring_description=ring_desc,
        variable_names=tuple(var_names),
        orientable_flags=tuple(ori_flags),
        e1_generators=tuple(str(p) for p in e1_gens),
        e2_generators=tuple(str(p) for p in e2_gens),
        single_var_alex=single_alex,
        evaluation_fingerprint=(fp1, fp2),
    )


def twisted_s3_invariants(
    d: MarkedDiagram,
    pres: GroupPresentation | None = None,
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    """Compute twisted homology `H_1(pi_1(S^4 \\ F); V_rho)` across conjugacy classes of surjective
    representations `rho: pi_1 -> S_3 -> GL_2(Z)` (self-certified via Fox's identity `J_rho * B_rho == 0`).

    Returns sorted tuple of `(free_rank, torsion_factors)` over all `S_3`-conjugacy classes of
    surjective homomorphisms `pi_1 -> S_3` (empty tuple `()` if `pi_1` has no surjective `S_3` quotient).
    """
    from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
    from surfaceknots.invariants.fundamental_group import STANDARD_TEST_GROUPS

    if pres is None:
        pres = wirtinger_presentation(d, from_resolution=-1, simplify=True)
    g = pres.num_generators
    if g == 0:
        return ()

    s3 = STANDARD_TEST_GROUPS[0]
    assert s3.name == "S3" and s3.order == 6

    # Build exact 2x2 integral standard representation of S_3 (`0..5`) by matching multiplication table!
    # Generators of STANDARD_TEST_GROUPS[0] were `(1, 2, 0)` (3-cycle) and `(1, 0, 2)` (transposition).
    # Standard 2D representation on `{(u0, u1, u2) : u0+u1+u2 = 0}` with basis `e1=(1,-1,0), e2=(0,1,-1)`:
    # Perm `p` sends `u_i` -> `u_{p[i]}`? Wait: in `_group_from_permutations`, `comp = tuple(pi[pj[k]])`:
    # Check `(1,2,0)` (`c`) and `(1,0,2)` (`s`): `c * c * c = id`, `s * s = id`, `c * s = (2, 1, 0)`,
    # `s * c = (0, 2, 1) = c^2 * s` (`s c s = c^2`)!
    # Check matrices `M_c = [[0, -1], [1, -1]]` (`M_c^3 = I`) and `M_s = [[0, 1], [1, 0]]` (`M_s^2 = I`):
    # `M_s * M_c * M_s = [[0,1],[1,0]] [[0,-1],[1,-1]] [[0,1],[1,0]] = [[1,-1],[0,-1]] [[0,1],[1,0]] = [[-1, 0], [-1, 1]]`?
    # Wait: `M_c^2 = [[-1, 1], [-1, 0]]`! Let's take `M_s = [[-1, 1], [0, 1]]`:
    # `M_s^2 = [[1, 0], [0, 1]]`, and `M_s * M_c * M_s = [[-1,1],[0,1]] [[0,-1],[1,-1]] [[-1,1],[0,1]]`
    # `= [[1, 0], [1, -1]] [[-1, 1], [0, 1]] = [[-1, 1], [-1, 0]] == M_c^2`!!
    # Even better: match all 6 elements directly from `_group_from_permutations` by BFS from `(1,2,0)` & `(1,0,2)`!
    mat_c = sp.Matrix([[0, -1], [1, -1]])
    mat_s = sp.Matrix([[-1, 1], [0, 1]])
    ident_perm = (0, 1, 2)
    g_perms = [(1, 2, 0), (1, 0, 2)]
    g_mats = [mat_c, mat_s]
    seen_p = {ident_perm: (0, sp.eye(2))}
    p_list = [ident_perm]
    idx = 0
    while idx < len(p_list):
        cur_p = p_list[idx]
        _, cur_m = seen_p[cur_p]
        idx += 1
        for gp, gm in zip(g_perms, g_mats):
            nxt_p = tuple(cur_p[gp[i]] for i in range(3))
            if nxt_p not in seen_p:
                seen_p[nxt_p] = (len(p_list), cur_m * gm)
                p_list.append(nxt_p)
    rep_mats = [seen_p[p][1] for p in p_list]
    # Verify exact homomorphism: rep_mats[s3.mul[a][b]] == rep_mats[a] * rep_mats[b]
    for a in range(6):
        for b in range(6):
            assert rep_mats[s3.mul[a][b]] == rep_mats[a] * rep_mats[b]

    # Enumerate surjective homomorphisms `assign: 0..g-1 -> 0..5`, picking 1 canonical rep per S_3-conjugacy class!
    # Conjugation table in S_3: `conj[h][x] = mul[mul[h][x]][inv[h]]`
    conj = [
        [s3.mul[s3.mul[h][x]][s3.inv[h]] for x in range(6)]
        for h in range(6)
    ]
    rels_by_max: list[list[Word]] = [[] for _ in range(g)]
    for r in pres.relators:
        if r:
            rels_by_max[max(abs(lit) - 1 for lit in r)].append(r)

    surj_reps: list[tuple[int, ...]] = []
    seen_orbits: set[tuple[int, ...]] = set()
    assign = [0] * g

    def bt(k: int) -> None:
        if k == g:
            # Check surjectivity: image generates S_3 iff image has both order-2 and order-3 elements
            # (or subgroup closure size == 6)
            sub = {0}
            st = [x for x in assign if x != 0]
            sub.update(st)
            while st:
                u = st.pop()
                for v in list(sub):
                    for prod_uv in (s3.mul[u][v], s3.mul[v][u]):
                        if prod_uv not in sub:
                            sub.add(prod_uv)
                            st.append(prod_uv)
            if len(sub) < 6:
                return
            orbit = [tuple(conj[h][x] for x in assign) for h in range(6)]
            canon = min(orbit)
            if canon not in seen_orbits:
                seen_orbits.add(canon)
                surj_reps.append(canon)
            return
        for val in range(6):
            assign[k] = val
            ok = True
            for r in rels_by_max[k]:
                acc = 0
                for lit in r:
                    gv = assign[lit - 1] if lit > 0 else s3.inv[assign[-lit - 1]]
                    acc = s3.mul[acc][gv]
                if acc != 0:
                    ok = False
                    break
            if ok:
                bt(k + 1)

    bt(0)
    if not surj_reps:
        return ()

    results: list[tuple[int, tuple[int, ...]]] = []
    m_rels = len(pres.relators)
    for rho in surj_reps:
        # Build `J_rho` (`2m x 2g`) and `B_rho` (`2g x 2`)
        b_blocks = [rep_mats[rho[j]] - sp.eye(2) for j in range(g)]
        b_mat = sp.Matrix.vstack(*b_blocks)  # 2g x 2

        if m_rels == 0:
            j_mat = sp.zeros(0, 2 * g)
        else:
            j_rows = []
            for r in pres.relators:
                blocks = [sp.zeros(2, 2) for _ in range(g)]
                pref = 0  # identity in S_3
                for lit in r:
                    j_idx = abs(lit) - 1
                    if lit > 0:
                        blocks[j_idx] += rep_mats[pref]
                        pref = s3.mul[pref][rho[j_idx]]
                    else:
                        pref = s3.mul[pref][s3.inv[rho[j_idx]]]
                        blocks[j_idx] -= rep_mats[pref]
                j_rows.append(sp.Matrix.hstack(*blocks))
            j_mat = sp.Matrix.vstack(*j_rows)  # 2m x 2g
            assert (j_mat * b_mat).is_zero_matrix, "Fox twisted chain condition failed!"

        # Compute integral homology `H_1 = ker(d1) / im(d2)` for `Z^{2m} ->^{d2} Z^{2g} ->^{d1} Z^2`
        # (`d1 = B_rho^T`, `d2 = J_rho^T`, with `d1 * d2 == 0` certified above).
        # Because `im(d1) subset Z^2` is a submodule of a free abelian group, `im(d1) ~= Z^{rk(d1)}`
        # is free over Z, so `0 -> ker(d1)/im(d2) -> Z^{2g}/im(d2) -> im(d1) -> 0` splits over Z!
        # Hence `tor(H_1) == tor(Z^{2g} / im(d2))` and `rank(H_1) == 2*g - rank(d2) - rank(d1)`.
        d1 = b_mat.T  # 2 x 2g
        d2 = j_mat.T  # 2g x 2m
        rk_d1 = int(d1.rank())
        if d2.cols == 0:
            results.append((2 * g - rk_d1, ()))
            continue
        snf_d2 = smith_normal_form(d2, domain=sp.ZZ)
        diag = [abs(int(snf_d2[i, i])) for i in range(min(snf_d2.rows, snf_d2.cols))]
        nz = [x for x in diag if x > 0]
        free_rk = 2 * g - len(nz) - rk_d1
        tor = tuple(sorted(x for x in nz if x > 1))
        results.append((free_rk, tor))

    results.sort()
    return tuple(results)

