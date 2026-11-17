"""Test candidate 4D geometric constructions against 12_2."""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path("/afs/csail.mit.edu/u/e/ejli/.gemini/antigravity-cli/worktrees/1dadd8f0-f304-455a-b682-76a831212ce6/subagent-4D-Topology---Mystery-Knot-Analyst-self-69be3150")
sys.path.insert(0, str(REPO / "code" / "src"))
sys.path.insert(0, "/afs/csail.mit.edu/u/e/ejli/.gemini/antigravity-cli/brain/ee536e5c-0e8d-469c-b4f6-08d6f2745f37/scratch")

from audit_12_2 import todd_coxeter_order
from surfaceknots.diagrams import MarkedDiagram

TARGET_ORDERS = {2: 2, 3: 3, 4: 108, 5: 5, 6: 720, 7: 2352}


def word_pow(w: tuple[int, ...], e: int) -> tuple[int, ...]:
    if e == 0:
        return ()
    if e > 0:
        return w * e
    inv = tuple(-x for x in reversed(w))
    return inv * (-e)


def word_inv(w: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(-x for x in reversed(w))


def test_symmetry_spins_composite_trefoil():
    print("=== 1. Testing Litherland Symmetry-Spins of 3_1 # 3_1 and 3_1 # -3_1 ===")
    # Gens: 1 = a, 2 = b0, 3 = b1
    # Trefoil T0: a b0 a = b0 a b0  => (1, 2, 1, -2, -1, -2)
    rel_t0 = (1, 2, 1, -2, -1, -2)
    for mirror in (False, True):
        # Mirror T1: a^-1 b1^-1 a^-1 = b1^-1 a^-1 b1^-1 (same relation!) or inverse crossing: a b1^-1 a^-1 ...
        # Wait: mirror of a b a = b a b with meridian a is a^-1 b a^-1 = b a^-1 b^-1? No:
        # Right-handed trefoil: b_out = a b_in a^-1 ... standard Wirtinger of mirror trefoil simply reverses all crossings,
        # which replaces a -> a^-1, b1 -> b1^-1, giving a^-1 b1^-1 a^-1 = b1^-1 a^-1 b1^-1 <=> b1 a b1 = a b1 a (IDENTICAL!).
        rel_t1 = (1, 3, 1, -3, -1, -3)
        # Half-turn symmetry swaps b0 <-> b1 (or b0 -> b1, b1 -> a^k b0 a^-k), composed with conjugation by u0 = a^c * (a*b0)^(3*d)
        for c in range(-6, 7):
            for d in range(-2, 3):
                u0 = word_pow((1,), c) + word_pow((1, 2, 1, 2, 1, 2), d)
                u1 = word_pow((1,), c) + word_pow((1, 3, 1, 3, 1, 3), d)
                # Monodromy sends b0 -> u0 b1 u0^-1 and b1 -> u1 b0 u1^-1 (or u0 b0 u0^-1)
                for mode in ("swap_u0_u1", "swap_u0_id"):
                    if mode == "swap_u0_u1":
                        r_sym1 = (-2,) + u0 + (3,) + word_inv(u0)
                        r_sym2 = (-3,) + u1 + (2,) + word_inv(u1)
                    else:
                        r_sym1 = (-2,) + u0 + (3,) + word_inv(u0)
                        r_sym2 = (-3, 2)  # if second half is untwisted, then b1 = b0 => u0 b0 u0^-1 = b0
                    rels = [rel_t0, rel_t1, r_sym1, r_sym2]
                    try:
                        o2 = todd_coxeter_order(3, rels + [(1, 1)], max_cosets=20000)
                        o3 = todd_coxeter_order(3, rels + [(1, 1, 1)], max_cosets=20000)
                        o4 = todd_coxeter_order(3, rels + [(1, 1, 1, 1)], max_cosets=20000)
                        if (o2, o3, o4) == (2, 3, 108):
                            full = {m: todd_coxeter_order(3, rels + [(1,) * m], max_cosets=100000) for m in range(2, 8)}
                            print(f"  MATCH! mirror={mirror}, mode={mode}, c={c}, d={d}: {full}")
                    except RuntimeError:
                        pass


def test_symmetry_spins_2bridge():
    print("\n=== 2. Testing Teragaito Twisting Symmetry-Spins of 2-Bridge Knots K(alpha, beta) ===")
    # Schubert 2-bridge knot K(alpha, beta), alpha odd, 0 < beta < alpha, gcd(alpha, beta) = 1:
    # Wait: a 2-bridge knot K(alpha, beta) has period-2 symmetry exchanging the two bridges x <-> z
    # WHEN beta^2 == +-1 (mod alpha)! (By Schubert / Burde-Zieschang: K(alpha, beta) is invertible always;
    # the orientation-preserving involution reversing bridges vs shifting along K!).
    # Wait: EVERY 2-bridge knot has an orientation-preserving involution h of (S^3, K) with Fix(h) cap K = empty
    # (period 2 free on K) WHEN beta^2 == 1 (mod alpha)? Wait: for 4_1 = K(5, 2), 2^2 = 4 == -1 (mod 5)!
    # Let's build exact Schubert 2-bridge presentation <x, z | w x = z w> (or <x, z | w x = x w> depending on parity)
    # and test ALL relations identifying z with x^c w^d ...!
    from math import gcd
    for alpha in range(3, 21, 2):
        for beta in range(1, alpha):
            if gcd(alpha, beta) != 1:
                continue
            # Schubert word w = z^{eps_1} x^{eps_2} ... z^{eps_{alpha-1}} (since alpha-1 is even!):
            eps = []
            for i in range(1, alpha):
                val = (i * beta) % (2 * alpha)
                eps.append(1 if (0 < val < alpha) else -1)
            # Build word w in gens 1=x, 2=z:
            # Standard Schubert: w = x^{e_1} z^{e_2} x^{e_3} ... z^{e_{alpha-1}}
            w = []
            for idx, e in enumerate(eps):
                gen = 1 if (idx % 2 == 0) else 2
                w.append(gen if e > 0 else -gen)
            w = tuple(w)
            # Knot relation: w * x * w^-1 * z^-1 = 1 (Wait: check abelianization: w has equal x and z? Let's check!)
            # Actually in Schubert presentation of 2-bridge knot: w x = z w where w = x^{e_1} z^{e_2} ... x^{e_{alpha-2}} z^{e_{alpha-1}}?
            # Wait: if rel is w * x * w^-1 * z^-1 = 1, abelianization gives x - z = 0 => H_1 = Z ALWAYS!
            # Wait: check for alpha=3, beta=1: eps for i=1,2: 1, 2 -> (+1, +1), w = x z.
            # w x = z w => x z x = z x z (3_1!!).
            # Check for alpha=5, beta=2: i*2 mod 10 for i=1..4: 2, 4, 6, 8 -> (+1, +1, -1, -1)?
            # Wait: for 4_1 = K(5, 3) or K(5, 2): let's check K(5, 3): 3, 6, 9, 12%10=2 -> (+1, -1, -1, +1):
            # w = x z^-1 x^-1 z => w x = z w <=> (x z^-1 x^-1 z) x = z (x z^-1 x^-1 z) (4_1!!).
            rel_knot = w + (1,) + word_inv(w) + (-2,)
            # Test adding ((z x^-1)^k * x^m) relations or symmetry-spin relations!
            # In Teragaito (1989) "Twisting symmetry-spins of 2-bridge knots":
            # What is the group presentation of the m-twisting symmetry-spin of K(alpha, beta) (when beta^2 == 1 mod alpha, or +-1)?
            # Wait: let's check (z x^-1)^k = 1 and also x^c = ...
            for k in (2, 3, 4):
                r2 = (2, -1) * k
                rels = [rel_knot, r2]
                try:
                    o2 = todd_coxeter_order(2, rels + [(1, 1)], max_cosets=10000)
                    o3 = todd_coxeter_order(2, rels + [(1, 1, 1)], max_cosets=10000)
                    o4 = todd_coxeter_order(2, rels + [(1, 1, 1, 1)], max_cosets=10000)
                    if (o2, o3, o4) == (2, 3, 108):
                        full = {m: todd_coxeter_order(2, rels + [(1,) * m], max_cosets=100000) for m in range(2, 8)}
                        print(f"  2-bridge quotient K({alpha},{beta}) / (z x^-1)^{k}=1 MATCHES: {full}")
                except RuntimeError:
                    pass


def inspect_12_2_movie():
    print("\n=== 3. Inspecting Exact Resolutions & Saddle Attachment of 12_2.a and 12_2.b ===")
    ch12 = json.loads((REPO / "code" / "data" / "expanded_ch12.json").read_text())
    e = next(x for x in ch12["entries"] if x["label"] == "12_2")
    for tag, epd in [("12_2.a", e["epd"]), ("12_2.b", e["alternate_epds"][0])]:
        d = MarkedDiagram.from_epd(epd)
        print(f"\n{tag}: {epd}")
        for sgn, name in [(-1, "L_-"), (+1, "L_+")]:
            res = d.resolve(sgn)
            print(f"  {name}: {res}")


if __name__ == "__main__":
    test_symmetry_spins_composite_trefoil()
    test_symmetry_spins_2bridge()
    inspect_12_2_movie()
