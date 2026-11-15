"""Independent helpers for audit D (free-group words, tangle kernel test, Kauffman bracket)."""
import itertools, pickle, os

def red(w):
    out = []
    for x in w:
        if out and out[-1] == -x:
            out.pop()
        else:
            out.append(x)
    return tuple(out)

def inv(w):
    return tuple(-x for x in reversed(w))

def artin_up(words, bw):
    """Artin action: +k: (u,v)->(u v u^-1, u); -k: (u,v)->(v, v^-1 u v)."""
    w = list(words)
    for g in bw:
        k = abs(g) - 1
        u, v = w[k], w[k + 1]
        if g > 0:
            w[k], w[k + 1] = red(u + v + inv(u)), u
        else:
            w[k], w[k + 1] = v, red(inv(v) + u + v)
    return w

def top_words(b, matching, bw):
    """phi_T: top port p -> word in F_b (arch r generator r+1)."""
    wb = [()] * (2 * b)
    for r, (u, v) in enumerate(matching):
        lo, hi = min(u, v), max(u, v)
        wb[lo] = (r + 1,)
        wb[hi] = (-(r + 1),)
    return artin_up(wb, bw)

def cap_words(b, matching, bw):
    """Normal generators of ker(phi_T) as words in top-port gens 1..2b."""
    wp = [(p + 1,) for p in range(2 * b)]
    wbot = artin_up(wp, tuple(-g for g in reversed(bw)))
    return [red(wbot[min(u, v)] + wbot[max(u, v)]) for u, v in matching]

def subst(word, images):
    out = []
    for lit in word:
        img = images[abs(lit) - 1]
        if lit < 0:
            img = inv(img)
        for x in img:
            if out and out[-1] == -x:
                out.pop()
            else:
                out.append(x)
    return tuple(out)

def isotopic(capsA, topB):
    """ker(phi_A) subset ker(phi_B)  (== equality by Hopficity of F_b)."""
    return all(len(subst(c, topB)) == 0 for c in capsA)
