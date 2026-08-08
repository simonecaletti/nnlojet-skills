#!/usr/bin/env python3
"""Reduced-Born rule as an algorithm: classify every infrared mode of a
channel as GENUINE or DEAD.

Usage:  python genuine_modes.py spec.json [--json]
        python genuine_modes.py --selftest

Spec:
  partons : {"<momentum index>": "<flavour>", ...} final-state partons.
            Flavours: "g" (gluon), "q<tag>" (quark), "qb<tag>"
            (antiquark) — <tag> is any flavour label ("1", "2", "u"...);
            same tag = same flavour.
  born    : list of legal Born final states (flavour lists), e.g.
            [["q1","qb1","g"]]. Tags are matched up to renaming.
  families: optional subset of
            ["ss","sco","ds","tc","sc","dc"]; default all.

The rule (empirically validated, 240/240 mode classifications): a
configuration is a genuine limit iff, after replacing every collinear
cluster by its parent parton (a cluster is valid iff its NET flavour is
a single parton) and deleting every soft parton (a gluon, or a
same-flavour q qb PAIR going soft together), what remains is a legal
radiative state of the process — implemented here as: reducible to a
Born state by further flavour-valid collapses/deletions. Consequences
the implementation reproduces by construction: a single quark never
goes soft; composites are NOT decomposed pairwise; two individually
genuine pairs can combine into a dead double limit if collapsing both
breaks the Born.
"""
import itertools
import json
import sys
from functools import lru_cache


# ---------------- flavour algebra ----------------

def parse(f):
    """flavour string -> ('g', None) | ('q', tag) | ('qb', tag)"""
    if f == "g":
        return ("g", None)
    if f.startswith("qb"):
        return ("qb", f[2:])
    if f.startswith("q"):
        return ("q", f[1:])
    raise ValueError(f"bad flavour {f!r} (want g / q<tag> / qb<tag>)")


def canonical(flavs):
    """Multiset signature, invariant under flavour-tag renaming:
    (ngluons, sorted per-tag (nq, nqb) tuples)."""
    ng = 0
    per = {}
    for f in flavs:
        k, t = parse(f)
        if k == "g":
            ng += 1
        else:
            a, b = per.get(t, (0, 0))
            per[t] = (a + 1, b) if k == "q" else (a, b + 1)
    return (ng, tuple(sorted(per.values())))


def cluster_parent(flavs):
    """Net flavour of a collinear cluster: 'g', 'q<tag>'/'qb<tag>', or
    None if the net is not a single parton (invalid cluster)."""
    net = {}
    for f in flavs:
        k, t = parse(f)
        if k == "q":
            net[t] = net.get(t, 0) + 1
        elif k == "qb":
            net[t] = net.get(t, 0) - 1
    nonzero = {t: n for t, n in net.items() if n != 0}
    total = sum(abs(n) for n in nonzero.values())
    if total == 0:
        return "g"
    if total == 1:
        (t, n), = nonzero.items()
        return f"q{t}" if n > 0 else f"qb{t}"
    return None


def soft_deletable(flavs):
    """May this set go soft TOGETHER? A single gluon, or a same-flavour
    q qb pair, or two gluons (as a double-soft pair)."""
    ks = sorted(parse(f)[0] for f in flavs)
    if ks == ["g"] or ks == ["g", "g"]:
        return True
    if ks == ["q", "qb"]:
        (t1,) = {parse(f)[1] for f in flavs if parse(f)[0] == "q"}
        (t2,) = {parse(f)[1] for f in flavs if parse(f)[0] == "qb"}
        return t1 == t2
    return False


def make_legal_checker(born):
    """legal(state): reducible to a Born multiset by further
    flavour-valid collapses / soft deletions (memoised closure)."""
    born_sigs = {canonical(b) for b in born}

    @lru_cache(maxsize=None)
    def legal(state):  # state: sorted tuple of flavour strings
        if canonical(state) in born_sigs:
            return True
        fl = list(state)
        n = len(fl)
        # delete a soft gluon
        for i in range(n):
            if fl[i] == "g":
                if legal(tuple(sorted(fl[:i] + fl[i + 1:]))):
                    return True
        # delete a same-flavour q qb pair
        for i, j in itertools.combinations(range(n), 2):
            if soft_deletable([fl[i], fl[j]]) and "g" not in (fl[i], fl[j]):
                rest = [f for k, f in enumerate(fl) if k not in (i, j)]
                if legal(tuple(sorted(rest))):
                    return True
        # merge a valid 2-cluster
        for i, j in itertools.combinations(range(n), 2):
            p = cluster_parent([fl[i], fl[j]])
            if p is not None:
                rest = [f for k, f in enumerate(fl) if k not in (i, j)]
                if legal(tuple(sorted(rest + [p]))):
                    return True
        return False

    return legal


# ---------------- mode classification ----------------

def classify(partons, born, families=None):
    """partons: {index(int): flavour}. -> list of mode dicts."""
    families = families or ["ss", "sco", "ds", "tc", "sc", "dc"]
    idx = sorted(partons)
    legal = make_legal_checker([tuple(sorted(b)) for b in born])

    def state_after(collapsed=(), deleted=()):
        rest = [partons[i] for i in idx
                if i not in set().union(*[set(c) for c in collapsed] or [set()])
                and i not in deleted]
        for c in collapsed:
            rest.append(cluster_parent([partons[i] for i in c]))
        return tuple(sorted(rest))

    out = []

    def add(fam, name, ok, reason):
        out.append({"family": fam, "name": name,
                    "genuine": bool(ok), "reason": reason})

    for fam in families:
        if fam == "ss":
            for i in idx:
                if partons[i] != "g":
                    add(fam, f"{i} soft", False, "single quark never soft")
                    continue
                ok = legal(state_after(deleted=(i,)))
                add(fam, f"{i} soft", ok,
                    "remaining state legal" if ok else "no Born reachable")
        elif fam == "sco":
            for i, j in itertools.combinations(idx, 2):
                p = cluster_parent([partons[i], partons[j]])
                if p is None:
                    add(fam, f"{i}||{j}", False, "cluster net not single parton")
                    continue
                ok = legal(state_after(collapsed=((i, j),)))
                add(fam, f"{i}||{j}", ok,
                    f"parent {p}" if ok else "collapsed state not legal")
        elif fam == "ds":
            for i, j in itertools.combinations(idx, 2):
                if not soft_deletable([partons[i], partons[j]]):
                    add(fam, f"{i},{j} soft", False,
                        "pair not soft-deletable (gg or same-flavour qqb)")
                    continue
                ok = legal(state_after(deleted=(i, j)))
                add(fam, f"{i},{j} soft", ok,
                    "remaining state legal" if ok else
                    "deleting pair breaks Born")
        elif fam == "tc":
            for i, j, k in itertools.combinations(idx, 3):
                p = cluster_parent([partons[i], partons[j], partons[k]])
                if p is None:
                    add(fam, f"{i}||{j}||{k}", False,
                        "cluster net not single parton")
                    continue
                ok = legal(state_after(collapsed=((i, j, k),)))
                add(fam, f"{i}||{j}||{k}", ok,
                    f"parent {p}" if ok else "collapsed state not legal")
        elif fam == "sc":
            for i in idx:
                for j, k in itertools.combinations([x for x in idx if x != i], 2):
                    nm = f"{i} soft + {j}||{k}"
                    if partons[i] != "g":
                        add(fam, nm, False, "single quark never soft")
                        continue
                    p = cluster_parent([partons[j], partons[k]])
                    if p is None:
                        add(fam, nm, False, "cluster net not single parton")
                        continue
                    ok = legal(state_after(collapsed=((j, k),), deleted=(i,)))
                    add(fam, nm, ok,
                        "state legal" if ok else "no Born reachable")
        elif fam == "dc":
            for (i, j), (k, l) in itertools.combinations(
                    itertools.combinations(idx, 2), 2):
                if {i, j} & {k, l}:
                    continue
                nm = f"{i}||{j} + {k}||{l}"
                p1 = cluster_parent([partons[i], partons[j]])
                p2 = cluster_parent([partons[k], partons[l]])
                if p1 is None or p2 is None:
                    add(fam, nm, False, "a cluster net is not single parton")
                    continue
                ok = legal(state_after(collapsed=((i, j), (k, l))))
                add(fam, nm, ok,
                    "state legal" if ok else
                    "collapsing both breaks Born (pairs may be individually genuine)")
    return out


# ---------------- selftest ----------------

def selftest():
    """Checks ALGORITHM invariants and the rule's defining clauses on
    synthetic content. Encodes no real process's mode answer."""
    partons = {3: "q1", 4: "qb1", 5: "g", 6: "q2", 7: "qb2"}
    born = [["q1", "qb1", "g"]]
    res = classify(partons, born)

    # 1. permutation covariance: relabel momenta, classification follows
    import re

    def canon_key(fam, name, perm=None):
        """(family, canonical index structure) — invariant under the
        name's own ordering conventions."""
        nums = [int(x) for x in re.findall(r"\d+", name)]
        if perm:
            nums = [perm[x] for x in nums]
        if fam in ("sco", "ds", "tc", "ss"):
            key = tuple(sorted(nums))
        elif fam == "sc":                       # soft index + collinear pair
            key = (nums[0], tuple(sorted(nums[1:])))
        elif fam == "dc":                       # unordered pair of pairs
            p1, p2 = sorted(nums[:2]), sorted(nums[2:])
            key = tuple(sorted([tuple(p1), tuple(p2)]))
        else:
            raise ValueError(fam)
        return (fam, key)

    perm = {3: 7, 4: 6, 5: 5, 6: 4, 7: 3}
    partons_p = {perm[i]: f for i, f in partons.items()}
    res_p = classify(partons_p, born)
    a = {canon_key(m["family"], m["name"], perm): m["genuine"] for m in res}
    b = {canon_key(m["family"], m["name"]): m["genuine"] for m in res_p}
    assert a == b, "not permutation-covariant"
    as_map = lambda r: {(m["family"], m["name"]): m["genuine"] for m in r}
    a = as_map(res)

    # 2. flavour-tag renaming invariance
    ren = {"q1": "qA", "qb1": "qbA", "q2": "qB", "qb2": "qbB", "g": "g"}
    res_r = classify({i: ren[f] for i, f in partons.items()},
                     [[ren[f] for f in born[0]]])
    assert as_map(res_r) == a, "not tag-renaming invariant"

    # 3. defining clauses of the rule (not process answers):
    for m in res:
        if m["family"] == "ss" and "quark" in m["reason"]:
            assert not m["genuine"]           # a single quark never soft
    # cross-flavour collinear pair is dead
    assert a[("sco", "4||6")] is False        # qb1 || q2: net two quarks
    # Born + one soft gluon: the soft-gluon modes are genuine by the
    # definition of the closure (not a process-specific answer)
    r2 = classify({1: "q1", 2: "qb1", 3: "g", 4: "g"}, born, ["ss"])
    gluon_soft = [m for m in r2 if m["reason"] != "single quark never soft"]
    assert gluon_soft and all(m["genuine"] for m in gluon_soft), \
        "Born+g: soft gluon must be genuine"

    # 4. mode-count combinatorics for 5 partons: 5/10/10/10/30/15
    from math import comb
    cnt = {}
    for m in res:
        cnt[m["family"]] = cnt.get(m["family"], 0) + 1
    n = 5
    assert cnt["ss"] == comb(n, 1) and cnt["sco"] == comb(n, 2)
    assert cnt["ds"] == comb(n, 2) and cnt["tc"] == comb(n, 3)
    assert cnt["sc"] == n * comb(n - 1, 2)
    assert cnt["dc"] == comb(n, 2) * comb(n - 2, 2) // 2
    print("genuine_modes selftest OK")


def main():
    if "--selftest" in sys.argv:
        selftest()
        return
    spec = json.load(open(sys.argv[1]))
    # keys are momentum labels: ints in check-program specs, but any
    # string label works (e.g. the l,k,i,j,m of a .map's FN line)
    partons = {int(k) if str(k).isdigit() else k: v
               for k, v in spec["partons"].items()}
    res = classify(partons, spec["born"], spec.get("families"))
    if "--json" in sys.argv:
        json.dump(res, sys.stdout, indent=1)
        return
    ng = sum(1 for m in res if m["genuine"])
    print(f"{len(res)} modes, {ng} genuine / {len(res) - ng} dead")
    for m in res:
        tag = "GENUINE" if m["genuine"] else "dead   "
        print(f"  [{tag}] {m['family']:3s}  {m['name']:22s} {m['reason']}")


if __name__ == "__main__":
    main()
