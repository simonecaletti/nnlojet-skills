#!/usr/bin/env python3
"""Predict the BLOCK SKELETON of an RR subtraction term from colour connection.

Rung 0 of the cost ladder: costs seconds, runs nothing from the repo,
and turns "which blocks does this term need?" from a build-cycle
question into a static one. Its output is a PREDICTION that the
measurement stack (probe-me-ir-structure) and the spike test then
adjudicate — a disagreement between prediction and measurement is
itself the highest-value diagnostic, because it means the block
STRUCTURE is wrong, not just a coefficient.

Usage:
  python predict_blocks.py spec.json                  # table
  python predict_blocks.py spec.json --json
  python predict_blocks.py spec.json --emit-markers   # map_blocks.py skeleton
  python predict_blocks.py spec.json --audit TERM.map # exits 1 on mismatch
  python predict_blocks.py --selftest
  (-o FILE writes instead of stdout, as in map_blocks.py)

Spec — a SUPERSET of genuine_modes.py's spec, so one file can feed both.
CAVEAT: the two read DIFFERENT keys (that script: `partons`, final state
only, plus `born`; this one: `chain`, the full colour ordering including
initial legs, plus `flavours`). The flavour information therefore appears
twice and is NOT cross-validated. Keep them in sync by hand.
  chain     : colour-ordered list of leg labels (ints), e.g. [1,3,4,5,2].
              This is ONE colour ordering. Colour-summed channels: run
              once per ordering and take the union (see the limitations
              at the end of this docstring).
  flavours  : {"<leg>": "g" | "q<tag>" | "qb<tag>"} for EVERY chain leg.
  initial   : optional list of chain legs that are initial-state
              (default []). Fixes the FF/IF/FI/II configuration.
  cyclic    : optional bool (default false). True for a pure-gluon
              (adjoint) string, where the two chain endpoints are also
              colour connected. False for a quark string.
  unresolvable : optional explicit list of legs that admit an unresolved
              limit. If absent it is derived, but the RELIABLE source is
              genuine_modes.py (write-spike-test):

                  python genuine_modes.py spec.json --json > modes.json
                  python predict_blocks.py spec.json --modes modes.json

  partons / born : optional — read by genuine_modes.py, ignored here, so
              one spec file can drive both scripts.

Predicts: which blocks exist and their line counts, every sign, each
antenna's FF/IF/FI/II configuration, the X30 family and coefficient, and
the colour-neighbouring sub-class. Does NOT predict: antenna letters,
slot conventions, which X40, X40 coefficients, or reduced-ME argument
orders — measure those (probe-me-ir-structure).

Two limitations that matter:
  - the prediction is per COLOUR ORDERING; check programs sum over
    orderings, so counts are per-ordering upper bounds. Adjacency is the
    RIGHT criterion for block assignment inside one ordered line and the
    WRONG one for deciding which limits are genuine — opposite uses.
  - the fallback unresolvable-leg rule over-counts; it warns on stderr.
"""
import argparse
import itertools
import json
import re
import sys
from fractions import Fraction

# MUST stay byte-identical to map_blocks.py's pair: that script WRITES the
# format this one AUDITS.
MARKER = re.compile(r"^\s*#\s*block:\s*(\S+)")
LABEL = re.compile(r"\*a(\d+)\b")

# Coefficient with which each tree 3-parton antenna enters its integrated
# dipole in maple/form/common/J21.map — verbatim, FF entries:
#     J21QQFF  = calA30FF          J21QGFF  = 1/2*calD30FF
#     J21hQGFF = 1/2*calE30FF      J21GGFF  = 1/3*calF30FF
#     J21hGGFF = calG30FF
# (the quark-loop / N_F family carries an infix `h`, NOT an `_NF` suffix —
# grep `J21hQGFF`, not `J21_NF`.) The value is 1/(number of distinct antenna
# configurations the parent process contains); an S,a line mirrors it.
X30_COEFF = {
    "A": Fraction(1),      # gamma* -> q qb g  : quark endpoints fixed, 1 config
    "D": Fraction(1, 2),   # chi -> gluino g g : 2 choices of outer gluon
    "E": Fraction(1, 2),   # chi -> gluino q' qb'
    "F": Fraction(1, 3),   # H -> g g g        : 3 cyclic configs
    "G": Fraction(1),      # H -> g q qb
}

# The X40 coefficients are NOT predicted here. maple/form/common/J22.map
# assembles integrated four-parton antennae with process- and colour-level
# dependent factors (FF block, lines 2-57): calA40FF 1, 1/2*calD40FF,
# calE40FF 1, calG40FF 1, 1/2*calH40FF, calB40FF 1, 2*calC40FF,
# 1/2*calAt40FF, 1/2*calEt40FF, 1/2*calGt40FF, 1/4*calF40FF. Read that file
# for the case at hand rather than extrapolating a rule.
X40_COEFF_SOURCE = "maple/form/common/J22.map"

# The discrete set observed across J21.map and J22.map. NOTE it is NOT
# all unit fractions: `2*calC40FF` and `-2/3*calG30FF*calF30FF` are both
# live in J22.map. Use it as a MEMBERSHIP test on fitted coefficients —
# a fit landing outside it is diagnostic — never as a prediction of which
# member applies.
RATIONAL_FAMILY = {Fraction(2), Fraction(1), Fraction(2, 3), Fraction(1, 2),
                   Fraction(1, 3), Fraction(1, 4), Fraction(1, 9)}

BLOCK_SIGN = {"Sa": "+", "Sb1": "+", "Sb2": "-", "Sc": "-", "Sd": "-"}
BLOCK_ORDER = ["Sa", "Sb1", "Sb2", "Sc", "Sd"]


# ---------------- chain algebra ----------------

def species(f):
    """flavour string -> 'g' or 'q' (quark and antiquark both 'q' for
    radiator-pair purposes)."""
    if f == "g":
        return "g"
    if f.startswith("q"):
        return "q"
    raise ValueError(f"bad flavour {f!r} (want g / q<tag> / qb<tag>)")


def check_spec(spec):
    chain = spec.get("chain")
    if not chain or len(chain) < 3:
        raise SystemExit("ERROR: spec needs a 'chain' of at least 3 legs")
    if len(set(chain)) != len(chain):
        raise SystemExit("ERROR: repeated leg in 'chain'")
    fl = spec.get("flavours", {})
    missing = [l for l in chain if str(l) not in fl]
    if missing:
        raise SystemExit(f"ERROR: no flavour given for chain leg(s) {missing}")
    bad = [l for l in spec.get("initial", []) if l not in chain]
    if bad:
        raise SystemExit(f"ERROR: 'initial' leg(s) {bad} are not in 'chain'")
    for l in chain:
        species(fl[str(l)])                       # raises on a bad flavour
    return chain, fl, set(spec.get("initial", [])), bool(spec.get("cyclic"))


def distance(chain, cyclic, a, b):
    """Colour-chain distance. 1 = adjacent = colour connected."""
    n = len(chain)
    d = abs(chain.index(a) - chain.index(b))
    return min(d, n - d) if cyclic else d


def neighbours(chain, cyclic, leg):
    """Chain neighbours of a leg — its candidate hard radiators, in chain
    order. A leg with fewer than two cannot sit in an antenna's middle
    slot."""
    n = len(chain)
    i = chain.index(leg)
    out = []
    if i > 0:
        out.append(chain[i - 1])
    elif cyclic:
        out.append(chain[-1])
    if i < n - 1:
        out.append(chain[i + 1])
    elif cyclic:
        out.append(chain[0])
    return out


def in_chain_order(chain, legs):
    """Order legs by colour-chain position — NOT by numeric label. The
    FF/IF/FI distinction depends on which radiator comes first. Only
    meaningful on a LINEAR chain; for a cyclic one use
    connected_radiators, which walks the chain instead."""
    return sorted(legs, key=chain.index)


def connected_radiators(chain, cyclic, j, k):
    """For an ADJACENT pair, return the two hard radiators in true chain
    order: the leg before the pair and the leg after it. Walking the
    chain (rather than sorting by index) is what makes this correct
    across the cyclic wrap, where index order is meaningless."""
    n = len(chain)
    ij = chain.index(j)
    first, second = (j, k) if chain[(ij + 1) % n] == k else (k, j)
    i1, i2 = chain.index(first), chain.index(second)
    before = chain[(i1 - 1) % n] if (cyclic or i1 > 0) else None
    after = chain[(i2 + 1) % n] if (cyclic or i2 < n - 1) else None
    return [x for x in (before, after) if x is not None]


def cluster_parent(flavs):
    """Net flavour of a collinear cluster: 'g', 'q<tag>'/'qb<tag>', or
    None if the net is not a single parton. Same rule as
    genuine_modes.py (write-spike-test) — kept local to keep this script
    dependency-free; if the two ever disagree, genuine_modes.py wins."""

    net = {}
    for f in flavs:
        if f == "g":
            continue
        if f.startswith("qb"):
            t, s = f[2:], -1
        else:
            t, s = f[1:], +1
        net[t] = net.get(t, 0) + s
    nz = {t: v for t, v in net.items() if v}
    total = sum(abs(v) for v in nz.values())
    if total == 0:
        return "g"
    if total == 1:
        (t, v), = nz.items()
        return f"q{t}" if v > 0 else f"qb{t}"
    return None


def config_of(chain_initial, legs):
    """FF / IF / FI / II from which radiators are initial-state."""
    flags = [l in chain_initial for l in legs]
    if not any(flags):
        return "FF"
    if all(flags):
        return "II"
    return "IF" if flags[0] else "FI"


def family_hint(fl, rad, unres):
    """Antenna letter family from the radiator pair and the unresolved
    content. A HINT, not a slot convention: run the pole scan before
    using an antenna you have not used before."""
    sp = tuple(sorted(species(fl[str(r)]) for r in rad))
    unres_q = bool(unres) and all(species(fl[str(u)]) == "q" for u in unres)
    if sp == ("q", "q"):
        return "A", X30_COEFF["A"]
    if sp == ("g", "q"):
        return ("E", X30_COEFF["E"]) if unres_q else ("D", X30_COEFF["D"])
    return ("G", X30_COEFF["G"]) if unres_q else ("F", X30_COEFF["F"])


# ---------------- prediction ----------------

def unresolvable_legs(spec, chain, cyclic, modes=None):
    """Legs that admit an unresolved limit, plus the provenance of that
    answer. Prefer genuine_modes.py output; the fallback CANNOT see the
    reduced-Born rule and will over-count."""
    if spec.get("unresolvable") is not None:
        return list(spec["unresolvable"]), "spec"
    if modes:
        live = set()
        for m in modes:
            if not m.get("genuine"):
                continue
            txt = str(m.get("name", ""))
            for ch in "|+,":
                txt = txt.replace(ch, " ")
            for tok in txt.split():
                if tok.isdigit() and int(tok) in chain:
                    live.add(int(tok))
        if live:
            return sorted(live), "genuine_modes"
    return [l for l in chain if len(neighbours(chain, cyclic, l)) == 2], \
        "fallback"


def predict(spec, modes=None):
    chain, fl, initial, cyclic = check_spec(spec)
    raw, src = unresolvable_legs(spec, chain, cyclic, modes)
    warnings = []
    if src == "fallback":
        warnings.append(
            "unresolvable legs came from the FALLBACK rule: it cannot see "
            "the reduced-Born rule and OVER-counts. Feed genuine_modes.py "
            "output via --modes before trusting any line count.")
    legs, dropped = [], []
    for l in raw:
        if l not in chain:
            dropped.append((l, "not in chain"))
        elif len(neighbours(chain, cyclic, l)) != 2:
            dropped.append((l, "chain endpoint — cannot sit in an antenna's "
                               "middle slot; if it should be unresolved, the "
                               "chain or `cyclic` is wrong"))
        else:
            legs.append(l)
    for l, why in dropped:
        warnings.append(f"leg {l} dropped from the unresolvable list: {why}")

    lines, pairs = [], []

    # --- S,a : one X30 per single-unresolved leg -------------------------
    for j in legs:
        rad = neighbours(chain, cyclic, j)
        fam, coeff = family_hint(fl, rad, [j])
        lines.append({
            "block": "Sa", "sign": "+", "unresolved": [j], "radiators": rad,
            "config": config_of(initial, rad), "family_hint": fam + "30",
            "coeff_hint": str(coeff),
            "structure": "X30 * M_(n-1) * JET",
            "note": "the NLO term with the jet function demoted"})

    # --- pairs -----------------------------------------------------------
    for j, k in itertools.combinations(sorted(legs), 2):
        d = distance(chain, cyclic, j, k)
        cls = ("colour-connected" if d == 1 else
               "almost-colour-connected" if d == 2 else "colour-unconnected")
        nb_j, nb_k = set(neighbours(chain, cyclic, j)), \
            set(neighbours(chain, cyclic, k))
        shared = in_chain_order(chain, nb_j & nb_k)
        outer = connected_radiators(chain, cyclic, j, k) if d == 1 else []
        # Colour-neighbouring is a property of a PAIR OF CLUSTERS, not of a
        # single unresolved pair: the two partons cluster OUTWARD onto their
        # respective radiators, giving a double collinear limit of two
        # chain-adjacent clusters. Only then do the Sb2 iterated products
        # fail to vanish — each equals the double unresolved limit and
        # cancels Sa's DOUBLED term as well as X40's spurious single, which
        # is the second job pole-sharing does not see
        # (references/block-counting.md). It requires the pair to be
        # adjacent, to have two DISTINCT outer radiators, and BOTH outward
        # clusters to be flavour-valid.
        nconf = None
        if d == 1 and len(outer) == 2 and outer[0] != outer[1]:
            ij = chain.index(j)
            first, second = ((j, k) if chain[(ij + 1) % len(chain)] == k
                             else (k, j))
            c1 = (outer[0], first)
            c2 = (second, outer[1])
            if all(cluster_parent([fl[str(a)], fl[str(b)]]) is not None
                   for a, b in (c1, c2)):
                nconf = [list(c1), list(c2)]
        pairs.append({"pair": [j, k], "distance": d, "class": cls,
                      "shared_radiator": shared,
                      "colour_neighbouring": nconf is not None,
                      "neighbouring_clusters": nconf})

        if d == 1:
            rad = [outer[0], outer[-1]]
            fam, _ = family_hint(fl, rad, [j, k])
            lines.append({
                "block": "Sb1", "sign": "+", "unresolved": [j, k],
                "radiators": rad, "config": config_of(initial, rad),
                "family_hint": fam + "40",
                "coeff_hint": f"not predicted — read {X40_COEFF_SOURCE}",
                "structure": "X40 * M_(n-2) * JET",
                "note": "genuine double-unresolved; choose the X40 by "
                        "MEASURED residue, never by colour adjacency"})
            for u in (j, k):
                lines.append({
                    "block": "Sb2", "sign": "-", "unresolved": [j, k],
                    "radiators": rad, "config": config_of(initial, rad),
                    "family_hint": "X30*X30",
                    "coeff_hint": "1, or 1/2 iff the two clusters are "
                                  "symmetric",
                    "structure": "- X30 * X30 * M_(n-2) * JET",
                    "note": f"removes the X40's spurious {u}-unresolved pole; "
                            "second antenna evaluated on the first's MAPPED "
                            "momenta (the two writings are NOT equivalent)"
                            + (f" | COLOUR-NEIGHBOURING via {nconf}: this "
                               "line ALSO cancels Sa's doubled double limit"
                               if nconf else "")})
        elif d == 2:
            for first, second in ((j, k), (k, j)):
                lines.append({
                    "block": "Sc", "sign": "-", "unresolved": [first, second],
                    "radiators": shared,
                    "config": config_of(initial, shared or [first]),
                    "family_hint": "SS-difference * X30",
                    "coeff_hint": "1/2 per sub-antenna (x30 yields HALF the "
                                  "soft eikonal)",
                    "structure": "- (sum of +-SS(..)) * X30 * M * JET",
                    "note": f"strong ordering {first} then {second}; use a "
                            "SUB-antenna, not a Full composite — a mapped "
                            "momentum must not be a collinear partner"})
        else:
            lines.append({
                "block": "Sd", "sign": "-", "unresolved": [j, k],
                "radiators": in_chain_order(chain, nb_j | nb_k),
                "config": "mixed",
                "family_hint": "X30*X30 (disjoint)",
                "coeff_hint": "product of the two X30 coefficients",
                "structure": "- X30 * X30 * M_(n-2) * JET",
                "note": "disjoint clusters; UNORDERED pairs only — no product "
                        "of two antenna configurations may appear twice"})

    counts = {}
    for ln in lines:
        counts[ln["block"]] = counts.get(ln["block"], 0) + 1
    return {"chain": chain, "cyclic": cyclic, "initial": sorted(initial),
            "unresolvable": legs, "unresolvable_source": src,
            "warnings": warnings,
            "pairs": pairs, "lines": lines, "counts": counts}


# ---------------- outputs ----------------

def render(pred):
    out = [f"chain {pred['chain']}  cyclic={pred['cyclic']}  "
           f"initial={pred['initial']}",
           f"unresolvable legs {pred['unresolvable']}  "
           f"(source: {pred['unresolvable_source']})"]
    for w in pred["warnings"]:
        out.append(f"  !! {w}")
    out += ["", "pair classification"]
    for p in pred["pairs"]:
        tag = (f"  [COLOUR-NEIGHBOURING {p['neighbouring_clusters']}]"
               if p["colour_neighbouring"] else "")
        out.append(f"  {str(p['pair']):10s} d={p['distance']}  "
                   f"{p['class']:24s}{tag}")
    out += ["", "predicted blocks"]
    for b in BLOCK_ORDER:
        if b in pred["counts"]:
            out.append(f"  {b:4s} {BLOCK_SIGN[b]}  {pred['counts'][b]:3d} "
                       "line(s)")
    out += ["", "predicted lines"]
    for ln in pred["lines"]:
        out.append(f"  {ln['block']:4s} {ln['sign']} unres={ln['unresolved']} "
                   f"rad={ln['radiators']} {ln['config']:5s} "
                   f"{ln['family_hint']}  coeff~{ln['coeff_hint']}")
        out.append(f"        {ln['structure']}")
        out.append(f"        {ln['note']}")
    out += ["", "NEXT: this is a PREDICTION. Confirm each line by measurement "
            "(probe-me-ir-structure) before composing. A disagreement means "
            "the block STRUCTURE is wrong, not the coefficient."]
    return "\n".join(out)


def emit_markers(pred, fn="FN_PLACEHOLDER"):
    """A map_blocks.py-composable skeleton: one placeholder term per
    predicted line, grouped under '# block:' markers, aN gap-free. The
    label is placed BEFORE the trailing comment so it survives as a
    label, not as commented-out text."""
    body, i = [], 0
    for b in BLOCK_ORDER:
        sel = [l for l in pred["lines"] if l["block"] == b]
        if not sel:
            continue
        body.append(f"# block: {b}")
        for ln in sel:
            i += 1
            stem = (ln["family_hint"].replace("*", "x").replace(" ", "")
                    .replace("-", "_"))
            tag = "_".join(map(str, ln["unresolved"]))
            body.append(f"{ln['sign']}TODO_{stem}_{tag}*a{i}"
                        f"    # rad={ln['radiators']} {ln['config']} "
                        f"coeff~{ln['coeff_hint']}")
    blocks = ",".join(b for b in BLOCK_ORDER if b in pred["counts"])
    return "\n".join(
        ["# PREDICTED SKELETON — placeholders, NOT a valid .map.",
         "# Replace each TODO_* with antenna(args)*reducedME(args)*JETnm(args),",
         f"# then: map_blocks.py compose <this> --blocks {blocks} -o TERM.map"]
        + [f"# !! {w}" for w in pred["warnings"]]
        + [f"FN:={fn}:", "XX:="] + body + [":"]) + "\n"


def audit(pred, text):
    """Compare a written .map against the prediction. STRUCTURAL ONLY: it
    checks which blocks exist and how many terms each carries. It cannot
    check antenna letters, arguments or mappings — those are measured."""
    actual, cur = {}, None
    for ln in text.splitlines():
        m = MARKER.match(ln)
        if m:
            cur = m.group(1)
            actual.setdefault(cur, 0)
        elif cur is not None:
            actual[cur] += len(LABEL.findall(ln))
    if not actual:
        return ["no '# block:' markers found — the file is not composed from "
                "a block master, so the audit cannot run."]
    stems = {k.split("_")[0] for k in actual}
    rep = []
    for b in BLOCK_ORDER:
        if b not in pred["counts"]:
            continue
        n = pred["counts"][b]
        if b not in stems:
            rep.append(f"MISSING  {b}: predicted {n} line(s), block absent")
            continue
        got = sum(v for k, v in actual.items() if k.split("_")[0] == b)
        if got != n:
            rep.append(f"COUNT    {b}: predicted {n}, file has {got} "
                       "(per-ordering counts are upper bounds — a colour-"
                       "summed term may legitimately differ)")
        else:
            rep.append(f"ok       {b}: {n}")
    for k in sorted(actual):
        if k != "_default" and k.split("_")[0] not in pred["counts"]:
            rep.append(f"EXTRA    {k}: present but not predicted — either a "
                       "second colour ordering, or a spurious block")
    return rep


# ---------------- selftest ----------------

def selftest():
    """Checks ALGORITHM invariants only. Encodes no real process's block
    answer, and no coefficient beyond the table's internal consistency."""
    spec = {"chain": [1, 3, 4, 5, 2],
            "flavours": {"1": "q1", "3": "g", "4": "g", "5": "g", "2": "qb1"},
            "initial": [1, 2],
            "unresolvable": [3, 4, 5]}
    p = predict(spec)

    # 1. every unresolved pair gets exactly one class, monotone in distance
    assert len({tuple(x["pair"]) for x in p["pairs"]}) == 3
    for x in p["pairs"]:
        d, c = x["distance"], x["class"]
        assert (d == 1) == (c == "colour-connected")
        assert (d == 2) == (c == "almost-colour-connected")
        assert (d >= 3) == (c == "colour-unconnected")

    # 2. block bookkeeping is fixed by the counting argument, not by content
    ncc = sum(1 for x in p["pairs"] if x["distance"] == 1)
    nac = sum(1 for x in p["pairs"] if x["distance"] == 2)
    nun = sum(1 for x in p["pairs"] if x["distance"] >= 3)
    c = p["counts"]
    assert c.get("Sa", 0) == len(p["unresolvable"]), c
    assert c.get("Sb1", 0) == ncc, c
    assert c.get("Sb2", 0) == 2 * ncc, c      # one per single limit of the X40
    assert c.get("Sc", 0) == 2 * nac, c       # both strong orderings
    assert c.get("Sd", 0) == nun, c           # unordered pairs, no duplicates

    # 3. signs come from the block, never from the content
    for ln in p["lines"]:
        assert ln["sign"] == BLOCK_SIGN[ln["block"]], ln

    # 4. chain-reversal invariance (an ordering may be read either way)
    pr = predict(dict(spec, chain=list(reversed(spec["chain"]))))
    assert pr["counts"] == p["counts"], (pr["counts"], p["counts"])
    assert {tuple(sorted(x["pair"])): x["class"] for x in pr["pairs"]} == \
           {tuple(sorted(x["pair"])): x["class"] for x in p["pairs"]}

    # 5. leg-relabelling invariance
    ren = {1: 11, 2: 12, 3: 13, 4: 14, 5: 15}
    p2 = predict({"chain": [ren[x] for x in spec["chain"]],
                  "flavours": {str(ren[int(k)]): v
                               for k, v in spec["flavours"].items()},
                  "initial": [ren[x] for x in spec["initial"]],
                  "unresolvable": [ren[x] for x in spec["unresolvable"]]})
    assert p2["counts"] == p["counts"]
    assert [x["class"] for x in p2["pairs"]] == [x["class"] for x in p["pairs"]]

    # 6. cyclic closure: endpoints adjacent iff cyclic
    gchain = [1, 2, 3, 4, 5]
    assert distance(gchain, True, 1, 5) == 1
    assert distance(gchain, False, 1, 5) == 4
    gl = predict({"chain": gchain, "cyclic": True,
                  "flavours": {str(i): "g" for i in gchain},
                  "unresolvable": [1, 5]})
    assert gl["pairs"][0]["distance"] == 1, "cyclic wrap failed"
    assert gl["counts"].get("Sb1") == 1, gl["counts"]
    # radiators must be read by WALKING the chain: across the wrap the pair
    # reads 4-5-1-2, so the radiators are [4,2] and NOT index-sorted [2,4]
    assert connected_radiators(gchain, True, 1, 5) == [4, 2], \
        connected_radiators(gchain, True, 1, 5)
    assert connected_radiators(gchain, False, 2, 3) == [1, 4]
    sb1 = [l for l in gl["lines"] if l["block"] == "Sb1"][0]
    assert sb1["radiators"] == [4, 2], sb1["radiators"]
    # the same legs are not even unresolvable on a linear chain (endpoints
    # have a single neighbour) — and the drop must be REPORTED, not silent
    lin = predict({"chain": gchain, "cyclic": False,
                   "flavours": {str(i): "g" for i in gchain},
                   "unresolvable": [1, 5]})
    assert lin["unresolvable"] == [] and lin["pairs"] == [], lin
    assert len(lin["warnings"]) == 2, lin["warnings"]
    assert all("dropped" in w for w in lin["warnings"]), lin["warnings"]

    # 7. configuration follows CHAIN position, not numeric label
    cfg = {tuple(l["unresolved"]): l["config"]
           for l in p["lines"] if l["block"] == "Sa"}
    assert cfg[(3,)] == "IF", cfg      # radiators 1(initial), 4 -> I then F
    assert cfg[(4,)] == "FF", cfg      # radiators 3, 5
    assert cfg[(5,)] == "FI", cfg      # radiators 4, 2(initial) -> F then I

    # 8. emitted skeleton is gap-free and audits clean against itself
    skel = emit_markers(p)
    labels = [int(x) for x in LABEL.findall(skel)]
    assert labels == list(range(1, len(p["lines"]) + 1)), labels
    # the label must not be swallowed by the trailing comment
    for ln in skel.splitlines():
        if "*a" in ln and "#" in ln:
            assert ln.index("*a") < ln.index("#"), ln
    assert all(r.startswith("ok") for r in audit(p, skel)), audit(p, skel)

    # 9. the audit detects under-count, missing blocks and extra blocks
    rep = audit(p, "# block: Sa\n+X*a1\n")
    assert any(r.startswith("COUNT") and " Sa:" in r for r in rep), rep
    assert any(r.startswith("MISSING") for r in rep), rep
    rep = audit(p, "# block: Sa\n+X*a1\n# block: Zz\n+Y*a2\n")
    assert any(r.startswith("EXTRA") for r in rep), rep
    assert audit(p, "no markers here")[0].startswith("no '# block:'")

    # 10. malformed specs are refused, never guessed
    for bad in ({"chain": [1, 2], "flavours": {"1": "g", "2": "g"}},
                {"chain": [1, 2, 3], "flavours": {"1": "g", "2": "g"}},
                {"chain": [1, 2, 3],
                 "flavours": {"1": "g", "2": "g", "3": "x"}},
                {"chain": [1, 2, 2],
                 "flavours": {"1": "g", "2": "g"}},
                {"chain": [1, 2, 3], "initial": [9],
                 "flavours": {"1": "g", "2": "g", "3": "g"}}):
        try:
            predict(bad)
            raise AssertionError(f"malformed spec accepted: {bad}")
        except (SystemExit, ValueError):
            pass

    # 11. the X30 table is a subset of the observed rational family, and the
    #     family is NOT all unit fractions (2 and 2/3 are live in J22.map —
    #     an assertion demanding numerator==1 would be wrong, not strict)
    assert set(X30_COEFF.values()) <= RATIONAL_FAMILY, X30_COEFF
    assert Fraction(2) in RATIONAL_FAMILY and \
        Fraction(2, 3) in RATIONAL_FAMILY, RATIONAL_FAMILY
    # no X40 coefficient is claimed anywhere in the output
    for ln in p["lines"]:
        if ln["block"] == "Sb1":
            assert "not predicted" in ln["coeff_hint"], ln

    # 11b. colour-neighbouring is a PAIR-OF-CLUSTERS property, so it must be
    #      falsifiable: with the outward clusters flavour-invalid it is off
    cn = {tuple(x["pair"]): x["colour_neighbouring"] for x in p["pairs"]}
    assert cn[(3, 4)] and cn[(4, 5)], cn        # q-g and g-qb clusters valid
    dead = predict({"chain": [1, 3, 4, 2],
                    "flavours": {"1": "q1", "3": "g", "4": "q2", "2": "q3"},
                    "unresolvable": [3, 4]})
    # the outward cluster (4,2) is q2 with q3 — net TWO quarks, not a single
    # parton, so no colour-neighbouring configuration exists here even
    # though the pair (3,4) is colour connected
    assert dead["pairs"][0]["distance"] == 1
    assert dead["pairs"][0]["colour_neighbouring"] is False, dead["pairs"][0]
    assert dead["pairs"][0]["neighbouring_clusters"] is None

    # 12. --modes wiring: genuine modes drive the leg list, dead ones do not
    modes = [{"family": "ss", "name": "4 soft", "genuine": True},
             {"family": "ss", "name": "3 soft", "genuine": False},
             {"family": "sco", "name": "4||5", "genuine": True}]
    pm = predict({k: v for k, v in spec.items() if k != "unresolvable"}, modes)
    assert pm["unresolvable_source"] == "genuine_modes", pm
    assert pm["unresolvable"] == [4, 5], pm["unresolvable"]

    print("predict_blocks selftest OK")


def main():
    if "--selftest" in sys.argv:
        selftest()
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("spec")
    ap.add_argument("--modes", help="genuine_modes.py --json output")
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--emit-markers", action="store_true")
    ap.add_argument("--audit", metavar="TERM.map")
    ap.add_argument("-o", "--output")
    args = ap.parse_args()
    spec = json.load(open(args.spec))
    modes = json.load(open(args.modes)) if args.modes else None
    pred = predict(spec, modes)

    # warnings go to stderr on EVERY path — a prediction built from the
    # fallback leg list must not look authoritative just because the caller
    # asked for --json or --audit
    for w in pred["warnings"]:
        print(f"WARNING: {w}", file=sys.stderr)

    if args.audit:
        rep = audit(pred, open(args.audit).read())
        for r in rep:
            print(r)
        # a gate that cannot fail is not a gate
        raise SystemExit(0 if all(r.startswith("ok") for r in rep) else 1)

    if args.json:
        out = json.dumps(pred, indent=1) + "\n"
    elif args.emit_markers:
        out = emit_markers(pred)
    else:
        out = render(pred) + "\n"
    if args.output:
        with open(args.output, "w") as fh:
            fh.write(out)
        print(f"wrote {args.output}: {len(pred['lines'])} predicted line(s)")
    else:
        sys.stdout.write(out)


if __name__ == "__main__":
    main()
