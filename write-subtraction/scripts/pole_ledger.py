#!/usr/bin/env python3
"""Pole ledger: static bookkeeping check of an RR subtraction term
against the antenna datasheet — the missing rung between
audit_blocks.py (counts blocks) and the spike test (costs a
regenerate+rebuild+run per hypothesis).

Given the .map, a flavour spec and the measured antenna datasheet
(antennae-naming-convention/scripts/antenna_datasheet.py), it checks —
in seconds, with no build:

  1. STRUCTURE  every antenna argument is a momentum the line has at
     that point (original leaves, or clusters produced by an earlier
     antenna of the same line); the reduced-ME and JET arguments are
     exactly the momentum set the line's mappings produce (typo and
     stale-cluster detector).
  2. FLAVOUR    slot species (quark-kind vs gluon) and same-flavour
     pair constraints per antenna; cluster flavours propagated with the
     datasheet's reduced-cluster semantics (NOT naive net-flavour
     arithmetic, which is wrong for the E/G families).
  3. SPLIT      a Full composite whose measured pole graph the generic
     cluster rule cannot support (requires_split, e.g. FullE40's
     endpoint quark||gluon pole) used on all-final clusters is an
     ERROR; the verified split identity is printed.
  4. COVERAGE   every genuine single (ss/sco) mode is covered by an
     S,a line singular there; every genuine double (ds/tc) mode by an
     X40 (or SS) line singular there.  sc/dc modes are only checkable
     when a single line carries both invariants — otherwise reported
     as unverified (their coverage also lives in reduced-ME poles,
     which a static check cannot see).
  5. PAIRING    every flavour-valid single-unresolved pole of every
     X40 line (its SPURIOUS singles) must be matched by at least one
     negative iterated (X30 x X30) line singular in the same
     invariant; orphaned iterated lines and sign-rule violations are
     reported.
  6. SOFT DIPOLES  every soft leg of every antenna emits an eikonal on
     a definite dipole; that dipole must be a COLOUR-CONNECTED pair of
     one of the declared chains.  A dipole present in no chain is a
     singularity the matrix element does not have, so no counterterm
     anywhere can cancel it -- the term is wrong however the rest of
     the bookkeeping looks.  This is the check that catches a wrong
     X40 SLOT ASSIGNMENT, which is otherwise invisible until the
     spike test: the empirical signature is the single-soft mode
     sitting at a stable fraction (0.3-0.9, wide spread) while every
     collinear mode reads 1.000000, and it costs a build cycle per
     hypothesis to localise.  Requires "chains" in the spec; without
     it the check is reported UNCHECKED, never silently skipped.
     Legs that are mapped momenta are skipped (the reduced chain is
     not tracked) and listed in the info block.

     The dipole is taken from the datasheet when recorded and
     otherwise DERIVED as "the soft leg's two collinear partners in
     the antenna's own pole graph" -- s(A,S) and s(S,B) singular =>
     eikonal on (A,B).  The selftest asserts that derivation against
     every dipole the datasheet does record, which is what licenses
     using it for the X40 entries whose dipole field is still empty.

Usage:
  pole_ledger.py TERM.map --spec spec.json [--datasheet FILE]
                 [--modes modes.json] [--verbose]
  pole_ledger.py --selftest

spec.json:
  { "flavours": {"l":"qb1","k":"g","i":"q2","j":"qb2","m":"q1"},
    "born": [["q1","qb1","g"]],
    "chains": [["l","m","i","j","k"], ["l","i","j","m","k"]] }
  Momentum labels are the FN argument names of the .map (letters or
  numbers); args of FN not in "flavours" are non-partons (leptons).
  "chains" is OPTIONAL but strongly recommended: one entry per colour
  ordering the matrix element sums over, each a list of parton labels
  in colour order.  For a colour-summed ME list them all -- the check
  accepts a dipole adjacent in ANY declared chain, so a missing chain
  produces false errors and a spurious one hides real ones.  For a
  4-quark ME remember that BOTH readings of the colour structure are
  chains (the pair that plays "secondary" swaps between them).

Exit status: 1 if any ERROR, else 0.  Warnings and unverified modes do
not fail the check.

SCOPE — existence, not sufficiency (measured on mutation tests): the
ledger detects MISSING structure (an invariant nothing touches, an
uncovered genuine mode, an unsplittable composite, a stale cluster) in
seconds without a build.  It cannot detect INSUFFICIENT structure — a
counterterm block whose lines touch the right invariants with the
wrong coefficients or the wrong mapping still passes (removing
alternative.map's A30xE30 block goes unflagged because the E30xd30
lines also touch s(l,k)).  Coefficient-level completeness remains the
spike test's job; run the ledger FIRST, the spike test after.
"""
import argparse
import itertools
import json
import os
import re
import subprocess
import sys
import tempfile

from _colour_algebra import fparse, kind_class, make_valid_cluster
from _map_parser import norm, leaves, show, parse_map

SKILLS = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
DEFAULT_DATASHEET = os.path.join(
    SKILLS, "antennae-naming-convention", "assets",
    "antenna_datasheet.json")
GENUINE_MODES = os.path.join(
    SKILLS, "write-spike-test", "scripts", "genuine_modes.py")


# Momentum terms (norm/leaves/show) and the .map term parser live in
# _map_parser.py; flavour algebra (fparse/kind_class/make_valid_cluster)
# in _colour_algebra.py — both shared libraries in this directory.

# ---------------- genuine modes ----------------

def genuine_modes(spec, modes_path=None):
    if modes_path:
        return json.load(open(modes_path))
    with tempfile.NamedTemporaryFile("w", suffix=".json",
                                     delete=False) as f:
        json.dump({"partons": spec["flavours"],
                   "born": spec["born"],
                   "families": spec.get("families")}, f)
        tmp = f.name
    try:
        out = subprocess.run(
            [sys.executable, GENUINE_MODES, tmp, "--json"],
            capture_output=True, text=True, check=True)
        return json.loads(out.stdout)
    finally:
        os.unlink(tmp)


def parse_mode_name(family, name):
    """mode name -> frozenset structure for matching."""
    if family == "ss":
        return frozenset([name.replace(" soft", "").strip()])
    if family == "sco":
        return frozenset(v.strip() for v in name.split("||"))
    if family == "ds":
        return frozenset(v.strip() for v in
                         name.replace(" soft", "").split(","))
    if family == "tc":
        return frozenset(v.strip() for v in name.split("||"))
    if family == "sc":
        soft, coll = name.split("+")
        return (frozenset([soft.replace(" soft", "").strip()]),
                frozenset(v.strip() for v in coll.split("||")))
    if family == "dc":
        a, b = name.split("+")
        return frozenset([frozenset(v.strip() for v in a.split("||")),
                          frozenset(v.strip() for v in b.split("||"))])
    return None


# ---------------- line analysis ----------------

class LineCheck:
    def __init__(self, aN, term, flavours, sheet, report,
                 subleading=False):
        self.aN = aN
        self.term = term
        self.subleading = subleading      # spec: colour == 'subleading'
        self.flavours = flavours          # leaf -> flavour string
        self.sheet = sheet
        self.rep = report
        self.antennae = []                # (token, args, entry)
        self.jet = None
        self.me = None
        for name, args in term["factors"]:
            if re.match(r"JET\d\d$", name):
                self.jet = (name, args)
            elif name in sheet:
                self.antennae.append((name, args, sheet[name]))
            else:
                self.me = (name, args)
        self.cluster_flav = {}            # norm(cluster) -> flavour
        self.produced = []                # clusters produced, in order
        self.final_set = None
        self.ok_structure = self.run_mapping()

    # -- mapping / structure --

    def err(self, msg):
        self.rep["errors"].append(f"a{self.aN}: {msg}")

    def warn(self, msg):
        self.rep["warnings"].append(f"a{self.aN}: {msg}")

    def flav_of(self, x):
        n = norm(x)
        if isinstance(n, str):
            return self.flavours.get(n)      # None for leptons
        return self.cluster_flav.get(n)

    def run_mapping(self):
        partons = [p for p in self.flavours]
        current = {norm(p) for p in partons}
        ok = True
        for token, args, entry in self.antennae:
            nslot = entry["arity"]
            if len(args) != nslot:
                self.err(f"{token} called with {len(args)} args, "
                         f"arity {nslot}")
                return False
            nargs = [norm(a) for a in args]
            for a in nargs:
                if a not in current:
                    self.err(f"{token} argument {show(a)} is not an "
                             f"available momentum at this point "
                             f"(stale or mistyped cluster)")
                    ok = False
            # generic cluster rule
            if nslot == 3:
                cl = [(args[0], args[1]), (args[1], args[2])]
            else:
                cl = [(args[0], args[1], args[2]),
                      (args[3], args[2], args[1])]
            red = entry.get("reduced")
            for ic, c in enumerate(cl):
                nc = norm(c)
                if red:
                    cls, slot = red[ic]
                    if cls == "g":
                        self.cluster_flav[nc] = "g"
                    else:
                        src = self.flav_of(args[slot - 1])
                        self.cluster_flav[nc] = src
                else:
                    self.cluster_flav[nc] = None
            current -= set(nargs)
            current |= {norm(c) for c in cl}
            self.produced.append((token, [norm(c) for c in cl]))
        self.final_set = current
        # reduced-ME / JET argument consistency
        for what, fac in (("reduced ME", self.me), ("JET", self.jet)):
            if not fac:
                continue
            name, args = fac
            got = {norm(a) for a in args
                   if norm(a) in current
                   or not isinstance(norm(a), str)
                   or norm(a) in self.flavours}
            pargs = [norm(a) for a in args
                     if not (isinstance(norm(a), str)
                             and norm(a) not in self.flavours)]
            if set(pargs) != current:
                miss = [show(c) for c in sorted(current - set(pargs),
                                                key=repr)]
                extra = [show(c) for c in sorted(set(pargs) - current,
                                                 key=repr)]
                self.err(f"{what} {name} momentum set mismatch: "
                         f"missing {miss or '-'} / not produced by the "
                         f"mappings {extra or '-'}")
                ok = False
        return ok

    # -- soft eikonal dipoles vs the colour chain --

    def soft_dipoles(self):
        """[(token, soft_arg, (legA, legB), how)] for every soft leg of
        every antenna on this line, expressed in .map arguments.

        The dipole is read from the datasheet when it is recorded, and
        otherwise DERIVED: a soft gluon's eikonal sits on the dipole
        formed by its two COLLINEAR partners in the antenna's own pole
        graph (s(A,S) and s(S,B) singular => eikonal on (A,B)).  That
        derivation reproduces every dipole the datasheet does record
        (asserted in the selftest), which is what makes it usable for
        the X40 entries whose dipole field is still empty.
        """
        out = []
        for token, args, entry in self.antennae:
            meas = entry.get("measured") or {}
            sco = set()
            for k in (meas.get("sco") or {}):
                try:
                    p, q = (int(t) for t in k.split(","))
                except ValueError:
                    continue
                sco.add(frozenset((p, q)))
            for b, d in sorted((meas.get("ss") or {}).items()):
                try:
                    ib = int(b)
                except ValueError:
                    continue
                rec = [k for k in (d.get("dipoles") or {})]
                if rec:
                    for k in rec:
                        try:
                            p, q = (int(t) for t in k.split(","))
                        except ValueError:
                            continue
                        out.append((token, args[ib - 1],
                                    (args[p - 1], args[q - 1]),
                                    "datasheet"))
                    continue
                partners = sorted(next(iter(pr - {ib})) for pr in sco
                                  if ib in pr)
                if len(partners) != 2:
                    continue        # not a plain eikonal: leave alone
                out.append((token, args[ib - 1],
                            (args[partners[0] - 1], args[partners[1] - 1]),
                            "derived"))
        return out

    def check_soft_dipoles(self, adjacency, unchecked):
        """Every soft-gluon eikonal a line emits must sit on a dipole the
        matrix element actually has, i.e. a colour-connected pair of one
        of the declared chains.  A dipole that exists in no chain is a
        singularity nothing can cancel: the classic symptom is the
        single-soft mode reading a stable fraction (~0.4) while every
        collinear mode is exact.

        Legs that are mapped momenta are skipped (the reduced chain is
        not tracked here) and counted in `unchecked`.
        """
        for token, soft, (a, b) in ((t, s, d)
                                    for t, s, d, _ in self.soft_dipoles()):
            legs = [norm(a), norm(b)]
            if any(not isinstance(x, str) or x not in self.flavours
                   for x in legs + [norm(soft)]):
                unchecked.append((self.aN, token))
                continue
            if frozenset(legs) not in adjacency:
                # Chain adjacency is a LEADING-COLOUR criterion. A
                # subleading-colour (tilde) ME is a different colour
                # projection of the same amplitudes and genuinely carries
                # soft correlations on non-adjacent (non-planar) dipoles —
                # measurable with probe-me-ir-structure mode 1. Declaring
                # "colour": "subleading" in the spec downgrades this to a
                # warning so such a term can be checked at all; it stays a
                # hard error by default, where it catches real slot errors.
                msg = (
                    f"{token}: soft {show(norm(soft))} puts its eikonal on "
                    f"the dipole ({show(legs[0])},{show(legs[1])}), which "
                    f"is NOT colour-connected in any declared chain. At "
                    f"LEADING colour the matrix element has no such "
                    f"singularity, so nothing can cancel it — check the "
                    f"antenna's slot assignment (permute the ARGUMENTS, "
                    f"never wrap the function).")
                if self.subleading:
                    self.warn(msg + " [spec says colour=subleading: "
                              "non-planar dipoles are expected here — "
                              "CONFIRM this one by a dipole fit before "
                              "trusting it]")
                else:
                    self.err(msg)

    # -- flavour / species --

    def check_species(self):
        for token, args, entry in self.antennae:
            species = entry.get("species")
            if species:
                for k, (sp, a) in enumerate(zip(species, args), 1):
                    fl = self.flav_of(a)
                    if fl is None:
                        continue
                    want = "g" if sp == "g" else "Q"
                    if kind_class(fl) != want:
                        self.err(
                            f"{token} slot {k} wants {want}, got "
                            f"{show(a)} = {fl}")
            for (p, q) in entry.get("slot_pairs", []):
                fa = self.flav_of(args[p - 1])
                fb = self.flav_of(args[q - 1])
                if fa is None or fb is None:
                    continue
                ka, ta = fparse(fa)
                kb, tb = fparse(fb)
                if ta != tb or {ka, kb} != {"q", "qb"}:
                    self.err(
                        f"{token} slots ({p},{q}) must hold a "
                        f"same-flavour quark/antiquark pair, got "
                        f"{fa},{fb}")

    def check_split(self):
        for token, args, entry in self.antennae:
            if not entry.get("requires_split"):
                continue
            allfinal = all(isinstance(norm(a), str) or True
                           for a in args)   # all-final process
            halves = entry.get("halves") or [
                t for t in self.sheet
                if t.startswith(token) and t != token]
            if entry.get("halves") or token.startswith("Full") \
                    or entry.get("fortran", "").startswith("Full"):
                ident = entry.get("split_identity")
                self.err(
                    f"{token} is a Full composite whose pole graph the "
                    f"generic cluster rule cannot support "
                    f"(unsupported poles at slots "
                    f"{entry.get('unsupported_poles')}): it MUST be "
                    f"split into halves with their own mappings"
                    + (f" — verified identity: Full = {ident} "
                       f"(halves {halves})" if ident else
                       f" — registered sub-antennae: {halves or '?'}"))
            else:
                self.warn(
                    f"{token} carries a pole the generic clusters do "
                    f"not support (slots "
                    f"{entry.get('unsupported_poles')}); its "
                    f"cancellation must come from another line's "
                    f"mapping — verify the pairing below")

    # -- pole instantiation --

    def primary(self):
        """The line's antennae acting on ORIGINAL leaves only."""
        return [(t, a, e) for (t, a, e) in self.antennae
                if all(isinstance(norm(x), str) for x in a)]

    def poles(self):
        """All (kind, frozenset-of-leaf-labels) invariants any antenna
        of this line is singular in, restricted to leaf-only slots.
        kind in {'sco','ss','tc','ds'}."""
        out = set()
        for token, args, entry in self.antennae:
            meas = entry.get("measured", {})
            for key in meas.get("sco", {}):
                p, q = (int(v) for v in key.split(","))
                a, b = norm(args[p - 1]), norm(args[q - 1])
                if isinstance(a, str) and isinstance(b, str):
                    out.add(("sco", frozenset([a, b])))
            for bslot, d in meas.get("ss", {}).items():
                a = norm(args[int(bslot) - 1])
                if isinstance(a, str):
                    out.add(("ss", frozenset([a])))
                elif len(leaves(a)) == 2:
                    # a soft-capable slot holding a 2-leaf cluster is a
                    # double-soft of the pair in original momenta
                    out.add(("ds", frozenset(leaves(a))))
            for key in meas.get("tc", {}):
                ids = [norm(args[int(v) - 1]) for v in key.split(",")]
                if all(isinstance(v, str) for v in ids):
                    out.add(("tc", frozenset(ids)))
            for key in meas.get("ds", {}):
                p, q = (int(v) for v in key.split(","))
                a, b = norm(args[p - 1]), norm(args[q - 1])
                if isinstance(a, str) and isinstance(b, str):
                    out.add(("ds", frozenset([a, b])))
        return out


# ---------------- the ledger ----------------

def run_ledger(mapfile, spec, sheet, modes, verbose=False, out=print):
    flavours = spec["flavours"]
    fn_name, fn_args, terms = parse_map(mapfile)
    unknown = [p for p in flavours if p not in fn_args]
    if unknown:
        raise ValueError(f"spec flavours {unknown} not FN arguments "
                         f"{fn_args}")
    rep = {"errors": [], "warnings": [], "info": []}
    lines = []
    for aN, term in terms:
        lc = LineCheck(aN, term, flavours, sheet, rep,
                       subleading=(str(spec.get('colour', ''))
                                   .lower() == 'subleading'))
        if not lc.antennae:
            rep["warnings"].append(
                f"a{aN}: no antenna factor recognized in the "
                f"datasheet — line skipped "
                f"(factors: {[f[0] for f in term['factors']]})")
            continue
        lc.check_species()
        lc.check_split()
        lines.append(lc)

    # soft-eikonal dipoles against the declared colour chain(s)
    chains = spec.get("chains")
    if chains:
        bad = [c for c in chains
               if any(leg not in flavours for leg in c)]
        if bad:
            raise ValueError(f"spec chains contain non-parton legs: {bad}")
        adjacency = {frozenset(p) for c in chains
                     for p in zip(c, c[1:])}
        unchecked = []
        for lc in lines:
            lc.check_soft_dipoles(adjacency, unchecked)
        if unchecked:
            rep["info"].append(
                "soft-dipole check skipped on mapped momenta: "
                + ", ".join(f"a{n}:{t}" for n, t in unchecked))
    else:
        rep["info"].append(
            "no 'chains' in spec -> soft-eikonal dipoles UNCHECKED. Add "
            "\"chains\": [[...], ...] (colour orderings as lists of FN "
            "argument names; list every ordering the ME sums over) to "
            "enable the check that catches a soft eikonal sitting on a "
            "dipole the matrix element does not have.")

    # sign rule (rung 0): single-antenna lines +, iterated lines -
    for lc in lines:
        n_ant = len(lc.antennae)
        if n_ant == 1 and lc.term["sign"] < 0:
            lc.warn("single-antenna line with NEGATIVE sign "
                    "(S,a/S,b1 are + by the rung-0 sign rule)")
        if n_ant >= 2 and lc.term["sign"] > 0:
            lc.warn("iterated line with POSITIVE sign "
                    "(S,b2/S,c/S,d are - by the rung-0 sign rule)")

    # coverage of genuine modes
    line_poles = [(lc, lc.poles()) for lc in lines]
    for m in modes:
        if not m["genuine"]:
            continue
        fam, name = m["family"], m["name"]
        key = parse_mode_name(fam, name)
        if fam in ("ss", "sco", "ds", "tc"):
            need = ("ss", key) if fam == "ss" else \
                   ("sco", key) if fam == "sco" else \
                   ("ds", key) if fam == "ds" else ("tc", key)
            cov = [lc for lc, pl in line_poles
                   if need in pl and lc.term["sign"] > 0]
            if fam in ("ss", "sco"):
                cov = [lc for lc in cov
                       if any(len(a) == 3 for _, a, _ in lc.primary())]
            else:
                cov = [lc for lc in cov
                       if any(len(a) == 4 for _, a, _ in lc.primary())]
            if not cov:
                rep["errors"].append(
                    f"genuine {fam} mode '{name}' has NO covering "
                    + ("S,a line (3-parton antenna singular there)"
                       if fam in ("ss", "sco") else
                       "S,b1 line (4-parton antenna singular there)"))
            elif verbose:
                rep["info"].append(
                    f"{fam} '{name}' covered by "
                    f"{['a%d' % lc.aN for lc in cov]}")
        elif fam == "sc":
            soft, coll = key
            cov = [lc for lc, pl in line_poles
                   if ("ss", soft) in pl and ("sco", coll) in pl]
            rep["info"].append(
                f"sc '{name}': "
                + (f"single-line carriers {['a%d' % c.aN for c in cov]}"
                   if cov else
                   "no single line carries both invariants — coverage "
                   "lives in reduced-ME poles, unverified statically"))
        elif fam == "dc":
            pa, pb = tuple(key)
            cov = [lc for lc, pl in line_poles
                   if ("sco", pa) in pl and ("sco", pb) in pl]
            rep["info"].append(
                f"dc '{name}': "
                + (f"single-line carriers {['a%d' % c.aN for c in cov]}"
                   if cov else
                   "no single line carries both pairs — coverage "
                   "lives in reduced-ME poles, unverified statically"))

    # A term may be PARTIAL by design: e.g. an identical-flavour channel
    # whose single-unresolved limits are carried by the non-identical
    # sibling term. Declaring `"families"` in the spec says so, and the
    # single-unresolved bookkeeping below then reports instead of fails.
    fams = set(spec.get("families")
               or ["ss", "sco", "ds", "tc", "sc", "dc"])
    singles_here = bool(fams & {"ss", "sco"})
    note = ("" if singles_here else
            "  [spec declares no ss/sco families: singles are another "
            "term's job, reported not failed]")

    def single_issue(msg):
        (rep["errors"] if singles_here else rep["info"]).append(
            msg + note)

    # spurious-single pairing for X40 lines
    valid_cluster = make_valid_cluster(flavours)
    iterated = [(lc, lc.poles()) for lc in lines
                if len(lc.antennae) >= 2 and lc.term["sign"] < 0]
    for lc in lines:
        prim4 = [(t, a, e) for (t, a, e) in lc.primary()
                 if len(a) == 4]
        if not prim4 or len(lc.antennae) != 1:
            continue
        token, args, entry = prim4[0]
        meas = entry.get("measured", {})
        singles = []
        for key in meas.get("sco", {}):
            p, q = (int(v) for v in key.split(","))
            pair = frozenset([norm(args[p - 1]), norm(args[q - 1])])
            if valid_cluster(pair):
                singles.append(("sco", pair))
        for bslot in meas.get("ss", {}):
            leaf = norm(args[int(bslot) - 1])
            if flavours.get(leaf) == "g":
                singles.append(("ss", frozenset([leaf])))
        for kind, inv in singles:
            partners = [ic for ic, pl in iterated if (kind, inv) in pl]
            nm = ("s(" + ",".join(sorted(inv)) + ")") \
                if kind == "sco" else f"soft {next(iter(inv))}"
            if not partners:
                single_issue(
                    f"a{lc.aN}: spurious single pole {nm} of {token} "
                    f"has NO negative iterated counterterm on the same "
                    f"invariant — it will survive every "
                    f"single-unresolved limit touching it")
            elif verbose:
                rep["info"].append(
                    f"a{lc.aN}: {token} {nm} paired with "
                    f"{['a%d' % ic.aN for ic in partners]}")

    # --- SPLIT HALVES: a registered half never appears alone ----------
    allow_half, allow_orphan = read_pragmas(mapfile)
    partners = half_partners(sheet)
    used = {}                                   # token -> [(aN, args)]
    for lc in lines:
        for token, args, _e in lc.antennae:
            used.setdefault(token, []).append((lc.aN, [norm(a)
                                                       for a in args]))
    for token, calls in sorted(used.items()):
        if token not in partners or token in allow_half:
            continue
        other, mine, theirs = partners[token]
        want = []                            # (aN, own args, partner args)
        for aN, args in calls:
            pa = partner_args(args, mine, theirs)
            if pa is not None:
                want.append((aN, args, pa))

        def callstr(tok, argv):
            return tok + "(" + ",".join(show(x) for x in argv) + ")"

        if other not in used:
            recipe = "; ".join(callstr(other, pa) for _n, _a, pa in want)
            rep["errors"].append(
                f"{token} is one HALF of a split antenna and its partner "
                f"{other} appears nowhere: the halves carry different "
                f"mappings and cover DIFFERENT limits, so using one alone "
                f"leaves the other's poles unsubtracted. Identity "
                f"{token}({mine})+{other}({theirs}). Add: "
                + (recipe or f"{other}(...)"))
        else:
            # Compare as a SET of momenta, not an ordered tuple: within a
            # half, slots that the antenna is symmetric in may be written
            # either way (both conventions occur in verified terms) and
            # only the MAPPING differs, not the value.  What must not
            # happen is a half with no partner on those momenta at all.
            have = {frozenset(map(repr, a)) for _n, a in used[other]}
            for aN, args, pa in want:
                if frozenset(map(repr, pa)) not in have:
                    rep["warnings"].append(
                        f"a{aN}: {callstr(token, args)} has no {other} "
                        f"line on the same four momenta — expected "
                        + callstr(other, pa))

    # --- PAIRING MATRIX: X40 spurious singles <-> iterated lines ------
    # Forward (X40 -> counterterm) was checked above.  This is the
    # REVERSE direction: an iterated line is singular at leading power
    # in every leaf invariant its FIRST antenna carries, so each of
    # those must be a spurious single of some X40 line in the term.
    # Sharing the invariant with an S,a line is NOT enough — S,a is
    # exact there, so an unmatched counterterm subtracts a limit that
    # nothing restores.  (Empirical signature: the single-unresolved
    # mode on that invariant lands at ME/(ME - c), c = O(1).)
    x40_single_map = {}
    valid_cluster2 = make_valid_cluster(flavours)
    for lc in lines:
        if len(lc.antennae) != 1:
            continue
        for token, args, entry in lc.primary():
            if len(args) != 4:
                continue
            meas = entry.get("measured", {})
            for key in meas.get("sco", {}):
                p, q = (int(v) for v in key.split(","))
                pair = frozenset([norm(args[p - 1]), norm(args[q - 1])])
                if valid_cluster2(pair):
                    x40_single_map.setdefault(("sco", pair),
                                              []).append(lc.aN)
            for bslot in meas.get("ss", {}):
                leaf = norm(args[int(bslot) - 1])
                if flavours.get(leaf) == "g":
                    x40_single_map.setdefault(
                        ("ss", frozenset([leaf])), []).append(lc.aN)
    for ic, pl in iterated:
        leafpoles = sorted(
            (p for p in pl
             if p[0] in ("sco", "ss")
             and all(isinstance(v, str) for v in p[1])),
            key=lambda p: (p[0], sorted(p[1])))
        if not leafpoles:
            continue
        unmatched = [p for p in leafpoles if p not in x40_single_map]
        cells = ", ".join(
            f"{polename(*p)}->"
            + (",".join("a%d" % n for n in x40_single_map[p])
               if p in x40_single_map else "NONE")
            for p in leafpoles)
        rep["info"].append(f"pairing a{ic.aN}: {cells}")
        if unmatched and ic.aN not in allow_orphan:
            single_issue(
                f"a{ic.aN}: iterated counterterm is singular at leading "
                f"power in {', '.join(polename(*p) for p in unmatched)} "
                f"but NO X40 line has that spurious single pole — it "
                f"will over-subtract that single-unresolved limit. "
                f"Either add the X40 half that carries it, or move this "
                f"overlap to the other writing.")

    return rep


def polename(kind, inv):
    if kind == "ss":
        return f"soft {next(iter(inv))}"
    return "s(" + ",".join(sorted(inv)) + ")"


def half_partners(sheet):
    """half token -> (partner token, letters_self, letters_partner).

    Derived from the datasheet's registered `halves` + `split_identity`
    (e.g. "E40a(ABCD)+E40b(ADCB)"); nothing antenna-specific is encoded
    here, so a newly measured split is picked up for free.
    """
    out = {}
    for _full, e in sheet.items():
        halves = e.get("halves") or []
        ident = e.get("split_identity") or ""
        if len(halves) != 2 or not ident:
            continue
        perms = dict(re.findall(r"([A-Za-z0-9]+)\(([A-Z]+)\)", ident))
        a, b = halves
        if a in perms and b in perms:
            out[a] = (b, perms[a], perms[b])
            out[b] = (a, perms[b], perms[a])
    return out


def partner_args(args, letters_self, letters_partner):
    """Re-order a half's call arguments into its partner's call."""
    if len(args) != len(letters_self) or \
            sorted(letters_self) != sorted(letters_partner):
        return None
    m = dict(zip(letters_self, args))
    return [m[l] for l in letters_partner]


def read_pragmas(mapfile):
    """`# ledger: allow-half <token>` / `allow-orphan <aN>` escapes."""
    allow_half, allow_orphan = set(), set()
    try:
        text = open(mapfile).read()
    except OSError:
        return allow_half, allow_orphan
    for m in re.finditer(r"#\s*ledger:\s*allow-half\s+(\S+)", text):
        allow_half.add(m.group(1))
    for m in re.finditer(r"#\s*ledger:\s*allow-orphan\s+a?(\d+)", text):
        allow_orphan.add(int(m.group(1)))
    return allow_half, allow_orphan


def report(rep, out=print):
    for e in rep["errors"]:
        out(f"ERROR   {e}")
    for w in rep["warnings"]:
        out(f"warning {w}")
    for i in rep["info"]:
        out(f"info    {i}")
    out(f"\n{len(rep['errors'])} error(s), "
        f"{len(rep['warnings'])} warning(s)")
    return 1 if rep["errors"] else 0


# ---------------- selftest ----------------

def selftest():
    """Synthetic datasheet + synthetic .map — encodes no real antenna's
    physics (fake tokens ZZ30/YY40 with made-up pole graphs)."""
    sheet = {
        "ZZ30": {"arity": 3, "species": ["Q", "g", "Q"],
                 "slot_pairs": [[1, 3]],
                 "reduced": [["Q", 1], ["Q", 3]],
                 "measured": {"sco": {"1,2": {"power": 1},
                                      "2,3": {"power": 1}},
                              "ss": {"2": {"alpha": 2.0}},
                              "tc": {}, "ds": {}},
                 "requires_split": False},
        "YY40": {"arity": 4, "species": ["Q", "g", "g", "Q"],
                 "slot_pairs": [[1, 4]],
                 "reduced": [["Q", 1], ["Q", 4]],
                 "measured": {"sco": {"1,2": {"power": 1},
                                      "3,4": {"power": 1}},
                              "ss": {"2": {"alpha": 2.0},
                                     "3": {"alpha": 2.0}},
                              "tc": {"1,2,3": {"alpha": 2.0}},
                              "ds": {"2,3": {"alpha": 4.0}}},
                 "requires_split": False},
        "WW40": {"arity": 4, "species": ["Q", "g", "g", "Q"],
                 "slot_pairs": [], "reduced": [["Q", 1], ["Q", 4]],
                 "measured": {"sco": {"1,4": {"power": 1}},
                              "ss": {}, "tc": {}, "ds": {}},
                 "requires_split": True, "unsupported_poles": ["1,4"],
                 "halves": ["WW40a", "WW40b"], "fortran": "FullWW40",
                 "split_identity": "WW40a(ABCD)+WW40b(ADCB)"},
    }
    spec = {"flavours": {"p": "qb1", "g1": "g", "g2": "g", "q": "q1"},
            "born": [["q1", "qb1"]]}
    modes = [
        {"family": "ss", "name": "g1 soft", "genuine": True},
        {"family": "sco", "name": "p||g1", "genuine": True},
        {"family": "ds", "name": "g1,g2 soft", "genuine": True},
        {"family": "tc", "name": "p||g1||g2", "genuine": True},
        {"family": "sco", "name": "p||q", "genuine": False},
    ]

    def write_map(body):
        f = tempfile.NamedTemporaryFile("w", suffix=".map",
                                        delete=False)
        f.write("FN:=TESTS(p,g1,g2,q,1,2):\nXX:=\n" + body + "\n:\n")
        f.close()
        return f.name

    # (1) a correct little term: S,a + S,b1 + two iterated
    good = "\n".join([
        "+ZZ30(p,g1,g2)*RED([p,g1],[g1,g2],q,1,2)"
        "*JET33([p,g1],[g1,g2],q)*a1",
        "+ZZ30(p,g2,g1)*RED([p,g2],[g2,g1],q,1,2)"
        "*JET33([p,g2],[g2,g1],q)*a2",
        "+YY40(p,g1,g2,q)*RED2([p,g1,g2],[q,g2,g1],1,2)"
        "*JET22([p,g1,g2],[q,g2,g1])*a3",
        "-ZZ30(p,g1,g2)*ZZ30([p,g1],[g1,g2],q)"
        "*RED2([[p,g1],[g1,g2]],[q,[g1,g2]],1,2)"
        "*JET22([[p,g1],[g1,g2]],[q,[g1,g2]])*a4",
        "-ZZ30(p,g2,g1)*ZZ30([p,g2],[g2,g1],q)"
        "*RED2([[p,g2],[g2,g1]],[q,[g2,g1]],1,2)"
        "*JET22([[p,g2],[g2,g1]],[q,[g2,g1]])*a5",
    ])
    rep = run_ledger(write_map(good), spec, sheet, modes)
    # YY40's soft legs g1,g2 and sco (3,4)={g2,q} poles: soft g1/g2
    # matched by a4/a5 ZZ30 poles; sco {g2,q} matched by a5 second
    # antenna? -> second antenna args are clusters, not leaves; a5's
    # first antenna has (p,g2),(g2,g1) -> {g2,q} NOT matched => 1 error
    assert any("spurious single pole s(g2,q)" in e
               for e in rep["errors"]), rep["errors"]
    assert not any("mismatch" in e for e in rep["errors"])
    assert not any("mode" in e for e in rep["errors"])

    # (1b) soft-eikonal dipoles vs the declared chain.  ZZ30's soft leg
    # is slot 2, whose collinear partners are slots 1,3 -> dipole on the
    # outer legs; on the chain p-g1-g2-q that is (p,g2) for ZZ30(p,g1,g2)
    # and (p,g1) for ZZ30(p,g2,g1): both colour-connected? (p,g1) is,
    # (p,g2) is not -- so the chain check is genuinely discriminating.
    cspec = dict(spec, chains=[["p", "g1", "g2", "q"]])
    rep = run_ledger(write_map(good), cspec, sheet, modes)
    assert any("eikonal on the dipole (p,g2)" in e
               for e in rep["errors"]), rep["errors"]
    assert not any("eikonal on the dipole (p,g1)" in e
                   for e in rep["errors"]), rep["errors"]
    # both orderings declared -> both dipoles legal, no eikonal error
    cspec2 = dict(spec, chains=[["p", "g1", "g2", "q"],
                                ["p", "g2", "g1", "q"]])
    rep = run_ledger(write_map(good), cspec2, sheet, modes)
    assert not any("eikonal" in e for e in rep["errors"]), rep["errors"]
    # without chains the check is announced as skipped, never silent
    rep = run_ledger(write_map(good), spec, sheet, modes)
    assert any("UNCHECKED" in i for i in rep["info"]), rep["info"]

    # (1c) the DERIVED dipole rule must reproduce every dipole the real
    # datasheet records, otherwise deriving it for the X40s is unsound.
    try:
        real = json.load(open(DEFAULT_DATASHEET))
    except OSError:                       # asset absent: skip, not fail
        real = {}
    nchecked = 0
    for tok, e in real.items():
        meas = e.get("measured") or {}
        sco = set()
        for k in (meas.get("sco") or {}):
            try:
                a, b = (int(t) for t in k.split(","))
            except ValueError:
                continue
            sco.add(frozenset((a, b)))
        for b, d in (meas.get("ss") or {}).items():
            for k in (d.get("dipoles") or {}):
                want = frozenset(int(t) for t in k.split(","))
                ib = int(b)
                got = frozenset(next(iter(pr - {ib})) for pr in sco
                                if ib in pr)
                assert got == want, (tok, b, sorted(got), sorted(want))
                nchecked += 1
    if nchecked:
        print(f"  derived-dipole rule agrees with {nchecked} recorded "
              f"datasheet dipole(s)")

    # (2) structural typo: stale cluster in the second antenna
    bad = good.replace("*ZZ30([p,g1],[g1,g2],q)",
                       "*ZZ30([p,g2],[g1,g2],q)", 1)
    rep = run_ledger(write_map(bad), spec, sheet, modes)
    assert any("not an available momentum" in e for e in rep["errors"])

    # (3) species violation: gluon slot fed a quark
    bad2 = good.replace("+ZZ30(p,g1,g2)", "+ZZ30(p,q,g2)", 1)
    rep = run_ledger(write_map(bad2), spec, sheet, modes)
    assert any("slot 2 wants g" in e for e in rep["errors"])

    # (4) requires_split composite used directly
    bad3 = good + "\n+WW40(p,g1,g2,q)*RED2([p,g1,g2],[q,g2,g1],1,2)" \
                  "*JET22([p,g1,g2],[q,g2,g1])*a6"
    rep = run_ledger(write_map(bad3.replace("*a6", "*a6")), spec,
                     sheet, modes)
    assert any("MUST be split" in e and "WW40a(ABCD)+WW40b(ADCB)" in e
               for e in rep["errors"])

    # (5) missing coverage: drop the S,b1 line -> ds/tc uncovered
    nob1 = "\n".join(good.splitlines()[:2] + good.splitlines()[3:])
    nob1 = re.sub(r"\*a\d+", lambda m, c=itertools.count(1):
                  f"*a{next(c)}", nob1)
    rep = run_ledger(write_map(nob1), spec, sheet, modes)
    assert any("ds mode" in e and "NO covering" in e
               for e in rep["errors"])
    assert any("tc mode" in e for e in rep["errors"])

    # (6) sign rule
    flip = good.replace("+YY40", "-YY40")
    rep = run_ledger(write_map(flip), spec, sheet, modes)
    assert any("NEGATIVE sign" in w for w in rep["warnings"])
    print("pole_ledger selftest OK")


def main():
    if "--selftest" in sys.argv:
        selftest()
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("mapfile")
    ap.add_argument("--spec", required=True)
    ap.add_argument("--datasheet", default=DEFAULT_DATASHEET)
    ap.add_argument("--modes")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    spec = json.load(open(args.spec))
    sheet = json.load(open(args.datasheet))
    modes = genuine_modes(spec, args.modes)
    rep = run_ledger(args.mapfile, spec, sheet, modes,
                     verbose=args.verbose)
    sys.exit(report(rep))


if __name__ == "__main__":
    main()
