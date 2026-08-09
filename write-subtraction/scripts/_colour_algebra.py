"""Flavour species and colour-chain geometry — shared library.

Imported by predict_blocks.py / emit_skeleton.py / audit_blocks.py (chain
geometry, cluster flavour) and pole_ledger.py (flavour parsing, valid
collinear clusters). Not a command: it has no CLI. Its behaviour is
exercised by the importing commands' --selftest.

Flavour strings are 'g', 'q<tag>' or 'qb<tag>'. cluster_parent applies
the same net-flavour rule as genuine_modes.py (write-spike-test) — kept
separate because the two skills must stay independently copyable; if the
two ever disagree, genuine_modes.py wins.
"""


# ---------------- flavour species ----------------

def species(f):
    """flavour string -> 'g' or 'q' (quark and antiquark both 'q' for
    radiator-pair purposes)."""
    if f == "g":
        return "g"
    if f.startswith("q"):
        return "q"
    raise ValueError(f"bad flavour {f!r} (want g / q<tag> / qb<tag>)")


def fparse(f):
    """flavour string -> (kind, tag): ('g', None), ('q', tag), ('qb', tag)."""
    if f == "g":
        return ("g", None)
    if f.startswith("qb"):
        return ("qb", f[2:])
    if f.startswith("q"):
        return ("q", f[1:])
    raise ValueError(f"bad flavour {f!r}")


def kind_class(f):
    """'g' for a gluon, 'Q' for any quark kind — the antenna slot classes."""
    return "g" if f == "g" else "Q"


def cluster_parent(flavs):
    """Net flavour of a collinear cluster: 'g', 'q<tag>'/'qb<tag>', or
    None if the net is not a single parton."""
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


def make_valid_cluster(flavours):
    """-> predicate(frozenset of leaf labels): is this pair a flavour-valid
    collinear cluster (q||g, g||g, or same-flavour q/qb)?"""
    def valid(pair):
        fl = [flavours.get(v) for v in pair]
        if any(f is None for f in fl):
            return False
        kinds = [fparse(f) for f in fl]
        ks = sorted(k for k, t in kinds)
        if "g" in [k for k, t in kinds]:
            return True                      # q||g or g||g
        return ks == ["q", "qb"] and kinds[0][1] == kinds[1][1]
    return valid


# ---------------- colour-chain geometry ----------------

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


def config_of(chain_initial, legs):
    """FF / IF / FI / II from which radiators are initial-state."""
    flags = [l in chain_initial for l in legs]
    if not any(flags):
        return "FF"
    if all(flags):
        return "II"
    return "IF" if flags[0] else "FI"
