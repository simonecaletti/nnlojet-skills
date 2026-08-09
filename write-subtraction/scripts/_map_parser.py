"""Term-level parser of a .map subtraction file — shared library.

Parses the FN:= line and every XX term into sign, coefficient and
factors, with momenta as leaves (strings) or clusters (nested tuples),
normalized direction-insensitively (bracket contents are equal up to
reversal, matching getpmapIK's convention). Imported by pole_ledger.py;
not a command — its behaviour is exercised by pole_ledger.py --selftest.
"""
import re


# ---------------- momentum terms ----------------

def norm(x):
    """Normalize a momentum: leaf -> str; cluster -> tuple, children
    normalized, direction-insensitive (contents equal up to reversal)."""
    if isinstance(x, str):
        return x
    t = tuple(norm(c) for c in x)
    return min(t, t[::-1], key=repr)


def leaves(x):
    if isinstance(x, str):
        return [x]
    out = []
    for c in x:
        out.extend(leaves(c))
    return out


def show(x):
    if isinstance(x, str):
        return x
    return "[" + ",".join(show(c) for c in x) + "]"


# ---------------- .map parsing ----------------

def split_args(s):
    out, depth, cur = [], 0, ""
    for ch in s:
        if ch in "[(":
            depth += 1
        elif ch in "])":
            depth -= 1
        if ch == "," and depth == 0:
            out.append(cur.strip())
            cur = ""
        else:
            cur += ch
    if cur.strip():
        out.append(cur.strip())
    return out


def parse_arg(s):
    s = s.strip()
    if s.startswith("["):
        assert s.endswith("]"), s
        return tuple(parse_arg(a) for a in split_args(s[1:-1]))
    return s


FACTOR = re.compile(r"([A-Za-z]\w*)\(")


def parse_term(text):
    """One XX term -> {sign, coeff, factors:[(name,[args])]}."""
    text = text.strip()
    sign = 1
    while text and text[0] in "+-":
        if text[0] == "-":
            sign = -sign
        text = text[1:].lstrip()
    m = re.match(r"(\d+(?:/\d+)?)\s*\*", text)
    coeff = 1.0
    if m:
        p = m.group(1).split("/")
        coeff = float(p[0]) / (float(p[1]) if len(p) > 1 else 1.0)
        text = text[m.end():]
    factors = []
    pos = 0
    while pos < len(text):
        m = FACTOR.search(text, pos)
        if not m:
            break
        depth = 1
        k = m.end()
        while k < len(text) and depth:
            if text[k] in "[(":
                depth += 1
            elif text[k] in "])":
                depth -= 1
            k += 1
        args = [parse_arg(a) for a in split_args(text[m.end():k - 1])]
        factors.append((m.group(1), args))
        pos = k
    return {"sign": sign, "coeff": coeff, "factors": factors}


def parse_map(path):
    """-> (fn_name, fn_args, [terms]) with terms carrying their aN."""
    raw = open(path, errors="replace").read()
    lines = [ln.split("#", 1)[0] for ln in raw.splitlines()]
    text = "\n".join(lines)
    m = re.search(r"FN\s*:=\s*(\w+)\(([^)]*)\)", text)
    if not m:
        raise ValueError("no FN:= line")
    fn_name = m.group(1)
    fn_args = [a.strip() for a in m.group(2).split(",")]
    m2 = re.search(r"XX\s*:=(.*?)\n\s*:", text, re.DOTALL)
    if not m2:
        raise ValueError("no XX:= ... : body")
    body = m2.group(1)
    terms = []
    # terms end at *aN labels
    pat = re.compile(r"\*\s*a(\d+)\b")
    pos = 0
    for m3 in pat.finditer(body):
        chunk = body[pos:m3.start()]
        terms.append((int(m3.group(1)),
                      parse_term(" ".join(chunk.split()))))
        pos = m3.end()
    return fn_name, fn_args, terms
