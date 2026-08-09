"""Block/marker/label/axis layer of a .map master file — shared library.

The ONE definition of the '# block:', '*aN' and '# axis:' conventions.
compose_blocks.py (the composer), audit_blocks.py and emit_skeleton.py
all import it, so the writer and the auditor of the format can never
drift apart. Not a command: it has no CLI; it is exercised by the
importing commands' --selftest.

Master-file convention: inside the XX:= body, group term lines under
maple comment markers

    # block: <name>

Everything before the first marker belongs to block "_default".
Everything outside the XX:= body (FN line, header comments, the
terminating ':') is structural and always preserved.

Axis convention (see compose_blocks.py for the workflow):

    # axis: <name> = <opt>:<blk>,<blk> | <opt>:<blk> | none:
"""
import re

MARKER = re.compile(r"^\s*#\s*block:\s*(\S+)")
AXIS = re.compile(r"^\s*#\s*axis:\s*(\w+)\s*=\s*(.+?)\s*$")
LABEL = re.compile(r"\*a(\d+)\b")


def split_map(text):
    """-> (head_lines, body_lines, tail_lines). Body = inside XX:= ... :"""
    lines = text.splitlines()
    try:
        istart = next(i for i, ln in enumerate(lines)
                      if ln.strip().replace(" ", "").startswith("XX:="))
    except StopIteration:
        raise SystemExit("ERROR: no 'XX:=' line found — not a term .map?")
    # terminating line: last line that is exactly ':' (possibly spaces)
    iend = None
    for i in range(len(lines) - 1, istart, -1):
        if lines[i].strip() == ":":
            iend = i
            break
    if iend is None:
        raise SystemExit("ERROR: no terminating ':' line after XX:=")
    return lines[:istart + 1], lines[istart + 1:iend], lines[iend:]


def blocks_of(body):
    """-> ordered list of (name, [lines])."""
    out = []
    cur_name, cur = "_default", []
    for ln in body:
        m = MARKER.match(ln)
        if m:
            out.append((cur_name, cur))
            cur_name, cur = m.group(1), [ln]
        else:
            cur.append(ln)
    out.append((cur_name, cur))
    # drop an empty _default
    return [(n, ls) for n, ls in out if n != "_default" or
            any(s.strip() for s in ls)]


def renumber(lines):
    """Renumber every *aN occurrence sequentially, gap-free from 1."""
    counter = [0]

    def sub(_m):
        counter[0] += 1
        return f"*a{counter[0]}"

    return [LABEL.sub(sub, ln) for ln in lines], counter[0]


def compose(text, selected):
    """Compose a .map from the named blocks, aN renumbered gap-free."""
    head, body, tail = split_map(text)
    blks = blocks_of(body)
    names = [n for n, _ in blks]
    unknown = [s for s in selected if s not in names]
    if unknown:
        raise SystemExit(f"ERROR: unknown block(s) {unknown}; "
                         f"available: {names}")
    kept = []
    for n, ls in blks:
        if n == "_default" or n in selected:
            kept.extend(ls)
    new_body, nterms = renumber(kept)
    if nterms == 0:
        raise SystemExit("ERROR: composition selects zero terms")
    return "\n".join(head + new_body + tail) + "\n", nterms


# ---------------- variant axes ----------------

def axes_of(text):
    """-> ordered dict {axis: [(option, [blocks]), ...]} from '# axis:' lines."""
    out = {}
    for ln in text.splitlines():
        m = AXIS.match(ln)
        if not m:
            continue
        name, rhs = m.group(1), m.group(2)
        opts = []
        for chunk in rhs.split("|"):
            chunk = chunk.strip()
            if not chunk:
                continue
            if ":" not in chunk:
                raise SystemExit(
                    f"ERROR: axis '{name}' option '{chunk}' is not <opt>:<blocks>")
            opt, blks = chunk.split(":", 1)
            opts.append((opt.strip(),
                         [b.strip() for b in blks.split(",") if b.strip()]))
        if len(opts) < 2:
            raise SystemExit(f"ERROR: axis '{name}' needs >= 2 options")
        if name in out:
            raise SystemExit(f"ERROR: axis '{name}' declared twice")
        out[name] = opts
    return out


def check_axes(text):
    """Every block an axis names must exist. Returns the axis table."""
    _, body, _ = split_map(text)
    known = {n for n, _ in blocks_of(body)}
    ax = axes_of(text)
    for name, opts in ax.items():
        for opt, blks in opts:
            missing = [b for b in blks if b not in known]
            if missing:
                raise SystemExit(
                    f"ERROR: axis '{name}' option '{opt}' names unknown "
                    f"block(s): {','.join(missing)}")
    return ax


def expand(text, fixed, selection):
    """fixed: [blocks]; selection: {axis: option} -> ordered block list."""
    ax = check_axes(text)
    out = list(fixed)
    for name, opt in selection.items():
        if name not in ax:
            raise SystemExit(f"ERROR: unknown axis '{name}'")
        match = [b for o, b in ax[name] if o == opt]
        if not match:
            have = ",".join(o for o, _ in ax[name])
            raise SystemExit(
                f"ERROR: axis '{name}' has no option '{opt}' (have: {have})")
        for b in match[0]:
            if b not in out:
                out.append(b)
    return out
