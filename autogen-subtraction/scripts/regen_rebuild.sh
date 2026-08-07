#!/usr/bin/env bash
# One-command regenerate + rebuild for an NNLOJET subtraction term under
# spike test: consistency check -> Fortran generation (with semantic
# restore of files that only changed in boilerplate) -> test build.
# Turns "does this .map edit work?" into a single short command — the
# iteration primitive behind block bisection (run-spike-test) and the
# block composer (write-subtraction).
#
# Usage:
#   regen_rebuild.sh -n <iprocess> -l <RR|RV> -s <src/process/DIR> \
#                    -t <test/process/PROC> -m <make target> [-j N]
#   regen_rebuild.sh --selftest
#
# All paths relative to the repo root (found via `hg root`, falling back
# to the current directory). Aborts on generator ERRORs (aN gaps,
# invalid ME argument, left-over list). No maple available -> aborts:
# NEVER hand-edit auto*.f to mimic the generator.

set -u

# --- semantic filter: the physics content of a generated auto*.f ------
semantic() {  # $1 = file
  grep -E 'set_map|wt\(|getqcdnorm|set_flav_perm|\* FF:' "$1" \
    | sed 's/^[[:space:]]*//'
}

selftest() {
  # Structural self-test: no maple, no physics answer. Verifies that the
  # semantic filter distinguishes a boilerplate-only change from a real
  # one — the property the restore step relies on.
  t=$(mktemp -d)
  cat > "$t/a.f" <<'EOF'
      double precision partons(10)
      call set_map(7,6, (/3,4,5/), (/2,1,3,5,6,7/), ipass)
      wt(1) = x*FullA30FF(3,4,5,7)
      call getqcdnorm(5,facnorm)
EOF
  # boilerplate-only change (declaration style), semantics identical
  sed 's/double precision partons(10)/      dimension partons(10)/' \
      "$t/a.f" > "$t/b.f"
  # semantic change (different mapping)
  sed 's/(\/3,4,5\/)/(\/3,4,6\/)/' "$t/a.f" > "$t/c.f"
  if ! diff <(semantic "$t/a.f") <(semantic "$t/b.f") >/dev/null; then
    echo "selftest FAIL: boilerplate change flagged as semantic"; exit 1
  fi
  if diff <(semantic "$t/a.f") <(semantic "$t/c.f") >/dev/null; then
    echo "selftest FAIL: semantic change not detected"; exit 1
  fi
  rm -rf "$t"
  echo "regen_rebuild selftest OK"
  exit 0
}

[ "${1:-}" = "--selftest" ] && selftest

NPROC="" LAYER="" SRCDIR="" TESTDIR="" TARGET="" JOBS=8
while getopts "n:l:s:t:m:j:" opt; do
  case $opt in
    n) NPROC=$OPTARG ;;
    l) LAYER=$OPTARG ;;
    s) SRCDIR=$OPTARG ;;
    t) TESTDIR=$OPTARG ;;
    m) TARGET=$OPTARG ;;
    j) JOBS=$OPTARG ;;
    *) exit 2 ;;
  esac
done
if [ -z "$NPROC" ] || [ -z "$LAYER" ] || [ -z "$SRCDIR" ] \
   || [ -z "$TESTDIR" ] || [ -z "$TARGET" ]; then
  sed -n '2,20p' "$0"; exit 2
fi
case $LAYER in RR|RV) ;; *) echo "ERROR: -l must be RR or RV"; exit 2;; esac

ROOT=$(hg root 2>/dev/null || pwd)
cd "$ROOT" || exit 1
command -v maple >/dev/null 2>&1 || {
  echo "ERROR: no maple on PATH. Do NOT hand-edit auto*.f; only"
  echo "already-checked-in auto*.f can be spike-tested."; exit 1; }

# --- 1. consistency check (catches aN gaps) ---------------------------
log=$(mktemp)
( cd maple && maple "make${LAYER}check" -Diprocess="$NPROC" ) \
  > "$log" 2>&1
if grep -qi "error" "$log"; then
  echo "make${LAYER}check FAILED:"; grep -i "error" "$log"; exit 1
fi

# --- 2. backup, generate, semantic restore ----------------------------
BK=$(mktemp -d)
cp -a "$SRCDIR"/auto*.f "$BK"/ 2>/dev/null
( cd maple && maple "makefort${LAYER}" -Diprocess="$NPROC" ) \
  > "$log" 2>&1
if grep -qiE "invalid ME argument|left-over list|error" "$log"; then
  echo "makefort${LAYER} FAILED (malformed .map? -> write-subtraction):"
  grep -iE "invalid ME argument|left-over list|error" "$log"
  exit 1
fi
nrestored=0
for f in "$SRCDIR"/auto*.f; do
  b="$BK/$(basename "$f")"
  [ -f "$b" ] || continue                      # brand-new file: keep
  if diff <(semantic "$f") <(semantic "$b") >/dev/null 2>&1; then
    cp "$b" "$f"                               # boilerplate-only churn
    nrestored=$((nrestored + 1))
  fi
done
echo "generated $(ls "$SRCDIR"/auto*.f | wc -l) file(s);" \
     "restored $nrestored boilerplate-only rewrite(s)"

# --- 3. rebuild the spike test ----------------------------------------
cd "$ROOT/$TESTDIR" || exit 1
if ! make "$TARGET" -j"$JOBS"; then
  echo "first make failed (module race?) — retrying with -j1"
  make "$TARGET" -j1 || exit 1
fi
echo "OK: $TESTDIR/$TARGET rebuilt — run it via run-spike-test"
