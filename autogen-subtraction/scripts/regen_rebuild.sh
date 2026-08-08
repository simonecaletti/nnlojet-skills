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
#                    -t <test/process/PROC> -m <make target> [-j N] \
#                    [--skip-ledger]
#   regen_rebuild.sh --selftest
#
# All paths relative to the repo root (found via `hg root`, falling back
# to the current directory). Aborts on generator ERRORs (aN gaps,
# invalid ME argument, left-over list). No maple available -> aborts:
# NEVER hand-edit auto*.f to mimic the generator.
#
# STEP 0 — the pole ledger is not optional. For every .map in
# maple/process/<DIR>/ carrying a sibling spec <TERM>.spec.json, the
# static ledger (write-subtraction/scripts/pole_ledger.py) runs BEFORE
# maple does, and a ledger ERROR aborts the run. It costs seconds and
# catches the classes that otherwise cost one build cycle each: an
# unpaired X40 spurious pole, a split half used without its partner, a
# stale cluster, an orphaned counterterm. A .map with no spec is
# reported UNCHECKED — the spec is four lines, write it. Escapes:
#   --skip-ledger      one-off bypass (prints a warning)
#   LEDGER_REQUIRED=1  make "no spec" a hard error too

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

  # Ledger gate: a .map that fails the static ledger must abort the run
  # BEFORE maple is invoked. Synthetic term — a real antenna token but a
  # stale cluster, so the failure is mechanical and encodes no physics.
  sdir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
  self="$sdir/$(basename "${BASH_SOURCE[0]}")"
  here=$(cd "$sdir/../../.." && pwd)          # the .claude directory
  g=$(mktemp -d)
  mkdir -p "$g/maple/process/FOO" "$g/.claude"
  ln -s "$here/skills" "$g/.claude/skills"
  cat > "$g/maple/process/FOO/T.map" <<'EOF'
FN:=T(l,m,k,1,2):
XX:=
+A30FF(l,m,k)*B1g0Z([l,zzz],[m,k],1,2)*JET33([l,zzz],[m,k])*a1
:
EOF
  cat > "$g/maple/process/FOO/T.spec.json" <<'EOF'
{"flavours": {"l":"qb1","m":"g","k":"q1"}, "born": [["q1","qb1","g"]],
 "partons": {"l":"qb1","m":"g","k":"q1"}}
EOF
  out=$( cd "$g" && bash "$self" -n 1 -l RR -s src/process/FOO \
         -t test/process/FOO -m dummy 2>&1 ); rc=$?
  rm -rf "$g"
  if [ $rc -eq 0 ] || ! printf '%s' "$out" | grep -q "LEDGER FAIL"; then
    echo "selftest FAIL: ledger gate did not abort on a bad .map"
    printf '%s\n' "$out" | head -5
    exit 1
  fi
  echo "regen_rebuild selftest OK"
  exit 0
}

[ "${1:-}" = "--selftest" ] && selftest

NPROC="" LAYER="" SRCDIR="" TESTDIR="" TARGET="" JOBS=8 SKIP_LEDGER=0
# getopts has no long options: pull the word flags out of $@ first.
ARGS=()
for a in "$@"; do
  case $a in
    --skip-ledger) SKIP_LEDGER=1 ;;
    *) ARGS+=("$a") ;;
  esac
done
set -- ${ARGS+"${ARGS[@]}"}
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

# --- 0. static pole ledger (seconds; not optional) --------------------
LEDGER="$ROOT/.claude/skills/write-subtraction/scripts/pole_ledger.py"
MAPDIR="maple/process/$(basename "$SRCDIR")"
if [ "$SKIP_LEDGER" = 1 ]; then
  echo "WARNING: --skip-ledger — the static pole ledger did NOT run."
  echo "         Structural errors will now cost a full build cycle each."
elif [ ! -f "$LEDGER" ] || ! command -v python3 >/dev/null 2>&1; then
  echo "WARNING: pole_ledger.py or python3 unavailable — ledger skipped."
elif [ ! -d "$MAPDIR" ]; then
  echo "WARNING: no $MAPDIR — ledger skipped."
else
  nunchk=0 nfail=0
  for mp in "$MAPDIR"/*.map; do
    [ -f "$mp" ] || continue
    base=$(basename "$mp" .map)
    case $base in auto*|*_master) continue ;; esac
    spec="$MAPDIR/$base.spec.json"
    if [ ! -f "$spec" ]; then
      echo "  UNCHECKED  $base.map   (no $base.spec.json)"
      nunchk=$((nunchk + 1)); continue
    fi
    out=$(python3 "$LEDGER" "$mp" --spec "$spec" 2>&1); rc=$?
    if [ $rc -ne 0 ]; then
      echo "  LEDGER FAIL  $base.map"
      printf '%s\n' "$out" | grep '^ERROR' | sed 's/^/    /'
      nfail=$((nfail + 1))
    else
      echo "  ledger ok    $base.map   ($(printf '%s\n' "$out" | tail -1))"
    fi
  done
  if [ "$nfail" -gt 0 ]; then
    echo
    echo "ABORT: $nfail .map file(s) fail the static pole ledger."
    echo "These are structure errors — fix them in the .map"
    echo "(write-subtraction skill) before spending a build cycle."
    echo "Deliberately testing a partial term? re-run --skip-ledger."
    exit 1
  fi
  if [ "$nunchk" -gt 0 ] && [ "${LEDGER_REQUIRED:-0}" = 1 ]; then
    echo "ABORT: LEDGER_REQUIRED=1 and $nunchk .map file(s) have no spec."
    exit 1
  fi
fi

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
