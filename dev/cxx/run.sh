#!/bin/sh
# Build and run the expression-template harness.
#
#   dev/cxx/run.sh              build, run, and run again under ASan/UBSan
#   dev/cxx/run.sh --record F   write the numeric output to F (reference run)
#   dev/cxx/run.sh --against F  diff this build's output against F
#
# The output is the assertion: two revisions that compute the same thing must
# produce byte-identical output.
set -eu

REPO=$(cd "$(dirname "$0")/../.." && pwd)
SRC="$REPO/dev/cxx/test_dual_expr.cpp"
OUT=${TMPDIR:-/tmp}/cppde_etest
RINC=$(Rscript -e 'cat(R.home("include"))')
RLIB=$(Rscript -e 'cat(R.home("lib"))')

CXX=${CXX:-g++}
STD=-std=gnu++17
INC="-I $REPO/inst/include -I $RINC"
# cppde.hpp declares the BLAS/LAPACK entry points R provides; link against R
# so any that get instantiated resolve.
LIBS="-L $RLIB -lR"

# -O2 matches how generated models are built.
$CXX $STD -O2 -DNDEBUG -Wall -Wextra $INC -o "$OUT" "$SRC" $LIBS

# -O1, not -O0: CPPDE_ET_INLINE is always_inline, which gcc can refuse to
# honour at -O0.
$CXX $STD -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
     $INC -o "$OUT.asan" "$SRC" $LIBS

export LD_LIBRARY_PATH="$RLIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

case "${1:-}" in
  --record)
    "$OUT" > "$2"
    echo "recorded -> $2"
    ;;
  --against)
    "$OUT" > "$OUT.txt"
    if diff -u "$2" "$OUT.txt"; then
      echo "output identical to $2"
    else
      echo "OUTPUT DIFFERS from $2" >&2
      exit 1
    fi
    ;;
  *)
    "$OUT"
    ;;
esac

ASAN_OPTIONS=detect_stack_use_after_scope=1 "$OUT.asan" > /dev/null
echo "asan/ubsan clean"
