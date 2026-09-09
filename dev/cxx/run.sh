#!/bin/sh
# Build and run a direct C++ harness.
#
#   dev/cxx/run.sh              build, run, and run again under ASan/UBSan
#   dev/cxx/run.sh --record F   write the numeric output to F (reference run)
#   dev/cxx/run.sh --against F  diff this build's output against F
#
# A leading --codual, --reverse-step, --reverse-step-rb4,
# --reverse-step-multistep, --reverse-trajectory, --reverse-trajectory-methods
# --reverse-events, --err-weights or --sparse-transpose selects a reverse-AD
# harness.
#
# The output is the assertion: two revisions that compute the same thing must
# produce byte-identical output.
set -eu

REPO=$(cd "$(dirname "$0")/../.." && pwd)
SRC="$REPO/dev/cxx/test_dual_expr.cpp"
OUT=${TMPDIR:-/tmp}/cppde_etest

case "${1:-}" in
  --codual)
    SRC="$REPO/dev/cxx/test_codual.cpp"
    OUT=${TMPDIR:-/tmp}/cppde_codual
    shift
    ;;
  --reverse-step)
    SRC="$REPO/dev/cxx/test_reverse_step.cpp"
    OUT=${TMPDIR:-/tmp}/cppde_reverse_step
    shift
    ;;
  --reverse-step-rb4)
    SRC="$REPO/dev/cxx/test_reverse_step_rb4.cpp"
    OUT=${TMPDIR:-/tmp}/cppde_reverse_step_rb4
    shift
    ;;
  --reverse-step-multistep)
    SRC="$REPO/dev/cxx/test_reverse_step_multistep.cpp"
    OUT=${TMPDIR:-/tmp}/cppde_reverse_step_multistep
    shift
    ;;
  --reverse-trajectory-methods)
    SRC="$REPO/dev/cxx/test_reverse_trajectory_methods.cpp"
    OUT=${TMPDIR:-/tmp}/cppde_reverse_trajectory_methods
    shift
    ;;
  --reverse-trajectory)
    SRC="$REPO/dev/cxx/test_reverse_trajectory.cpp"
    OUT=${TMPDIR:-/tmp}/cppde_reverse_trajectory
    shift
    ;;
  --reverse-events)
    SRC="$REPO/dev/cxx/test_reverse_events.cpp"
    OUT=${TMPDIR:-/tmp}/cppde_reverse_events
    shift
    ;;
  --err-weights)
    SRC="$REPO/dev/cxx/test_err_weights.cpp"
    OUT=${TMPDIR:-/tmp}/cppde_err_weights
    shift
    ;;
  --sparse-transpose)
    SRC="$REPO/dev/cxx/test_sparse_transpose.cpp"
    OUT=${TMPDIR:-/tmp}/cppde_sparse_transpose
    # KLU comes from the install-time probe, so the harness links what a
    # generated sparse model links. An uninstalled package leaves it empty and
    # the test reports itself skipped rather than failing to build.
    KLU=$(Rscript -e 'cfg <- try(get("cvodeConfig", envir = asNamespace("cppDE")), silent = TRUE); if (!inherits(cfg, "try-error") && isTRUE(cfg$klu_available)) cat("-DKLU", cfg$klu_cflags, cfg$klu_libs)' 2>/dev/null || true)
    shift
    ;;
esac
RINC=$(Rscript -e 'cat(R.home("include"))')
# Windows keeps no import libraries under R.home("lib"); bin/<arch> holds the
# DLLs and mingw links straight against those.
RLIB=$(Rscript -e 'cat(if (.Platform$OS.type == "windows") R.home(file.path("bin", .Platform$r_arch)) else R.home("lib"))')

CXX=${CXX:-g++}
STD=-std=gnu++17
INC="-I $REPO/inst/include -I $RINC"
# cppde.hpp declares the BLAS/LAPACK entry points R provides; link against R
# so any that get instantiated resolve.
# libRlapack is absent where R takes LAPACK from the BLAS it links, as a
# FlexiBLAS build does; there the symbols come out of libRblas.
LIBS="-L $RLIB -lR -lRblas"
for f in "$RLIB"/libRlapack.* "$RLIB"/Rlapack.*; do
  if [ -e "$f" ]; then LIBS="$LIBS -lRlapack"; break; fi
done

# -O2 matches how generated models are built.
$CXX $STD -O2 -DNDEBUG -Wall -Wextra $INC ${KLU:-} -o "$OUT" "$SRC" $LIBS

# Windows resolves the DLLs off PATH, which needs the POSIX spelling of RLIB.
case $(uname -s) in
  MINGW*|MSYS*|CYGWIN*) PATH=$(cygpath -u "$RLIB"):$PATH; export PATH ;;
  *) export LD_LIBRARY_PATH="$RLIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" ;;
esac

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

# -O1, not -O0: CPPDE_ET_INLINE is always_inline, which gcc can refuse to
# honour at -O0. Skipped where the sanitizer runtimes are not installed.
if $CXX $STD -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
        $INC ${KLU:-} -o "$OUT.asan" "$SRC" $LIBS 2>/dev/null; then
  ASAN_OPTIONS=detect_stack_use_after_scope=1 "$OUT.asan" > /dev/null
  echo "asan/ubsan clean"
else
  echo "asan/ubsan skipped: no sanitizer runtime"
fi
