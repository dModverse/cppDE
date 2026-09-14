#!/bin/sh
# The root jump's adjoint against the nested-dual sandwich, with no solver, no
# store and no event engine around it. Seconds to build, milliseconds to run,
# so a term can be derived and checked without a package rebuild.
#
#   sh dev/run-jump-check.sh
set -e
# The toolchain writes its intermediates where TMPDIR points, and on Windows
# that can inherit a directory the user cannot write.
: "${TMPDIR:=${HOME:-.}/.cache/cppde-jumpcheck}"
mkdir -p "$TMPDIR"
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"
REPO=$(cd "$(dirname "$0")/.." && pwd)
RINC=$(Rscript -e 'cat(R.home("include"))')
# Windows keeps no import libraries under R.home("lib"); bin/<arch> holds the
# DLLs and mingw links straight against those.
RLIB=$(Rscript -e 'cat(if (.Platform$OS.type == "windows") R.home(file.path("bin", .Platform$r_arch)) else R.home("lib"))')

CXX=${CXX:-g++}
OUT="$TMPDIR"
# libRlapack is absent where R takes LAPACK from the BLAS it links, as a
# FlexiBLAS build does; there the symbols come out of libRblas.
LIBS="-L $RLIB -lR -lRblas"
for f in "$RLIB"/libRlapack.* "$RLIB"/Rlapack.*; do
  if [ -e "$f" ]; then LIBS="$LIBS -lRlapack"; break; fi
done

"$CXX" -std=c++17 -O1 \
  -I "$REPO/inst/include" -I "$RINC" \
  "$REPO/dev/jump-adjoint-check.cpp" -o "$OUT/jumpcheck" \
  $LIBS

case $(uname -s) in
  MINGW*|MSYS*|CYGWIN*) PATH=$(cygpath -u "$RLIB"):$PATH; export PATH ;;
  *) export LD_LIBRARY_PATH="$RLIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" ;;
esac

"$OUT/jumpcheck"
