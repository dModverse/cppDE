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
R_HOME=${R_HOME:-"C:/Program Files/R/R-4.6.1"}
CXX=${CXX:-g++}
OUT="$TMPDIR"
"$CXX" -std=c++17 -O1 \
  -I inst/include -I "$R_HOME/include" \
  dev/jump-adjoint-check.cpp -o "$OUT/jumpcheck" \
  -L "$R_HOME/bin/x64" -lRblas -lRlapack -lR
"$OUT/jumpcheck"
