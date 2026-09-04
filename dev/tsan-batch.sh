#!/usr/bin/env bash
# ThreadSanitizer over solveODEBatch().
#
# Two things make this awkward, both handled here:
#
#   * libgomp is not instrumented, so TSan cannot see the barriers around a
#     parallel region and reports every value written inside it and read after
#     it.  run_batch() annotates both edges (see CPPDE_TSAN_HB/HA in
#     cppde_r_batch.hpp), which removes the bulk of that noise.  What remains
#     are reads of the outlined function's argument block on the main thread's
#     stack at region entry, before any statement we could annotate.  Those
#     are expected; a report is only interesting if BOTH accesses are worker
#     threads.
#
#   * Codegen runs Python through reticulate, and LD_PRELOAD would instrument
#     `uv` too, drowning the output.  So generate and compile first, without
#     the preload, and only run the solves under it.
#
# Needs libtsan.  Without root:
#   dnf download libtsan && rpm2cpio libtsan-*.rpm | cpio -idm
#   printf 'INPUT ( %s/usr/lib64/libtsan.so.2.0.0 )\n' "$PWD" > <dir>/libtsan.so
# and point TSAN_LIBDIR / TSAN_SCRIPTDIR at them.
set -euo pipefail

PKG="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${TMPDIR:-/tmp}/cppde-tsan"
TSAN_LIBDIR="${TSAN_LIBDIR:-/usr/lib64}"
TSAN_SCRIPTDIR="${TSAN_SCRIPTDIR:-}"
TSAN_SO="$TSAN_LIBDIR/libtsan.so.2.0.0"

[ -f "$TSAN_SO" ] || { echo "libtsan not found at $TSAN_SO (set TSAN_LIBDIR)"; exit 1; }

rm -rf "$WORK"; mkdir -p "$WORK"

# --- 1. generate + compile, no preload -------------------------------------
Rscript -e "
  suppressMessages(devtools::load_all('$PKG', quiet = TRUE))
  for (cfg in list(list(nm='t_d1', d2=FALSE, ns=NULL),
                   list(nm='t_d2', d2=TRUE,  ns=4L),
                   list(nm='t_hp', d2=FALSE, ns=Inf))) {
    m <- cppODE(c(A='-k*A', B='k*A', C='0.3*B-0.1*C'), modelname=cfg\$nm,
                deriv=TRUE, deriv2=cfg\$d2, nStack=cfg\$ns,
                outdir='$WORK', compile=FALSE, verbose=FALSE)
    saveRDS(m, file.path('$WORK', paste0(cfg\$nm, '.rds')))
  }" >/dev/null

LIBS="-fopenmp -fsanitize=thread"
[ -n "$TSAN_SCRIPTDIR" ] && LIBS="$LIBS -L$TSAN_SCRIPTDIR"
LIBS="$LIBS -Wl,-rpath,$TSAN_LIBDIR"

for f in t_d1 t_d2 t_hp; do
  ( cd "$WORK" && \
    PKG_CXXFLAGS="-O1 -g -w -fPIC -fopenmp -fsanitize=thread" \
    PKG_CPPFLAGS="-I$PKG/inst/include" \
    PKG_LIBS="$LIBS" \
    R CMD SHLIB "$f.cpp" >/dev/null 2>&1 ) || { echo "build failed: $f"; exit 1; }
done

# --- 2. run the solves under TSan ------------------------------------------
cat > "$WORK/run.R" <<'RRR'
suppressMessages(devtools::load_all(Sys.getenv("PKG"), quiet = TRUE))
w <- Sys.getenv("WORK"); ext <- .Platform$dynlib.ext
tt <- seq(0, 4, 0.2)
cs <- lapply(seq(0.3, 2.0, length.out = 16),
             function(k) list(parms = c(A = 1, B = 0, C = 0, k = k)))
for (nm in c("t_d1", "t_d2", "t_hp")) {
  m <- readRDS(file.path(w, paste0(nm, ".rds")))
  dyn.load(file.path(w, paste0(nm, ext))); cppDE::clearNativeSymbols()
  invisible(solveODE(m, times = tt, parms = c(A = 1, B = 0, C = 0, k = 0.7)))
  a <- solveODEBatch(m, cs, times = tt, cores = 8)
  b <- solveODEBatch(m, cs, times = tt, cores = 1)
  stopifnot(all(mapply(function(x, y) identical(x$variable, y$variable), a, b)))
  cat("ok:", nm, "\n")
}
RRR

PKG="$PKG" WORK="$WORK" LD_PRELOAD="$TSAN_SO" \
TSAN_OPTIONS="halt_on_error=0:report_bugs=1:history_size=7" \
  Rscript "$WORK/run.R" > "$WORK/out.txt" 2>&1 || true

# --- 3. verdict: only worker-vs-worker reports matter ----------------------
python3 - "$WORK/out.txt" <<'PYY'
import re, sys
txt = open(sys.argv[1], errors="replace").read()
reps = [r for r in re.split(r"(?=WARNING: ThreadSanitizer: data race)", txt) if r.startswith("WARNING")]
ww = [r for r in reps
      if all(w.startswith("thread T")
             for w in re.findall(r"of size \d+ at 0x[0-9a-f]+ by (main thread|thread T\d+)", r)[:2])]
print(f"reports: {len(reps)}   worker-vs-worker: {len(ww)}")
for line in re.findall(r"^ok: .*$", txt, re.M): print(line)
if ww:
    print("\n--- worker-vs-worker (these are real) ---")
    print("\n".join(ww[0].splitlines()[:25]))
    sys.exit(1)
print("no worker-vs-worker races")
PYY
