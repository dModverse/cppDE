## =====================================================================
##  diagnose-cvode.R -- why does the CVODE backend fail on *this* host?
##
##  Run it on the machine in question, over ssh and with the same R that
##  a runbg job would get:
##
##      ssh <host> 'R --vanilla --no-echo' < dev/diagnose-cvode.R
##
##  It prints everything static first and attempts a solve last, because
##  some failures (an MPI abort, a segfault in a linked library) kill the
##  process rather than raise an R condition: if the output stops before
##  the last section, the solve took the process down with it, and that
##  is itself the answer.
## =====================================================================

hr <- function(s) cat("\n== ", s, " ", strrep("=", max(0, 60 - nchar(s))), "\n", sep = "")
kv <- function(k, v) cat(sprintf("  %-22s %s\n", k, paste(v, collapse = " ")))

hr("host and R")
kv("node", Sys.info()[["nodename"]])
kv("R", paste(R.version$major, R.version$minor, sep = "."))
kv("platform", R.version$platform)
kv("R_HOME", R.home())
kv("BLAS_LIBS", system2(file.path(R.home("bin"), "R"), c("CMD", "config", "BLAS_LIBS"),
                        stdout = TRUE, stderr = FALSE))
kv("LAPACK_LIBS", system2(file.path(R.home("bin"), "R"), c("CMD", "config", "LAPACK_LIBS"),
                          stdout = TRUE, stderr = FALSE))

hr("cppDE configuration")
suppressPackageStartupMessages(library(cppDE))
kv("cppDE", as.character(utils::packageVersion("cppDE")))
kv("libPath", dirname(find.package("cppDE")))
for (k in c("available", "cflags", "libs", "klu_available", "klu_cflags",
            "klu_libs", "cvode_lapack_available", "cvode_lapack_libs"))
  kv(k, format(cppDE:::cvodeConfig[[k]]))

## The prefix cppDE was *installed against*, which is not necessarily
## the one install_libs() last built: the paths are frozen into
## cvodeConfig.dcf at install time.
hr("SUNDIALS actually linked")
prefix <- sub(".*-I *([^ ]*)/include.*", "\\1", cppDE:::cvodeConfig$cflags)
kv("prefix", prefix)
cfg <- file.path(prefix, "include", "sundials", "sundials_config.h")
if (file.exists(cfg)) {
  hits <- grep("MPI|PROFILING|VERSION ", readLines(cfg, warn = FALSE), value = TRUE)
  for (h in hits) cat("   ", h, "\n")
} else cat("    no sundials_config.h at that prefix\n")

for (lib in Sys.glob(file.path(prefix, "lib", "libsundials_c*"))) {
  ld <- suppressWarnings(system2("ldd", shQuote(lib), stdout = TRUE, stderr = TRUE))
  kv(basename(lib), c(grep("mpi", ld, ignore.case = TRUE, value = TRUE), "")[1L])
}

hr("compiling a one-state CVODE model")
m <- cvode(c(A = "-k*A"), modelname = "cvode_probe", compile = TRUE, verbose = TRUE)
so <- sub("\\.cpp$", .Platform$dynlib.ext, attr(m, "srcfile"))
kv("shared object", so)
## RPATH outranks LD_LIBRARY_PATH, RUNPATH does not.
dt <- suppressWarnings(system2("readelf", c("-d", shQuote(so)), stdout = TRUE,
                               stderr = FALSE))
kv("rpath kind", c(regmatches(dt, regexpr("RUNPATH|RPATH", dt))[1L], "unknown")[1L])
ld <- suppressWarnings(system2("ldd", shQuote(so), stdout = TRUE, stderr = TRUE))
for (l in grep("mpi|lapack|blas|sundials|klu", ld, ignore.case = TRUE, value = TRUE))
  cat("   ", trimws(l), "\n")

## The decisive evidence: what is mapped into this process now that the
## model is loaded.  A library that arrives indirectly -- pulled in by
## BLAS, by R itself, by an LD_PRELOAD -- shows up here and nowhere in
## the ldd output above.
hr("MPI mapped into the process?")
maps <- tryCatch(readLines("/proc/self/maps", warn = FALSE), error = function(e) character(0))
## Match the library, not the string: R's tempdir is named RtmpXXXXXX and
## a random suffix beginning with "I" makes "Rtmp..." match a
## case-insensitive "mpi".
mpi <- unique(sub(".*\\s", "", grep("libmpi", maps, value = TRUE)))
if (length(mpi)) for (p in mpi) cat("    ", p, "\n") else cat("    none\n")
kv("LD_PRELOAD", Sys.getenv("LD_PRELOAD", "<unset>"))
kv("LD_LIBRARY_PATH", Sys.getenv("LD_LIBRARY_PATH", "<unset>"))

hr("solving -- serial")
print(solveODE(m, seq(0, 5, 1), c(A = 1, k = 0.5))$variable)

## Only reached if the serial solve survived.  The benchmark runs its
## problems under mclapply, so a failure that appears only after a fork
## has to be told apart from one that was there all along.
hr("solving -- forked (mclapply, 2 workers)")
if (.Platform$OS.type == "unix") {
  r <- parallel::mclapply(1:2, function(i)
    tryCatch(solveODE(m, seq(0, 5, 1), c(A = 1, k = 0.5))$variable[2L, ],
             error = function(e) conditionMessage(e)), mc.cores = 2L)
  print(r)
} else cat("    not unix\n")

hr("done -- CVODE works on this host")
