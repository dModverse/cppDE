# compile() on hand-written sources: what reaches the object compile of a
# combined output, and how a failed compile ends.

src_dir <- function() {
  d <- tempfile("cppDE_compile_")
  dir.create(d)
  d
}

# A compile() handle for a source of R entry points.
handle <- function(dir, name, body, compileArgs = NULL) {
  f <- file.path(dir, paste0(name, ".cpp"))
  writeLines(c("#include <R.h>", "#include <Rinternals.h>", body), f)
  structure(list(), srcfile = f, compileArgs = compileArgs)
}

test_that("a user Makevars reaches the object compile of a combined output", {
  d <- src_dir()
  inc <- file.path(d, "inc")
  dir.create(inc)
  writeLines("#define CPPDE_TEST_HEADER 3", file.path(inc, "cppde_test_header.h"))
  mv <- file.path(d, "user.mk")
  writeLines(paste0("CPPFLAGS += -I", normalizePath(inc, winslash = "/"), " -DCPPDE_TEST_MACRO=4"), mv)
  withr::local_envvar(R_MAKEVARS_USER = mv)

  obj <- handle(d, "cppde_test_mv", c(
    '#include "cppde_test_header.h"',
    'extern "C" SEXP cppde_test_mv() {',
    '  return Rf_ScalarInteger(CPPDE_TEST_HEADER * CPPDE_TEST_MACRO);',
    '}'))
  compile(obj, output = "cppde_test_mv_out")
  expect_identical(.Call("cppde_test_mv", PACKAGE = "cppde_test_mv_out"), 12L)
})

test_that("the compileArgs of a source reach that source only", {
  d <- src_dir()
  body <- function(own, other) c(
    sprintf("#ifdef CPPDE_TEST_%s", other), "#error foreign compileArgs", "#endif",
    sprintf('extern "C" SEXP cppde_test_%s() {', tolower(own)),
    sprintf("  return Rf_ScalarInteger(CPPDE_TEST_%s);", own), "}")
  a <- handle(d, "cppde_test_a", body("A", "B"), "-DCPPDE_TEST_A=1")
  b <- handle(d, "cppde_test_b", body("B", "A"), "-DCPPDE_TEST_B=2")
  compile(a, b, output = "cppde_test_ab_out", cores = 2)
  expect_identical(.Call("cppde_test_a", PACKAGE = "cppde_test_ab_out"), 1L)
  expect_identical(.Call("cppde_test_b", PACKAGE = "cppde_test_ab_out"), 2L)
})

test_that("compile() rebuilds the object of a source compiled before", {
  d <- src_dir()
  obj <- handle(d, "cppde_test_re", c(
    'extern "C" SEXP cppde_test_re() { return Rf_ScalarInteger(CPPDE_TEST_V); }'),
    "-DCPPDE_TEST_V=1")
  compile(obj, output = "cppde_test_re_out1")
  attr(obj, "compileArgs") <- "-DCPPDE_TEST_V=2"
  compile(obj, output = "cppde_test_re_out2")
  expect_identical(.Call("cppde_test_re", PACKAGE = "cppde_test_re_out2"), 2L)
})

test_that("a source that fails to compile is an error naming it, also under cores = 2", {
  d <- src_dir()
  good <- handle(d, "cppde_test_good",
                 'extern "C" SEXP cppde_test_good() { return R_NilValue; }')
  bad <- handle(d, "cppde_test_bad",
                'extern "C" SEXP cppde_test_bad() { return cppde_undeclared; }')
  warned <- FALSE
  expect_error(
    withCallingHandlers(
      compile(good, bad, output = "cppde_test_fail_out", cores = 2),
      warning = function(w) {
        warned <<- TRUE
        invokeRestart("muffleWarning")
      }),
    "cppde_test_bad\\.cpp[\\s\\S]*cppde_undeclared", perl = TRUE)
  expect_false(warned)
})
