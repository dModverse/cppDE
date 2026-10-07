# Compile jobs per test file: CPPDE_TEST_CORES, one by default.
test_cores <- function() {
  n <- suppressWarnings(as.integer(Sys.getenv("CPPDE_TEST_CORES", "1")))
  if (is.na(n) || n < 1L) 1L else n
}
