#!/usr/bin/env Rscript

# Local R CMD check mimicking the GitHub Actions check-r-package step
# (build_args + --as-cran). Aborts on warnings, not just on errors.

# Needs vignettes/Methods.pdf (else run dev/render-methods.R first). Without
# ghostscript, --compact-vignettes=gs+qpdf skips compaction instead of failing.

# On Windows without pandoc on PATH: try the usual RStudio bundle locations.
if (!nzchar(Sys.which("pandoc")) && !nzchar(Sys.getenv("RSTUDIO_PANDOC"))) {
  rstudio_pandoc_candidates <- c(
    "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
    "C:/Program Files/RStudio/bin/pandoc",
    "C:/Program Files/Quarto/bin/tools"
  )
  hit <- rstudio_pandoc_candidates[file.exists(file.path(rstudio_pandoc_candidates, "pandoc.exe"))][1]
  if (!is.na(hit)) {
    Sys.setenv(RSTUDIO_PANDOC = hit)
    message("Setting RSTUDIO_PANDOC = ", hit)
  }
}

# Tests run in parallel (Config/testthat/parallel in DESCRIPTION). testthat uses
# min(TESTTHAT_CPUS, number of test files) workers, so this is an upper bound.
# Each worker shells out to a compiler, so the ceiling is the slowest file.
if (!nzchar(Sys.getenv("TESTTHAT_CPUS"))) {
  Sys.setenv(TESTTHAT_CPUS = max(1L, parallel::detectCores() - 2L))
}

devtools::check(
  error_on   = "warning",
  args       = c("--no-manual", "--as-cran"),
  build_args = c("--no-manual", "--compact-vignettes=gs+qpdf")
)
