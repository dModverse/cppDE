#!/usr/bin/env Rscript
# Render demo/<name>.R -> examples/<name>.md, figures in examples/figures/<name>/.
#
# The demos run against the installed package, so install it first. examples/ is
# .Rbuildignore'd and the README links the rendered pages. Top-level comments
# become prose, `# Title ----` rulers become headings. Arguments name the demos
# to render, by default all of them.
#
# Workflow: edit demo/<name>.R, install the package, run this script from the
# package root, commit the regenerated examples/<name>.md and its figures.

if (!dir.exists("demo") || !dir.exists("examples"))
  stop("Run this from the package root.")

demos <- sub("\\.R$", "", list.files("demo", pattern = "\\.R$"))
if (length(commandArgs(TRUE))) demos <- intersect(demos, commandArgs(TRUE))

# a top-level comment block after a blank line that opens with a capital is
# prose; a lowercase one or one right below code stays a code comment
.spinnable <- function(src) {
  src <- sub("^# (.*\\S) -{3,}$", "#' ## \\1", src, perl = TRUE)
  src <- sub("^## --- (.*) ---$", "#' ### \\1", src, perl = TRUE)
  prose <- FALSE
  for (i in seq_along(src)) {
    above <- if (i == 1) "" else src[i - 1]
    opens <- !nzchar(trimws(above)) && grepl("^# [A-Z`]", src[i])
    if (grepl("^#( |$)", src[i]) && (opens || prose))
      src[i] <- sub("^#( |$)", "#'\\1", src[i])
    prose <- grepl("^#'", src[i])
  }
  src
}

for (d in demos) {
  rmd <- file.path(tempdir(), paste0(d, ".Rmd"))
  src <- file.path(tempdir(), paste0(d, ".R"))
  head <- c(sprintf("#' # `demo(\"%s\")`", d), "#'",
            sprintf("#' Rendered by `dev/render-demos.R`; run it with `demo(\"%s\", package = \"cppDE\")`.", d),
            "#'", "")
  writeLines(c(head, .spinnable(readLines(file.path("demo", paste0(d, ".R"))))), src)
  knitr::spin(src, knit = FALSE, format = "Rmd")
  unlink(file.path("examples", "figures", d), recursive = TRUE)
  knitr::opts_knit$set(root.dir = normalizePath("examples"))
  knitr::opts_chunk$set(fig.path = file.path("figures", d, ""), fig.width = 10,
                        fig.height = 6.5, dpi = 96, comment = "#>", message = FALSE,
                        warning = FALSE, error = FALSE)
  owd <- setwd("examples")
  knitr::knit(rmd, output = paste0(d, ".md"), envir = new.env(), quiet = TRUE)
  setwd(owd)
  cat("rendered", d, "\n")
}
