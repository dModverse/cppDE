# The pkgdown site with the reference topics and the articles built in parallel R processes, then
# the index pages, news and search. Run from the package root: Rscript dev/build-site.R [workers]
# The site is built from scratch: docs/ is removed first, since pkgdown keeps pages of removed
# articles.
workers <- as.integer(commandArgs(TRUE)[1])
if (is.na(workers)) workers <- 8L

unlink("docs", recursive = TRUE)
pkg <- pkgdown::as_pkgdown(".")
pkgdown::init_site(pkg)

# Articles that read the same example project compile into its build/ and run one after the other
# in one job; every job is a separate R process with a separate tempdir(). Reference topics go in
# `workers` groups; examples run against the installed package (devel = FALSE).
files <- pkg$vignettes$file_in
project <- vapply(files, function(f) {
  src <- readLines(f, warn = FALSE)
  hit <- regmatches(src, regexpr("examples/[A-Za-z0-9]+", src))
  if (length(hit)) hit[1] else NA_character_
}, "")
key <- ifelse(is.na(project), pkg$vignettes$name, project)
groups <- split(pkg$vignettes$name, key)
cost <- vapply(split(file.size(files), key), sum, 0)
groups <- groups[order(cost[names(groups)], decreasing = TRUE)]
topics <- split(pkg$topics$name, seq_along(pkg$topics$name) %% workers)
jobs <- c(lapply(groups, function(a) list(kind = "articles", what = a)),
          lapply(topics, function(t) list(kind = "reference", what = t)))

build <- function(job) {
  t0 <- Sys.time()
  err <- tryCatch({
    callr::r(function(kind, what) {
      pkg <- pkgdown::as_pkgdown(".")
      if (kind == "reference")
        pkgdown::build_reference(pkg, lazy = FALSE, preview = FALSE, devel = FALSE, topics = what)
      else for (a in what) pkgdown::build_article(a, pkg, new_process = FALSE, quiet = TRUE)
    }, args = list(kind = job$kind, what = job$what), show = FALSE)
    NA_character_
  }, error = function(e) conditionMessage(if (!is.null(e$parent)) e$parent else e))
  label <- if (job$kind == "reference") sprintf("reference (%d topics)", length(job$what))
           else paste(basename(job$what), collapse = ", ")
  data.frame(job = label, secs = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
             error = err)
}
done <- do.call(rbind, parallel::mclapply(jobs, build, mc.cores = workers, mc.preschedule = FALSE))
# pkgdown collects the assets of an article through one scratch file in the article directory, so
# two articles finishing at once can collide there; failed jobs run once more, one after the other
again <- which(!is.na(done$error))
for (i in again) done[i, ] <- build(jobs[[i]])
print(done[order(-done$secs), c("job", "secs")], row.names = FALSE)

pkgdown::build_reference_index(pkg)
pkgdown::build_articles_index(pkg)
pkgdown::build_home(pkg, preview = FALSE)
# release dates from CRAN need the network, which build hosts may not have
pkgdown::build_news(pkg, preview = FALSE, override = list(news = list(cran_dates = FALSE)))
pkgdown::build_redirects(pkg)
pkgdown::build_search(pkg)

failed <- done[!is.na(done$error), ]
if (nrow(failed)) {
  for (i in seq_len(nrow(failed))) message("failed: ", failed$job[i], "\n  ", failed$error[i])
  quit(status = 1)
}
