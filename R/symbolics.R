#' Extract Symbol Names from R Expressions
#'
#' Returns the unique names of all symbols occurring in a character
#' vector of R expressions.
#'
#' @param expr Character vector of R expressions.
#' @param omit Optional character vector of symbol names to remove.
#'
#' @return Character vector of unique symbol names.
#'
#' @keywords internal
getSymbols <- function(expr, omit = NULL) {
  if (is.null(expr)) return(character(0))

  expr <- expr[expr != "0"]
  if (!length(expr)) return(character(0))

  parsed <- tryCatch(parse(text = expr, keep.source = TRUE), error = function(e) NULL)
  if (is.null(parsed)) return(character(0))

  pd <- utils::getParseData(parsed)
  syms <- unique(pd[pd$token == "SYMBOL", "text"])

  if (!is.null(omit)) syms <- setdiff(syms, omit)
  syms
}
## Python 3 `keyword.kwlist`.
.pythonKeywords <- c(
  "False","None","True","and","as","assert","async","await","break","class",
  "continue","def","del","elif","else","except","finally","for","from","global",
  "if","import","in","is","lambda","nonlocal","not","or","pass","raise","return",
  "try","while","with","yield"
)

#' Reject Symbol Names the Backend Cannot Parse
#'
#' Stops when an expression or a symbol name holds a Python keyword.
#'
#' @param ... Character vectors of expressions or symbol names, optionally
#'   named. Names are checked along with the values. `NULL` is ignored.
#' @return `NULL`, invisibly. Called for the error.
#'
#' @keywords internal
checkSymbolNames <- function(...) {
  ## Python's parser rejects a keyword and reads True, False and None as
  ## constants.
  ## C++ keywords need no check: the printer emits a slot, not the name.
  text <- unlist(list(...), use.names = TRUE)
  if (!length(text)) return(invisible(NULL))
  text <- c(as.character(text), names(text))

  ## One pass over all expressions.
  pat <- paste0("\\b(?:", paste(.pythonKeywords, collapse = "|"), ")\\b")
  found <- regmatches(text, gregexpr(pat, text, perl = TRUE))
  hit <- intersect(.pythonKeywords, unlist(found))
  if (!length(hit)) return(invisible(NULL))

  stop("Python keyword used as a symbol name: ",
       paste0("'", hit, "'", collapse = ", "),
       ". The code generator reads the equations with Python's parser and ",
       "cannot read such a name as a symbol. Rename it in the model definition.",
       call. = FALSE)
}
