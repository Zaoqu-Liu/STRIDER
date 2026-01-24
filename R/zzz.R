#' @useDynLib STRIDER, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom utils head
NULL

# Suppress R CMD check notes for NSE variables used in ggplot2 and data.table
utils::globalVariables(c(
  # ggplot2 aesthetics
  "x", "y", "x_aligned", "y_aligned", 
  "celltype", "cluster", "fraction", "topic", "value",
  "spot", "n_topics", "accuracy", "method", "coherence", "nmi",
  # data.table operators
  ":="
))

.onLoad <- function(libname, pkgname) {
  # Set default options
  op <- options()
  op.strider <- list(
    strider.verbose = TRUE,
    strider.seed = 42
  )
  toset <- !(names(op.strider) %in% names(op))
  if (any(toset)) options(op.strider[toset])
  
  invisible()
}

.onAttach <- function(libname, pkgname) {
  version <- utils::packageVersion(pkgname)
  packageStartupMessage(
    "STRIDER version ", version, "\n",
    "Spatial Transcriptomics Deconvolution and Integration in R\n",
    "Author: Zaoqu Liu (liuzaoqu@163.com)\n",
    "GitHub: https://github.com/Zaoqu-Liu/STRIDER"
  )
}

.onUnload <- function(libpath) {
  library.dynam.unload("STRIDER", libpath)
}
