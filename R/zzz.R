#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib STRIDER, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

.onLoad <- function(libname, pkgname) {
  # Set default options
  op <- options()
  op.strider <- list(
    strider.verbose = TRUE,
    strider.seed = 42,
    strider.ncores = 1L
  )
  toset <- !(names(op.strider) %in% names(op))
  if (any(toset)) options(op.strider[toset])
  
  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "STRIDER v", utils::packageVersion("STRIDER"), "\n",
    "Spatial TRanscrIptomics DEconvolution in R\n",
    "Author: Zaoqu Liu (liuzaoqu@163.com)\n",
    "GitHub: https://github.com/Zaoqu-Liu/STRIDER"
  )
}
