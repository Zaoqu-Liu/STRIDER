#' @title Utility Functions for STRIDER
#' @description Internal utility functions for the STRIDER package
#' @name utils
#' @keywords internal
NULL

#' Default Color Palette for STRIDER
#'
#' @description A carefully curated color palette for visualizing cell types
#'   and clusters in spatial transcriptomics data.
#'
#' @return A character vector of hex color codes
#' @export
#'
#' @examples
#' colors <- get_default_colors()
#' scales::show_col(colors)
get_default_colors <- function() {
  c(
    "#E58606", "#5D69B1", "#52BCA3", "#99C945", "#CC61B0", "#24796C",
    "#DAA51B", "#2F8AC4", "#764E9F", "#ED645A", "#9EB4C7", "#BCBD22",
    "#B279A2", "#EECA3B", "#17BECF", "#FF9DA6", "#778AAE", "#1B9E77",
    "#A6761D", "#526A83", "#B82E2E", "#80B1D3", "#68855C", "#D95F02",
    "#BEBADA", "#AF6458", "#D9AF6B", "#9C9C5E", "#625377", "#8C785D",
    "#88CCEE", "#E73F74", "#FFFFB3", "#CCEBC5", "#332288", "#A65628"
  )
}

#' Extended Color Palette Generator
#'
#' @description Generate a color palette with a specified number of colors.
#'   Uses default colors first, then extends with viridis if needed.
#'
#' @param n Number of colors needed
#' @param palette Base palette to use ("default", "viridis", "spectral")
#'
#' @return A character vector of hex color codes
#' @export
#'
#' @examples
#' get_extended_colors(10)
#' get_extended_colors(50, palette = "viridis")
get_extended_colors <- function(n, palette = "default") {
  if (palette == "default") {
    default_cols <- get_default_colors()
    if (n <= length(default_cols)) {
      return(default_cols[seq_len(n)])
    } else {
      # Extend with viridis
      extra_n <- n - length(default_cols)
      extra_cols <- viridis::viridis(extra_n)
      return(c(default_cols, extra_cols))
    }
  } else if (palette == "viridis") {
    return(viridis::viridis(n))
  } else if (palette == "spectral") {
    return(grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(11, "Spectral")
    )(n))
  } else {
    rlang::abort(paste0("Unknown palette: ", palette))
  }
}

#' Check if Object is Sparse Matrix
#'
#' @param x Object to check
#' @return Logical
#' @keywords internal
is_sparse <- function(x) {
  inherits(x, "sparseMatrix") || inherits(x, "dgCMatrix") || 
    inherits(x, "dgRMatrix") || inherits(x, "dgTMatrix")
}

#' Convert to Sparse Matrix
#'
#' @param x Matrix to convert
#' @return dgCMatrix
#' @keywords internal
as_sparse <- function(x) {

  if (is_sparse(x)) {
    return(methods::as(x, "dgCMatrix"))
  }
  Matrix::Matrix(x, sparse = TRUE)
}

#' Safe Log Transformation
#'
#' @description Perform log transformation with a small offset to avoid log(0)
#'
#' @param x Numeric vector or matrix
#' @param base Log base (default: natural log)
#' @param pseudocount Pseudocount to add (default: 1)
#'
#' @return Log-transformed values
#' @keywords internal
safe_log <- function(x, base = exp(1), pseudocount = 1) {
  log(x + pseudocount, base = base)
}

#' Normalize to Sum to 1
#'
#' @description Normalize a vector or matrix rows to sum to 1
#'
#' @param x Numeric vector or matrix
#' @param margin For matrices: 1 = rows, 2 = columns
#'
#' @return Normalized values
#' @keywords internal
normalize_to_one <- function(x, margin = 1) {
  if (is.vector(x)) {
    x_sum <- sum(x)
    if (x_sum == 0) return(x)
    return(x / x_sum)
  }
  
  if (margin == 1) {
    row_sums <- Matrix::rowSums(x)
    row_sums[row_sums == 0] <- 1
    return(x / row_sums)
  } else {
    col_sums <- Matrix::colSums(x)
    col_sums[col_sums == 0] <- 1
    return(Matrix::t(Matrix::t(x) / col_sums))
  }
}

#' Message with Timestamp
#'
#' @param ... Message components
#' @param verbose Whether to print
#' @keywords internal
vmessage <- function(..., verbose = getOption("strider.verbose", TRUE)) {
  if (verbose) {
    timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
    message(timestamp, " ", ...)
  }
}

#' Calculate Quantile Normalization Factor
#'
#' @param x Numeric vector
#' @param q Quantile (default: 0.75)
#' @return Normalization factor
#' @keywords internal
quantile_norm_factor <- function(x, q = 0.75) {
  qval <- stats::quantile(x, probs = q, na.rm = TRUE)
  round(qval / 1000, 0) * 1000
}

#' Fast Row/Column Operations on Sparse Matrices
#'
#' @name sparse_ops
#' @keywords internal
NULL

#' @rdname sparse_ops
sparse_row_scale <- function(x, scale_factor) {
  # Efficient row scaling for sparse matrices
  if (!is_sparse(x)) {
    x <- as_sparse(x)
  }
  
  # Get dimensions
  nr <- nrow(x)
  nc <- ncol(x)
  
  # Column sums for normalization
  col_sums <- Matrix::colSums(x)
  col_sums[col_sums == 0] <- 1
  
  # Scale each column
  x_scaled <- x %*% Matrix::Diagonal(x = scale_factor / col_sums)
  
  return(x_scaled)
}

#' @rdname sparse_ops
sparse_col_scale <- function(x, scale_factor) {
  Matrix::t(sparse_row_scale(Matrix::t(x), scale_factor))
}

#' Check Required Packages
#'
#' @param pkgs Character vector of package names
#' @param error If TRUE, throw error; if FALSE, return logical
#' @return Logical indicating if all packages are available
#' @keywords internal
check_packages <- function(pkgs, error = TRUE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  
  if (length(missing) > 0) {
    msg <- paste0(
      "The following packages are required but not installed: ",
      paste(missing, collapse = ", "),
      "\nPlease install them with: install.packages(c(",
      paste0('"', missing, '"', collapse = ", "), "))"
    )
    if (error) {
      rlang::abort(msg)
    } else {
      rlang::warn(msg)
      return(FALSE)
    }
  }
  return(TRUE)
}

#' Validate Input Data
#'
#' @param x Data to validate
#' @param name Name for error messages
#' @param type Expected type
#' @param allow_null Allow NULL values
#' @keywords internal
validate_input <- function(x, name, type = NULL, allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) return(invisible(NULL))
    rlang::abort(paste0("`", name, "` cannot be NULL"))
  }
  
  if (!is.null(type)) {
    type_check <- switch(
      type,
      "matrix" = is.matrix(x) || is_sparse(x),
      "data.frame" = is.data.frame(x) || data.table::is.data.table(x),
      "numeric" = is.numeric(x),
      "character" = is.character(x),
      "logical" = is.logical(x),
      "list" = is.list(x),
      TRUE
    )
    
    if (!type_check) {
      rlang::abort(paste0(
        "`", name, "` must be of type ", type, 
        ", not ", class(x)[1]
      ))
    }
  }
  
  invisible(NULL)
}

#' Set Random Seed
#'
#' @param seed Seed value
#' @keywords internal
set_seed <- function(seed = NULL) {
  if (is.null(seed)) {
    seed <- getOption("strider.seed", 42)
  }
  set.seed(seed)
}

#' Create Progress Bar
#'
#' @param n Total iterations
#' @param name Task name
#' @return Progress bar object
#' @keywords internal
create_progress <- function(n, name = "Processing") {
  if (requireNamespace("progressr", quietly = TRUE)) {
    progressr::progressor(steps = n)
  } else {
    function(...) invisible(NULL)
  }
}

#' Safe Division
#'
#' @description Perform division while handling zero denominators
#'
#' @param x Numerator
#' @param y Denominator
#' @param replace Value to use when y is zero (default: 0)
#'
#' @return Result of division
#' @keywords internal
safe_divide <- function(x, y, replace = 0) {
  result <- x / y
  result[!is.finite(result)] <- replace
  result
}
