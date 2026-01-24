#' @title Utility Functions for STRIDER
#' @description Internal and exported utility functions for the STRIDER package
#' @name utils
#' @keywords internal
NULL

#' Default Color Palette
#'
#' @description A carefully curated color palette for visualizing cell types
#'   and clusters. Designed to be colorblind-friendly and distinguishable.
#'
#' @return A character vector of 36 hex color codes
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

#' Generate Color Palette
#'
#' @description Generate a color palette with a specified number of colors.
#'
#' @param n Number of colors needed
#' @param palette Base palette: "default" or "viridis"
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
    }
    # Extend with viridis
    extra_n <- n - length(default_cols)
    extra_cols <- viridis::viridis(extra_n)
    return(c(default_cols, extra_cols))
  } else if (palette == "viridis") {
    return(viridis::viridis(n))
  }
  rlang::abort(paste0("Unknown palette: ", palette))
}

#' Check if Object is Sparse Matrix
#'
#' @param x Object to check
#' @return Logical
#' @keywords internal
is_sparse <- function(x) {
  inherits(x, c("sparseMatrix", "dgCMatrix", "dgRMatrix", "dgTMatrix"))
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
#' @param x Numeric vector or matrix
#' @param base Log base (default: e)
#' @param pseudocount Pseudocount (default: 1)
#' @return Log-transformed values
#' @keywords internal
safe_log <- function(x, base = exp(1), pseudocount = 1) {
  log(x + pseudocount, base = base)
}

#' Normalize Rows/Columns to Sum to 1
#'
#' @param x Numeric vector or matrix
#' @param margin For matrices: 1 = rows, 2 = columns
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
  }
  
  col_sums <- Matrix::colSums(x)
  col_sums[col_sums == 0] <- 1
  return(Matrix::t(Matrix::t(x) / col_sums))
}

#' Print Message with Timestamp
#'
#' @param ... Message components
#' @param verbose Whether to print (default: from options)
#' @keywords internal
vmessage <- function(..., verbose = getOption("strider.verbose", TRUE)) {
  if (verbose) {
    timestamp <- format(Sys.time(), "[%H:%M:%S]")
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

#' Scale Sparse Matrix by Column
#'
#' @param x Sparse matrix
#' @param scale_factor Scale factor
#' @return Scaled matrix
#' @keywords internal
sparse_row_scale <- function(x, scale_factor) {
  if (!is_sparse(x)) {
    x <- as_sparse(x)
  }
  
  col_sums <- Matrix::colSums(x)
  col_sums[col_sums == 0] <- 1
  
  x %*% Matrix::Diagonal(x = scale_factor / col_sums)
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
      "Required packages not installed: ",
      paste(missing, collapse = ", "),
      "\nInstall with: install.packages(c(",
      paste0('"', missing, '"', collapse = ", "), "))"
    )
    if (error) rlang::abort(msg)
    rlang::warn(msg)
    return(FALSE)
  }
  TRUE
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
      rlang::abort(paste0("`", name, "` must be ", type, ", not ", class(x)[1]))
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

#' Safe Division
#'
#' @param x Numerator
#' @param y Denominator
#' @param replace Value when y is zero (default: 0)
#' @return Result of division
#' @keywords internal
safe_divide <- function(x, y, replace = 0) {
  result <- x / y
  result[!is.finite(result)] <- replace
  result
}

#' Convert STRIDER Object to Seurat
#'
#' @description Add STRIDER deconvolution results to a Seurat object.
#'
#' @param strider STRIDER result object
#' @param seurat Seurat object
#' @param assay.name Name for new assay (default: "STRIDER")
#'
#' @return Seurat object with deconvolution results
#' @export
as_seurat <- function(strider, seurat, assay.name = "STRIDER") {
  if (!inherits(strider, "STRIDER")) {
    rlang::abort("strider must be a STRIDER object")
  }
  if (!inherits(seurat, "Seurat")) {
    rlang::abort("seurat must be a Seurat object")
  }
  
  # Match spots
  common <- intersect(rownames(strider$deconv), colnames(seurat))
  if (length(common) == 0) {
    rlang::abort("No matching spots between STRIDER and Seurat objects")
  }
  
  # Add to metadata
  deconv_df <- as.data.frame(strider$deconv[common, , drop = FALSE])
  for (ct in colnames(deconv_df)) {
    seurat@meta.data[[paste0("strider_", ct)]] <- deconv_df[[ct]]
  }
  
  seurat
}

#' Create STRIDER Object from Seurat
#'
#' @description Extract data from Seurat object for STRIDER analysis.
#'
#' @param seurat Seurat object
#' @param celltype.col Column name in metadata containing cell types
#' @param assay Assay to use (default: "RNA")
#'
#' @return List with matrix and celltypes
#' @export
from_seurat <- function(seurat, celltype.col, assay = "RNA") {
  if (!inherits(seurat, "Seurat")) {
    rlang::abort("seurat must be a Seurat object")
  }
  
  # Extract counts
  counts <- .extract_seurat_counts(seurat)
  
  # Extract celltypes
  celltypes <- .extract_seurat_metadata(seurat, celltype.col)
  
  list(
    matrix = counts,
    celltypes = celltypes
  )
}

# Re-export the internal function for as_seurat
#' @keywords internal
.extract_seurat_counts <- function(seurat_obj) {
  counts <- NULL
  
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    tryCatch({
      counts <- SeuratObject::GetAssayData(seurat_obj, layer = "counts")
    }, error = function(e) {
      tryCatch({
        counts <<- SeuratObject::GetAssayData(seurat_obj, slot = "counts")
      }, error = function(e2) NULL)
    })
  }
  
  if (is.null(counts)) {
    tryCatch({
      counts <- seurat_obj@assays$RNA@counts
    }, error = function(e) {
      tryCatch({
        counts <- seurat_obj@assays$Spatial@counts
      }, error = function(e2) {
        rlang::abort("Could not extract counts from Seurat object")
      })
    })
  }
  
  as_sparse(counts)
}

#' @keywords internal
.extract_seurat_metadata <- function(seurat_obj, column) {
  if (!column %in% colnames(seurat_obj@meta.data)) {
    rlang::abort(paste0("Column '", column, "' not found in Seurat metadata"))
  }
  
  celltypes <- seurat_obj@meta.data[[column]]
  names(celltypes) <- rownames(seurat_obj@meta.data)
  celltypes
}
