#' @title Preprocessing Functions for STRIDER
#' @description Functions for preprocessing single-cell and spatial transcriptomics data
#' @name preprocess
NULL

#' Filter Genes
#'
#' @description Filter genes based on minimum expression criteria.
#'
#' @param matrix Count matrix (genes x cells/spots)
#' @param min.cells Minimum number of cells expressing the gene (default: 10)
#' @param min.counts Minimum total counts for a gene (default: 0)
#'
#' @return Filtered matrix
#' @export
#'
#' @examples
#' \dontrun{
#' filtered <- filter_genes(mat, min.cells = 10)
#' }
filter_genes <- function(matrix, min.cells = 10, min.counts = 0) {
  vmessage("Filtering genes...")
  
  n_genes_before <- nrow(matrix)
  
  # Calculate gene statistics
  if (is_sparse(matrix)) {
    n_cells_per_gene <- Matrix::rowSums(matrix > 0)
    total_counts_per_gene <- Matrix::rowSums(matrix)
  } else {
    n_cells_per_gene <- rowSums(matrix > 0)
    total_counts_per_gene <- rowSums(matrix)
  }
  
  # Filter
  keep <- (n_cells_per_gene >= min.cells) & (total_counts_per_gene >= min.counts)
  matrix <- matrix[keep, , drop = FALSE]
  
  n_genes_after <- nrow(matrix)
  vmessage("  Kept ", n_genes_after, " of ", n_genes_before, " genes")
  
  matrix
}

#' Filter Cells/Spots
#'
#' @description Filter cells or spots based on minimum expression criteria.
#'
#' @param matrix Count matrix (genes x cells/spots)
#' @param min.genes Minimum number of genes detected (default: 200)
#' @param min.counts Minimum total counts (default: 0)
#' @param max.counts Maximum total counts (default: Inf)
#'
#' @return Filtered matrix
#' @export
#'
#' @examples
#' \dontrun{
#' filtered <- filter_cells(mat, min.genes = 200)
#' }
filter_cells <- function(matrix, min.genes = 200, min.counts = 0, max.counts = Inf) {
  vmessage("Filtering cells/spots...")
  
  n_cells_before <- ncol(matrix)
  
  # Calculate cell statistics
  if (is_sparse(matrix)) {
    n_genes_per_cell <- Matrix::colSums(matrix > 0)
    total_counts_per_cell <- Matrix::colSums(matrix)
  } else {
    n_genes_per_cell <- colSums(matrix > 0)
    total_counts_per_cell <- colSums(matrix)
  }
  
  # Filter
  keep <- (n_genes_per_cell >= min.genes) & 
          (total_counts_per_cell >= min.counts) &
          (total_counts_per_cell <= max.counts)
  matrix <- matrix[, keep, drop = FALSE]
  
  n_cells_after <- ncol(matrix)
  vmessage("  Kept ", n_cells_after, " of ", n_cells_before, " cells/spots")
  
  matrix
}

#' Normalize Counts
#'
#' @description Normalize count matrix by total counts per cell/spot.
#'
#' @param matrix Count matrix (genes x cells/spots)
#' @param scale.factor Target sum for normalization (default: NULL, auto-detect)
#' @param log.transform Apply log1p transformation (default: FALSE)
#'
#' @return Normalized matrix
#' @export
#'
#' @examples
#' \dontrun{
#' normalized <- normalize_counts(mat, scale.factor = 10000)
#' }
normalize_counts <- function(matrix, scale.factor = NULL, log.transform = FALSE) {
  vmessage("Normalizing counts...")
  
  # Calculate total counts per cell
  if (is_sparse(matrix)) {
    total_counts <- Matrix::colSums(matrix)
  } else {
    total_counts <- colSums(matrix)
  }
  
  # Auto-detect scale factor if not provided
  if (is.null(scale.factor)) {
    scale.factor <- quantile_norm_factor(total_counts, q = 0.75)
    if (scale.factor == 0) scale.factor <- 10000
    vmessage("  Using scale factor: ", scale.factor)
  }
  
  # Normalize
  # Avoid division by zero
  total_counts[total_counts == 0] <- 1
  
  if (is_sparse(matrix)) {
    # Efficient sparse matrix normalization
    norm_factors <- scale.factor / total_counts
    matrix <- Matrix::t(Matrix::t(matrix) * norm_factors)
  } else {
    matrix <- sweep(matrix, 2, total_counts, "/") * scale.factor
  }
  
  # Log transform if requested
  if (log.transform) {
    vmessage("  Applying log1p transformation")
    matrix <- log1p(matrix)
  }
  
  matrix
}

#' Scale Data by Standard Deviation
#'
#' @description Scale data by dividing by standard deviation (without centering).
#'
#' @param matrix Normalized matrix (genes x cells/spots)
#' @param center Center data (subtract mean) (default: FALSE)
#' @param scale Scale data (divide by SD) (default: TRUE)
#'
#' @return Scaled matrix
#' @export
#'
#' @examples
#' \dontrun{
#' scaled <- scale_data(mat)
#' }
scale_data <- function(matrix, center = FALSE, scale = TRUE) {
  vmessage("Scaling data...")
  
  if (!scale && !center) {
    return(matrix)
  }
  
  # Convert sparse to dense for scaling (necessary for proper computation)
  was_sparse <- is_sparse(matrix)
  if (was_sparse && center) {
    vmessage("  Converting to dense matrix for centered scaling")
    matrix <- as.matrix(matrix)
  }
  
  if (scale) {
    # Calculate row-wise SD
    if (is_sparse(matrix)) {
      # For sparse matrices, compute SD without converting to dense
      # SD = sqrt(E[X^2] - E[X]^2)
      n <- ncol(matrix)
      row_means <- Matrix::rowMeans(matrix)
      row_sq_means <- Matrix::rowMeans(matrix^2)
      row_sds <- sqrt(row_sq_means - row_means^2)
      row_sds <- sqrt(row_sds^2 * n / (n - 1))  # Correct for sample SD
    } else {
      row_sds <- apply(matrix, 1, stats::sd)
    }
    
    # Avoid division by zero
    row_sds[row_sds == 0] <- 1
    
    if (center) {
      row_means <- if(is_sparse(matrix)) Matrix::rowMeans(matrix) else rowMeans(matrix)
      matrix <- (matrix - row_means) / row_sds
    } else {
      if (is_sparse(matrix)) {
        # Efficient sparse scaling
        matrix <- matrix / row_sds
      } else {
        matrix <- sweep(matrix, 1, row_sds, "/")
      }
    }
  } else if (center) {
    row_means <- if(is_sparse(matrix)) Matrix::rowMeans(matrix) else rowMeans(matrix)
    matrix <- matrix - row_means
  }
  
  matrix
}

#' Calculate Quality Control Metrics
#'
#' @description Calculate quality control metrics for cells/spots.
#'
#' @param matrix Count matrix (genes x cells/spots)
#' @param mito.pattern Pattern to match mitochondrial genes (default: "^MT-|^mt-")
#'
#' @return Data frame with QC metrics
#' @export
#'
#' @examples
#' \dontrun{
#' qc <- calculate_qc_metrics(mat)
#' }
calculate_qc_metrics <- function(matrix, mito.pattern = "^MT-|^mt-") {
  vmessage("Calculating QC metrics...")
  
  if (is_sparse(matrix)) {
    total_counts <- Matrix::colSums(matrix)
    n_genes <- Matrix::colSums(matrix > 0)
  } else {
    total_counts <- colSums(matrix)
    n_genes <- colSums(matrix > 0)
  }
  
  # Calculate mitochondrial percentage
  mito_genes <- grep(mito.pattern, rownames(matrix), value = TRUE)
  if (length(mito_genes) > 0) {
    if (is_sparse(matrix)) {
      mito_counts <- Matrix::colSums(matrix[mito_genes, , drop = FALSE])
    } else {
      mito_counts <- colSums(matrix[mito_genes, , drop = FALSE])
    }
    pct_mito <- mito_counts / total_counts * 100
  } else {
    pct_mito <- rep(0, ncol(matrix))
  }
  
  data.frame(
    barcode = colnames(matrix),
    nCount = total_counts,
    nFeature = n_genes,
    pct_mito = pct_mito,
    row.names = colnames(matrix)
  )
}

#' Preprocess Single-Cell Data
#'
#' @description Complete preprocessing pipeline for single-cell RNA-seq data.
#'
#' @param matrix Count matrix (genes x cells)
#' @param celltypes Named vector of cell type annotations
#' @param scale.factor Target sum for normalization (default: NULL)
#' @param min.cells Minimum cells per gene (default: 10)
#' @param min.genes Minimum genes per cell (default: 200)
#' @param normalize Whether to normalize (default: TRUE)
#' @param scale Whether to scale by SD (default: FALSE)
#'
#' @return A list containing preprocessed data
#' @export
#'
#' @examples
#' \dontrun{
#' sc_data <- preprocess_sc(mat, celltypes)
#' }
preprocess_sc <- function(matrix, celltypes, scale.factor = NULL,
                          min.cells = 10, min.genes = 200,
                          normalize = TRUE, scale = FALSE) {
  
  vmessage("Preprocessing single-cell data...")
  
  # Ensure sparse
  matrix <- as_sparse(matrix)
  
  # Filter cells present in annotations
  common_cells <- intersect(colnames(matrix), names(celltypes))
  if (length(common_cells) < ncol(matrix)) {
    vmessage("  Subsetting to ", length(common_cells), " annotated cells")
    matrix <- matrix[, common_cells]
    celltypes <- celltypes[common_cells]
  }
  
  # Store raw matrix
  raw_matrix <- matrix
  
  # Filter genes
  matrix <- filter_genes(matrix, min.cells = min.cells)
  
  # Filter cells
  matrix <- filter_cells(matrix, min.genes = min.genes)
  
  # Update celltypes to match filtered cells
  celltypes <- celltypes[colnames(matrix)]
  raw_matrix <- raw_matrix[, colnames(matrix)]
  
  # Normalize
  if (normalize) {
    norm_matrix <- normalize_counts(matrix, scale.factor = scale.factor)
  } else {
    norm_matrix <- matrix
  }
  
  # Scale if requested
  if (scale) {
    scale_matrix <- scale_data(norm_matrix, center = FALSE, scale = TRUE)
  } else {
    scale_matrix <- norm_matrix
  }
  
  list(
    scale_matrix = scale_matrix,
    raw_matrix = raw_matrix[rownames(scale_matrix), ],
    genes = rownames(scale_matrix),
    cells = colnames(scale_matrix),
    celltypes = celltypes
  )
}

#' Preprocess Spatial Data
#'
#' @description Complete preprocessing pipeline for spatial transcriptomics data.
#'
#' @param matrix Count matrix (genes x spots)
#' @param scale.factor Target sum for normalization (default: NULL)
#' @param normalize Whether to normalize (default: TRUE)
#' @param scale Whether to scale by SD (default: FALSE)
#'
#' @return A list containing preprocessed data
#' @export
#'
#' @examples
#' \dontrun{
#' st_data <- preprocess_st(mat)
#' }
preprocess_st <- function(matrix, scale.factor = NULL, 
                          normalize = TRUE, scale = FALSE) {
  
  vmessage("Preprocessing spatial data...")
  
  # Ensure sparse
  matrix <- as_sparse(matrix)
  
  # Store raw matrix
  raw_matrix <- matrix
  
  # Normalize
  if (normalize) {
    scale_matrix <- normalize_counts(matrix, scale.factor = scale.factor)
  } else {
    scale_matrix <- matrix
  }
  
  # Scale if requested
  if (scale) {
    scale_matrix <- scale_data(scale_matrix, center = FALSE, scale = TRUE)
  }
  
  list(
    scale_matrix = scale_matrix,
    raw_matrix = raw_matrix,
    genes = rownames(scale_matrix),
    spots = colnames(scale_matrix)
  )
}
