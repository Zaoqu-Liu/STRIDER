#' @title Marker Gene Identification for STRIDER
#' @description Functions for identifying cell-type specific marker genes
#' @name marker_find
NULL

#' Find Marker Genes for Each Cell Type
#'
#' @description Identify differentially expressed marker genes for each cell type
#'   using Wilcoxon rank-sum test. This is crucial for training the LDA model.
#'
#' @param matrix Count matrix (genes x cells), should be normalized
#' @param celltypes Named vector of cell type annotations
#' @param n.top Number of top markers per cell type (default: 200)
#' @param min.pct Minimum fraction of cells expressing the gene (default: 0.1)
#' @param logfc.threshold Minimum log fold change threshold (default: 0.25)
#' @param only.pos Only return positive markers (default: TRUE)
#' @param test.use Statistical test to use: "wilcox" or "t" (default: "wilcox")
#' @param adjust.method P-value adjustment method (default: "BH")
#' @param ncores Number of cores for parallel computation (default: 1)
#'
#' @return A list containing:
#'   \itemize{
#'     \item{markers}{Data frame of all markers with statistics}
#'     \item{top_markers}{Character vector of top marker genes}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' markers <- find_markers(norm_mat, celltypes, n.top = 200)
#' }
find_markers <- function(matrix, celltypes, n.top = 200, min.pct = 0.1,
                         logfc.threshold = 0.25, only.pos = TRUE,
                         test.use = "wilcox", adjust.method = "BH",
                         ncores = 1) {
  
  vmessage("Finding marker genes...")
  
  # Ensure celltypes is aligned with matrix columns
  common_cells <- intersect(colnames(matrix), names(celltypes))
  if (length(common_cells) == 0) {
    rlang::abort("No matching cells between matrix and celltypes")
  }
  
  matrix <- matrix[, common_cells]
  celltypes <- celltypes[common_cells]
  
  # Get unique cell types
  celltype_levels <- sort(unique(celltypes))
  n_celltypes <- length(celltype_levels)
  
  vmessage("  Testing ", nrow(matrix), " genes across ", n_celltypes, " cell types")
  
  # Pre-compute log-transformed data for fold change calculation
  if (is_sparse(matrix)) {
    log_data <- matrix
    log_data@x <- log1p(log_data@x)
  } else {
    log_data <- log1p(matrix)
  }
  
  # Function to test one cell type vs rest
  test_one_celltype <- function(ct) {
    cells_in <- names(celltypes)[celltypes == ct]
    cells_out <- names(celltypes)[celltypes != ct]
    
    n_in <- length(cells_in)
    n_out <- length(cells_out)
    
    results <- lapply(seq_len(nrow(matrix)), function(i) {
      gene <- rownames(matrix)[i]
      
      # Get expression values
      if (is_sparse(matrix)) {
        expr_in <- as.numeric(matrix[i, cells_in])
        expr_out <- as.numeric(matrix[i, cells_out])
        log_in <- as.numeric(log_data[i, cells_in])
        log_out <- as.numeric(log_data[i, cells_out])
      } else {
        expr_in <- matrix[i, cells_in]
        expr_out <- matrix[i, cells_out]
        log_in <- log_data[i, cells_in]
        log_out <- log_data[i, cells_out]
      }
      
      # Calculate percentage expressing
      pct_in <- sum(expr_in > 0) / n_in
      pct_out <- sum(expr_out > 0) / n_out
      pct_diff <- pct_in - pct_out
      
      # Skip if below minimum percentage
      if (pct_in < min.pct && pct_out < min.pct) {
        return(NULL)
      }
      
      # Calculate log fold change
      mean_in <- mean(log_in)
      mean_out <- mean(log_out)
      avg_logfc <- mean_in - mean_out
      
      # Skip if below fold change threshold
      if (only.pos && avg_logfc < logfc.threshold) {
        return(NULL)
      }
      if (!only.pos && abs(avg_logfc) < logfc.threshold) {
        return(NULL)
      }
      
      # Statistical test
      if (test.use == "wilcox") {
        # Wilcoxon rank-sum test
        tryCatch({
          test_result <- stats::wilcox.test(
            expr_in, expr_out,
            alternative = if (only.pos) "greater" else "two.sided",
            exact = FALSE
          )
          pval <- test_result$p.value
        }, error = function(e) {
          pval <<- 1
        })
      } else if (test.use == "t") {
        # t-test
        tryCatch({
          test_result <- stats::t.test(
            expr_in, expr_out,
            alternative = if (only.pos) "greater" else "two.sided"
          )
          pval <- test_result$p.value
        }, error = function(e) {
          pval <<- 1
        })
      }
      
      data.frame(
        gene = gene,
        cluster = ct,
        avg_logfc = avg_logfc,
        pct.1 = pct_in,
        pct.2 = pct_out,
        pct_diff = pct_diff,
        p_val = pval,
        stringsAsFactors = FALSE
      )
    })
    
    do.call(rbind, results)
  }
  
  # Run tests for all cell types
  if (ncores > 1 && requireNamespace("future.apply", quietly = TRUE)) {
    # Parallel execution
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    
    future::plan(future::multisession, workers = ncores)
    all_results <- future.apply::future_lapply(
      celltype_levels,
      test_one_celltype,
      future.seed = TRUE
    )
  } else {
    # Sequential execution
    all_results <- lapply(celltype_levels, test_one_celltype)
  }
  
  # Combine results
  markers_df <- do.call(rbind, all_results)
  
  if (is.null(markers_df) || nrow(markers_df) == 0) {
    rlang::warn("No markers found with current thresholds")
    return(list(markers = data.frame(), top_markers = character()))
  }
  
  # Adjust p-values
  markers_df$p_val_adj <- stats::p.adjust(markers_df$p_val, method = adjust.method)
  
  # Sort by cell type and log fold change
  markers_df <- markers_df[order(markers_df$cluster, -markers_df$avg_logfc), ]
  
  # Get top markers per cell type
  top_markers <- do.call(c, lapply(celltype_levels, function(ct) {
    ct_markers <- markers_df[markers_df$cluster == ct, ]
    head(ct_markers$gene, n.top)
  }))
  
  top_markers <- unique(top_markers)
  
  vmessage("  Found ", length(top_markers), " unique marker genes")
  
  list(
    markers = markers_df,
    top_markers = top_markers
  )
}

#' Find Highly Variable Genes
#'
#' @description Identify highly variable genes using the variance stabilization
#'   approach (similar to Seurat's VST method).
#'
#' @param matrix Count matrix (genes x cells)
#' @param n.top Number of top variable genes (default: 2000)
#' @param min.mean Minimum mean expression (default: 0.0125)
#' @param max.mean Maximum mean expression (default: 3)
#' @param min.var Minimum variance (default: 0.5)
#' @param span Loess span parameter (default: 0.3)
#'
#' @return A list containing:
#'   \itemize{
#'     \item{hvg}{Character vector of highly variable genes}
#'     \item{stats}{Data frame with gene statistics}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' hvg <- find_variable_genes(mat, n.top = 2000)
#' }
find_variable_genes <- function(matrix, n.top = 2000, min.mean = 0.0125,
                                 max.mean = 3, min.var = 0.5, span = 0.3) {
  
  vmessage("Finding highly variable genes...")
  
  # Normalize if needed
  if (max(matrix) > 100) {
    vmessage("  Normalizing data for HVG selection")
    matrix <- normalize_counts(matrix, log.transform = TRUE)
  }
  
  # Calculate mean and variance for each gene
  if (is_sparse(matrix)) {
    gene_mean <- Matrix::rowMeans(matrix)
    gene_var <- Matrix::rowMeans(matrix^2) - gene_mean^2
    gene_var <- gene_var * ncol(matrix) / (ncol(matrix) - 1)  # Sample variance
  } else {
    gene_mean <- rowMeans(matrix)
    gene_var <- apply(matrix, 1, stats::var)
  }
  
  # Create stats data frame
  stats_df <- data.frame(
    gene = rownames(matrix),
    mean = gene_mean,
    variance = gene_var,
    stringsAsFactors = FALSE
  )
  
  # Filter by mean expression
  stats_df <- stats_df[stats_df$mean >= min.mean & stats_df$mean <= max.mean, ]
  
  if (nrow(stats_df) == 0) {
    rlang::warn("No genes pass the mean expression filters")
    return(list(hvg = character(), stats = data.frame()))
  }
  
  # Fit loess to get expected variance
  loess_fit <- stats::loess(variance ~ mean, data = stats_df, span = span)
  stats_df$expected_var <- stats::predict(loess_fit)
  stats_df$expected_var[stats_df$expected_var < 0] <- 0
  
  # Calculate standardized variance
  stats_df$std_var <- (stats_df$variance - stats_df$expected_var) / 
                       sqrt(stats_df$expected_var + 1e-10)
  
  # Rank genes
  stats_df <- stats_df[order(-stats_df$std_var), ]
  stats_df$rank <- seq_len(nrow(stats_df))
  stats_df$variable <- stats_df$rank <= n.top
  
  hvg <- head(stats_df$gene, n.top)
  
  vmessage("  Selected ", length(hvg), " highly variable genes")
  
  list(
    hvg = hvg,
    stats = stats_df
  )
}

#' Get Shared Genes Between Datasets
#'
#' @description Find genes shared between single-cell and spatial datasets,
#'   optionally filtered by a provided gene list.
#'
#' @param sc_genes Character vector of single-cell genes
#' @param st_genes Character vector of spatial genes
#' @param gene_use Optional character vector of genes to use, or "All"
#'
#' @return Character vector of shared genes
#' @export
#'
#' @examples
#' \dontrun{
#' shared <- get_shared_genes(sc_genes, st_genes, marker_genes)
#' }
get_shared_genes <- function(sc_genes, st_genes, gene_use = NULL) {
  
  vmessage("Finding shared genes...")
  
  if (is.null(gene_use) || identical(gene_use, "All")) {
    shared <- intersect(sc_genes, st_genes)
  } else {
    shared <- intersect(intersect(sc_genes, st_genes), gene_use)
  }
  
  shared <- sort(shared)
  
  vmessage("  Found ", length(shared), " shared genes")
  
  shared
}

#' Subset Matrix by Genes
#'
#' @description Subset a matrix to include only specified genes.
#'
#' @param matrix Count matrix (genes x cells)
#' @param genes Character vector of genes to keep
#'
#' @return Subsetted matrix
#' @keywords internal
subset_by_genes <- function(matrix, genes) {
  genes_in_matrix <- intersect(genes, rownames(matrix))
  if (length(genes_in_matrix) == 0) {
    rlang::abort("No specified genes found in matrix")
  }
  matrix[genes_in_matrix, , drop = FALSE]
}
