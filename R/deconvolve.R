#' @title Cell Type Deconvolution Functions for STRIDER
#' @description Core functions for deconvolving cell type proportions in spatial transcriptomics
#' @name deconvolve
NULL

#' Run STRIDER Deconvolution Pipeline
#'
#' @description Main entry point for STRIDER cell type deconvolution. This function
#'   performs the complete pipeline including marker identification, LDA training,
#'   model selection, and spatial deconvolution.
#'
#' @param sc.data Single-cell data. Can be:
#'   \itemize{
#'     \item A matrix (genes x cells)
#'     \item A Seurat object (V4 or V5)
#'     \item A file path to count matrix (.h5 or .txt)
#'   }
#' @param st.data Spatial transcriptomics data. Same formats as sc.data.
#' @param sc.celltype Cell type annotations. Can be:
#'   \itemize{
#'     \item A named vector (cell -> celltype)
#'     \item Column name in Seurat metadata
#'     \item A file path to annotation file
#'   }
#' @param gene.use Genes to use for model training:
#'   \itemize{
#'     \item NULL: auto-detect markers (default)
#'     \item "All": use all shared genes
#'     \item Character vector: specific genes
#'     \item File path: gene list file
#'   }
#' @param n.topics Number of topics, or vector of numbers to test.
#'   Default: seq(n_celltypes, 3*n_celltypes)
#' @param method Deconvolution method: "auto", "Raw", "Norm", "NormBySD", "Bayes", "BayesNorm"
#' @param sc.scale.factor Scale factor for single-cell normalization (default: auto)
#' @param st.scale.factor Scale factor for spatial normalization (default: auto)
#' @param normalize Whether to scale by SD (default: FALSE)
#' @param n.markers Number of marker genes per cell type (default: 200)
#' @param lda.iter LDA iterations (default: 500)
#' @param seed Random seed (default: 42)
#' @param out.dir Output directory (default: ".")
#' @param out.prefix Output file prefix (default: "STRIDER")
#' @param save.model Save trained model (default: TRUE)
#'
#' @return A STRIDER object containing:
#'   \itemize{
#'     \item{deconv}{Cell type fractions matrix (spots x celltypes)}
#'     \item{topic_spot}{Topic distributions for spots (topics x spots)}
#'     \item{model}{Trained LDA model}
#'     \item{method}{Selected deconvolution method}
#'     \item{n_topics}{Selected number of topics}
#'     \item{metrics}{Model evaluation metrics}
#'   }
#'
#' @details
#' The STRIDER pipeline consists of the following steps:
#' \enumerate{
#'   \item Load and preprocess single-cell and spatial data
#'   \item Identify marker genes for each cell type
#'   \item Train LDA topic model on single-cell data
#'   \item Evaluate models with different topic numbers
#'   \item Select optimal model and deconvolution method
#'   \item Apply model to spatial data for deconvolution
#' }
#'
#' @references
#' Sun, D., et al. (2022). STRIDE: accurately decomposing and integrating
#' spatial transcriptomics using single-cell RNA sequencing.
#' Nucleic Acids Research.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # From matrices
#' result <- run_strider(
#'   sc.data = sc_matrix,
#'   st.data = st_matrix,
#'   sc.celltype = celltypes
#' )
#'
#' # From Seurat objects
#' result <- run_strider(
#'   sc.data = seurat_sc,
#'   st.data = seurat_st,
#'   sc.celltype = "celltype"
#' )
#' }
run_strider <- function(sc.data, st.data, sc.celltype,
                        gene.use = NULL, n.topics = NULL,
                        method = "auto",
                        sc.scale.factor = NULL, st.scale.factor = NULL,
                        normalize = FALSE, n.markers = 200,
                        lda.iter = 500, seed = 42,
                        out.dir = ".", out.prefix = "STRIDER",
                        save.model = TRUE) {
  
  vmessage("=== STRIDER Cell Type Deconvolution ===")
  vmessage("Author: Zaoqu Liu (liuzaoqu@163.com)")
  vmessage("")
  
  set_seed(seed)
  
  # Create output directory
  if (!dir.exists(out.dir)) {
    dir.create(out.dir, recursive = TRUE)
  }
  
  # ===== Step 1: Load single-cell data =====
  vmessage("Step 1: Loading single-cell data...")
  sc_data <- .load_sc_data(sc.data, sc.celltype)
  sc_matrix <- sc_data$matrix
  celltypes <- sc_data$celltypes
  
  vmessage("  ", nrow(sc_matrix), " genes x ", ncol(sc_matrix), " cells")
  vmessage("  ", length(unique(celltypes)), " cell types: ", 
           paste(head(unique(celltypes), 5), collapse = ", "),
           if(length(unique(celltypes)) > 5) "..." else "")
  
  # ===== Step 2: Load spatial data =====
  vmessage("\nStep 2: Loading spatial data...")
  st_data <- .load_st_data(st.data)
  st_matrix <- st_data$matrix
  
  vmessage("  ", nrow(st_matrix), " genes x ", ncol(st_matrix), " spots")
  
  # ===== Step 3: Preprocess single-cell data =====
  vmessage("\nStep 3: Preprocessing single-cell data...")
  sc_processed <- preprocess_sc(
    sc_matrix, celltypes,
    scale.factor = sc.scale.factor,
    normalize = TRUE, scale = normalize
  )
  
  # ===== Step 4: Preprocess spatial data =====
  vmessage("\nStep 4: Preprocessing spatial data...")
  st_processed <- preprocess_st(
    st_matrix,
    scale.factor = st.scale.factor,
    normalize = TRUE, scale = normalize
  )
  
  # ===== Step 5: Select genes =====
  vmessage("\nStep 5: Selecting genes for model training...")
  
  if (is.null(gene.use)) {
    # Find marker genes using differential expression
    marker_result <- find_markers(
      sc_processed$raw_matrix,
      sc_processed$celltypes,
      n.top = n.markers
    )
    genes_use <- marker_result$top_markers
    vmessage("  Identified ", length(genes_use), " marker genes")
  } else if (identical(gene.use, "All")) {
    genes_use <- "All"
    vmessage("  Using all shared genes")
  } else if (is.character(gene.use) && length(gene.use) == 1 && file.exists(gene.use)) {
    genes_use <- read_gene_list(gene.use)
    vmessage("  Loaded ", length(genes_use), " genes from file")
  } else {
    genes_use <- gene.use
    vmessage("  Using ", length(genes_use), " provided genes")
  }
  
  # Get shared genes
  shared_genes <- get_shared_genes(
    sc_processed$genes,
    st_processed$genes,
    genes_use
  )
  
  if (length(shared_genes) < 10) {
    rlang::abort("Too few shared genes. Check gene name format (e.g., symbol vs ENSEMBL).")
  }
  
  # Subset matrices to shared genes
  sc_matrix_use <- sc_processed$scale_matrix[shared_genes, ]
  st_matrix_use <- st_processed$scale_matrix[shared_genes, ]
  
  # ===== Step 6: Determine topic range =====
  if (is.null(n.topics)) {
    n_celltypes <- length(unique(sc_processed$celltypes))
    n.topics <- seq(n_celltypes, min(3 * n_celltypes, n_celltypes + 10))
    vmessage("\nAuto-selected topic range: ", min(n.topics), " to ", max(n.topics))
  }
  
  # ===== Step 7: Train and evaluate models =====
  vmessage("\nStep 6: Training and evaluating LDA models...")
  
  if (length(n.topics) > 1) {
    eval_results <- evaluate_topic_range(
      sc_matrix_use,
      sc_processed$celltypes,
      topic_range = n.topics,
      n.iter = lda.iter,
      seed = seed,
      save.models = save.model,
      out.dir = out.dir
    )
    
    best_model_info <- select_best_model(eval_results)
    model <- best_model_info$model
    topic_celltype <- best_model_info$topic_celltype
    celltype_topic_bayes <- best_model_info$celltype_topic_bayes
    selected_method <- if (method == "auto") best_model_info$method else method
    selected_n_topics <- best_model_info$n_topics
    metrics <- eval_results$metrics
    
  } else {
    # Train single model
    model <- train_lda(sc_matrix_use, n.topics = n.topics, n.iter = lda.iter, seed = seed)
    
    topic_celltype <- calculate_topic_celltype(model$topic_cell, sc_processed$celltypes)
    
    celltype_counts <- table(sc_processed$celltypes)
    celltype_topic_bayes <- calculate_celltype_topic_bayes(
      topic_celltype, celltype_counts, length(sc_processed$celltypes)
    )
    
    selected_method <- if (method == "auto") "Bayes" else method
    selected_n_topics <- n.topics
    metrics <- NULL
  }
  
  vmessage("\n  Selected: ", selected_n_topics, " topics, method = ", selected_method)
  
  # ===== Step 8: Deconvolve spatial data =====
  vmessage("\nStep 7: Deconvolving spatial transcriptomics...")
  
  deconv_result <- deconvolve(
    st_matrix_use,
    model,
    topic_celltype = topic_celltype,
    celltype_topic_bayes = celltype_topic_bayes,
    method = selected_method
  )
  
  # ===== Create output object =====
  result <- create_strider_object(
    deconv = deconv_result$celltype_frac,
    topic_spot = deconv_result$topic_spot,
    model = model,
    method = selected_method,
    n_topics = selected_n_topics,
    metrics = metrics,
    genes_used = shared_genes,
    celltypes = sort(unique(sc_processed$celltypes))
  )
  
  # ===== Step 9: Save results =====
  vmessage("\nStep 8: Saving results...")
  
  # Save deconvolution result
  deconv_file <- file.path(out.dir, paste0(out.prefix, "_spot_celltype_frac.txt"))
  write_results(deconv_result$celltype_frac, deconv_file)
  
  # Save topic distribution
  topic_file <- file.path(out.dir, paste0(out.prefix, "_topic_spot_mat_", selected_n_topics, ".txt"))
  write_results(deconv_result$topic_spot, topic_file)
  
  # Save model selection metrics
  if (!is.null(metrics)) {
    metrics_file <- file.path(out.dir, paste0(out.prefix, "_Model_selection.txt"))
    write_results(metrics, metrics_file)
  }
  
  # Save model and associated matrices
  if (save.model) {
    model_dir <- file.path(out.dir, "model")
    if (!dir.exists(model_dir)) dir.create(model_dir, recursive = TRUE)
    
    saveRDS(model, file.path(model_dir, paste0("lda_model_", selected_n_topics, ".rds")))
    
    # Save topic-celltype matrix
    tc_file <- file.path(model_dir, paste0("topic_celltype_mat_", selected_n_topics, ".txt"))
    write_results(topic_celltype, tc_file)
    
    # Save Bayesian celltype-topic matrix
    ct_bayes_file <- file.path(model_dir, paste0("celltype_topic_mat_bayes_", selected_n_topics, ".txt"))
    write_results(celltype_topic_bayes, ct_bayes_file)
    
    # Save gene dictionary
    gene_dict_file <- file.path(out.dir, "Gene_dict.txt")
    writeLines(model$vocab, gene_dict_file)
  }
  
  vmessage("\n=== STRIDER Complete ===")
  vmessage("Results saved to: ", out.dir)
  
  result
}

#' Deconvolve Spatial Spots
#'
#' @description Perform cell type deconvolution on spatial transcriptomics data
#'   using a trained LDA model.
#'
#' @param st.matrix Spatial count matrix (genes x spots)
#' @param model Trained LDA model from train_lda()
#' @param topic_celltype Topic-celltype association matrix
#' @param celltype_topic_bayes Bayesian celltype-topic matrix (required for Bayes methods)
#' @param method Deconvolution method:
#'   \itemize{
#'     \item "Raw": Direct matrix multiplication
#'     \item "Norm": Column-normalized weights
#'     \item "NormBySD": SD-normalized weights
#'     \item "Bayes": Bayesian posterior (default)
#'     \item "BayesNorm": Normalized Bayesian posterior
#'   }
#' @param n.iter Inference iterations (default: 100)
#'
#' @return A list containing:
#'   \itemize{
#'     \item{celltype_frac}{Cell type fractions (spots x celltypes)}
#'     \item{topic_spot}{Topic distributions (topics x spots)}
#'   }
#'
#' @details
#' The deconvolution computes cell type fractions for each spot by:
#' \enumerate{
#'   \item Inferring topic distributions for spatial spots
#'   \item Applying celltype-topic weights to convert topics to cell types
#'   \item Normalizing to obtain fractions that sum to 1
#' }
#'
#' @export
deconvolve <- function(st.matrix, model, topic_celltype,
                       celltype_topic_bayes = NULL,
                       method = "Bayes", n.iter = 100) {
  
  vmessage("Deconvolving with method: ", method)
  
  # Step 1: Transform spatial data to topic space
  topic_spot <- transform_lda(model, st.matrix, n.iter = n.iter)
  
  # Step 2: Get celltype-topic weights based on method
  celltype_topic <- t(topic_celltype)  # celltypes x topics
  
  ct_weights <- switch(
    method,
    "Raw" = celltype_topic,
    "Norm" = {
      col_sums <- colSums(celltype_topic)
      col_sums[col_sums == 0] <- 1
      sweep(celltype_topic, 2, col_sums, "/")
    },
    "NormBySD" = {
      # Use population SD (N), not sample SD (N-1), to match sklearn StandardScaler
      col_sds <- apply(celltype_topic, 2, function(x) {
        sqrt(mean((x - mean(x))^2))
      })
      col_sds[col_sds == 0] <- 1
      sweep(celltype_topic, 2, col_sds, "/")
    },
    "Bayes" = {
      if (is.null(celltype_topic_bayes)) {
        rlang::abort("celltype_topic_bayes required for Bayes method")
      }
      celltype_topic_bayes
    },
    "BayesNorm" = {
      if (is.null(celltype_topic_bayes)) {
        rlang::abort("celltype_topic_bayes required for BayesNorm method")
      }
      col_sums <- colSums(celltype_topic_bayes)
      col_sums[col_sums == 0] <- 1
      sweep(celltype_topic_bayes, 2, col_sums, "/")
    },
    rlang::abort(paste0("Unknown method: ", method))
  )
  
  # Ensure dimension compatibility
  n_topics <- nrow(topic_spot)
  if (ncol(ct_weights) > n_topics) {
    ct_weights <- ct_weights[, seq_len(n_topics), drop = FALSE]
  }
  
  # Step 3: Deconvolve - matrix multiplication
  # (celltypes x topics) %*% (topics x spots) = (celltypes x spots)
  celltype_spot <- ct_weights %*% topic_spot
  
  # Step 4: Normalize to fractions (each column sums to 1)
  col_sums <- colSums(celltype_spot)
  col_sums[col_sums == 0] <- 1
  celltype_frac <- sweep(celltype_spot, 2, col_sums, "/")
  
  # Transpose to spots x celltypes (conventional format)
  celltype_frac <- t(celltype_frac)
  
  vmessage("  Deconvolved ", ncol(topic_spot), " spots into ", 
           nrow(ct_weights), " cell types")
  
  list(
    celltype_frac = celltype_frac,
    topic_spot = topic_spot
  )
}

#' Train LDA Model (Convenience Wrapper)
#'
#' @description Train an LDA model from single-cell data with automatic
#'   topic-celltype association calculation.
#'
#' @param sc.matrix Single-cell count matrix (genes x cells)
#' @param celltypes Named vector of cell type annotations
#' @param gene.use Genes to use (optional)
#' @param n.topics Number of topics
#' @param n.iter LDA iterations
#' @param seed Random seed
#'
#' @return Trained LDA model object with topic-celltype matrices
#' @export
train_model <- function(sc.matrix, celltypes, gene.use = NULL,
                        n.topics = 10, n.iter = 500, seed = 42) {
  
  # Subset genes if specified
  if (!is.null(gene.use) && !identical(gene.use, "All")) {
    shared <- intersect(gene.use, rownames(sc.matrix))
    if (length(shared) == 0) {
      rlang::abort("No specified genes found in matrix")
    }
    sc.matrix <- sc.matrix[shared, , drop = FALSE]
  }
  
  # Train model
  model <- train_lda(sc.matrix, n.topics = n.topics, n.iter = n.iter, seed = seed)
  
  # Calculate topic-celltype associations
  topic_celltype <- calculate_topic_celltype(model$topic_cell, celltypes)
  
  # Calculate Bayesian associations
  celltype_counts <- table(celltypes)
  celltype_topic_bayes <- calculate_celltype_topic_bayes(
    topic_celltype, celltype_counts, length(celltypes)
  )
  
  model$topic_celltype <- topic_celltype
  model$celltype_topic_bayes <- celltype_topic_bayes
  
  model
}

#' Create STRIDER Object
#'
#' @description Create a STRIDER result object.
#'
#' @param deconv Deconvolution results matrix (spots x celltypes)
#' @param topic_spot Topic-spot distribution matrix
#' @param model Trained LDA model
#' @param method Selected deconvolution method
#' @param n_topics Number of topics
#' @param metrics Model evaluation metrics
#' @param genes_used Genes used in model training
#' @param celltypes Cell types in the model
#'
#' @return A STRIDER object
#' @export
create_strider_object <- function(deconv, topic_spot = NULL, model = NULL,
                                   method = NULL, n_topics = NULL,
                                   metrics = NULL, genes_used = NULL,
                                   celltypes = NULL) {
  
  obj <- list(
    deconv = deconv,
    topic_spot = topic_spot,
    model = model,
    method = method,
    n_topics = n_topics,
    metrics = metrics,
    genes_used = genes_used,
    celltypes = celltypes
  )
  
  class(obj) <- "STRIDER"
  obj
}

#' Print STRIDER Object
#' @param x A STRIDER object
#' @param ... Additional arguments (unused)
#' @export
print.STRIDER <- function(x, ...) {
  cat("STRIDER Deconvolution Result\n")
  cat("============================\n")
  cat("Spots:", nrow(x$deconv), "\n")
  cat("Cell types:", ncol(x$deconv), "\n")
  if (!is.null(x$n_topics)) cat("Topics:", x$n_topics, "\n")
  if (!is.null(x$method)) cat("Method:", x$method, "\n")
  if (!is.null(x$genes_used)) cat("Genes:", length(x$genes_used), "\n")
  invisible(x)
}

#' Summary of STRIDER Object
#' @param object A STRIDER object
#' @param ... Additional arguments (unused)
#' @export
summary.STRIDER <- function(object, ...) {
  cat("STRIDER Deconvolution Summary\n")
  cat("=============================\n\n")
  
  cat("Data:\n")
  cat("  Spots:", nrow(object$deconv), "\n")
  cat("  Cell types:", paste(colnames(object$deconv), collapse = ", "), "\n\n")
  
  cat("Model:\n")
  cat("  Topics:", object$n_topics, "\n")
  cat("  Method:", object$method, "\n")
  cat("  Genes:", length(object$genes_used), "\n\n")
  
  cat("Cell type proportions:\n")
  ct_means <- colMeans(object$deconv)
  for (ct in names(sort(ct_means, decreasing = TRUE))) {
    cat("  ", ct, ": ", sprintf("%.1f%%", ct_means[ct] * 100), "\n", sep = "")
  }
  
  invisible(object)
}

# ============================================================================
# Internal Helper Functions
# ============================================================================

#' @keywords internal
.load_sc_data <- function(sc.data, sc.celltype) {
  # Handle matrix input
  if (is.matrix(sc.data) || is_sparse(sc.data)) {
    sc_matrix <- as_sparse(sc.data)
  }
  # Handle file path input
  else if (is.character(sc.data) && length(sc.data) == 1 && file.exists(sc.data)) {
    if (grepl("\\.h5$", sc.data, ignore.case = TRUE)) {
      data <- read_10x_h5(sc.data)
      sc_matrix <- data$matrix
    } else {
      data <- read_count_matrix(sc.data)
      sc_matrix <- data$matrix
    }
  }
  # Handle Seurat object input
  else if (inherits(sc.data, "Seurat")) {
    sc_matrix <- .extract_seurat_counts(sc.data)
  }
  else {
    rlang::abort("sc.data must be a matrix, file path, or Seurat object")
  }
  
  # Handle cell type annotations
  if (is.character(sc.celltype) && length(sc.celltype) == 1) {
    if (file.exists(sc.celltype)) {
      celltypes <- read_celltype_annotations(sc.celltype)
    } else if (inherits(sc.data, "Seurat")) {
      celltypes <- .extract_seurat_metadata(sc.data, sc.celltype)
    } else {
      rlang::abort(paste0("Cell type file not found: ", sc.celltype))
    }
  } else if (is.vector(sc.celltype) && !is.null(names(sc.celltype))) {
    celltypes <- sc.celltype
  } else if (is.factor(sc.celltype) || is.character(sc.celltype)) {
    if (length(sc.celltype) == ncol(sc_matrix)) {
      celltypes <- sc.celltype
      names(celltypes) <- colnames(sc_matrix)
    } else {
      rlang::abort("Length of sc.celltype must match number of cells")
    }
  } else {
    rlang::abort("sc.celltype must be a named vector, file path, or Seurat column name")
  }
  
  list(matrix = sc_matrix, celltypes = celltypes)
}

#' @keywords internal
.load_st_data <- function(st.data) {
  if (is.matrix(st.data) || is_sparse(st.data)) {
    st_matrix <- as_sparse(st.data)
  } else if (is.character(st.data) && length(st.data) == 1 && file.exists(st.data)) {
    if (grepl("\\.h5$", st.data, ignore.case = TRUE)) {
      data <- read_10x_h5(st.data)
      st_matrix <- data$matrix
    } else {
      data <- read_count_matrix(st.data)
      st_matrix <- data$matrix
    }
  } else if (inherits(st.data, "Seurat")) {
    st_matrix <- .extract_seurat_counts(st.data)
  } else {
    rlang::abort("st.data must be a matrix, file path, or Seurat object")
  }
  
  list(matrix = st_matrix)
}

#' @keywords internal
.extract_seurat_counts <- function(seurat_obj) {
  # Support both Seurat v4 and v5
  counts <- NULL
  
  # Try SeuratObject functions first (v5 compatible)
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    tryCatch({
      # Try GetAssayData with layer argument (v5)
      counts <- SeuratObject::GetAssayData(seurat_obj, layer = "counts")
    }, error = function(e) {
      tryCatch({
        # Try slot argument (v4)
        counts <<- SeuratObject::GetAssayData(seurat_obj, slot = "counts")
      }, error = function(e2) {
        NULL
      })
    })
  }
  
  # Fallback to direct slot access (v4)
  if (is.null(counts)) {
    tryCatch({
      counts <- seurat_obj@assays$RNA@counts
    }, error = function(e) {
      # Try Spatial assay
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
    rlang::abort(paste0(
      "Column '", column, "' not found in Seurat metadata. ",
      "Available columns: ", paste(colnames(seurat_obj@meta.data), collapse = ", ")
    ))
  }
  
  celltypes <- seurat_obj@meta.data[[column]]
  names(celltypes) <- rownames(seurat_obj@meta.data)
  
  celltypes
}
