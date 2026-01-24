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
#'   - A matrix (genes x cells)
#'   - A Seurat object (V4 or V5)
#'   - A file path to count matrix
#' @param st.data Spatial transcriptomics data. Can be:
#'   - A matrix (genes x spots)
#'   - A Seurat object
#'   - A file path to count matrix
#' @param sc.celltype Cell type annotations. Can be:
#'   - A named vector (cell -> celltype)
#'   - Column name in Seurat metadata
#'   - A file path to annotation file
#' @param gene.use Genes to use for model training:
#'   - NULL: auto-detect markers (default)
#'   - "All": use all shared genes
#'   - Character vector: specific genes
#'   - File path: gene list file
#' @param n.topics Number of topics, or vector of numbers to test
#' @param method Deconvolution method: "Raw", "Norm", "NormBySD", "Bayes", "BayesNorm"
#' @param sc.scale.factor Scale factor for single-cell normalization
#' @param st.scale.factor Scale factor for spatial normalization
#' @param normalize Whether to scale by SD (default: FALSE)
#' @param n.markers Number of markers per cell type (default: 200)
#' @param lda.iter LDA iterations (default: 500)
#' @param seed Random seed (default: 42)
#' @param out.dir Output directory (default: current directory)
#' @param out.prefix Output file prefix (default: "STRIDER")
#' @param save.model Save trained model (default: TRUE)
#'
#' @return A STRIDER object containing:
#'   \itemize{
#'     \item{deconv}{Cell type fractions matrix (spots x celltypes)}
#'     \item{topic_spot}{Topic distributions for spots}
#'     \item{model}{Trained LDA model}
#'     \item{method}{Selected method}
#'     \item{n_topics}{Selected number of topics}
#'     \item{metrics}{Model evaluation metrics}
#'   }
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
#'
#' # From files
#' result <- run_strider(
#'   sc.data = "sc_count.h5",
#'   st.data = "st_count.txt",
#'   sc.celltype = "celltype.txt"
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
  
  # ===== Load single-cell data =====
  vmessage("Step 1: Loading single-cell data...")
  sc_data <- .load_sc_data(sc.data, sc.celltype)
  sc_matrix <- sc_data$matrix
  celltypes <- sc_data$celltypes
  
  vmessage("  ", nrow(sc_matrix), " genes x ", ncol(sc_matrix), " cells")
  vmessage("  ", length(unique(celltypes)), " cell types")
  
  # ===== Load spatial data =====
  vmessage("\nStep 2: Loading spatial data...")
  st_data <- .load_st_data(st.data)
  st_matrix <- st_data$matrix
  
  vmessage("  ", nrow(st_matrix), " genes x ", ncol(st_matrix), " spots")
  
  # ===== Preprocess single-cell data =====
  vmessage("\nStep 3: Preprocessing single-cell data...")
  sc_processed <- preprocess_sc(
    sc_matrix, celltypes,
    scale.factor = sc.scale.factor,
    normalize = TRUE, scale = normalize
  )
  
  # ===== Preprocess spatial data =====
  vmessage("\nStep 4: Preprocessing spatial data...")
  st_processed <- preprocess_st(
    st_matrix,
    scale.factor = st.scale.factor,
    normalize = TRUE, scale = normalize
  )
  
  # ===== Get genes to use =====
  vmessage("\nStep 5: Selecting genes for model training...")
  
  if (is.null(gene.use)) {
    # Find marker genes
    marker_result <- find_markers(
      sc_processed$raw_matrix,
      sc_processed$celltypes,
      n.top = n.markers
    )
    genes_use <- marker_result$top_markers
  } else if (identical(gene.use, "All")) {
    genes_use <- "All"
  } else if (is.character(gene.use) && length(gene.use) == 1 && file.exists(gene.use)) {
    genes_use <- read_gene_list(gene.use)
  } else {
    genes_use <- gene.use
  }
  
  # Get shared genes
  shared_genes <- get_shared_genes(
    sc_processed$genes,
    st_processed$genes,
    genes_use
  )
  
  if (length(shared_genes) < 10) {
    rlang::abort("Too few shared genes. Check gene name format.")
  }
  
  # Subset matrices
  sc_matrix_use <- sc_processed$scale_matrix[shared_genes, ]
  st_matrix_use <- st_processed$scale_matrix[shared_genes, ]
  
  vmessage("  Using ", length(shared_genes), " genes")
  
  # ===== Determine topic range =====
  if (is.null(n.topics)) {
    n_celltypes <- length(unique(celltypes))
    n.topics <- seq(n_celltypes, 3 * n_celltypes)
  }
  
  # ===== Train and evaluate models =====
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
  
  # ===== Deconvolve spatial data =====
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
    celltypes = unique(sc_processed$celltypes)
  )
  
  # ===== Save results =====
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
  
  # Save model
  if (save.model) {
    model_dir <- file.path(out.dir, "model")
    if (!dir.exists(model_dir)) dir.create(model_dir, recursive = TRUE)
    
    saveRDS(model, file.path(model_dir, paste0("lda_model_", selected_n_topics, ".rds")))
    
    # Save topic-celltype matrix
    tc_file <- file.path(model_dir, paste0("topic_celltype_mat_", selected_n_topics, ".txt"))
    write_results(topic_celltype, tc_file)
    
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
#' @param celltype_topic_bayes Bayesian celltype-topic matrix (optional)
#' @param method Deconvolution method
#' @param n.iter Inference iterations (default: 100)
#'
#' @return A list containing:
#'   \itemize{
#'     \item{celltype_frac}{Cell type fractions (spots x celltypes)}
#'     \item{topic_spot}{Topic distributions (topics x spots)}
#'   }
#'
#' @export
deconvolve <- function(st.matrix, model, topic_celltype,
                       celltype_topic_bayes = NULL,
                       method = "Bayes", n.iter = 100) {
  
  vmessage("Deconvolving with method: ", method)
  
  # Transform spatial data to topic space
  topic_spot <- transform_lda(model, st.matrix, n.iter = n.iter)
  
  # Get celltype-topic weights based on method
  celltype_topic <- t(topic_celltype)
  
  if (method == "Raw") {
    ct_weights <- celltype_topic
  } else if (method == "Norm") {
    col_sums <- colSums(celltype_topic)
    col_sums[col_sums == 0] <- 1
    ct_weights <- sweep(celltype_topic, 2, col_sums, "/")
  } else if (method == "NormBySD") {
    col_sds <- apply(celltype_topic, 2, stats::sd)
    col_sds[col_sds == 0] <- 1
    ct_weights <- sweep(celltype_topic, 2, col_sds, "/")
  } else if (method == "Bayes") {
    if (is.null(celltype_topic_bayes)) {
      rlang::abort("celltype_topic_bayes required for Bayes method")
    }
    ct_weights <- celltype_topic_bayes
  } else if (method == "BayesNorm") {
    if (is.null(celltype_topic_bayes)) {
      rlang::abort("celltype_topic_bayes required for BayesNorm method")
    }
    col_sums <- colSums(celltype_topic_bayes)
    col_sums[col_sums == 0] <- 1
    ct_weights <- sweep(celltype_topic_bayes, 2, col_sums, "/")
  }
  
  # Ensure dimensions match
  n_topics <- nrow(topic_spot)
  if (ncol(ct_weights) > n_topics) {
    ct_weights <- ct_weights[, seq_len(n_topics), drop = FALSE]
  }
  
  # Deconvolve: (celltypes x topics) %*% (topics x spots) = (celltypes x spots)
  celltype_spot <- ct_weights %*% topic_spot
  
  # Normalize to fractions
  col_sums <- colSums(celltype_spot)
  col_sums[col_sums == 0] <- 1
  celltype_frac <- sweep(celltype_spot, 2, col_sums, "/")
  
  # Transpose to spots x celltypes
  celltype_frac <- t(celltype_frac)
  
  vmessage("  Deconvolved ", ncol(topic_spot), " spots into ", nrow(ct_weights), " cell types")
  
  list(
    celltype_frac = celltype_frac,
    topic_spot = topic_spot
  )
}

#' Train LDA Model
#'
#' @description Train an LDA model from single-cell data. This is a wrapper
#'   around the core train_lda function with additional preprocessing.
#'
#' @param sc.matrix Single-cell count matrix (genes x cells)
#' @param celltypes Named vector of cell type annotations
#' @param gene.use Genes to use (optional)
#' @param n.topics Number of topics
#' @param n.iter LDA iterations
#' @param seed Random seed
#'
#' @return Trained LDA model object
#' @export
train_model <- function(sc.matrix, celltypes, gene.use = NULL,
                        n.topics = 10, n.iter = 500, seed = 42) {
  
  # Subset genes if specified
  if (!is.null(gene.use) && !identical(gene.use, "All")) {
    shared <- intersect(gene.use, rownames(sc.matrix))
    sc.matrix <- sc.matrix[shared, ]
  }
  
  # Train model
  model <- train_lda(sc.matrix, n.topics = n.topics, n.iter = n.iter, seed = seed)
  
  # Calculate associations
  topic_celltype <- calculate_topic_celltype(model$topic_cell, celltypes)
  
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
#' @param deconv Deconvolution results matrix
#' @param topic_spot Topic-spot matrix
#' @param model LDA model
#' @param method Selected method
#' @param n_topics Number of topics
#' @param metrics Evaluation metrics
#' @param genes_used Genes used in model
#' @param celltypes Cell types in model
#'
#' @return STRIDER object
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
#' @export
print.STRIDER <- function(x, ...) {
  cat("STRIDER Deconvolution Result\n")
  cat("============================\n")
  cat("Spots:", nrow(x$deconv), "\n")
  cat("Cell types:", ncol(x$deconv), "\n")
  if (!is.null(x$n_topics)) cat("Topics:", x$n_topics, "\n")
  if (!is.null(x$method)) cat("Method:", x$method, "\n")
  if (!is.null(x$genes_used)) cat("Genes used:", length(x$genes_used), "\n")
  invisible(x)
}

#' Summary of STRIDER Object
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
  
  cat("Cell type proportions summary:\n")
  print(summary(object$deconv))
  
  invisible(object)
}

# Internal helper functions

#' @keywords internal
.load_sc_data <- function(sc.data, sc.celltype) {
  # Handle different input types
  
  if (is.character(sc.data) && length(sc.data) == 1 && file.exists(sc.data)) {
    # File path
    if (grepl("\\.h5$", sc.data)) {
      data <- read_10x_h5(sc.data)
      sc_matrix <- data$matrix
    } else {
      data <- read_count_matrix(sc.data)
      sc_matrix <- data$matrix
    }
  } else if (is.matrix(sc.data) || is_sparse(sc.data)) {
    sc_matrix <- as_sparse(sc.data)
  } else if (inherits(sc.data, "Seurat")) {
    sc_matrix <- .extract_seurat_counts(sc.data)
  } else {
    rlang::abort("sc.data must be a matrix, file path, or Seurat object")
  }
  
  # Handle cell types
  if (is.character(sc.celltype) && length(sc.celltype) == 1) {
    if (file.exists(sc.celltype)) {
      celltypes <- read_celltype_annotations(sc.celltype)
    } else if (inherits(sc.data, "Seurat")) {
      # Column name in Seurat metadata
      celltypes <- .extract_seurat_metadata(sc.data, sc.celltype)
    } else {
      rlang::abort("sc.celltype file not found")
    }
  } else if (is.vector(sc.celltype) && !is.null(names(sc.celltype))) {
    celltypes <- sc.celltype
  } else {
    rlang::abort("sc.celltype must be a named vector or file path")
  }
  
  list(matrix = sc_matrix, celltypes = celltypes)
}

#' @keywords internal
.load_st_data <- function(st.data) {
  if (is.character(st.data) && length(st.data) == 1 && file.exists(st.data)) {
    if (grepl("\\.h5$", st.data)) {
      data <- read_10x_h5(st.data)
      st_matrix <- data$matrix
    } else {
      data <- read_count_matrix(st.data)
      st_matrix <- data$matrix
    }
  } else if (is.matrix(st.data) || is_sparse(st.data)) {
    st_matrix <- as_sparse(st.data)
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
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    # Try v5 style first
    tryCatch({
      counts <- SeuratObject::GetAssayData(seurat_obj, slot = "counts")
    }, error = function(e) {
      # Fall back to v4 style
      counts <<- seurat_obj@assays$RNA@counts
    })
  } else {
    # Direct access (v4)
    counts <- seurat_obj@assays$RNA@counts
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
