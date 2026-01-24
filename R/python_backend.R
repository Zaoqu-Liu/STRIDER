#' @title Python Backend for STRIDER
#' @description Use Python (gensim, POT, scanpy) via reticulate for exact equivalence
#'   with the original STRIDE implementation.
#' @name python_backend
NULL

# ============================================================================
# Python Environment Setup
# ============================================================================

#' Check if Python Backend is Available
#'
#' @description Check if required Python packages are available.
#'
#' @param packages Character vector of package names to check
#' @return Logical indicating if all packages are available
#' @export
check_python_backend <- function(packages = c("gensim", "numpy", "scipy", "sklearn")) {
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(FALSE)
  }
  
  tryCatch({
    for (pkg in packages) {
      if (!reticulate::py_module_available(pkg)) {
        return(FALSE)
      }
    }
    TRUE
  }, error = function(e) {
    FALSE
  })
}

#' Setup Python Environment for STRIDER
#'
#' @description Initialize Python environment with required packages.
#'
#' @param conda_env Name of conda environment (optional)
#' @param virtualenv Path to virtualenv (optional)
#' @param required_packages Packages to install if missing
#' @param verbose Print messages
#'
#' @return Invisible TRUE on success
#' @export
setup_python_backend <- function(conda_env = NULL, 
                                  virtualenv = NULL,
                                  required_packages = c("gensim", "numpy", "scipy", 
                                                        "scikit-learn", "POT", "scanpy"),
                                  verbose = TRUE) {
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    rlang::abort("Package 'reticulate' is required. Install with: install.packages('reticulate')")
  }
  
  # Configure Python
  if (!is.null(conda_env)) {
    reticulate::use_condaenv(conda_env, required = TRUE)
  } else if (!is.null(virtualenv)) {
    reticulate::use_virtualenv(virtualenv, required = TRUE)
  }
  
  # Check and install packages
  for (pkg in required_packages) {
    if (!reticulate::py_module_available(pkg)) {
      if (verbose) vmessage("Installing Python package: ", pkg)
      reticulate::py_install(pkg)
    }
  }
  
  if (verbose) vmessage("Python backend ready")
  invisible(TRUE)
}

# ============================================================================
# LDA Training via Gensim (Exact Match with Python STRIDE)
# ============================================================================

#' Train LDA Model Using Gensim (Python)
#'
#' @description Train LDA model using gensim for exact equivalence with Python STRIDE.
#'
#' @param matrix Normalized count matrix (genes x cells)
#' @param n.topics Number of topics
#' @param n.iter Number of iterations (passes)
#' @param alpha Dirichlet prior (default: 'symmetric')
#' @param eta Topic-word prior (default: 'auto')
#' @param seed Random seed
#' @param verbose Print progress
#'
#' @return LDA model object compatible with STRIDER
#' @export
train_lda_gensim <- function(matrix, n.topics = 10, n.iter = 50,
                              alpha = "symmetric", eta = "auto",
                              seed = 42, verbose = TRUE) {
  
  if (!check_python_backend(c("gensim", "numpy"))) {
    rlang::abort("Python backend not available. Run setup_python_backend() first.")
  }
  
  if (verbose) vmessage("Training LDA using gensim (Python backend)...")
  
  # Import Python modules
  gensim <- reticulate::import("gensim")
  np <- reticulate::import("numpy")
  scipy_sparse <- reticulate::import("scipy.sparse")
  
  # Prepare data - transpose so cells are rows
  if (is_sparse(matrix)) {
    mat_t <- Matrix::t(matrix)
    # Convert to scipy sparse
    mat_scipy <- scipy_sparse$csc_matrix(
      reticulate::tuple(mat_t@x, mat_t@i, mat_t@p),
      shape = reticulate::tuple(as.integer(nrow(mat_t)), as.integer(ncol(mat_t)))
    )
  } else {
    mat_t <- t(matrix)
    mat_scipy <- scipy_sparse$csc_matrix(np$array(mat_t))
  }
  
  # Store vocabulary (genes)
  vocab <- colnames(mat_t)
  if (is.null(vocab)) {
    vocab <- paste0("Gene_", seq_len(ncol(mat_t)))
  }
  cell_ids <- rownames(mat_t)
  if (is.null(cell_ids)) {
    cell_ids <- paste0("Cell_", seq_len(nrow(mat_t)))
  }
  
  # Create gensim dictionary
  genes_dict <- gensim$corpora$Dictionary(list(vocab))
  
  # Convert to gensim corpus format using Sparse2Corpus
  corpus <- gensim$matutils$Sparse2Corpus(mat_scipy)
  
  # Train LDA model - exactly matching STRIDE
  if (verbose) vmessage("  Training with ", n.topics, " topics, ", n.iter, " passes...")
  
  lda <- gensim$models$LdaModel(
    corpus = corpus,
    num_topics = as.integer(n.topics),
    id2word = genes_dict,
    passes = as.integer(n.iter),
    alpha = alpha,
    eta = eta,
    random_state = as.integer(seed)
  )
  
  # Extract topic-cell matrix
  topic_cell_corpus <- lda$get_document_topics(corpus)
  topic_cell_mat <- gensim$matutils$corpus2csc(topic_cell_corpus)
  topic_cell <- as.matrix(topic_cell_mat$todense())
  rownames(topic_cell) <- paste0("Topic_", seq_len(n.topics))
  colnames(topic_cell) <- cell_ids
  
  # Extract topic-gene matrix
  topic_gene_list <- lda$get_topics()
  topic_gene <- np$array(topic_gene_list)
  topic_gene <- as.matrix(topic_gene)
  rownames(topic_gene) <- paste0("Topic_", seq_len(n.topics))
  colnames(topic_gene) <- vocab
  
  if (verbose) vmessage("  LDA training complete (gensim)")
  
  structure(
    list(
      topic_gene = topic_gene,
      topic_cell = topic_cell,
      vocab = vocab,
      n_topics = n.topics,
      n_iter = n.iter,
      alpha = alpha,
      eta = eta,
      gensim_model = lda,
      genes_dict = genes_dict,
      backend = "gensim"
    ),
    class = c("strider_lda_gensim", "strider_lda")
  )
}

#' Transform Data Using Gensim LDA Model
#'
#' @description Apply gensim LDA model to new data.
#'
#' @param model Trained gensim LDA model
#' @param newdata New data matrix (genes x spots)
#' @param verbose Print progress
#'
#' @return Topic distribution matrix (topics x spots)
#' @export
transform_lda_gensim <- function(model, newdata, verbose = TRUE) {
  
  if (!inherits(model, "strider_lda_gensim")) {
    rlang::abort("model must be a strider_lda_gensim object")
  }
  
  if (verbose) vmessage("Transforming data using gensim...")
  
  gensim <- reticulate::import("gensim")
  np <- reticulate::import("numpy")
  scipy_sparse <- reticulate::import("scipy.sparse")
  
  # Get shared genes
  shared_genes <- intersect(model$vocab, rownames(newdata))
  gene_idx <- match(shared_genes, model$vocab)
  
  if (verbose) vmessage("  Using ", length(shared_genes), " shared genes")
  
  # Subset and reorder
  newdata_subset <- newdata[shared_genes, , drop = FALSE]
  
  # Create full matrix with correct vocabulary order
  full_mat <- matrix(0, nrow = ncol(newdata_subset), ncol = length(model$vocab))
  colnames(full_mat) <- model$vocab
  rownames(full_mat) <- colnames(newdata_subset)
  full_mat[, gene_idx] <- t(as.matrix(newdata_subset))
  
  # Convert to corpus
  spot_ids <- rownames(full_mat)
  mat_scipy <- scipy_sparse$csc_matrix(np$array(full_mat))
  corpus <- gensim$matutils$Sparse2Corpus(mat_scipy)
  
  # Transform
  topic_spot_corpus <- model$gensim_model$get_document_topics(corpus)
  topic_spot_mat <- gensim$matutils$corpus2csc(topic_spot_corpus)
  topic_spot <- as.matrix(topic_spot_mat$todense())
  
  rownames(topic_spot) <- paste0("Topic_", seq_len(model$n_topics))
  colnames(topic_spot) <- spot_ids
  
  if (verbose) vmessage("  Transformation complete")
  
  topic_spot
}

# ============================================================================
# Coherence Calculation (UMass and C_v) via Gensim
# ============================================================================

#' Calculate Coherence Using Gensim
#'
#' @description Calculate both UMass and C_v coherence using gensim.
#'
#' @param model LDA model (gensim backend)
#' @param matrix Original count matrix
#' @param cell_gene_list List of genes expressed per cell (for c_v)
#' @param coherence_type Type: "u_mass", "c_v", or "both"
#'
#' @return List with coherence scores
#' @export
calculate_coherence_gensim <- function(model, matrix = NULL, cell_gene_list = NULL,
                                        coherence_type = "both") {
  
  if (!inherits(model, "strider_lda_gensim")) {
    rlang::abort("model must be a strider_lda_gensim object")
  }
  
  vmessage("Calculating coherence using gensim...")
  
  gensim <- reticulate::import("gensim")
  
  # Get corpus from model
  lda <- model$gensim_model
  genes_dict <- model$genes_dict
  
  # We need the corpus - reconstruct if needed
  if (is.null(matrix)) {
    rlang::abort("matrix is required for coherence calculation")
  }
  
  np <- reticulate::import("numpy")
  scipy_sparse <- reticulate::import("scipy.sparse")
  
  # Reconstruct corpus
  mat_t <- t(as.matrix(matrix[model$vocab, , drop = FALSE]))
  mat_scipy <- scipy_sparse$csc_matrix(np$array(mat_t))
  corpus <- gensim$matutils$Sparse2Corpus(mat_scipy)
  
  result <- list()
  
  # UMass coherence
  if (coherence_type %in% c("u_mass", "both")) {
    vmessage("  Calculating UMass coherence...")
    cm_umass <- gensim$models$CoherenceModel(
      model = lda,
      corpus = corpus,
      coherence = "u_mass"
    )
    result$umass <- cm_umass$get_coherence()
    vmessage("    UMass: ", round(result$umass, 4))
  }
  
  # C_v coherence
  if (coherence_type %in% c("c_v", "both")) {
    vmessage("  Calculating C_v coherence...")
    
    # Create cell_gene_list if not provided
    if (is.null(cell_gene_list)) {
      vmessage("    Building gene lists per cell...")
      cell_gene_list <- lapply(seq_len(ncol(matrix)), function(i) {
        expressed <- matrix[, i] > 0
        model$vocab[expressed[model$vocab]]
      })
    }
    
    # Convert to Python list
    texts <- reticulate::r_to_py(cell_gene_list)
    
    cm_cv <- gensim$models$CoherenceModel(
      model = lda,
      corpus = corpus,
      texts = texts,
      coherence = "c_v"
    )
    result$cv <- cm_cv$get_coherence()
    vmessage("    C_v: ", round(result$cv, 4))
  }
  
  result
}

# ============================================================================
# Fused Gromov-Wasserstein via POT Library
# ============================================================================

#' Compute FGW Distance Using POT (Python Optimal Transport)
#'
#' @description Compute Fused Gromov-Wasserstein distance using the POT library
#'   for exact equivalence with Python STRIDE.
#'
#' @param M Feature cost matrix (n x m)
#' @param C1 Source structure matrix (n x n)
#' @param C2 Target structure matrix (m x m)
#' @param p Source distribution (n)
#' @param q Target distribution (m)
#' @param alpha Trade-off parameter (0 = OT, 1 = GW)
#' @param loss_fun Loss function: "square_loss" or "kl_loss"
#' @param verbose Print progress
#'
#' @return List with transport plan T and FGW distance
#' @export
fused_gromov_wasserstein_pot <- function(M, C1, C2, p, q, 
                                          alpha = 0.5, 
                                          loss_fun = "square_loss",
                                          verbose = TRUE) {
  
  if (!check_python_backend("ot")) {
    rlang::abort("POT library not available. Run setup_python_backend() first.")
  }
  
  if (verbose) vmessage("Computing FGW using POT (Python backend)...")
  
  ot <- reticulate::import("ot")
  np <- reticulate::import("numpy")
  
  # Convert to numpy arrays
  M_np <- np$array(as.matrix(M), dtype = "float64")
  C1_np <- np$array(as.matrix(C1), dtype = "float64")
  C2_np <- np$array(as.matrix(C2), dtype = "float64")
  p_np <- np$array(as.numeric(p), dtype = "float64")
  q_np <- np$array(as.numeric(q), dtype = "float64")
  
  # Call POT's fused_gromov_wasserstein - exactly like STRIDE
  result <- ot$gromov$fused_gromov_wasserstein(
    M = M_np,
    C1 = C1_np,
    C2 = C2_np,
    p = p_np,
    q = q_np,
    loss_fun = loss_fun,
    alpha = alpha
  )
  
  # Convert result back to R
  T_mat <- as.matrix(result)
  
  # Compute FGW distance
  # FGW = (1-alpha) * <M, T> + alpha * GW(C1, C2, T)
  ot_cost <- sum(M * T_mat)
  
  # GW cost: sum_{ijkl} L(C1_ik, C2_jl) * T_ij * T_kl
  # For square loss: sum_{ijkl} (C1_ik - C2_jl)^2 * T_ij * T_kl
  gw_cost <- 0
  for (i in seq_len(nrow(T_mat))) {
    for (j in seq_len(ncol(T_mat))) {
      for (k in seq_len(nrow(T_mat))) {
        for (l in seq_len(ncol(T_mat))) {
          gw_cost <- gw_cost + (C1[i, k] - C2[j, l])^2 * T_mat[i, j] * T_mat[k, l]
        }
      }
    }
  }
  
  fgw_dist <- (1 - alpha) * ot_cost + alpha * gw_cost
  
  if (verbose) vmessage("  FGW distance: ", round(fgw_dist, 6))
  
  list(
    T = T_mat,
    fgw_dist = fgw_dist,
    ot_cost = ot_cost,
    gw_cost = gw_cost
  )
}

#' Compute KL Divergence Using Python
#'
#' @description Compute KL divergence matrix exactly like Python STRIDE.
#'
#' @param sample1_topic Topic distribution for sample 1 (spots x topics)
#' @param sample2_topic Topic distribution for sample 2 (spots x topics)
#'
#' @return KL divergence matrix
#' @export
compute_kl_divergence_python <- function(sample1_topic, sample2_topic) {
  
  np <- reticulate::import("numpy")
  
  # Normalize rows to sum to 1
  sample1_norm <- sample1_topic / rowSums(sample1_topic)
  sample2_norm <- sample2_topic / rowSums(sample2_topic)
  
  # Exactly like Python STRIDE's KLDivergency function
  sample1_log <- log(sample1_norm)
  sample2_log <- log(sample2_norm)
  
  # sample1_log_sample1[i] = sum_k p[i,k] * log(p[i,k])
  sample1_log_sample1 <- rowSums(sample1_norm * sample1_log)
  
  # sample1_log_sample2[i,j] = sum_k p[i,k] * log(q[j,k])
  sample1_log_sample2 <- as.matrix(sample1_norm) %*% t(sample2_log)
  
  # KL(P||Q) = sum_k p_k * log(p_k) - sum_k p_k * log(q_k)
  kl_divergence <- matrix(sample1_log_sample1, nrow = length(sample1_log_sample1), 
                          ncol = nrow(sample2_norm), byrow = FALSE) - sample1_log_sample2
  
  kl_divergence
}

# ============================================================================
# Marker Gene Identification via Scanpy
# ============================================================================

#' Find Marker Genes Using Scanpy
#'
#' @description Find marker genes using scanpy for exact equivalence with Python STRIDE.
#'
#' @param matrix Count matrix (genes x cells)
#' @param celltypes Named vector of cell type annotations
#' @param n.top Number of top markers per cell type
#' @param method Test method (default: "wilcoxon")
#' @param tie_correct Use tie correction (default: TRUE)
#' @param pts Calculate percentage of cells expressing (default: TRUE)
#'
#' @return List with markers data frame and top marker list
#' @export
find_markers_scanpy <- function(matrix, celltypes, n.top = 200,
                                 method = "wilcoxon", 
                                 tie_correct = TRUE,
                                 pts = TRUE) {
  
  if (!check_python_backend("scanpy")) {
    rlang::abort("scanpy not available. Run setup_python_backend() first.")
  }
  
  vmessage("Finding markers using scanpy (Python backend)...")
  
  sc <- reticulate::import("scanpy")
  np <- reticulate::import("numpy")
  pd <- reticulate::import("pandas")
  anndata <- reticulate::import("anndata")
  
  # Create AnnData object - cells as rows, genes as columns
  if (is_sparse(matrix)) {
    mat_t <- as.matrix(Matrix::t(matrix))
  } else {
    mat_t <- t(matrix)
  }
  
  # Ensure cell names match celltypes
  common_cells <- intersect(rownames(mat_t), names(celltypes))
  mat_t <- mat_t[common_cells, , drop = FALSE]
  celltypes <- celltypes[common_cells]
  
  # Create AnnData
  adata <- anndata$AnnData(
    X = np$array(mat_t, dtype = "float32"),
    obs = pd$DataFrame(list(Celltype = as.character(celltypes)), 
                       index = names(celltypes)),
    var = pd$DataFrame(index = colnames(mat_t))
  )
  
  # Preprocess exactly like STRIDE
  vmessage("  Preprocessing...")
  sc$pp$filter_cells(adata, min_genes = as.integer(200))
  sc$pp$filter_genes(adata, min_cells = as.integer(10))
  sc$pp$calculate_qc_metrics(adata, inplace = TRUE)
  
  # Normalize
  sc$pp$normalize_total(adata, target_sum = 1e4)
  sc$pp$log1p(adata)
  
  # Find markers - exactly like STRIDE
  vmessage("  Running rank_genes_groups...")
  sc$tl$rank_genes_groups(
    adata,
    groupby = "Celltype",
    method = method,
    pts = pts,
    use_raw = FALSE,
    tie_correct = tie_correct
  )
  
  # Extract results
  marker_names <- adata$uns["rank_genes_groups"]["names"]
  marker_df <- pd$DataFrame(marker_names)
  marker_df_r <- reticulate::py_to_r(marker_df)
  
  # Get top markers - exactly like STRIDE
  top_markers <- unique(as.vector(t(as.matrix(marker_df_r[seq_len(n.top), ]))))
  top_markers <- top_markers[!is.na(top_markers)]
  
  vmessage("  Found ", length(top_markers), " unique marker genes")
  
  list(
    markers_df = marker_df_r,
    top_markers = top_markers
  )
}

# ============================================================================
# High-Level Integration Functions
# ============================================================================
#' Align Two Samples Using POT FGW
#'
#' @description Align two samples using POT's FGW implementation.
#'
#' @param topic1 Topic distribution for sample 1 (topics x spots)
#' @param topic2 Topic distribution for sample 2 (topics x spots)
#' @param coords1 Coordinates for sample 1
#' @param coords2 Coordinates for sample 2
#' @param alpha Trade-off parameter (default: 0.5)
#'
#' @return List with transport plan, rotation matrix, and FGW distance
#' @export
align_samples_pot <- function(topic1, topic2, coords1, coords2, alpha = 0.5) {
  
  vmessage("Aligning samples using POT...")
  
  # Transpose topics to spots x topics
  topic1_t <- t(topic1)
  topic2_t <- t(topic2)
  
  n1 <- nrow(topic1_t)
  n2 <- nrow(topic2_t)
  
  # Uniform distributions
  p <- rep(1/n1, n1)
  q <- rep(1/n2, n2)
  
  # KL divergence matrix (exactly like Python)
  M <- compute_kl_divergence_python(topic1_t + 1e-10, topic2_t + 1e-10)
  
  # Distance matrices
  scipy <- reticulate::import("scipy.spatial")
  C1 <- as.matrix(scipy$distance_matrix(as.matrix(coords1[, 1:2]), 
                                         as.matrix(coords1[, 1:2])))
  C2 <- as.matrix(scipy$distance_matrix(as.matrix(coords2[, 1:2]), 
                                         as.matrix(coords2[, 1:2])))
  
  # FGW using POT
  fgw_result <- fused_gromov_wasserstein_pot(M, C1, C2, p, q, alpha = alpha)
  
  # Compute rotation
  c1_centered <- center_coords(coords1, p)
  c2_centered <- center_coords(coords2, q)
  
  # SVD for rotation - exactly like Python
  np <- reticulate::import("numpy")
  H <- t(as.matrix(c2_centered)) %*% t(fgw_result$T) %*% as.matrix(c1_centered)
  svd_result <- svd(H)
  R <- svd_result$u %*% t(svd_result$v)
  
  # Ensure proper rotation
  if (det(R) < 0) {
    svd_result$v[, 2] <- -svd_result$v[, 2]
    R <- svd_result$u %*% t(svd_result$v)
  }
  
  list(
    transport = fgw_result$T,
    rotation = R,
    fgw_dist = fgw_result$fgw_dist
  )
}

# ============================================================================
# Unified Interface
# ============================================================================

#' Train LDA Model (Unified Interface)
#'
#' @description Train LDA model with option to use Python backend.
#'
#' @param matrix Count matrix (genes x cells)
#' @param n.topics Number of topics
#' @param backend Backend to use: "auto", "python", or "r"
#' @param ... Additional arguments passed to backend
#'
#' @return LDA model object
#' @export
train_lda_unified <- function(matrix, n.topics = 10, backend = "auto", ...) {
  
  if (backend == "auto") {
    backend <- if (check_python_backend(c("gensim", "numpy"))) "python" else "r"
    vmessage("Using backend: ", backend)
  }
  
  if (backend == "python") {
    train_lda_gensim(matrix, n.topics = n.topics, ...)
  } else {
    train_lda(matrix, n.topics = n.topics, ...)
  }
}

#' Transform Data with LDA (Unified Interface)
#'
#' @description Transform data using trained LDA model.
#'
#' @param model Trained LDA model
#' @param newdata New data matrix
#' @param ... Additional arguments
#'
#' @return Topic distribution matrix
#' @export
transform_lda_unified <- function(model, newdata, ...) {
  
  if (inherits(model, "strider_lda_gensim")) {
    transform_lda_gensim(model, newdata, ...)
  } else {
    transform_lda(model, newdata, ...)
  }
}

#' Calculate Coherence (Unified Interface)
#'
#' @description Calculate coherence with option to use Python backend.
#'
#' @param model LDA model
#' @param matrix Count matrix
#' @param backend Backend: "auto", "python", or "r"
#' @param ... Additional arguments
#'
#' @return List with coherence scores
#' @export
calculate_coherence_unified <- function(model, matrix, backend = "auto", ...) {
  
  if (backend == "auto") {
    backend <- if (inherits(model, "strider_lda_gensim")) "python" else "r"
  }
  
  if (backend == "python" && inherits(model, "strider_lda_gensim")) {
    calculate_coherence_gensim(model, matrix, ...)
  } else {
    # R backend only has UMass
    result <- calculate_coherence(model, matrix, ...)
    list(umass = result$mean, cv = NA)
  }
}

#' Find Markers (Unified Interface)
#'
#' @description Find marker genes with option to use Python backend.
#'
#' @param matrix Count matrix
#' @param celltypes Cell type annotations
#' @param backend Backend: "auto", "python", or "r"
#' @param ... Additional arguments
#'
#' @return List with markers
#' @export
find_markers_unified <- function(matrix, celltypes, backend = "auto", ...) {
  
  if (backend == "auto") {
    backend <- if (check_python_backend("scanpy")) "python" else "r"
    vmessage("Using backend: ", backend)
  }
  
  if (backend == "python") {
    find_markers_scanpy(matrix, celltypes, ...)
  } else {
    find_markers(matrix, celltypes, ...)
  }
}

#' Align Samples (Unified Interface)
#'
#' @description Align two samples with option to use Python backend.
#'
#' @param topic1 Topic distribution for sample 1
#' @param topic2 Topic distribution for sample 2
#' @param coords1 Coordinates for sample 1
#' @param coords2 Coordinates for sample 2
#' @param backend Backend: "auto", "python", or "r"
#' @param ... Additional arguments
#'
#' @return Alignment result
#' @export
align_samples_unified <- function(topic1, topic2, coords1, coords2, 
                                   backend = "auto", ...) {
  
  if (backend == "auto") {
    backend <- if (check_python_backend("ot")) "python" else "r"
    vmessage("Using backend: ", backend)
  }
  
  if (backend == "python") {
    align_samples_pot(topic1, topic2, coords1, coords2, ...)
  } else {
    align_samples(topic1, topic2, coords1, coords2, ...)
  }
}
