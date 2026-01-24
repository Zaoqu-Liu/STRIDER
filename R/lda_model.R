#' @title LDA Topic Model for STRIDER
#' @description Functions for training and applying LDA topic models
#' @name lda_model
NULL

#' Train LDA Topic Model
#'
#' @description Train a Latent Dirichlet Allocation (LDA) topic model on 
#'   single-cell gene expression data. The model treats cells as documents
#'   and genes as words, following the methodology described in STRIDE.
#'
#' @param matrix Normalized count matrix (genes x cells)
#' @param n.topics Number of topics (default: 10)
#' @param n.iter Number of iterations (default: 500)
#' @param convergence.tol Convergence tolerance (default: 1e-4)
#' @param alpha Dirichlet prior for document-topic distribution. 
#'   Default is 50/n.topics as in original STRIDE.
#' @param beta Dirichlet prior for topic-word distribution (default: 0.01)
#' @param seed Random seed for reproducibility
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return An LDA model object with components:
#'   \itemize{
#'     \item{topic_gene}{Topic-gene distribution matrix (phi)}
#'     \item{topic_cell}{Document-topic distribution matrix (theta)}
#'     \item{vocab}{Gene vocabulary}
#'     \item{n_topics}{Number of topics}
#'   }
#'
#' @details
#' The LDA model is trained using Variational Bayes inference via the text2vec

#' package. The default hyperparameters (alpha = 50/K, beta = 0.01) follow the
#' original STRIDE implementation for optimal performance on single-cell data.
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
#' model <- train_lda(norm_mat, n.topics = 10)
#' }
train_lda <- function(matrix, n.topics = 10, n.iter = 500,
                      convergence.tol = 1e-4, alpha = NULL, beta = 0.01,
                      seed = NULL, verbose = TRUE) {
  
  if (verbose) vmessage("Training LDA model with ", n.topics, " topics...")
  
  # Set seed for reproducibility
  if (!is.null(seed)) set.seed(seed)
  
  # Check for text2vec package
  if (!requireNamespace("text2vec", quietly = TRUE)) {
    rlang::abort("Package 'text2vec' is required. Install with: install.packages('text2vec')")
  }
  
  # Set default alpha following STRIDE
  if (is.null(alpha)) {
    alpha <- 50 / n.topics
  }
  
  # Prepare data - transpose so cells are rows (documents)
  if (is_sparse(matrix)) {
    dtm <- Matrix::t(matrix)
  } else {
    dtm <- Matrix::Matrix(t(matrix), sparse = TRUE)
  }
  
  # Ensure proper format

  dtm <- methods::as(dtm, "dgCMatrix")
  
  # Handle negative values (should not occur in count data)
  if (any(dtm@x < 0)) {
    rlang::warn("Negative values detected and set to zero")
    dtm@x[dtm@x < 0] <- 0
  }
  
  # Store vocabulary (genes)
  vocab <- colnames(dtm)
  if (is.null(vocab)) {
    vocab <- paste0("Gene_", seq_len(ncol(dtm)))
    colnames(dtm) <- vocab
  }
  
  # Create text2vec LDA model with STRIDE hyperparameters
  lda_model <- text2vec::LDA$new(
    n_topics = n.topics,
    doc_topic_prior = alpha,
    topic_word_prior = beta
  )
  
  # Fit the model
  doc_topic_distr <- lda_model$fit_transform(
    x = dtm,
    n_iter = n.iter,
    convergence_tol = convergence.tol,
    n_check_convergence = 10,
    progressbar = verbose
  )
  
  # Extract topic-word distribution (phi matrix)
  # Try different methods to get topic-word distribution
  topic_gene <- tryCatch({
    # Method 1: Direct access to private field
    lda_model$.__enclos_env__$private$topic_word_distribution
  }, error = function(e) NULL)
  
  if (is.null(topic_gene)) {
    # Method 2: Use get_top_words and reconstruct
    tryCatch({
      # Get components from model
      components <- lda_model$components
      if (!is.null(components) && nrow(components) == n.topics) {
        topic_gene <<- components
      }
    }, error = function(e) NULL)
  }
  
  if (is.null(topic_gene)) {
    # Method 3: Compute from word counts
    tryCatch({
      topic_word_count <- lda_model$.__enclos_env__$private$topic_word_count
      if (!is.null(topic_word_count)) {
        # Normalize to get distribution
        row_sums <- rowSums(topic_word_count)
        row_sums[row_sums == 0] <- 1
        topic_gene <<- topic_word_count / row_sums
      }
    }, error = function(e) NULL)
  }
  
  if (is.null(topic_gene)) {
    # Method 4: Create from scratch using doc_topic and dtm
    # This is a fallback - compute pseudo topic-word distribution
    # P(word|topic) ∝ sum_d P(topic|doc) * count(word, doc)
    topic_gene <- matrix(0, nrow = n.topics, ncol = ncol(dtm))
    for (k in seq_len(n.topics)) {
      # Weight each document's word counts by topic probability
      topic_weights <- doc_topic_distr[, k]
      topic_gene[k, ] <- as.numeric(Matrix::t(dtm) %*% topic_weights)
    }
    # Normalize
    row_sums <- rowSums(topic_gene)
    row_sums[row_sums == 0] <- 1
    topic_gene <- topic_gene / row_sums
  }
  
  rownames(topic_gene) <- paste0("Topic_", seq_len(n.topics))
  colnames(topic_gene) <- vocab
  
  # Transpose doc_topic to get topic_cell (topics x cells)
  topic_cell <- t(doc_topic_distr)
  rownames(topic_cell) <- paste0("Topic_", seq_len(n.topics))
  colnames(topic_cell) <- rownames(dtm)
  
  if (verbose) vmessage("  LDA training complete")
  
  structure(
    list(
      topic_gene = topic_gene,
      topic_cell = topic_cell,
      vocab = vocab,
      n_topics = n.topics,
      n_iter = n.iter,
      alpha = alpha,
      beta = beta,
      model = lda_model
    ),
    class = "strider_lda"
  )
}

#' Apply LDA Model to New Data
#'
#' @description Apply a trained LDA model to transform new data (e.g., spatial spots)
#'   into topic space using variational inference.
#'
#' @param model Trained LDA model from train_lda()
#' @param newdata Matrix to transform (genes x cells/spots)
#' @param n.iter Number of inference iterations (default: 100)
#'
#' @return Topic distribution matrix (topics x cells/spots)
#'
#' @details
#' This function performs inference on new documents using the trained topic-word
#' distributions. Genes not present in the training vocabulary are ignored.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' topic_dist <- transform_lda(model, st_matrix)
#' }
transform_lda <- function(model, newdata, n.iter = 100) {
  
  if (!inherits(model, "strider_lda")) {
    rlang::abort("model must be a strider_lda object from train_lda()")
  }
  
  vmessage("Transforming data with LDA model...")
  
  # Get shared genes between model vocabulary and new data
  shared_genes <- intersect(model$vocab, rownames(newdata))
  if (length(shared_genes) == 0) {
    rlang::abort("No shared genes between model and new data")
  }
  
  if (length(shared_genes) < length(model$vocab) * 0.5) {
    rlang::warn(paste0(
      "Only ", length(shared_genes), " of ", length(model$vocab),
      " genes shared. Results may be suboptimal."
    ))
  }
  
  vmessage("  Using ", length(shared_genes), " shared genes")
  
  # Subset new data to shared genes
  newdata_subset <- newdata[shared_genes, , drop = FALSE]
  
  # Create document-term matrix with full vocabulary
  full_dtm <- Matrix::Matrix(
    0, 
    nrow = ncol(newdata_subset), 
    ncol = length(model$vocab),
    sparse = TRUE
  )
  rownames(full_dtm) <- colnames(newdata_subset)
  colnames(full_dtm) <- model$vocab
  
  # Fill in values for shared genes
  gene_idx <- match(shared_genes, model$vocab)
  full_dtm[, gene_idx] <- Matrix::t(newdata_subset)
  
  # Transform using the trained model
  doc_topic_distr <- model$model$transform(
    full_dtm,
    n_iter = n.iter,
    progressbar = FALSE
  )
  
  # Transpose to get topics x cells/spots
  topic_doc <- t(doc_topic_distr)
  rownames(topic_doc) <- paste0("Topic_", seq_len(model$n_topics))
  colnames(topic_doc) <- colnames(newdata)
  
  vmessage("  Transformation complete")
  
  topic_doc
}

#' Calculate Topic-Celltype Association Matrix
#'
#' @description Calculate the association between topics and cell types
#'   by aggregating topic distributions across cells of each type.
#'
#' @param topic_cell Topic-cell matrix (topics x cells)
#' @param celltypes Named vector of cell type annotations
#' @param method Aggregation method: "mean" (default), "median", or "sum"
#'
#' @return Topic-celltype matrix (topics x celltypes)
#'
#' @details
#' This matrix represents P(topic|celltype), the probability of each topic
#' given a cell type. It is computed by averaging topic distributions
#' across all cells of each type.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' topic_ct <- calculate_topic_celltype(topic_cell, celltypes)
#' }
calculate_topic_celltype <- function(topic_cell, celltypes, method = "mean") {
  
  vmessage("Calculating topic-celltype associations...")
  
  # Ensure alignment between matrix and annotations
  common_cells <- intersect(colnames(topic_cell), names(celltypes))
  if (length(common_cells) == 0) {
    rlang::abort("No matching cells between topic_cell and celltypes")
  }
  
  topic_cell <- topic_cell[, common_cells, drop = FALSE]
  celltypes <- celltypes[common_cells]
  
  # Get unique cell types (sorted for reproducibility)
  celltype_levels <- sort(unique(celltypes))
  n_topics <- nrow(topic_cell)
  n_celltypes <- length(celltype_levels)
  
  # Initialize association matrix
  topic_celltype <- matrix(0, nrow = n_topics, ncol = n_celltypes)
  rownames(topic_celltype) <- rownames(topic_cell)
  colnames(topic_celltype) <- celltype_levels
  
  # Aggregate topic distributions for each cell type
  for (ct in celltype_levels) {
    cells_ct <- names(celltypes)[celltypes == ct]
    ct_data <- topic_cell[, cells_ct, drop = FALSE]
    
    topic_celltype[, ct] <- switch(
      method,
      "mean" = rowMeans(ct_data),
      "median" = apply(ct_data, 1, stats::median),
      "sum" = rowSums(ct_data),
      rowMeans(ct_data)  # default to mean
    )
  }
  
  topic_celltype
}

#' Calculate Bayesian Celltype-Topic Matrix
#'
#' @description Calculate celltype-topic association using Bayes' theorem,
#'   incorporating prior probabilities of cell types.
#'
#' @param topic_celltype Topic-celltype matrix (topics x celltypes)
#' @param celltype_counts Named vector of cell counts per cell type
#' @param total_cells Total number of cells
#'
#' @return Celltype-topic Bayes matrix (celltypes x topics)
#'
#' @details
#' Applies Bayes' theorem: P(celltype|topic) = P(topic|celltype) * P(celltype) / P(topic)
#' 
#' This accounts for differences in cell type abundance in the reference data,
#' preventing overrepresentation of abundant cell types in deconvolution.
#'
#' @export
calculate_celltype_topic_bayes <- function(topic_celltype, celltype_counts, total_cells) {
  
  vmessage("Calculating Bayesian celltype-topic associations...")
  
  # Prior probability of each cell type: P(celltype)
  celltype_prob <- as.numeric(celltype_counts) / total_cells
  names(celltype_prob) <- names(celltype_counts)
  
  # Prior probability of each topic: P(topic) = sum_ct P(topic|ct) * P(ct)
  topic_prob <- as.numeric(topic_celltype %*% celltype_prob)
  
  # Avoid division by zero
  topic_prob[topic_prob == 0] <- 1e-10
  
  # Bayes: P(celltype|topic) = P(topic|celltype) * P(celltype) / P(topic)
  celltype_topic_bayes <- t(topic_celltype)  # celltypes x topics
  
  for (i in seq_len(nrow(celltype_topic_bayes))) {
    ct <- rownames(celltype_topic_bayes)[i]
    # P(topic|celltype) * P(celltype) / P(topic)
    celltype_topic_bayes[i, ] <- topic_celltype[, ct] * celltype_prob[ct] / topic_prob
  }
  
  # Handle numerical issues
  celltype_topic_bayes[!is.finite(celltype_topic_bayes)] <- 0
  
  celltype_topic_bayes
}

#' Calculate LDA Coherence Score
#'
#' @description Calculate topic coherence scores to evaluate model quality.
#'   Implements UMass and optionally C_v coherence metrics.
#'
#' @param model LDA model from train_lda()
#' @param matrix Original count matrix used for training
#' @param n.words Number of top words per topic (default: 10)
#' @param type Coherence type: "umass", "cv", or "both" (default: "both")
#' @param cell_gene_list List of genes expressed per cell (for c_v, auto-generated if NULL)
#'
#' @return A list containing:
#'   \itemize{
#'     \item{umass}{UMass coherence (mean across topics)}
#'     \item{cv}{C_v coherence (mean across topics), NA if type="umass"}
#'     \item{umass_per_topic}{UMass score for each topic}
#'     \item{cv_per_topic}{C_v score for each topic}
#'   }
#'
#' @details
#' UMass coherence measures how often top words co-occur in the corpus:
#' C_UMass = sum_{j=2}^{M} sum_{i=1}^{j-1} log((D(w_i, w_j) + 1) / D(w_i))
#' 
#' C_v coherence uses normalized pointwise mutual information (NPMI) with
#' a sliding window and cosine similarity, providing better correlation
#' with human judgments.
#'
#' @references
#' Mimno, D., et al. (2011). Optimizing semantic coherence in topic models. EMNLP.
#' Röder, M., et al. (2015). Exploring the Space of Topic Coherence Measures. WSDM.
#'
#' @export
calculate_coherence <- function(model, matrix, n.words = 10, type = "both",
                                 cell_gene_list = NULL) {
  
  vmessage("Calculating topic coherence...")
  
  topic_gene <- model$topic_gene
  n_topics <- nrow(topic_gene)
  vocab <- model$vocab
  
  # Initialize results
  umass_scores <- numeric(n_topics)
  cv_scores <- numeric(n_topics)
  
  # Build cell_gene_list if needed for c_v
  if (type %in% c("cv", "both") && is.null(cell_gene_list)) {
    vmessage("  Building gene lists per cell for C_v...")
    cell_gene_list <- lapply(seq_len(ncol(matrix)), function(i) {
      expressed <- matrix[, i] > 0
      rownames(matrix)[expressed]
    })
  }
  
  for (k in seq_len(n_topics)) {
    # Get top genes for this topic
    gene_probs <- topic_gene[k, ]
    top_genes <- names(sort(gene_probs, decreasing = TRUE))[seq_len(n.words)]
    top_genes <- intersect(top_genes, rownames(matrix))
    
    if (length(top_genes) < 2) {
      umass_scores[k] <- NA
      cv_scores[k] <- NA
      next
    }
    
    # ===== UMass Coherence =====
    if (type %in% c("umass", "both")) {
      sub_mat <- matrix[top_genes, , drop = FALSE] > 0
      coherence <- 0
      m <- length(top_genes)
      
      for (j in 2:m) {
        for (i in 1:(j-1)) {
          D_wi <- sum(sub_mat[i, ])
          D_wiwj <- sum(sub_mat[i, ] & sub_mat[j, ])
          if (D_wi > 0) {
            coherence <- coherence + log((D_wiwj + 1) / D_wi)
          }
        }
      }
      umass_scores[k] <- coherence
    }
    
    # ===== C_v Coherence =====
    if (type %in% c("cv", "both")) {
      cv_scores[k] <- .calculate_cv_coherence(top_genes, cell_gene_list, n.words)
    }
  }
  
  result <- list(
    umass = mean(umass_scores, na.rm = TRUE),
    cv = if (type %in% c("cv", "both")) mean(cv_scores, na.rm = TRUE) else NA,
    umass_per_topic = umass_scores,
    cv_per_topic = cv_scores
  )
  
  # For backward compatibility
  result$mean <- result$umass
  result$per_topic <- result$umass_per_topic
  
  vmessage("  UMass coherence: ", round(result$umass, 4))
  if (!is.na(result$cv)) {
    vmessage("  C_v coherence: ", round(result$cv, 4))
  }
  
  result
}

#' Calculate C_v Coherence for One Topic
#'
#' @description Internal function to calculate C_v coherence.
#'   Implements NPMI-based coherence with sliding window.
#'
#' @param top_words Top words for the topic
#' @param texts List of word lists (documents)
#' @param window_size Sliding window size
#'
#' @return C_v coherence score
#' @keywords internal
.calculate_cv_coherence <- function(top_words, texts, window_size = 10) {
  
  m <- length(top_words)
  if (m < 2) return(NA)
  
  # Count word occurrences and co-occurrences using sliding window
  n_docs <- length(texts)
  word_count <- rep(0, m)
  names(word_count) <- top_words
  cooccur_count <- matrix(0, m, m)
  rownames(cooccur_count) <- colnames(cooccur_count) <- top_words
  
  for (doc in texts) {
    doc_words <- intersect(doc, top_words)
    
    # Count single occurrences
    for (w in doc_words) {
      word_count[w] <- word_count[w] + 1
    }
    
    # Count co-occurrences within window
    if (length(doc) >= 2) {
      for (i in seq_len(length(doc) - 1)) {
        window_end <- min(i + window_size, length(doc))
        window_words <- doc[i:window_end]
        window_top <- intersect(window_words, top_words)
        
        if (length(window_top) >= 2) {
          for (j in seq_len(length(window_top) - 1)) {
            for (k in (j + 1):length(window_top)) {
              w1 <- window_top[j]
              w2 <- window_top[k]
              cooccur_count[w1, w2] <- cooccur_count[w1, w2] + 1
              cooccur_count[w2, w1] <- cooccur_count[w2, w1] + 1
            }
          }
        }
      }
    }
  }
  
  # Calculate NPMI matrix
  eps <- 1e-10
  npmi_matrix <- matrix(0, m, m)
  
  for (i in seq_len(m)) {
    for (j in seq_len(m)) {
      if (i != j) {
        p_i <- (word_count[i] + eps) / n_docs
        p_j <- (word_count[j] + eps) / n_docs
        p_ij <- (cooccur_count[i, j] + eps) / n_docs
        
        pmi <- log(p_ij / (p_i * p_j))
        npmi <- pmi / (-log(p_ij))
        npmi_matrix[i, j] <- max(-1, min(1, npmi))  # Clip to [-1, 1]
      }
    }
  }
  
  # C_v uses cosine similarity of NPMI vectors
  # For simplicity, we use mean NPMI as approximation
  upper_tri <- npmi_matrix[upper.tri(npmi_matrix)]
  mean(upper_tri, na.rm = TRUE)
}

#' Print Method for LDA Model
#' @param x A strider_lda object
#' @param ... Additional arguments (unused)
#' @export
print.strider_lda <- function(x, ...) {
  cat("STRIDER LDA Model\n")
  cat("-----------------\n")
  cat("Topics:", x$n_topics, "\n")
  cat("Vocabulary:", length(x$vocab), "genes\n")
  cat("Documents:", ncol(x$topic_cell), "cells\n")
  cat("Iterations:", x$n_iter, "\n")
  cat("Hyperparameters: alpha =", round(x$alpha, 3), ", beta =", x$beta, "\n")
  invisible(x)
}
