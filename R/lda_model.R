#' @title LDA Topic Model for STRIDER
#' @description Functions for training and applying LDA topic models
#' @name lda_model
NULL

#' Train LDA Topic Model
#'
#' @description Train a Latent Dirichlet Allocation (LDA) topic model on 
#'   single-cell gene expression data. The model treats cells as documents
#'   and genes as words.
#'
#' @param matrix Normalized count matrix (genes x cells)
#' @param n.topics Number of topics (default: 10)
#' @param n.iter Number of iterations (default: 500)
#' @param convergence.tol Convergence tolerance (default: 1e-4)
#' @param seed Random seed for reproducibility
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return An LDA model object with components:
#'   \itemize{
#'     \item{topic_gene}{Topic-gene distribution matrix}
#'     \item{doc_topic}{Document(cell)-topic distribution matrix}
#'     \item{vocab}{Gene vocabulary}
#'     \item{n_topics}{Number of topics}
#'     \item{coherence}{Model coherence scores}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' model <- train_lda(norm_mat, n.topics = 10)
#' }
train_lda <- function(matrix, n.topics = 10, n.iter = 500,
                      convergence.tol = 1e-4, seed = NULL, verbose = TRUE) {
  
  if (verbose) vmessage("Training LDA model with ", n.topics, " topics...")
  
  # Set seed
  if (!is.null(seed)) set.seed(seed)
  
  # Check for text2vec package
  if (!requireNamespace("text2vec", quietly = TRUE)) {
    rlang::abort("Package 'text2vec' is required. Install with: install.packages('text2vec')")
  }
  
  # Prepare data - transpose so cells are rows (documents)
  if (is_sparse(matrix)) {
    dtm <- Matrix::t(matrix)
  } else {
    dtm <- Matrix::Matrix(t(matrix), sparse = TRUE)
  }
  
  # Ensure non-negative integers (scale and round)
  dtm <- methods::as(dtm, "dgCMatrix")
  if (any(dtm@x < 0)) {
    dtm@x[dtm@x < 0] <- 0
  }
  
  # Store vocabulary (genes)
  vocab <- colnames(dtm)
  
  # Create text2vec LDA model
  lda_model <- text2vec::LDA$new(
    n_topics = n.topics,
    doc_topic_prior = 50 / n.topics,  # alpha
    topic_word_prior = 0.01           # beta
  )
  
  # Fit the model
  doc_topic_distr <- lda_model$fit_transform(
    x = dtm,
    n_iter = n.iter,
    convergence_tol = convergence.tol,
    n_check_convergence = 10,
    progressbar = verbose
  )
  
  # Get topic-word distribution
  topic_word_distr <- lda_model$get_top_words(n = length(vocab), lambda = 1)
  
  # Convert to proper matrix format (topics x genes)
  topic_gene <- matrix(0, nrow = n.topics, ncol = length(vocab))
  rownames(topic_gene) <- paste0("Topic_", seq_len(n.topics))
  colnames(topic_gene) <- vocab
  
  # Fill topic-gene matrix from model components
  tw <- lda_model$.__enclos_env__$private$topic_word_distribution
  if (!is.null(tw)) {
    topic_gene <- tw
    rownames(topic_gene) <- paste0("Topic_", seq_len(n.topics))
    colnames(topic_gene) <- vocab
  }
  
  # Transpose doc_topic to get topic_doc (topics x cells)
  topic_doc <- t(doc_topic_distr)
  rownames(topic_doc) <- paste0("Topic_", seq_len(n.topics))
  colnames(topic_doc) <- rownames(dtm)
  
  if (verbose) vmessage("  LDA model training complete")
  
  structure(
    list(
      topic_gene = topic_gene,
      topic_cell = topic_doc,
      vocab = vocab,
      n_topics = n.topics,
      n_iter = n.iter,
      model = lda_model
    ),
    class = "strider_lda"
  )
}

#' Apply LDA Model to New Data
#'
#' @description Apply a trained LDA model to transform new data (e.g., spatial spots)
#'   into topic space.
#'
#' @param model Trained LDA model from train_lda()
#' @param newdata Matrix to transform (genes x cells/spots)
#' @param n.iter Number of inference iterations (default: 100)
#'
#' @return Topic distribution matrix (topics x cells/spots)
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
  
  # Get shared genes
  shared_genes <- intersect(model$vocab, rownames(newdata))
  if (length(shared_genes) == 0) {
    rlang::abort("No shared genes between model and new data")
  }
  
  vmessage("  Using ", length(shared_genes), " shared genes")
  
  # Subset and reorder to match model vocabulary
  newdata_subset <- newdata[shared_genes, , drop = FALSE]
  
  # Create full vocabulary matrix with zeros
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
  
  # Transform using the model
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
#'   based on the topic distribution of cells with known cell types.
#'
#' @param topic_cell Topic-cell matrix (topics x cells)
#' @param celltypes Named vector of cell type annotations
#' @param method Method for calculating association: "mean", "median", or "sum"
#'
#' @return Topic-celltype matrix (topics x celltypes)
#' @export
#'
#' @examples
#' \dontrun{
#' topic_ct <- calculate_topic_celltype(topic_cell, celltypes)
#' }
calculate_topic_celltype <- function(topic_cell, celltypes, method = "mean") {
  
  vmessage("Calculating topic-celltype associations...")
  
  # Ensure alignment
  common_cells <- intersect(colnames(topic_cell), names(celltypes))
  topic_cell <- topic_cell[, common_cells]
  celltypes <- celltypes[common_cells]
  
  # Get unique cell types
  celltype_levels <- sort(unique(celltypes))
  n_topics <- nrow(topic_cell)
  n_celltypes <- length(celltype_levels)
  
  # Calculate association matrix
  topic_celltype <- matrix(0, nrow = n_topics, ncol = n_celltypes)
  rownames(topic_celltype) <- rownames(topic_cell)
  colnames(topic_celltype) <- celltype_levels
  
  for (ct in celltype_levels) {
    cells_ct <- names(celltypes)[celltypes == ct]
    ct_data <- topic_cell[, cells_ct, drop = FALSE]
    
    if (method == "mean") {
      topic_celltype[, ct] <- rowMeans(ct_data)
    } else if (method == "median") {
      topic_celltype[, ct] <- apply(ct_data, 1, stats::median)
    } else if (method == "sum") {
      topic_celltype[, ct] <- rowSums(ct_data)
    }
  }
  
  topic_celltype
}

#' Calculate Celltype-Topic Bayes Matrix
#'
#' @description Calculate Bayesian celltype-topic association using prior
#'   probabilities of cell types.
#'
#' @param topic_celltype Topic-celltype matrix (topics x celltypes)
#' @param celltype_counts Named vector of cell counts per cell type
#' @param total_cells Total number of cells
#'
#' @return Celltype-topic Bayes matrix (celltypes x topics)
#' @export
calculate_celltype_topic_bayes <- function(topic_celltype, celltype_counts, total_cells) {
  
  vmessage("Calculating Bayesian celltype-topic associations...")
  
  # Prior probability of each cell type
  celltype_prob <- celltype_counts / total_cells
  
  # Prior probability of each topic (sum across cell types, normalized)
  topic_prob <- rowSums(topic_celltype)
  topic_prob <- topic_prob / sum(topic_prob)
  
  # Bayes: P(celltype|topic) = P(topic|celltype) * P(celltype) / P(topic)
  # topic_celltype is P(topic|celltype)
  
  celltype_topic_bayes <- t(topic_celltype)  # celltypes x topics
  
  for (i in seq_len(nrow(celltype_topic_bayes))) {
    ct <- rownames(celltype_topic_bayes)[i]
    celltype_topic_bayes[i, ] <- topic_celltype[, ct] * celltype_prob[ct] / topic_prob
  }
  
  # Handle NaN/Inf
  celltype_topic_bayes[!is.finite(celltype_topic_bayes)] <- 0
  
  celltype_topic_bayes
}

#' Calculate LDA Coherence Score
#'
#' @description Calculate topic coherence scores to evaluate model quality.
#'   Implements UMass coherence metric.
#'
#' @param model LDA model from train_lda()
#' @param matrix Original count matrix used for training
#' @param n.words Number of top words per topic for coherence calculation
#'
#' @return Coherence score (higher is better, typically negative)
#' @export
calculate_coherence <- function(model, matrix, n.words = 10) {
  
  vmessage("Calculating topic coherence...")
  
  # Get top words for each topic
  topic_gene <- model$topic_gene
  
  coherence_scores <- numeric(nrow(topic_gene))
  
  for (k in seq_len(nrow(topic_gene))) {
    # Get top genes for this topic
    top_genes <- names(sort(topic_gene[k, ], decreasing = TRUE))[seq_len(n.words)]
    top_genes <- top_genes[top_genes %in% rownames(matrix)]
    
    if (length(top_genes) < 2) {
      coherence_scores[k] <- NA
      next
    }
    
    # Calculate UMass coherence
    # C = sum_{j=2}^M sum_{i=1}^{j-1} log((D(wi, wj) + 1) / D(wi))
    # where D(wi, wj) = count of documents containing both words
    # D(wi) = count of documents containing word wi
    
    # Binarize: present or not
    sub_mat <- matrix[top_genes, , drop = FALSE] > 0
    
    coherence <- 0
    for (j in 2:length(top_genes)) {
      for (i in 1:(j-1)) {
        D_wi <- sum(sub_mat[i, ])
        D_wiwj <- sum(sub_mat[i, ] & sub_mat[j, ])
        if (D_wi > 0) {
          coherence <- coherence + log((D_wiwj + 1) / D_wi)
        }
      }
    }
    
    coherence_scores[k] <- coherence
  }
  
  # Average coherence across topics
  mean_coherence <- mean(coherence_scores, na.rm = TRUE)
  
  vmessage("  Mean coherence (UMass): ", round(mean_coherence, 3))
  
  list(
    mean = mean_coherence,
    per_topic = coherence_scores
  )
}

#' Print Method for LDA Model
#' @export
print.strider_lda <- function(x, ...) {
  cat("STRIDER LDA Model\n")
  cat("-----------------\n")
  cat("Number of topics:", x$n_topics, "\n")
  cat("Vocabulary size:", length(x$vocab), "\n")
  cat("Number of documents:", ncol(x$topic_cell), "\n")
  cat("Training iterations:", x$n_iter, "\n")
  invisible(x)
}
