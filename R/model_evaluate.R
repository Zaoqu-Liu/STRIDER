#' @title Model Evaluation Functions for STRIDER
#' @description Functions for evaluating and selecting the best LDA model
#' @name model_evaluate
NULL

#' Evaluate Model Using Cell Type Prediction
#'
#' @description Evaluate LDA model quality by predicting cell types and 
#'   comparing with known annotations.
#'
#' @param topic_cell Topic-cell matrix (topics x cells)
#' @param topic_celltype Topic-celltype matrix (topics x celltypes)
#' @param true_celltypes Named vector of true cell type annotations
#' @param method Deconvolution method: "Raw", "Norm", "NormBySD", "Bayes", "BayesNorm"
#' @param celltype_topic_bayes Bayesian celltype-topic matrix (for Bayes methods)
#'
#' @return A list containing:
#'   \itemize{
#'     \item{accuracy}{Prediction accuracy}
#'     \item{predictions}{Predicted cell types}
#'     \item{confusion}{Confusion matrix}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' eval <- evaluate_model_prediction(topic_cell, topic_ct, true_ct)
#' }
evaluate_model_prediction <- function(topic_cell, topic_celltype, true_celltypes,
                                       method = "Raw", celltype_topic_bayes = NULL) {
  
  vmessage("Evaluating model with method: ", method)
  
  # Ensure alignment
  common_cells <- intersect(colnames(topic_cell), names(true_celltypes))
  topic_cell <- topic_cell[, common_cells]
  true_celltypes <- true_celltypes[common_cells]
  
  # Transpose topic_celltype to get celltype_topic
  celltype_topic <- t(topic_celltype)
  
  # Apply normalization based on method
  if (method == "Raw") {
    ct_weights <- celltype_topic
  } else if (method == "Norm") {
    # Normalize by column (topic) sums
    col_sums <- colSums(celltype_topic)
    col_sums[col_sums == 0] <- 1
    ct_weights <- sweep(celltype_topic, 2, col_sums, "/")
  } else if (method == "NormBySD") {
    # Normalize by standard deviation
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
    # Normalize Bayes matrix by column sums
    col_sums <- colSums(celltype_topic_bayes)
    col_sums[col_sums == 0] <- 1
    ct_weights <- sweep(celltype_topic_bayes, 2, col_sums, "/")
  } else {
    rlang::abort(paste0("Unknown method: ", method))
  }
  
  # Ensure ct_weights has correct dimensions
  n_topics <- nrow(topic_cell)
  if (ncol(ct_weights) > n_topics) {
    ct_weights <- ct_weights[, seq_len(n_topics), drop = FALSE]
  }
  
  # Calculate cell-celltype scores: (celltypes x topics) %*% (topics x cells)
  celltype_cell <- ct_weights %*% topic_cell
  
  # Normalize to fractions
  col_sums <- colSums(celltype_cell)
  col_sums[col_sums == 0] <- 1
  celltype_cell_frac <- sweep(celltype_cell, 2, col_sums, "/")
  
  # Predict cell types (maximum score)
  predictions <- rownames(celltype_cell_frac)[apply(celltype_cell_frac, 2, which.max)]
  names(predictions) <- colnames(celltype_cell_frac)
  
  # Calculate accuracy
  correct <- predictions == true_celltypes
  accuracy <- sum(correct) / length(correct)
  
  # Create confusion matrix
  confusion <- table(Prediction = predictions, Truth = true_celltypes)
  
  vmessage("  Accuracy: ", round(accuracy * 100, 2), "%")
  
  list(
    accuracy = accuracy,
    predictions = predictions,
    celltype_frac = celltype_cell_frac,
    confusion = confusion
  )
}

#' Calculate Normalized Mutual Information
#'
#' @description Calculate NMI between predicted clusters and true cell types.
#'
#' @param pred Predicted labels
#' @param truth True labels
#'
#' @return NMI score between 0 and 1
#' @export
calculate_nmi <- function(pred, truth) {
  # Ensure same length
  if (length(pred) != length(truth)) {
    rlang::abort("pred and truth must have same length")
  }
  
  # Convert to factors
  pred <- as.factor(pred)
  truth <- as.factor(truth)
  
  n <- length(pred)
  
  # Contingency table
  cont_table <- table(pred, truth)
  
  # Marginal probabilities
  p_pred <- rowSums(cont_table) / n
  p_truth <- colSums(cont_table) / n
  
  # Joint probabilities
  p_joint <- cont_table / n
  
  # Entropy
  H_pred <- -sum(p_pred[p_pred > 0] * log(p_pred[p_pred > 0]))
  H_truth <- -sum(p_truth[p_truth > 0] * log(p_truth[p_truth > 0]))
  
  # Mutual information
  MI <- 0
  for (i in seq_len(nrow(cont_table))) {
    for (j in seq_len(ncol(cont_table))) {
      if (p_joint[i, j] > 0) {
        MI <- MI + p_joint[i, j] * log(p_joint[i, j] / (p_pred[i] * p_truth[j]))
      }
    }
  }
  
  # Normalized MI
  if (H_pred + H_truth == 0) {
    return(0)
  }
  
  NMI <- 2 * MI / (H_pred + H_truth)
  
  NMI
}

#' Perform K-Means Clustering on Topic Space
#'
#' @description Cluster cells in topic space and evaluate clustering quality.
#'
#' @param topic_cell Topic-cell matrix (topics x cells)
#' @param true_celltypes Named vector of true cell type annotations
#' @param n_clusters Number of clusters (default: auto from celltypes)
#' @param n_iter Maximum k-means iterations (default: 1000)
#' @param seed Random seed
#'
#' @return A list containing cluster assignments and NMI score
#' @export
cluster_in_topic_space <- function(topic_cell, true_celltypes, 
                                    n_clusters = NULL, n_iter = 1000, seed = 42) {
  
  # Determine number of clusters
  if (is.null(n_clusters)) {
    n_clusters <- length(unique(true_celltypes))
  }
  
  # Transpose to cells x topics
  cell_topic <- t(topic_cell)
  
  # K-means clustering
  set.seed(seed)
  km <- stats::kmeans(cell_topic, centers = n_clusters, iter.max = n_iter, nstart = 10)
  
  clusters <- km$cluster[names(true_celltypes)]
  
  # Calculate NMI
  nmi <- calculate_nmi(clusters, true_celltypes)
  
  list(
    clusters = clusters,
    nmi = nmi,
    kmeans = km
  )
}

#' Evaluate Multiple Topic Numbers
#'
#' @description Train and evaluate LDA models with different numbers of topics
#'   to find the optimal configuration.
#'
#' @param matrix Normalized count matrix (genes x cells)
#' @param celltypes Named vector of cell type annotations
#' @param topic_range Vector of topic numbers to test
#' @param methods Vector of evaluation methods
#' @param n.iter LDA iterations (default: 500)
#' @param seed Random seed
#' @param save.models Save intermediate models (default: FALSE)
#' @param out.dir Output directory for models
#' @param backend Computational backend: "auto", "python", or "r" (default: "auto").
#'   When "auto", uses Python gensim if available for exact equivalence
#'   with the original STRIDE implementation.
#'
#' @return A list containing:
#'   \itemize{
#'     \item{metrics}{Data frame of evaluation metrics}
#'     \item{best_model}{Best LDA model}
#'     \item{best_n_topics}{Optimal number of topics}
#'     \item{best_method}{Best evaluation method}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' results <- evaluate_topic_range(mat, celltypes, topic_range = 5:15)
#' }
evaluate_topic_range <- function(matrix, celltypes, topic_range,
                                  methods = c("Raw", "Norm", "NormBySD", "Bayes", "BayesNorm"),
                                  n.iter = 500, seed = 42,
                                  save.models = FALSE, out.dir = NULL,
                                  backend = "auto") {
  
  vmessage("Evaluating ", length(topic_range), " topic configurations...")
  
  # Select backend
  if (backend == "auto") {
    backend <- if (check_python_backend(c("gensim", "numpy"))) "python" else "r"
  }
  
  # Select training function based on backend
  train_lda_fn <- if (backend == "python" && check_python_backend(c("gensim", "numpy"))) {
    train_lda_gensim
  } else {
    train_lda
  }
  
  # Store results
  all_metrics <- list()
  all_models <- list()
  all_predictions <- list()
  
  # Cell type counts for Bayesian methods
  celltype_counts <- table(celltypes)
  total_cells <- length(celltypes)
  
  for (n_topics in topic_range) {
    vmessage("\n=== Number of topics: ", n_topics, " ===")
    
    # Train model
    model <- train_lda_fn(matrix, n.topics = n_topics, n.iter = n.iter, seed = seed)
    
    # Calculate topic-celltype association
    topic_celltype <- calculate_topic_celltype(model$topic_cell, celltypes)
    
    # Calculate Bayesian association
    celltype_topic_bayes <- calculate_celltype_topic_bayes(
      topic_celltype, celltype_counts, total_cells
    )
    
    # Calculate coherence (use gensim if available for both u_mass and c_v)
    if (backend == "python" && inherits(model, "strider_lda_gensim")) {
      coherence <- calculate_coherence_gensim(model, matrix, coherence_type = "both")
    } else {
      coherence <- calculate_coherence(model, matrix, type = "both")
    }
    
    # Evaluate with different methods
    method_accuracies <- numeric(length(methods))
    names(method_accuracies) <- methods
    
    predictions_by_method <- list()
    
    for (method in methods) {
      eval_result <- evaluate_model_prediction(
        model$topic_cell, topic_celltype, celltypes,
        method = method, celltype_topic_bayes = celltype_topic_bayes
      )
      method_accuracies[method] <- eval_result$accuracy
      predictions_by_method[[method]] <- eval_result$predictions
    }
    
    # K-means clustering NMI
    cluster_result <- cluster_in_topic_space(model$topic_cell, celltypes)
    
    # Store metrics - include both umass and cv coherence
    metrics <- data.frame(
      n_topics = n_topics,
      umass_coherence = coherence$umass,
      cv_coherence = if (!is.null(coherence$cv) && !is.na(coherence$cv)) coherence$cv else NA,
      nmi = cluster_result$nmi,
      t(method_accuracies),
      stringsAsFactors = FALSE
    )
    colnames(metrics)[5:ncol(metrics)] <- paste0(methods, "_accuracy")
    
    # Keep backward compatibility
    metrics$coherence <- metrics$umass_coherence
    
    all_metrics[[as.character(n_topics)]] <- metrics
    all_models[[as.character(n_topics)]] <- list(
      model = model,
      topic_celltype = topic_celltype,
      celltype_topic_bayes = celltype_topic_bayes
    )
    all_predictions[[as.character(n_topics)]] <- predictions_by_method
    
    # Save model if requested
    if (save.models && !is.null(out.dir)) {
      model_dir <- file.path(out.dir, "models")
      if (!dir.exists(model_dir)) dir.create(model_dir, recursive = TRUE)
      saveRDS(model, file.path(model_dir, paste0("lda_model_", n_topics, ".rds")))
    }
  }
  
  # Combine metrics
  metrics_df <- do.call(rbind, all_metrics)
  rownames(metrics_df) <- NULL
  
  # Find best model (highest mean accuracy across methods)
  accuracy_cols <- grep("_accuracy$", colnames(metrics_df), value = TRUE)
  metrics_df$mean_accuracy <- rowMeans(metrics_df[, accuracy_cols])
  
  best_idx <- which.max(metrics_df$mean_accuracy)
  best_n_topics <- metrics_df$n_topics[best_idx]
  
  # Find best method (highest mean accuracy across topic numbers)
  method_means <- colMeans(metrics_df[, accuracy_cols])
  best_method <- gsub("_accuracy$", "", names(which.max(method_means)))
  
  vmessage("\n=== Model Selection Results ===")
  vmessage("Best number of topics: ", best_n_topics)
  vmessage("Best method: ", best_method)
  vmessage("Best accuracy: ", round(max(metrics_df$mean_accuracy) * 100, 2), "%")
  
  list(
    metrics = metrics_df,
    best_model = all_models[[as.character(best_n_topics)]]$model,
    best_topic_celltype = all_models[[as.character(best_n_topics)]]$topic_celltype,
    best_celltype_topic_bayes = all_models[[as.character(best_n_topics)]]$celltype_topic_bayes,
    best_n_topics = best_n_topics,
    best_method = best_method,
    all_models = all_models,
    predictions = all_predictions
  )
}

#' Select Best Model
#'
#' @description Automatically select the best model from evaluation results.
#'
#' @param eval_results Results from evaluate_topic_range()
#' @param criterion Selection criterion: "accuracy", "coherence", or "balanced"
#'
#' @return Selected model components
#' @export
select_best_model <- function(eval_results, criterion = "accuracy") {
  
  metrics <- eval_results$metrics
  
  if (criterion == "accuracy") {
    best_idx <- which.max(metrics$mean_accuracy)
  } else if (criterion == "coherence") {
    best_idx <- which.max(metrics$coherence)
  } else if (criterion == "balanced") {
    # Normalize both metrics and combine
    norm_acc <- (metrics$mean_accuracy - min(metrics$mean_accuracy)) / 
                (max(metrics$mean_accuracy) - min(metrics$mean_accuracy) + 1e-10)
    norm_coh <- (metrics$coherence - min(metrics$coherence)) / 
                (max(metrics$coherence) - min(metrics$coherence) + 1e-10)
    combined <- norm_acc + norm_coh
    best_idx <- which.max(combined)
  }
  
  best_n_topics <- metrics$n_topics[best_idx]
  
  list(
    n_topics = best_n_topics,
    model = eval_results$all_models[[as.character(best_n_topics)]]$model,
    topic_celltype = eval_results$all_models[[as.character(best_n_topics)]]$topic_celltype,
    celltype_topic_bayes = eval_results$all_models[[as.character(best_n_topics)]]$celltype_topic_bayes,
    method = eval_results$best_method
  )
}
