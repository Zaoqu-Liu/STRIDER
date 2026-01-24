#' @title Scientific Validation Functions for STRIDER
#' @description Functions for validating STRIDER results against benchmarks
#'   and assessing scientific equivalence with Python STRIDE.
#' @name validate
NULL

#' Validate STRIDER Results Against Reference
#'
#' @description Compare STRIDER deconvolution results with reference results
#'   (e.g., from Python STRIDE) to assess scientific equivalence.
#'
#' @param strider_result STRIDER result object or deconvolution matrix (spots x celltypes)
#' @param reference_result Reference result (matrix, data frame, or file path)
#' @param metrics Which metrics to compute (default: all)
#' @param tolerance Numerical tolerance for "exact" comparisons (default: 0.01)
#'
#' @return A validation report object with:
#'   \itemize{
#'     \item{pearson_r}{Pearson correlation of flattened matrices}
#'     \item{spearman_r}{Spearman correlation of flattened matrices}
#'     \item{rmse}{Root mean squared error}
#'     \item{mae}{Mean absolute error}
#'     \item{dominant_match}{Fraction of spots with same dominant cell type}
#'     \item{rank_correlation}{Mean rank correlation per spot}
#'     \item{verdict}{"EQUIVALENT" if metrics pass thresholds}
#'   }
#'
#' @details
#' The validation assesses scientific equivalence rather than numerical identity.
#' Two results are considered equivalent if:
#' \itemize{
#'   \item Pearson correlation > 0.95
#'   \item Dominant cell type match > 0.90
#'   \item RMSE < 0.05
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' validation <- validate_strider(strider_result, "stride_result.txt")
#' print(validation)
#' }
validate_strider <- function(strider_result, reference_result,
                              metrics = c("correlation", "dominant", "error", "rank"),
                              tolerance = 0.01) {
  
  vmessage("=== Validating STRIDER Results ===")
  
  # Extract deconvolution matrix from STRIDER object
  if (inherits(strider_result, "STRIDER")) {
    strider_mat <- as.matrix(strider_result$deconv)
  } else {
    strider_mat <- as.matrix(strider_result)
  }
  
  # Load reference result
  if (is.character(reference_result) && length(reference_result) == 1) {
    if (file.exists(reference_result)) {
      ref_mat <- as.matrix(data.table::fread(reference_result, data.table = FALSE, header = TRUE, row.names = 1))
    } else {
      rlang::abort(paste0("Reference file not found: ", reference_result))
    }
  } else if (is.data.frame(reference_result)) {
    ref_mat <- as.matrix(reference_result)
  } else {
    ref_mat <- as.matrix(reference_result)
  }
  
  # Align spots and cell types
  common_spots <- intersect(rownames(strider_mat), rownames(ref_mat))
  common_ct <- intersect(colnames(strider_mat), colnames(ref_mat))
  
  if (length(common_spots) == 0) {
    rlang::abort("No common spots found between results")
  }
  if (length(common_ct) == 0) {
    rlang::abort("No common cell types found between results")
  }
  
  vmessage("  Common spots: ", length(common_spots), " / ", 
           nrow(strider_mat), " (STRIDER), ", nrow(ref_mat), " (reference)")
  vmessage("  Common cell types: ", length(common_ct))
  
  strider_mat <- strider_mat[common_spots, common_ct]
  ref_mat <- ref_mat[common_spots, common_ct]
  
  # Initialize result
  result <- list()
  
  # ===== Correlation metrics =====
  if ("correlation" %in% metrics) {
    strider_vec <- as.vector(strider_mat)
    ref_vec <- as.vector(ref_mat)
    
    result$pearson_r <- stats::cor(strider_vec, ref_vec, method = "pearson")
    result$spearman_r <- stats::cor(strider_vec, ref_vec, method = "spearman")
    
    vmessage("  Pearson r: ", round(result$pearson_r, 4))
    vmessage("  Spearman r: ", round(result$spearman_r, 4))
  }
  
  # ===== Error metrics =====
  if ("error" %in% metrics) {
    diff <- strider_mat - ref_mat
    result$rmse <- sqrt(mean(diff^2))
    result$mae <- mean(abs(diff))
    result$max_error <- max(abs(diff))
    
    vmessage("  RMSE: ", round(result$rmse, 4))
    vmessage("  MAE: ", round(result$mae, 4))
  }
  
  # ===== Dominant cell type match =====
  if ("dominant" %in% metrics) {
    strider_dominant <- colnames(strider_mat)[apply(strider_mat, 1, which.max)]
    ref_dominant <- colnames(ref_mat)[apply(ref_mat, 1, which.max)]
    
    result$dominant_match <- mean(strider_dominant == ref_dominant)
    
    vmessage("  Dominant match: ", round(result$dominant_match * 100, 1), "%")
  }
  
  # ===== Per-spot rank correlation =====
  if ("rank" %in% metrics) {
    rank_cors <- sapply(seq_len(nrow(strider_mat)), function(i) {
      stats::cor(strider_mat[i, ], ref_mat[i, ], method = "spearman")
    })
    result$rank_correlation <- mean(rank_cors, na.rm = TRUE)
    result$rank_correlation_sd <- stats::sd(rank_cors, na.rm = TRUE)
    
    vmessage("  Mean rank correlation: ", round(result$rank_correlation, 4))
  }
  
  # ===== Verdict =====
  result$verdict <- "EQUIVALENT"
  
  if (!is.null(result$pearson_r) && result$pearson_r < 0.95) {
    result$verdict <- "NEEDS REVIEW"
  }
  if (!is.null(result$dominant_match) && result$dominant_match < 0.90) {
    result$verdict <- "NEEDS REVIEW"
  }
  if (!is.null(result$rmse) && result$rmse > 0.05) {
    result$verdict <- "NEEDS REVIEW"
  }
  
  vmessage("\n  Verdict: ", result$verdict)
  
  # Store metadata
  result$n_spots <- length(common_spots)
  result$n_celltypes <- length(common_ct)
  result$celltypes <- common_ct
  
  class(result) <- "strider_validation"
  result
}

#' Print Validation Report
#' @param x A strider_validation object
#' @param ... Additional arguments (unused)
#' @export
print.strider_validation <- function(x, ...) {
  cat("\n")
  cat("STRIDER Validation Report\n")
  cat("=========================\n\n")
  
  cat("Data:\n")
  cat("  Spots:", x$n_spots, "\n")
  cat("  Cell types:", x$n_celltypes, "\n\n")
  
  cat("Metrics:\n")
  if (!is.null(x$pearson_r)) {
    cat("  Pearson correlation:", sprintf("%.4f", x$pearson_r), "\n")
    cat("  Spearman correlation:", sprintf("%.4f", x$spearman_r), "\n")
  }
  if (!is.null(x$rmse)) {
    cat("  RMSE:", sprintf("%.4f", x$rmse), "\n")
    cat("  MAE:", sprintf("%.4f", x$mae), "\n")
  }
  if (!is.null(x$dominant_match)) {
    cat("  Dominant cell type match:", sprintf("%.1f%%", x$dominant_match * 100), "\n")
  }
  if (!is.null(x$rank_correlation)) {
    cat("  Mean rank correlation:", sprintf("%.4f", x$rank_correlation), "\n")
  }
  
  cat("\n")
  if (x$verdict == "EQUIVALENT") {
    cat("Verdict: ", x$verdict, " ✓\n")
    cat("Results are scientifically equivalent to the reference.\n")
  } else {
    cat("Verdict: ", x$verdict, "\n")
    cat("Results may differ from the reference. Review recommended.\n")
  }
  
  invisible(x)
}

#' Compare Cell Type Proportions
#'
#' @description Visualize comparison between STRIDER and reference results.
#'
#' @param strider_result STRIDER result
#' @param reference_result Reference result
#' @param plot_type Type of plot: "scatter", "bland_altman", "correlation"
#'
#' @return A ggplot object
#' @export
plot_validation <- function(strider_result, reference_result, 
                            plot_type = "scatter") {
  
  # Extract matrices
  if (inherits(strider_result, "STRIDER")) {
    strider_mat <- as.matrix(strider_result$deconv)
  } else {
    strider_mat <- as.matrix(strider_result)
  }
  
  if (is.character(reference_result)) {
    ref_mat <- as.matrix(data.table::fread(reference_result, data.table = FALSE, 
                                            header = TRUE, row.names = 1))
  } else {
    ref_mat <- as.matrix(reference_result)
  }
  
  # Align
  common_spots <- intersect(rownames(strider_mat), rownames(ref_mat))
  common_ct <- intersect(colnames(strider_mat), colnames(ref_mat))
  
  strider_mat <- strider_mat[common_spots, common_ct]
  ref_mat <- ref_mat[common_spots, common_ct]
  
  # Create data frame for plotting
  plot_df <- data.frame(
    STRIDER = as.vector(strider_mat),
    Reference = as.vector(ref_mat),
    CellType = rep(colnames(strider_mat), each = nrow(strider_mat)),
    Spot = rep(rownames(strider_mat), ncol(strider_mat))
  )
  
  if (plot_type == "scatter") {
    # Scatter plot
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Reference, y = STRIDER)) +
      ggplot2::geom_point(alpha = 0.3, size = 0.5) +
      ggplot2::geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      ggplot2::facet_wrap(~ CellType, scales = "free") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "STRIDER vs Reference",
        x = "Reference Fraction",
        y = "STRIDER Fraction"
      )
    
  } else if (plot_type == "bland_altman") {
    # Bland-Altman plot
    plot_df$Mean <- (plot_df$STRIDER + plot_df$Reference) / 2
    plot_df$Diff <- plot_df$STRIDER - plot_df$Reference
    
    mean_diff <- mean(plot_df$Diff)
    sd_diff <- stats::sd(plot_df$Diff)
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Mean, y = Diff)) +
      ggplot2::geom_point(alpha = 0.3, size = 0.5) +
      ggplot2::geom_hline(yintercept = mean_diff, color = "blue") +
      ggplot2::geom_hline(yintercept = mean_diff + 1.96 * sd_diff, 
                          color = "red", linetype = "dashed") +
      ggplot2::geom_hline(yintercept = mean_diff - 1.96 * sd_diff, 
                          color = "red", linetype = "dashed") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Bland-Altman Plot",
        x = "Mean Fraction",
        y = "Difference (STRIDER - Reference)"
      )
    
  } else if (plot_type == "correlation") {
    # Per-cell-type correlation
    cor_df <- data.frame(
      CellType = colnames(strider_mat),
      Correlation = sapply(seq_len(ncol(strider_mat)), function(i) {
        stats::cor(strider_mat[, i], ref_mat[, i])
      })
    )
    cor_df <- cor_df[order(-cor_df$Correlation), ]
    cor_df$CellType <- factor(cor_df$CellType, levels = cor_df$CellType)
    
    p <- ggplot2::ggplot(cor_df, ggplot2::aes(x = CellType, y = Correlation)) +
      ggplot2::geom_bar(stat = "identity", fill = "#3C5488FF") +
      ggplot2::geom_hline(yintercept = 0.95, color = "red", linetype = "dashed") +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Per-Cell-Type Correlation",
        x = NULL,
        y = "Pearson Correlation"
      )
  }
  
  p
}

#' Benchmark STRIDER Against Known Ground Truth
#'
#' @description Evaluate STRIDER performance using simulated data with known
#'   cell type proportions.
#'
#' @param strider_result STRIDER result
#' @param true_proportions True cell type proportions (spots x celltypes)
#'
#' @return Benchmark metrics
#' @export
benchmark_strider <- function(strider_result, true_proportions) {
  
  vmessage("Benchmarking against ground truth...")
  
  if (inherits(strider_result, "STRIDER")) {
    pred <- as.matrix(strider_result$deconv)
  } else {
    pred <- as.matrix(strider_result)
  }
  true <- as.matrix(true_proportions)
  
  # Align
  common_spots <- intersect(rownames(pred), rownames(true))
  common_ct <- intersect(colnames(pred), colnames(true))
  
  pred <- pred[common_spots, common_ct]
  true <- true[common_spots, common_ct]
  
  # Calculate metrics
  results <- list()
  
  # Overall correlation
  results$pearson_overall <- stats::cor(as.vector(pred), as.vector(true))
  results$spearman_overall <- stats::cor(as.vector(pred), as.vector(true), method = "spearman")
  
  # Per-cell-type correlation
  results$pearson_per_ct <- sapply(seq_len(ncol(pred)), function(i) {
    stats::cor(pred[, i], true[, i])
  })
  names(results$pearson_per_ct) <- common_ct
  
  # RMSE
  results$rmse_overall <- sqrt(mean((pred - true)^2))
  results$rmse_per_ct <- sqrt(colMeans((pred - true)^2))
  
  # Dominant cell type accuracy
  pred_dominant <- colnames(pred)[apply(pred, 1, which.max)]
  true_dominant <- colnames(true)[apply(true, 1, which.max)]
  results$dominant_accuracy <- mean(pred_dominant == true_dominant)
  
  # Jensen-Shannon divergence (per spot)
  results$jsd_per_spot <- sapply(seq_len(nrow(pred)), function(i) {
    p <- pred[i, ]
    q <- true[i, ]
    m <- (p + q) / 2
    jsd <- 0.5 * sum(p * log2(p / m + 1e-10)) + 0.5 * sum(q * log2(q / m + 1e-10))
    jsd
  })
  results$jsd_mean <- mean(results$jsd_per_spot)
  
  vmessage("  Overall Pearson r: ", round(results$pearson_overall, 4))
  vmessage("  Overall RMSE: ", round(results$rmse_overall, 4))
  vmessage("  Dominant accuracy: ", round(results$dominant_accuracy * 100, 1), "%")
  vmessage("  Mean JSD: ", round(results$jsd_mean, 4))
  
  class(results) <- "strider_benchmark"
  results
}

#' Print Benchmark Results
#' @param x A strider_benchmark object
#' @param ... Additional arguments (unused)
#' @export
print.strider_benchmark <- function(x, ...) {
  cat("\n")
  cat("STRIDER Benchmark Results\n")
  cat("=========================\n\n")
  
  cat("Overall Metrics:\n")
  cat("  Pearson correlation:", sprintf("%.4f", x$pearson_overall), "\n")
  cat("  Spearman correlation:", sprintf("%.4f", x$spearman_overall), "\n")
  cat("  RMSE:", sprintf("%.4f", x$rmse_overall), "\n")
  cat("  Dominant cell type accuracy:", sprintf("%.1f%%", x$dominant_accuracy * 100), "\n")
  cat("  Mean Jensen-Shannon divergence:", sprintf("%.4f", x$jsd_mean), "\n\n")
  
  cat("Per-Cell-Type Pearson r:\n")
  for (ct in names(x$pearson_per_ct)) {
    cat("  ", ct, ":", sprintf("%.4f", x$pearson_per_ct[ct]), "\n")
  }
  
  invisible(x)
}
