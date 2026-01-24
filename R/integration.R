#' @title Multi-Sample Integration Functions for STRIDER
#' @description Functions for integrating multiple spatial transcriptomics samples
#'   using Fused Gromov-Wasserstein optimal transport
#' @name integration
NULL

#' Integrate Multiple Spatial Transcriptomics Samples
#'
#' @description Integrate multiple spatial transcriptomics samples using 
#'   Fused Gromov-Wasserstein (FGW) optimal transport. This aligns samples
#'   based on both topic distribution similarity and spatial structure.
#'
#' @param deconv.list List of deconvolution result matrices (spots x celltypes)
#' @param topic.list List of topic distribution matrices (topics x spots)
#' @param coords.list List of spatial coordinate data frames
#' @param sample.ids Sample identifiers (default: S1, S2, ...)
#' @param alpha Trade-off parameter between feature and spatial similarity (default: 0.5)
#'   - 0: pure OT (only topic similarity)
#'   - 1: pure GW (only spatial structure)
#' @param reg Sinkhorn regularization parameter (default: 0.1)
#' @param max.iter Maximum iterations for FGW (default: 100)
#' @param seed Random seed
#'
#' @return A list containing:
#'   \itemize{
#'     \item{aligned_coords}{Data frame of aligned coordinates for all samples}
#'     \item{transport_plans}{List of optimal transport plans between consecutive samples}
#'     \item{deconv}{Combined deconvolution results}
#'     \item{rotation_matrices}{List of rotation matrices}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' integrated <- integrate_samples(
#'   deconv.list = list(deconv1, deconv2, deconv3),
#'   topic.list = list(topic1, topic2, topic3),
#'   coords.list = list(coords1, coords2, coords3)
#' )
#' }
integrate_samples <- function(deconv.list, topic.list, coords.list,
                               sample.ids = NULL, alpha = 0.5, reg = 0.1,
                               max.iter = 100, seed = 42) {
  
  vmessage("=== Integrating Multiple Samples ===")
  set.seed(seed)
  
  n_samples <- length(deconv.list)
  
  if (n_samples < 2) {
    rlang::abort("At least 2 samples are required for integration")
  }
  
  # Generate sample IDs if not provided
  if (is.null(sample.ids)) {
    sample.ids <- paste0("S", seq_len(n_samples))
  }
  
  vmessage("Number of samples: ", n_samples)
  vmessage("Sample IDs: ", paste(sample.ids, collapse = ", "))
  vmessage("Alpha (structure vs feature): ", alpha)
  
  # ===== Process first sample =====
  vmessage("\nProcessing sample 1: ", sample.ids[1])
  
  # Get data for first sample
  s1_deconv <- deconv.list[[1]]
  s1_topic <- topic.list[[1]]
  s1_coords <- coords.list[[1]]
  
  # Center first sample
  s1_weights <- rep(1/nrow(s1_coords), nrow(s1_coords))
  s1_centered <- center_coords(s1_coords, s1_weights)
  
  # Initialize results
  all_coords <- data.frame(
    spot = rownames(s1_coords),
    x_original = s1_coords[, 1],
    y_original = s1_coords[, 2],
    x_aligned = s1_centered[, 1],
    y_aligned = s1_centered[, 2],
    sample = sample.ids[1],
    stringsAsFactors = FALSE
  )
  rownames(all_coords) <- paste0(sample.ids[1], "_", rownames(s1_coords))
  
  transport_plans <- list()
  rotation_matrices <- list()
  
  # Cumulative rotation matrix
  R_cumulative <- diag(2)
  
  # ===== Align subsequent samples =====
  for (i in 2:n_samples) {
    vmessage("\nAligning sample ", i, ": ", sample.ids[i], " to ", sample.ids[i-1])
    
    # Get data for current sample
    s2_deconv <- deconv.list[[i]]
    s2_topic <- topic.list[[i]]
    s2_coords <- coords.list[[i]]
    
    # Get data for previous sample
    s1_deconv <- deconv.list[[i-1]]
    s1_topic <- topic.list[[i-1]]
    s1_coords <- coords.list[[i-1]]
    
    # Align samples
    alignment <- align_samples(
      topic1 = s1_topic,
      topic2 = s2_topic,
      coords1 = s1_coords,
      coords2 = s2_coords,
      alpha = alpha,
      reg = reg,
      max.iter = max.iter
    )
    
    # Update cumulative rotation
    R_cumulative <- R_cumulative %*% alignment$rotation
    
    # Apply cumulative rotation to center current sample
    s2_weights <- rep(1/nrow(s2_coords), nrow(s2_coords))
    s2_centered <- center_coords(s2_coords, s2_weights)
    s2_aligned <- as.matrix(s2_centered) %*% t(R_cumulative)
    
    # Add to results
    sample_coords <- data.frame(
      spot = rownames(s2_coords),
      x_original = s2_coords[, 1],
      y_original = s2_coords[, 2],
      x_aligned = s2_aligned[, 1],
      y_aligned = s2_aligned[, 2],
      sample = sample.ids[i],
      stringsAsFactors = FALSE
    )
    rownames(sample_coords) <- paste0(sample.ids[i], "_", rownames(s2_coords))
    
    all_coords <- rbind(all_coords, sample_coords)
    
    # Store transport plan and rotation
    transport_plans[[paste0(sample.ids[i-1], "_", sample.ids[i])]] <- alignment$transport
    rotation_matrices[[sample.ids[i]]] <- alignment$rotation
    
    vmessage("  FGW distance: ", round(alignment$fgw_dist, 4))
  }
  
  # ===== Combine deconvolution results =====
  vmessage("\nCombining deconvolution results...")
  
  combined_deconv <- do.call(rbind, lapply(seq_along(deconv.list), function(i) {
    df <- as.data.frame(deconv.list[[i]])
    rownames(df) <- paste0(sample.ids[i], "_", rownames(deconv.list[[i]]))
    df
  }))
  
  vmessage("\n=== Integration Complete ===")
  vmessage("Total spots: ", nrow(all_coords))
  
  list(
    aligned_coords = all_coords,
    transport_plans = transport_plans,
    rotation_matrices = rotation_matrices,
    deconv = combined_deconv,
    sample_ids = sample.ids
  )
}

#' Align Two Samples Using FGW
#'
#' @description Align two spatial transcriptomics samples using Fused 
#'   Gromov-Wasserstein optimal transport.
#'
#' @param topic1 Topic distribution for sample 1 (topics x spots)
#' @param topic2 Topic distribution for sample 2 (topics x spots)
#' @param coords1 Spatial coordinates for sample 1
#' @param coords2 Spatial coordinates for sample 2
#' @param alpha Trade-off parameter (default: 0.5)
#' @param reg Regularization parameter (default: 0.1)
#' @param max.iter Maximum iterations (default: 100)
#'
#' @return A list containing transport plan, rotation matrix, and FGW distance
#' @export
#'
#' @examples
#' \dontrun{
#' alignment <- align_samples(topic1, topic2, coords1, coords2)
#' }
align_samples <- function(topic1, topic2, coords1, coords2,
                          alpha = 0.5, reg = 0.1, max.iter = 100) {
  
  # Transpose topics to spots x topics
  topic1_t <- t(topic1)
  topic2_t <- t(topic2)
  
  n1 <- nrow(topic1_t)
  n2 <- nrow(topic2_t)
  
  # Uniform weights
  p <- rep(1/n1, n1)
  q <- rep(1/n2, n2)
  
  # Compute KL divergence matrix for topic similarity
  M <- compute_kl_divergence(topic1_t, topic2_t)
  
  # Compute spatial distance matrices
  C1 <- compute_distance_matrix(coords1)
  C2 <- compute_distance_matrix(coords2)
  
  # Run FGW
  fgw_result <- fused_gromov_wasserstein_cpp(
    M = M,
    C1 = C1,
    C2 = C2,
    p = p,
    q = q,
    alpha = alpha,
    reg = reg,
    max_iter = max.iter,
    tol = 1e-7
  )
  
  # Compute optimal rotation
  rotation <- compute_rotation(coords1, coords2, fgw_result$T, p, q)
  
  list(
    transport = fgw_result$T,
    rotation = rotation,
    fgw_dist = fgw_result$fgw_dist
  )
}

#' Center Coordinates
#'
#' @description Center coordinates by weighted mean.
#'
#' @param coords Coordinate data frame or matrix
#' @param weights Point weights
#'
#' @return Centered coordinates
#' @keywords internal
center_coords <- function(coords, weights = NULL) {
  coord_mat <- as.matrix(coords[, 1:2])
  
  if (is.null(weights)) {
    weights <- rep(1/nrow(coord_mat), nrow(coord_mat))
  }
  
  # Use C++ implementation
  centered <- center_points_cpp(coord_mat, weights)
  
  colnames(centered) <- c("x", "y")
  rownames(centered) <- rownames(coords)
  
  as.data.frame(centered)
}

#' Compute KL Divergence Matrix
#'
#' @description Compute pairwise KL divergence between topic distributions.
#'
#' @param P Topic distribution matrix 1 (spots x topics)
#' @param Q Topic distribution matrix 2 (spots x topics)
#'
#' @return KL divergence matrix (n1 x n2)
#' @keywords internal
compute_kl_divergence <- function(P, Q) {
  # Use C++ implementation
  kl_divergence_matrix_cpp(as.matrix(P), as.matrix(Q))
}

#' Compute Spatial Distance Matrix
#'
#' @description Compute pairwise Euclidean distance matrix.
#'
#' @param coords Coordinate data frame or matrix
#'
#' @return Distance matrix
#' @keywords internal
compute_distance_matrix <- function(coords) {
  coord_mat <- as.matrix(coords[, 1:2])
  euclidean_distance_cpp(coord_mat, NULL)
}

#' Compute Rotation Matrix
#'
#' @description Compute optimal rotation matrix using SVD.
#'
#' @param coords1 Coordinates of sample 1
#' @param coords2 Coordinates of sample 2
#' @param transport Optimal transport plan
#' @param p Source weights
#' @param q Target weights
#'
#' @return Rotation matrix (2 x 2)
#' @keywords internal
compute_rotation <- function(coords1, coords2, transport, p, q) {
  
  # Center coordinates
  c1 <- as.matrix(coords1[, 1:2])
  c2 <- as.matrix(coords2[, 1:2])
  
  c1_centered <- center_points_cpp(c1, p)
  c2_centered <- center_points_cpp(c2, q)
  
  # Compute optimal rotation
  # Object matrix: X_2^T * Pi^T * X_1
  obj_mat <- t(c2_centered) %*% t(transport) %*% c1_centered
  
  # SVD
  svd_result <- svd(obj_mat)
  
  # Rotation matrix
  R <- svd_result$u %*% t(svd_result$v)
  
  # Ensure proper rotation (det = 1)
  if (det(R) < 0) {
    svd_result$v[, 2] <- -svd_result$v[, 2]
    R <- svd_result$u %*% t(svd_result$v)
  }
  
  R
}

#' Plot Integration Results
#'
#' @description Visualize integrated samples in 3D or 2D.
#'
#' @param integration_result Result from integrate_samples()
#' @param deconv Deconvolution results for pie chart coloring
#' @param plot.3d Create 3D plot (default: FALSE)
#' @param pt.size Point size (default: 2)
#' @param sample.spacing Spacing between samples in 3D (default: 1)
#'
#' @return A ggplot object
#' @export
plot_integration <- function(integration_result, deconv = NULL, 
                             plot.3d = FALSE, pt.size = 2,
                             sample.spacing = 1) {
  
  coords <- integration_result$aligned_coords
  sample_ids <- unique(coords$sample)
  
  if (plot.3d) {
    rlang::warn("3D plotting requires additional setup. Returning 2D faceted plot.")
    plot.3d <- FALSE
  }
  
  # Create 2D faceted plot
  if (!is.null(deconv)) {
    # Color by dominant cell type
    aligned_spots <- rownames(coords)
    common_spots <- intersect(aligned_spots, rownames(deconv))
    
    deconv_subset <- deconv[common_spots, , drop = FALSE]
    dominant_ct <- colnames(deconv_subset)[apply(deconv_subset, 1, which.max)]
    coords$celltype <- dominant_ct[match(rownames(coords), common_spots)]
    
    colors <- get_extended_colors(length(unique(dominant_ct)))
    
    p <- ggplot2::ggplot(coords, ggplot2::aes(x = x_aligned, y = y_aligned, color = celltype)) +
      ggplot2::geom_point(size = pt.size) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::facet_wrap(~ sample, scales = "free") +
      ggplot2::coord_equal() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Integrated Samples",
        x = "Aligned X",
        y = "Aligned Y",
        color = "Cell Type"
      )
  } else {
    # Color by sample
    colors <- get_extended_colors(length(sample_ids))
    
    p <- ggplot2::ggplot(coords, ggplot2::aes(x = x_aligned, y = y_aligned, color = sample)) +
      ggplot2::geom_point(size = pt.size) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::coord_equal() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Integrated Samples",
        x = "Aligned X",
        y = "Aligned Y",
        color = "Sample"
      )
  }
  
  p
}

#' Calculate Integration Quality Metrics
#'
#' @description Calculate metrics to assess integration quality.
#'
#' @param integration_result Result from integrate_samples()
#' @param topic.list List of topic distributions
#'
#' @return Data frame with quality metrics
#' @export
integration_quality <- function(integration_result, topic.list) {
  
  transport_plans <- integration_result$transport_plans
  
  metrics <- lapply(names(transport_plans), function(name) {
    T_mat <- transport_plans[[name]]
    
    # Transport entropy (higher = more diffuse mapping)
    T_flat <- as.vector(T_mat)
    T_flat <- T_flat[T_flat > 1e-10]
    entropy <- -sum(T_flat * log(T_flat))
    
    # Sparsity (fraction of near-zero entries)
    sparsity <- mean(T_mat < 1e-6)
    
    data.frame(
      pair = name,
      entropy = entropy,
      sparsity = sparsity,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, metrics)
}
