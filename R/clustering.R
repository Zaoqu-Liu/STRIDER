#' @title Spatial Clustering Functions for STRIDER
#' @description Functions for neighborhood-aware clustering of spatial spots
#' @name clustering
NULL

#' Find Spatial Neighbors
#'
#' @description Find neighbors for each spot based on spatial coordinates
#'   using k-nearest neighbors or distance threshold.
#'
#' @param coords Data frame with spatial coordinates (columns: x, y)
#' @param k Number of nearest neighbors (default: 6)
#' @param radius Distance radius for neighbors (alternative to k)
#' @param include.self Include the spot itself as a neighbor (default: FALSE)
#'
#' @return A list where each element contains indices of neighbors for each spot
#'
#' @details
#' For Visium data, k=6 corresponds to the hexagonal grid structure.
#' For other technologies, adjust k based on expected spatial relationships.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' neighbors <- find_neighbors(spatial_coords, k = 6)
#' }
find_neighbors <- function(coords, k = 6, radius = NULL, include.self = FALSE) {
  
  vmessage("Finding spatial neighbors...")
  
  # Extract coordinate columns
  if ("x" %in% colnames(coords) && "y" %in% colnames(coords)) {
    coord_matrix <- as.matrix(coords[, c("x", "y")])
  } else if (ncol(coords) >= 2) {
    coord_matrix <- as.matrix(coords[, 1:2])
  } else {
    rlang::abort("coords must have at least 2 columns for x and y coordinates")
  }
  
  n_spots <- nrow(coord_matrix)
  rownames(coord_matrix) <- rownames(coords)
  
  if (!is.null(radius)) {
    # Radius-based neighbors
    dist_mat <- as.matrix(stats::dist(coord_matrix, method = "euclidean"))
    
    neighbor_list <- lapply(seq_len(n_spots), function(i) {
      neighbors <- which(dist_mat[i, ] <= radius)
      if (!include.self) {
        neighbors <- neighbors[neighbors != i]
      }
      neighbors
    })
  } else {
    # K-nearest neighbors
    if (!requireNamespace("FNN", quietly = TRUE)) {
      rlang::abort("Package 'FNN' is required. Install with: install.packages('FNN')")
    }
    
    # FNN includes the point itself, so adjust k
    k_search <- if (include.self) k else k + 1
    
    nn_result <- FNN::get.knn(coord_matrix, k = k_search)
    
    neighbor_list <- lapply(seq_len(n_spots), function(i) {
      neighbors <- nn_result$nn.index[i, ]
      if (!include.self) {
        neighbors <- neighbors[neighbors != i]
      }
      neighbors[seq_len(min(k, length(neighbors)))]
    })
  }
  
  names(neighbor_list) <- rownames(coords)
  
  vmessage("  Spots: ", n_spots)
  vmessage("  Mean neighbors: ", round(mean(lengths(neighbor_list)), 1))
  
  neighbor_list
}

#' Calculate Neighbor Average Fractions
#'
#' @description Calculate the average cell type fractions from neighboring spots.
#'
#' @param deconv Cell type fraction matrix (spots x celltypes)
#' @param neighbors List of neighbor indices from find_neighbors()
#'
#' @return Matrix of neighbor-averaged fractions (spots x celltypes)
#' @export
calculate_neighbor_fractions <- function(deconv, neighbors) {
  
  n_spots <- nrow(deconv)
  n_celltypes <- ncol(deconv)
  
  neighbor_frac <- matrix(0, nrow = n_spots, ncol = n_celltypes)
  rownames(neighbor_frac) <- rownames(deconv)
  colnames(neighbor_frac) <- colnames(deconv)
  
  for (i in seq_len(n_spots)) {
    neighbor_idx <- neighbors[[i]]
    if (length(neighbor_idx) > 0) {
      neighbor_frac[i, ] <- colMeans(deconv[neighbor_idx, , drop = FALSE])
    } else {
      # If no neighbors, use own values
      neighbor_frac[i, ] <- deconv[i, ]
    }
  }
  
  neighbor_frac
}

#' Spatial Clustering Based on Cell Type Composition
#'
#' @description Perform clustering of spatial spots based on cell type
#'   composition with optional neighborhood information.
#'
#' @param deconv Deconvolution result (STRIDER object or matrix)
#' @param coords Spatial coordinates data frame
#' @param n.clusters Number of clusters (default: 5)
#' @param weight Weight for spot's own composition vs neighbors (default: 0.5)
#'   \itemize{
#'     \item 1.0: only own composition
#'     \item 0.0: only neighbor composition
#'     \item 0.5: equal weight (default)
#'   }
#' @param k Number of neighbors for local composition (default: 6)
#' @param method Clustering method: "kmeans" or "hierarchical" (default: "kmeans")
#' @param seed Random seed
#'
#' @return A data frame with spot coordinates and cluster assignments
#'
#' @details
#' The clustering combines each spot's own cell type composition with the
#' average composition of its spatial neighbors. This produces spatially
#' coherent clusters that reflect local tissue microenvironments.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' clusters <- spatial_cluster(strider_result, coords, n.clusters = 5)
#' }
spatial_cluster <- function(deconv, coords, n.clusters = 5, weight = 0.5,
                            k = 6, method = "kmeans", seed = 42) {
  
  vmessage("Performing spatial clustering...")
  
  # Extract deconv matrix
  if (inherits(deconv, "STRIDER")) {
    deconv_mat <- deconv$deconv
  } else {
    deconv_mat <- as.matrix(deconv)
  }
  
  # Ensure coords and deconv are aligned
  common_spots <- intersect(rownames(deconv_mat), rownames(coords))
  if (length(common_spots) == 0) {
    rlang::abort("No matching spots between deconv and coords")
  }
  
  deconv_mat <- deconv_mat[common_spots, , drop = FALSE]
  coords <- coords[common_spots, , drop = FALSE]
  
  # Find neighbors
  neighbors <- find_neighbors(coords, k = k)
  
  # Calculate neighbor fractions
  neighbor_frac <- calculate_neighbor_fractions(deconv_mat, neighbors)
  
  # Combine own and neighbor fractions
  combined_frac <- weight * deconv_mat + (1 - weight) * neighbor_frac
  
  # Clustering
  set.seed(seed)
  
  if (method == "kmeans") {
    km <- stats::kmeans(combined_frac, centers = n.clusters, iter.max = 1000, nstart = 25)
    clusters <- km$cluster
  } else if (method == "hierarchical") {
    dist_mat <- stats::dist(combined_frac)
    hc <- stats::hclust(dist_mat, method = "ward.D2")
    clusters <- stats::cutree(hc, k = n.clusters)
  } else {
    rlang::abort(paste0("Unknown clustering method: ", method))
  }
  
  # Create result data frame
  result <- data.frame(
    spot = rownames(coords),
    x = as.numeric(coords[, 1]),
    y = as.numeric(coords[, 2]),
    cluster = factor(clusters, levels = seq_len(n.clusters)),
    row.names = rownames(coords)
  )
  
  vmessage("  Clustered ", nrow(result), " spots into ", n.clusters, " clusters")
  
  # Report cluster sizes
  cluster_sizes <- table(result$cluster)
  vmessage("  Cluster sizes: ", paste(names(cluster_sizes), "=", cluster_sizes, collapse = ", "))
  
  result
}

#' Calculate Cluster Statistics
#'
#' @description Calculate cell type composition statistics for each cluster.
#'
#' @param deconv Deconvolution matrix (spots x celltypes)
#' @param clusters Cluster assignments (vector or data frame with 'cluster' column)
#'
#' @return Data frame with mean composition per cluster
#' @export
cluster_composition <- function(deconv, clusters) {
  
  # Extract deconv matrix
  if (inherits(deconv, "STRIDER")) {
    deconv_mat <- deconv$deconv
  } else {
    deconv_mat <- as.matrix(deconv)
  }
  
  # Extract cluster vector
  if (is.data.frame(clusters)) {
    cluster_vec <- clusters$cluster
    names(cluster_vec) <- rownames(clusters)
  } else {
    cluster_vec <- clusters
  }
  
  # Align
  common <- intersect(names(cluster_vec), rownames(deconv_mat))
  if (length(common) == 0) {
    rlang::abort("No matching spots between deconv and clusters")
  }
  
  deconv_mat <- deconv_mat[common, , drop = FALSE]
  cluster_vec <- cluster_vec[common]
  
  # Calculate mean composition per cluster
  cluster_levels <- sort(unique(cluster_vec))
  
  composition <- t(sapply(cluster_levels, function(cl) {
    idx <- cluster_vec == cl
    colMeans(deconv_mat[idx, , drop = FALSE])
  }))
  
  rownames(composition) <- paste0("Cluster_", cluster_levels)
  
  as.data.frame(composition)
}

#' Find Enriched Cell Types per Cluster
#'
#' @description Identify cell types that are enriched in each cluster compared
#'   to the overall distribution using statistical testing.
#'
#' @param deconv Deconvolution matrix (spots x celltypes)
#' @param clusters Cluster assignments
#' @param method Statistical test: "wilcox" or "t" (default: "wilcox")
#' @param p.adjust.method P-value adjustment method (default: "BH")
#'
#' @return Data frame with enrichment statistics per cluster-celltype pair
#' @export
find_cluster_markers <- function(deconv, clusters, method = "wilcox",
                                  p.adjust.method = "BH") {
  
  vmessage("Finding enriched cell types per cluster...")
  
  # Extract data
  if (inherits(deconv, "STRIDER")) {
    deconv_mat <- deconv$deconv
  } else {
    deconv_mat <- as.matrix(deconv)
  }
  
  if (is.data.frame(clusters)) {
    cluster_vec <- clusters$cluster
    names(cluster_vec) <- rownames(clusters)
  } else {
    cluster_vec <- clusters
  }
  
  # Align data
  common <- intersect(names(cluster_vec), rownames(deconv_mat))
  deconv_mat <- deconv_mat[common, , drop = FALSE]
  cluster_vec <- cluster_vec[common]
  
  cluster_levels <- sort(unique(cluster_vec))
  celltypes <- colnames(deconv_mat)
  
  results <- list()
  
  for (cl in cluster_levels) {
    cl_results <- lapply(celltypes, function(ct) {
      in_cluster <- deconv_mat[cluster_vec == cl, ct]
      out_cluster <- deconv_mat[cluster_vec != cl, ct]
      
      # Effect size: mean difference
      mean_diff <- mean(in_cluster) - mean(out_cluster)
      
      # Statistical test
      test_result <- tryCatch({
        if (method == "wilcox") {
          stats::wilcox.test(in_cluster, out_cluster, alternative = "greater")
        } else {
          stats::t.test(in_cluster, out_cluster, alternative = "greater")
        }
      }, error = function(e) list(p.value = 1))
      
      data.frame(
        cluster = cl,
        celltype = ct,
        mean_in = mean(in_cluster),
        mean_out = mean(out_cluster),
        diff = mean_diff,
        p_val = test_result$p.value,
        stringsAsFactors = FALSE
      )
    })
    
    results[[as.character(cl)]] <- do.call(rbind, cl_results)
  }
  
  result_df <- do.call(rbind, results)
  result_df$p_adj <- stats::p.adjust(result_df$p_val, method = p.adjust.method)
  
  # Sort by cluster and difference
  result_df <- result_df[order(result_df$cluster, -result_df$diff), ]
  rownames(result_df) <- NULL
  
  result_df
}

#' Calculate Silhouette Score for Clustering
#'
#' @description Calculate silhouette score to assess clustering quality.
#'
#' @param deconv Deconvolution matrix
#' @param clusters Cluster assignments
#'
#' @return Average silhouette score
#' @export
cluster_silhouette <- function(deconv, clusters) {
  
  if (!requireNamespace("cluster", quietly = TRUE)) {
    rlang::abort("Package 'cluster' is required. Install with: install.packages('cluster')")
  }
  
  # Extract data
  if (inherits(deconv, "STRIDER")) {
    deconv_mat <- deconv$deconv
  } else {
    deconv_mat <- as.matrix(deconv)
  }
  
  if (is.data.frame(clusters)) {
    cluster_vec <- as.integer(clusters$cluster)
    names(cluster_vec) <- rownames(clusters)
  } else {
    cluster_vec <- as.integer(clusters)
  }
  
  # Align
  common <- intersect(names(cluster_vec), rownames(deconv_mat))
  deconv_mat <- deconv_mat[common, , drop = FALSE]
  cluster_vec <- cluster_vec[common]
  
  # Calculate silhouette
  dist_mat <- stats::dist(deconv_mat)
  sil <- cluster::silhouette(cluster_vec, dist_mat)
  
  mean(sil[, "sil_width"])
}
