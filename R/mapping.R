#' @title Spot-to-Cell Mapping Functions for STRIDER
#' @description Functions for mapping spatial spots to similar single cells
#' @name mapping
NULL

#' Map Spatial Spots to Similar Single Cells
#'
#' @description Identify the most similar single cells for each spatial spot
#'   based on topic distribution similarity.
#'
#' @param topic.spot Topic distribution for spots (topics x spots)
#' @param topic.cell Topic distribution for cells (topics x cells)
#' @param n.top Number of top similar cells per spot (default: 10)
#' @param method Distance metric: "euclidean", "cosine", or "correlation" (default: "euclidean")
#'
#' @return A data frame with spot-cell mappings and similarity scores
#' @export
#'
#' @examples
#' \dontrun{
#' mapping <- map_spots_to_cells(topic_spot, topic_cell, n.top = 10)
#' }
map_spots_to_cells <- function(topic.spot, topic.cell, n.top = 10,
                                method = "euclidean") {
  
  vmessage("Mapping spots to similar cells...")
  
  # Transpose to spots/cells x topics
  spot_topic <- t(topic.spot)
  cell_topic <- t(topic.cell)
  
  n_spots <- nrow(spot_topic)
  n_cells <- nrow(cell_topic)
  
  vmessage("  ", n_spots, " spots, ", n_cells, " cells")
  vmessage("  Finding top ", n.top, " similar cells per spot")
  
  # Compute distance matrix
  if (method == "euclidean") {
    dist_mat <- euclidean_distance_cpp(as.matrix(spot_topic), as.matrix(cell_topic))
  } else if (method == "cosine") {
    dist_mat <- compute_cosine_distance(spot_topic, cell_topic)
  } else if (method == "correlation") {
    dist_mat <- compute_correlation_distance(spot_topic, cell_topic)
  } else {
    rlang::abort(paste0("Unknown distance method: ", method))
  }
  
  # Find top matches for each spot
  result_list <- lapply(seq_len(n_spots), function(i) {
    distances <- dist_mat[i, ]
    
    # Get top n.top closest cells
    top_idx <- order(distances)[seq_len(min(n.top, n_cells))]
    
    data.frame(
      spot = rownames(spot_topic)[i],
      rank = seq_along(top_idx),
      cell = colnames(topic.cell)[top_idx],
      distance = distances[top_idx],
      stringsAsFactors = FALSE
    )
  })
  
  result <- do.call(rbind, result_list)
  
  vmessage("  Mapping complete")
  
  result
}

#' Get Mapped Cell Information
#'
#' @description Extract information about cells mapped to each spot.
#'
#' @param mapping Mapping result from map_spots_to_cells()
#' @param celltypes Named vector of cell type annotations
#' @param aggregate How to aggregate cell types: "majority", "all", or "weighted"
#'
#' @return Data frame with spot-level cell type information
#' @export
get_mapped_celltypes <- function(mapping, celltypes, aggregate = "majority") {
  
  vmessage("Extracting mapped cell type information...")
  
  # Add cell types to mapping
  mapping$celltype <- celltypes[mapping$cell]
  
  # Aggregate by spot
  spots <- unique(mapping$spot)
  
  result <- lapply(spots, function(s) {
    spot_data <- mapping[mapping$spot == s, ]
    
    if (aggregate == "majority") {
      # Most common cell type
      ct_table <- table(spot_data$celltype)
      majority_ct <- names(ct_table)[which.max(ct_table)]
      
      data.frame(
        spot = s,
        predicted_celltype = majority_ct,
        n_cells = nrow(spot_data),
        stringsAsFactors = FALSE
      )
      
    } else if (aggregate == "all") {
      # All mapped cell types with counts
      ct_table <- table(spot_data$celltype)
      
      data.frame(
        spot = s,
        celltype = names(ct_table),
        count = as.numeric(ct_table),
        stringsAsFactors = FALSE
      )
      
    } else if (aggregate == "weighted") {
      # Weighted by inverse distance
      spot_data$weight <- 1 / (spot_data$distance + 1e-10)
      spot_data$weight <- spot_data$weight / sum(spot_data$weight)
      
      ct_weights <- tapply(spot_data$weight, spot_data$celltype, sum)
      majority_ct <- names(ct_weights)[which.max(ct_weights)]
      
      data.frame(
        spot = s,
        predicted_celltype = majority_ct,
        confidence = max(ct_weights),
        stringsAsFactors = FALSE
      )
    }
  })
  
  do.call(rbind, result)
}

#' Compute Cosine Distance Matrix
#'
#' @keywords internal
compute_cosine_distance <- function(X, Y) {
  # Normalize rows to unit length
  X_norm <- X / sqrt(rowSums(X^2) + 1e-10)
  Y_norm <- Y / sqrt(rowSums(Y^2) + 1e-10)
  
  # Cosine similarity
  sim <- as.matrix(X_norm) %*% t(as.matrix(Y_norm))
  
  # Convert to distance
  1 - sim
}

#' Compute Correlation Distance Matrix
#'
#' @keywords internal
compute_correlation_distance <- function(X, Y) {
  # Center rows
  X_centered <- X - rowMeans(X)
  Y_centered <- Y - rowMeans(Y)
  
  # Correlation
  X_sd <- sqrt(rowSums(X_centered^2) + 1e-10)
  Y_sd <- sqrt(rowSums(Y_centered^2) + 1e-10)
  
  cor_mat <- (as.matrix(X_centered) %*% t(as.matrix(Y_centered))) / 
             (X_sd %*% t(Y_sd))
  
  # Convert to distance
  1 - cor_mat
}

#' Transfer Labels from Single Cells to Spots
#'
#' @description Transfer cell type labels from single-cell reference to spatial spots
#'   using mapping results.
#'
#' @param mapping Mapping result from map_spots_to_cells()
#' @param celltypes Named vector of cell type annotations
#' @param method Transfer method: "knn" (k-nearest neighbor voting) or "weighted"
#' @param k Number of neighbors for kNN (default: all mapped cells)
#'
#' @return Named vector of transferred labels
#' @export
transfer_labels <- function(mapping, celltypes, method = "knn", k = NULL) {
  
  vmessage("Transferring labels from cells to spots...")
  
  spots <- unique(mapping$spot)
  
  if (is.null(k)) {
    k <- max(mapping$rank)
  }
  
  labels <- sapply(spots, function(s) {
    spot_data <- mapping[mapping$spot == s & mapping$rank <= k, ]
    spot_data$celltype <- celltypes[spot_data$cell]
    
    if (method == "knn") {
      # Simple majority vote
      ct_table <- table(spot_data$celltype)
      names(ct_table)[which.max(ct_table)]
      
    } else if (method == "weighted") {
      # Distance-weighted vote
      spot_data$weight <- 1 / (spot_data$distance + 1e-10)
      ct_weights <- tapply(spot_data$weight, spot_data$celltype, sum)
      names(ct_weights)[which.max(ct_weights)]
    }
  })
  
  names(labels) <- spots
  labels
}

#' Create Spot-Cell Neighborhood Graph
#'
#' @description Create a bipartite graph connecting spots to their most similar cells.
#'
#' @param mapping Mapping result from map_spots_to_cells()
#' @param threshold Distance threshold for edges (optional)
#'
#' @return An igraph object (if igraph is available) or edge list
#' @export
create_mapping_graph <- function(mapping, threshold = NULL) {
  
  # Filter by threshold if specified
  if (!is.null(threshold)) {
    mapping <- mapping[mapping$distance <= threshold, ]
  }
  
  # Create edge list
  edges <- data.frame(
    from = paste0("spot_", mapping$spot),
    to = paste0("cell_", mapping$cell),
    weight = 1 / (mapping$distance + 1e-10)
  )
  
  # Try to create igraph
  if (requireNamespace("igraph", quietly = TRUE)) {
    g <- igraph::graph_from_data_frame(edges, directed = FALSE)
    igraph::V(g)$type <- grepl("^spot_", igraph::V(g)$name)
    return(g)
  } else {
    rlang::warn("Package 'igraph' not available. Returning edge list.")
    return(edges)
  }
}

#' Visualize Spot-Cell Mapping
#'
#' @description Create visualization of spot-cell mapping relationships.
#'
#' @param mapping Mapping result from map_spots_to_cells()
#' @param celltypes Named vector of cell type annotations
#' @param spots Specific spots to visualize (default: first 10)
#'
#' @return A ggplot object
#' @export
plot_mapping <- function(mapping, celltypes, spots = NULL) {
  
  if (is.null(spots)) {
    spots <- unique(mapping$spot)[1:min(10, length(unique(mapping$spot)))]
  }
  
  # Subset mapping
  mapping_sub <- mapping[mapping$spot %in% spots, ]
  mapping_sub$celltype <- celltypes[mapping_sub$cell]
  
  # Create plot
  colors <- get_extended_colors(length(unique(mapping_sub$celltype)))
  
  ggplot2::ggplot(mapping_sub, ggplot2::aes(x = spot, y = rank, fill = celltype)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title = "Spot-Cell Mapping",
      x = "Spot",
      y = "Rank",
      fill = "Cell Type"
    )
}
