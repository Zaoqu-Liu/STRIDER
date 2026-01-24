#' @title Visualization Functions for STRIDER
#' @description Functions for visualizing deconvolution and analysis results
#' @name plot
NULL

#' Plot Deconvolution Results
#'
#' @description Main plotting function for STRIDER deconvolution results.
#'   Supports scatter pie plots, scatter plots, and heatmaps.
#'
#' @param deconv Deconvolution result (STRIDER object or matrix)
#' @param coords Spatial coordinates data frame
#' @param plot.type Type of plot: "scatterpie", "scatter", or "heatmap"
#' @param sample.id Sample ID to plot (if multiple samples)
#' @param pt.size Point size (default: auto)
#' @param colors Color palette (default: STRIDER palette)
#' @param title Plot title (default: auto)
#' @param coord.flip Flip x and y coordinates (default: FALSE)
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' plot_deconv(strider_result, coords, plot.type = "scatterpie")
#' }
plot_deconv <- function(deconv, coords, plot.type = "scatterpie",
                        sample.id = NULL, pt.size = NULL, colors = NULL,
                        title = NULL, coord.flip = FALSE) {
  
  # Extract deconv matrix
  if (inherits(deconv, "STRIDER")) {
    deconv_mat <- deconv$deconv
  } else {
    deconv_mat <- as.matrix(deconv)
  }
  
  # Handle sample filtering
  if (!is.null(sample.id) && "sample" %in% colnames(coords)) {
    coords <- coords[coords$sample == sample.id, ]
  }
  
  # Align data
  common_spots <- intersect(rownames(deconv_mat), rownames(coords))
  if (length(common_spots) == 0) {
    rlang::abort("No matching spots between deconv and coords")
  }
  
  deconv_mat <- deconv_mat[common_spots, , drop = FALSE]
  coords <- coords[common_spots, , drop = FALSE]
  
  # Set default point size
  if (is.null(pt.size)) {
    pt.size <- if (plot.type == "scatterpie") 8 else 2
  }
  
  # Set default colors
  if (is.null(colors)) {
    colors <- get_extended_colors(ncol(deconv_mat))
  }
  
  # Set title
  if (is.null(title)) {
    title <- if (!is.null(sample.id)) sample.id else "Cell Type Distribution"
  }
  
  # Create plot based on type
  p <- switch(
    plot.type,
    "scatterpie" = plot_deconv_pie(deconv_mat, coords, pt.size, colors, title, coord.flip),
    "scatter" = plot_deconv_scatter(deconv_mat, coords, pt.size, colors, title, coord.flip),
    "heatmap" = plot_celltype_heatmap(deconv_mat, title),
    rlang::abort(paste0("Unknown plot type: ", plot.type))
  )
  
  p
}

#' Scatter Pie Plot for Deconvolution
#'
#' @description Create scatter pie plot showing cell type composition at each spot.
#'
#' @param deconv_mat Deconvolution matrix (spots x celltypes)
#' @param coords Spatial coordinates
#' @param pt.size Point/pie size
#' @param colors Color palette
#' @param title Plot title
#' @param coord.flip Flip coordinates
#'
#' @return A ggplot object
#' @export
plot_deconv_pie <- function(deconv_mat, coords, pt.size = 8, colors = NULL,
                            title = "Cell Type Distribution", coord.flip = FALSE) {
  
  # Check for scatterpie package
  if (!requireNamespace("scatterpie", quietly = TRUE)) {
    rlang::abort("Package 'scatterpie' is required. Install with: install.packages('scatterpie')")
  }
  
  # Prepare data
  celltypes <- colnames(deconv_mat)
  
  plot_data <- data.frame(
    x = as.numeric(coords[, 1]),
    y = as.numeric(coords[, 2]),
    deconv_mat
  )
  rownames(plot_data) <- rownames(coords)
  
  # Ensure numeric
  for (ct in celltypes) {
    plot_data[[ct]] <- as.numeric(plot_data[[ct]])
  }
  
  if (coord.flip) {
    tmp <- plot_data$x
    plot_data$x <- plot_data$y
    plot_data$y <- tmp
  }
  
  # Set colors
  if (is.null(colors)) {
    colors <- get_extended_colors(length(celltypes))
  }
  names(colors) <- celltypes
  
  # Create plot
  p <- ggplot2::ggplot() +
    scatterpie::geom_scatterpie(
      data = plot_data,
      ggplot2::aes(x = x, y = y),
      cols = celltypes,
      pie_scale = pt.size / 4,
      color = "grey30",
      linewidth = 0.1
    ) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 10),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title = title,
      x = if(coord.flip) colnames(coords)[2] else colnames(coords)[1],
      y = if(coord.flip) colnames(coords)[1] else colnames(coords)[2],
      fill = "Cell Type"
    )
  
  p
}

#' Scatter Plot for Deconvolution
#'
#' @description Create scatter plot colored by dominant cell type.
#'
#' @param deconv_mat Deconvolution matrix (spots x celltypes)
#' @param coords Spatial coordinates
#' @param pt.size Point size
#' @param colors Color palette
#' @param title Plot title
#' @param coord.flip Flip coordinates
#'
#' @return A ggplot object
#' @export
plot_deconv_scatter <- function(deconv_mat, coords, pt.size = 2, colors = NULL,
                                 title = "Dominant Cell Type", coord.flip = FALSE) {
  
  # Get dominant cell type for each spot
  dominant_ct <- colnames(deconv_mat)[apply(deconv_mat, 1, which.max)]
  
  # Prepare data
  plot_data <- data.frame(
    x = as.numeric(coords[, 1]),
    y = as.numeric(coords[, 2]),
    celltype = factor(dominant_ct)
  )
  
  if (coord.flip) {
    tmp <- plot_data$x
    plot_data$x <- plot_data$y
    plot_data$y <- tmp
  }
  
  # Set colors
  celltypes <- levels(plot_data$celltype)
  if (is.null(colors)) {
    colors <- get_extended_colors(length(celltypes))
  }
  names(colors) <- celltypes
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = celltype)) +
    ggplot2::geom_point(size = pt.size) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title = title,
      x = if(coord.flip) colnames(coords)[2] else colnames(coords)[1],
      y = if(coord.flip) colnames(coords)[1] else colnames(coords)[2],
      color = "Cell Type"
    )
  
  p
}

#' Plot Cluster Results
#'
#' @description Visualize spatial clustering results.
#'
#' @param cluster_result Clustering result from spatial_cluster()
#' @param pt.size Point size (default: 4)
#' @param colors Color palette
#' @param title Plot title
#'
#' @return A ggplot object
#' @export
plot_cluster <- function(cluster_result, pt.size = 4, colors = NULL, 
                         title = NULL) {
  
  # Validate input
  if (!is.data.frame(cluster_result) || !"cluster" %in% colnames(cluster_result)) {
    rlang::abort("cluster_result must be a data frame with 'cluster' column")
  }
  
  plot_data <- cluster_result
  
  # Set colors
  n_clusters <- length(unique(plot_data$cluster))
  if (is.null(colors)) {
    colors <- get_extended_colors(n_clusters)
  }
  
  # Set title
  if (is.null(title)) {
    title <- paste0("Spatial Clusters (k=", n_clusters, ")")
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = factor(cluster))) +
    ggplot2::geom_point(size = pt.size) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title = title,
      x = "X",
      y = "Y",
      color = "Cluster"
    )
  
  p
}

#' Plot Cell Type Composition Heatmap
#'
#' @description Create heatmap of cell type composition across spots.
#'
#' @param deconv_mat Deconvolution matrix (spots x celltypes)
#' @param title Plot title
#' @param scale.by Scale values: "row", "column", or "none"
#'
#' @return A ggplot object
#' @export
plot_celltype_heatmap <- function(deconv_mat, title = "Cell Type Composition",
                                   scale.by = "none") {
  
  # Apply scaling if requested
  if (scale.by == "row") {
    deconv_mat <- t(scale(t(deconv_mat)))
  } else if (scale.by == "column") {
    deconv_mat <- scale(deconv_mat)
  }
  
  # Handle NaN from scaling
  deconv_mat[is.nan(deconv_mat)] <- 0
  
  # Convert to long format
  df <- as.data.frame(as.matrix(deconv_mat))
  df$spot <- factor(rownames(df), levels = rownames(df))
  
  celltypes <- setdiff(colnames(df), "spot")
  
  df_long <- stats::reshape(
    df,
    direction = "long",
    varying = celltypes,
    v.names = "fraction",
    timevar = "celltype",
    times = celltypes
  )
  
  # Cluster spots for better visualization
  if (nrow(deconv_mat) > 1) {
    spot_order <- stats::hclust(stats::dist(deconv_mat))$order
    df_long$spot <- factor(df_long$spot, levels = rownames(deconv_mat)[spot_order])
  }
  
  # Create heatmap
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = spot, y = celltype, fill = fraction)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(option = "plasma") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = title,
      x = "Spots",
      y = "Cell Type",
      fill = "Fraction"
    )
  
  p
}

#' Plot Topic Distribution Heatmap
#'
#' @description Create heatmap of topic distributions across samples.
#'
#' @param topic_mat Topic distribution matrix (topics x cells/spots)
#' @param celltypes Named vector of cell type annotations (optional)
#' @param title Plot title
#' @param cluster.rows Cluster rows (default: TRUE)
#' @param cluster.cols Cluster columns (default: FALSE)
#'
#' @return A ggplot object
#' @export
plot_topic_heatmap <- function(topic_mat, celltypes = NULL, title = "Topic Distribution",
                                cluster.rows = TRUE, cluster.cols = FALSE) {
  
  # Ensure topics are rows
  if (!is.null(rownames(topic_mat)) && grepl("Topic", rownames(topic_mat)[1], ignore.case = TRUE)) {
    # Already in correct orientation
  } else if (!is.null(colnames(topic_mat)) && grepl("Topic", colnames(topic_mat)[1], ignore.case = TRUE)) {
    topic_mat <- t(topic_mat)
  }
  
  # Convert to data frame
  df <- as.data.frame(as.matrix(topic_mat))
  df$topic <- rownames(df)
  
  samples <- setdiff(colnames(df), "topic")
  
  df_long <- stats::reshape(
    df,
    direction = "long",
    varying = samples,
    v.names = "value",
    timevar = "sample",
    times = samples
  )
  
  # Clustering
  if (cluster.rows && nrow(topic_mat) > 1) {
    row_order <- stats::hclust(stats::dist(topic_mat))$order
    df_long$topic <- factor(df_long$topic, levels = rownames(topic_mat)[row_order])
  }
  
  if (cluster.cols && ncol(topic_mat) > 1) {
    col_order <- stats::hclust(stats::dist(t(topic_mat)))$order
    df_long$sample <- factor(df_long$sample, levels = colnames(topic_mat)[col_order])
  }
  
  # Create heatmap
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = sample, y = topic, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(option = "viridis") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = title,
      x = "Samples",
      y = "Topics",
      fill = "Probability"
    )
  
  p
}

#' Plot Model Selection Results
#'
#' @description Visualize model evaluation metrics across topic numbers.
#'
#' @param metrics Metrics data frame from evaluate_topic_range()
#' @param metric Which metric to plot: "accuracy", "coherence", or "nmi"
#'
#' @return A ggplot object
#' @export
plot_model_selection <- function(metrics, metric = "accuracy") {
  
  if (metric == "accuracy") {
    # Get accuracy columns
    acc_cols <- grep("_accuracy$", colnames(metrics), value = TRUE)
    
    if (length(acc_cols) == 0) {
      rlang::abort("No accuracy columns found in metrics")
    }
    
    df_long <- stats::reshape(
      metrics[, c("n_topics", acc_cols)],
      direction = "long",
      varying = acc_cols,
      v.names = "accuracy",
      timevar = "method",
      times = gsub("_accuracy$", "", acc_cols)
    )
    
    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = factor(n_topics), y = accuracy, fill = method)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::scale_fill_viridis_d() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Model Accuracy by Topic Number",
        x = "Number of Topics",
        y = "Accuracy",
        fill = "Method"
      )
    
  } else if (metric == "coherence") {
    p <- ggplot2::ggplot(metrics, ggplot2::aes(x = factor(n_topics), y = coherence)) +
      ggplot2::geom_bar(stat = "identity", fill = "#3C5488") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Model Coherence by Topic Number",
        x = "Number of Topics",
        y = "Coherence (UMass)"
      )
    
  } else if (metric == "nmi") {
    if (!"nmi" %in% colnames(metrics)) {
      rlang::abort("NMI column not found in metrics")
    }
    
    p <- ggplot2::ggplot(metrics, ggplot2::aes(x = factor(n_topics), y = nmi)) +
      ggplot2::geom_bar(stat = "identity", fill = "#E64B35") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Normalized Mutual Information by Topic Number",
        x = "Number of Topics",
        y = "NMI"
      )
    
  } else {
    rlang::abort(paste0("Unknown metric: ", metric))
  }
  
  p
}

#' Combine Multiple Plots
#'
#' @description Combine multiple ggplot objects into a single figure.
#'
#' @param plots List of ggplot objects
#' @param ncol Number of columns
#' @param nrow Number of rows (optional)
#' @param title Overall title
#'
#' @return A combined plot
#' @export
combine_plots <- function(plots, ncol = 2, nrow = NULL, title = NULL) {
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    rlang::warn("Package 'patchwork' recommended. Returning first plot.")
    return(plots[[1]])
  }
  
  combined <- patchwork::wrap_plots(plots, ncol = ncol, nrow = nrow)
  
  if (!is.null(title)) {
    combined <- combined + patchwork::plot_annotation(title = title)
  }
  
  combined
}

#' Save Plot to File
#'
#' @description Save a ggplot object to file with sensible defaults.
#'
#' @param plot A ggplot object
#' @param filename Output filename
#' @param width Plot width in inches (default: 8)
#' @param height Plot height in inches (default: 6)
#' @param dpi Resolution (default: 300)
#'
#' @export
save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
  vmessage("Plot saved to: ", filename)
}
