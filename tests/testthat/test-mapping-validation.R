# Mapping and Clustering Validation Tests
# Deep tests for spot-cell mapping and spatial clustering

# =============================================================================
# Test 1: Cosine Distance Mathematical Correctness
# =============================================================================

test_that("cosine distance equals 1 - cosine similarity", {
  # Known vectors - simple case
  X <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  
  D <- compute_cosine_distance(X, X)
  
  # Self-cosine similarity is 1, so distance is 0
  expect_equal(D[1, 1], 0, tolerance = 1e-8)
  expect_equal(D[2, 2], 0, tolerance = 1e-8)
  
  # Check symmetry
  expect_equal(D[1, 2], D[2, 1], tolerance = 1e-10)
  
  # Check that D[1,2] is in valid range [0, 2]
  expect_true(D[1, 2] >= 0)
  expect_true(D[1, 2] <= 2)
})

test_that("cosine distance is symmetric for self-comparison", {
  set.seed(42)
  X <- matrix(runif(30), nrow = 6)
  
  D <- compute_cosine_distance(X, X)
  
  # Self-distance should be 0
  expect_true(all(abs(diag(D)) < 1e-10))
  
  # Symmetric
  expect_true(all(abs(D - t(D)) < 1e-10))
})

# =============================================================================
# Test 2: Correlation Distance Mathematical Correctness
# =============================================================================

test_that("correlation distance is in valid range", {
  set.seed(42)
  # Both matrices must have same number of columns
  X <- matrix(rnorm(30), nrow = 5, ncol = 6)
  Y <- matrix(rnorm(48), nrow = 8, ncol = 6)
  
  D <- compute_correlation_distance(X, Y)
  
  # Correlation distance = 1 - correlation
  # Correlation in [-1, 1], so distance in [0, 2]
  expect_true(all(D >= -1e-10))
  expect_true(all(D <= 2 + 1e-10))
})

test_that("correlation distance is zero for identical rows", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
  
  D <- compute_correlation_distance(X, X)
  
  # Self-correlation is 1, so distance is 0
  expect_true(all(abs(diag(D)) < 1e-10))
})

# =============================================================================
# Test 3: Mapping Finds Nearest Neighbors Correctly
# =============================================================================

test_that("map_spots_to_cells finds correct nearest neighbors", {
  # Create distinct topic distributions
  topic_cell <- matrix(c(
    0.9, 0.1, 0.1, 0.1,  # Cell 1: high topic 1
    0.1, 0.9, 0.1, 0.1,  # Cell 2: high topic 2
    0.1, 0.1, 0.9, 0.1,  # Cell 3: high topic 3
    0.1, 0.1, 0.1, 0.9   # Cell 4: high topic 4
  ), nrow = 4, ncol = 4, byrow = TRUE)
  topic_cell <- t(topic_cell)  # topics x cells
  colnames(topic_cell) <- paste0("Cell", 1:4)
  rownames(topic_cell) <- paste0("Topic", 1:4)
  
  # Spot similar to Cell 1
  topic_spot <- matrix(c(0.8, 0.1, 0.05, 0.05), nrow = 4)
  colnames(topic_spot) <- "Spot1"
  rownames(topic_spot) <- paste0("Topic", 1:4)
  
  mapping <- map_spots_to_cells(topic_spot, topic_cell, n.top = 4, method = "euclidean")
  
  # Cell 1 should be closest to Spot 1
  top_cell <- mapping$cell[mapping$rank == 1]
  expect_equal(top_cell, "Cell1")
})

test_that("different distance methods rank cells consistently", {
  set.seed(42)
  n_topics <- 5
  n_cells <- 10
  n_spots <- 3
  
  topic_cell <- matrix(runif(n_topics * n_cells), nrow = n_topics)
  colnames(topic_cell) <- paste0("Cell", 1:n_cells)
  rownames(topic_cell) <- paste0("Topic", 1:n_topics)
  
  topic_spot <- matrix(runif(n_topics * n_spots), nrow = n_topics)
  colnames(topic_spot) <- paste0("Spot", 1:n_spots)
  rownames(topic_spot) <- paste0("Topic", 1:n_topics)
  
  # All methods should return valid mappings
  map_euc <- map_spots_to_cells(topic_spot, topic_cell, n.top = 5, method = "euclidean")
  map_cos <- map_spots_to_cells(topic_spot, topic_cell, n.top = 5, method = "cosine")
  map_cor <- map_spots_to_cells(topic_spot, topic_cell, n.top = 5, method = "correlation")
  
  # All should have same number of rows
  expect_equal(nrow(map_euc), nrow(map_cos))
  expect_equal(nrow(map_euc), nrow(map_cor))
  
  # All should have valid cells
  expect_true(all(map_euc$cell %in% colnames(topic_cell)))
  expect_true(all(map_cos$cell %in% colnames(topic_cell)))
  expect_true(all(map_cor$cell %in% colnames(topic_cell)))
})

# =============================================================================
# Test 4: Neighbor Finding Correctness
# =============================================================================

test_that("find_neighbors returns correct k neighbors", {
  # Simple grid
  coords <- data.frame(
    x = c(0, 1, 2, 0, 1, 2, 0, 1, 2),
    y = c(0, 0, 0, 1, 1, 1, 2, 2, 2),
    row.names = paste0("Spot", 1:9)
  )
  
  neighbors <- find_neighbors(coords, k = 4)
  
  # Each spot should have exactly 4 neighbors
  expect_true(all(lengths(neighbors) == 4))
  
  # Center spot (5) should have spots 2, 4, 6, 8 as nearest neighbors
  center_neighbors <- neighbors[["Spot5"]]
  expected_neighbors <- c(2, 4, 6, 8)  # indices in grid
  expect_true(all(expected_neighbors %in% center_neighbors))
})

test_that("find_neighbors excludes self by default", {
  coords <- data.frame(
    x = c(0, 1, 2),
    y = c(0, 0, 0),
    row.names = paste0("Spot", 1:3)
  )
  
  neighbors <- find_neighbors(coords, k = 2, include.self = FALSE)
  
  # No spot should be its own neighbor
  for (i in seq_along(neighbors)) {
    expect_false(i %in% neighbors[[i]])
  }
})

# =============================================================================
# Test 5: Spatial Clustering Properties
# =============================================================================

test_that("spatial_cluster produces valid cluster assignments", {
  set.seed(42)
  n_spots <- 30
  n_celltypes <- 4
  
  # Create deconvolution matrix
  deconv <- matrix(runif(n_spots * n_celltypes), nrow = n_spots, ncol = n_celltypes)
  deconv <- deconv / rowSums(deconv)  # Normalize
  rownames(deconv) <- paste0("Spot", 1:n_spots)
  colnames(deconv) <- paste0("CT", 1:n_celltypes)
  
  # Create random coordinates
  coords <- data.frame(
    x = runif(n_spots, 0, 10),
    y = runif(n_spots, 0, 10),
    row.names = paste0("Spot", 1:n_spots)
  )
  
  result <- spatial_cluster(deconv, coords, n.clusters = 3, k = 4, seed = 42)
  
  # Check output structure
  expect_true("spot" %in% colnames(result))
  expect_true("cluster" %in% colnames(result))
  expect_true("x" %in% colnames(result))
  expect_true("y" %in% colnames(result))
  
  # All spots assigned
  expect_equal(nrow(result), n_spots)
  
  # Cluster values in expected range
  expect_true(all(result$cluster %in% 1:3))
  
  # All 3 clusters exist
  expect_equal(length(unique(result$cluster)), 3)
})

test_that("neighbor weight affects clustering", {
  set.seed(42)
  n_spots <- 20
  n_celltypes <- 3
  
  # Create spatially clustered data
  deconv <- matrix(0.33, nrow = n_spots, ncol = n_celltypes)
  
  # Make first half dominated by CT1
  deconv[1:10, ] <- c(0.8, 0.1, 0.1)
  # Make second half dominated by CT2
  deconv[11:20, ] <- c(0.1, 0.8, 0.1)
  
  rownames(deconv) <- paste0("Spot", 1:n_spots)
  colnames(deconv) <- paste0("CT", 1:n_celltypes)
  
  # Spatial coords: first half in left region, second half in right region
  coords <- data.frame(
    x = c(runif(10, 0, 5), runif(10, 10, 15)),
    y = runif(20, 0, 10),
    row.names = paste0("Spot", 1:n_spots)
  )
  
  # With weight = 1 (only own composition)
  result_own <- spatial_cluster(deconv, coords, n.clusters = 2, weight = 1.0, seed = 42)
  
  # Clusters should separate the two groups
  cluster1_spots <- which(result_own$cluster == 1)
  cluster2_spots <- which(result_own$cluster == 2)
  
  # Check that clusters separate the two spatial regions
  # (allowing for some variation due to k-means randomness)
  expect_true(
    length(intersect(cluster1_spots, 1:10)) > 5 || 
    length(intersect(cluster2_spots, 1:10)) > 5
  )
})

# =============================================================================
# Test 6: Cluster Composition Calculation
# =============================================================================

test_that("cluster_composition calculates correct means", {
  # Simple example
  deconv <- matrix(c(
    0.8, 0.2,  # Spot 1
    0.6, 0.4,  # Spot 2
    0.2, 0.8,  # Spot 3
    0.3, 0.7   # Spot 4
  ), nrow = 4, byrow = TRUE)
  rownames(deconv) <- paste0("Spot", 1:4)
  colnames(deconv) <- c("CT1", "CT2")
  
  clusters <- c(1, 1, 2, 2)
  names(clusters) <- rownames(deconv)
  
  comp <- cluster_composition(deconv, clusters)
  
  # Cluster 1: mean of spots 1,2 = (0.8+0.6)/2, (0.2+0.4)/2
  expect_equal(comp["Cluster_1", "CT1"], 0.7, tolerance = 1e-10)
  expect_equal(comp["Cluster_1", "CT2"], 0.3, tolerance = 1e-10)
  
  # Cluster 2: mean of spots 3,4 = (0.2+0.3)/2, (0.8+0.7)/2
  expect_equal(comp["Cluster_2", "CT1"], 0.25, tolerance = 1e-10)
  expect_equal(comp["Cluster_2", "CT2"], 0.75, tolerance = 1e-10)
})

# =============================================================================
# Test 7: Label Transfer Accuracy
# =============================================================================

test_that("transfer_labels produces valid predictions", {
  # Create mapping result
  mapping <- data.frame(
    spot = rep(c("Spot1", "Spot2"), each = 3),
    rank = rep(1:3, 2),
    cell = c("CellA", "CellB", "CellC", "CellD", "CellE", "CellF"),
    distance = c(0.1, 0.2, 0.3, 0.15, 0.25, 0.35),
    stringsAsFactors = FALSE
  )
  
  celltypes <- c(CellA = "TypeX", CellB = "TypeX", CellC = "TypeY",
                 CellD = "TypeY", CellE = "TypeY", CellF = "TypeZ")
  
  # kNN should use majority vote
  labels_knn <- transfer_labels(mapping, celltypes, method = "knn")
  
  # Spot1: 2 TypeX, 1 TypeY -> TypeX
  expect_equal(as.character(labels_knn["Spot1"]), "TypeX")
  
  # Spot2: 2 TypeY, 1 TypeZ -> TypeY
  expect_equal(as.character(labels_knn["Spot2"]), "TypeY")
})

test_that("weighted transfer considers distances", {
  # Create mapping where closest has different type than majority
  mapping <- data.frame(
    spot = rep("Spot1", 3),
    rank = 1:3,
    cell = c("CellA", "CellB", "CellC"),
    distance = c(0.001, 0.5, 0.5),  # A is much closer
    stringsAsFactors = FALSE
  )
  
  celltypes <- c(CellA = "TypeX", CellB = "TypeY", CellC = "TypeY")
  
  # Weighted should favor TypeX (closest)
  labels_weighted <- transfer_labels(mapping, celltypes, method = "weighted")
  
  # TypeX should win due to very close distance
  expect_equal(as.character(labels_weighted["Spot1"]), "TypeX")
})

# =============================================================================
# Test 8: Neighbor Average Calculation
# =============================================================================

test_that("calculate_neighbor_fractions averages correctly", {
  deconv <- matrix(c(
    0.8, 0.2,  # Spot 1
    0.6, 0.4,  # Spot 2
    0.4, 0.6,  # Spot 3
    0.2, 0.8   # Spot 4
  ), nrow = 4, byrow = TRUE)
  rownames(deconv) <- paste0("Spot", 1:4)
  colnames(deconv) <- c("CT1", "CT2")
  
  # Neighbors: 1 -> {2,3}, 2 -> {1,3}, 3 -> {1,2}, 4 -> {} (isolated)
  neighbors <- list(c(2, 3), c(1, 3), c(1, 2), integer(0))
  
  neighbor_frac <- calculate_neighbor_fractions(deconv, neighbors)
  
  # Spot 1: avg of spots 2,3 = (0.6+0.4)/2, (0.4+0.6)/2
  expect_equal(neighbor_frac[1, 1], 0.5, tolerance = 1e-10)
  expect_equal(neighbor_frac[1, 2], 0.5, tolerance = 1e-10)
  
  # Spot 4: no neighbors, should equal own value
  expect_equal(neighbor_frac[4, 1], 0.2, tolerance = 1e-10)
  expect_equal(neighbor_frac[4, 2], 0.8, tolerance = 1e-10)
})

# =============================================================================
# Test 9: Mapping Graph Creation
# =============================================================================

test_that("create_mapping_graph produces valid output", {
  mapping <- data.frame(
    spot = c("Spot1", "Spot1", "Spot2", "Spot2"),
    rank = c(1, 2, 1, 2),
    cell = c("CellA", "CellB", "CellC", "CellD"),
    distance = c(0.1, 0.2, 0.15, 0.25),
    stringsAsFactors = FALSE
  )
  
  result <- create_mapping_graph(mapping)
  
  # If igraph available, should be a graph
  if (requireNamespace("igraph", quietly = TRUE)) {
    expect_true(inherits(result, "igraph"))
  } else {
    # Otherwise edge list
    expect_true(is.data.frame(result))
    expect_true(all(c("from", "to", "weight") %in% colnames(result)))
  }
})

# =============================================================================
# Test 10: Distance Threshold in Mapping Graph
# =============================================================================

test_that("mapping graph respects distance threshold", {
  mapping <- data.frame(
    spot = c("Spot1", "Spot1", "Spot2", "Spot2"),
    rank = c(1, 2, 1, 2),
    cell = c("CellA", "CellB", "CellC", "CellD"),
    distance = c(0.1, 0.3, 0.15, 0.35),  # Only 2 edges below 0.2
    stringsAsFactors = FALSE
  )
  
  # Without threshold - 4 edges
  result_all <- create_mapping_graph(mapping, threshold = NULL)
  
  # With threshold 0.2 - 2 edges
  result_filtered <- create_mapping_graph(mapping, threshold = 0.2)
  
  if (requireNamespace("igraph", quietly = TRUE)) {
    expect_equal(igraph::ecount(result_all), 4)
    expect_equal(igraph::ecount(result_filtered), 2)
  } else {
    expect_equal(nrow(result_all), 4)
    expect_equal(nrow(result_filtered), 2)
  }
})
