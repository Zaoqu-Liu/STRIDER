# Tests for clustering functions

test_that("find_neighbors works correctly", {
  coords <- data.frame(
    x = c(0, 1, 2, 0, 1, 2, 3, 4, 5),
    y = c(0, 0, 0, 1, 1, 1, 2, 2, 2),
    row.names = paste0("Spot", 1:9)
  )
  
  neighbors <- find_neighbors(coords, k = 2)
  
  expect_equal(length(neighbors), 9)
  expect_true(all(lengths(neighbors) >= 1))
})

test_that("calculate_neighbor_fractions preserves dimensions", {
  deconv <- matrix(runif(36), nrow = 9, ncol = 4)
  rownames(deconv) <- paste0("Spot", 1:9)
  colnames(deconv) <- paste0("CT", 1:4)
  
  coords <- data.frame(
    x = c(0, 1, 2, 0, 1, 2, 3, 4, 5),
    y = c(0, 0, 0, 1, 1, 1, 2, 2, 2),
    row.names = paste0("Spot", 1:9)
  )
  
  neighbors <- find_neighbors(coords, k = 2)
  neighbor_frac <- calculate_neighbor_fractions(deconv, neighbors)
  
  expect_equal(dim(neighbor_frac), dim(deconv))
})

test_that("spatial_cluster returns correct structure", {
  # Create enough data points for clustering
  set.seed(42)
  n_spots <- 20
  
  # Create deconvolution matrix with some pattern
  deconv <- matrix(0, nrow = n_spots, ncol = 4)
  deconv[1:10, 1:2] <- runif(20, 0.3, 0.7)
  deconv[11:20, 3:4] <- runif(20, 0.3, 0.7)
  deconv <- deconv / rowSums(deconv)  # Normalize
  rownames(deconv) <- paste0("Spot", 1:n_spots)
  colnames(deconv) <- paste0("CT", 1:4)
  
  # Create coordinates on a grid
  coords <- data.frame(
    x = rep(1:5, 4),
    y = rep(1:4, each = 5),
    row.names = paste0("Spot", 1:n_spots)
  )
  
  result <- spatial_cluster(deconv, coords, n.clusters = 2, seed = 42)
  
  expect_true("cluster" %in% colnames(result))
  expect_equal(nrow(result), n_spots)
  expect_equal(length(unique(result$cluster)), 2)
})

test_that("cluster_composition calculates means correctly", {
  deconv <- matrix(c(0.8, 0.2, 0.7, 0.3, 0.2, 0.8, 0.3, 0.7), 
                   nrow = 4, ncol = 2, byrow = TRUE)
  rownames(deconv) <- paste0("Spot", 1:4)
  colnames(deconv) <- c("CT1", "CT2")
  
  clusters <- c(1, 1, 2, 2)
  names(clusters) <- paste0("Spot", 1:4)
  
  comp <- cluster_composition(deconv, clusters)
  
  expect_equal(nrow(comp), 2)
  expect_equal(ncol(comp), 2)
  # Cluster 1 should have higher CT1
  expect_true(comp["Cluster_1", "CT1"] > comp["Cluster_1", "CT2"])
  # Cluster 2 should have higher CT2
  expect_true(comp["Cluster_2", "CT2"] > comp["Cluster_2", "CT1"])
})
