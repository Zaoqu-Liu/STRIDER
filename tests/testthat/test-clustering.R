# Tests for clustering functions

test_that("find_neighbors works correctly", {
  coords <- data.frame(
    x = c(0, 1, 2, 0, 1, 2),
    y = c(0, 0, 0, 1, 1, 1),
    row.names = paste0("Spot", 1:6)
  )
  
  neighbors <- find_neighbors(coords, k = 2)
  
  expect_equal(length(neighbors), 6)
  expect_true(all(lengths(neighbors) >= 1))
})

test_that("calculate_neighbor_fractions preserves dimensions", {
  deconv <- matrix(runif(24), nrow = 6, ncol = 4)
  rownames(deconv) <- paste0("Spot", 1:6)
  colnames(deconv) <- paste0("CT", 1:4)
  
  coords <- data.frame(
    x = c(0, 1, 2, 0, 1, 2),
    y = c(0, 0, 0, 1, 1, 1),
    row.names = paste0("Spot", 1:6)
  )
  
  neighbors <- find_neighbors(coords, k = 2)
  neighbor_frac <- calculate_neighbor_fractions(deconv, neighbors)
  
  expect_equal(dim(neighbor_frac), dim(deconv))
})

test_that("spatial_cluster returns correct structure", {
  deconv <- matrix(runif(24), nrow = 6, ncol = 4)
  rownames(deconv) <- paste0("Spot", 1:6)
  colnames(deconv) <- paste0("CT", 1:4)
  
  coords <- data.frame(
    x = c(0, 1, 2, 0, 1, 2),
    y = c(0, 0, 0, 1, 1, 1),
    row.names = paste0("Spot", 1:6)
  )
  
  result <- spatial_cluster(deconv, coords, n.clusters = 2)
  
  expect_true("cluster" %in% colnames(result))
  expect_equal(nrow(result), 6)
  expect_equal(length(unique(result$cluster)), 2)
})
