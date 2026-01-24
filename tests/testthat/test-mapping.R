# Tests for spot-cell mapping functions

test_that("map_spots_to_cells finds similar cells", {
  set.seed(42)
  n_topics <- 5
  n_spots <- 10
  n_cells <- 30
  
  # Create topic distributions
  topic_spot <- matrix(runif(n_topics * n_spots), nrow = n_topics)
  topic_spot <- topic_spot / rep(colSums(topic_spot), each = n_topics)
  rownames(topic_spot) <- paste0("Topic", 1:n_topics)
  colnames(topic_spot) <- paste0("Spot", 1:n_spots)
  
  topic_cell <- matrix(runif(n_topics * n_cells), nrow = n_topics)
  topic_cell <- topic_cell / rep(colSums(topic_cell), each = n_topics)
  rownames(topic_cell) <- paste0("Topic", 1:n_topics)
  colnames(topic_cell) <- paste0("Cell", 1:n_cells)
  
  # Make first 3 cells very similar to first spot
  topic_cell[, 1:3] <- topic_spot[, 1] + matrix(rnorm(n_topics * 3, sd = 0.01), 
                                                 nrow = n_topics)
  topic_cell[, 1:3] <- topic_cell[, 1:3] / rep(colSums(topic_cell[, 1:3]), 
                                                each = n_topics)
  
  result <- map_spots_to_cells(topic_spot, topic_cell, n.top = 5, method = "euclidean")
  
  # Check structure
  expect_true(is.data.frame(result))
  expect_true(all(c("spot", "rank", "cell", "distance") %in% colnames(result)))
  
  # Each spot should have n.top mappings
  expect_equal(nrow(result), n_spots * 5)
  
  # Distances should be non-negative
  expect_true(all(result$distance >= 0))
  
  # First spot's top cells should include Cell1-3
  spot1_cells <- result[result$spot == "Spot1" & result$rank <= 3, "cell"]
  expect_true(any(spot1_cells %in% c("Cell1", "Cell2", "Cell3")))
})

test_that("map_spots_to_cells works with different distance methods", {
  set.seed(42)
  n_topics <- 3
  n_spots <- 5
  n_cells <- 10
  
  topic_spot <- matrix(runif(n_topics * n_spots), nrow = n_topics)
  topic_spot <- topic_spot / rep(colSums(topic_spot), each = n_topics)
  rownames(topic_spot) <- paste0("Topic", 1:n_topics)
  colnames(topic_spot) <- paste0("Spot", 1:n_spots)
  
  topic_cell <- matrix(runif(n_topics * n_cells), nrow = n_topics)
  topic_cell <- topic_cell / rep(colSums(topic_cell), each = n_topics)
  rownames(topic_cell) <- paste0("Topic", 1:n_topics)
  colnames(topic_cell) <- paste0("Cell", 1:n_cells)
  
  methods <- c("euclidean", "cosine", "correlation")
  
  for (method in methods) {
    result <- map_spots_to_cells(topic_spot, topic_cell, n.top = 3, method = method)
    
    expect_equal(nrow(result), n_spots * 3,
                 label = paste("Method:", method))
    expect_true(all(result$distance >= 0),
                label = paste("Method:", method))
  }
})

test_that("get_mapped_celltypes aggregates correctly", {
  # Create mapping result
  mapping <- data.frame(
    spot = rep(c("Spot1", "Spot2"), each = 5),
    rank = rep(1:5, 2),
    cell = c("Cell1", "Cell2", "Cell3", "Cell4", "Cell5",
             "Cell6", "Cell7", "Cell8", "Cell9", "Cell10"),
    distance = runif(10, 0.1, 0.5)
  )
  
  # Cell type annotations
  celltypes <- c(Cell1 = "TypeA", Cell2 = "TypeA", Cell3 = "TypeA",
                 Cell4 = "TypeB", Cell5 = "TypeB",
                 Cell6 = "TypeB", Cell7 = "TypeB", Cell8 = "TypeB",
                 Cell9 = "TypeC", Cell10 = "TypeC")
  
  # Test majority vote
  result_majority <- get_mapped_celltypes(mapping, celltypes, aggregate = "majority")
  
  expect_true(is.data.frame(result_majority))
  expect_equal(nrow(result_majority), 2)  # 2 spots
  expect_true("predicted_celltype" %in% colnames(result_majority))
  
  # Spot1: 3 TypeA, 2 TypeB -> TypeA
  expect_equal(result_majority[result_majority$spot == "Spot1", "predicted_celltype"], "TypeA")
  
  # Spot2: 3 TypeB, 2 TypeC -> TypeB
  expect_equal(result_majority[result_majority$spot == "Spot2", "predicted_celltype"], "TypeB")
})

test_that("transfer_labels produces valid labels", {
  mapping <- data.frame(
    spot = rep(c("Spot1", "Spot2", "Spot3"), each = 3),
    rank = rep(1:3, 3),
    cell = paste0("Cell", 1:9),
    distance = runif(9, 0.1, 0.5)
  )
  
  celltypes <- setNames(
    c(rep("TypeA", 3), rep("TypeB", 3), rep("TypeC", 3)),
    paste0("Cell", 1:9)
  )
  
  labels <- transfer_labels(mapping, celltypes, method = "knn", k = 3)
  
  expect_equal(length(labels), 3)
  expect_true(all(names(labels) %in% c("Spot1", "Spot2", "Spot3")))
  expect_true(all(labels %in% c("TypeA", "TypeB", "TypeC")))
})

test_that("compute_cosine_distance is valid", {
  X <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  Y <- matrix(c(1, 0, 0, 1, 1, 1), nrow = 3, byrow = TRUE)
  
  D <- STRIDER:::compute_cosine_distance(X, Y)
  
  expect_equal(dim(D), c(2, 3))
  
  # Identical vectors should have distance ~0
  expect_equal(D[1, 1], 0, tolerance = 1e-6)  # (1,0) vs (1,0)
  expect_equal(D[2, 2], 0, tolerance = 1e-6)  # (0,1) vs (0,1)
  
  # Orthogonal vectors should have distance 1
  expect_equal(D[1, 2], 1, tolerance = 1e-6)  # (1,0) vs (0,1)
  expect_equal(D[2, 1], 1, tolerance = 1e-6)  # (0,1) vs (1,0)
})

test_that("compute_correlation_distance is valid", {
  X <- matrix(c(1, 2, 3, 3, 2, 1), nrow = 2, byrow = TRUE)
  
  D <- STRIDER:::compute_correlation_distance(X, X)
  
  expect_equal(dim(D), c(2, 2))
  
  # Perfect correlation with self
  expect_equal(D[1, 1], 0, tolerance = 1e-10)
  expect_equal(D[2, 2], 0, tolerance = 1e-10)
  
  # Perfect negative correlation -> distance = 2
  expect_equal(D[1, 2], 2, tolerance = 1e-10)
})
