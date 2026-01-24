# Tests for integration functions

test_that("center_coords correctly centers coordinates", {
  coords <- data.frame(
    x = c(10, 20, 30),
    y = c(10, 20, 30),
    row.names = c("S1", "S2", "S3")
  )
  
  # Uniform weights
  weights <- rep(1/3, 3)
  centered <- center_coords(coords, weights)
  
  # Centroid should be at (20, 20), so centered should be (-10, 0, 10)
  expect_equal(mean(centered$x), 0, tolerance = 1e-10)
  expect_equal(mean(centered$y), 0, tolerance = 1e-10)
})

test_that("compute_distance_matrix returns valid distances", {
  coords <- data.frame(
    x = c(0, 3, 0),
    y = c(0, 4, 4),
    row.names = c("S1", "S2", "S3")
  )
  
  D <- compute_distance_matrix(coords)
  
  expect_equal(dim(D), c(3, 3))
  expect_equal(D[1, 1], 0)  # Self-distance
  expect_equal(D[1, 2], 5)  # 3-4-5 triangle
  expect_equal(D[1, 3], 4)  # Vertical distance
  expect_equal(D[2, 3], 3)  # Horizontal distance
  
  # Should be symmetric
  expect_equal(D[1, 2], D[2, 1])
})

test_that("compute_kl_divergence produces valid KL matrix", {
  P <- matrix(c(0.5, 0.5, 0.9, 0.1), nrow = 2, byrow = TRUE)
  Q <- matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2, byrow = TRUE)
  
  KL <- compute_kl_divergence(P, Q)
  
  expect_equal(dim(KL), c(2, 2))
  expect_true(all(KL >= 0))
  
  # KL divergence of uniform from uniform should be 0
  expect_equal(KL[1, 1], 0, tolerance = 1e-10)
})

test_that("align_samples returns valid alignment", {
  set.seed(42)
  n1 <- 10
  n2 <- 12
  n_topics <- 3
  
  # Topic distributions
  topic1 <- matrix(runif(n_topics * n1), nrow = n_topics, ncol = n1)
  topic1 <- topic1 / rep(colSums(topic1), each = n_topics)
  rownames(topic1) <- paste0("Topic", 1:n_topics)
  colnames(topic1) <- paste0("S1_", 1:n1)
  
  topic2 <- matrix(runif(n_topics * n2), nrow = n_topics, ncol = n2)
  topic2 <- topic2 / rep(colSums(topic2), each = n_topics)
  rownames(topic2) <- paste0("Topic", 1:n_topics)
  colnames(topic2) <- paste0("S2_", 1:n2)
  
  # Coordinates
  coords1 <- data.frame(
    x = runif(n1),
    y = runif(n1),
    row.names = paste0("S1_", 1:n1)
  )
  
  coords2 <- data.frame(
    x = runif(n2),
    y = runif(n2),
    row.names = paste0("S2_", 1:n2)
  )
  
  result <- align_samples(topic1, topic2, coords1, coords2, 
                          alpha = 0.5, reg = 0.1, max.iter = 30)
  
  # Check structure
  expect_true(is.list(result))
  expect_true("transport" %in% names(result))
  expect_true("rotation" %in% names(result))
  expect_true("fgw_dist" %in% names(result))
  
  # Check transport plan
  T_mat <- result$transport
  expect_equal(dim(T_mat), c(n1, n2))
  expect_true(all(T_mat >= 0))
  
  # Check rotation matrix is orthogonal
  R <- result$rotation
  expect_equal(dim(R), c(2, 2))
  expect_true(all(abs(t(R) %*% R - diag(2)) < 1e-6))
  
  # Distance should be finite
  expect_true(is.finite(result$fgw_dist))
})

test_that("integrate_samples combines multiple samples", {
  set.seed(42)
  n_spots <- 10
  n_celltypes <- 3
  n_topics <- 3
  n_samples <- 3
  
  # Create sample data
  deconv_list <- lapply(1:n_samples, function(i) {
    mat <- matrix(runif(n_spots * n_celltypes), nrow = n_spots)
    mat <- mat / rowSums(mat)
    rownames(mat) <- paste0("S", i, "_Spot", 1:n_spots)
    colnames(mat) <- paste0("CT", 1:n_celltypes)
    mat
  })
  
  topic_list <- lapply(1:n_samples, function(i) {
    mat <- matrix(runif(n_topics * n_spots), nrow = n_topics)
    mat <- mat / rep(colSums(mat), each = n_topics)
    rownames(mat) <- paste0("Topic", 1:n_topics)
    colnames(mat) <- paste0("S", i, "_Spot", 1:n_spots)
    mat
  })
  
  coords_list <- lapply(1:n_samples, function(i) {
    data.frame(
      x = runif(n_spots) * 10 + (i - 1) * 20,  # Offset by sample
      y = runif(n_spots) * 10,
      row.names = paste0("S", i, "_Spot", 1:n_spots)
    )
  })
  
  result <- integrate_samples(
    deconv_list, topic_list, coords_list,
    sample.ids = paste0("Sample", 1:n_samples),
    alpha = 0.5, reg = 0.1, max.iter = 30, seed = 42
  )
  
  # Check structure
  expect_s3_class(result, "strider_integration")
  expect_true("aligned_coords" %in% names(result))
  expect_true("transport_plans" %in% names(result))
  expect_true("deconv" %in% names(result))
  
  # Check aligned coordinates
  expect_equal(nrow(result$aligned_coords), n_spots * n_samples)
  expect_true("x_aligned" %in% colnames(result$aligned_coords))
  expect_true("y_aligned" %in% colnames(result$aligned_coords))
  expect_true("sample" %in% colnames(result$aligned_coords))
  
  # Check transport plans (n_samples - 1 pairs)
  expect_equal(length(result$transport_plans), n_samples - 1)
  
  # Check combined deconv
  expect_equal(nrow(result$deconv), n_spots * n_samples)
  expect_equal(ncol(result$deconv), n_celltypes)
})

test_that("integration_quality computes valid metrics", {
  set.seed(42)
  
  # Create mock integration result
  n1 <- 5
  n2 <- 6
  
  T_mat <- matrix(runif(n1 * n2), nrow = n1)
  T_mat <- T_mat / sum(T_mat)  # Normalize
  
  result <- list(
    transport_plans = list("S1_to_S2" = T_mat),
    class = "strider_integration"
  )
  class(result) <- "strider_integration"
  
  metrics <- integration_quality(result)
  
  expect_true(is.data.frame(metrics))
  expect_true("pair" %in% colnames(metrics))
  expect_true("entropy" %in% colnames(metrics))
  expect_true("sparsity" %in% colnames(metrics))
  
  # Entropy should be non-negative
  expect_true(metrics$entropy >= 0)
  
  # Sparsity should be between 0 and 1
  expect_true(metrics$sparsity >= 0)
  expect_true(metrics$sparsity <= 1)
})
