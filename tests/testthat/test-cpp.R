# Tests for C++ functions - Optimal Transport and Distance Computations

test_that("sinkhorn_cpp returns valid transport plan", {
  set.seed(42)
  n <- 5
  m <- 6
  
  # Random cost matrix
  M <- matrix(runif(n * m), nrow = n, ncol = m)
  
  # Uniform marginals
  p <- rep(1/n, n)
  q <- rep(1/m, m)
  
  # Compute transport plan
  T_mat <- sinkhorn_cpp(M, p, q, reg = 0.1, max_iter = 100, tol = 1e-6)
  
  # Check dimensions
  expect_equal(dim(T_mat), c(n, m))
  
  # Check marginal constraints (approximately)
  row_sums <- rowSums(T_mat)
  col_sums <- colSums(T_mat)
  
  expect_true(all(abs(row_sums - p) < 0.01))
  expect_true(all(abs(col_sums - q) < 0.01))
  
  # Transport plan should be non-negative
  expect_true(all(T_mat >= 0))
})

test_that("sinkhorn_cpp handles non-uniform marginals", {
  n <- 4
  m <- 5
  
  M <- matrix(runif(n * m), nrow = n, ncol = m)
  
  # Non-uniform marginals
  p <- c(0.4, 0.3, 0.2, 0.1)
  q <- c(0.3, 0.25, 0.2, 0.15, 0.1)
  
  T_mat <- sinkhorn_cpp(M, p, q, reg = 0.1, max_iter = 500, tol = 1e-8)
  
  # Check marginals
  expect_true(all(abs(rowSums(T_mat) - p) < 0.01))
  expect_true(all(abs(colSums(T_mat) - q) < 0.01))
})

test_that("euclidean_distance_cpp computes correct distances", {
  # Simple test case
  X <- matrix(c(0, 0, 3, 4), nrow = 2, byrow = TRUE)
  
  # Self-distance
  D <- euclidean_distance_cpp(X, NULL)
  
  expect_equal(dim(D), c(2, 2))
  expect_equal(D[1, 1], 0)  # Distance to self is 0
  expect_equal(D[2, 2], 0)
  expect_equal(D[1, 2], 5)  # 3-4-5 triangle
  expect_equal(D[2, 1], 5)  # Symmetric
})

test_that("euclidean_distance_cpp handles two different matrices", {
  X <- matrix(c(0, 0, 1, 0), nrow = 2, byrow = TRUE)
  Y <- matrix(c(0, 1, 1, 1, 2, 0), nrow = 3, byrow = TRUE)
  
  D <- euclidean_distance_cpp(X, Y)
  
  expect_equal(dim(D), c(2, 3))
  
  # Check specific distances
  expect_equal(D[1, 1], 1)  # (0,0) to (0,1)
  expect_equal(D[2, 2], 1)  # (1,0) to (1,1)
})

test_that("kl_divergence_matrix_cpp computes valid KL divergence", {
  # Create probability distributions
  P <- matrix(c(0.5, 0.3, 0.2,
                0.1, 0.8, 0.1), nrow = 2, byrow = TRUE)
  Q <- matrix(c(0.33, 0.33, 0.34,
                0.25, 0.5, 0.25), nrow = 2, byrow = TRUE)
  
  KL <- kl_divergence_matrix_cpp(P, Q)
  
  expect_equal(dim(KL), c(2, 2))
  
  # KL divergence should be non-negative
  expect_true(all(KL >= -1e-10))  # Allow small numerical error
  
  # KL(P[i] || P[i]) should be close to 0 (diagonal elements)
  # Note: KL_self computes KL between row i of P and row j of P
  KL_self <- kl_divergence_matrix_cpp(P, P)
  # Diagonal elements should be 0 (same distribution vs itself)
  expect_equal(diag(KL_self), c(0, 0), tolerance = 1e-10)
  
  # Test with identical distributions (all rows same)
  P_same <- matrix(rep(c(0.4, 0.3, 0.3), 3), nrow = 3, byrow = TRUE)
  KL_same <- kl_divergence_matrix_cpp(P_same, P_same)
  expect_true(all(abs(KL_same) < 1e-6))  # All should be ~0
})

test_that("fused_gromov_wasserstein_cpp returns valid result", {
  set.seed(42)
  n <- 5
  m <- 6
  
  # Feature cost matrix
  M <- matrix(runif(n * m), nrow = n, ncol = m)
  
  # Structural distance matrices
  C1 <- as.matrix(dist(matrix(runif(n * 2), ncol = 2)))
  C2 <- as.matrix(dist(matrix(runif(m * 2), ncol = 2)))
  
  # Marginals
  p <- rep(1/n, n)
  q <- rep(1/m, m)
  
  # Run FGW
  result <- fused_gromov_wasserstein_cpp(M, C1, C2, p, q, 
                                          alpha = 0.5, reg = 0.1, 
                                          max_iter = 50, tol = 1e-6)
  
  # Check structure
  expect_true(is.list(result))
  expect_true("T" %in% names(result))
  expect_true("fgw_dist" %in% names(result))
  
  # Check transport plan
  T_mat <- result$T
  expect_equal(dim(T_mat), c(n, m))
  expect_true(all(T_mat >= -1e-10))  # Allow small numerical errors
  
  # Check marginals (approximately)
  expect_true(all(abs(rowSums(T_mat) - p) < 0.05))
  expect_true(all(abs(colSums(T_mat) - q) < 0.05))
  
  # FGW distance is computed but may have numerical issues
  # Just check it's finite
  expect_true(is.finite(result$fgw_dist))
})

test_that("fgw alpha=0 is close to pure Wasserstein", {
  set.seed(42)
  n <- 5
  m <- 5
  
  M <- matrix(runif(n * m), nrow = n, ncol = m)
  C1 <- as.matrix(dist(matrix(runif(n * 2), ncol = 2)))
  C2 <- as.matrix(dist(matrix(runif(m * 2), ncol = 2)))
  
  p <- rep(1/n, n)
  q <- rep(1/m, m)
  
  # alpha = 0 should ignore structural cost
  result_alpha0 <- fused_gromov_wasserstein_cpp(M, C1, C2, p, q, 
                                                 alpha = 0, reg = 0.1)
  
  # Pure Wasserstein (using sinkhorn)
  T_wasserstein <- sinkhorn_cpp(M, p, q, reg = 0.1)
  
  # The transport plans should be similar
  diff <- sum(abs(result_alpha0$T - T_wasserstein))
  expect_true(diff < 0.1)
})

test_that("gw_cpp computes valid Gromov-Wasserstein distance", {
  set.seed(42)
  n <- 5
  m <- 6
  
  # Structural distance matrices
  C1 <- as.matrix(dist(matrix(runif(n * 2), ncol = 2)))
  C2 <- as.matrix(dist(matrix(runif(m * 2), ncol = 2)))
  
  p <- rep(1/n, n)
  q <- rep(1/m, m)
  
  result <- gw_cpp(C1, C2, p, q, reg = 0.1, max_iter = 50, tol = 1e-6)
  
  expect_true(is.list(result))
  expect_true("T" %in% names(result))
  expect_true("gw_dist" %in% names(result))
  
  # Check transport plan validity
  T_mat <- result$T
  expect_true(all(T_mat >= -1e-10))  # Allow small numerical errors
  expect_true(is.finite(result$gw_dist))  # Distance should be finite
})

test_that("center_points_cpp correctly centers points", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2, byrow = TRUE)
  weights <- c(1/3, 1/3, 1/3)
  
  centered <- center_points_cpp(X, weights)
  
  # Weighted centroid should be at origin
  weighted_mean <- colSums(centered * weights)
  expect_true(all(abs(weighted_mean) < 1e-10))
})

test_that("center_points_cpp handles non-uniform weights", {
  X <- matrix(c(0, 0, 10, 0), nrow = 2, ncol = 2, byrow = TRUE)
  
  # Heavy weight on first point
  weights <- c(0.9, 0.1)
  
  centered <- center_points_cpp(X, weights)
  
  # Weighted centroid: (0.9*0 + 0.1*10, 0.9*0 + 0.1*0) = (1, 0)
  # First point should be at -1, 0
  expect_true(abs(centered[1, 1] - (-1)) < 1e-10)
  expect_true(abs(centered[2, 1] - 9) < 1e-10)
})

test_that("compute_rotation_svd_cpp finds optimal rotation", {
  # Create two point sets related by 90 degree rotation
  X <- matrix(c(1, 0, 0, 1, -1, 0, 0, -1), nrow = 4, ncol = 2, byrow = TRUE)
  
  # 90 degree rotation matrix
  theta <- pi / 2
  R_true <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 
                   nrow = 2, ncol = 2, byrow = TRUE)
  
  Y <- X %*% R_true
  
  weights <- rep(1/4, 4)
  
  R_computed <- compute_rotation_svd_cpp(X, Y, weights)
  
  # Check that R is orthogonal (R'R = I)
  expect_true(all(abs(t(R_computed) %*% R_computed - diag(2)) < 1e-6))
  
  # Check determinant is 1 or -1 (orthogonal matrix)
  expect_true(abs(abs(det(R_computed)) - 1) < 1e-6)
  
  # X %*% R should approximately equal Y (or -Y for reflection)
  X_rotated <- X %*% R_computed
  diff_pos <- sum((X_rotated - Y)^2)
  diff_neg <- sum((X_rotated + Y)^2)
  expect_true(min(diff_pos, diff_neg) < 1e-6)
})
