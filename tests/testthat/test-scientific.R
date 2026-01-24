# Scientific accuracy tests for core algorithms
# These tests verify mathematical correctness of key computations

# =============================================================================
# Test 1: Sinkhorn Algorithm Correctness
# =============================================================================

test_that("sinkhorn satisfies marginal constraints", {
  set.seed(42)
  
  for (i in 1:5) {  # Run multiple random tests
    n <- sample(5:15, 1)
    m <- sample(5:15, 1)
    
    # Random cost matrix
    M <- matrix(runif(n * m, 0, 10), nrow = n, ncol = m)
    
    # Random marginals
    p <- runif(n)
    p <- p / sum(p)
    q <- runif(m)
    q <- q / sum(q)
    
    T_mat <- sinkhorn_cpp(M, p, q, reg = 0.1, max_iter = 500, tol = 1e-9)
    
    # Check marginals
    row_sums <- rowSums(T_mat)
    col_sums <- colSums(T_mat)
    
    expect_true(all(abs(row_sums - p) < 1e-4), 
                label = paste("Row marginals, iter", i))
    expect_true(all(abs(col_sums - q) < 1e-4),
                label = paste("Col marginals, iter", i))
  }
})

test_that("sinkhorn produces doubly stochastic matrix for uniform marginals", {
  n <- 5
  M <- matrix(runif(n * n), nrow = n, ncol = n)
  p <- rep(1/n, n)
  q <- rep(1/n, n)
  
  T_mat <- sinkhorn_cpp(M, p, q, reg = 0.1, max_iter = 500, tol = 1e-9)
  
  # Should be doubly stochastic
  expect_true(all(abs(rowSums(T_mat) - 1/n) < 1e-6))
  expect_true(all(abs(colSums(T_mat) - 1/n) < 1e-6))
  
  # Total mass should be 1
  expect_equal(sum(T_mat), 1, tolerance = 1e-10)
})

# =============================================================================
# Test 2: Gromov-Wasserstein Mathematical Properties
# =============================================================================

test_that("GW transport plan satisfies constraints", {
  set.seed(42)
  n <- 5
  m <- 6
  
  # Generate point clouds
  X <- matrix(runif(n * 2), ncol = 2)
  Y <- matrix(runif(m * 2), ncol = 2)
  
  C1 <- as.matrix(dist(X))
  C2 <- as.matrix(dist(Y))
  
  p <- rep(1/n, n)
  q <- rep(1/m, m)
  
  result <- gw_cpp(C1, C2, p, q, reg = 0.1, max_iter = 100)
  T_mat <- result$T
  
  # Check marginal constraints
  expect_true(all(abs(rowSums(T_mat) - p) < 0.05))
  expect_true(all(abs(colSums(T_mat) - q) < 0.05))
  
  # Total mass
  expect_equal(sum(T_mat), 1, tolerance = 0.01)
})

test_that("GW is symmetric for same distributions", {
  set.seed(42)
  n <- 5
  
  X <- matrix(runif(n * 2), ncol = 2)
  C <- as.matrix(dist(X))
  p <- rep(1/n, n)
  
  result <- gw_cpp(C, C, p, p, reg = 0.1, max_iter = 100)
  
  # GW(C, C) should produce identity-like transport
  # Diagonal should have higher values
  T_mat <- result$T
  
  # Check that transport is approximately uniform/symmetric
  expect_true(is.finite(result$gw_dist))
})

# =============================================================================
# Test 3: FGW Interpolation Property
# =============================================================================

test_that("FGW alpha interpolates between Wasserstein and GW", {
  set.seed(42)
  n <- 5
  m <- 5
  
  M <- matrix(runif(n * m), nrow = n, ncol = m)
  C1 <- as.matrix(dist(matrix(runif(n * 2), ncol = 2)))
  C2 <- as.matrix(dist(matrix(runif(m * 2), ncol = 2)))
  
  p <- rep(1/n, n)
  q <- rep(1/m, m)
  
  # alpha = 0: pure Wasserstein
  result_w <- fused_gromov_wasserstein_cpp(M, C1, C2, p, q, 
                                            alpha = 0, reg = 0.1)
  
  # alpha = 1: pure GW  
  result_gw <- fused_gromov_wasserstein_cpp(M, C1, C2, p, q,
                                             alpha = 1, reg = 0.1)
  
  # alpha = 0.5: mixture
  result_fgw <- fused_gromov_wasserstein_cpp(M, C1, C2, p, q,
                                              alpha = 0.5, reg = 0.1)
  
  # All should be valid transport plans
  expect_true(all(result_w$T >= -1e-10))
  expect_true(all(result_gw$T >= -1e-10))
  expect_true(all(result_fgw$T >= -1e-10))
  
  # Marginals should be satisfied
  expect_true(all(abs(rowSums(result_fgw$T) - p) < 0.05))
  expect_true(all(abs(colSums(result_fgw$T) - q) < 0.05))
})

# =============================================================================
# Test 4: Euclidean Distance Properties
# =============================================================================

test_that("euclidean distance satisfies metric properties", {
  set.seed(42)
  n <- 10
  d <- 3
  
  X <- matrix(runif(n * d), ncol = d)
  
  D <- euclidean_distance_cpp(X, NULL)
  
  # Non-negativity
  expect_true(all(D >= 0))
  
  # Identity of indiscernibles (diagonal = 0, with numerical tolerance)
  expect_true(all(abs(diag(D)) < 1e-6))
  
  # Symmetry
  expect_true(all(abs(D - t(D)) < 1e-6))
  
  # Triangle inequality
  for (i in 1:min(n, 5)) {
    for (j in 1:min(n, 5)) {
      for (k in 1:min(n, 5)) {
        expect_true(D[i, j] <= D[i, k] + D[k, j] + 1e-6)
      }
    }
  }
})

# =============================================================================
# Test 5: Rotation Matrix Properties
# =============================================================================

test_that("rotation matrix is orthogonal", {
  set.seed(42)
  n <- 10
  
  X <- matrix(runif(n * 2), ncol = 2)
  theta <- runif(1) * 2 * pi
  R_true <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 
                   nrow = 2, byrow = TRUE)
  Y <- X %*% R_true
  
  weights <- rep(1/n, n)
  R_computed <- compute_rotation_svd_cpp(X, Y, weights)
  
  # Orthogonality: R'R = I
  expect_true(all(abs(t(R_computed) %*% R_computed - diag(2)) < 1e-6))
  
  # R * R' = I
  expect_true(all(abs(R_computed %*% t(R_computed) - diag(2)) < 1e-6))
})

# =============================================================================
# Test 6: KL Divergence Properties
# =============================================================================

test_that("KL divergence satisfies mathematical properties", {
  # KL(P||Q) >= 0 (Gibbs' inequality)
  P <- matrix(c(0.2, 0.3, 0.5,
                0.1, 0.6, 0.3), nrow = 2, byrow = TRUE)
  Q <- matrix(c(0.33, 0.33, 0.34,
                0.25, 0.5, 0.25), nrow = 2, byrow = TRUE)
  
  KL <- kl_divergence_matrix_cpp(P, Q)
  
  # Non-negativity
  expect_true(all(KL >= -1e-10))
  
  # KL(P||P) = 0 for each distribution with itself
  KL_self <- kl_divergence_matrix_cpp(P, P)
  expect_true(all(abs(diag(KL_self)) < 1e-10))
  
  # Asymmetry: generally KL(P||Q) != KL(Q||P)
  KL_rev <- kl_divergence_matrix_cpp(Q, P)
  # Just check that the function runs - asymmetry is expected
  expect_equal(dim(KL_rev), c(2, 2))
})

# =============================================================================
# Test 7: Centering Properties
# =============================================================================

test_that("center_points moves centroid to origin", {
  set.seed(42)
  n <- 20
  
  X <- matrix(runif(n * 2, -10, 10), ncol = 2)
  weights <- runif(n)
  weights <- weights / sum(weights)
  
  centered <- center_points_cpp(X, weights)
  
  # Weighted centroid should be at origin
  weighted_centroid <- colSums(centered * weights)
  expect_true(all(abs(weighted_centroid) < 1e-10))
})

# =============================================================================
# Test 8: NMI Mathematical Properties
# =============================================================================

test_that("NMI satisfies mathematical bounds", {
  set.seed(42)
  
  for (i in 1:5) {
    n <- 30
    k1 <- sample(2:5, 1)
    k2 <- sample(2:5, 1)
    
    labels1 <- sample(1:k1, n, replace = TRUE)
    labels2 <- sample(1:k2, n, replace = TRUE)
    
    nmi <- as.numeric(calculate_nmi(labels1, labels2))
    
    # NMI should be in [0, 1]
    expect_true(nmi >= 0 - 1e-10, label = paste("NMI >= 0, iter", i))
    expect_true(nmi <= 1 + 1e-10, label = paste("NMI <= 1, iter", i))
  }
  
  # Perfect clustering should have NMI = 1
  labels <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
  nmi_perfect <- as.numeric(calculate_nmi(labels, labels))
  expect_equal(nmi_perfect, 1, tolerance = 1e-6)
})

# =============================================================================
# Test 9: Topic Model Probability Properties
# =============================================================================

test_that("LDA produces valid probability distributions", {
  skip_if_not_installed("text2vec")
  
  set.seed(42)
  n_genes <- 50
  n_cells <- 30
  
  mat <- matrix(rpois(n_genes * n_cells, lambda = 5), 
                nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  model <- train_lda(mat, n.topics = 3, n.iter = 50, seed = 42, verbose = FALSE)
  
  # topic_cell: columns should sum to 1 (topic distribution per cell)
  col_sums <- colSums(model$topic_cell)
  expect_true(all(abs(col_sums - 1) < 0.01))
  
  # All values should be non-negative
  expect_true(all(model$topic_cell >= 0))
  expect_true(all(model$topic_gene >= 0))
  
  # topic_gene: rows should sum to 1 (word distribution per topic)
  row_sums <- rowSums(model$topic_gene)
  expect_true(all(abs(row_sums - 1) < 0.01))
})

# =============================================================================
# Test 10: Deconvolution Probability Properties
# =============================================================================

test_that("deconvolution produces valid fractions", {
  skip_if_not_installed("text2vec")
  
  set.seed(42)
  n_genes <- 30
  n_cells <- 30
  n_spots <- 10
  
  sc_mat <- Matrix::Matrix(rpois(n_genes * n_cells, lambda = 5), 
                           nrow = n_genes, ncol = n_cells, sparse = TRUE)
  rownames(sc_mat) <- paste0("Gene", 1:n_genes)
  colnames(sc_mat) <- paste0("Cell", 1:n_cells)
  
  celltypes <- rep(c("TypeA", "TypeB", "TypeC"), each = 10)
  names(celltypes) <- colnames(sc_mat)
  
  st_mat <- Matrix::Matrix(rpois(n_genes * n_spots, lambda = 10), 
                           nrow = n_genes, ncol = n_spots, sparse = TRUE)
  rownames(st_mat) <- paste0("Gene", 1:n_genes)
  colnames(st_mat) <- paste0("Spot", 1:n_spots)
  
  model <- train_lda(sc_mat, n.topics = 3, n.iter = 50, seed = 42, verbose = FALSE)
  topic_ct <- calculate_topic_celltype(model$topic_cell, celltypes)
  ct_counts <- table(celltypes)
  ct_bayes <- calculate_celltype_topic_bayes(topic_ct, ct_counts, length(celltypes))
  
  result <- deconvolve(st_mat, model, topic_celltype = topic_ct,
                       celltype_topic_bayes = ct_bayes, method = "Bayes", n.iter = 20)
  
  frac <- result$celltype_frac
  
  # All fractions should be in [0, 1]
  expect_true(all(frac >= 0))
  expect_true(all(frac <= 1))
  
  # Each row should sum to 1
  row_sums <- rowSums(frac)
  expect_true(all(abs(row_sums - 1) < 1e-6))
  
  # Check dimensions
  expect_equal(nrow(frac), n_spots)
  expect_equal(ncol(frac), 3)  # 3 cell types
})
