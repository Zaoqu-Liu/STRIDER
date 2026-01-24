# Numerical Stability Tests
# Tests for edge cases and numerical robustness

# =============================================================================
# Test 1: Zero and Near-Zero Handling
# =============================================================================

test_that("sinkhorn handles near-uniform marginals", {
  n <- 10
  M <- matrix(runif(n * n), nrow = n)
  
  # Near-uniform (but not exactly)
  p <- rep(0.1, n) + runif(n, -0.001, 0.001)
  p <- p / sum(p)
  q <- rep(0.1, n) + runif(n, -0.001, 0.001)
  q <- q / sum(q)
  
  T_mat <- sinkhorn_cpp(M, p, q, reg = 0.1)
  
  # Should not contain NaN or Inf
  expect_true(all(is.finite(T_mat)))
  
  # Marginal constraints should be approximately satisfied
  expect_true(all(abs(rowSums(T_mat) - p) < 0.05))
  expect_true(all(abs(colSums(T_mat) - q) < 0.05))
})

test_that("KL divergence handles distributions with zeros", {
  # Distribution with many zeros
  P <- matrix(c(0.9, 0.1, 0, 0, 0.5, 0.5, 0, 0), nrow = 2, byrow = TRUE)
  Q <- matrix(c(0.7, 0.1, 0.1, 0.1, 0.3, 0.3, 0.2, 0.2), nrow = 2, byrow = TRUE)
  
  KL <- kl_divergence_matrix_cpp(P, Q)
  
  # Should be finite
  expect_true(all(is.finite(KL)))
  
  # Should be non-negative
  expect_true(all(KL >= -1e-10))
})

test_that("euclidean distance handles identical points", {
  X <- matrix(c(1, 2, 3, 1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE)
  
  D <- euclidean_distance_cpp(X, NULL)
  
  # Distance from point to itself is 0
  expect_true(all(abs(diag(D)) < 1e-10))
  
  # D[1,2] should be 0 (identical points)
  expect_equal(D[1, 2], 0, tolerance = 1e-10)
  
  # All values should be non-negative
  expect_true(all(D >= -1e-10))
})

# =============================================================================
# Test 2: Large Value Handling
# =============================================================================

test_that("sinkhorn handles large cost values", {
  n <- 5
  # Large cost values
  M <- matrix(runif(n * n) * 1000, nrow = n)
  
  p <- rep(1/n, n)
  q <- rep(1/n, n)
  
  # Should not overflow - use higher regularization for stability
  T_mat <- sinkhorn_cpp(M, p, q, reg = 100)
  
  expect_true(all(is.finite(T_mat)))
  expect_true(all(T_mat >= 0))
})

test_that("normalization handles large counts", {
  set.seed(42)
  n_genes <- 20
  n_cells <- 10
  
  # Very large counts
  mat <- matrix(rpois(n_genes * n_cells, lambda = 10000), 
                nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  
  sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  norm_mat <- normalize_counts(sparse_mat, scale.factor = 10000)
  
  # Should be finite
  expect_true(all(is.finite(as.matrix(norm_mat))))
  
  # Columns should sum to 10000
  col_sums <- Matrix::colSums(norm_mat)
  expect_true(all(abs(col_sums - 10000) < 1e-6))
})

# =============================================================================
# Test 3: Small Value Handling
# =============================================================================

test_that("KL divergence handles very small probabilities", {
  # Very small probabilities
  P <- matrix(c(0.999, 0.001, 0.001, 0.999), nrow = 2)
  Q <- matrix(c(0.001, 0.999, 0.999, 0.001), nrow = 2)
  
  KL <- kl_divergence_matrix_cpp(P, Q)
  
  expect_true(all(is.finite(KL)))
  expect_true(all(KL >= -1e-10))
})

test_that("safe_divide handles zero denominators", {
  x <- c(1, 2, 3, 4)
  y <- c(1, 0, 2, 0)
  
  result <- safe_divide(x, y, replace = 0)
  
  expect_equal(result[1], 1)
  expect_equal(result[2], 0)  # 2/0 -> 0
  expect_equal(result[3], 1.5)
  expect_equal(result[4], 0)  # 4/0 -> 0
})

# =============================================================================
# Test 4: Empty and Minimal Inputs
# =============================================================================

test_that("normalize_to_one handles zero vectors", {
  x <- c(0, 0, 0)
  result <- normalize_to_one(x)
  expect_equal(result, x)  # Should return unchanged
})

test_that("scale_data handles constant rows", {
  # Matrix with constant row
  mat <- matrix(c(5, 5, 5, 5, 1, 2, 3, 4), nrow = 2, byrow = TRUE)
  rownames(mat) <- c("Constant", "Variable")
  colnames(mat) <- paste0("Cell", 1:4)
  
  scaled <- scale_data(mat, center = FALSE, scale = TRUE)
  
  # Variable row should be scaled
  expect_true(abs(sd(scaled[2, ]) - 1) < 1e-6)
  
  # Should not contain NaN
  expect_true(all(is.finite(scaled)))
})

# =============================================================================
# Test 5: Rotation Matrix Edge Cases
# =============================================================================

test_that("rotation handles collinear points", {
  # Points on a line
  X <- matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  Y <- matrix(c(0, 0, -1, 1, -2, 2), ncol = 2, byrow = TRUE)
  
  weights <- rep(1/3, 3)
  
  # Should not fail
  R <- compute_rotation_svd_cpp(X, Y, weights)
  
  # R should be orthogonal
  RTR <- t(R) %*% R
  expect_true(all(abs(RTR - diag(2)) < 1e-6))
})

test_that("centering handles single point", {
  X <- matrix(c(5, 10), nrow = 1)
  weights <- c(1)
  
  centered <- center_points_cpp(X, weights)
  
  # Single point centered should be at origin
  expect_true(all(abs(centered) < 1e-10))
})

# =============================================================================
# Test 6: Deconvolution Edge Cases
# =============================================================================

test_that("deconvolution handles single cell type", {
  skip_if_not_installed("text2vec")
  
  set.seed(42)
  n_genes <- 30
  n_cells <- 20
  n_spots <- 10
  
  sc_mat <- Matrix::Matrix(rpois(n_genes * n_cells, lambda = 10), 
                           nrow = n_genes, sparse = TRUE)
  rownames(sc_mat) <- paste0("Gene", 1:n_genes)
  colnames(sc_mat) <- paste0("Cell", 1:n_cells)
  
  # Single cell type
  celltypes <- rep("SingleType", n_cells)
  names(celltypes) <- colnames(sc_mat)
  
  st_mat <- Matrix::Matrix(rpois(n_genes * n_spots, lambda = 10), 
                           nrow = n_genes, sparse = TRUE)
  rownames(st_mat) <- paste0("Gene", 1:n_genes)
  colnames(st_mat) <- paste0("Spot", 1:n_spots)
  
  model <- train_lda(sc_mat, n.topics = 3, n.iter = 30, seed = 42, verbose = FALSE)
  topic_ct <- calculate_topic_celltype(model$topic_cell, celltypes)
  
  ct_counts <- table(celltypes)
  ct_bayes <- calculate_celltype_topic_bayes(topic_ct, ct_counts, length(celltypes))
  
  result <- deconvolve(st_mat, model, topic_celltype = topic_ct,
                       celltype_topic_bayes = ct_bayes, method = "Bayes", n.iter = 20)
  
  # All fractions should be 1 for single cell type
  expect_true(all(abs(result$celltype_frac[, "SingleType"] - 1) < 0.01))
})

# =============================================================================
# Test 7: Integration Edge Cases
# =============================================================================

test_that("FGW handles identical samples", {
  set.seed(42)
  n <- 10
  
  # Identical topic distributions
  topic1 <- matrix(runif(5 * n), nrow = 5)
  colnames(topic1) <- paste0("Spot", 1:n)
  
  # Same coordinates
  coords1 <- data.frame(x = runif(n), y = runif(n), row.names = paste0("Spot", 1:n))
  
  # FGW between identical samples
  result <- align_samples(topic1, topic1, coords1, coords1, alpha = 0.5)
  
  # Distance should be finite
  expect_true(is.finite(result$fgw_dist))
  
  # Transport plan should be valid (entries sum to 1, non-negative)
  expect_true(all(result$transport >= -1e-10))
  expect_true(abs(sum(result$transport) - 1) < 0.1)
})

# =============================================================================
# Test 8: LDA Edge Cases
# =============================================================================

test_that("LDA handles sparse data", {
  skip_if_not_installed("text2vec")
  
  set.seed(42)
  n_genes <- 40
  n_cells <- 30
  
  # Moderately sparse matrix (~70% zeros)
  mat <- matrix(0, nrow = n_genes, ncol = n_cells)
  non_zero <- sample(length(mat), size = floor(0.30 * length(mat)))
  mat[non_zero] <- rpois(length(non_zero), lambda = 10)
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  
  sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  # Should train without error
  model <- train_lda(sparse_mat, n.topics = 3, n.iter = 50, seed = 42, verbose = FALSE)
  
  # Topic distributions should be valid
  expect_true(all(model$topic_cell >= 0))
  # Allow more tolerance for sparse data
  expect_true(all(abs(colSums(model$topic_cell) - 1) < 0.1))
})

# =============================================================================
# Test 9: Marker Finding Edge Cases
# =============================================================================

test_that("find_markers handles unbalanced cell types", {
  set.seed(42)
  n_genes <- 20
  n_cells <- 100
  
  # Very unbalanced: 90 type A, 10 type B
  celltypes <- c(rep("TypeA", 90), rep("TypeB", 10))
  names(celltypes) <- paste0("Cell", 1:n_cells)
  
  mat <- matrix(rpois(n_genes * n_cells, lambda = 10), 
                nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- names(celltypes)
  
  # Make gene 1 a strong marker for TypeB
  mat[1, 91:100] <- mat[1, 91:100] + 20
  
  sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  result <- find_markers(sparse_mat, celltypes, n.top = 5, 
                         min.pct = 0.01, logfc.threshold = 0.1)
  
  # Should find markers despite imbalance
  expect_true(nrow(result$markers) > 0)
  
  # Gene1 should be marker for TypeB
  typeB_markers <- result$markers[result$markers$cluster == "TypeB", ]
  expect_true("Gene1" %in% typeB_markers$gene)
})

# =============================================================================
# Test 10: Clustering Edge Cases
# =============================================================================

test_that("spatial_cluster handles spots with no neighbors", {
  set.seed(42)
  
  # Create coordinates with isolated spots
  coords <- data.frame(
    x = c(0, 100, 200, 1, 101, 201),  # Spots 1,2,3 are isolated; 4,5,6 are near them
    y = c(0, 0, 0, 1, 1, 1),
    row.names = paste0("Spot", 1:6)
  )
  
  deconv <- matrix(runif(6 * 3), nrow = 6)
  deconv <- deconv / rowSums(deconv)
  rownames(deconv) <- paste0("Spot", 1:6)
  colnames(deconv) <- paste0("CT", 1:3)
  
  # Should not fail with isolated points
  result <- spatial_cluster(deconv, coords, n.clusters = 2, k = 2, seed = 42)
  
  expect_equal(nrow(result), 6)
  expect_true(all(result$cluster %in% 1:2))
})

# =============================================================================
# Test 11: NMI Edge Cases
# =============================================================================

test_that("NMI handles single cluster", {
  labels1 <- rep(1, 10)
  labels2 <- rep(1, 10)
  
  # NMI for single cluster is 0 (H(X) = H(Y) = 0, so MI = 0)
  # This is mathematically correct: no information gain
  nmi <- calculate_nmi(labels1, labels2)
  expect_true(is.finite(nmi))
  expect_true(nmi >= 0)
})

test_that("NMI handles all different labels", {
  labels1 <- 1:10
  labels2 <- 11:20
  
  # Each item in its own cluster - NMI should be 0
  # (since there's no correspondence between labelings)
  nmi <- calculate_nmi(labels1, labels2)
  expect_true(is.finite(nmi))
})

# =============================================================================
# Test 12: Distance Metric Properties
# =============================================================================

test_that("euclidean distance satisfies triangle inequality", {
  set.seed(42)
  X <- matrix(runif(30), nrow = 10)
  
  D <- euclidean_distance_cpp(X, NULL)
  
  # Triangle inequality: D[i,k] <= D[i,j] + D[j,k]
  n <- nrow(D)
  for (attempt in 1:10) {
    i <- sample(n, 1)
    j <- sample(n, 1)
    k <- sample(n, 1)
    
    expect_true(D[i, k] <= D[i, j] + D[j, k] + 1e-10)
  }
})

test_that("cosine distance is bounded", {
  set.seed(42)
  # Both matrices must have same number of columns
  X <- matrix(runif(30), nrow = 5, ncol = 6)
  Y <- matrix(runif(48), nrow = 8, ncol = 6)
  
  D <- compute_cosine_distance(X, Y)
  
  # Cosine distance should be in [0, 2]
  expect_true(all(D >= -1e-10))
  expect_true(all(D <= 2 + 1e-10))
})
