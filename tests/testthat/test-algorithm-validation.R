# Algorithm Validation Tests
# These tests verify mathematical correctness of core algorithms against known results

# =============================================================================
# Test 1: Sinkhorn Algorithm - Known Solution Verification
# =============================================================================

test_that("sinkhorn recovers known optimal transport solution", {
  # For uniform marginals and symmetric cost, transport should be uniform
  n <- 4
  M <- matrix(1, nrow = n, ncol = n)
  diag(M) <- 0  # Zero cost for diagonal
  
  p <- rep(1/n, n)
  q <- rep(1/n, n)
  
  T_mat <- sinkhorn_cpp(M, p, q, reg = 0.01, max_iter = 2000, tol = 1e-10)
  
  # Should concentrate mass on diagonal (minimum cost)
  diag_mass <- sum(diag(T_mat))
  expect_true(diag_mass > 0.9)  # Most mass on diagonal
})

test_that("sinkhorn preserves total mass", {
  set.seed(42)
  n <- 10
  m <- 8
  
  M <- matrix(runif(n * m), nrow = n, ncol = m)
  p <- runif(n)
  p <- p / sum(p)
  q <- runif(m)
  q <- q / sum(q)
  
  T_mat <- sinkhorn_cpp(M, p, q, reg = 0.1)
  
  # Total mass should be 1
  expect_equal(sum(T_mat), 1, tolerance = 1e-6)
})

# =============================================================================
# Test 2: GW Gradient - Mathematical Verification
# =============================================================================

test_that("GW gradient has correct structure", {
  # For squared loss, gradient should be -4 * C1 * T * C2
  n <- 3
  m <- 4
  
  C1 <- matrix(c(0, 1, 2, 1, 0, 1, 2, 1, 0), nrow = n)
  C2 <- matrix(c(0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0), nrow = m)
  T_mat <- matrix(1/(n*m), nrow = n, ncol = m)  # Uniform transport
  
  grad <- gwgrad_cpp(C1, C2, T_mat)
  
  # Manual calculation
  expected <- -4 * C1 %*% T_mat %*% C2
  
  expect_equal(as.vector(grad), as.vector(expected), tolerance = 1e-10)
})

# =============================================================================
# Test 3: KL Divergence - Mathematical Properties
# =============================================================================

test_that("KL divergence satisfies fundamental properties", {
  # Property 1: KL(P||P) = 0
  P <- matrix(c(0.5, 0.3, 0.2), nrow = 1)
  KL_self <- kl_divergence_matrix_cpp(P, P)
  expect_equal(KL_self[1, 1], 0, tolerance = 1e-10)
  
  # Property 2: KL(P||Q) >= 0 (Gibbs' inequality)
  P <- matrix(c(0.5, 0.3, 0.2, 0.1, 0.8, 0.1), nrow = 2, byrow = TRUE)
  Q <- matrix(c(0.33, 0.33, 0.34, 0.25, 0.5, 0.25), nrow = 2, byrow = TRUE)
  KL <- kl_divergence_matrix_cpp(P, Q)
  expect_true(all(KL >= -1e-10))
  
  # Property 3: KL is asymmetric
  KL_PQ <- kl_divergence_matrix_cpp(P, Q)
  KL_QP <- kl_divergence_matrix_cpp(Q, P)
  # These should NOT be equal in general
  expect_true(!all(abs(KL_PQ - t(KL_QP)) < 1e-10))
})

test_that("KL divergence equals analytic formula for simple case", {
  # KL(p||q) = sum p_i * log(p_i/q_i)
  p <- c(0.6, 0.4)
  q <- c(0.5, 0.5)
  
  P <- matrix(p, nrow = 1)
  Q <- matrix(q, nrow = 1)
  
  KL_cpp <- kl_divergence_matrix_cpp(P, Q)[1, 1]
  KL_analytic <- sum(p * log(p / q))
  
  expect_equal(KL_cpp, KL_analytic, tolerance = 1e-10)
})

# =============================================================================
# Test 4: Euclidean Distance - Known Values
# =============================================================================

test_that("euclidean distance computes correct values", {
  # Known distances
  X <- matrix(c(0, 0, 3, 4, 1, 0), nrow = 3, byrow = TRUE)
  
  D <- euclidean_distance_cpp(X, NULL)
  
  # Distance from (0,0) to (3,4) = 5 (3-4-5 triangle)
  expect_equal(D[1, 2], 5, tolerance = 1e-10)
  
  # Distance from (0,0) to (1,0) = 1
  expect_equal(D[1, 3], 1, tolerance = 1e-10)
  
  # Distance from (3,4) to (1,0) = sqrt((3-1)^2 + (4-0)^2) = sqrt(20)
  expect_equal(D[2, 3], sqrt(20), tolerance = 1e-10)
})

# =============================================================================
# Test 5: Rotation Matrix - Orthogonality Verification
# =============================================================================

test_that("rotation matrix preserves orthogonality", {
  set.seed(42)
  n <- 20
  
  # Random points
  X <- matrix(rnorm(n * 2), ncol = 2)
  
  # Apply known rotation (45 degrees)
  theta <- pi / 4
  R_true <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 
                   nrow = 2, byrow = TRUE)
  Y <- X %*% R_true
  
  weights <- rep(1/n, n)
  R_est <- compute_rotation_svd_cpp(X, Y, weights)
  
  # R^T * R should be identity
  RTR <- t(R_est) %*% R_est
  expect_true(all(abs(RTR - diag(2)) < 1e-6))
  
  # R * R^T should be identity
  RRT <- R_est %*% t(R_est)
  expect_true(all(abs(RRT - diag(2)) < 1e-6))
})

# =============================================================================
# Test 6: NMI - Known Clustering Results
# =============================================================================

test_that("NMI returns correct values for known cases", {
  # Perfect match: NMI = 1
  labels <- c(1, 1, 2, 2, 3, 3)
  nmi_perfect <- as.numeric(calculate_nmi(labels, labels))
  expect_equal(nmi_perfect, 1, tolerance = 1e-6)
  
  # Random permutation should still give NMI = 1
  # (NMI is invariant to label permutation)
  labels_permuted <- c("A", "A", "B", "B", "C", "C")
  labels_original <- c(1, 1, 2, 2, 3, 3)
  nmi_perm <- as.numeric(calculate_nmi(labels_permuted, labels_original))
  expect_equal(nmi_perm, 1, tolerance = 1e-6)
})

test_that("NMI is bounded between 0 and 1", {
  set.seed(42)
  
  for (i in 1:10) {
    n <- 50
    labels1 <- sample(1:5, n, replace = TRUE)
    labels2 <- sample(1:5, n, replace = TRUE)
    
    nmi <- as.numeric(calculate_nmi(labels1, labels2))
    
    expect_true(nmi >= 0 - 1e-10)
    expect_true(nmi <= 1 + 1e-10)
  }
})

# =============================================================================
# Test 7: Bayes Theorem Application
# =============================================================================

test_that("Bayesian celltype-topic matrix satisfies Bayes theorem", {
  skip_if_not_installed("text2vec")
  
  # Create simple topic-celltype matrix: P(topic|celltype)
  topic_ct <- matrix(c(
    0.7, 0.2, 0.1,  # Topic 1: high for CT1
    0.2, 0.7, 0.1,  # Topic 2: high for CT2
    0.1, 0.1, 0.8   # Topic 3: high for CT3
  ), nrow = 3, byrow = TRUE)
  rownames(topic_ct) <- paste0("Topic", 1:3)
  colnames(topic_ct) <- c("CT1", "CT2", "CT3")
  
  # Unequal celltype counts
  ct_counts <- c(CT1 = 100, CT2 = 200, CT3 = 300)
  total <- sum(ct_counts)
  
  # Calculate Bayesian matrix
  bayes_ct <- calculate_celltype_topic_bayes(topic_ct, ct_counts, total)
  
  # Verify dimensions
  expect_equal(dim(bayes_ct), c(3, 3))
  
  # Verify Bayes: P(ct|topic) = P(topic|ct) * P(ct) / P(topic)
  # P(ct) = ct_counts / total
  p_ct <- ct_counts / total
  
  # P(topic) = sum_ct P(topic|ct) * P(ct)
  p_topic <- as.vector(topic_ct %*% p_ct)
  
  # Manual calculation for CT1, Topic1
  # P(CT1|Topic1) = P(Topic1|CT1) * P(CT1) / P(Topic1)
  expected_ct1_t1 <- topic_ct[1, "CT1"] * p_ct["CT1"] / p_topic[1]
  # Remove names for comparison
  expect_equal(as.numeric(bayes_ct["CT1", 1]), as.numeric(expected_ct1_t1), tolerance = 1e-6)
})

# =============================================================================
# Test 8: LDA Probability Constraints
# =============================================================================

test_that("LDA produces valid probability distributions", {
  skip_if_not_installed("text2vec")
  
  set.seed(42)
  n_genes <- 50
  n_cells <- 40
  n_topics <- 5
  
  # Create count matrix
  mat <- matrix(rpois(n_genes * n_cells, lambda = 5), 
                nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  model <- train_lda(mat, n.topics = n_topics, n.iter = 50, seed = 42, verbose = FALSE)
  
  # topic_cell: each column should sum to 1
  col_sums <- colSums(model$topic_cell)
  expect_true(all(abs(col_sums - 1) < 0.01))
  
  # All values non-negative
  expect_true(all(model$topic_cell >= 0))
  
  # topic_gene: each row should sum to 1
  row_sums <- rowSums(model$topic_gene)
  expect_true(all(abs(row_sums - 1) < 0.01))
  
  # All values non-negative
  expect_true(all(model$topic_gene >= 0))
})

# =============================================================================
# Test 9: Deconvolution Method Consistency
# =============================================================================

test_that("deconvolution methods produce valid fractions", {
  skip_if_not_installed("text2vec")
  
  set.seed(42)
  n_genes <- 40
  n_cells <- 30
  n_spots <- 15
  n_celltypes <- 3
  
  # Create data
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
  
  # Train model
  model <- train_lda(sc_mat, n.topics = 5, n.iter = 50, seed = 42, verbose = FALSE)
  topic_ct <- calculate_topic_celltype(model$topic_cell, celltypes)
  ct_counts <- table(celltypes)
  ct_bayes <- calculate_celltype_topic_bayes(topic_ct, ct_counts, length(celltypes))
  
  # Test all methods
  methods <- c("Raw", "Norm", "NormBySD", "Bayes", "BayesNorm")
  
  for (method in methods) {
    result <- deconvolve(st_mat, model, topic_celltype = topic_ct,
                         celltype_topic_bayes = ct_bayes, method = method, n.iter = 20)
    
    frac <- result$celltype_frac
    
    # All fractions in [0, 1]
    expect_true(all(frac >= 0), label = paste("Method", method, "non-negative"))
    expect_true(all(frac <= 1), label = paste("Method", method, "<= 1"))
    
    # Rows sum to 1
    row_sums <- rowSums(frac)
    expect_true(all(abs(row_sums - 1) < 1e-6), 
                label = paste("Method", method, "sum to 1"))
    
    # Correct dimensions
    expect_equal(nrow(frac), n_spots, label = paste("Method", method, "rows"))
    expect_equal(ncol(frac), n_celltypes, label = paste("Method", method, "cols"))
  }
})

# =============================================================================
# Test 10: FGW Distance Properties
# =============================================================================

test_that("FGW converges and produces valid transport", {
  set.seed(42)
  n <- 8
  m <- 10
  
  # Feature cost
  M <- matrix(runif(n * m), nrow = n, ncol = m)
  
  # Structure costs
  C1 <- as.matrix(dist(matrix(runif(n * 2), ncol = 2)))
  C2 <- as.matrix(dist(matrix(runif(m * 2), ncol = 2)))
  
  p <- rep(1/n, n)
  q <- rep(1/m, m)
  
  result <- fused_gromov_wasserstein_cpp(M, C1, C2, p, q, 
                                          alpha = 0.5, reg = 0.1, 
                                          max_iter = 100, tol = 1e-7)
  
  T_mat <- result$T
  
  # Transport plan is non-negative
  expect_true(all(T_mat >= -1e-10))
  
  # Total mass is 1
  expect_equal(sum(T_mat), 1, tolerance = 0.01)
  
  # Marginal constraints (approximately)
  expect_true(all(abs(rowSums(T_mat) - p) < 0.05))
  expect_true(all(abs(colSums(T_mat) - q) < 0.05))
  
  # Distance is finite
  expect_true(is.finite(result$fgw_dist))
})

# =============================================================================
# Test 11: Topic-Celltype Association Logic
# =============================================================================

test_that("topic-celltype association captures cell type patterns", {
  skip_if_not_installed("text2vec")
  
  set.seed(42)
  n_topics <- 4
  n_cells <- 60
  n_celltypes <- 3
  
  # Create topic-cell matrix where each cell type has dominant topic
  # Use Topic_ prefix to match what calculate_topic_celltype expects
  topic_cell <- matrix(0.1, nrow = n_topics, ncol = n_cells)
  rownames(topic_cell) <- paste0("Topic_", 1:n_topics)
  colnames(topic_cell) <- paste0("Cell", 1:n_cells)
  
  # CT1 (cells 1-20): Topic 1 dominant
  topic_cell[1, 1:20] <- 0.7
  topic_cell[2:4, 1:20] <- 0.1
  
  # CT2 (cells 21-40): Topic 2 dominant
  topic_cell[2, 21:40] <- 0.7
  topic_cell[c(1, 3, 4), 21:40] <- 0.1
  
  # CT3 (cells 41-60): Topic 3 dominant
  topic_cell[3, 41:60] <- 0.7
  topic_cell[c(1, 2, 4), 41:60] <- 0.1
  
  # Normalize columns
  topic_cell <- sweep(topic_cell, 2, colSums(topic_cell), "/")
  
  celltypes <- c(rep("CT1", 20), rep("CT2", 20), rep("CT3", 20))
  names(celltypes) <- colnames(topic_cell)
  
  topic_ct <- calculate_topic_celltype(topic_cell, celltypes)
  
  # Topic 1 should be highest for CT1
  expect_true(topic_ct["Topic_1", "CT1"] > topic_ct["Topic_1", "CT2"])
  expect_true(topic_ct["Topic_1", "CT1"] > topic_ct["Topic_1", "CT3"])
  
  # Topic 2 should be highest for CT2
  expect_true(topic_ct["Topic_2", "CT2"] > topic_ct["Topic_2", "CT1"])
  expect_true(topic_ct["Topic_2", "CT2"] > topic_ct["Topic_2", "CT3"])
  
  # Topic 3 should be highest for CT3
  expect_true(topic_ct["Topic_3", "CT3"] > topic_ct["Topic_3", "CT1"])
  expect_true(topic_ct["Topic_3", "CT3"] > topic_ct["Topic_3", "CT2"])
})

# =============================================================================
# Test 12: Centering Preserves Relative Positions
# =============================================================================

test_that("centering preserves relative distances", {
  set.seed(42)
  n <- 20
  
  coords <- matrix(runif(n * 2, 0, 100), ncol = 2)
  weights <- rep(1/n, n)
  
  centered <- center_points_cpp(coords, weights)
  
  # Compute pairwise distances before and after centering
  D_before <- as.matrix(dist(coords))
  D_after <- as.matrix(dist(centered))
  
  # Distances should be identical
  expect_true(all(abs(D_before - D_after) < 1e-10))
})
