# Tests for LDA model functions

test_that("train_lda creates valid model structure", {
  skip_if_not_installed("text2vec")
  
  # Create test matrix with realistic structure
  set.seed(42)
  n_genes <- 100
  n_cells <- 50
  
  mat <- matrix(0, nrow = n_genes, ncol = n_cells)
  # Create some pattern - different gene sets active in different cells
  mat[1:30, 1:25] <- rpois(30 * 25, lambda = 5)
  mat[31:60, 26:50] <- rpois(30 * 25, lambda = 5)
  mat[61:100, ] <- rpois(40 * 50, lambda = 2)
  
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  # Train LDA with small number of iterations for speed
  model <- train_lda(mat, n.topics = 3, n.iter = 50, seed = 42, verbose = FALSE)
  
  # Check structure
  expect_s3_class(model, "strider_lda")
  expect_equal(model$n_topics, 3)
  expect_equal(length(model$vocab), n_genes)
  expect_equal(nrow(model$topic_gene), 3)
  expect_equal(ncol(model$topic_gene), n_genes)
  expect_equal(nrow(model$topic_cell), 3)
  expect_equal(ncol(model$topic_cell), n_cells)
  
  # Check that probabilities sum to ~1 per column
  topic_sums <- colSums(model$topic_cell)
  expect_true(all(abs(topic_sums - 1) < 0.01))
})

test_that("calculate_topic_celltype produces valid associations", {
  skip_if_not_installed("text2vec")
  
  set.seed(42)
  n_genes <- 50
  n_cells <- 30
  
  # Create structured matrix
  mat <- matrix(rpois(n_genes * n_cells, lambda = 3), nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  # Cell type annotations
  celltypes <- rep(c("TypeA", "TypeB", "TypeC"), each = 10)
  names(celltypes) <- colnames(mat)
  
  # Train model
  model <- train_lda(mat, n.topics = 3, n.iter = 50, seed = 42, verbose = FALSE)
  
  # Calculate topic-celltype association
  topic_ct <- calculate_topic_celltype(model$topic_cell, celltypes)
  
  # Check dimensions
  expect_equal(nrow(topic_ct), 3)  # n_topics
  expect_equal(ncol(topic_ct), 3)  # n_celltypes
  expect_equal(colnames(topic_ct), sort(unique(celltypes)))
  
  # Values should be non-negative
  expect_true(all(topic_ct >= 0))
})

test_that("calculate_celltype_topic_bayes applies Bayes correctly", {
  # Create test topic-celltype matrix
  topic_ct <- matrix(
    c(0.6, 0.2, 0.2,
      0.3, 0.5, 0.2,
      0.1, 0.3, 0.6),
    nrow = 3, ncol = 3, byrow = TRUE
  )
  rownames(topic_ct) <- paste0("Topic", 1:3)
  colnames(topic_ct) <- c("TypeA", "TypeB", "TypeC")
  
  # Cell type counts (unequal)
  ct_counts <- c(TypeA = 100, TypeB = 50, TypeC = 150)
  total <- sum(ct_counts)
  
  bayes_ct <- calculate_celltype_topic_bayes(topic_ct, ct_counts, total)
  
  # Check dimensions
  expect_equal(nrow(bayes_ct), 3)  # n_celltypes
  expect_equal(ncol(bayes_ct), 3)  # n_topics
  
  # Values should be non-negative
  expect_true(all(bayes_ct >= 0))
  expect_true(all(is.finite(bayes_ct)))
})

test_that("transform_lda produces valid topic distributions", {
  skip_if_not_installed("text2vec")
  
  set.seed(42)
  n_genes <- 50
  n_train_cells <- 30
  n_test_cells <- 10
  
  # Training matrix
  train_mat <- matrix(rpois(n_genes * n_train_cells, lambda = 3), 
                      nrow = n_genes, ncol = n_train_cells)
  rownames(train_mat) <- paste0("Gene", 1:n_genes)
  colnames(train_mat) <- paste0("TrainCell", 1:n_train_cells)
  train_mat <- Matrix::Matrix(train_mat, sparse = TRUE)
  
  # Test matrix (same genes)
  test_mat <- matrix(rpois(n_genes * n_test_cells, lambda = 3), 
                     nrow = n_genes, ncol = n_test_cells)
  rownames(test_mat) <- paste0("Gene", 1:n_genes)
  colnames(test_mat) <- paste0("TestCell", 1:n_test_cells)
  test_mat <- Matrix::Matrix(test_mat, sparse = TRUE)
  
  # Train model
  model <- train_lda(train_mat, n.topics = 3, n.iter = 50, seed = 42, verbose = FALSE)
  
  # Transform test data
  topic_test <- transform_lda(model, test_mat, n.iter = 20)
  
  # Check dimensions
  expect_equal(nrow(topic_test), 3)  # n_topics
  expect_equal(ncol(topic_test), n_test_cells)
  
  # Probabilities should sum to ~1
  topic_sums <- colSums(topic_test)
  expect_true(all(abs(topic_sums - 1) < 0.01))
})

test_that("calculate_coherence returns valid scores", {
  skip_if_not_installed("text2vec")
  
  set.seed(42)
  n_genes <- 100
  n_cells <- 50
  
  mat <- matrix(rpois(n_genes * n_cells, lambda = 3), 
                nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  model <- train_lda(mat, n.topics = 3, n.iter = 50, seed = 42, verbose = FALSE)
  
  coherence <- calculate_coherence(model, mat, n.words = 5)
  
  # Check structure
  expect_true(is.list(coherence))
  expect_true("mean" %in% names(coherence))
  expect_true("per_topic" %in% names(coherence))
  expect_equal(length(coherence$per_topic), 3)
  
  # UMass coherence is typically negative
  expect_true(is.finite(coherence$mean))
})
