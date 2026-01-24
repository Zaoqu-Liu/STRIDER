# Tests for deconvolution functions

test_that("deconvolve produces valid cell type fractions", {
  skip_if_not_installed("text2vec")
  
  set.seed(42)
  n_genes <- 50
  n_cells <- 30
  n_spots <- 20
  
  # Create training data
  sc_mat <- matrix(rpois(n_genes * n_cells, lambda = 5), 
                   nrow = n_genes, ncol = n_cells)
  rownames(sc_mat) <- paste0("Gene", 1:n_genes)
  colnames(sc_mat) <- paste0("Cell", 1:n_cells)
  sc_mat <- Matrix::Matrix(sc_mat, sparse = TRUE)
  
  # Create cell type annotations
  celltypes <- rep(c("TypeA", "TypeB", "TypeC"), each = 10)
  names(celltypes) <- colnames(sc_mat)
  
  # Create spatial data (same genes)
  st_mat <- matrix(rpois(n_genes * n_spots, lambda = 10), 
                   nrow = n_genes, ncol = n_spots)
  rownames(st_mat) <- paste0("Gene", 1:n_genes)
  colnames(st_mat) <- paste0("Spot", 1:n_spots)
  st_mat <- Matrix::Matrix(st_mat, sparse = TRUE)
  
  # Train model
  model <- train_lda(sc_mat, n.topics = 3, n.iter = 50, seed = 42, verbose = FALSE)
  
  # Calculate associations
  topic_ct <- calculate_topic_celltype(model$topic_cell, celltypes)
  
  ct_counts <- table(celltypes)
  ct_bayes <- calculate_celltype_topic_bayes(topic_ct, ct_counts, length(celltypes))
  
  # Deconvolve
  result <- deconvolve(st_mat, model, topic_celltype = topic_ct,
                       celltype_topic_bayes = ct_bayes, method = "Bayes", n.iter = 20)
  
  # Check structure
  expect_true(is.list(result))
  expect_true("celltype_frac" %in% names(result))
  expect_true("topic_spot" %in% names(result))
  
  # Check cell type fractions
  frac <- result$celltype_frac
  expect_equal(nrow(frac), n_spots)
  expect_equal(ncol(frac), 3)  # 3 cell types
  
  # Fractions should sum to 1
  row_sums <- rowSums(frac)
  expect_true(all(abs(row_sums - 1) < 1e-6))
  
  # All fractions should be non-negative
  expect_true(all(frac >= 0))
  expect_true(all(frac <= 1))
})

test_that("deconvolve works with different methods", {
  skip_if_not_installed("text2vec")
  
  set.seed(42)
  n_genes <- 30
  n_cells <- 20
  n_spots <- 10
  
  sc_mat <- Matrix::Matrix(rpois(n_genes * n_cells, lambda = 3), 
                           nrow = n_genes, ncol = n_cells, sparse = TRUE)
  rownames(sc_mat) <- paste0("Gene", 1:n_genes)
  colnames(sc_mat) <- paste0("Cell", 1:n_cells)
  
  celltypes <- rep(c("TypeA", "TypeB"), each = 10)
  names(celltypes) <- colnames(sc_mat)
  
  st_mat <- Matrix::Matrix(rpois(n_genes * n_spots, lambda = 5), 
                           nrow = n_genes, ncol = n_spots, sparse = TRUE)
  rownames(st_mat) <- paste0("Gene", 1:n_genes)
  colnames(st_mat) <- paste0("Spot", 1:n_spots)
  
  model <- train_lda(sc_mat, n.topics = 2, n.iter = 30, seed = 42, verbose = FALSE)
  topic_ct <- calculate_topic_celltype(model$topic_cell, celltypes)
  ct_counts <- table(celltypes)
  ct_bayes <- calculate_celltype_topic_bayes(topic_ct, ct_counts, length(celltypes))
  
  methods <- c("Raw", "Norm", "NormBySD", "Bayes", "BayesNorm")
  
  for (method in methods) {
    result <- deconvolve(st_mat, model, topic_celltype = topic_ct,
                         celltype_topic_bayes = ct_bayes, method = method, n.iter = 20)
    
    frac <- result$celltype_frac
    
    # Basic validity checks
    expect_equal(nrow(frac), n_spots, 
                 label = paste("Method:", method, "- rows"))
    expect_true(all(frac >= 0), 
                label = paste("Method:", method, "- non-negative"))
    expect_true(all(abs(rowSums(frac) - 1) < 1e-6), 
                label = paste("Method:", method, "- sum to 1"))
  }
})

test_that("train_model combines training and association calculation", {
  skip_if_not_installed("text2vec")
  
  set.seed(42)
  n_genes <- 30
  n_cells <- 30
  
  sc_mat <- Matrix::Matrix(rpois(n_genes * n_cells, lambda = 3), 
                           nrow = n_genes, ncol = n_cells, sparse = TRUE)
  rownames(sc_mat) <- paste0("Gene", 1:n_genes)
  colnames(sc_mat) <- paste0("Cell", 1:n_cells)
  
  celltypes <- rep(c("TypeA", "TypeB", "TypeC"), each = 10)
  names(celltypes) <- colnames(sc_mat)
  
  model <- train_model(sc_mat, celltypes, n.topics = 3, n.iter = 30, seed = 42)
  
  # Check that additional matrices are computed
  expect_true("topic_celltype" %in% names(model))
  expect_true("celltype_topic_bayes" %in% names(model))
  
  expect_equal(nrow(model$topic_celltype), 3)
  expect_equal(ncol(model$topic_celltype), 3)
})

test_that("create_strider_object creates valid object", {
  deconv <- matrix(c(0.5, 0.3, 0.2, 0.2, 0.6, 0.2), nrow = 2, byrow = TRUE)
  rownames(deconv) <- c("Spot1", "Spot2")
  colnames(deconv) <- c("TypeA", "TypeB", "TypeC")
  
  obj <- create_strider_object(
    deconv = deconv,
    n_topics = 5,
    method = "Bayes",
    celltypes = c("TypeA", "TypeB", "TypeC")
  )
  
  expect_s3_class(obj, "STRIDER")
  expect_equal(nrow(obj$deconv), 2)
  expect_equal(obj$n_topics, 5)
  expect_equal(obj$method, "Bayes")
})
