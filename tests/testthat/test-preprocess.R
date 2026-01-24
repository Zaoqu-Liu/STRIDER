# Tests for preprocessing functions

test_that("filter_genes works correctly", {
  # Create test matrix
  mat <- matrix(c(0, 1, 2, 3, 4, 5, 0, 0, 0, 1, 2, 3), nrow = 3)
  rownames(mat) <- c("Gene1", "Gene2", "Gene3")
  colnames(mat) <- c("Cell1", "Cell2", "Cell3", "Cell4")
  mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  # Filter
  filtered <- filter_genes(mat, min.cells = 2, min.counts = 0)
  
  expect_true(nrow(filtered) <= nrow(mat))
})

test_that("normalize_counts preserves dimensions", {
  mat <- matrix(runif(100), nrow = 10)
  rownames(mat) <- paste0("Gene", 1:10)
  colnames(mat) <- paste0("Cell", 1:10)
  mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  normalized <- normalize_counts(mat, scale.factor = 10000)
  
  expect_equal(dim(normalized), dim(mat))
})

test_that("scale_data works correctly", {
  mat <- matrix(runif(100), nrow = 10)
  rownames(mat) <- paste0("Gene", 1:10)
  colnames(mat) <- paste0("Cell", 1:10)
  
  scaled <- scale_data(mat, center = FALSE, scale = TRUE)
  
  expect_equal(dim(scaled), dim(mat))
})

test_that("calculate_qc_metrics returns correct columns", {
  mat <- matrix(rpois(100, lambda = 10), nrow = 10)
  rownames(mat) <- paste0("Gene", 1:10)
  colnames(mat) <- paste0("Cell", 1:10)
  mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  qc <- calculate_qc_metrics(mat)
  
  expect_true("nCount" %in% colnames(qc))
  expect_true("nFeature" %in% colnames(qc))
  expect_equal(nrow(qc), ncol(mat))
})
