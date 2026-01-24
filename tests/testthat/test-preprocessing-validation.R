# Preprocessing and Marker Finding Validation Tests
# Deep tests for data preprocessing pipeline

# =============================================================================
# Test 1: Normalization Correctness
# =============================================================================

test_that("normalize_counts produces correct scale factors", {
  set.seed(42)
  n_genes <- 100
  n_cells <- 50
  
  # Create count matrix with known total counts
  mat <- matrix(rpois(n_genes * n_cells, lambda = 10), 
                nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  
  # Known total counts per cell
  orig_totals <- colSums(mat)
  
  # Normalize to scale factor 10000
  norm_mat <- normalize_counts(Matrix::Matrix(mat, sparse = TRUE), 
                               scale.factor = 10000)
  
  # After normalization, all columns should sum to 10000
  norm_totals <- Matrix::colSums(norm_mat)
  expect_true(all(abs(norm_totals - 10000) < 1e-6))
})

test_that("normalization preserves zero entries", {
  set.seed(42)
  n_genes <- 50
  n_cells <- 30
  
  # Create sparse matrix with ~70% zeros
  mat <- matrix(0, nrow = n_genes, ncol = n_cells)
  non_zero <- sample(length(mat), size = floor(0.3 * length(mat)))
  mat[non_zero] <- rpois(length(non_zero), lambda = 5)
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  
  sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  # Normalize
  norm_mat <- normalize_counts(sparse_mat, scale.factor = 10000)
  
  # Zero entries should remain zero
  expect_equal(sum(mat == 0), sum(norm_mat == 0))
})

test_that("normalization with log transform is monotonic", {
  set.seed(42)
  n_genes <- 30
  n_cells <- 20
  
  mat <- matrix(rpois(n_genes * n_cells, lambda = 10), 
                nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  
  sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  # Normalize without log
  norm_mat <- normalize_counts(sparse_mat, scale.factor = 10000, log.transform = FALSE)
  
  # Normalize with log
  norm_log_mat <- normalize_counts(sparse_mat, scale.factor = 10000, log.transform = TRUE)
  
  # log1p is monotonic, so if a > b in norm_mat, then log1p(a) > log1p(b)
  for (i in 1:10) {
    idx1 <- sample(length(norm_mat), 1)
    idx2 <- sample(length(norm_mat), 1)
    if (norm_mat[idx1] > norm_mat[idx2]) {
      expect_true(norm_log_mat[idx1] > norm_log_mat[idx2])
    }
  }
})

# =============================================================================
# Test 2: Scaling Validation
# =============================================================================

test_that("scale_data produces unit standard deviation", {
  set.seed(42)
  n_genes <- 50
  n_cells <- 100
  
  mat <- matrix(rnorm(n_genes * n_cells, mean = 10, sd = 3), 
                nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  
  scaled <- scale_data(mat, center = FALSE, scale = TRUE)
  
  # Each row should have POPULATION SD = 1 (N, not N-1 to match sklearn)
  row_sds <- apply(scaled, 1, function(x) sqrt(mean((x - mean(x))^2)))
  expect_true(all(abs(row_sds - 1) < 1e-6))
})

test_that("scale_data with centering produces zero mean", {
  set.seed(42)
  n_genes <- 30
  n_cells <- 50
  
  mat <- matrix(rnorm(n_genes * n_cells, mean = 10, sd = 3), 
                nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  
  scaled <- scale_data(mat, center = TRUE, scale = TRUE)
  
  # Each row should have mean = 0
  row_means <- rowMeans(scaled)
  expect_true(all(abs(row_means) < 1e-10))
  
  # Each row should have POPULATION SD = 1 (N, not N-1 to match sklearn)
  row_sds <- apply(scaled, 1, function(x) sqrt(mean((x - mean(x))^2)))
  expect_true(all(abs(row_sds - 1) < 1e-6))
})

# =============================================================================
# Test 3: Filtering Validation
# =============================================================================

test_that("filter_genes correctly applies thresholds", {
  set.seed(42)
  n_genes <- 100
  n_cells <- 50
  
  # Create matrix where first 20 genes are expressed in < 5 cells
  mat <- matrix(0, nrow = n_genes, ncol = n_cells)
  
  # Low expression genes (expressed in 3 cells each)
  for (i in 1:20) {
    cells <- sample(n_cells, 3)
    mat[i, cells] <- rpois(3, lambda = 5)
  }
  
  # High expression genes (expressed in 30 cells each)
  for (i in 21:n_genes) {
    cells <- sample(n_cells, 30)
    mat[i, cells] <- rpois(30, lambda = 10)
  }
  
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  
  sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  # Filter with min.cells = 5
  filtered <- filter_genes(sparse_mat, min.cells = 5)
  
  # Should remove the first 20 low-expression genes
  expect_true(nrow(filtered) <= 80)
  
  # All kept genes should be expressed in >= 5 cells
  cells_per_gene <- Matrix::rowSums(filtered > 0)
  expect_true(all(cells_per_gene >= 5))
})

test_that("filter_cells correctly applies thresholds", {
  set.seed(42)
  n_genes <- 100
  n_cells <- 50
  
  # Create matrix where first 10 cells have very few genes detected
  mat <- matrix(rpois(n_genes * n_cells, lambda = 5), 
                nrow = n_genes, ncol = n_cells)
  
  # Make first 10 cells have < 50 genes detected
  for (j in 1:10) {
    zeros <- sample(n_genes, 60)
    mat[zeros, j] <- 0
  }
  
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  
  sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  # Filter with min.genes = 50
  filtered <- filter_cells(sparse_mat, min.genes = 50)
  
  # Should remove most of the first 10 cells
  expect_true(ncol(filtered) <= 45)
  
  # All kept cells should have >= 50 genes
  genes_per_cell <- Matrix::colSums(filtered > 0)
  expect_true(all(genes_per_cell >= 50))
})

# =============================================================================
# Test 4: Marker Finding Statistical Validity
# =============================================================================

test_that("find_markers identifies known differentially expressed genes", {
  set.seed(42)
  n_genes <- 50
  n_cells <- 100
  
  # Create two cell types
  celltypes <- c(rep("TypeA", 50), rep("TypeB", 50))
  names(celltypes) <- paste0("Cell", 1:n_cells)
  
  # Base expression
  mat <- matrix(rnorm(n_genes * n_cells, mean = 5, sd = 1), 
                nrow = n_genes, ncol = n_cells)
  mat[mat < 0] <- 0  # No negative values
  
  # Make genes 1-10 highly expressed in TypeA
  mat[1:10, 1:50] <- mat[1:10, 1:50] + 3
  
  # Make genes 11-20 highly expressed in TypeB
  mat[11:20, 51:100] <- mat[11:20, 51:100] + 3
  
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- names(celltypes)
  
  sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  # Find markers
  result <- find_markers(sparse_mat, celltypes, n.top = 15, 
                         min.pct = 0.1, logfc.threshold = 0.1)
  
  # Check TypeA markers
  typeA_markers <- result$markers[result$markers$cluster == "TypeA", ]
  expect_true(any(paste0("Gene", 1:10) %in% typeA_markers$gene))
  
  # Check TypeB markers
  typeB_markers <- result$markers[result$markers$cluster == "TypeB", ]
  expect_true(any(paste0("Gene", 11:20) %in% typeB_markers$gene))
})

test_that("find_markers returns valid p-values", {
  set.seed(42)
  n_genes <- 30
  n_cells <- 60
  
  celltypes <- c(rep("A", 30), rep("B", 30))
  names(celltypes) <- paste0("Cell", 1:n_cells)
  
  mat <- matrix(rpois(n_genes * n_cells, lambda = 10), 
                nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- names(celltypes)
  
  # Add differential expression
  mat[1:5, 1:30] <- mat[1:5, 1:30] + 10
  
  sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  result <- find_markers(sparse_mat, celltypes, n.top = 10, 
                         min.pct = 0.01, logfc.threshold = 0.1)
  
  # All p-values should be in [0, 1]
  expect_true(all(result$markers$p_val >= 0))
  expect_true(all(result$markers$p_val <= 1))
  
  # Adjusted p-values should also be in [0, 1]
  expect_true(all(result$markers$p_val_adj >= 0))
  expect_true(all(result$markers$p_val_adj <= 1))
  
  # Adjusted p-values should be >= raw p-values
  expect_true(all(result$markers$p_val_adj >= result$markers$p_val - 1e-10))
})

# =============================================================================
# Test 5: QC Metrics Calculation
# =============================================================================

test_that("calculate_qc_metrics computes correct values", {
  set.seed(42)
  n_genes <- 50
  n_cells <- 20
  
  mat <- matrix(rpois(n_genes * n_cells, lambda = 5), 
                nrow = n_genes, ncol = n_cells)
  
  # Add mitochondrial genes
  rownames(mat) <- c(paste0("MT-", 1:5), paste0("Gene", 6:n_genes))
  colnames(mat) <- paste0("Cell", 1:n_cells)
  
  sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  qc <- calculate_qc_metrics(sparse_mat)
  
  # Manual calculation
  expected_counts <- as.numeric(colSums(mat))
  expected_genes <- as.numeric(colSums(mat > 0))
  expected_mito <- as.numeric(colSums(mat[1:5, ]) / colSums(mat) * 100)
  
  expect_equal(qc$nCount, expected_counts)
  expect_equal(qc$nFeature, expected_genes)
  expect_equal(qc$pct_mito, expected_mito, tolerance = 1e-10)
})

# =============================================================================
# Test 6: Variable Gene Selection
# =============================================================================

test_that("find_variable_genes selects high-variance genes", {
  set.seed(42)
  n_genes <- 200
  n_cells <- 100
  
  # Create matrix with different variance levels
  mat <- matrix(0, nrow = n_genes, ncol = n_cells)
  
  # Low variance genes (first 150)
  for (i in 1:150) {
    mat[i, ] <- rnorm(n_cells, mean = 5, sd = 0.5)
  }
  
  # High variance genes (last 50)
  for (i in 151:200) {
    mat[i, ] <- rnorm(n_cells, mean = 5, sd = 3)
  }
  
  mat[mat < 0] <- 0
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  
  sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  # Normalize first (required for HVG selection)
  norm_mat <- normalize_counts(sparse_mat, scale.factor = 10000, log.transform = TRUE)
  
  result <- find_variable_genes(norm_mat, n.top = 50, min.mean = 0.01, max.mean = 10)
  
  # High variance genes (151-200) should be preferentially selected
  high_var_selected <- sum(paste0("Gene", 151:200) %in% result$hvg)
  expect_true(high_var_selected >= 20)  # At least 40% of high-var genes selected
})

# =============================================================================
# Test 7: Preprocessing Pipeline Integration
# =============================================================================

test_that("preprocess_sc produces consistent output", {
  set.seed(42)
  n_genes <- 100
  n_cells <- 80
  
  mat <- matrix(rpois(n_genes * n_cells, lambda = 10), 
                nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  
  celltypes <- rep(c("A", "B", "C", "D"), each = 20)
  names(celltypes) <- colnames(mat)
  
  sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  result <- preprocess_sc(sparse_mat, celltypes, 
                          min.cells = 5, min.genes = 10,
                          normalize = TRUE, scale = FALSE)
  
  # Check output structure
  expect_true("scale_matrix" %in% names(result))
  expect_true("raw_matrix" %in% names(result))
  expect_true("genes" %in% names(result))
  expect_true("cells" %in% names(result))
  expect_true("celltypes" %in% names(result))
  
  # Dimensions should match
  expect_equal(nrow(result$scale_matrix), length(result$genes))
  expect_equal(ncol(result$scale_matrix), length(result$cells))
  expect_equal(length(result$celltypes), length(result$cells))
  
  # Cell types should be aligned with cells
  expect_true(all(names(result$celltypes) == result$cells))
})

# =============================================================================
# Test 8: Sparse Matrix Handling
# =============================================================================

test_that("preprocessing handles sparse matrices correctly", {
  set.seed(42)
  n_genes <- 50
  n_cells <- 40
  
  # Create very sparse matrix (~90% zeros)
  mat <- matrix(0, nrow = n_genes, ncol = n_cells)
  non_zero <- sample(length(mat), size = floor(0.1 * length(mat)))
  mat[non_zero] <- rpois(length(non_zero), lambda = 20)
  
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  
  celltypes <- rep(c("X", "Y"), each = 20)
  names(celltypes) <- colnames(mat)
  
  # Test with sparse matrix
  sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  result_sparse <- preprocess_sc(sparse_mat, celltypes, 
                                  min.cells = 1, min.genes = 1,
                                  normalize = TRUE, scale = FALSE)
  
  # Test with dense matrix
  result_dense <- preprocess_sc(mat, celltypes,
                                 min.cells = 1, min.genes = 1, 
                                 normalize = TRUE, scale = FALSE)
  
  # Results should be approximately equal
  expect_equal(dim(result_sparse$scale_matrix), dim(result_dense$scale_matrix))
  
  # Values should be close (after converting both to dense for comparison)
  sparse_dense <- as.matrix(result_sparse$scale_matrix)
  dense_result <- as.matrix(result_dense$scale_matrix)
  
  # Check non-zero elements are similar
  expect_true(cor(as.vector(sparse_dense), as.vector(dense_result)) > 0.99)
})

# =============================================================================
# Test 9: Edge Cases
# =============================================================================

test_that("preprocessing handles single cell type", {
  set.seed(42)
  n_genes <- 30
  n_cells <- 20
  
  mat <- matrix(rpois(n_genes * n_cells, lambda = 10), 
                nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  
  # Single cell type
  celltypes <- rep("SingleType", n_cells)
  names(celltypes) <- colnames(mat)
  
  sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  result <- preprocess_sc(sparse_mat, celltypes, 
                          min.cells = 1, min.genes = 1,
                          normalize = TRUE)
  
  expect_equal(length(unique(result$celltypes)), 1)
})

test_that("normalization handles zero-sum columns gracefully", {
  set.seed(42)
  n_genes <- 20
  n_cells <- 10
  
  mat <- matrix(rpois(n_genes * n_cells, lambda = 5), 
                nrow = n_genes, ncol = n_cells)
  
  # Make one column all zeros
  mat[, 5] <- 0
  
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  
  sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  # Should not fail
  norm_mat <- normalize_counts(sparse_mat, scale.factor = 10000)
  
  expect_true(is.finite(sum(norm_mat)))
  expect_false(any(is.nan(as.matrix(norm_mat))))
})

# =============================================================================
# Test 10: Wilcoxon Test Mathematical Validity
# =============================================================================

test_that("Wilcoxon test in find_markers agrees with stats::wilcox.test", {
  set.seed(42)
  n_genes <- 10
  n_cells <- 40
  
  celltypes <- c(rep("A", 20), rep("B", 20))
  names(celltypes) <- paste0("Cell", 1:n_cells)
  
  mat <- matrix(rpois(n_genes * n_cells, lambda = 10), 
                nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- names(celltypes)
  
  # Add clear differential expression to gene 1
  mat[1, 1:20] <- mat[1, 1:20] + 20
  
  sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  # Run find_markers
  result <- find_markers(sparse_mat, celltypes, n.top = 10,
                         min.pct = 0.01, logfc.threshold = 0,
                         only.pos = FALSE, test.use = "wilcox")
  
  # Get p-value for gene 1, cell type A
  gene1_pval <- result$markers[result$markers$gene == "Gene1" & 
                                result$markers$cluster == "A", "p_val"]
  
  # Manual Wilcoxon test
  expr_A <- as.numeric(mat[1, 1:20])
  expr_B <- as.numeric(mat[1, 21:40])
  manual_test <- wilcox.test(expr_A, expr_B, alternative = "two.sided", exact = FALSE)
  
  # P-values should be similar
  expect_true(length(gene1_pval) > 0)
  expect_equal(gene1_pval[1], manual_test$p.value, tolerance = 0.01)
})
