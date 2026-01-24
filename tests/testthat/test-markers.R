# Tests for marker gene identification functions

test_that("find_markers identifies differentially expressed genes", {
  set.seed(42)
  n_genes <- 50
  n_cells <- 60
  
  # Create matrix with clear marker genes
  mat <- matrix(rpois(n_genes * n_cells, lambda = 2), 
                nrow = n_genes, ncol = n_cells)
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  
  # Make some genes specific to each cell type
  # Genes 1-5: markers for TypeA
  mat[1:5, 1:20] <- mat[1:5, 1:20] + rpois(5 * 20, lambda = 15)
  # Genes 6-10: markers for TypeB
  mat[6:10, 21:40] <- mat[6:10, 21:40] + rpois(5 * 20, lambda = 15)
  # Genes 11-15: markers for TypeC
  mat[11:15, 41:60] <- mat[11:15, 41:60] + rpois(5 * 20, lambda = 15)
  
  mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  celltypes <- c(rep("TypeA", 20), rep("TypeB", 20), rep("TypeC", 20))
  names(celltypes) <- colnames(mat)
  
  result <- find_markers(mat, celltypes, n.top = 10, min.pct = 0.1,
                         logfc.threshold = 0.1, test.use = "wilcox")
  
  # Check structure
  expect_true(is.list(result))
  expect_true("markers" %in% names(result))
  expect_true("top_markers" %in% names(result))
  
  # Should find some markers
  expect_true(nrow(result$markers) > 0)
  expect_true(length(result$top_markers) > 0)
  
  # Check that expected marker genes are found
  typeA_markers <- result$markers[result$markers$cluster == "TypeA", ]
  expect_true(any(typeA_markers$gene %in% paste0("Gene", 1:5)))
})

test_that("find_markers returns correct columns", {
  set.seed(42)
  
  mat <- Matrix::Matrix(rpois(100 * 30, lambda = 3), 
                        nrow = 100, ncol = 30, sparse = TRUE)
  rownames(mat) <- paste0("Gene", 1:100)
  colnames(mat) <- paste0("Cell", 1:30)
  
  celltypes <- rep(c("TypeA", "TypeB", "TypeC"), each = 10)
  names(celltypes) <- colnames(mat)
  
  result <- find_markers(mat, celltypes, n.top = 5, min.pct = 0.05,
                         logfc.threshold = 0.1)
  
  expected_cols <- c("gene", "cluster", "avg_logfc", "pct.1", "pct.2", 
                     "pct_diff", "p_val", "p_val_adj")
  expect_true(all(expected_cols %in% colnames(result$markers)))
})

test_that("find_variable_genes identifies highly variable genes", {
  set.seed(42)
  n_genes <- 100
  n_cells <- 50
  
  # Create matrix with varying levels of variance
  mat <- matrix(0, nrow = n_genes, ncol = n_cells)
  
  # Low variance genes (1-30) - constant expression
  mat[1:30, ] <- matrix(rpois(30 * n_cells, lambda = 10), nrow = 30)
  
  # High variance genes (31-50) - very clear bimodal expression
  for (i in 31:50) {
    mat[i, 1:25] <- rpois(25, lambda = 1)
    mat[i, 26:50] <- rpois(25, lambda = 50)
  }
  
  # Medium variance genes (51-100)
  mat[51:100, ] <- matrix(rpois(50 * n_cells, lambda = 5), nrow = 50)
  
  rownames(mat) <- paste0("Gene", 1:n_genes)
  colnames(mat) <- paste0("Cell", 1:n_cells)
  mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  # Use wider mean range to include high-variance genes
  result <- find_variable_genes(mat, n.top = 20, min.mean = 0.01, max.mean = 100)
  
  # Check structure
  expect_true(is.list(result))
  expect_true("hvg" %in% names(result))
  expect_true("stats" %in% names(result))
  
  # Should return requested number of HVGs (or less if not enough genes pass filters)
  expect_true(length(result$hvg) <= 20)
  
  # High variance genes (31-50) should be well-represented (at least 2)
  n_highvar_in_hvg <- sum(result$hvg %in% paste0("Gene", 31:50))
  expect_true(n_highvar_in_hvg >= 2)
})

test_that("get_shared_genes finds intersection correctly", {
  sc_genes <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")
  st_genes <- c("Gene2", "Gene3", "Gene4", "Gene6", "Gene7")
  
  # Without gene_use filter
  shared <- get_shared_genes(sc_genes, st_genes)
  expect_equal(sort(shared), c("Gene2", "Gene3", "Gene4"))
  
  # With gene_use filter
  gene_use <- c("Gene2", "Gene3", "Gene5", "Gene6")
  shared_filtered <- get_shared_genes(sc_genes, st_genes, gene_use)
  expect_equal(sort(shared_filtered), c("Gene2", "Gene3"))
  
  # With "All"
  shared_all <- get_shared_genes(sc_genes, st_genes, "All")
  expect_equal(sort(shared_all), c("Gene2", "Gene3", "Gene4"))
})
