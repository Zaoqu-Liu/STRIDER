# End-to-end workflow tests - complete pipeline validation

# =============================================================================
# Test 1: Complete deconvolution workflow
# =============================================================================

test_that("complete deconvolution workflow produces valid results", {
  skip_if_not_installed("text2vec")
  
  set.seed(42)
  
  # Create realistic scRNA-seq data
  n_genes <- 100
  n_cells <- 90  # 30 cells per type
  n_spots <- 30
  n_celltypes <- 3
  
  # scRNA-seq matrix with cell type specific patterns
  sc_mat <- matrix(rpois(n_genes * n_cells, lambda = 3), 
                   nrow = n_genes, ncol = n_cells)
  
  # Add cell type specific markers
  cells_per_type <- 30
  for (ct in 1:n_celltypes) {
    start_cell <- (ct - 1) * cells_per_type + 1
    end_cell <- ct * cells_per_type
    start_gene <- (ct - 1) * 10 + 1
    end_gene <- ct * 10
    sc_mat[start_gene:end_gene, start_cell:end_cell] <- 
      sc_mat[start_gene:end_gene, start_cell:end_cell] + rpois(10 * cells_per_type, lambda = 15)
  }
  
  rownames(sc_mat) <- paste0("Gene", 1:n_genes)
  colnames(sc_mat) <- paste0("Cell", 1:n_cells)
  sc_mat <- Matrix::Matrix(sc_mat, sparse = TRUE)
  
  # Cell type annotations
  celltypes <- rep(paste0("Type", 1:n_celltypes), each = cells_per_type)
  names(celltypes) <- colnames(sc_mat)
  
  # Spatial data (mixture of cell types)
  st_mat <- matrix(rpois(n_genes * n_spots, lambda = 8), 
                   nrow = n_genes, ncol = n_spots)
  rownames(st_mat) <- paste0("Gene", 1:n_genes)
  colnames(st_mat) <- paste0("Spot", 1:n_spots)
  st_mat <- Matrix::Matrix(st_mat, sparse = TRUE)
  
  # Workflow steps
  # 1. Filter and normalize
  sc_filtered <- filter_genes(sc_mat, min.cells = 3, min.counts = 1)
  sc_norm <- normalize_counts(sc_filtered, scale.factor = 10000)
  
  st_filtered <- filter_genes(st_mat, min.cells = 1, min.counts = 1)
  st_norm <- normalize_counts(st_filtered, scale.factor = 10000)
  
  # 2. Find shared genes
  shared_genes <- get_shared_genes(rownames(sc_norm), rownames(st_norm))
  expect_true(length(shared_genes) > 30)
  
  # 3. Subset to shared genes
  sc_shared <- sc_norm[shared_genes, ]
  st_shared <- st_norm[shared_genes, ]
  
  # 4. Train LDA model
  model <- train_lda(sc_shared, n.topics = n_celltypes + 2, 
                     n.iter = 50, seed = 42, verbose = FALSE)
  
  expect_s3_class(model, "strider_lda")
  expect_equal(model$n_topics, n_celltypes + 2)
  
  # 5. Calculate topic-celltype associations
  topic_ct <- calculate_topic_celltype(model$topic_cell, celltypes)
  expect_equal(nrow(topic_ct), n_celltypes + 2)
  expect_equal(ncol(topic_ct), n_celltypes)
  
  # 6. Calculate Bayesian associations
  ct_counts <- table(celltypes)
  ct_bayes <- calculate_celltype_topic_bayes(topic_ct, ct_counts, length(celltypes))
  expect_equal(nrow(ct_bayes), n_celltypes)
  
  # 7. Deconvolve
  result <- deconvolve(st_shared, model, 
                       topic_celltype = topic_ct,
                       celltype_topic_bayes = ct_bayes,
                       method = "Bayes", n.iter = 20)
  
  # Validate results
  expect_true(is.list(result))
  frac <- result$celltype_frac
  
  # All fractions in [0, 1]
  expect_true(all(frac >= 0))
  expect_true(all(frac <= 1))
  
  # Rows sum to 1
  expect_true(all(abs(rowSums(frac) - 1) < 1e-6))
  
  # Correct dimensions
  expect_equal(nrow(frac), n_spots)
  expect_equal(ncol(frac), n_celltypes)
})

# =============================================================================
# Test 2: Spatial clustering workflow
# =============================================================================

test_that("spatial clustering workflow produces meaningful clusters", {
  set.seed(42)
  n_spots <- 50
  n_celltypes <- 4
  
  # Create spatially organized deconvolution results
  # Top half dominated by Type1/2, bottom half by Type3/4
  deconv <- matrix(0, nrow = n_spots, ncol = n_celltypes)
  
  coords <- data.frame(
    x = rep(1:10, 5),
    y = rep(1:5, each = 10),
    row.names = paste0("Spot", 1:n_spots)
  )
  
  for (i in 1:n_spots) {
    if (coords$y[i] <= 2.5) {
      # Top region
      deconv[i, ] <- c(0.4 + runif(1, -0.1, 0.1), 
                       0.4 + runif(1, -0.1, 0.1),
                       0.1 + runif(1, -0.05, 0.05), 
                       0.1 + runif(1, -0.05, 0.05))
    } else {
      # Bottom region
      deconv[i, ] <- c(0.1 + runif(1, -0.05, 0.05), 
                       0.1 + runif(1, -0.05, 0.05),
                       0.4 + runif(1, -0.1, 0.1), 
                       0.4 + runif(1, -0.1, 0.1))
    }
    deconv[i, ] <- deconv[i, ] / sum(deconv[i, ])
  }
  
  rownames(deconv) <- paste0("Spot", 1:n_spots)
  colnames(deconv) <- paste0("Type", 1:n_celltypes)
  
  # Cluster
  result <- spatial_cluster(deconv, coords, n.clusters = 2, k = 4)
  
  expect_equal(nrow(result), n_spots)
  expect_true("cluster" %in% colnames(result))
  
  # Should identify ~2 main regions
  n_clusters <- length(unique(result$cluster))
  expect_equal(n_clusters, 2)
})

# =============================================================================
# Test 3: Spot-cell mapping workflow
# =============================================================================

test_that("spot-cell mapping identifies similar cells", {
  skip_if_not_installed("text2vec")
  
  set.seed(42)
  n_topics <- 5
  n_spots <- 15
  n_cells <- 50
  
  # Create topic distributions
  # Some spots should clearly match certain cells
  topic_spot <- matrix(0, nrow = n_topics, ncol = n_spots)
  topic_cell <- matrix(0, nrow = n_topics, ncol = n_cells)
  
  # First 5 spots similar to first 10 cells (Type A pattern)
  for (i in 1:5) {
    topic_spot[, i] <- c(0.7, 0.1, 0.1, 0.05, 0.05) + rnorm(n_topics, sd = 0.02)
  }
  for (i in 1:10) {
    topic_cell[, i] <- c(0.7, 0.1, 0.1, 0.05, 0.05) + rnorm(n_topics, sd = 0.02)
  }
  
  # Next 5 spots similar to cells 11-25 (Type B pattern)
  for (i in 6:10) {
    topic_spot[, i] <- c(0.1, 0.7, 0.1, 0.05, 0.05) + rnorm(n_topics, sd = 0.02)
  }
  for (i in 11:25) {
    topic_cell[, i] <- c(0.1, 0.7, 0.1, 0.05, 0.05) + rnorm(n_topics, sd = 0.02)
  }
  
  # Last 5 spots similar to cells 26-50 (Type C pattern)
  for (i in 11:15) {
    topic_spot[, i] <- c(0.1, 0.1, 0.7, 0.05, 0.05) + rnorm(n_topics, sd = 0.02)
  }
  for (i in 26:50) {
    topic_cell[, i] <- c(0.1, 0.1, 0.7, 0.05, 0.05) + rnorm(n_topics, sd = 0.02)
  }
  
  # Ensure valid probabilities
  topic_spot[topic_spot < 0] <- 0.01
  topic_cell[topic_cell < 0] <- 0.01
  topic_spot <- topic_spot / rep(colSums(topic_spot), each = n_topics)
  topic_cell <- topic_cell / rep(colSums(topic_cell), each = n_topics)
  
  rownames(topic_spot) <- paste0("Topic", 1:n_topics)
  colnames(topic_spot) <- paste0("Spot", 1:n_spots)
  rownames(topic_cell) <- paste0("Topic", 1:n_topics)
  colnames(topic_cell) <- paste0("Cell", 1:n_cells)
  
  # Map spots to cells
  mapping <- map_spots_to_cells(topic_spot, topic_cell, n.top = 5, method = "euclidean")
  
  expect_equal(nrow(mapping), n_spots * 5)
  
  # Check that Spot1's top matches are from Cell1-10 (mostly)
  spot1_cells <- mapping[mapping$spot == "Spot1" & mapping$rank <= 3, "cell"]
  spot1_cell_nums <- as.numeric(gsub("Cell", "", spot1_cells))
  expect_true(sum(spot1_cell_nums <= 10) >= 2)
  
  # Cell type annotations
  celltypes <- c(rep("TypeA", 10), rep("TypeB", 15), rep("TypeC", 25))
  names(celltypes) <- paste0("Cell", 1:n_cells)
  
  # Get mapped cell types
  mapped_ct <- get_mapped_celltypes(mapping, celltypes, aggregate = "majority")
  
  expect_equal(nrow(mapped_ct), n_spots)
  
  # First 5 spots should be TypeA
  spot1_5_ct <- mapped_ct[mapped_ct$spot %in% paste0("Spot", 1:5), "predicted_celltype"]
  expect_true(sum(spot1_5_ct == "TypeA") >= 3)
})

# =============================================================================
# Test 4: Multi-sample integration workflow
# =============================================================================

test_that("multi-sample integration aligns samples correctly", {
  set.seed(42)
  n_samples <- 3
  n_spots_per_sample <- 20
  n_topics <- 4
  n_celltypes <- 3
  
  # Create sample data
  deconv_list <- list()
  topic_list <- list()
  coords_list <- list()
  
  for (s in 1:n_samples) {
    # Deconvolution results
    deconv <- matrix(runif(n_spots_per_sample * n_celltypes), 
                     nrow = n_spots_per_sample)
    deconv <- deconv / rowSums(deconv)
    rownames(deconv) <- paste0("S", s, "_Spot", 1:n_spots_per_sample)
    colnames(deconv) <- paste0("CT", 1:n_celltypes)
    deconv_list[[s]] <- deconv
    
    # Topic distributions
    topic <- matrix(runif(n_topics * n_spots_per_sample), nrow = n_topics)
    topic <- topic / rep(colSums(topic), each = n_topics)
    rownames(topic) <- paste0("Topic", 1:n_topics)
    colnames(topic) <- paste0("S", s, "_Spot", 1:n_spots_per_sample)
    topic_list[[s]] <- topic
    
    # Coordinates - offset by sample
    coords <- data.frame(
      x = runif(n_spots_per_sample) * 100 + (s - 1) * 200,
      y = runif(n_spots_per_sample) * 100,
      row.names = paste0("S", s, "_Spot", 1:n_spots_per_sample)
    )
    coords_list[[s]] <- coords
  }
  
  # Integrate
  result <- integrate_samples(
    deconv_list, topic_list, coords_list,
    sample.ids = paste0("Sample", 1:n_samples),
    alpha = 0.5, reg = 0.1, max.iter = 30, seed = 42
  )
  
  expect_s3_class(result, "strider_integration")
  
  # Check combined coordinates
  expect_equal(nrow(result$aligned_coords), n_spots_per_sample * n_samples)
  expect_true("x_aligned" %in% colnames(result$aligned_coords))
  expect_true("sample" %in% colnames(result$aligned_coords))
  
  # Check transport plans
  expect_equal(length(result$transport_plans), n_samples - 1)
  
  # Check combined deconv
  expect_equal(nrow(result$deconv), n_spots_per_sample * n_samples)
  
  # Quality metrics
  quality <- integration_quality(result)
  expect_true(is.data.frame(quality))
  expect_true(all(quality$entropy >= 0))
})

# =============================================================================
# Test 5: Simple pipeline test
# =============================================================================

test_that("basic STRIDER pipeline works", {
  skip_if_not_installed("text2vec")
  
  set.seed(42)
  
  # Simulate data
  n_genes <- 50
  n_cells <- 60
  n_spots <- 20
  n_celltypes <- 3
  
  # scRNA-seq data
  sc_mat <- matrix(rpois(n_genes * n_cells, lambda = 5), 
                   nrow = n_genes, ncol = n_cells)
  
  # Add marker genes
  cells_per_type <- 20
  for (ct in 1:n_celltypes) {
    cell_range <- ((ct - 1) * cells_per_type + 1):(ct * cells_per_type)
    gene_range <- ((ct - 1) * 5 + 1):(ct * 5)
    sc_mat[gene_range, cell_range] <- sc_mat[gene_range, cell_range] + 
                                       rpois(5 * cells_per_type, lambda = 15)
  }
  
  rownames(sc_mat) <- paste0("Gene", 1:n_genes)
  colnames(sc_mat) <- paste0("Cell", 1:n_cells)
  sc_mat <- Matrix::Matrix(sc_mat, sparse = TRUE)
  
  celltypes <- rep(paste0("Type", 1:n_celltypes), each = cells_per_type)
  names(celltypes) <- colnames(sc_mat)
  
  # Spatial data
  st_mat <- matrix(rpois(n_genes * n_spots, lambda = 10), 
                   nrow = n_genes, ncol = n_spots)
  rownames(st_mat) <- paste0("Gene", 1:n_genes)
  colnames(st_mat) <- paste0("Spot", 1:n_spots)
  st_mat <- Matrix::Matrix(st_mat, sparse = TRUE)
  
  coords <- data.frame(
    x = rep(1:5, 4),
    y = rep(1:4, each = 5),
    row.names = paste0("Spot", 1:n_spots)
  )
  
  # Train model with associations
  model <- train_model(sc_mat, celltypes, n.topics = 5, n.iter = 40, seed = 42)
  
  expect_s3_class(model, "strider_lda")
  expect_true("topic_celltype" %in% names(model))
  
  # Deconvolve
  result <- deconvolve(st_mat, model, 
                       topic_celltype = model$topic_celltype,
                       celltype_topic_bayes = model$celltype_topic_bayes,
                       method = "Bayes", n.iter = 20)
  
  frac <- result$celltype_frac
  
  # Validate
  expect_true(all(frac >= 0))
  expect_true(all(frac <= 1))
  expect_true(all(abs(rowSums(frac) - 1) < 1e-6))
  
  # Spatial clustering
  clusters <- spatial_cluster(frac, coords, n.clusters = 3, k = 4)
  
  expect_equal(nrow(clusters), n_spots)
  expect_true("cluster" %in% colnames(clusters))
  
  # Map to cells
  mapping <- map_spots_to_cells(result$topic_spot, model$topic_cell, 
                                n.top = 5, method = "euclidean")
  
  expect_true(nrow(mapping) > 0)
  
  # Pipeline completed successfully
  expect_true(TRUE)
})
