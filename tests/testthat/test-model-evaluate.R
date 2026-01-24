# Tests for model evaluation functions

test_that("calculate_nmi returns valid NMI score", {
  # Perfect match should give NMI = 1
  labels_true <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
  labels_pred <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
  
  nmi <- calculate_nmi(labels_pred, labels_true)
  
  expect_equal(as.numeric(nmi), 1, tolerance = 1e-6)
})

test_that("calculate_nmi handles different label values", {
  # Same clustering but different labels
  labels_true <- c("A", "A", "A", "B", "B", "B")
  labels_pred <- c("X", "X", "X", "Y", "Y", "Y")
  
  nmi <- calculate_nmi(labels_pred, labels_true)
  
  expect_equal(as.numeric(nmi), 1, tolerance = 1e-6)
})

test_that("calculate_nmi is symmetric", {
  set.seed(42)
  labels1 <- sample(1:3, 30, replace = TRUE)
  labels2 <- sample(1:3, 30, replace = TRUE)
  
  nmi_12 <- calculate_nmi(labels1, labels2)
  nmi_21 <- calculate_nmi(labels2, labels1)
  
  expect_equal(nmi_12, nmi_21, tolerance = 1e-10)
})

test_that("calculate_nmi is between 0 and 1", {
  set.seed(42)
  
  for (i in 1:5) {
    labels1 <- sample(1:4, 50, replace = TRUE)
    labels2 <- sample(1:4, 50, replace = TRUE)
    
    nmi <- calculate_nmi(labels1, labels2)
    
    expect_true(nmi >= 0)
    expect_true(nmi <= 1)
  }
})

test_that("evaluate_model_prediction calculates accuracy correctly", {
  skip_if_not_installed("text2vec")
  
  set.seed(42)
  n_topics <- 3
  n_cells <- 30
  n_celltypes <- 3
  
  # Create topic-cell distribution
  topic_cell <- matrix(0, nrow = n_topics, ncol = n_cells)
  rownames(topic_cell) <- paste0("Topic", 1:n_topics)
  colnames(topic_cell) <- paste0("Cell", 1:n_cells)
  
  # Make topic 1 dominant for cells 1-10 (TypeA)
  topic_cell[1, 1:10] <- runif(10, 0.6, 0.8)
  topic_cell[2, 1:10] <- runif(10, 0.1, 0.2)
  topic_cell[3, 1:10] <- 1 - topic_cell[1, 1:10] - topic_cell[2, 1:10]
  
  # Topic 2 dominant for cells 11-20 (TypeB)
  topic_cell[2, 11:20] <- runif(10, 0.6, 0.8)
  topic_cell[1, 11:20] <- runif(10, 0.1, 0.2)
  topic_cell[3, 11:20] <- 1 - topic_cell[1, 11:20] - topic_cell[2, 11:20]
  
  # Topic 3 dominant for cells 21-30 (TypeC)
  topic_cell[3, 21:30] <- runif(10, 0.6, 0.8)
  topic_cell[1, 21:30] <- runif(10, 0.1, 0.2)
  topic_cell[2, 21:30] <- 1 - topic_cell[1, 21:30] - topic_cell[3, 21:30]
  
  # Topic-celltype association
  topic_ct <- matrix(
    c(0.8, 0.1, 0.1,
      0.1, 0.8, 0.1,
      0.1, 0.1, 0.8),
    nrow = 3, byrow = TRUE
  )
  rownames(topic_ct) <- paste0("Topic", 1:3)
  colnames(topic_ct) <- c("TypeA", "TypeB", "TypeC")
  
  # True cell types
  true_ct <- c(rep("TypeA", 10), rep("TypeB", 10), rep("TypeC", 10))
  names(true_ct) <- colnames(topic_cell)
  
  result <- evaluate_model_prediction(topic_cell, topic_ct, true_ct, method = "Raw")
  
  # Check structure
  expect_true(is.list(result))
  expect_true("accuracy" %in% names(result))
  expect_true("predictions" %in% names(result))
  
  # Accuracy should be reasonable with this setup
  expect_true(result$accuracy > 0.5)
  expect_true(result$accuracy <= 1)
})

test_that("cluster_in_topic_space returns valid clusters", {
  skip_if_not_installed("text2vec")
  
  set.seed(42)
  n_topics <- 3
  n_cells <- 30
  
  # Create topic-cell distribution with clear clusters
  topic_cell <- matrix(0, nrow = n_topics, ncol = n_cells)
  rownames(topic_cell) <- paste0("Topic", 1:n_topics)
  colnames(topic_cell) <- paste0("Cell", 1:n_cells)
  
  # Clear separation
  topic_cell[1, 1:10] <- 0.8
  topic_cell[2, 1:10] <- 0.1
  topic_cell[3, 1:10] <- 0.1
  
  topic_cell[1, 11:20] <- 0.1
  topic_cell[2, 11:20] <- 0.8
  topic_cell[3, 11:20] <- 0.1
  
  topic_cell[1, 21:30] <- 0.1
  topic_cell[2, 21:30] <- 0.1
  topic_cell[3, 21:30] <- 0.8
  
  # True cell types
  true_ct <- c(rep("TypeA", 10), rep("TypeB", 10), rep("TypeC", 10))
  names(true_ct) <- colnames(topic_cell)
  
  result <- cluster_in_topic_space(topic_cell, true_ct, n_clusters = 3, seed = 42)
  
  expect_true(is.list(result))
  expect_true("clusters" %in% names(result))
  expect_true("nmi" %in% names(result))
  
  # With clear separation, NMI should be high
  expect_true(result$nmi > 0.7)
  
  # Cluster assignments should have correct length
  expect_equal(length(result$clusters), n_cells)
})
