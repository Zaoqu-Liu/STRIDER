# Tests for utility functions

test_that("get_default_colors returns correct length", {
  colors <- get_default_colors()
  expect_true(length(colors) >= 30)
  expect_true(all(grepl("^#", colors)))
})

test_that("get_extended_colors works correctly", {
  colors_10 <- get_extended_colors(10)
  colors_50 <- get_extended_colors(50)
  
  expect_equal(length(colors_10), 10)
  expect_equal(length(colors_50), 50)
})

test_that("normalize_to_one works for vectors", {
  x <- c(1, 2, 3, 4)
  result <- normalize_to_one(x)
  
  expect_equal(sum(result), 1)
  expect_equal(length(result), 4)
})

test_that("normalize_to_one works for matrices", {
  mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)
  
  result_row <- normalize_to_one(mat, margin = 1)
  result_col <- normalize_to_one(mat, margin = 2)
  
  expect_equal(nrow(result_row), 2)
  expect_equal(ncol(result_row), 3)
})

test_that("safe_divide handles zeros", {
  result <- safe_divide(c(1, 2, 3), c(1, 0, 3))
  
  expect_equal(result[1], 1)
  expect_equal(result[2], 0)  # Zero division returns 0
  expect_equal(result[3], 1)
})
