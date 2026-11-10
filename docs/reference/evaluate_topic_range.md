# Evaluate Multiple Topic Numbers

Train and evaluate LDA models with different numbers of topics to find
the optimal configuration.

## Usage

``` r
evaluate_topic_range(
  matrix,
  celltypes,
  topic_range,
  methods = c("Raw", "Norm", "NormBySD", "Bayes", "BayesNorm"),
  n.iter = 500,
  seed = 42,
  save.models = FALSE,
  out.dir = NULL
)
```

## Arguments

- matrix:

  Normalized count matrix (genes x cells)

- celltypes:

  Named vector of cell type annotations

- topic_range:

  Vector of topic numbers to test

- methods:

  Vector of evaluation methods

- n.iter:

  LDA iterations (default: 500)

- seed:

  Random seed

- save.models:

  Save intermediate models (default: FALSE)

- out.dir:

  Output directory for models

## Value

A list containing:

- metricsData frame of evaluation metrics

- best_modelBest LDA model

- best_n_topicsOptimal number of topics

- best_methodBest evaluation method

## Examples

``` r
if (FALSE) { # \dontrun{
results <- evaluate_topic_range(mat, celltypes, topic_range = 5:15)
} # }
```
