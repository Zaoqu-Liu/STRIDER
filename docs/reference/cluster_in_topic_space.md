# Perform K-Means Clustering on Topic Space

Cluster cells in topic space and evaluate clustering quality.

## Usage

``` r
cluster_in_topic_space(
  topic_cell,
  true_celltypes,
  n_clusters = NULL,
  n_iter = 1000,
  seed = 42
)
```

## Arguments

- topic_cell:

  Topic-cell matrix (topics x cells)

- true_celltypes:

  Named vector of true cell type annotations

- n_clusters:

  Number of clusters (default: auto from celltypes)

- n_iter:

  Maximum k-means iterations (default: 1000)

- seed:

  Random seed

## Value

A list containing cluster assignments and NMI score
