# Plot Topic Distribution Heatmap

Create heatmap of topic distributions across samples.

## Usage

``` r
plot_topic_heatmap(
  topic_mat,
  celltypes = NULL,
  title = "Topic Distribution",
  cluster.rows = TRUE,
  cluster.cols = FALSE
)
```

## Arguments

- topic_mat:

  Topic distribution matrix (topics x cells/spots)

- celltypes:

  Named vector of cell type annotations (optional)

- title:

  Plot title

- cluster.rows:

  Cluster rows (default: TRUE)

- cluster.cols:

  Cluster columns (default: FALSE)

## Value

A ggplot object
