# Calculate Cluster Statistics

Calculate cell type composition statistics for each cluster.

## Usage

``` r
cluster_composition(deconv, clusters)
```

## Arguments

- deconv:

  Deconvolution matrix (spots x celltypes)

- clusters:

  Cluster assignments (vector or data frame with 'cluster' column)

## Value

Data frame with mean composition per cluster
