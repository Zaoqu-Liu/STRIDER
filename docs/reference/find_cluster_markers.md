# Find Enriched Cell Types per Cluster

Identify cell types that are enriched in each cluster compared to the
overall distribution using statistical testing.

## Usage

``` r
find_cluster_markers(
  deconv,
  clusters,
  method = "wilcox",
  p.adjust.method = "BH"
)
```

## Arguments

- deconv:

  Deconvolution matrix (spots x celltypes)

- clusters:

  Cluster assignments

- method:

  Statistical test: "wilcox" or "t" (default: "wilcox")

- p.adjust.method:

  P-value adjustment method (default: "BH")

## Value

Data frame with enrichment statistics per cluster-celltype pair
