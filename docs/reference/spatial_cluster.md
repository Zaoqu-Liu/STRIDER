# Spatial Clustering Based on Cell Type Composition

Perform clustering of spatial spots based on cell type composition with
optional neighborhood information.

## Usage

``` r
spatial_cluster(
  deconv,
  coords,
  n.clusters = 5,
  weight = 0.5,
  k = 6,
  method = "kmeans",
  seed = 42
)
```

## Arguments

- deconv:

  Deconvolution result (STRIDER object or matrix)

- coords:

  Spatial coordinates data frame

- n.clusters:

  Number of clusters (default: 5)

- weight:

  Weight for spot's own composition vs neighbors (default: 0.5)

  - 1.0: only own composition

  - 0.0: only neighbor composition

  - 0.5: equal weight (default)

- k:

  Number of neighbors for local composition (default: 6)

- method:

  Clustering method: "kmeans" or "hierarchical" (default: "kmeans")

- seed:

  Random seed

## Value

A data frame with spot coordinates and cluster assignments

## Details

The clustering combines each spot's own cell type composition with the
average composition of its spatial neighbors. This produces spatially
coherent clusters that reflect local tissue microenvironments.

## Examples

``` r
if (FALSE) { # \dontrun{
clusters <- spatial_cluster(strider_result, coords, n.clusters = 5)
} # }
```
