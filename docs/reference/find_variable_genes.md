# Find Highly Variable Genes

Identify highly variable genes using the variance stabilization approach
(similar to Seurat's VST method).

## Usage

``` r
find_variable_genes(
  matrix,
  n.top = 2000,
  min.mean = 0.0125,
  max.mean = 3,
  min.var = 0.5,
  span = 0.3
)
```

## Arguments

- matrix:

  Count matrix (genes x cells)

- n.top:

  Number of top variable genes (default: 2000)

- min.mean:

  Minimum mean expression (default: 0.0125)

- max.mean:

  Maximum mean expression (default: 3)

- min.var:

  Minimum variance (default: 0.5)

- span:

  Loess span parameter (default: 0.3)

## Value

A list containing:

- hvgCharacter vector of highly variable genes

- statsData frame with gene statistics

## Examples

``` r
if (FALSE) { # \dontrun{
hvg <- find_variable_genes(mat, n.top = 2000)
} # }
```
