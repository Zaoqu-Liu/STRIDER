# Find Highly Variable Genes (scanpy-style)

Identify highly variable genes using the scanpy algorithm. This matches
sc.pp.highly_variable_genes() with flavor='seurat'.

## Usage

``` r
find_variable_genes_scanpy(
  matrix,
  n.top = 2000,
  min.mean = 0.0125,
  max.mean = 3,
  min.disp = 0.5
)
```

## Arguments

- matrix:

  Log-normalized count matrix (genes x cells)

- n.top:

  Number of top variable genes (default: 2000)

- min.mean:

  Minimum mean expression (default: 0.0125)

- max.mean:

  Maximum mean expression (default: 3)

- min.disp:

  Minimum dispersion (default: 0.5)

## Value

A list with hvg vector and stats data frame
