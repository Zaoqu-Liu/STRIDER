# Preprocess Single-Cell Data

Complete preprocessing pipeline for single-cell RNA-seq data.

## Usage

``` r
preprocess_sc(
  matrix,
  celltypes,
  scale.factor = NULL,
  min.cells = 10,
  min.genes = 200,
  normalize = TRUE,
  scale = FALSE
)
```

## Arguments

- matrix:

  Count matrix (genes x cells)

- celltypes:

  Named vector of cell type annotations

- scale.factor:

  Target sum for normalization (default: NULL)

- min.cells:

  Minimum cells per gene (default: 10)

- min.genes:

  Minimum genes per cell (default: 200)

- normalize:

  Whether to normalize (default: TRUE)

- scale:

  Whether to scale by SD (default: FALSE)

## Value

A list containing preprocessed data

## Examples

``` r
if (FALSE) { # \dontrun{
sc_data <- preprocess_sc(mat, celltypes)
} # }
```
