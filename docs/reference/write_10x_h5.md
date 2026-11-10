# Write 10X Genomics HDF5 File

Write gene expression data to 10X Genomics HDF5 format.

## Usage

``` r
write_10x_h5(matrix, features, barcodes, filename, overwrite = FALSE)
```

## Arguments

- matrix:

  Sparse matrix (genes x cells/spots)

- features:

  Gene names

- barcodes:

  Cell/spot barcodes

- filename:

  Output file path

- overwrite:

  Overwrite existing file (default: FALSE)

## Value

Invisible NULL

## Examples

``` r
if (FALSE) { # \dontrun{
write_10x_h5(mat, genes, cells, "output.h5")
} # }
```
