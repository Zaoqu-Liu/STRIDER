# Read 10X Genomics HDF5 File

Read gene expression data from 10X Genomics HDF5 format. Supports both
Cell Ranger v2 and v3 formats.

## Usage

``` r
read_10x_h5(filename, use.names = TRUE, unique.features = TRUE)
```

## Arguments

- filename:

  Path to the HDF5 file

- use.names:

  Use gene names instead of IDs (default: TRUE)

- unique.features:

  Make feature names unique (default: TRUE)

## Value

A list containing:

- matrixSparse dgCMatrix of counts

- featuresGene names/IDs

- barcodesCell/spot barcodes

## Examples

``` r
if (FALSE) { # \dontrun{
data <- read_10x_h5("filtered_feature_bc_matrix.h5")
} # }
```
