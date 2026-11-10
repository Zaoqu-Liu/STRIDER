# Calculate Quality Control Metrics

Calculate quality control metrics for cells/spots.

## Usage

``` r
calculate_qc_metrics(matrix, mito.pattern = "^MT-|^mt-")
```

## Arguments

- matrix:

  Count matrix (genes x cells/spots)

- mito.pattern:

  Pattern to match mitochondrial genes (default: "^MT-\|^mt-")

## Value

Data frame with QC metrics

## Examples

``` r
if (FALSE) { # \dontrun{
qc <- calculate_qc_metrics(mat)
} # }
```
