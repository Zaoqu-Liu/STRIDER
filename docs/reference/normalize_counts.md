# Normalize Counts

Normalize count matrix by total counts per cell/spot.

## Usage

``` r
normalize_counts(matrix, scale.factor = NULL, log.transform = FALSE)
```

## Arguments

- matrix:

  Count matrix (genes x cells/spots)

- scale.factor:

  Target sum for normalization (default: NULL, auto-detect)

- log.transform:

  Apply log1p transformation (default: FALSE)

## Value

Normalized matrix

## Examples

``` r
if (FALSE) { # \dontrun{
normalized <- normalize_counts(mat, scale.factor = 10000)
} # }
```
