# Scale Data by Standard Deviation

Scale data by dividing by standard deviation (without centering). Uses
population SD (N) instead of sample SD (N-1) to match sklearn's
StandardScaler(with_mean=False) used in Python STRIDE.

## Usage

``` r
scale_data(matrix, center = FALSE, scale = TRUE)
```

## Arguments

- matrix:

  Normalized matrix (genes x cells/spots)

- center:

  Center data (subtract mean) (default: FALSE)

- scale:

  Scale data (divide by SD) (default: TRUE)

## Value

Scaled matrix

## Examples

``` r
if (FALSE) { # \dontrun{
scaled <- scale_data(mat)
} # }
```
