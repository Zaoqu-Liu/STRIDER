# Filter Cells/Spots

Filter cells or spots based on minimum expression criteria.

## Usage

``` r
filter_cells(matrix, min.genes = 200, min.counts = 0, max.counts = Inf)
```

## Arguments

- matrix:

  Count matrix (genes x cells/spots)

- min.genes:

  Minimum number of genes detected (default: 200)

- min.counts:

  Minimum total counts (default: 0)

- max.counts:

  Maximum total counts (default: Inf)

## Value

Filtered matrix

## Examples

``` r
if (FALSE) { # \dontrun{
filtered <- filter_cells(mat, min.genes = 200)
} # }
```
