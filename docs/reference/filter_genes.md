# Filter Genes

Filter genes based on minimum expression criteria.

## Usage

``` r
filter_genes(matrix, min.cells = 10, min.counts = 0)
```

## Arguments

- matrix:

  Count matrix (genes x cells/spots)

- min.cells:

  Minimum number of cells expressing the gene (default: 10)

- min.counts:

  Minimum total counts for a gene (default: 0)

## Value

Filtered matrix

## Examples

``` r
if (FALSE) { # \dontrun{
filtered <- filter_genes(mat, min.cells = 10)
} # }
```
