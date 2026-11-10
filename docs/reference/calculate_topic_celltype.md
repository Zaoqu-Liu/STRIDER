# Calculate Topic-Celltype Association Matrix

Calculate the association between topics and cell types by aggregating
topic distributions across cells of each type.

## Usage

``` r
calculate_topic_celltype(topic_cell, celltypes, method = "mean")
```

## Arguments

- topic_cell:

  Topic-cell matrix (topics x cells)

- celltypes:

  Named vector of cell type annotations

- method:

  Aggregation method: "mean" (default), "median", or "sum"

## Value

Topic-celltype matrix (topics x celltypes)

## Details

This matrix represents P(topic\|celltype), the probability of each topic
given a cell type. It is computed by averaging topic distributions
across all cells of each type.

## Examples

``` r
if (FALSE) { # \dontrun{
topic_ct <- calculate_topic_celltype(topic_cell, celltypes)
} # }
```
