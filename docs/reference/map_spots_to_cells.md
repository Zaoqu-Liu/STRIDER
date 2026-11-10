# Map Spatial Spots to Similar Single Cells

Identify the most similar single cells for each spatial spot based on
topic distribution similarity.

## Usage

``` r
map_spots_to_cells(topic.spot, topic.cell, n.top = 10, method = "euclidean")
```

## Arguments

- topic.spot:

  Topic distribution for spots (topics x spots)

- topic.cell:

  Topic distribution for cells (topics x cells)

- n.top:

  Number of top similar cells per spot (default: 10)

- method:

  Distance metric: "euclidean", "cosine", or "correlation" (default:
  "euclidean")

## Value

A data frame with spot-cell mappings and similarity scores

## Examples

``` r
if (FALSE) { # \dontrun{
mapping <- map_spots_to_cells(topic_spot, topic_cell, n.top = 10)
} # }
```
