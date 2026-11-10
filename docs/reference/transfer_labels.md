# Transfer Labels from Single Cells to Spots

Transfer cell type labels from single-cell reference to spatial spots
using mapping results.

## Usage

``` r
transfer_labels(mapping, celltypes, method = "knn", k = NULL)
```

## Arguments

- mapping:

  Mapping result from map_spots_to_cells()

- celltypes:

  Named vector of cell type annotations

- method:

  Transfer method: "knn" (k-nearest neighbor voting) or "weighted"

- k:

  Number of neighbors for kNN (default: all mapped cells)

## Value

Named vector of transferred labels
