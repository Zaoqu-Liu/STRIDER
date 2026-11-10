# Get Mapped Cell Information

Extract information about cells mapped to each spot.

## Usage

``` r
get_mapped_celltypes(mapping, celltypes, aggregate = "majority")
```

## Arguments

- mapping:

  Mapping result from map_spots_to_cells()

- celltypes:

  Named vector of cell type annotations

- aggregate:

  How to aggregate cell types: "majority", "all", or "weighted"

## Value

Data frame with spot-level cell type information
