# Visualize Spot-Cell Mapping

Create visualization of spot-cell mapping relationships.

## Usage

``` r
plot_mapping(mapping, celltypes, spots = NULL)
```

## Arguments

- mapping:

  Mapping result from map_spots_to_cells()

- celltypes:

  Named vector of cell type annotations

- spots:

  Specific spots to visualize (default: first 10)

## Value

A ggplot object
