# Find Spatial Neighbors

Find neighbors for each spot based on spatial coordinates using
k-nearest neighbors or distance threshold.

## Usage

``` r
find_neighbors(coords, k = 6, radius = NULL, include.self = FALSE)
```

## Arguments

- coords:

  Data frame with spatial coordinates (columns: x, y)

- k:

  Number of nearest neighbors (default: 6)

- radius:

  Distance radius for neighbors (alternative to k)

- include.self:

  Include the spot itself as a neighbor (default: FALSE)

## Value

A list where each element contains indices of neighbors for each spot

## Details

For Visium data, k=6 corresponds to the hexagonal grid structure. For
other technologies, adjust k based on expected spatial relationships.

## Examples

``` r
if (FALSE) { # \dontrun{
neighbors <- find_neighbors(spatial_coords, k = 6)
} # }
```
