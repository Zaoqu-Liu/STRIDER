# Create Spot-Cell Neighborhood Graph

Create a bipartite graph connecting spots to their most similar cells.

## Usage

``` r
create_mapping_graph(mapping, threshold = NULL)
```

## Arguments

- mapping:

  Mapping result from map_spots_to_cells()

- threshold:

  Distance threshold for edges (optional)

## Value

An igraph object (if igraph is available) or edge list
