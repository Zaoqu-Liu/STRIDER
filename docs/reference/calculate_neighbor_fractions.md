# Calculate Neighbor Average Fractions

Calculate the average cell type fractions from neighboring spots.

## Usage

``` r
calculate_neighbor_fractions(deconv, neighbors)
```

## Arguments

- deconv:

  Cell type fraction matrix (spots x celltypes)

- neighbors:

  List of neighbor indices from find_neighbors()

## Value

Matrix of neighbor-averaged fractions (spots x celltypes)
