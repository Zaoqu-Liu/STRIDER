# Compute Rotation Matrix

Compute optimal rotation matrix to align two point sets using the
optimal transport plan.

## Usage

``` r
compute_rotation(coords1, coords2, transport, p, q)
```

## Arguments

- coords1:

  Coordinates of sample 1

- coords2:

  Coordinates of sample 2

- transport:

  Optimal transport plan

- p:

  Source weights

- q:

  Target weights

## Value

Rotation matrix (2 x 2)
