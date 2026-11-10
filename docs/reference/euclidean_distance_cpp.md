# Compute Euclidean Distance Matrix

Compute pairwise Euclidean distances between points.

## Usage

``` r
euclidean_distance_cpp(X, Y_nullable = NULL)
```

## Arguments

- X:

  First point set (n x d)

- Y_nullable:

  Second point set (m x d), optional. If NULL, computes pairwise
  distances within X.

## Value

Distance matrix (n x m or n x n)
