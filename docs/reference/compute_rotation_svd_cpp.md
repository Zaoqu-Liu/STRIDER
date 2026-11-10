# Compute Optimal Rotation via SVD

Find the optimal rotation matrix R to align X to Y using Procrustes
analysis (SVD of cross-covariance).

## Usage

``` r
compute_rotation_svd_cpp(X, Y, weights)
```

## Arguments

- X:

  Source points (n x d)

- Y:

  Target points (n x d)

- weights:

  Point weights (n)

## Value

Rotation matrix R (d x d) such that X \* R aligns with Y

## Details

Minimizes weighted sum of squared alignment errors. Solution uses SVD: R
= U \* V' where H = Y' \* W \* X = U \* S \* V'
