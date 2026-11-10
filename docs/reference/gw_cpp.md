# Gromov-Wasserstein Distance

Compute Gromov-Wasserstein distance between two metric spaces using
conditional gradient descent with Sinkhorn projections.

## Usage

``` r
gw_cpp(C1, C2, p, q, reg = 0.1, max_iter = 100L, tol = 1e-07)
```

## Arguments

- C1:

  Source distance matrix (n x n)

- C2:

  Target distance matrix (m x m)

- p:

  Source distribution (n)

- q:

  Target distribution (m)

- reg:

  Entropic regularization (default: 0.1)

- max_iter:

  Maximum iterations (default: 100)

- tol:

  Convergence tolerance (default: 1e-7)

## Value

List with transport plan T and GW distance
