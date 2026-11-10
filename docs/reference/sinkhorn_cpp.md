# Sinkhorn-Knopp Algorithm

Compute optimal transport plan using the Sinkhorn-Knopp algorithm with
entropic regularization.

## Usage

``` r
sinkhorn_cpp(M, p, q, reg = 0.1, max_iter = 1000L, tol = 1e-09)
```

## Arguments

- M:

  Cost matrix (n x m)

- p:

  Source distribution (n), must sum to 1

- q:

  Target distribution (m), must sum to 1

- reg:

  Entropic regularization parameter (default: 0.1)

- max_iter:

  Maximum iterations (default: 1000)

- tol:

  Convergence tolerance (default: 1e-9)

## Value

Transport plan matrix (n x m)
