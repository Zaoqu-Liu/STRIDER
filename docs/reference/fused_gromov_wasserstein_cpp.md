# Fused Gromov-Wasserstein Distance

Compute Fused Gromov-Wasserstein (FGW) distance, which combines
feature-based Wasserstein distance and structure-based
Gromov-Wasserstein distance.

## Usage

``` r
fused_gromov_wasserstein_cpp(
  M,
  C1,
  C2,
  p,
  q,
  alpha = 0.5,
  reg = 0.1,
  max_iter = 100L,
  tol = 1e-07
)
```

## Arguments

- M:

  Feature distance/cost matrix (n x m)

- C1:

  Source structural distance matrix (n x n)

- C2:

  Target structural distance matrix (m x m)

- p:

  Source distribution (n), must sum to 1

- q:

  Target distribution (m), must sum to 1

- alpha:

  Trade-off parameter between 0 and 1 (default: 0.5)

  - alpha = 0: pure Wasserstein (only feature cost)

  - alpha = 1: pure GW (only structural cost)

- reg:

  Entropic regularization (default: 0.1)

- max_iter:

  Maximum iterations (default: 100)

- tol:

  Convergence tolerance (default: 1e-7)

## Value

List with:

- T: Optimal transport plan (n x m)

- fgw_dist: FGW distance value

## Details

The FGW distance is defined as: FGW(p,q) = min_T (1-alpha) \* \<M, T\> +
alpha \* GW(C1, C2, T) subject to T \* 1 = p, T^T \* 1 = q

## References

Vayer, T., et al. (2020). Fused Gromov-Wasserstein distance for
structured objects. Algorithms.
