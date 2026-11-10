# Compute KL Divergence Matrix

Compute pairwise KL divergence between probability distributions.
Matches the exact implementation in Python STRIDE's KLDivergency
function.

## Usage

``` r
compute_kl_divergence(P, Q, eps = 1e-10)
```

## Arguments

- P:

  First matrix (n x k), rows are distributions

- Q:

  Second matrix (m x k), rows are distributions

- eps:

  Small constant to avoid log(0) (default: 1e-10)

## Value

KL divergence matrix (n x m)
