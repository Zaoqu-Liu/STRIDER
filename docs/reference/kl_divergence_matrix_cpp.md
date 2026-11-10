# Compute KL Divergence Matrix

Compute pairwise KL divergence between probability distributions (rows
of P and Q).

## Usage

``` r
kl_divergence_matrix_cpp(P, Q)
```

## Arguments

- P:

  First matrix (n x k), rows are distributions

- Q:

  Second matrix (m x k), rows are distributions

## Value

KL divergence matrix (n x m) where entry at row i, col j is KL(P_i \|\|
Q_j)

## Details

KL(p\|\|q) = sum_k p_k \* log(p_k / q_k)
