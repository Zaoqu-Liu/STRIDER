# Compute GW Gradient (Tensor Product)

Compute the gradient of GW loss w.r.t. transport plan. For squared loss:
L(C1, C2, T) = sum_ijkl (C1_ik - C2_jl)^2 \* T_ij \* T_kl Gradient:
dL/dT_ij = -4 \* sum_kl C1_ik \* T_kl \* C2_jl = -4 \* (C1 \* T \*
C2)\_ij

## Usage

``` r
gwgrad_cpp(C1, C2, T)
```

## Arguments

- C1:

  Source distance matrix (n x n)

- C2:

  Target distance matrix (m x m)

- T:

  Current transport plan (n x m)

## Value

Gradient matrix (n x m)
