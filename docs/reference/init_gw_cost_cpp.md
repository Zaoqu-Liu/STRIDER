# Initialize Constant Term for GW

Initialize the constant term in GW loss for squared loss. constC =
f1(C1) \* p \* 1_m^T + 1_n \* q^T \* f2(C2)^T where f1(a) = f2(a) = a^2
for squared loss

## Usage

``` r
init_gw_cost_cpp(C1, C2, p, q)
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

## Value

Constant cost matrix (n x m)
