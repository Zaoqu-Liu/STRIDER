# Align Two Samples Using FGW

Align two spatial transcriptomics samples using Fused Gromov-Wasserstein
optimal transport.

## Usage

``` r
align_samples(
  topic1,
  topic2,
  coords1,
  coords2,
  alpha = 0.5,
  reg = 0.1,
  max.iter = 100
)
```

## Arguments

- topic1:

  Topic distribution for sample 1 (topics x spots)

- topic2:

  Topic distribution for sample 2 (topics x spots)

- coords1:

  Spatial coordinates for sample 1

- coords2:

  Spatial coordinates for sample 2

- alpha:

  Trade-off parameter (default: 0.5)

- reg:

  Regularization parameter (default: 0.1)

- max.iter:

  Maximum iterations (default: 100)

## Value

A list containing:

- transportOptimal transport plan matrix

- rotationRotation matrix for alignment

- fgw_distFGW distance value

## Examples

``` r
if (FALSE) { # \dontrun{
alignment <- align_samples(topic1, topic2, coords1, coords2)
} # }
```
