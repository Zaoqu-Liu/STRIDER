# Integrate Multiple Spatial Transcriptomics Samples

Integrate multiple spatial transcriptomics samples using Fused
Gromov-Wasserstein (FGW) optimal transport. This aligns samples based on
both topic distribution similarity and spatial structure.

## Usage

``` r
integrate_samples(
  deconv.list,
  topic.list,
  coords.list,
  sample.ids = NULL,
  alpha = 0.5,
  reg = 0.1,
  max.iter = 100,
  seed = 42
)
```

## Arguments

- deconv.list:

  List of deconvolution result matrices (spots x celltypes)

- topic.list:

  List of topic distribution matrices (topics x spots)

- coords.list:

  List of spatial coordinate data frames

- sample.ids:

  Sample identifiers (default: S1, S2, ...)

- alpha:

  Trade-off parameter between feature and spatial similarity (default:
  0.5)

  - 0: pure OT (only topic similarity)

  - 1: pure GW (only spatial structure)

- reg:

  Sinkhorn regularization parameter (default: 0.1)

- max.iter:

  Maximum iterations for FGW (default: 100)

- seed:

  Random seed

## Value

A list containing:

- aligned_coordsData frame of aligned coordinates for all samples

- transport_plansList of optimal transport plans between samples

- deconvCombined deconvolution results

- rotation_matricesList of rotation matrices

## Details

The integration is performed sequentially, aligning each sample to the
previous one. The FGW distance balances:

- Feature similarity: KL divergence between topic distributions

- Structural similarity: Gromov-Wasserstein distance on spatial
  structure

## References

Vayer, T., et al. (2020). Fused Gromov-Wasserstein distance for
structured objects. Algorithms.

## Examples

``` r
if (FALSE) { # \dontrun{
integrated <- integrate_samples(
  deconv.list = list(deconv1, deconv2, deconv3),
  topic.list = list(topic1, topic2, topic3),
  coords.list = list(coords1, coords2, coords3)
)
} # }
```
