# Plot Integration Results

Visualize integrated samples.

## Usage

``` r
plot_integration(
  integration_result,
  color.by = "sample",
  pt.size = 2,
  deconv = NULL
)
```

## Arguments

- integration_result:

  Result from integrate_samples()

- color.by:

  Color spots by: "sample", "celltype", or NULL

- pt.size:

  Point size (default: 2)

- deconv:

  Deconvolution results for cell type coloring

## Value

A ggplot object
