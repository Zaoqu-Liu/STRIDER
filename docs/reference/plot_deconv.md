# Plot Deconvolution Results

Main plotting function for STRIDER deconvolution results. Supports
scatter pie plots, scatter plots, and heatmaps.

## Usage

``` r
plot_deconv(
  deconv,
  coords,
  plot.type = "scatterpie",
  sample.id = NULL,
  pt.size = NULL,
  colors = NULL,
  title = NULL,
  coord.flip = FALSE
)
```

## Arguments

- deconv:

  Deconvolution result (STRIDER object or matrix)

- coords:

  Spatial coordinates data frame

- plot.type:

  Type of plot: "scatterpie", "scatter", or "heatmap"

- sample.id:

  Sample ID to plot (if multiple samples)

- pt.size:

  Point size (default: auto)

- colors:

  Color palette (default: STRIDER palette)

- title:

  Plot title (default: auto)

- coord.flip:

  Flip x and y coordinates (default: FALSE)

## Value

A ggplot object

## Examples

``` r
if (FALSE) { # \dontrun{
plot_deconv(strider_result, coords, plot.type = "scatterpie")
} # }
```
