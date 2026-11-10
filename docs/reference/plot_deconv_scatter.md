# Scatter Plot for Deconvolution

Create scatter plot colored by dominant cell type.

## Usage

``` r
plot_deconv_scatter(
  deconv_mat,
  coords,
  pt.size = 2,
  colors = NULL,
  title = "Dominant Cell Type",
  coord.flip = FALSE
)
```

## Arguments

- deconv_mat:

  Deconvolution matrix (spots x celltypes)

- coords:

  Spatial coordinates

- pt.size:

  Point size

- colors:

  Color palette

- title:

  Plot title

- coord.flip:

  Flip coordinates

## Value

A ggplot object
