# Scatter Pie Plot for Deconvolution

Create scatter pie plot showing cell type composition at each spot.

## Usage

``` r
plot_deconv_pie(
  deconv_mat,
  coords,
  pt.size = 8,
  colors = NULL,
  title = "Cell Type Distribution",
  coord.flip = FALSE
)
```

## Arguments

- deconv_mat:

  Deconvolution matrix (spots x celltypes)

- coords:

  Spatial coordinates

- pt.size:

  Point/pie size

- colors:

  Color palette

- title:

  Plot title

- coord.flip:

  Flip coordinates

## Value

A ggplot object
