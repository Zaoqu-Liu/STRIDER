# Plot Cell Type Composition Heatmap

Create heatmap of cell type composition across spots.

## Usage

``` r
plot_celltype_heatmap(
  deconv_mat,
  title = "Cell Type Composition",
  scale.by = "none"
)
```

## Arguments

- deconv_mat:

  Deconvolution matrix (spots x celltypes)

- title:

  Plot title

- scale.by:

  Scale values: "row", "column", or "none"

## Value

A ggplot object
