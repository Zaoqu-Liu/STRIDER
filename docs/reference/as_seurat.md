# Convert STRIDER Object to Seurat

Add STRIDER deconvolution results to a Seurat object.

## Usage

``` r
as_seurat(strider, seurat, assay.name = "STRIDER")
```

## Arguments

- strider:

  STRIDER result object

- seurat:

  Seurat object

- assay.name:

  Name for new assay (default: "STRIDER")

## Value

Seurat object with deconvolution results
