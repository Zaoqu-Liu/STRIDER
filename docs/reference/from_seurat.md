# Create STRIDER Object from Seurat

Extract data from Seurat object for STRIDER analysis.

## Usage

``` r
from_seurat(seurat, celltype.col, assay = "RNA")
```

## Arguments

- seurat:

  Seurat object

- celltype.col:

  Column name in metadata containing cell types

- assay:

  Assay to use (default: "RNA")

## Value

List with matrix and celltypes
