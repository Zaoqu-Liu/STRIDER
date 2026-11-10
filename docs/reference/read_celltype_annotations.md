# Read Cell Type Annotations

Read cell type annotations from tab-separated text file. File should
have two columns: cell ID and cell type.

## Usage

``` r
read_celltype_annotations(filename, sep = "\t", header = FALSE)
```

## Arguments

- filename:

  Path to the annotation file

- sep:

  Field separator (default: tab)

- header:

  Whether file has header (default: FALSE)

## Value

A named vector mapping cell IDs to cell types

## Examples

``` r
if (FALSE) { # \dontrun{
celltypes <- read_celltype_annotations("celltype.txt")
} # }
```
