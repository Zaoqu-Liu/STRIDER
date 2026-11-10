# Get Shared Genes Between Datasets

Find genes shared between single-cell and spatial datasets, optionally
filtered by a provided gene list.

## Usage

``` r
get_shared_genes(sc_genes, st_genes, gene_use = NULL)
```

## Arguments

- sc_genes:

  Character vector of single-cell genes

- st_genes:

  Character vector of spatial genes

- gene_use:

  Optional character vector of genes to use, or "All"

## Value

Character vector of shared genes

## Examples

``` r
if (FALSE) { # \dontrun{
shared <- get_shared_genes(sc_genes, st_genes, marker_genes)
} # }
```
