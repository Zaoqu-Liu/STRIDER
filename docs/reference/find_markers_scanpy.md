# Find Marker Genes (scanpy-style workflow)

Identify differentially expressed marker genes for each cell type using
the exact scanpy/STRIDE workflow including preprocessing steps.

## Usage

``` r
find_markers_scanpy(
  matrix,
  celltypes,
  n.top = 200,
  min.genes = 200,
  min.cells = 10,
  target.sum = 10000,
  n.hvg = NULL,
  ncores = 1
)
```

## Arguments

- matrix:

  Raw count matrix (genes x cells)

- celltypes:

  Named vector of cell type annotations

- n.top:

  Number of top markers per cell type (default: 200)

- min.genes:

  Minimum genes per cell for filtering (default: 200)

- min.cells:

  Minimum cells per gene for filtering (default: 10)

- target.sum:

  Target sum for normalization (default: 1e4)

- n.hvg:

  Number of highly variable genes to use (default: NULL, use all)

- ncores:

  Number of cores for parallel computation (default: 1)

## Value

A list containing:

- markersData frame of all markers with statistics

- top_markersCharacter vector of top marker genes

## Details

This function replicates the scanpy workflow used in Python STRIDE:

1.  Filter cells with \< min.genes expressed genes

2.  Filter genes expressed in \< min.cells cells

3.  Normalize counts to target.sum per cell

4.  Log-transform: log1p(x)

5.  Optionally select highly variable genes

6.  Run rank_genes_groups with Wilcoxon test and tie correction

## Examples

``` r
if (FALSE) { # \dontrun{
markers <- find_markers_scanpy(raw_count_mat, celltypes, n.top = 200)
} # }
```
