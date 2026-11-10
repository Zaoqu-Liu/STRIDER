# Find Marker Genes for Each Cell Type

Identify differentially expressed marker genes for each cell type using
Wilcoxon rank-sum test. This is crucial for training the LDA model.

## Usage

``` r
find_markers(
  matrix,
  celltypes,
  n.top = 200,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  only.pos = TRUE,
  test.use = "wilcox",
  adjust.method = "BH",
  ncores = 1
)
```

## Arguments

- matrix:

  Count matrix (genes x cells), should be normalized

- celltypes:

  Named vector of cell type annotations

- n.top:

  Number of top markers per cell type (default: 200)

- min.pct:

  Minimum fraction of cells expressing the gene (default: 0.1)

- logfc.threshold:

  Minimum log fold change threshold (default: 0.25)

- only.pos:

  Only return positive markers (default: TRUE)

- test.use:

  Statistical test to use: "wilcox" or "t" (default: "wilcox")

- adjust.method:

  P-value adjustment method (default: "BH")

- ncores:

  Number of cores for parallel computation (default: 1)

## Value

A list containing:

- markersData frame of all markers with statistics

- top_markersCharacter vector of top marker genes

## Examples

``` r
if (FALSE) { # \dontrun{
markers <- find_markers(norm_mat, celltypes, n.top = 200)
} # }
```
