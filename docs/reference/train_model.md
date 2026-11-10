# Train LDA Model (Convenience Wrapper)

Train an LDA model from single-cell data with automatic topic-celltype
association calculation.

## Usage

``` r
train_model(
  sc.matrix,
  celltypes,
  gene.use = NULL,
  n.topics = 10,
  n.iter = 500,
  seed = 42
)
```

## Arguments

- sc.matrix:

  Single-cell count matrix (genes x cells)

- celltypes:

  Named vector of cell type annotations

- gene.use:

  Genes to use (optional)

- n.topics:

  Number of topics

- n.iter:

  LDA iterations

- seed:

  Random seed

## Value

Trained LDA model object with topic-celltype matrices
