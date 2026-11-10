# Calculate Bayesian Celltype-Topic Matrix

Calculate celltype-topic association using Bayes' theorem, incorporating
prior probabilities of cell types.

## Usage

``` r
calculate_celltype_topic_bayes(topic_celltype, celltype_counts, total_cells)
```

## Arguments

- topic_celltype:

  Topic-celltype matrix (topics x celltypes)

- celltype_counts:

  Named vector of cell counts per cell type

- total_cells:

  Total number of cells

## Value

Celltype-topic Bayes matrix (celltypes x topics)

## Details

Applies Bayes' theorem: P(celltype\|topic) = P(topic\|celltype) \*
P(celltype) / P(topic)

This accounts for differences in cell type abundance in the reference
data, preventing overrepresentation of abundant cell types in
deconvolution.
