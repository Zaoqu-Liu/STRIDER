# Run STRIDER Deconvolution Pipeline

Main entry point for STRIDER cell type deconvolution. This function
performs the complete pipeline including marker identification, LDA
training, model selection, and spatial deconvolution.

## Usage

``` r
run_strider(
  sc.data,
  st.data,
  sc.celltype,
  gene.use = NULL,
  n.topics = NULL,
  method = "auto",
  sc.scale.factor = NULL,
  st.scale.factor = NULL,
  normalize = FALSE,
  n.markers = 200,
  lda.iter = 500,
  seed = 42,
  out.dir = ".",
  out.prefix = "STRIDER",
  save.model = TRUE
)
```

## Arguments

- sc.data:

  Single-cell data. Can be:

  - A matrix (genes x cells)

  - A Seurat object (V4 or V5)

  - A file path to count matrix (.h5 or .txt)

- st.data:

  Spatial transcriptomics data. Same formats as sc.data.

- sc.celltype:

  Cell type annotations. Can be:

  - A named vector (cell -\> celltype)

  - Column name in Seurat metadata

  - A file path to annotation file

- gene.use:

  Genes to use for model training:

  - NULL: auto-detect markers (default)

  - "All": use all shared genes

  - Character vector: specific genes

  - File path: gene list file

- n.topics:

  Number of topics, or vector of numbers to test. Default:
  seq(n_celltypes, 3\*n_celltypes)

- method:

  Deconvolution method: "auto", "Raw", "Norm", "NormBySD", "Bayes",
  "BayesNorm"

- sc.scale.factor:

  Scale factor for single-cell normalization (default: auto)

- st.scale.factor:

  Scale factor for spatial normalization (default: auto)

- normalize:

  Whether to scale by SD (default: FALSE)

- n.markers:

  Number of marker genes per cell type (default: 200)

- lda.iter:

  LDA iterations (default: 500)

- seed:

  Random seed (default: 42)

- out.dir:

  Output directory (default: ".")

- out.prefix:

  Output file prefix (default: "STRIDER")

- save.model:

  Save trained model (default: TRUE)

## Value

A STRIDER object containing:

- deconvCell type fractions matrix (spots x celltypes)

- topic_spotTopic distributions for spots (topics x spots)

- modelTrained LDA model

- methodSelected deconvolution method

- n_topicsSelected number of topics

- metricsModel evaluation metrics

## Details

The STRIDER pipeline consists of the following steps:

1.  Load and preprocess single-cell and spatial data

2.  Identify marker genes for each cell type

3.  Train LDA topic model on single-cell data

4.  Evaluate models with different topic numbers

5.  Select optimal model and deconvolution method

6.  Apply model to spatial data for deconvolution

## References

Sun, D., et al. (2022). STRIDE: accurately decomposing and integrating
spatial transcriptomics using single-cell RNA sequencing. Nucleic Acids
Research.

## Examples

``` r
if (FALSE) { # \dontrun{
# From matrices
result <- run_strider(
  sc.data = sc_matrix,
  st.data = st_matrix,
  sc.celltype = celltypes
)

# From Seurat objects
result <- run_strider(
  sc.data = seurat_sc,
  st.data = seurat_st,
  sc.celltype = "celltype"
)
} # }
```
