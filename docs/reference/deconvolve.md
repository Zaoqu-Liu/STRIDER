# Cell Type Deconvolution Functions for STRIDER

Core functions for deconvolving cell type proportions in spatial
transcriptomics

Perform cell type deconvolution on spatial transcriptomics data using a
trained LDA model.

## Usage

``` r
deconvolve(
  st.matrix,
  model,
  topic_celltype,
  celltype_topic_bayes = NULL,
  method = "Bayes",
  n.iter = 100
)
```

## Arguments

- st.matrix:

  Spatial count matrix (genes x spots)

- model:

  Trained LDA model from train_lda()

- topic_celltype:

  Topic-celltype association matrix

- celltype_topic_bayes:

  Bayesian celltype-topic matrix (required for Bayes methods)

- method:

  Deconvolution method:

  - "Raw": Direct matrix multiplication

  - "Norm": Column-normalized weights

  - "NormBySD": SD-normalized weights

  - "Bayes": Bayesian posterior (default)

  - "BayesNorm": Normalized Bayesian posterior

- n.iter:

  Inference iterations (default: 100)

## Value

A list containing:

- celltype_fracCell type fractions (spots x celltypes)

- topic_spotTopic distributions (topics x spots)

## Details

The deconvolution computes cell type fractions for each spot by:

1.  Inferring topic distributions for spatial spots

2.  Applying celltype-topic weights to convert topics to cell types

3.  Normalizing to obtain fractions that sum to 1
