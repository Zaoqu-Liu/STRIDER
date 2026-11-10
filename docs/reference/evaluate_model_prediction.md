# Evaluate Model Using Cell Type Prediction

Evaluate LDA model quality by predicting cell types and comparing with
known annotations.

## Usage

``` r
evaluate_model_prediction(
  topic_cell,
  topic_celltype,
  true_celltypes,
  method = "Raw",
  celltype_topic_bayes = NULL
)
```

## Arguments

- topic_cell:

  Topic-cell matrix (topics x cells)

- topic_celltype:

  Topic-celltype matrix (topics x celltypes)

- true_celltypes:

  Named vector of true cell type annotations

- method:

  Deconvolution method: "Raw", "Norm", "NormBySD", "Bayes", "BayesNorm"

- celltype_topic_bayes:

  Bayesian celltype-topic matrix (for Bayes methods)

## Value

A list containing:

- accuracyPrediction accuracy

- predictionsPredicted cell types

- confusionConfusion matrix

## Examples

``` r
if (FALSE) { # \dontrun{
eval <- evaluate_model_prediction(topic_cell, topic_ct, true_ct)
} # }
```
