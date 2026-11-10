# Apply LDA Model to New Data

Apply a trained LDA model to transform new data (e.g., spatial spots)
into topic space using variational inference.

## Usage

``` r
transform_lda(model, newdata, n.iter = 100)
```

## Arguments

- model:

  Trained LDA model from train_lda()

- newdata:

  Matrix to transform (genes x cells/spots)

- n.iter:

  Number of inference iterations (default: 100)

## Value

Topic distribution matrix (topics x cells/spots)

## Details

This function performs inference on new documents using the trained
topic-word distributions. Genes not present in the training vocabulary
are ignored.

## Examples

``` r
if (FALSE) { # \dontrun{
topic_dist <- transform_lda(model, st_matrix)
} # }
```
