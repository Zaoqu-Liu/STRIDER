# Create STRIDER Object

Create a STRIDER result object.

## Usage

``` r
create_strider_object(
  deconv,
  topic_spot = NULL,
  model = NULL,
  method = NULL,
  n_topics = NULL,
  metrics = NULL,
  genes_used = NULL,
  celltypes = NULL
)
```

## Arguments

- deconv:

  Deconvolution results matrix (spots x celltypes)

- topic_spot:

  Topic-spot distribution matrix

- model:

  Trained LDA model

- method:

  Selected deconvolution method

- n_topics:

  Number of topics

- metrics:

  Model evaluation metrics

- genes_used:

  Genes used in model training

- celltypes:

  Cell types in the model

## Value

A STRIDER object
