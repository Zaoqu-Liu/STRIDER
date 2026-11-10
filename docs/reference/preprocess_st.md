# Preprocess Spatial Data

Complete preprocessing pipeline for spatial transcriptomics data.

## Usage

``` r
preprocess_st(matrix, scale.factor = NULL, normalize = TRUE, scale = FALSE)
```

## Arguments

- matrix:

  Count matrix (genes x spots)

- scale.factor:

  Target sum for normalization (default: NULL)

- normalize:

  Whether to normalize (default: TRUE)

- scale:

  Whether to scale by SD (default: FALSE)

## Value

A list containing preprocessed data

## Examples

``` r
if (FALSE) { # \dontrun{
st_data <- preprocess_st(mat)
} # }
```
