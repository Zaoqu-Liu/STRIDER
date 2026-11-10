# Validate STRIDER Results Against Reference

Compare STRIDER deconvolution results with reference results (e.g., from
Python STRIDE) to assess scientific equivalence.

## Usage

``` r
validate_strider(
  strider_result,
  reference_result,
  metrics = c("correlation", "dominant", "error", "rank"),
  tolerance = 0.01
)
```

## Arguments

- strider_result:

  STRIDER result object or deconvolution matrix (spots x celltypes)

- reference_result:

  Reference result (matrix, data frame, or file path)

- metrics:

  Which metrics to compute (default: all)

- tolerance:

  Numerical tolerance for "exact" comparisons (default: 0.01)

## Value

A validation report object with:

- pearson_rPearson correlation of flattened matrices

- spearman_rSpearman correlation of flattened matrices

- rmseRoot mean squared error

- maeMean absolute error

- dominant_matchFraction of spots with same dominant cell type

- rank_correlationMean rank correlation per spot

- verdict"EQUIVALENT" if metrics pass thresholds

## Details

The validation assesses scientific equivalence rather than numerical
identity. Two results are considered equivalent if:

- Pearson correlation \> 0.95

- Dominant cell type match \> 0.90

- RMSE \< 0.05

## Examples

``` r
if (FALSE) { # \dontrun{
validation <- validate_strider(strider_result, "stride_result.txt")
print(validation)
} # }
```
