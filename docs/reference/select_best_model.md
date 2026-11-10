# Select Best Model

Automatically select the best model from evaluation results.

## Usage

``` r
select_best_model(eval_results, criterion = "accuracy")
```

## Arguments

- eval_results:

  Results from evaluate_topic_range()

- criterion:

  Selection criterion: "accuracy", "coherence", or "balanced"

## Value

Selected model components
