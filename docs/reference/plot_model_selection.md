# Plot Model Selection Results

Visualize model evaluation metrics across topic numbers.

## Usage

``` r
plot_model_selection(metrics, metric = "accuracy")
```

## Arguments

- metrics:

  Metrics data frame from evaluate_topic_range()

- metric:

  Which metric to plot: "accuracy", "coherence", or "nmi"

## Value

A ggplot object
