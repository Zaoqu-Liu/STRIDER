# Calculate C_v Coherence for One Topic

Internal function to calculate C_v coherence. Implements NPMI-based
coherence with sliding window.

## Usage

``` r
.calculate_cv_coherence(top_words, texts, window_size = 10)
```

## Arguments

- top_words:

  Top words for the topic

- texts:

  List of word lists (documents)

- window_size:

  Sliding window size

## Value

C_v coherence score
