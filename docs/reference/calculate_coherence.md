# Calculate LDA Coherence Score

Calculate topic coherence scores to evaluate model quality. Implements
UMass and optionally C_v coherence metrics.

## Usage

``` r
calculate_coherence(
  model,
  matrix,
  n.words = 10,
  type = "both",
  cell_gene_list = NULL
)
```

## Arguments

- model:

  LDA model from train_lda()

- matrix:

  Original count matrix used for training

- n.words:

  Number of top words per topic (default: 10)

- type:

  Coherence type: "umass", "cv", or "both" (default: "both")

- cell_gene_list:

  List of genes expressed per cell (for c_v, auto-generated if NULL)

## Value

A list containing:

- umassUMass coherence (mean across topics)

- cvC_v coherence (mean across topics), NA if type="umass"

- umass_per_topicUMass score for each topic

- cv_per_topicC_v score for each topic

## Details

UMass coherence measures how often top words co-occur in the corpus:
C_UMass = sum_j=2^M sum_i=1^j-1 log((D(w_i, w_j) + 1) / D(w_i))

C_v coherence uses normalized pointwise mutual information (NPMI) with a
sliding window and cosine similarity, providing better correlation with
human judgments.

## References

Mimno, D., et al. (2011). Optimizing semantic coherence in topic models.
EMNLP. Röder, M., et al. (2015). Exploring the Space of Topic Coherence
Measures. WSDM.
