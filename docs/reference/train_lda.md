# Train LDA Topic Model

Train a Latent Dirichlet Allocation (LDA) topic model on single-cell
gene expression data. The model treats cells as documents and genes as
words, following the methodology described in STRIDE.

## Usage

``` r
train_lda(
  matrix,
  n.topics = 10,
  n.iter = 500,
  n.passes = 1,
  convergence.tol = 1e-04,
  alpha = NULL,
  beta = 0.01,
  gamma.threshold = 0.001,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- matrix:

  Normalized count matrix (genes x cells)

- n.topics:

  Number of topics (default: 10)

- n.iter:

  Number of iterations per pass (default: 500)

- n.passes:

  Number of passes over the corpus (default: 1)

- convergence.tol:

  Convergence tolerance (default: 1e-4)

- alpha:

  Dirichlet prior for document-topic distribution. Default is
  50/n.topics as in original STRIDE (matching gensim).

- beta:

  Dirichlet prior for topic-word distribution (default: 0.01)

- gamma.threshold:

  Minimum change in gamma for document convergence (default: 0.001)

- seed:

  Random seed for reproducibility

- verbose:

  Print progress messages (default: TRUE)

## Value

An LDA model object with components:

- topic_geneTopic-gene distribution matrix (phi)

- topic_cellDocument-topic distribution matrix (theta)

- vocabGene vocabulary

- n_topicsNumber of topics

## Details

The LDA model is trained using Variational Bayes inference via the
text2vec package, matching gensim's LdaModel behavior as closely as
possible. The default hyperparameters (alpha = 50/K, beta = 0.01) follow
the original STRIDE implementation for optimal performance on
single-cell data.

Note: gensim uses Online Variational Bayes by default, while text2vec
uses standard Variational Bayes. For single-cell data with moderate
corpus size, results are highly correlated.

## References

Sun, D., et al. (2022). STRIDE: accurately decomposing and integrating
spatial transcriptomics using single-cell RNA sequencing. Nucleic Acids
Research.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- train_lda(norm_mat, n.topics = 10)
} # }
```
