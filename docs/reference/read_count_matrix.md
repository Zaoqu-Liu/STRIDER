# Read Count Matrix from Text File

Read gene expression count matrix from tab-separated text file. First
column should contain gene names, first row should contain cell/spot
IDs.

## Usage

``` r
read_count_matrix(
  filename,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  sparse = TRUE
)
```

## Arguments

- filename:

  Path to the text file

- sep:

  Field separator (default: tab)

- header:

  Whether file has header (default: TRUE)

- row.names:

  Column containing row names (default: 1)

- sparse:

  Convert to sparse matrix (default: TRUE)

## Value

A list containing:

- matrixCount matrix (sparse or dense)

- featuresGene names

- barcodesCell/spot barcodes

## Examples

``` r
if (FALSE) { # \dontrun{
data <- read_count_matrix("gene_count.txt")
} # }
```
