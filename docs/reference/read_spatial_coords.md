# Read Spatial Coordinates

Read spatial coordinates from tab-separated text file. First column
should contain spot IDs, followed by coordinates.

## Usage

``` r
read_spatial_coords(
  filename,
  sep = "\t",
  header = TRUE,
  spot.col = 1,
  coord.cols = c(2, 3),
  sample.col = NULL
)
```

## Arguments

- filename:

  Path to the text file

- sep:

  Field separator (default: tab)

- header:

  Whether file has header (default: TRUE)

- spot.col:

  Column containing spot IDs (default: 1)

- coord.cols:

  Columns containing coordinates (default: c(2, 3))

- sample.col:

  Column containing sample IDs (optional)

## Value

A data.frame with columns: spot, x, y (and optionally sample)

## Examples

``` r
if (FALSE) { # \dontrun{
coords <- read_spatial_coords("location.txt")
} # }
```
