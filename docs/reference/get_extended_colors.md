# Generate Color Palette

Generate a color palette with a specified number of colors.

## Usage

``` r
get_extended_colors(n, palette = "default")
```

## Arguments

- n:

  Number of colors needed

- palette:

  Base palette: "default" or "viridis"

## Value

A character vector of hex color codes

## Examples

``` r
get_extended_colors(10)
get_extended_colors(50, palette = "viridis")
```
