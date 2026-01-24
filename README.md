# STRIDER <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/Zaoqu-Liu/STRIDER/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Zaoqu-Liu/STRIDER/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**STRIDER** (**S**patial **T**ranscriptomics deconvolu**T**ion and integ**R**ation **I**n **D**eveloped **E**nvironment for **R**) is an R package for comprehensive analysis of spatial transcriptomics data, including cell type deconvolution, multi-sample integration, and spatial clustering.

## Overview

STRIDER implements:

- **Cell Type Deconvolution**: Topic modeling (LDA) based decomposition of spatial spots into cell type proportions
- **Multi-Sample Integration**: Fused Gromov-Wasserstein (FGW) optimal transport for aligning multiple spatial samples
- **Spatial Clustering**: Neighborhood-aware clustering of spatial domains
- **Visualization**: Publication-ready scatter pie plots, heatmaps, and integration visualizations

## Installation

### From GitHub

```r
# Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

# Install STRIDER
devtools::install_github("Zaoqu-Liu/STRIDER")
```

### Dependencies

STRIDER requires several packages. Install them with:

```r
# Required packages
install.packages(c("Matrix", "data.table", "ggplot2", "viridis", "Rcpp"))

# For LDA topic modeling
install.packages("text2vec")

# For spatial neighbor analysis
install.packages("FNN")

# For scatter pie plots
install.packages("scatterpie")

# For HDF5 file support (choose one)
install.packages("hdf5r")
# or
BiocManager::install("rhdf5")

# For Seurat integration
install.packages("SeuratObject")
```

## Quick Start

### Basic Deconvolution

```r
library(STRIDER)

# Run deconvolution pipeline
result <- run_strider(
  sc.data = sc_matrix,       # Single-cell count matrix (genes x cells)
  st.data = st_matrix,       # Spatial count matrix (genes x spots)
  sc.celltype = celltypes,   # Named vector of cell type annotations
  n.topics = NULL,           # Auto-select optimal topic number
  out.dir = "strider_output"
)

# View results
print(result)
summary(result)
```

### From Seurat Objects

```r
# From Seurat objects (v4 or v5)
result <- run_strider(
  sc.data = seurat_sc,
  st.data = seurat_st,
  sc.celltype = "celltype"   # Column name in metadata
)
```

### Visualization

```r
# Scatter pie plot
plot_deconv(result, coords, plot.type = "scatterpie")

# Dominant cell type scatter plot
plot_deconv(result, coords, plot.type = "scatter")

# Cell type composition heatmap
plot_deconv(result, coords, plot.type = "heatmap")
```

### Spatial Clustering

```r
# Neighborhood-aware clustering
clusters <- spatial_cluster(
  result,
  coords,
  n.clusters = 5,
  weight = 0.5    # Balance between spot and neighbor composition
)

# Visualize clusters
plot_cluster(clusters)

# Find enriched cell types per cluster
markers <- find_cluster_markers(result, clusters)
```

### Multi-Sample Integration

```r
# Integrate multiple samples using FGW
integrated <- integrate_samples(
  deconv.list = list(deconv1, deconv2, deconv3),
  topic.list = list(topic1, topic2, topic3),
  coords.list = list(coords1, coords2, coords3),
  alpha = 0.5    # Trade-off: 0=feature, 1=structure
)

# Visualize integration
plot_integration(integrated)
```

## Key Functions

| Function | Description |
|----------|-------------|
| `run_strider()` | Main deconvolution pipeline |
| `train_lda()` | Train LDA topic model |
| `deconvolve()` | Deconvolve spatial spots |
| `spatial_cluster()` | Neighborhood-aware clustering |
| `integrate_samples()` | Multi-sample integration |
| `find_markers()` | Identify cell type markers |
| `plot_deconv()` | Visualize deconvolution |
| `plot_cluster()` | Visualize spatial clusters |

## Input Formats

STRIDER supports multiple input formats:

- **Matrices**: Dense or sparse gene x cell/spot matrices
- **Seurat objects**: v4 and v5 compatible
- **10X HDF5**: Cell Ranger output files
- **Text files**: Tab-separated count matrices

## Performance

The package is optimized for performance:

- Sparse matrix operations throughout
- C++ implementation of FGW optimal transport
- Efficient data I/O via data.table
- Optional parallel processing

## Citation

If you use STRIDER in your research, please cite:

```
Liu, Z. (2026). STRIDER: Spatial Transcriptomics Deconvolution and Integration in R.
GitHub repository: https://github.com/Zaoqu-Liu/STRIDER
```

## License

MIT License

## Contact

- **Author**: Zaoqu Liu
- **Email**: liuzaoqu@163.com
- **GitHub**: https://github.com/Zaoqu-Liu/STRIDER

## Acknowledgments

STRIDER is inspired by the STRIDE Python package (Sun et al., 2022). The Fused Gromov-Wasserstein implementation follows Vayer et al. (2020).
