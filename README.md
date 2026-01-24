# STRIDER <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/Zaoqu-Liu/STRIDER/workflows/R-CMD-check/badge.svg)](https://github.com/Zaoqu-Liu/STRIDER/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/STRIDER)](https://CRAN.R-project.org/package=STRIDER)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

**STRIDER** (**S**patial **TR**anscriptomics **I**ntegration and **DE**convolution in **R**) is a comprehensive R package for cell-type deconvolution of spatial transcriptomics data using topic modeling (LDA).

## 🌟 Features

- **Cell-type Deconvolution**: Decompose cell-type proportions in spatial spots using single-cell reference
- **Multi-sample Integration**: Align multiple samples using Fused Gromov-Wasserstein optimal transport
- **Spatial Clustering**: Neighborhood-aware clustering based on cell-type composition
- **Spot-to-Cell Mapping**: Identify most similar single cells for each spatial spot
- **Publication-ready Visualization**: Beautiful scatter pie plots, heatmaps, and more

## 📦 Installation

### From GitHub (recommended)

```r
# Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Install STRIDER
devtools::install_github("Zaoqu-Liu/STRIDER")
```

### Dependencies

STRIDER requires the following packages:

```r
# Core dependencies
install.packages(c("Matrix", "data.table", "ggplot2", "Rcpp", "RcppArmadillo"))

# For visualization
install.packages(c("scatterpie", "viridis", "patchwork"))

# For LDA modeling
install.packages("text2vec")

# Optional: For reading HDF5 files
install.packages("hdf5r")
# or
BiocManager::install("rhdf5")
```

## 🚀 Quick Start

### Basic Usage

```r
library(STRIDER)

# Run STRIDER deconvolution
result <- run_strider(
  sc.data = "sc_count.h5",           # Single-cell count matrix
  st.data = "st_count.h5",           # Spatial count matrix
  sc.celltype = "celltype.txt",      # Cell type annotations
  out.dir = "output",                 # Output directory
  out.prefix = "my_sample"            # Output file prefix
)

# Visualize results
plot_deconv(result, coords, plot.type = "scatterpie")
```

### From Seurat Objects (V4 Compatible)

```r
library(STRIDER)
library(Seurat)

# Load Seurat objects
seurat_sc <- readRDS("sc_seurat.rds")
seurat_st <- readRDS("st_seurat.rds")

# Run STRIDER
result <- run_strider(
  sc.data = seurat_sc,
  st.data = seurat_st,
  sc.celltype = "celltype"  # Column name in metadata
)
```

## 📖 Workflow

### 1. Deconvolution

```r
# Full pipeline
result <- run_strider(
  sc.data = sc_matrix,
  st.data = st_matrix,
  sc.celltype = celltypes,
  n.markers = 200,        # Markers per cell type
  n.topics = 5:15,        # Topic range to test
  method = "auto"         # Auto-select best method
)

# Access results
deconv_frac <- result$deconv       # Cell type fractions
topic_dist <- result$topic_spot    # Topic distributions
```

### 2. Spatial Clustering

```r
# Neighborhood-aware clustering
clusters <- spatial_cluster(
  result, 
  coords,
  n.clusters = 5,
  weight = 0.5  # Balance between own vs neighbor composition
)

# Visualize
plot_cluster(clusters, pt.size = 4)
```

### 3. Multi-sample Integration

```r
# Integrate multiple samples using FGW optimal transport
integrated <- integrate_samples(
  deconv.list = list(deconv1, deconv2, deconv3),
  topic.list = list(topic1, topic2, topic3),
  coords.list = list(coords1, coords2, coords3),
  alpha = 0.5  # Balance structure vs features
)

# Visualize
plot_integration(integrated)
```

### 4. Spot-to-Cell Mapping

```r
# Map spots to similar cells
mapping <- map_spots_to_cells(
  topic.spot = topic_spot,
  topic.cell = topic_cell,
  n.top = 10
)

# Transfer labels
labels <- transfer_labels(mapping, celltypes)
```

## 📊 Visualization

```r
# Scatter pie plot
plot_deconv_pie(result$deconv, coords, pt.size = 8)

# Scatter plot (dominant cell type)
plot_deconv_scatter(result$deconv, coords, pt.size = 2)

# Heatmap
plot_celltype_heatmap(result$deconv)

# Model selection
plot_model_selection(result$metrics)
```

## ⚡ Performance

STRIDER is optimized for performance:

- **Vectorized operations** using Matrix package
- **C++ backend** (RcppArmadillo) for Fused Gromov-Wasserstein
- **Efficient sparse matrix** handling throughout
- **Optional parallelization** via future package

## 📚 Vignettes

- [Getting Started](vignettes/getting_started.html)
- [Multi-sample Integration](vignettes/integration.html)
- [Working with Seurat](vignettes/seurat_workflow.html)

## 🔬 Citation

If you use STRIDER in your research, please cite:

```
Liu Z, et al. (2024). STRIDER: Spatial Transcriptomics Integration and 
Deconvolution in R. GitHub: https://github.com/Zaoqu-Liu/STRIDER
```

Original Python implementation (STRIDE):
```
Sun D, et al. (2022). STRIDE: accurately decomposing and integrating 
spatial transcriptomics using single-cell RNA sequencing. 
Nucleic Acids Research. doi: 10.1093/nar/gkac150
```

## 📧 Contact

- **Author**: Zaoqu Liu
- **Email**: liuzaoqu@163.com
- **GitHub**: [https://github.com/Zaoqu-Liu/STRIDER](https://github.com/Zaoqu-Liu/STRIDER)

## 📄 License

MIT License © 2024 Zaoqu Liu

## 🙏 Acknowledgments

This package is an R implementation of the original Python STRIDE tool by Dongqing Sun.
Special thanks to the original authors for their excellent work.
