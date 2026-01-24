# STRIDER

<!-- badges: start -->
[![R-universe](https://zaoqu-liu.r-universe.dev/badges/STRIDER)](https://zaoqu-liu.r-universe.dev/STRIDER)
[![R-CMD-check](https://github.com/Zaoqu-Liu/STRIDER/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Zaoqu-Liu/STRIDER/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub issues](https://img.shields.io/github/issues/Zaoqu-Liu/STRIDER)](https://github.com/Zaoqu-Liu/STRIDER/issues)
<!-- badges: end -->

## Overview

**STRIDER** (**S**patial **T**ranscriptomics deconvolu**T**ion and integ**R**ation **I**n **D**eveloped **E**nvironment for **R**) is a comprehensive R package for the analysis of spatial transcriptomics data. It provides a unified framework for cell type deconvolution, multi-sample integration, and spatial domain identification.

Spatial transcriptomics technologies enable the measurement of gene expression while preserving spatial context, but most platforms capture mixtures of cells within each measurement spot. STRIDER addresses this challenge by leveraging topic modeling approaches to decompose spatial spots into constituent cell type proportions, enabling downstream biological interpretation of tissue microenvironments.

## Key Features

- **Cell Type Deconvolution**: Latent Dirichlet Allocation (LDA)-based decomposition of spatial spots into cell type proportions using single-cell RNA-seq reference data
- **Multi-Sample Integration**: Fused Gromov-Wasserstein (FGW) optimal transport for alignment and integration of multiple spatial transcriptomics samples
- **Spatial Domain Identification**: Neighborhood-aware clustering that incorporates both cell type composition and spatial proximity
- **Spot-Cell Mapping**: Identification of transcriptionally similar single cells for each spatial spot
- **Visualization**: Publication-ready visualizations including scatter pie plots, spatial heatmaps, and integration plots

## Installation

### From R-universe (Recommended)

```r
install.packages("STRIDER", repos = "https://zaoqu-liu.r-universe.dev")
```

### From GitHub

```r
# Install remotes if needed
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

# Install STRIDER
remotes::install_github("Zaoqu-Liu/STRIDER")
```

### Dependencies

STRIDER requires several packages for full functionality:

```r
# Core dependencies
install.packages(c("Matrix", "data.table", "ggplot2", "viridis", "Rcpp", "RcppArmadillo"))

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

# For Seurat integration (optional)
install.packages("SeuratObject")
```

## Methodology

### Cell Type Deconvolution

STRIDER employs a topic modeling framework where:
- **Documents** correspond to spatial spots
- **Words** correspond to genes
- **Topics** represent latent cell type signatures

The deconvolution process consists of:

1. **Model Training**: An LDA model is trained on single-cell reference data to learn topic-gene distributions (φ) and cell-topic distributions (θ)
2. **Topic-Celltype Association**: Topic distributions are aggregated across cells of each type to establish topic-celltype relationships
3. **Inference**: The trained model infers topic distributions for spatial spots
4. **Deconvolution**: Cell type fractions are computed using Bayesian inference:

$$P(\text{celltype}|\text{topic}) = \frac{P(\text{topic}|\text{celltype}) \cdot P(\text{celltype})}{P(\text{topic})}$$

### Multi-Sample Integration

For integrating multiple spatial samples, STRIDER implements Fused Gromov-Wasserstein (FGW) optimal transport, which jointly optimizes:

- **Feature similarity**: KL divergence between topic distributions
- **Structural similarity**: Gromov-Wasserstein distance preserving spatial relationships

The FGW distance is defined as:

$$\text{FGW}(\mu, \nu) = \min_{\pi} (1-\alpha) \langle M, \pi \rangle + \alpha \cdot \text{GW}(C_1, C_2, \pi)$$

where α balances feature and structural costs, and the optimization is solved via conditional gradient descent with Sinkhorn projections.

## Quick Start

### Basic Deconvolution Pipeline

```r
library(STRIDER)

# Run the complete deconvolution pipeline
result <- run_strider(
  sc.data = sc_matrix,       # Single-cell count matrix (genes × cells)
  st.data = st_matrix,       # Spatial count matrix (genes × spots)
  sc.celltype = celltypes,   # Named vector of cell type annotations
  n.topics = NULL,           # Auto-select optimal number of topics
  method = "auto",           # Automatic method selection
  out.dir = "strider_output"
)

# Examine results
print(result)
summary(result)
```

### Working with Seurat Objects

```r
# Direct input from Seurat objects (v4 and v5 compatible)
result <- run_strider(
  sc.data = seurat_sc,
  st.data = seurat_st,
  sc.celltype = "celltype"   # Column name in metadata
)
```

### Visualization

```r
# Scatter pie plot showing cell type composition
plot_deconv(result, coords, plot.type = "scatterpie")

# Dominant cell type visualization
plot_deconv(result, coords, plot.type = "scatter")

# Cell type composition heatmap
plot_celltype_heatmap(result)
```

### Spatial Domain Identification

```r
# Neighborhood-aware spatial clustering
clusters <- spatial_cluster(
  result,
  coords,
  n.clusters = 5,
  weight = 0.5,    # Balance between spot and neighbor composition
  k = 6            # Number of spatial neighbors
)

# Visualize spatial domains
plot_cluster(clusters)

# Identify enriched cell types per domain
enrichment <- find_cluster_markers(result, clusters)
```

### Multi-Sample Integration

```r
# Integrate multiple spatial samples
integrated <- integrate_samples(
  deconv.list = list(result1$deconv, result2$deconv, result3$deconv),
  topic.list = list(result1$topic_spot, result2$topic_spot, result3$topic_spot),
  coords.list = list(coords1, coords2, coords3),
  alpha = 0.5      # Trade-off: 0 = feature only, 1 = structure only
)

# Visualize aligned samples
plot_integration(integrated)

# Assess integration quality
quality <- integration_quality(integrated)
```

## Function Reference

| Function | Description |
|----------|-------------|
| `run_strider()` | Complete deconvolution pipeline |
| `train_lda()` | Train LDA topic model on single-cell data |
| `transform_lda()` | Apply trained model to spatial data |
| `deconvolve()` | Compute cell type fractions |
| `spatial_cluster()` | Neighborhood-aware spatial clustering |
| `integrate_samples()` | Multi-sample integration via FGW |
| `map_spots_to_cells()` | Find similar cells for each spot |
| `find_markers()` | Identify cell type marker genes |
| `plot_deconv()` | Visualize deconvolution results |
| `plot_cluster()` | Visualize spatial clusters |
| `plot_integration()` | Visualize integrated samples |

## Input Data Formats

STRIDER supports multiple input formats:

| Format | Description |
|--------|-------------|
| **Matrix** | Dense or sparse gene × cell/spot matrices |
| **Seurat** | Seurat objects (v4 and v5 compatible) |
| **10X HDF5** | Cell Ranger output files (.h5) |
| **Text** | Tab-separated count matrices |

## Performance Considerations

STRIDER is optimized for computational efficiency:

- **Sparse matrix operations** throughout the pipeline
- **C++ implementation** of FGW optimal transport via RcppArmadillo
- **Efficient I/O** using data.table for large files
- **Memory-efficient** sparse matrix storage

Typical runtime for a dataset with 10,000 spots and 20,000 genes: ~2-5 minutes on a standard laptop.

## References

If you use STRIDER in your research, please cite:

```bibtex
@software{liu2026strider,
  author = {Liu, Zaoqu},
  title = {STRIDER: Spatial Transcriptomics Deconvolution and Integration in R},
  year = {2026},
  url = {https://github.com/Zaoqu-Liu/STRIDER}
}
```

STRIDER builds upon methodological foundations from:

- **STRIDE**: Sun, D., et al. (2022). STRIDE: accurately decomposing and integrating spatial transcriptomics using single-cell RNA sequencing. *Nucleic Acids Research*, 50(7), e42.
- **Fused Gromov-Wasserstein**: Vayer, T., et al. (2020). Fused Gromov-Wasserstein distance for structured objects. *Algorithms*, 13(9), 212.
- **LDA**: Blei, D.M., et al. (2003). Latent Dirichlet Allocation. *Journal of Machine Learning Research*, 3, 993-1022.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

- **Author**: Zaoqu Liu
- **Email**: liuzaoqu@163.com
- **GitHub**: [https://github.com/Zaoqu-Liu/STRIDER](https://github.com/Zaoqu-Liu/STRIDER)
- **Issues**: [https://github.com/Zaoqu-Liu/STRIDER/issues](https://github.com/Zaoqu-Liu/STRIDER/issues)

## Acknowledgments

We thank the developers of the original STRIDE Python package for their pioneering work in spatial transcriptomics deconvolution. The optimal transport implementations are based on the Python Optimal Transport (POT) library.
