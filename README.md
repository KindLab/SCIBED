
# SCIBED: Single Cell Imputation Benchmarking on Epigenomics Datasets

<!-- badges: start -->
<!-- badges: end -->

**Welcome to SCIBED**, a resource and pipeline designed for benchmarking imputation methods for single-cell epigenomic data. This repository hosts the package used in our research to evaluate the applicability of imputation algorithms originally developed for single-cell RNA (scRNA) and single-cell ATAC (scATAC) sequencing data to single-cell histone post-translational modification (scHPTM) data.

## Overview

In the last decade, there have been significant advances in the field of profiling histone modifications at the single-cell level (scHPTM). Similar to single-cell RNA profiling, scHPTM techniques hold the potential to reveal cell-to-cell heterogeneity among various cell subpopulations. However, scHPTM data is often plagued by high levels of noise due to non-specific antibody binding and other technical challenges. Additionally, scHPTM data is characterized by its sparsity, typically containing only a few hundred to a few thousand reads per cell.

This issue of noisy and sparse data is not unique to scHPTM. In the more mature fields of scRNA and scATAC data analysis, several algorithms have been developed to computationally mitigate noise and impute dropouts. While most of these algorithms were initially developed with scRNA data in mind, recent years have seen the adaptation of these methods for scATAC data. To our knowledge, no imputation methods have been developed explicitly for scHPTM data, but we hypothesize that some of the existing methods from scRNA and scATAC may be applicable to scHPTM data.

## Getting Started

### Quickest option: Docker 🐋

A ready-built Docker image lives here: [docker hub](docker.io/robinhweide/scibed-scibed-rstudio).

### Prerequisites

To get started with SCIBED, you'll need R (>= 3.5.0) with the following packages:

- GenomicRanges
- GenomeInfoDb
- IRanges
- Matrix
- Seurat

## Installation

You can install the development version of SCIBED like so:

``` r
remotes::install_github("kindlab/SCIBED")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SCIBED)
## basic example code
```

## Contributing

We welcome contributions from the community! If you'd like to contribute, please fork the repository and submit a pull request. For major changes, please open an issue to discuss your proposed changes.

## Citation

If you use SCIBED in your research, please cite our paper:

> Your citation details here.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

We would like to thank all contributors and collaborators who have made this project possible. Especially Peter Zeller for assistance with the SortChiC datasets and the Colomé-Tatché lab for fruitful discussions.
