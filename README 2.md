# SCIBER: Single-Cell Imputation Benchmark for Epigenomics Research

**Welcome to SCIBER**, a resource and pipeline designed for benchmarking imputation methods for single-cell epigenomic data. This repository hosts the code, pipelines, and datasets used in our research to evaluate the applicability of imputation algorithms originally developed for single-cell RNA (scRNA) and single-cell ATAC (scATAC) sequencing data to single-cell histone post-translational modification (scHPTM) data.

## Overview

In the last decade, there have been significant advances in the field of profiling histone modifications at the single-cell level (scHPTM). Similar to single-cell RNA profiling, scHPTM techniques hold the potential to reveal cell-to-cell heterogeneity among various cell subpopulations. However, scHPTM data is often plagued by high levels of noise due to non-specific antibody binding and other technical challenges. Additionally, scHPTM data is characterized by its sparsity, typically containing only a few hundred to a few thousand reads per cell.

This issue of noisy and sparse data is not unique to scHPTM. In the more mature fields of scRNA and scATAC data analysis, several algorithms have been developed to computationally mitigate noise and impute dropouts. While most of these algorithms were initially developed with scRNA data in mind, recent years have seen the adaptation of these methods for scATAC data. To our knowledge, no imputation methods have been developed explicitly for scHPTM data, but we hypothesize that some of the existing methods from scRNA and scATAC may be applicable to scHPTM data.

## Key Features

- **Imputation Benchmarking**: Evaluate and compare the performance of various imputation methods on scHPTM data.
- **Pipeline Automation**: Streamlined pipelines for processing scHPTM datasets, running imputation algorithms, and generating comparative analyses.
- **Reproducibility**: Full reproducibility of all analyses performed in our study, including pre-processing, imputation, and downstream analysis.
- **Extensibility**: Easily extend the framework to include new datasets or additional imputation methods.

## Repository Structure

- **`data/`**: Contains sample scHPTM datasets used in our benchmarks.
- **`scripts/`**: Includes scripts for data pre-processing, running imputation algorithms, and post-imputation analysis.
- **`pipelines/`**: Workflow files for automating the benchmarking process using various computational resources.
- **`results/`**: Output directory where benchmark results and visualizations are stored.
- **`notebooks/`**: Jupyter and Rstudio notebooks for exploratory data analysis and result interpretation.

## Getting Started

### Prerequisites

To get started with SCIBER, you'll need:

- R (version 4.4.1 )
- Python (version X.X.X or later)
- The following R libraries: `proxyC`, `Matrix`, `tidyverse`, `ComplexHeatmap`, and `reshape2`
- Additional libraries listed in `requirements.txt`

### Installation

Clone the SCIBER repository to your local machine:

```bash
git clone https://github.com/robinweide/SCIBER.git
cd SCIBER
```

Install the required R packages:

```r
install.packages(c("proxyC", "Matrix", "tidyverse", "ComplexHeatmap", "reshape2"))
```

Install the required Python packages:

```bash
pip install -r requirements.txt
```

## Contributing

We welcome contributions from the community! If you'd like to contribute, please fork the repository and submit a pull request. For major changes, please open an issue to discuss your proposed changes.

## Citation

If you use SCIBER in your research, please cite our paper:

> Your citation details here.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

We would like to thank all contributors and collaborators who have made this project possible. Especially Peter Zeller for assistance with the SortChiC dataset and the XXX lab for fruitful discussions.
