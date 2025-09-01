## Issues

### 1. Add specific noise vector ✅

> Include tutorials for custom noise models (e.g., open-chromatin bias) or 
broader histone marks. --- R2.m.3

### 2. Add second example with H3K9me3 ✅

> Include tutorials for custom noise models (e.g., open-chromatin bias) or 
broader histone marks. --- R2.m.3

### 3. Expand tutorials ✅

> The example provided in the GitHub repository 
(https://github.com/KindLab/SCIBED) should be expanded to include a 
more complete usage workflow.

### 4. Reproducibility ✅

> Please specify the required package versions for dependencies such as 
GenomeInfoDb, IRanges, and GenomicRanges. These packages require 
updated versions of Bioconductor, which in turn depend on R ≥ 4.3. 
Resolving these dependencies can be nontrivial, so providing a version-
controlled environment—ideally as a Docker container with all required 
packages pre-installed—would greatly improve reproducibility and ease of 
use for other researchers --- R1.m.6

## Changelog

### 20250828

1. The `DESCIPTION` file in this repository already has all (versions of) dependencies listed. [4]
2. A dockerfile has been added/ [4]
3. Added `noise_vector` parameter to `generate_in_silico()` [1]
4. Added explainer on the `noise_vector` parameter of `generate_in_silico()` to the SCIBED_tutorial [1]

### 20250829

1. Docker image was built and uploaded to Docker hub [4]
2. Added H3K9me3-datasets and referred to it in the tutorial [2]

### 20250901

1. Added second vignette (ranking_vignette) to show complete workflow using H3K9me3 [2,3]
