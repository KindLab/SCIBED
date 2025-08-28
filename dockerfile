# R + RStudio Server (based on R 4.4.1)
FROM rocker/rstudio:4.4.1

# Install system dependencies
RUN apt-get update && apt-get install -y \
libxml2-dev \
libssl-dev \
libcurl4-openssl-dev \
libfontconfig1-dev \
libharfbuzz-dev \
libfribidi-dev \
libfreetype6-dev \
libpng-dev \
libtiff5-dev \
libjpeg-dev \
&& rm -rf /var/lib/apt/lists/*

  # Install pak + BiocManager
  RUN R -e "install.packages(c('pak','BiocManager'), repos='https://cloud.r-project.org')"

# Pre-install Bioconductor deps (from DESCRIPTION)
RUN R -e "BiocManager::install(c('GenomicRanges','IRanges','GenomeInfoDb','S4Vectors'), ask=FALSE, update=TRUE)"

# Pre-install heavy CRAN deps
RUN R -e "pak::pak(c('Matrix','Seurat','proxyC','RcppAnnoy','sparseMatrixStats'))"

# Install SCIBED from GitHub
RUN R -e "pak::pak('KindLab/SCIBED')"

# Set working directory (mounted from host if using docker-compose)
WORKDIR /home/rstudio

# Default CMD is RStudio (from rocker/rstudio)
