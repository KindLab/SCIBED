#' H3K9me3 sortChIC data
#'
#' These data sets contain a subset of the cells from the H3K9me3 dataset of Zeller et al. (2023)
#'
#' Contains a sparse matrix of 100 B-cells and 100 Erythrocytes, binned at 50kb.
#'
#' @format A bins by cells sparse matrix
#' @examples
#' test_matrix <- zeller_H3K9me3_matrix
#' test_matrix_clusters <- zeller_H3K9me3_metadata
"zeller_H3K9me3_matrix"


#' H3K9me3 sortChIC metadata
#'
#' These datasets contain a subset of the cells from the H3K9me3 dataset of Zeller et al. (2023)
#'
#' Contains a vector of the FACS-based celltype per column in `zeller_H3K9me3_matrix`.
#'
#' @format A vector of `length(ncol(zeller_H3K9me3_matrix))`
#' @examples
#' test_matrix <- zeller_H3K9me3_matrix
#' test_matrix_clusters <- zeller_H3K9me3_metadata
"zeller_H3K9me3_metadata"

#' H3K9me3 sortChIC peaks
#'
#' These datasets contain filtered peaks for two cell types from the H3K9me3 dataset of Zeller et al. (2023)
#'
#' Contains a `GRangesList` for peaks from the Bcells and Eryths clusters
#'
#' @format A `GRangesList` object with two entries
#' @examples
#' zeller_H3K9me3_peaks <- zeller_H3K9me3_peaks
"zeller_H3K9me3_peaks"



