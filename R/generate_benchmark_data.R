# tools to transform a bin by cell matrix to a set of matrices with different levels of noise and depth

################################### internals ##################################


#' Calculate Fraction of Zero-Bins
#'
#' Computes the fraction of elements in a numeric vector that are zero.
#'
#' @param x A numeric vector for which the fraction of zero elements is calculated.
#'
#' @return A numeric value between 0 and 1 representing the fraction of zero values in `x`.
#' @examples
#' \dontrun{
#' calculate_fzc(c(0, 1, 0, 2, 3))
#' calculate_fzc(c(0, 0, 0, 0))
#' calculate_fzc(c(1, 2, 3, 4))
#' }
#'
calculate_fzc <- function(x) {
  if (!is.vector(x)) {
    stop("`x` must be a vector.")
  }
  if (!is.numeric(x)) {
    stop("`x` must be numeric.")
  }

  mean(x == 0)
}

#' Filter Matrix by Fraction of Zero Counts (FZC)
#'
#' Filters columns of a sparse matrix based on the fraction of zero values in each column.
#' Columns with an FZC below a specified threshold (mean FZC minus `sd_magnitude` standard deviations)
#' are excluded from the returned matrix.
#'
#' @param x A sparse matrix of class `Matrix` or `dgCMatrix`.
#' @param sd_magnitude A numeric value indicating the number of standard deviations below the mean FZC
#' used to set the threshold. Default is 2.
#'
#' @return A sparse matrix containing only the columns with an FZC above the calculated threshold.
#' @importFrom stats sd
#' @examples
#' \dontrun{
#' # Assuming `sparse_mat` is a dgCMatrix
#' filter_by_fzc(sparse_mat, sd_magnitude = 2)
#' }
#'
filter_by_fzc <- function(x, sd_magnitude = 2) {
  if (!inherits(x, "dgCMatrix")) {
    stop("`x` must be a sparse matrix of class `dgCMatrix`.")
  }

  fzc_data <- apply(x, MARGIN = 2, calculate_fzc)
  fzc_threshold <- mean(fzc_data) - (sd_magnitude * stats::sd(fzc_data))

  x[, fzc_data >= fzc_threshold]
}


#' Add Noise to Sparse Matrix
#'
#' Adds noise to a sparse matrix by increasing the original noise level by a specified factor.
#'
#' @param x A sparse matrix of class `dgCMatrix`.
#' @param increase_factor A numeric value indicating the factor by which to increase the original noise level. Default is 0.
#'
#' @return A sparse matrix with added noise.
#' @examples
#' \dontrun{
#' # Assuming `sparse_mat` is a dgCMatrix
#' increase_noise(sparse_mat, increase_factor = 2)
#'}
#'
increase_noise <- function(x, increase_factor = 0) {
  if (!inherits(x, "dgCMatrix")) {
    stop("`x` must be a sparse matrix of class `dgCMatrix`.")
  }

  orig_noise <- sum(x) / (nrow(x) * ncol(x))

  # round is there to get to whole integers again
  output <- x + (round(orig_noise) * increase_factor)

  output
}


#' Calculate Pseudobulk from Sparse Matrix
#'
#' Aggregates columns of a sparse matrix into pseudobulk profiles based on a grouping vector.
#'
#' @param x A sparse matrix of class `dgCMatrix`.
#' @param grouping_vector A vector indicating group assignments for each column in `x`. If `NULL`,
#' all columns will be grouped into a single pseudobulk profile.
#'
#' @return A sparse matrix of class `dgCMatrix` where each column represents the pseudobulk profile
#' of a group specified in `grouping_vector`.
#' @import Matrix
#' @examples
#' \dontrun{
#' # Assuming `sparse_mat` is a dgCMatrix and `groups` is a vector of group labels
#' calculate_pseudobulk(sparse_mat, grouping_vector = groups)
#' }
calculate_pseudobulk <- function(x, grouping_vector = NULL) {
  if (!inherits(x, "dgCMatrix")) {
    stop("`x` must be a sparse matrix of class `dgCMatrix`.")
  }

  if (!is.null(grouping_vector)) {
    if (!is.vector(grouping_vector)) {
      stop("`grouping_vector` must be a vector.")
    }
    if (length(grouping_vector) != ncol(x)) {
      stop("`grouping_vector` must have the same length as the number of columns in `x`.")
    }
  } else {
    grouping_vector <- rep("group", ncol(x))
  }

  pseudobulks <- Matrix(nrow = nrow(x), ncol = 0, sparse = TRUE)

  for (group_name in unique(grouping_vector)) {
    x_sub <- x[, grouping_vector == group_name]
    x_ps <- Matrix(rowSums(x_sub), ncol = 1, dimnames = list(rownames(x), group_name), sparse = TRUE)
    pseudobulks <- cbind(pseudobulks, x_ps)
  }

  pseudobulks
}

#' Fill Missing Bins with Zeros in a Sparse Matrix
#'
#' Expands a sparse matrix to include all expected genomic bins, filling missing bins with zeros.
#' Handles cases where the last bin of each chromosome may be truncated, and detects 0-based or 1-based
#' binning formats automatically.
#'
#' @param mat A sparse matrix of class `dgCMatrix` with genomic bins as row names in the format `chr:start-end`.
#' @param bin_size An integer indicating the expected bin size in base pairs.
#'
#' @return A sparse matrix of class `dgCMatrix` with all expected bins, including zeros for missing bins.
#' @importFrom utils strcapture
#' @import Matrix
#' @examples
#' \dontrun{
#' # Assuming `sparse_mat` is a dgCMatrix with row names like "chr1:0-50000" or "chr1:1-50000"
#' fill_missing_bins(sparse_mat, bin_size = 50000)
#' }
#'
fill_missing_bins <- function(mat, bin_size) {
  # Extract chromosome, start, and end positions without assuming "chr" prefix
  bin_info <- utils::strcapture(
    pattern = "(.*):(\\d+)-(\\d+)",
    rownames(mat),
    data.frame(chr = character(), start = integer(), end = integer())
  )

  # Determine if the matrix is 1-based by checking the minimum start position
  is_1_based <- any(bin_info$start == 1)

  # Convert to 0-based if necessary for internal processing
  if (is_1_based) {
    bin_info$start <- bin_info$start - 1
  }

  # Create the full range of expected bins for each chromosome
  full_bins <- do.call(rbind, lapply(unique(bin_info$chr), function(chr) {
    chr_bins <- bin_info[bin_info$chr == chr, ]
    max_pos <- max(as.numeric(chr_bins$end))
    seq_start <- seq(0, max_pos, by = bin_size)

    # Adjust last bin to end at max_pos if truncated
    end_positions <- ifelse(seq_start + bin_size > max_pos, max_pos, seq_start + bin_size)

    data.frame(
      chr = chr,
      start = seq_start,
      end = end_positions
    )
  }))

  # Convert full_bins to bin names, adjusting back to 1-based if necessary
  full_bins$bin_name <- if (is_1_based) {
    with(full_bins, sprintf("%s:%d-%d", chr, start + 1, end))
  } else {
    with(full_bins, sprintf("%s:%d-%d", chr, start, end))
  }

  # Identify and create zero-filled rows for missing bins
  missing_bins <- setdiff(full_bins$bin_name, rownames(mat))
  zero_matrix <- Matrix(0, nrow = length(missing_bins), ncol = ncol(mat), sparse = TRUE)
  rownames(zero_matrix) <- missing_bins

  # Combine original matrix with zero-filled matrix for missing bins
  complete_mat <- rbind(mat, zero_matrix)
  sorted_mat <- complete_mat[order(match(rownames(complete_mat), full_bins$bin_name)), , drop = FALSE]

  sorted_mat
}



################################### externals ##################################

#' Generate Ground Truth Pseudobulks from Sparse Matrix
#'
#' This function processes a sparse matrix to generate ground truth pseudobulks,
#' filtering by sparsity and minimum read depth. It expands the matrix to include
#' all expected genomic bins and calculates pseudobulks based on a provided grouping vector.
#'
#' @param x A sparse matrix of class `dgCMatrix`.
#' @param grouping_vector A vector indicating the grouping for each column in `x`.
#' @param bin_size An integer specifying the expected bin size in base pairs.
#' @param min_RD An integer specifying the minimum read depth required for inclusion.
#' @param fzc_sd An optional numeric value to filter by fraction of zero bins (FZC).
#'
#' @return A sparse matrix of class `dgCMatrix` containing the ground truth pseudobulks.
#' @examples
#' data(zeller_H3K4me3_matrix)
#' data(zeller_H3K4me3_metadata)
#'
#' zeller_ground_truth <- generate_ground_truth(
#'   zeller_H3K4me3_matrix,
#'   grouping_vector = zeller_H3K4me3_metadata,
#'   bin_size = 50000,
#'   min_RD = 10
#'   )
#'
#' @export
generate_ground_truth <- function(x, grouping_vector, bin_size, min_RD = 0, fzc_sd = NULL) {
  if (!inherits(x, "dgCMatrix")) {
    stop("`x` must be a sparse matrix of class `dgCMatrix`.")
  }

  if (!is.vector(grouping_vector)) {
    stop("`grouping_vector` must be a vector.")
  }

  if (length(grouping_vector) != ncol(x)) {
    stop("`grouping_vector` must have the same length as the number of columns in `x`.")
  }

  # Check if row names have expected chr:start-end notation
  row_names <- rownames(x)

  if (is.null(row_names) || any(!grepl("(.*):(\\d+)-(\\d+)", row_names))) {
    stop("Row names must follow the format 'chr:start-end'.")
  }

  # Filter by sparsity (fraction of bins with 0)
  mat <- x
  if (!is.null(fzc_sd)) {
    mat <- filter_by_fzc(mat, sd_magnitude = fzc_sd)
  }

  # Filter by read depth
  mat <- mat[, colSums(mat) >= min_RD]

  # Expand matrix to full genome size of available chromosomes
  mat <- fill_missing_bins(mat, bin_size)

  # Generate ground truth by calculating pseudobulks
  # this is not pretty: best would be if the filterings give us indexes of the good cells, so we can subset the grouping_vector
  ground_truth_pseudobulks <- calculate_pseudobulk(mat, grouping_vector = grouping_vector[colnames(x) %in% colnames(mat)])

  ground_truth_pseudobulks
}



#' Generate In Silico Datasets with Noise and Downsampling
#'
#' This function generates an in silico dataset by creating pseudobulks from a sparse matrix,
#' optionally adding noise and downsampling to achieve a target read depth.
#'
#' @param x A sparse matrix of class `dgCMatrix`.
#' @param grouping_vector A vector indicating the grouping for each column in `x`.
#' @param bin_size An integer specifying the expected bin size in base pairs.
#' @param target_RD A numeric value specifying the target read depth for downsampling.
#' @param noise_level A numeric value indicating the noise level to be added; defaults to 0.
#' @param min_RD An integer specifying the minimum read depth required for inclusion.
#' @param fzc_sd An optional numeric value for filtering by fraction of zero bins (FZC).
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{ground_truth}{A sparse matrix of class `dgCMatrix` containing the ground truth pseudobulks.}
#'     \item{in_silico_matrix}{A sparse matrix with added noise and downsampling, representing the in silico dataset.}
#'   }
#'
#' @importFrom Seurat SampleUMI
#' @importFrom stats complete.cases
#' @examples
#' data(zeller_H3K4me3_matrix)
#' data(zeller_H3K4me3_metadata)
#'
#' in_silico_output <- generate_in_silico(zeller_H3K4me3_matrix,
#'   grouping_vector = zeller_H3K4me3_metadata,
#'   bin_size = 50000,
#'   target_RD = 1000,
#'   noise_level = 0.1)
#'
#' @export
generate_in_silico <- function(x, grouping_vector, bin_size, target_RD, noise_level = 0, min_RD = 0, fzc_sd = NULL) {

  # Check that grouping_vector contains no NA values and only strings
  if (!all(stats::complete.cases(grouping_vector))) {
    stop("`grouping_vector` must have no NA values.")
  }

  # Generate ground truth pseudobulks
  ground_truth_pseudobulks <- generate_ground_truth(
    x = x,
    grouping_vector = grouping_vector,
    min_RD = min_RD,
    fzc_sd = fzc_sd,
    bin_size = bin_size
  )

  # Create ground truth matrix with replicates of group columns
  ground_truth_matrix <- do.call(cbind, lapply(unique(grouping_vector), function(group_name) {
    gtmat <- ground_truth_pseudobulks[, rep(group_name, each = sum(grouping_vector == group_name))]
    colnames(gtmat) <- paste0(group_name, "_", seq_along(colnames(gtmat)))
    gtmat
  }))

  # Add noise if specified
  if (noise_level > 0) {
    ground_truth_matrix <- increase_noise(x = ground_truth_matrix, increase_factor = noise_level)
  }

  # Downsample to introduce sparsity
  ground_truth_matrix <- Seurat::SampleUMI(ground_truth_matrix, max.umi = as.numeric(target_RD))

  list(
    ground_truth = ground_truth_pseudobulks,
    in_silico_matrix = ground_truth_matrix
  )
}



