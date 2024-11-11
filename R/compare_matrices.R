

#' Calculate Correlation per cell to Ground Truth
#'
#' This function calculates the correlation between each cell in a matrix and a ground truth pseudobulk,
#' grouped by cell type, using a specified correlation method.
#'
#' @param x A sparse matrix of class `dgCMatrix` with cells to correlate.
#' @param grouping_vector A vector indicating the group (cell type) for each column in `x`.
#' @param bin_size An integer specifying the bin size in base pairs.
#' @param ground_truth A sparse matrix containing the ground truth pseudobulk per cell type.
#' @param method A character string specifying the correlation method, either 'pearson' or 'spearman'; defaults to 'pearson'.
#'
#' @return A data frame with summary statistics of correlation values for each cell type, including
#'   \describe{
#'     \item{group}{The cell type or "complete" for overall statistics.}
#'     \item{cor_method}{The correlation method used.}
#'     \item{mu}{Mean correlation value.}
#'     \item{sd}{Standard deviation of the correlation values.}
#'     \item{N}{Number of cells in the group.}
#'   }
#' @importFrom stats cor
#' @importFrom stats sd
#' @examples
#' data(zeller_H3K4me3_matrix)
#' data(zeller_H3K4me3_metadata)
#'
#' ground_truth_mat <- generate_ground_truth(
#'   x = zeller_H3K4me3_matrix,
#'   grouping_vector = zeller_H3K4me3_metadata,
#'   bin_size = 50000,
#'   fzc_sd = 0)
#'
#' correlation_stats <- calculate_correlation(zeller_H3K4me3_matrix,
#'   grouping_vector = zeller_H3K4me3_metadata,
#'   bin_size = 50000,
#'   ground_truth = ground_truth_mat,
#'   method = "pearson")
#'
#' @export
calculate_correlation <- function(x, grouping_vector, bin_size, ground_truth, method = "pearson") {

  # Check if all expected clusters are in ground truth
  if (any(is.na(match(grouping_vector, colnames(ground_truth))))) {
    stop("Not all groups in `grouping_vector` are found as column names in `ground_truth`.")
  }

  # Expand ground truth and input matrix with missing bins
  data_gt <- fill_missing_bins(ground_truth, bin_size)
  data_imp <- fill_missing_bins(x, bin_size)

  # Ensure rows are sorted to match ground truth
  data_imp <- data_imp[rownames(data_gt), ]

  # Calculate correlation per cell type
  per_cell_stats <- do.call(rbind, lapply(unique(grouping_vector), function(group_name) {
    sub_mat <- data_imp[, grouping_vector == group_name, drop = FALSE]
    cor_vector <- as.vector(cor(method = method, as.matrix(sub_mat), data_gt[, group_name]))
    data.frame(group = group_name, cor_method = method, value = cor_vector)
  }))

  # Overall summary statistics
  overall_stats <- data.frame(
    group = "complete",
    cor_method = method,
    mu = mean(per_cell_stats$value),
    sd = sd(per_cell_stats$value),
    N = nrow(per_cell_stats)
  )

  # Summary statistics per group
  group_stats <- do.call(rbind, lapply(unique(per_cell_stats$group), function(group_name) {
    tmp <- per_cell_stats[per_cell_stats$group == group_name, ]
    data.frame(
      group = group_name,
      cor_method = method,
      mu = mean(tmp$value),
      sd = sd(tmp$value),
      N = nrow(tmp)
    )
  }))

  # Combine overall and group-specific statistics
  result <- rbind(overall_stats, group_stats)
  rownames(result) <- NULL
  return(result)
}

#' Min-Max Scaling
#'
#' This function applies min-max scaling to a numeric vector, rescaling values to a [0, 1] range.
#'
#' @param x A numeric vector to be rescaled.
#' @param na.rm A logical value indicating whether `NA` values should be removed before scaling.
#'   Defaults to `TRUE`.
#'
#' @return A numeric vector with values scaled to the [0, 1] range.
#' @examples
#' minmax(c(1, 2, 3, 4, 5))
#' minmax(c(1, NA, 3, 4, 5), na.rm = TRUE)
#'
#' @export
minmax <- function(x, na.rm = TRUE) {
  (x - min(x, na.rm = na.rm)) / (max(x, na.rm = na.rm) - min(x, na.rm = na.rm))
}

#' Calculate Enrichment for Regions
#'
#' This function calculates enrichment for specified genomic regions, using either a `GRanges`,
#' `GRangesList`, or `data.frame` of regions, and supports optional permutation testing.
#'
#' @param x A sparse matrix representing data to be transformed to pseudobulk for enrichment calculation.
#' @param grouping_vector A vector indicating group assignments for columns in `x`.
#' @param regions A `GRanges`, `GRangesList`, or `data.frame` specifying the regions for enrichment analysis.
#'   If a `GRangesList` is provided, it must be named to match entries in `grouping_vector`.
#' @param maxgap Integer specifying maximum gap between bins (default: -1).
#' @param minoverlap Integer specifying minimum overlap between bins (default: 0).
#' @param perm_times Integer specifying the number of permutations for permutation testing (optional).
#' @param perm_cores Integer specifying the number of cores for parallel computation in permutation testing (optional).
#'
#' @return A data frame with enrichment statistics, including optional permutation test results.
#' @importFrom GenomicRanges GRanges GRangesList reduce
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom S4Vectors mcols `mcols<-`
#' @importFrom methods is
#' @importFrom stats quantile
#'
#' @note This quantification is done on pseudobulks; `calculate_SIMIC` and
#' `calculate_correlation` are single-cell based. In the future, we want
#' to add a single-cell variant of this function.
#'
#' @examples
#' data(zeller_H3K4me3_matrix)
#' data(zeller_H3K4me3_metadata)
#' data(zeller_H3K4me3_peaks)
#'
#' enrichment_stats <- calculate_enrichment(
#'   x = zeller_H3K4me3_matrix,
#'   grouping_vector = zeller_H3K4me3_metadata,
#'   regions = zeller_H3K4me3_peaks)
#'
#' @export
calculate_enrichment <- function(x, grouping_vector, regions, maxgap = -1L, minoverlap = 0L, perm_times = NULL, perm_cores = 2) {

  # Check the format of the `regions` input and convert if necessary
  if (is(regions, "GRangesList")) {
    if (any(is.na(match(grouping_vector, names(regions))))) {
      stop("Names in `regions` (GRangesList) must be found in `grouping_vector`.")
    }
    regions <- regions[unique(grouping_vector)]

  } else if (is(regions, "GRanges")) {
    regions <- GRangesList(input_regions = regions)

  } else if (is(regions, "data.frame")) {
    if (ncol(regions) < 3) {
      stop("The `regions` data frame must have at least three columns.")
    }
    regions <- GRanges(seqnames = regions[, 1], ranges = IRanges(start = regions[, 2], end = regions[, 3]))
    regions <- GRangesList(input_regions = regions)

  } else {
    stop("`regions` must be a `GRanges`, `GRangesList`, or `data.frame`.")
  }

  # Transform input matrix `x` to pseudobulk GRanges list
  pb <- calculate_pseudobulk(x, grouping_vector)
  x_gr <- GRanges(rownames(pb))
  S4Vectors::mcols(x_gr) <- as.data.frame(pb)
  x_depth <- colSums(as.matrix(S4Vectors::mcols(x_gr)))

  # Ensure consistent chromosome naming
  regions <- lapply(regions, function(gr) {
    GenomeInfoDb::seqlevelsStyle(gr) <- GenomeInfoDb::seqlevelsStyle(x_gr)
    gr
  })

  # Calculate enrichment per region
  per_region_stats <- do.call(rbind, lapply(seq_along(regions), function(i) {

    gr <- regions[[i]]
    y <- colSums(as.matrix(S4Vectors::mcols(IRanges::subsetByOverlaps(x_gr, gr, maxgap = maxgap, minoverlap = minoverlap))))
    partial_result <- data.frame(
      group = names(y),
      region = names(regions)[i],
      SIP = y / x_depth[names(y)]
    )

    # Perform permutation testing if requested
    if (!is.null(perm_times)) {
      perm_genome <- GenomicRanges::reduce(x_gr)
      perm_results <- parallel::mclapply(seq_len(perm_times), mc.cores = perm_cores, function(p_i) {
        set.seed(p_i)
        gr_p <- regioneR::resampleGenome(
          A = gr,
          genome = perm_genome,
          min.tile.width = pmax(pmin(1e3, quantile(gr@ranges@width, .1)),100)
        )
        y <- colSums(as.matrix(S4Vectors::mcols(IRanges::subsetByOverlaps(x_gr, gr_p, maxgap = maxgap, minoverlap = minoverlap))))
        data.frame(
          group = names(y),
          SIP = y / x_depth[names(y)]
        )
      })

      # Aggregate permutation results
      perm_out <- do.call(rbind, perm_results)
      perm_stats <- do.call(rbind, lapply(split.data.frame(perm_out, perm_out$group), function(e) {
        data.frame(
          group = e$group[1],
          SIPperm_mu = mean(e$SIP),
          SIPperm_sd = sd(e$SIP)
        )
      }))

      # Merge and calculate Z-scores for enrichment
      partial_result <- merge(partial_result, perm_stats, by = "group", all.x = TRUE)
      partial_result$SIPperm_times <- perm_times
      partial_result$SIPperm_Z <- (partial_result$SIP - partial_result$SIPperm_mu) / partial_result$SIPperm_sd
    }

    partial_result
  }))

  rownames(per_region_stats) <- NULL
  per_region_stats
}

#' Generate Similarity Matrix Using Cosine Distance
#'
#' This function calculates the cosine similarity matrix for pseudobulks.
#'
#' @param x A sparse matrix containing data to be grouped by the `grouping_vector`.
#' @param grouping_vector A vector specifying group assignments for each column (cell) in `x`.
#'
#' @return A similarity matrix calculated using cosine distance on pseudobulk profiles.
#' @importFrom proxyC simil
#' @examples
#' data(zeller_H3K4me3_matrix)
#' data(zeller_H3K4me3_metadata)
#'
#' sim_mat <- generate_ab_initio(
#'   x = zeller_H3K4me3_matrix,
#'   grouping_vector = zeller_H3K4me3_metadata)
#'
#' @export
generate_ab_initio <- function(x, grouping_vector) {

  # Calculate pseudobulk matrix based on the grouping vector
  pb_mat <- calculate_pseudobulk(x = x, grouping_vector = grouping_vector)

  # Compute cosine similarity on the pseudobulk matrix
  similarity_matrix <- proxyC::simil(pb_mat, margin = 2, method = "cosine")

  similarity_matrix
}


#' Calculate SIMIC: Similarity Index for Clustering
#'
#' This function calculates the SIMIC score, which compares observed similarity to expected similarity
#' for clustering cell types. The user can specify quantiles for high-variance bins (HVB) to focus on
#' the most variable features.
#'
#' @param x A sparse matrix of data for similarity calculation.
#' @param grouping_vector A vector indicating group assignments for each cell in `x`.
#' @param ab_initio Optional precomputed similarity matrix. If `NULL`, it will be calculated.
#' @param HVB_q Quantile of most variable bins to consider (e.g., 0.95 means top-5 percent variable bins).
#' @param K Integer, number of nearest neighbors to consider in calculating similarity.
#' @param N_cores Integer, number of cores for parallel processing.
#' @return A data frame with observed similarity, expected similarity, and their ratio (delta).
#' @importFrom proxyC simil
#' @importFrom sparseMatrixStats rowVars
#' @importFrom parallel mclapply
#' @importFrom methods new
#' @importFrom stats setNames
#' @importFrom RcppAnnoy AnnoyEuclidean
#' @examples
#' data(zeller_H3K4me3_matrix)
#' data(zeller_H3K4me3_metadata)
#'
#' sim_mat <- generate_ab_initio(
#'   x = zeller_H3K4me3_matrix,
#'   grouping_vector = zeller_H3K4me3_metadata)
#'
#' simicstats <- calculate_simic(
#'   x = zeller_H3K4me3_matrix,
#'   grouping_vector = zeller_H3K4me3_metadata,
#'   ab_initio = sim_mat,
#'   HVB_q = 0.95,
#'   K = 10,
#'   N_cores = 2)
#'
#' @export
calculate_simic <- function(x, grouping_vector, ab_initio = NULL, HVB_q = 0.95, K = 10, N_cores = 2) {

  # Calculate the ab initio similarity matrix if not provided
  if (is.null(ab_initio)) {
    ab_initio <- generate_ab_initio(x = x, grouping_vector = grouping_vector)
  }
  similarity_matrix <- minmax(ab_initio)

  # Calculate expected similarity (ES_ij)
  cell_type_proportions <- prop.table(table(grouping_vector))
  cell_type_proportions <- cell_type_proportions[colnames(similarity_matrix)]

  cell_type_names <- names(cell_type_proportions)
  expected_similarity <- sum(sapply(cell_type_names, function(i) {
    sum(similarity_matrix[cell_type_names, i] * cell_type_proportions[i] * cell_type_proportions)
  }))

  # Identify high-variance bins (HVB) based on the specified quantile
  bin_variance <- sparseMatrixStats::rowVars(x, useNames = TRUE)
  HVB <- names(bin_variance[bin_variance >= quantile(bin_variance, HVB_q)])
  x_HVB <- x[HVB, ]

  # Compute cosine distance matrix for high-variance bins
  dmat <- as.matrix(1 - proxyC::simil(x_HVB, method = "correlation", margin = 2))

  # Set up clustering vector based on cell types
  cluster_vector <- setNames(as.character(grouping_vector), colnames(x))[colnames(dmat)]

  # Initialize and build Annoy index for Euclidean distance with 1000 trees
  a <- new(RcppAnnoy::AnnoyEuclidean, ncol(dmat))
  a$setSeed(42)
  for (i in seq_len(nrow(dmat))) {
    a$addItem(i - 1, dmat[i, ])  # Annoy uses zero indexing
  }
  a$build(1000)

  # Calculate observed similarity by finding K nearest neighbors and computing mean similarity
  observed_similarity <- mean(unlist(parallel::mclapply(1:a$getNItems(), function(i) {
    nn <- a$getNNsByItem(i - 1, K + 1) + 1  # Adjust for 0-based index
    ct_e <- cluster_vector[nn]
    mean(similarity_matrix[ct_e[1], ct_e[-1]])
  }, mc.cores = N_cores)))

  # Return observed, expected, and delta similarity values
  output <- data.frame(
    S_ij = observed_similarity,
    ES_ij = expected_similarity,
    delta = observed_similarity / expected_similarity
  )

  return(output)
}





