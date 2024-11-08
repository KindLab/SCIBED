test_that("correlation works", {

  data(zeller_H3K4me3_matrix)
  data(zeller_H3K4me3_metadata)

  ground_truth_mat <- generate_ground_truth(
    x = zeller_H3K4me3_matrix,
    grouping_vector = zeller_H3K4me3_metadata,
    bin_size = 50000,
    fzc_sd = 0)

  correlation_stats <- calculate_correlation(zeller_H3K4me3_matrix,
                                             grouping_vector = zeller_H3K4me3_metadata,
                                             bin_size = 50000,
                                             ground_truth = ground_truth_mat,
                                             method = "pearson")

  expect_equal(round(correlation_stats$mu[correlation_stats$group == 'Eryths'],2), 0.49)
})


test_that("enrichment works", {

  data(zeller_H3K4me3_matrix)
  data(zeller_H3K4me3_metadata)
  data(zeller_H3K4me3_peaks)

  enrichment_stats <- calculate_enrichment(
    x = zeller_H3K4me3_matrix,
    grouping_vector = zeller_H3K4me3_metadata,
    regions = zeller_H3K4me3_peaks)

  expect_equal(round(sum(enrichment_stats$SIP),4), 2.6814)
})



