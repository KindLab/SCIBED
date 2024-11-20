

# matrix

dat <- readRDS("~/Projects/imputation_benchmark/sortChIC/fulll_matrices/sortChIC_H3K4me3_50kb_groundTruth_renamed.rds")

## get 100 cells of B-cells and Eryths
set.seed(42); idx_b <- sample(which(dat$anno == "Bcells"), size = 100)
set.seed(42); idx_e <- sample(which(dat$anno == "Eryths"), size = 100)
idx <- unname(sort(c(idx_e,idx_b)))

## subsample
zeller_H3K4me3_matrix <- dat$mat[,idx]
zeller_H3K4me3_metadata <- dat$anno[idx]
names(zeller_H3K4me3_metadata) <- colnames(zeller_H3K4me3_matrix)

## save
save(zeller_H3K4me3_matrix, file = "data/zeller_H3K4me3_matrix.rda")
save(zeller_H3K4me3_metadata, file = "data/zeller_H3K4me3_metadata.rda")

# peaks
zeller_H3K4me3_peaks <- GRangesList("Bcells" = import.bed("~/Downloads/H3K4me3_SortChIC_filtered_peaks/Bcells_peaks.bed"),
                                    "Eryths" = import.bed("~/Downloads/H3K4me3_SortChIC_filtered_peaks/Eryths_peaks.bed"))
save(zeller_H3K4me3_peaks, file = "data/zeller_H3K4me3_peaks.rda")

