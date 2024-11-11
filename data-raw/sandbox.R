
data("zeller_H3K4me3_matrix")
data("zeller_H3K4me3_metadata")
library(Matrix)
x = zeller_H3K4me3_matrix
grouping_vector = zeller_H3K4me3_metadata
bin_size = 5e4
ground_truth <- generate_ground_truth(x = zeller_H3K4me3_matrix,
                                      grouping_vector = zeller_H3K4me3_metadata, bin_size = bin_size, fzc_sd = 0)



data("zeller_H3K4me3_peaks")
regions <- zeller_H3K4me3_peaks


### Functions for package
### 2024.10.23
#
# ########################### IN SILICO DATASET GENERATION ###########################
#
# # Calculate fraction of zero-bins
# calculate_FZC <- function(data) {
#   n_zero <- sum(data == 0)
#   FZC = n_zero/length(data)
#   return(FZC)
# }
#
# # Filter by fraction of zero-bins
# filter_by_ZFC <- function(data) {
#   FZC_data <- apply(data, MARGIN = 2, calculate_FZC)
#   min_FZC <- mean(FZC_data) - 2*sd(FZC_data)
#   filtered_data <- data[, FZC_data >= min_FZC]
#   return(filtered_data)
# }
#
# # Add noise
# increase_noise <- function(data, increase_factor) {
#   orig_noise <- sum(data) / (nrow(data)*ncol(data))
#   output <- data + round(orig_noise) * increase_factor
#   return(output)
# }
#
# # Calculate pseudobulks
# calculate_pseudobulk <- function(data, metadata) {
#   metadataSubset <- metadata %>% filter(cell %in% colnames(data))
#   pseudobulks <- matrix(nrow = nrow(data), ncol = 0)
#   for (ct in unique(metadataSubset$ctype.from.LL)) {
#     cellIDs <- metadataSubset %>% filter(ctype.from.LL == ct) %>% select(cell)
#     tmp <- data[,colnames(data) %in% cellIDs$cell] %>% rowSums() %>% as.matrix()
#     colnames(tmp) <- ct
#     pseudobulks <- cbind(pseudobulks, tmp)
#   }
#   return(pseudobulks)
# }
#
# # Generate Ground Truth datasets
# generateGroundTruth <- function(data, metadata, genome, batch="all", min_RD, FZC=FALSE, binsize) {
#     genome=BSgenome::getBSgenome(genome_name)
#     # Filter dataset by batch (if applicable)
#     if (batch != "all") {
#         message(paste0("\tSelected batch: ", batch))
#         data_newCells <- metadata$cell[metadata$batch == batch]
#         data <- data[,colnames(data) %in% data_newCells] %>% as.data.frame()
#     } else {
#         message("\tNo batch selected")
#     }
#     # By read depth (RD)
#     message(paste0("\tMinimum read depth: ", min_RD))
#     data <- data[,colSums(data) > as.numeric(min_RD)]
#     #By fraction of bins with 0 cuts (FZC)
#     if (FZC == TRUE) {
#         message("\tFiltering by ZFC")
#         data <- filter_by_ZFC(data)
#     }
#     message("\tExpanding matrix...")
#     if (genome_name=="BSgenome.Mmusculus.UCSC.mm10") {
#         chroms_to_use <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY")
#     } else {
#        stop("Are you sure this is the right genome?")
#     }
#     this_genome <- GRanges(names(seqlengths(genome)), IRanges(1,seqlengths(genome)), seqinfo = seqinfo(genome))
#     this_genome <- keepStandardChromosomes(this_genome, pruning.mode = 'coarse')
#     this_genome <- keepSeqlevels(this_genome,chroms_to_use,pruning.mode = 'coarse')
#
#     this_genome <- slidingWindows(this_genome, width = as.numeric(binsize), step=as.numeric(binsize)) %>%
#       as.data.frame() %>%
#       select(seqnames,start,end) %>%
#       mutate(start=str_trim(format(start-1, scientific = FALSE)),
#              end=str_trim(format(end, scientific=FALSE)),
#              region = paste0(seqnames,":",start,"-",end)) %>%
#       select(region) %>%
#       column_to_rownames("region")
#
#     data <- merge(data, this_genome, by="row.names", all = TRUE) %>% column_to_rownames("Row.names")
#     data[is.na(data)] <- 0
#     message(paste0("\tNumber bins final filtered matrix: ", dim(data)[1]))
#     message(paste0("\tNumber cells final filtered matrix: ", dim(data)[2]))
#
#     #Generate ground truth (By calculating pseudobulks)
#     message("Calculating ground truth...")
#     metadata <- metadata[metadata$cell %in% colnames(data),]
#     groundTruth_pseudobulks <- calculate_pseudobulk(data, metadata=metadata)
#
#     return(groundTruth_pseudobulks)
# }
#
# # Generate In silico datasets
# generateInSilico <- function(data, metadata, genome, batch="all", min_RD, FZC=FALSE, RD_is, noise, binsize) {
#     #Generate ground truth
#     groundTruth_pseudobulks <- generateGroundTruth(data, metadata, genome, batch, min_RD, FZC, binsize)
#     # Create ground truth matrix
#     groundTruth_matrix <- sapply(unique(metadata$ctype.from.LL), function(ct) {
#         n_ct <- metadata %>% group_by(ctype.from.LL) %>% tally() %>% as.data.frame()
#         output <- groundTruth_pseudobulks[,rep(grep(colnames(groundTruth_pseudobulks), pattern=paste0("^",ct)), each=n_ct$n[n_ct$ctype.from.LL==ct])]
#         colnames(output) <- paste0(colnames(output), "_", seq.int(ncol(output)))
#         return(output)
#     }) %>%
#         do.call(what="cbind")
#
#     ### Add noise
#     message("Adding noise...")
#     data_noised <- increase_noise(increase_factor=as.numeric(noise), data=data)
#
#     ### Downsample (to introduce sparcity)
#     message("Downsampling...")
#     data_is <- Seurat::SampleUMI(data_noised, max.umi=as.numeric(RD_is))
#
#     return(list("GroundTruth"=groundTruth_pseudobulks, "InSilico"=data_is))
# }






########################### METRICS ###########################

# # Pearson correlations per cell
# calculatePearsonCorr <- function(data, groundTruth, genome, binsize, name) {
#     genome <- BSgenome::getBSgenome(genome_name)
#
#     ### Expanding matrix
#     message("Expanding matrix...")
#     if (genome_name=="BSgenome.Mmusculus.UCSC.mm10") {
#         chroms_to_use <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY")
#     } else {
#        stop("Are you sure this is the right genome?")
#     }
#     this_genome <- GRanges(names(seqlengths(genome)), IRanges(1,seqlengths(genome)), seqinfo = seqinfo(genome))
#     this_genome <- keepStandardChromosomes(this_genome, pruning.mode = 'coarse')
#     this_genome <- keepSeqlevels(this_genome,chroms_to_use,pruning.mode = 'coarse')
#
#     this_genome <- slidingWindows(this_genome, width = as.numeric(binsize), step=as.numeric(binsize)) %>%
#         as.data.frame() %>%
#         select(seqnames,start,end) %>%
#         mutate(start=str_trim(format(start-1, scientific = FALSE)),
#             end=str_trim(format(end, scientific=FALSE)),
#             region = paste0(seqnames,":",start,"-",end)) %>%
#         select(region) %>%
#         column_to_rownames("region")
#
#     data_expand <- merge(data, this_genome, by="row.names", all.y=T) %>% column_to_rownames("Row.names")
#     groundTruth <- merge(groundTruth, this_genome, by="row.names") %>% column_to_rownames("Row.names")
#     data_expand[is.na(data_expand)] <- 0
#
#     ### Calculate correlations
#     message("Calculating correlations to ground truth...")
#     message("\tCalculating correlations...")
#     metadata <- data.frame('cell'=colnames(data_expand)) %>% mutate(tmp=cell) %>% separate(tmp, into=c("ctype.from.LL", "number"))
#     correlationsToGT <- lapply(metadata$cell, function(cell) {
#         ctype <- stringr::str_split_1(cell, "_")[1]
#         corr <- data.frame("cell"=cell,"corr"=cor(data_expand[,cell], groundTruth[,ctype], method = "pearson"))
#     }) %>% dplyr::bind_rows() %>%
#         mutate('dataset'=name)
#
#     return(correlationsToGT)
# }

# Similarity score
calculateSIMIC <- function(data, cos_dist, mark, name) {
    cell_metadata <- data.frame(sample=colnames(data)) %>%
    separate(sample, into = c("cell_type","cellNr"), remove = F)
    ctype_eq <- read.csv("/hpc/hub_kind/rvanderweide/projects/imputation/sortChIC/full_dataset_20240814/naming_conventions.csv", sep = ";")
    cell_metadata <- merge(cell_metadata, ctype_eq, by="cell_type") %>% column_to_rownames('sample')
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts=data), colData = cell_metadata[match(colnames(data), rownames(cell_metadata)),])

    similarity_matrix <- TSPred::minmax(1 - cos_dist)
    # Compute proportions of clusters
    cell_type_proportions <- prop.table(table(as.character(cell_metadata[[mark]])))
    cell_type_proportions <- cell_type_proportions[colnames(cos_dist)]
    # Get ES_ij (expected similarity)
    cell_type_names <- names(cell_type_proportions)
    expected_similarity <- sum(sapply(cell_type_names, function(i) {
        sum(similarity_matrix[cell_type_names, i] * cell_type_proportions[i] * cell_type_proportions)
    }))

    message("Calculate dataset similarity")
    # Compute similarity for observed dataset
    # Calculate variance for each row (bin) in the assay matrix
    bin_variance <- sparseMatrixStats::rowVars(assay(sce), useNames=T)
    # Select top 95% high variance bins
    HVB <- names(which(bin_variance >= quantile(bin_variance, .95)))
    # Subset the SCE object with high variance bins
    sce_sub <- sce[HVB, ]
    rm(sce); gc(verbose = F)
    # Calculate correlation-based similarity matrix
    dmat <- proxyC::simil(assay(sce_sub), method = "correlation", margin = 2)
    # Convert similarity matrix to distance matrix
    dmat <- as.matrix(1 - dmat)
    # Create a cluster vector for cell types based on column names
    cluster_vector <- setNames(sce_sub[[mark]], colnames(sce_sub))[rownames(dmat)]
    cluster_vector <- as.character(cluster_vector)
    # Initialize Annoy index for Euclidean distances
    a <- new(RcppAnnoy::AnnoyEuclidean, ncol(dmat))
    a$setSeed(42)
    a$setVerbose(0)
    # Add each row in the distance matrix as an item to the Annoy index
    for (i in 1:nrow(dmat)) {
        a$addItem(i - 1, dmat[i, ]) # Annoy uses zero indexing
    }
    # Build the Annoy index with 1000 trees
    a$build(1000)
    rm(dmat); gc(verbose = F)
    # Compute observed similarity by finding nearest neighbors and calculating mean similarity
    observed_similarity <- mean(unlist(parallel::mclapply(1:a$getNItems(), function(i) {
          # Get 10 nearest neighbors for the current item
          nn10 <- a$getNNsByItem(i - 1, 11) + 1 # Annoy uses zero indexing
        # Get the cell types of the nearest neighbors
        ct_e <- cluster_vector[nn10]
        # Calculate mean similarity for the first cell type compared to the others
        mean(similarity_matrix[ct_e[1], ct_e[-1]])
    }, mc.cores = 7)))
    # Return the results as a data frame
    output <- data.frame(
        dataset = name,
        S_ij = observed_similarity,
        ES_ij = expected_similarity,
        delta = observed_similarity / expected_similarity)

    return(output)
}


#
# # FRiP score
# calculateFRiP <- function(data, ct_lin, group_var, regions, binsize, name) {
#     metadata <- data.frame('cell'=colnames(data)) %>% mutate(tmp=cell) %>% separate(tmp, into=c("ctype.from.LL", "number"))
#     metadata <- merge(metadata, ct_lin, by="ctype.from.LL")
#     metadata <- metadata[match(colnames(data), metadata$cell),]
#     chroms_to_use <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY")
#
#     windows <- tileGenome(seqinfo(BSgenome.Mmusculus.UCSC.mm10), tilewidth = as.numeric(binsize), cut.last.tile.in.chrom = TRUE)
#     regions_rebinned <- lapply(regions, function(x) {
#         ctype <- unique(x$ctype.from.LL)
#         output <- subsetByOverlaps(windows,x)
#         output$ctype.from.LL <- ctype
#         return(output)
#     })
#
#     metrics <- lapply(na.omit(unique(regions$ctype.from.LL)), function(g){
#         # subset peaks
#         peaks_all <- regions[regions$ctype.from.LL == g]
#         peaks_all_rebin <- regions_rebinned[regions_rebinned$ctype.from.LL == g]
#         # subset count table to relevant cells
#         cells <- na.omit(metadata[metadata[,variable] == g,])$cell
#         subset_data <- data[,colnames(data) %in% cells] %>% as.data.frame() %>% rownames_to_column("regions") %>% separate(regions, into=c("chrom","coordinates"), sep=":") %>%
#              separate(coordinates, into=c("start","end"), sep="-") %>% makeGRangesFromDataFrame(keep.extra.columns=T)
#         subset_data <- keepSeqlevels(subset_data,chroms_to_use,pruning.mode = 'coarse')
#         subset_data_chrPeaks <- subset_data %>% filter(seqnames %in% unique(seqnames(peaks_all)))
#         # Calculate observed FRiP
#         message(paste0("\t",as.character(g)," - Creating FRiP-related metrics..."))
#         data_obsPeaks = subsetByOverlaps(subset_data, peaks_all)
#         data_obsPeaks_chrPeaks = subsetByOverlaps(subset_data_chrPeaks, peaks_all)
#         FRiP_obs_all = colSums(as.data.frame(data_obsPeaks@elementMetadata))/colSums(as.data.frame(subset_data@elementMetadata))
#         obs_FRiP_chrPeaks = colSums(as.data.frame(data_obsPeaks_chrPeaks@elementMetadata))/colSums(as.data.frame(subset_data_chrPeaks@elementMetadata))
#         # return dataframe with relevant information
#         FRiP_stats = data.frame("cell"=names(FRiP_obs_all),
#             "obs_FRiP_all"=FRiP_obs_all,
#             "obs_FRiP_chrPeaks"=obs_FRiP_chrPeaks) %>%
#         mutate(ctype.from.LL = g,
#             width_all_peaks = sum(width(reduce(peaks_all_rebin))),
#             perc_all_peaks = sum(width(reduce(peaks_all_rebin)))/sum(width(reduce(windows))))
#         rownames(FRiP_stats) <- NULL
#         return(FRiP_stats)
#     }) %>% do.call(rbind, .) %>%
#         mutate('dataset'=gsub(name, pattern=".rds", replacement=""))
#
#     return(metrics)
# }
#





########################### RANKING ###########################

rankImpMethods <- function(correlations, similarity, FRiP_peaks, fail="None") {
    read_depths <- c("100","1000","10000")
    noise_levels <- c("0","0.25","0.35","0.5","0.75","1","1.25","2","3.5","5","10","50")
    algorithm_datasets <- c("noImputation","MAGIC","knn3","knn7","knn15","SAVER","scImpute.Lineage","scImpute.Group","scImpute.CellType","scOpen","ALRA","cisTopic","DCA","SCALE","SCALEX.RNA","SCALEX.ATAC")
    imputed_dataset <- c("MAGIC","knn3","knn7","knn15","SAVER","scImpute.Lineage","scImpute.Group","scImpute.CellType","scOpen","ALRA","cisTopic","DCA","SCALE","SCALEX.RNA","SCALEX.ATAC")

    approachClass_levels <- c("Autoencoder", "Dictionary learning", "Matrix theory","Data smoothing","Statistical Model","noImputation")
    designedOn_levels <- c("scRNA","scATAC","scATAC+scRNA","noImputation")
    alg_classification <- read.csv("./raw_data/algorithmClassification_v2.csv", sep=";") %>%
        mutate(designedOn=factor(designedOn, levels=designedOn_levels),
            approachClass=factor(approachClass, levels=approachClass_levels))

    # Correlations
    correlations <- correlations %>%
        separate(col=dataset, into=c("dataType","tech","tissue", "mark","version","noise","RD","algorithm","countType"), sep="_") %>%
        mutate(algorithm = case_when(algorithm=="scImpute.k4"~"scImpute.Lineage",
                                    algorithm=="scImpute.k8"~"scImpute.Group",
                                    algorithm=="scImpute.k14"~"scImpute.CellType",
                                    TRUE~algorithm)) %>%
        filter(!(algorithm=="noImputation" & countType=="RAW") & !(algorithm=="scImpute.Lineage" & countType=="CPM") &
            !(algorithm=="scImpute.Group" & countType=="CPM") & !(algorithm=="scImpute.CellType" & countType=="CPM")) %>%
        mutate( noise=factor(gsub(noise, pattern="noise",replacement=""), levels=noise_levels),
            RD=factor(gsub(RD, pattern="RD", replacement = ""), levels=read_depths),
            algorithm=factor(algorithm, levels = algorithm_datasets),
            RD_levels=case_when(RD=="100"~"Low RD", RD=="1000"~"Mid RD", RD=="10000"~"High RD"),
            noise_levels=case_when(noise == "0"~"None", noise %in% c("0.25","0.35","0.5") ~"Low", noise %in% c("0.75","1","1.25")~"Mid", noise %in% c("2","3.5","5")~"High", noise %in% c("10","50")~"Extremely high")) %>%
        filter(RD %in% c("10000","1000","100") & noise %in% c("0","0.25","0.35","0.5","0.75","1","1.25","2","3.5","5","10","50")) %>%
        mutate(noise_levels = factor(noise_levels, levels=c("None","Low", "Mid","High", "Extremely high")),
            signalNoiseRatio = paste0(as.character(noise),"X"),
            signalNoiseRatio_levels=case_when(signalNoiseRatio == "0X"~"No noise",
                signalNoiseRatio %in% c("0.25X","0.35X","0.5X") ~"High",
                signalNoiseRatio %in% c("0.75X","1X","1.25X")~"Mid",
                signalNoiseRatio %in% c("2X","3.5X","5X")~"Low",
                signalNoiseRatio %in% c("10X","50X")~"Extremely Low"),
            RD_levels = factor(RD_levels, levels=c("High RD","Mid RD","Low RD"))) %>%
        mutate(signalNoiseRatio_levels = factor(signalNoiseRatio_levels, levels=c("No noise", "High", "Mid", "Low", "Extremely Low"))) %>%
        select(!countType)

    correlations_summarized <- correlations_data %>% group_by(noise_levels, RD_levels,algorithm,signalNoiseRatio_levels) %>%
        dplyr::summarise(medianCorr=median(corr))

    correlations_summarized2 <- correlations_data %>% group_by(noise, noise_levels, RD_levels,algorithm,signalNoiseRatio_levels, signalNoiseRatio) %>%
        dplyr::summarise(medianCorr=median(corr)) %>%
        mutate(x_labels = paste0(noise, "&", noise_levels),
        x_labels_SNratio = paste0(signalNoiseRatio, "&", signalNoiseRatio_levels))
    correlations_summarized2 <- correlations_summarized2 %>%
        mutate(x_labels = factor(x_labels, levels = unique(correlations_summarized2$x_labels)),
            x_labels_SNratio = factor(x_labels_SNratio, levels = unique(correlations_summarized2$x_labels_SNratio)))

    rm(correlations_data)

    # Similarity
    similarity <- similarity %>%
            as.data.frame()  %>%
        separate(col=dataset, into=c("dataType","tech","tissue", "mark","version","noise","RD","algorithm","countType"), sep="_") %>%
        mutate(algorithm = case_when(algorithm=="scImpute.k4"~"scImpute.Lineage",
                                algorithm=="scImpute.k8"~"scImpute.Group",
                                algorithm=="scImpute.k14"~"scImpute.CellType",
                                TRUE~algorithm)) %>%
        filter(!(algorithm=="noImputation" & countType=="RAW") & !(algorithm=="scImpute.Lineage" & countType=="CPM") &
            !(algorithm=="scImpute.Group" & countType=="CPM") & !(algorithm=="scImpute.CellType" & countType=="CPM")) %>%
        mutate( noise=factor(gsub(noise, pattern="noise",replacement=""), levels=noise_levels),
            RD=factor(gsub(RD, pattern="RD", replacement = ""), levels=read_depths),
            algorithm=factor(algorithm, levels = algorithm_datasets),
            RD_levels=case_when(RD=="100"~"Low RD", RD=="1000"~"Mid RD", RD=="10000"~"High RD"),
            noise_levels=case_when(noise == "0"~"None", noise %in% c("0.25","0.35","0.5") ~"Low", noise %in% c("0.75","1","1.25")~"Mid", noise %in% c("2","3.5","5")~"High", noise %in% c("10","50")~"Extremely high")) %>%
        filter(RD %in% c("10000","1000","100") & noise %in% c("0","0.25","0.35","0.5","0.75","1","1.25","2","3.5","5","10","50")) %>%
        mutate(noise_levels = factor(noise_levels, levels=c("None","Low", "Mid","High", "Extremely high")),
            signalNoiseRatio = paste0(as.character(noise),"X"),
            signalNoiseRatio_levels=case_when(signalNoiseRatio == "0X"~"No noise",
                signalNoiseRatio %in% c("0.25X","0.35X","0.5X") ~"High",
                signalNoiseRatio %in% c("0.75X","1X","1.25X")~"Mid",
                signalNoiseRatio %in% c("2X","3.5X","5X")~"Low",
                signalNoiseRatio %in% c("10X","50X")~"Extremely Low"),
            RD_levels = factor(RD_levels, levels=c("High RD","Mid RD","Low RD"))) %>%
        mutate(signalNoiseRatio_levels = factor(signalNoiseRatio_levels, levels=c("No noise", "High", "Mid", "Low", "Extremely Low"))) %>%
        select(!countType)

    similarity_data_summarized <- similarity_data %>%
        group_by(noise, noise_levels, RD_levels,algorithm,signalNoiseRatio_levels, signalNoiseRatio) %>%
        mutate(x_labels = paste0(noise, "&", noise_levels),
        x_labels_SNratio = paste0(signalNoiseRatio, "&", signalNoiseRatio_levels))
    similarity_data_summarized <- similarity_data_summarized[order(similarity_data_summarized$noise),]
    similarity_data_summarized <- similarity_data_summarized %>%
        mutate(x_labels = factor(x_labels, levels = unique(similarity_data_summarized$x_labels)),
        x_labels_SNratio = factor(x_labels_SNratio, levels = unique(similarity_data_summarized$x_labels_SNratio))) %>%
        select(!ES_ij & !S_ij)

    # FRiP
    FRiP_peaks <- FRiP_peaks %>%
        separate(col=dataset, into=c("dataType","tech","tissue", "mark","version","noise","RD","algorithm","countType"), sep="_") %>%
        mutate(algorithm = case_when(algorithm=="scImpute.k4"~"scImpute.Lineage",
                                    algorithm=="scImpute.k8"~"scImpute.Group",
                                    algorithm=="scImpute.k14"~"scImpute.CellType",
                                    TRUE~algorithm)) %>%
        filter(!(algorithm=="noImputation" & countType=="RAW") & !(algorithm=="scImpute.Lineage" & countType=="CPM") &
            !(algorithm=="scImpute.Group" & countType=="CPM") & !(algorithm=="scImpute.CellType" & countType=="CPM")) %>%
        mutate( noise=factor(gsub(noise, pattern="noise",replacement=""), levels=noise_levels),
                RD=factor(gsub(RD, pattern="RD", replacement = ""), levels=read_depths),
                algorithm=factor(algorithm, levels = algorithm_datasets),
                RD_levels=case_when(RD=="100"~"Low RD", RD=="1000"~"Mid RD", RD=="10000"~"High RD"),
                noise_levels=case_when(noise == "0"~"None", noise %in% c("0.25","0.35","0.5") ~"Low", noise %in% c("0.75","1","1.25")~"Mid", noise %in% c("2","3.5","5")~"High", noise %in% c("10","50")~"Extremely high")) %>%
        filter(RD %in% c("10000","1000","100") & noise %in% c("0","0.25","0.35","0.5","0.75","1","1.25","2","3.5","5","10","50")) %>%
        mutate(noise_levels = factor(noise_levels, levels=c("None","Low", "Mid","High","Extremely high")),
            signalNoiseRatio = paste0(as.character(noise),"X"),
            signalNoiseRatio_levels=case_when(signalNoiseRatio == "0X"~"No noise",
                signalNoiseRatio %in% c("0.25X","0.35X","0.5X") ~"High",
                signalNoiseRatio %in% c("0.75X","1X","1.25X")~"Mid",
                signalNoiseRatio %in% c("2X","3.5X","5X")~"Low",
                signalNoiseRatio %in% c("10X","50X")~"Extremely Low"),
            RD_levels = factor(RD_levels, levels=c("High RD","Mid RD","Low RD"))) %>%
        select(!countType)

    FRiP_peaks_summarized <- FRiP_peaks_data %>%
        group_by(noise_levels, RD_levels,algorithm,signalNoiseRatio_levels) %>%
        dplyr::summarise(medianObsFRiP_peaks=median(obs_FRiP_all))

    FRiP_peaks_summarized2 <- FRiP_peaks_data %>%
        group_by(noise, noise_levels, RD_levels,algorithm,signalNoiseRatio_levels, signalNoiseRatio) %>%
        dplyr::summarise(medianObsFRiP_peaks=median(obs_FRiP_all)) %>%
        mutate(x_labels = paste0(noise, "&", noise_levels),
            x_labels_SNratio = paste0(signalNoiseRatio, "&", signalNoiseRatio_levels))
    FRiP_peaks_summarized2 <- FRiP_peaks_summarized2 %>%
        mutate(x_labels = factor(x_labels, levels = unique(FRiP_peaks_summarized2$x_labels)),
        x_labels_SNratio = factor(x_labels_SNratio, levels = unique(FRiP_peaks_summarized2$x_labels_SNratio)))
    rm(FRiP_peaks_data)


    # merge all scores into one
    data <- expand.grid(algorithm=unique(correlations_summarized2$algorithm),
        RD_levels=unique(correlations_summarized2$RD_levels),
        noise=unique(correlations_summarized2$noise)) %>%
        mutate(algorithm=factor(algorithm, levels = algorithm_datasets),
           RD_level = factor(RD_levels, levels=c("High RD","Mid RD","Low RD")))%>%
        mutate(noise_levels=case_when(noise == "0"~"None", noise %in% c("0.25","0.35","0.5") ~"Low", noise %in% c("0.75","1","1.25")~"Mid", noise %in% c("2","3.5","5")~"High",noise %in% c("10","50")~"Extremely high"),
           noise_levels = factor(noise_levels, levels=c("None","Low", "Mid","High","Extremely high")),
            signalNoiseRatio = paste0(as.character(noise),"X"),
            signalNoiseRatio_levels=case_when(signalNoiseRatio == "0X"~"No noise",
                signalNoiseRatio %in% c("0.25X","0.35X","0.5X") ~"High",
                signalNoiseRatio %in% c("0.75X","1X","1.25X")~"Mid",
                signalNoiseRatio %in% c("2X","3.5X","5X")~"Low",
                signalNoiseRatio %in% c("10X","50X")~"Extremely Low"))

    data <- merge(data,correlations_summarized2, all.x=T)
    data <- merge(data,FRiP_peaks_summarized2, all.x=T)
    data <- merge(data,similarity_data_summarized, all.x=T)

    if (fail != "None") {
        runFails <- fail %>%
            mutate(file=gsub(file, pattern=".rds", replacement="")) %>%
            separate(col=file, into=c("dataType","tech","tissue", "mark","version","noise","RD","algorithm","countType"), sep="_") %>%
            mutate( noise=factor(gsub(noise, pattern="noise",replacement=""), levels=noise_levels),
                RD=factor(gsub(RD, pattern="RD", replacement = ""), levels=read_depths),
                algorithm=factor(algorithm, levels = algorithm_datasets),
                RD_levels=case_when(RD=="100"~"Low RD", RD=="1000"~"Mid RD", RD=="10000"~"High RD"),
                noise_levels=case_when(noise == "0"~"None", noise %in% c("0.25","0.35","0.5") ~"Low", noise %in% c("0.75","1","1.25")~"Mid", noise %in% c("2","3.5","5")~"High", noise %in% c("10","50")~"Extremely high")) %>%
            filter(RD %in% c("10000","1000","100") & noise %in% c("0","0.25","0.35","0.5","0.75","1","1.25","2","3.5","5","10","50")) %>%
            mutate(noise_levels = factor(noise_levels, levels=c("None","Low", "Mid","High","Extremely high")),
            signalNoiseRatio = paste0(as.character(noise),"X"),
            signalNoiseRatio_levels=case_when(signalNoiseRatio == "0X"~"No noise",
                signalNoiseRatio %in% c("0.25X","0.35X","0.5X") ~"High",
                signalNoiseRatio %in% c("0.75X","1X","1.25X")~"Mid",
                signalNoiseRatio %in% c("2X","3.5X","5X")~"Low",
                signalNoiseRatio %in% c("10X","50X")~"Extremely Low"),
                RD_levels = factor(RD_levels, levels=c("High RD","Mid RD","Low RD"))) %>%
            select(noise,algorithm,cause,RD_levels,noise_levels,signalNoiseRatio,signalNoiseRatio_levels)

    data <- merge(data, runFails, all.x=TRUE) %>%
    mutate(run=case_when(is.na(cause)==TRUE ~ TRUE,
        is.na(cause)==FALSE ~ FALSE,))
    }

    message("Fill NAs with the minimum scores & scale...")
    data2 <- data %>%
        filter(signalNoiseRatio_levels != "Extremely Low") %>%
        mutate(scenario = paste0(noise,"_",RD_levels), run = case_when(is.na(medianCorr)==T ~ F, T~T)) %>%
        select(!dataType & !tech &!tissue & !mark & !version & !RD)

    data3 <- lapply(unique(data2$scenario), function(s) {
        subset <- data2 %>% filter(scenario == s)

        subset[,10:(ncol(subset)-3)] <- apply(subset[,10:(ncol(subset)-3)],2,function(x) {
            y <- x
            y[is.na(y)] <- min(na.omit(y))
            if (min(y) == max(y)) {
                y <- 1
            } else {
                y <- (y-min(y))/(max(y)-min(y))
            }
            return(y)
        })
        return(subset)
    }) %>% bind_rows()

    data3 %>%
        mutate(
            Correlation=medianCorr,
            simScore = delta,
            FRiPscore=if("medianObsFRiP_genes" %in% colnames(.)) (medianObsFRiP_peaks + medianObsFRiP_genes)/2 else medianObsFRiP_peaks,
            score = (Correlation + simScore + FRiPscore)/3)


    message("Calculate ranking...")
    ranks <- data3 %>%
        mutate(
            Correlation=medianCorr,
            simScore = delta,
            FRiPscore=if("medianObsFRiP_genes" %in% colnames(.)) (medianObsFRiP_peaks + medianObsFRiP_genes)/2 else medianObsFRiP_peaks,
            score = (Correlation + simScore + FRiPscore)/3) %>%
        group_by(algorithm) %>%
        dplyr::reframe(
            similarity=mean(simScore),correlation=mean(medianCorr),signalEnrichment=mean(FRiPscore),Final_score = mean(score),
            FRiP_peaks=mean(medianObsFRiP_peaks), FRiP_genes=if("medianObsFRiP_genes" %in% colnames(.)) mean(medianObsFRiP_genes) else NA) %>%
        ungroup() %>%
        mutate(alg_rank=rank(-Final_score)) %>%
        arrange(alg_rank)

    data_per_scenario <- data3 %>%
        mutate(
            Correlation=medianCorr,
            simScore = delta,
            FRiPscore=if("medianObsFRiP_genes" %in% colnames(.)) (medianObsFRiP_peaks + medianObsFRiP_genes)/2 else medianObsFRiP_peaks,
            score = (Correlation + simScore + FRiPscore)/3) %>%
        group_by(noise_levels, RD_levels, algorithm) %>%
        dplyr::reframe(scenarioScore = mean(score), run=sum(run)) %>%
        group_by(noise_levels, RD_levels) %>%
        mutate(scenario_rank=rank(-na.omit(scenarioScore)))

    test <- data3 %>%
        mutate(cause = case_when(is.na(cause) == TRUE ~ "", TRUE ~ cause)) %>%
        group_by(noise_levels, RD_levels, algorithm) %>%
        dplyr::reframe(run=sum(run),
            cause=paste0(unique(cause), collapse=""))

    data_per_scenario <- merge(data_per_scenario, test, all.x=TRUE) %>%
        mutate(scenario_rank=case_when(run == 0 ~ NA,
            TRUE ~ scenario_rank)) %>%
        mutate(scenario_rank_text = case_when(is.na(scenario_rank) == FALSE ~ as.character(scenario_rank),
            is.na(scenario_rank) == TRUE ~ "")) %>%
        mutate(scenario_rank_text = case_when(
            cause=="Technical failure" ~ paste0(scenario_rank_text, "*"),
            cause=="Too many resources" ~ paste0(scenario_rank_text, "\u2020"),
            T ~ as.character(scenario_rank_text)
            ))

    data_per_scenario <- merge(data_per_scenario, ranks)
    data_per_scenario <- merge(data_per_scenario, alg_classification)

    return(list("scores" = data, "scoresScaled"=data3, "rankingPerScenario"=data_per_scenario, "rankingOverall"=ranks))
}
