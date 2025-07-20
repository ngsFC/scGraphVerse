library(scGraphVerse)
library(TENxPBMCData)
library(scater)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(SingleR)
library(celldex)

sce <- TENxPBMCData("pbmc3k")

sce <- logNormCounts(pbmc_obj)
symbols_tenx <- rowData(sce)$Symbol_TENx
valid <- !is.na(symbols_tenx) & symbols_tenx != ""
sce <- sce[valid, ]
rownames(sce) <- make.unique(symbols_tenx[valid])
logcounts(sce) <- as.matrix(logcounts(sce))
colnames(sce) <- paste0("cell_", seq_len(ncol(sce)))

ref <- celldex::HumanPrimaryCellAtlasData()
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
colData(sce)$predicted_celltype <- pred$labels

genes <- selgene(
    object = sce,
    top_n = 100,
    cell_type = "T_cells",
    cell_type_col = "predicted_celltype",
    remove_rib = TRUE,
    remove_mt = TRUE
)

str_res <- stringdb_adjacency(
    genes = genes,
    species = 9606,
    required_score = 900,
    keep_all_genes = FALSE
)

wadj_truth <- str_res$weighted
adj_truth <- str_res$binary

common <- intersect(rownames(adj_truth), colnames(adj_truth))
adj_truth <- adj_truth[common, common]
adj_truth <- adj_truth[order(rownames(adj_truth)), order(colnames(adj_truth))]

nodes <- nrow(adj_truth)
sims <- zinb_simdata(
    n = 40,
    p = nodes,
    B = adj_truth,
    mu_range = list(c(1, 4), c(1, 7), c(1, 10)),
    mu_noise = c(1, 3, 5),
    theta = c(1, 0.7, 0.5),
    pi = c(0.2, 0.2, 0.2),
    kmat = 3,
    depth_range = c(0.8 * nodes * 3, 1.2 * nodes * 3)
)

count_matrices <- lapply(sims, t)
count_matrices <- lapply(count_matrices, function(mat) {
    col_data <- DataFrame(CELL_TYPE = rep("T_cells", ncol(mat)))
    SingleCellExperiment(assays = list(counts = mat), colData = col_data)
})

usethis::use_data(count_matrices, overwrite = TRUE)
usethis::use_data(adj_truth, overwrite = TRUE)
