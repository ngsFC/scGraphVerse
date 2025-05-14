library(scGraphVerse)

# 1. Download PBMC data
url <- paste0(
    "https://www.dropbox.com/s/r8qwsng79rhp9gf/",
    "SCA_scRNASEQ_TISSUE_WHOLE_BLOOD.RDS?dl=1"
)
seu <- download_Atlas(file_url = url)

# 2. Select top 500 T-cell genes
genes <- selgene(
    object = seu,
    top_n = 100,
    cell_type = "T_cells",
    cell_type_col = "CELL_TYPE",
    remove_rib = TRUE,
    remove_mt = TRUE
)

# 3. Retrieve STRINGdb adjacency
str_res <- stringdb_adjacency(
    genes = genes,
    species = 9606,
    required_score = 900,
    keep_all_genes = FALSE
)

wadj_truth <- str_res$weighted
adj_truth <- str_res$binary

# 4. Symmetrize and sort
common <- intersect(rownames(adj_truth), colnames(adj_truth))
adj_truth <- adj_truth[common, common]
adj_truth <- adj_truth[order(rownames(adj_truth)), order(colnames(adj_truth))]

# 2. Simulating Zero-Inflated Count Data
# Simulation parameters
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
# Transpose to cells × genes
count_matrices <- lapply(sims, t)

usethis::use_data(count_matrices, overwrite = TRUE)
usethis::use_data(adj_truth, overwrite = TRUE)
