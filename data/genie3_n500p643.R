# ---- Load library ----

setwd("/home/francescoc/Desktop/GRN_project/data")

library(tidyverse)
library(GENIE3)
library(igraph)
library(ggraph)
library(pROC)
library(RColorBrewer)
library(scales)
library(gridExtra)
library(INetTool)

source("./../R/earlyj.R")
source("./../R/infer_networks.R")
source("./../R/generate_adjacency.R")
source("./../R/symmetrize.R")
source("./../R/plotROC.R")
source("./../R/cutoff_adjacency.R")
source("./../R/pscores.R")
source("./../R/plotg.R")
source("./../R/create_consensus.R")
source("./../R/compare_consensus.R")

time <- list()

# ---- Load Count matrices and adjm ----

adjm <- as.matrix(read.table("./../analysis/adjm_n500p643.txt"))
count_matrices <- readRDS("./../analysis/count_matrices_n500xp643.RDS")
dim(count_matrices[[1]])

# ---- GENIE3 late integration ----

set.seed(1234)
time[["GENIE3_late_15Cores"]] <- system.time(
  genie3_late <- infer_networks(count_matrices, 
                                method="GENIE3",
                                nCores = 15)
)

saveRDS(genie3_late, "./../analysis/genie3_late_n500p643.RDS")

# ---- Symmetrize and ROC ----

genie3_late_wadj <- generate_adjacency(genie3_late, ground.truth = adjm)
sgenie3_late_wadj <- symmetrize(genie3_late_wadj, weight_function = "mean")

png("./../analysis/plots/genie3_late_n500p643_auc.png", width = 2400, height = 1800, res = 300)
genie3_late_auc <- plotROC(sgenie3_late_wadj, adjm, plot_title = "ROC curve - GENIE3 Late Integration")
dev.off()

# ---- Cutoff ----

sgenie3_late_adj <- cutoff_adjacency(count_matrices = count_matrices,
                                     weighted_adjm_list = sgenie3_late_wadj, 
                                     ground.truth = adjm,
                                     n = 3,
                                     method = "GENIE3",
                                     nCores = 15)

png("./../analysis/plots/genie3_late_n500p643_scores.png", width = 2400, height = 1800, res = 300)
scores.genie3.late.all <- pscores(adjm, sgenie3_late_adj)
dev.off()

png("./../analysis/plots/genie3_late_n500p643_mplots.png", width = 2400, height = 2400, res = 300)
plots <- plotg(sgenie3_late_adj)
dev.off()

# ---- Consensus ----

consesusm <- create_consensus(sgenie3_late_adj, method="vote")
consesusu <- create_consensus(sgenie3_late_adj, method="union")
consesunet <- create_consensus(adj_matrix_list = sgenie3_late_adj, weighted_list = sgenie3_late_wadj, method = "INet", threshold = 0.05, ncores = 15)

png("./../analysis/plots/genie3_late_n500p643_vote_score.png", width = 2400, height = 2400, res = 300)
scores.genie3.late <- pscores(adjm, list(consesusm))
dev.off()
png("./../analysis/plots/genie3_late_n500p643_union_score.png", width = 2400, height = 2400, res = 300)
scoresu.genie3.late <- pscores(adjm, list(consesusu))
dev.off()
png("./../analysis/plots/genie3_late_n500p643_inet_score.png", width = 2400, height = 2400, res = 300)
scoresnet.genie3.late <- pscores(adjm, list(consesunet))
dev.off()

# ---- Plot comparison ----

png("./../analysis/plots/genie3_late_n500p643_vote_plot.png", width = 4000, height = 2400, res = 300)
ajm_compared <- compare_consensus(consesusm, adjm)
dev.off()
png("./../analysis/plots/genie3_late_n500p643_union_plot.png", width = 4000, height = 2400, res = 300)
ajm_compared <- compare_consensus(consesusu, adjm)
dev.off()
png("./../analysis/plots/genie3_late_n500p643_inet_plot.png", width = 4000, height = 2400, res = 300)
ajm_compared <- compare_consensus(consesunet, adjm)
dev.off()

# ---- GENIE3 early integration ----

early_matrix <- list(earlyj(count_matrices))

set.seed(1234)
time[["GENIE3_early_15Cores"]] <- system.time(
  genie3_early <- infer_networks(early_matrix, method="GENIE3", nCores = 15)
)

saveRDS(genie3_early, "./../analysis/genie3_early_n500p643.RDS")

# ---- Symmetrize and ROC ----

genie3_early_wadj <- generate_adjacency(genie3_early, ground.truth = adjm)
sgenie3_early_wadj <- symmetrize(genie3_early_wadj, weight_function = "mean")
png("./../analysis/plots/genie3_early_n500p643_auc.png", width = 2400, height = 1800, res = 300)
genie3_early_auc <- plotROC(sgenie3_early_wadj, adjm, plot_title = "ROC curve - GENIE3 Early Integration")
dev.off()

# ---- Cutoff ----

sgenie3_early_adj <- cutoff_adjacency(count_matrices = early_matrix,
                                      weighted_adjm_list = sgenie3_early_wadj, 
                                      ground.truth = adjm,
                                      n = 2,
                                      method = "GENIE3",
                                      nCores = 15)

png("./../analysis/plots/genie3_early_n500p643_scores.png", width = 2400, height = 1800, res = 300)
scores.genie3.early <- pscores(adjm, sgenie3_early_adj)
dev.off()

png("./../analysis/plots/genie3_early_n500p643_mplots.png", width = 2400, height = 2400, res = 300)
plots <- plotg(sgenie3_early_adj)
dev.off()

# ---- Plot comparison ----

png("./../analysis/plots/genie3_early_n500p643_plot.png", width = 4000, height = 2400, res = 300)
ajm_compared <- compare_consensus(sgenie3_early_adj[[1]], adjm)
dev.off()

