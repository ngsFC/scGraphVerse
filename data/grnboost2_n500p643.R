# ---- Load library ----

setwd("/home/francescoc/Desktop/GRN_project/data")

library(tidyverse)
library(igraph)
library(ggraph)
library(pROC)
library(RColorBrewer)
library(scales)
library(gridExtra)
library(INetTool)
library(reticulate)

use_python("/usr/bin/python3", required = TRUE)
arboreto <- import("arboreto.algo")
pandas <- import("pandas")
numpy <- import("numpy")

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

ddir <- "/home/francescoc/Desktop/GRN_project/analysis/GRNBoost2/data"
pdir <- "/home/francescoc/Desktop/GRN_project/analysis/GRNBoost2/plots"

# ---- Load Count matrices and adjm ----

adjm <- as.matrix(read.table("./../analysis/adjm_n500p643.txt"))
count_matrices <- readRDS("./../analysis/count_matrices_n500xp643.RDS")
dim(count_matrices[[1]])

# ---- GRNBoost2 late integration ----

set.seed(1234)
time[["grnboost_late"]] <- system.time(
  grnboost_late <- infer_networks(count_matrices, 
                                method="GRNBoost2",
                                nCores = 15)
)

saveRDS(grnboost_late, paste(ddir, "grnboost_late_n500p643.RDS", sep = "/"))

# ---- Symmetrize and ROC ----

grnboost_late_wadj <- generate_adjacency(grnboost_late, ground.truth = adjm)
sgrnboost_late_wadj <- symmetrize(grnboost_late_wadj, weight_function = "mean")

png(paste(pdir, "grnboost_late_n500p643_auc.png", sep = "/"), width = 2400, height = 1800, res = 300)
grnboost_late_auc <- plotROC(sgrnboost_late_wadj, adjm, plot_title = "ROC curve - grnboost Late Integration")
dev.off()

# ---- Cutoff ----

sgrnboost_late_adj <- cutoff_adjacency(count_matrices = count_matrices,
                                     weighted_adjm_list = sgrnboost_late_wadj, 
                                     ground.truth = adjm,
                                     n = 3,
                                     method = "GRNBoost2",
                                     nCores = 15)

png(paste(pdir, "grnboost_late_n500p643_scores.png", sep = "/"), width = 2400, height = 1800, res = 300)
scores.grnboost.late.all <- pscores(adjm, sgrnboost_late_adj)
dev.off()

png(paste(pdir, "grnboost_late_n500p643_mplots.png", sep = "/"), width = 2400, height = 2400, res = 300)
plots <- plotg(sgrnboost_late_adj)
dev.off()

cbind(scores.grnboost.late.all$Statistics, grnboost_late_auc) %>% write.table(., paste(ddir, "grnboost_late_n500p643_perf.txt", sep = "/"), sep = "\t", quote = F, col.names = T, row.names = F)

# ---- Consensus ----

consesusm <- create_consensus(sgrnboost_late_adj, method="vote")
consesusu <- create_consensus(sgrnboost_late_adj, method="union")
consesunet <- create_consensus(adj_matrix_list = sgrnboost_late_adj, weighted_list = sgrnboost_late_wadj, method = "INet", threshold = 0.05)

png(paste(pdir, "grnboost_late_n500p643_vote_score.png", sep = "/"), width = 2400, height = 2400, res = 300)
scores.grnboost.late <- pscores(adjm, list(consesusm))
dev.off()
png(paste(pdir, "grnboost_late_n500p643_union_score.png", sep = "/"), width = 2400, height = 2400, res = 300)
scoresu.grnboost.late <- pscores(adjm, list(consesusu))
dev.off()
png(paste(pdir, "grnboost_late_n500p643_inet_score.png", sep = "/"), width = 2400, height = 2400, res = 300)
scoresnet.grnboost.late <- pscores(adjm, list(consesunet))
dev.off()

# ---- Plot comparison ----

png(paste(pdir, "grnboost_late_n500p643_vote_plot.png", sep = "/"), width = 4000, height = 2400, res = 300)
ajm_compared <- compare_consensus(consesusm, adjm)
dev.off()
png(paste(pdir, "grnboost_late_n500p643_union_plot.png", sep = "/"), width = 4000, height = 2400, res = 300)
ajm_compared <- compare_consensus(consesusu, adjm)
dev.off()
png(paste(pdir, "grnboost_late_n500p643_inet_plot.png", sep = "/"), width = 4000, height = 2400, res = 300)
ajm_compared <- compare_consensus(consesunet, adjm)
dev.off()

# ---- GRNBoost2 early integration ----

early_matrix <- list(earlyj(count_matrices))

set.seed(1234)
time[["grnboost_early"]] <- system.time(
  grnboost_early <- infer_networks(early_matrix, method="GRNBoost2")
)

saveRDS(grnboost_early, paste(ddir, "grnboost_early_n500p643.RDS", sep = "/"))

# ---- Symmetrize and ROC ----

grnboost_early_wadj <- generate_adjacency(grnboost_early, ground.truth = adjm)
sgrnboost_early_wadj <- symmetrize(grnboost_early_wadj, weight_function = "mean")

png(paste(pdir, "grnboost_early_n500p643_auc.png", sep = "/"), width = 2400, height = 1800, res = 300)
grnboost_early_auc <- plotROC(sgrnboost_early_wadj, adjm, plot_title = "ROC curve - grnboost Early Integration")
dev.off()

# ---- Cutoff ----

sgrnboost_early_adj <- cutoff_adjacency(count_matrices = early_matrix,
                                      weighted_adjm_list = sgrnboost_early_wadj, 
                                      ground.truth = adjm,
                                      n = 2,
                                      method = "GRNBoost2",
                                      nCores = 15)

png(paste(pdir, "grnboost_early_n500p643_scores.png", sep = "/"), width = 2400, height = 1800, res = 300)
scores.grnboost.early <- pscores(adjm, sgrnboost_early_adj)
dev.off()

png(paste(pdir, "grnboost_early_n500p643_mplots.png", sep = "/"), width = 2400, height = 2400, res = 300)
plots <- plotg(sgrnboost_early_adj)
dev.off()

cbind(scores.grnboost.early$Statistics, grnboost_early_auc) %>% write.table(., paste(ddir, "grnboost_early_n500p643_perf.txt", sep = "/"), sep = "\t", quote = F, col.names = T, row.names = F)

# ---- Plot comparison ----

png(paste(pdir, "grnboost_early_n500p643_plot.png", sep = "/"), width = 4000, height = 2400, res = 300)
ajm_compared <- compare_consensus(sgrnboost_early_adj[[1]], adjm)
dev.off()

time_data <- data.frame(
  Method = names(time),
  Time_in_Seconds = sapply(time, function(x) {
    if ("elapsed" %in% names(x)) x["elapsed"] else NA
  })
)

time_data$Time_in_Minutes <- as.numeric(time_data$Time_in_Seconds) / 60
time_data$Time_in_Hours <- as.numeric(time_data$Time_in_Seconds) / 3600

time_data <- time_data[order(time_data$Time_in_Hours), ]
time_data$Method <- factor(time_data$Method, levels = time_data$Method)

write.table(time_data, paste(ddir, "grnboost_n500p643_timeres.txt", sep = "/"), sep = "\t", quote = F, col.names = T, row.names = F)

