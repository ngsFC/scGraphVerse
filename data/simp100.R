setwd("/home/francescoc/Desktop/scGraphVerse/data/")
ddir <- "/home/francescoc/Desktop/scGraphVerse/analysis/simulation/results/"
pdir <- "/home/francescoc/Desktop/scGraphVerse/analysis/simulation/plot/"

run_single_simulation <- function(run_id = 1, seed_base = 1234) {

library(tidyverse)
library(scGraphVerse)
modules <- init_py(python_path ="/usr/bin/python3", required = TRUE)

time <- list()

ddir <- "/home/francescoc/Desktop/scGraphVerse/analysis/simulation/results/"
pdir <- "/home/francescoc/Desktop/scGraphVerse/analysis/simulation/plot/"

# Adjacency and Count matrices

adjm <- as.matrix(read.table("./../analysis/simulation/adjacency/adjm_top200_p100.txt"))
colnames(adjm) <- rownames(adjm)

count_matrices <- readRDS("./../analysis/simulation/simdata/sim_n100p100.RDS")
count_matrices <- lapply(count_matrices, t)

# GENIE3
## Late integration

time[["GENIE3_late_15Cores"]] <- system.time(
  late <- infer_networks(count_matrices, 
                         method="GENIE3",
                         nCores = 15,
                         seed=seed_base+run_id)
)

### Symmetrize and ROC

late_wadj <- generate_adjacency(late)
slate_wadj <- symmetrize(late_wadj, weight_function = "mean")
late_auc <- plotROC(slate_wadj, adjm, plot_title = "ROC curve - GENIE3 Late Integration", is_binary = F)

### Cutoff

slate_adj <- cutoff_adjacency(count_matrices = count_matrices,
                              weighted_adjm_list = slate_wadj, 
                              n = 3,
                              method = "GENIE3",
                              nCores = 15)

scores.late.all <- pscores(adjm, slate_adj)
plotg(slate_adj)

### Consensus

consesusm <- create_consensus(slate_adj, method="vote")
consesusu <- create_consensus(slate_adj, method="union")
consesunet <- create_consensus(adj_matrix_list = slate_adj, weighted_list = slate_wadj, method = "INet", threshold = 0.05, ncores = 15)

scores.late <- pscores(adjm, list(consesusm,consesusu,consesunet))

keepscores <- scores.late$Statistics %>% mutate(Predicted_Matrix=c("vote", "union", "inet"), Method=c("GENIE3"), Ratio=nrow(adjm)/ncol(count_matrices[[1]]), p=nrow(adjm))

### Plot comparison

compare_consensus(consensus_matrix = consesusm, reference_matrix = adjm, false_plot = F)
compare_consensus(consensus_matrix = consesusu, reference_matrix = adjm, false_plot = F)
compare_consensus(consensus_matrix = consesunet, reference_matrix = adjm, false_plot = F)

### Community detection

adj_comm <- community_path(adjm)
comm_consesusm <- community_path(consesusm)
comm_consesusu <- community_path(consesusu)
comm_consesunet <- community_path(consesunet)

topscore <- community_similarity(adj_comm,list(comm_consesusm, comm_consesusu, comm_consesunet))

keepscores <- topscore$topology_similarity %>% rownames_to_column("Predicted_Matrix") %>% mutate(Predicted_Matrix=c("vote", "union", "inet"), Method=c("GENIE3"), Ratio=nrow(adjm)/ncol(count_matrices[[1]]), p=nrow(adjm)) %>%
  full_join(keepscores)
keepscores <- topscore$community_metrics %>% rownames_to_column("Predicted_Matrix") %>% mutate(Predicted_Matrix=c("vote", "union", "inet"), Method=c("GENIE3"), Ratio=nrow(adjm)/ncol(count_matrices[[1]]), p=nrow(adjm)) %>%
  full_join(keepscores)

## Early integration

count_matrices <- lapply(count_matrices, as.matrix)
early_matrix <- list(earlyj(count_matrices, rowg = T))

time[["GENIE3_early_15Cores"]] <- system.time(
  early <- infer_networks(early_matrix, method="GENIE3", nCores = 15,
                           seed=seed_base+run_id)
)

### Symmetrize and ROC

early_wadj <- generate_adjacency(early)
searly_wadj <- symmetrize(early_wadj, weight_function = "mean")

early_auc <- plotROC(searly_wadj, adjm, plot_title = "ROC curve - GENIE3 Early Integration", is_binary = F)

### Cutoff

searly_adj <- cutoff_adjacency(count_matrices = early_matrix,
                               weighted_adjm_list = searly_wadj, 
                               n = 2,
                               method = "GENIE3",
                               nCores = 15)

scores.early <- pscores(adjm, searly_adj)

keepscores <- scores.early$Statistics %>% mutate(Predicted_Matrix=c("early"), Method=c("GENIE3"), Ratio=nrow(adjm)/ncol(count_matrices[[1]]), p=nrow(adjm)) %>%
  full_join(keepscores)

plotg(searly_adj)

### Plot comparison

compare_consensus(consensus_matrix = searly_adj[[1]], reference_matrix = adjm, false_plot = F)

### Community detection

comm_consesusm <- community_path(searly_adj[[1]])
topscore <- community_similarity(adj_comm,list(comm_consesusm))

community_metrics <- topscore$community_metrics["Predicted_1", ]
topology_metrics <- topscore$topology_similarity["Predicted_1", ]

early_row_index <- which(keepscores$Predicted_Matrix == "early" & keepscores$Method == "GENIE3")

keepscores$VI[early_row_index] <- community_metrics$VI
keepscores$NMI[early_row_index] <- community_metrics$NMI
keepscores$ARI[early_row_index] <- community_metrics$ARI

keepscores$Modularity[early_row_index]   <- topology_metrics$Modularity
keepscores$Communities[early_row_index]  <- topology_metrics$Communities
keepscores$Density[early_row_index]      <- topology_metrics$Density
keepscores$Transitivity[early_row_index] <- topology_metrics$Transitivity

# GRNBoost2
## Late integration
time[["GRNBoost2_late"]] <- system.time(
  late <- infer_networks(count_matrices, 
                         method="GRNBoost2",
                         grnboost_modules = modules,
                         seed=seed_base+run_id)
)

### Symmetrize and ROC

late_wadj <- generate_adjacency(late)
slate_wadj <- symmetrize(late_wadj, weight_function = "mean")

late_auc <- plotROC(slate_wadj, adjm, plot_title = "ROC curve - grnboost Late Integration")

### Cutoff

slate_adj <- cutoff_adjacency(count_matrices = count_matrices,
                              weighted_adjm_list = slate_wadj, 
                              n = 3,
                              method = "GRNBoost2",
                              grnboost_modules = modules)

scores.late.all <- pscores(adjm, slate_adj)

plotg(slate_adj)

### Consensus

consesusm <- create_consensus(slate_adj, method="vote")
consesusu <- create_consensus(slate_adj, method="union")
consesunet <- create_consensus(adj_matrix_list = slate_adj, weighted_list = slate_wadj, method = "INet", threshold = 0.05)
scores.late <- pscores(adjm, list(consesusm,consesusu,consesunet))

keepscores <- scores.late$Statistics %>% mutate(Predicted_Matrix=c("vote", "union", "inet"), Method=c("GRNBoost2"), Ratio=nrow(adjm)/ncol(count_matrices[[1]]), p=nrow(adjm)) %>%
  full_join(keepscores)

### Plot comparison

compare_consensus(consensus_matrix = consesusm, reference_matrix = adjm, false_plot = F)
compare_consensus(consensus_matrix = consesusu, reference_matrix = adjm, false_plot = F)
compare_consensus(consensus_matrix = consesunet, reference_matrix = adjm, false_plot = F)

### Community detection

adj_comm <- community_path(adjm)
comm_consesusm <- community_path(consesusm)
comm_consesusu <- community_path(consesusu)
comm_consesunet <- community_path(consesunet)

topscore <- community_similarity(adj_comm,list(comm_consesusm, comm_consesusu, comm_consesunet))

check  <- topscore$topology_similarity %>% rownames_to_column("Predicted_Matrix") %>% mutate(Predicted_Matrix=c("vote", "union", "inet"), Method=c("GRNBoost2"), Ratio=nrow(adjm)/ncol(count_matrices[[1]]), p=nrow(adjm))

# Join check to keepscores by Predicted_Matrix and Method (to be safe)
keepscores <- keepscores %>%
  left_join(check, by = c("Predicted_Matrix", "Method", "Ratio", "p"), suffix = c("", "_check")) %>%
  mutate(
    Modularity = ifelse(is.na(Modularity), Modularity_check, Modularity),
    Communities = ifelse(is.na(Communities), Communities_check, Communities),
    Density = ifelse(is.na(Density), Density_check, Density),
    Transitivity = ifelse(is.na(Transitivity), Transitivity_check, Transitivity)
  ) %>%
  dplyr::select(-ends_with("_check"))

check  <- topscore$community_metrics %>% rownames_to_column("Predicted_Matrix") %>% mutate(Predicted_Matrix=c("vote", "union", "inet"), Method=c("GRNBoost2"), Ratio=nrow(adjm)/ncol(count_matrices[[1]]), p=nrow(adjm))

# Join check to keepscores by Predicted_Matrix and Method (to be safe)
keepscores <- keepscores %>%
  left_join(check, by = c("Predicted_Matrix", "Method", "Ratio", "p"), suffix = c("", "_check")) %>%
  mutate(
    VI = ifelse(is.na(VI), VI_check, VI),
    NMI = ifelse(is.na(NMI), NMI_check, NMI),
    ARI = ifelse(is.na(ARI), ARI_check, ARI)
  ) %>%
  dplyr::select(-ends_with("_check"))

## Early integration

early_matrix <- list(earlyj(count_matrices))

time[["GRNBoost2_early"]] <- system.time(
  early <- infer_networks(early_matrix, 
                          method="GRNBoost2", 
                          grnboost_modules = modules,
                          seed=seed_base+run_id)
)

### Symmetrize and ROC

early_wadj <- generate_adjacency(early)
searly_wadj <- symmetrize(early_wadj, weight_function = "mean")

early_auc <- plotROC(searly_wadj, adjm, plot_title = "ROC curve - grnboost Early Integration")

### Cutoff

searly_adj <- cutoff_adjacency(count_matrices = early_matrix,
                               weighted_adjm_list = searly_wadj, 
                               n = 2,
                               method = "GRNBoost2",
                               grnboost_modules = modules
)

scores.early <- pscores(adjm, searly_adj)

keepscores <- scores.early$Statistics %>% mutate(Predicted_Matrix=c("early"), Method=c("GRNBoost2"), Ratio=nrow(adjm)/ncol(count_matrices[[1]]), p=nrow(adjm)) %>%
  full_join(keepscores)

plotg(searly_adj)

### Plot comparison

compare_consensus(consensus_matrix = searly_adj[[1]], reference_matrix = adjm, false_plot = F)

### Community detection

comm_consesusm <- community_path(searly_adj[[1]])
topscore <- community_similarity(adj_comm,list(comm_consesusm))

community_metrics <- topscore$community_metrics["Predicted_1", ]
topology_metrics <- topscore$topology_similarity["Predicted_1", ]

early_row_index <- which(keepscores$Predicted_Matrix == "early" & keepscores$Method == "GRNBoost2")

keepscores$VI[early_row_index] <- community_metrics$VI
keepscores$NMI[early_row_index] <- community_metrics$NMI
keepscores$ARI[early_row_index] <- community_metrics$ARI

keepscores$Modularity[early_row_index]   <- topology_metrics$Modularity
keepscores$Communities[early_row_index]  <- topology_metrics$Communities
keepscores$Density[early_row_index]      <- topology_metrics$Density
keepscores$Transitivity[early_row_index] <- topology_metrics$Transitivity

# Joint Integration

## Joint Random Forest

#https://cran.r-project.org/src/contrib/Archive/JRF/
#install.packages("/home/francescoc/Downloads/JRF_0.1-4.tar.gz", repos = NULL, type = "source")
time[["JRF_15Cores"]] <- system.time(
  jrf_mat <- infer_networks(count_matrices, method="JRF", nCores = 15,
                           seed=seed_base+run_id)
)

### Prepare the output

jrf_list <- list()

importance_columns <- grep("importance", names(jrf_mat[[1]]), value = TRUE)

for (i in seq_along(importance_columns)) {
  # Select the 'gene1', 'gene2', and the current 'importance' column
  df <- jrf_mat[[1]][, c("gene1", "gene2", importance_columns[i])]
  
  # Rename the importance column to its original name (e.g., importance1, importance2, etc.)
  names(df)[3] <- importance_columns[i]
  
  # Add the data frame to the output list
  jrf_list[[i]] <- df
}

### symmetrize Output and ROC

jrf_wadj <- generate_adjacency(jrf_list)
sjrf_wadj <- symmetrize(jrf_wadj, weight_function = "mean")
jrf_auc_mine <- plotROC(sjrf_wadj, adjm, plot_title = "ROC curve - JRF Late Integration", is_binary = F)

### Generate Adjacency and Apply Cutoff

sjrf_adj <- cutoff_adjacency(count_matrices = count_matrices,
                             weighted_adjm_list = sjrf_wadj, 
                             n = 3,
                             method = "JRF")

### Comparison with the Ground Truth

scores.jrf.all <- pscores(adjm, sjrf_adj)
plotg(sjrf_adj)

consesusm <- create_consensus(sjrf_adj, method="vote")
consesusu <- create_consensus(sjrf_adj, method="union")
consesunet <- create_consensus(adj_matrix_list = sjrf_adj, weighted_list = sjrf_wadj, method = "INet", threshold = 0.1, ncores = 15)

scores.jrf <- pscores(adjm, list(consesusm, consesusu, consesunet))

keepscores <- scores.jrf$Statistics %>% mutate(Predicted_Matrix=c("vote", "union", "inet"), Method=c("JRF"), Ratio=nrow(adjm)/ncol(count_matrices[[1]]), p=nrow(adjm)) %>%
  full_join(keepscores)

compare_consensus(consensus_matrix = consesusm, reference_matrix = adjm, false_plot = F)
compare_consensus(consensus_matrix = consesusu, reference_matrix = adjm, false_plot = F)
compare_consensus(consensus_matrix = consesunet, reference_matrix = adjm, false_plot = F)

### Community detection

comm_consesusm <- community_path(consesusm)
comm_consesusu <- community_path(consesusu)
comm_consesunet <- community_path(consesunet)

community_similarity(adj_comm,list(comm_consesusm, comm_consesusu, comm_consesunet))

topscore <- community_similarity(adj_comm,list(comm_consesusm, comm_consesusu, comm_consesunet))

check  <- topscore$topology_similarity %>% rownames_to_column("Predicted_Matrix") %>% mutate(Predicted_Matrix=c("vote", "union", "inet"), Method=c("JRF"), Ratio=nrow(adjm)/ncol(count_matrices[[1]]), p=nrow(adjm))

# Join check to keepscores by Predicted_Matrix and Method (to be safe)
keepscores <- keepscores %>%
  left_join(check, by = c("Predicted_Matrix", "Method", "Ratio", "p"), suffix = c("", "_check")) %>%
  mutate(
    Modularity = ifelse(is.na(Modularity), Modularity_check, Modularity),
    Communities = ifelse(is.na(Communities), Communities_check, Communities),
    Density = ifelse(is.na(Density), Density_check, Density),
    Transitivity = ifelse(is.na(Transitivity), Transitivity_check, Transitivity)
  ) %>%
  dplyr::select(-ends_with("_check"))


check  <- topscore$community_metrics %>% rownames_to_column("Predicted_Matrix") %>% mutate(Predicted_Matrix=c("vote", "union", "inet"), Method=c("JRF"), Ratio=nrow(adjm)/ncol(count_matrices[[1]]), p=nrow(adjm))

# Join check to keepscores by Predicted_Matrix and Method (to be safe)
keepscores <- keepscores %>%
  left_join(check, by = c("Predicted_Matrix", "Method", "Ratio", "p"), suffix = c("", "_check")) %>%
  mutate(
    VI = ifelse(is.na(VI), VI_check, VI),
    NMI = ifelse(is.na(NMI), NMI_check, NMI),
    ARI = ifelse(is.na(ARI), ARI_check, ARI)
  ) %>%
  dplyr::select(-ends_with("_check"))

# Method Comparison

time_data <- data.frame(
  Method = names(time),
  Time_in_Seconds = sapply(time, function(x) if ("elapsed" %in% names(x)) x["elapsed"] else NA)
) %>%
  mutate(
    Time_in_Minutes = Time_in_Seconds / 60,
    Time_in_Hours = Time_in_Seconds / 3600,
    Ratio = nrow(adjm) / ncol(count_matrices[[1]]),
    p = nrow(adjm)
  ) %>%
  arrange(Time_in_Hours) %>%
  mutate(Method = factor(Method, levels = Method)) %>%
  separate(Method, into = c("Method", "Predicted_Matrix", "Cores"), sep = "_", fill = "right") %>%
  mutate(
    Predicted_Matrix = ifelse(Method == "JRF", "joint", Predicted_Matrix),
    Cores = case_when(
      Method %in% c("JRF", "GRNBoost2") ~ "15Cores",
      TRUE ~ Cores
    )
  )

df2 <- keepscores %>%
  mutate(
    Predicted_Matrix_Mapped = case_when(
      Method == "JRF" ~ "joint",
      Predicted_Matrix == "early" ~ "early",
      Predicted_Matrix %in% c("vote", "union", "inet") ~ "late",
      TRUE ~ Predicted_Matrix
    )
  )

df1_unique <- time_data %>%
  dplyr::select(Method, Predicted_Matrix, Time_in_Seconds, Time_in_Minutes, Time_in_Hours) %>%
  distinct()

df2 <- df2 %>%
  left_join(df1_unique, by = c("Method", "Predicted_Matrix_Mapped" = "Predicted_Matrix")) %>%
  dplyr::select(-Predicted_Matrix_Mapped)

return(df2)
}

all_runs <- lapply(1:10, function(i) run_single_simulation(run_id = i))
all_df <- bind_rows(all_runs, .id = "Run")

summary_df <- all_df %>%
  group_by(Method, Predicted_Matrix) %>%
  summarise(across(where(is.numeric), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"), .groups = "drop")

write.table(summary_df, file = file.path(ddir, "simp100_10runs.txt"), sep = "\t", quote = FALSE, row.names = FALSE)



