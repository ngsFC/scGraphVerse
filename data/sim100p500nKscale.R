run_single_simulation <- function(run_id = 1, count_matrices, seed_base = 1234, k = NA) {
  set.seed(seed_base + run_id)
  
  library(tidyverse)
  library(scGraphVerse)
  modules <- init_py(python_path = "/usr/bin/python3", required = TRUE)
  
  time <- list()
  
  adjm <- as.matrix(read.table("./../analysis/simulation/adjacency/adjm_top1200_p500.txt"))
  colnames(adjm) <- rownames(adjm)
  
  adj_comm <- community_path(adjm)
  
  if (k == 1) {
    message("Only 1 matrix: skipping consensus and JRF. Early = Late.")
    single_matrix <- list(count_matrices[[1]])
    
    # GENIE3
    time[["GENIE3_late_15Cores"]] <- system.time(
      net <- infer_networks(single_matrix, method = "GENIE3", nCores = 15)
    )
    wadj <- generate_adjacency(net)
    swadj <- symmetrize(wadj, weight_function = "mean")
    adj <- cutoff_adjacency(single_matrix, swadj, n = 1, method = "GENIE3", nCores = 15)
    stats <- pscores(adjm, adj)$Statistics %>% mutate(Method = "GENIE3")
    comm <- community_path(adj[[1]])
    cs <- community_similarity(adj_comm, list(comm))
    stats[, c("VI", "NMI", "ARI")] <- cs$community_metrics["Predicted_1", c("VI", "NMI", "ARI")]
    stats[, c("Modularity", "Communities", "Density", "Transitivity")] <- cs$topology_similarity["Predicted_1", c("Modularity", "Communities", "Density", "Transitivity")]
    stats_genie3 <- bind_rows(
      stats %>% mutate(Predicted_Matrix = "late"),
      stats %>% mutate(Predicted_Matrix = "early")
    )
    time[["GENIE3_early_15Cores"]] <- time[["GENIE3_late_15Cores"]]
    
    # GRNBoost2
    time[["GRNBoost2_late"]] <- system.time(
      net <- infer_networks(single_matrix, method = "GRNBoost2", nCores = 15, grnboost_modules = modules)
    )
    wadj <- generate_adjacency(net)
    swadj <- symmetrize(wadj, weight_function = "mean")
    adj <- cutoff_adjacency(single_matrix, swadj, n = 1, method = "GRNBoost2", nCores = 15, grnboost_modules = modules)
    stats <- pscores(adjm, adj)$Statistics %>% mutate(Method = "GRNBoost2")
    comm <- community_path(adj[[1]])
    cs <- community_similarity(adj_comm, list(comm))
    stats[, c("VI", "NMI", "ARI")] <- cs$community_metrics["Predicted_1", c("VI", "NMI", "ARI")]
    stats[, c("Modularity", "Communities", "Density", "Transitivity")] <- cs$topology_similarity["Predicted_1", c("Modularity", "Communities", "Density", "Transitivity")]
    stats_grnboost2 <- bind_rows(
      stats %>% mutate(Predicted_Matrix = "late"),
      stats %>% mutate(Predicted_Matrix = "early")
    )
    time[["GRNBoost2_early"]] <- time[["GRNBoost2_late"]]
    
    keepscores <- bind_rows(stats_genie3, stats_grnboost2)
    
  } else {
    early_matrix <- list(earlyj(count_matrices, rowg = TRUE))
    late_matrix <- count_matrices
    
    ########## GENIE3 ##########
    time[["GENIE3_late_15Cores"]] <- system.time(
      late <- infer_networks(late_matrix, method = "GENIE3", nCores = 15)
    )
    late_wadj <- generate_adjacency(late)
    slate_wadj <- symmetrize(late_wadj, weight_function = "mean")
    slate_adj <- cutoff_adjacency(count_matrices = late_matrix, weighted_adjm_list = slate_wadj, n = 3, method = "GENIE3", nCores = 15)
    stats_late <- pscores(adjm, list(
      vote = create_consensus(slate_adj, method = "vote"),
      union = create_consensus(slate_adj, method = "union"),
      inet = create_consensus(slate_adj, slate_wadj, method = "INet", threshold = 0.05, ncores = 15)
    ))$Statistics %>%
      mutate(Predicted_Matrix = c("vote", "union", "inet"), Method = "GENIE3")
    
    time[["GENIE3_early_15Cores"]] <- system.time(
      early <- infer_networks(early_matrix, method = "GENIE3", nCores = 15)
    )
    early_wadj <- generate_adjacency(early)
    searly_wadj <- symmetrize(early_wadj, weight_function = "mean")
    searly_adj <- cutoff_adjacency(count_matrices = early_matrix, weighted_adjm_list = searly_wadj, n = 2, method = "GENIE3", nCores = 15)
    stats_early <- pscores(adjm, searly_adj)$Statistics %>%
      mutate(Predicted_Matrix = "early", Method = "GENIE3")
    
    keepscores <- bind_rows(stats_late, stats_early)
    
    ########## GRNBoost2 ##########
    time[["GRNBoost2_late"]] <- system.time(
      late <- infer_networks(late_matrix, method = "GRNBoost2", nCores = 15, grnboost_modules = modules)
    )
    late_wadj <- generate_adjacency(late)
    slate_wadj <- symmetrize(late_wadj, weight_function = "mean")
    slate_adj <- cutoff_adjacency(late_matrix, slate_wadj, n = 3, method = "GRNBoost2", nCores = 15, grnboost_modules = modules)
    stats_late <- pscores(adjm, list(
      vote = create_consensus(slate_adj, method = "vote"),
      union = create_consensus(slate_adj, method = "union"),
      inet  = create_consensus(slate_adj, slate_wadj, method = "INet", threshold = 0.05)
    ))$Statistics %>%
      mutate(Predicted_Matrix = c("vote", "union", "inet"), Method = "GRNBoost2")
    
    time[["GRNBoost2_early"]] <- system.time(
      early <- infer_networks(early_matrix, method = "GRNBoost2", nCores = 15, grnboost_modules = modules)
    )
    early_wadj <- generate_adjacency(early)
    searly_wadj <- symmetrize(early_wadj, weight_function = "mean")
    searly_adj <- cutoff_adjacency(early_matrix, searly_wadj, n = 2, method = "GRNBoost2", nCores = 15, grnboost_modules = modules)
    stats_early <- pscores(adjm, searly_adj)$Statistics %>%
      mutate(Predicted_Matrix = "early", Method = "GRNBoost2")
    
    keepscores <- bind_rows(keepscores, stats_late, stats_early)
    
    ########## JRF ##########
    time[["JRF_15Cores"]] <- system.time(
      jrf_mat <- infer_networks(count_matrices, method = "JRF", nCores = 15)
    )
    jrf_list <- lapply(grep("importance", names(jrf_mat[[1]]), value = TRUE), function(col) {
      df <- jrf_mat[[1]][, c("gene1", "gene2", col)]
      colnames(df)[3] <- col
      df
    })
    jrf_wadj <- generate_adjacency(jrf_list)
    sjrf_wadj <- symmetrize(jrf_wadj, weight_function = "mean")
    sjrf_adj <- cutoff_adjacency(count_matrices, sjrf_wadj, n = 3, method = "JRF")
    stats_jrf <- pscores(adjm, list(
      vote  = create_consensus(sjrf_adj, method = "vote"),
      union = create_consensus(sjrf_adj, method = "union"),
      inet  = create_consensus(sjrf_adj, sjrf_wadj, method = "INet", threshold = 0.1, ncores = 15)
    ))$Statistics %>%
      mutate(Predicted_Matrix = c("vote", "union", "inet"), Method = "JRF")
    
    keepscores <- bind_rows(keepscores, stats_jrf)
  }
  
  timing_df <- data.frame(
    Method = names(time),
    Time_in_Seconds = sapply(time, function(x) x["elapsed"])
  ) %>%
    mutate(
      Time_in_Minutes = Time_in_Seconds / 60,
      Time_in_Hours = Time_in_Seconds / 3600,
      Method = gsub("_15Cores", "", Method),
      Predicted_Matrix = case_when(
        grepl("late", names(time)) ~ "late",
        grepl("early", names(time)) ~ "early",
        grepl("JRF", names(time)) ~ "joint",
        TRUE ~ "unknown"
      )
    )
  
  keepscores <- keepscores %>%
    left_join(timing_df, by = c("Method", "Predicted_Matrix")) %>%
    mutate(n_matrices = k)
  
  return(keepscores)
}

library(tidyverse)
library(glue)

all_counts <- readRDS("/home/francescoc/Desktop/scGraphVerse/analysis/simulation/simdata/sim_n100p500k5.RDS")
all_counts <- lapply(all_counts, t)

all_results <- list()

for (k in 1:5) {
  count_subset <- all_counts[1:k]
  runs_k <- lapply(1:10, function(i) {
    message(glue("Running k = {k}, run = {i}"))
    run_single_simulation(run_id = i, count_matrices = count_subset, k = k)
  })
  all_results[[as.character(k)]] <- bind_rows(runs_k, .id = "Run")
}

final_df <- bind_rows(all_results, .id = "Subset_Size")

write.table(final_df, file = file.path("/home/francescoc/Desktop/scGraphVerse/analysis/simulation/results/", "sim_n100p500k5_k1to5_10runs.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)


for (k in 1:2) {
  count_subset <- all_counts[1:k]
  runs_k <- lapply(1:2, function(i) {
    message(glue("Running k = {k}, run = {i}"))
    run_single_simulation(run_id = i, count_matrices = count_subset, k = k)
  })
  all_results[[as.character(k)]] <- bind_rows(runs_k, .id = "Run")
}

final_df <- bind_rows(all_results, .id = "Subset_Size")
