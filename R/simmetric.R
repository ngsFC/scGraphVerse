simmetric <- function(df_list, weight_function = "mean") {
  
  process_dataframe <- function(df, weight_function) {
    if (ncol(df) < 3) {
      stop("The dataframe must have three columns. Gene1-Gene2-weight")
    }
    
    colnames(df)[1:3] <- c("G1", "G2", "w")
    df$pair <- apply(df[, c("G1", "G2")], 1, function(x) paste(sort(x), collapse = "-"))
    
    aggregated_df <- df %>%
      group_by(pair) %>%
      summarise(
        Gene1 = unique(sort(c(G1, G2)))[1],  # Get the first gene in alphabetical order
        Gene2 = unique(sort(c(G1, G2)))[2],  # Get the second gene in alphabetical order
        weight = match.fun(weight_function)(w)  # Apply the specified aggregation function (mean, max, min)
      ) %>%
      ungroup()
    
    result_df <- aggregated_df %>%
      select(Gene1, Gene2, weight)
    
    return(result_df)
  }
  
  processed_df_list <- lapply(df_list, process_dataframe, weight_function = weight_function)
  
  return(processed_df_list)
}