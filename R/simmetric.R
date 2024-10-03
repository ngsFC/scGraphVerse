simmetric <- function(df_list, weight_function = "mean") {
  
  # Internal function to process each dataframe
  process_dataframe <- function(df, weight_function) {
    # Ensure the dataframe has the expected column names
    if (!all(c("regulatoryGene", "targetGene", "weight") %in% colnames(df))) {
      stop("The dataframe must have columns: 'regulatoryGene', 'targetGene', 'weight'")
    }
    
    # Create a unique key for each gene pair, ignoring order (e.g., "Gene1-Gene2" and "Gene2-Gene1" are the same)
    df$pair <- apply(df[, c("regulatoryGene", "targetGene")], 1, function(x) paste(sort(x), collapse = "-"))
    
    # Aggregate the weight for each unique gene pair
    aggregated_df <- df %>%
      group_by(pair) %>%
      summarise(
        Gene1 = unique(sort(c(regulatoryGene, targetGene)))[1],  # Get the first gene in alphabetical order
        Gene2 = unique(sort(c(regulatoryGene, targetGene)))[2],  # Get the second gene in alphabetical order
        weight = match.fun(weight_function)(weight)  # Apply the specified aggregation function (mean, max, min)
      ) %>%
      ungroup()
    
    # Select the relevant columns for the final result
    result_df <- aggregated_df %>%
      select(Gene1, Gene2, weight)  # Keep only the relevant columns

    return(result_df)
  }
  
  # Apply the process_dataframe function to each dataframe in the list
  processed_df_list <- lapply(df_list, process_dataframe, weight_function = weight_function)
  
  return(processed_df_list)
}
