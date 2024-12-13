exploreCells <- function(seurat_object, label_column = "cell_type") {
  if (!label_column %in% colnames(seurat_object@meta.data)) {
    stop("The specified label column does not exist in the metadata.")
  }
  
  cell_types <- seurat_object@meta.data[[label_column]]
  
  cell_type_counts <- as.data.frame(table(cell_types), stringsAsFactors = FALSE)
  colnames(cell_type_counts) <- c("Cell_Type", "Cell_Count")
  cell_type_counts <- cell_type_counts[order(-cell_type_counts$Cell_Count), ]
  
  message("Summary of cell types and their counts:")
  print(cell_type_counts)
  
  return(cell_type_counts)
}

