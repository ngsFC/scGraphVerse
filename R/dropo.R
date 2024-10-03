dropo <- function(matrix, dropout_rate) {
  total_elements <- length(matrix)
  num_dropout <- round(dropout_rate * total_elements)
  dropout_indices <- sample(1:total_elements, num_dropout)
  matrix[dropout_indices] <- 0
  return(matrix)
}
