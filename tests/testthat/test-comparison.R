test_that("cutoff_adjacency and create_consensus work correctly", {
  set.seed(42)
  mat1 <- matrix(rpois(100, lambda = 5), nrow = 10)
  mat2 <- matrix(rpois(100, lambda = 5), nrow = 10)
  rownames(mat1) <- rownames(mat2) <- paste0("Gene", 1:10)

  adjm1 <- matrix(runif(100), nrow = 10)
  adjm2 <- matrix(runif(100), nrow = 10)
  rownames(adjm1) <- colnames(adjm1) <- paste0("Gene", 1:10)
  rownames(adjm2) <- colnames(adjm2) <- paste0("Gene", 1:10)
  symm <- symmetrize(list(adjm1, adjm2), weight_function = "mean")
  bin_adj_list <- cutoff_adjacency(
    count_matrices = list(mat1, mat2),
    weighted_adjm_list = symm,
    n = 2,
    method = "GENIE3",
    quantile_threshold = 0.9,
    nCores = 1
  )

  expect_true(is.list(bin_adj_list))
  expect_equal(length(bin_adj_list), 2)
  expect_equal(dim(bin_adj_list[[1]]), c(10, 10))

  consensus <- create_consensus(bin_adj_list, method = "vote", threshold = 0.5)
  expect_true(is.matrix(consensus))
  expect_equal(dim(consensus), c(10, 10))

  comparison <- compare_consensus(consensus_matrix = consensus, reference_matrix = consensus, false_plot = FALSE)
  expect_true(inherits(comparison, "gg"))
})
