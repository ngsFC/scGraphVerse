test_that("infer_networks and generate_adjacency produce valid matrices", {
  set.seed(123)
  mat <- matrix(rpois(100, lambda = 5), nrow = 10)
  rownames(mat) <- paste0("Gene", 1:10)

  inferred <- infer_networks(count_matrices_list = list(mat), method = "GENIE3", nCores = 1)
  adjm <- generate_adjacency(inferred)
  symm <- symmetrize(adjm, weight_function = "mean")[[1]]

  expect_true(is.matrix(symm))
  expect_equal(dim(symm), c(10, 10))
  expect_true(all(symm >= 0))
  expect_true(all(abs(symm - t(symm)) < 1e-6)) # Symmetric
})
