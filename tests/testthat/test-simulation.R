test_that("zinb_simdata generates count matrices with correct dimensions", {
  B <- diag(5)
  rownames(B) <- colnames(B) <- paste0("Gene", 1:5)

  sim_matrices <- zinb_simdata(
    n = 100, p = 5, B = B,
    mu_range = list(c(1, 5)),
    mu_noise = 0.5,
    theta = 1,
    pi = 0.2,
    kmat = 1,
    depth_range = c(1000, 5000),
    seed = 42
  )

  expect_true(is.list(sim_matrices))
  expect_equal(length(sim_matrices), 1)
  expect_equal(dim(sim_matrices[[1]]), c(100, 5))
  expect_true(all(sim_matrices[[1]] >= 0))
})

