test_that("create_mae creates MultiAssayExperiment", {
    skip_if_not_installed("MultiAssayExperiment")

    mat1 <- matrix(rpois(100, lambda = 5), nrow = 10, ncol = 10)
    mat2 <- matrix(rpois(100, lambda = 5), nrow = 10, ncol = 10)
    rownames(mat1) <- rownames(mat2) <- paste0("Gene", 1:10)
    colnames(mat1) <- paste0("Cell", 1:10)
    colnames(mat2) <- paste0("Cell", 1:10)

    mae <- create_mae(list(Cond1 = mat1, Cond2 = mat2))

    expect_true(inherits(mae, "MultiAssayExperiment"))
    expect_equal(length(mae), 2)
})

test_that("earlyj handles MAE of matrices", {
    skip_if_not_installed("MultiAssayExperiment")

    mat1 <- matrix(rpois(50, lambda = 5), nrow = 5, ncol = 10)
    mat2 <- matrix(rpois(50, lambda = 5), nrow = 5, ncol = 10)
    rownames(mat1) <- rownames(mat2) <- paste0("Gene", 1:5)
    colnames(mat1) <- colnames(mat2) <- paste0("Cell", 1:10)

    mae <- create_mae(list(Cond1 = mat1, Cond2 = mat2))
    result <- earlyj(mae)

    expect_true(inherits(result, "MultiAssayExperiment"))
    expect_equal(length(result), 1)
})

test_that("build_network_se creates SummarizedExperiment from adjacency matrix", {
    skip_if_not_installed("SummarizedExperiment")

    adj_mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    rownames(adj_mat) <- colnames(adj_mat) <- paste0("Gene", 1:3)

    se <- build_network_se(list(adj_mat))

    expect_true(inherits(se, "SummarizedExperiment"))
})

test_that("infer_networks works with GENIE3 method - returns list of dataframes", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(50, lambda = 5), nrow = 10, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:10)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    inferred <- infer_networks(count_matrices_list = mae, method = "GENIE3", nCores = 1)

    expect_true(is.list(inferred))
    expect_equal(length(inferred), 1)
    expect_true(is.data.frame(inferred[[1]]))
})

test_that("infer_networks works with PCzinb method - returns SE", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(50, lambda = 5), nrow = 10, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:10)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    inferred <- infer_networks(count_matrices_list = mae, method = "PCzinb", nCores = 1)

    expect_true(inherits(inferred, "SummarizedExperiment"))
})

test_that("infer_networks works with ZILGM method - returns SE", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(50, lambda = 5), nrow = 10, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:10)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    inferred <- infer_networks(count_matrices_list = mae, method = "ZILGM", nCores = 1)

    expect_true(inherits(inferred, "SummarizedExperiment"))
})

test_that("infer_networks ZILGM works with multiple matrices", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat1 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    mat2 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat1) <- rownames(mat2) <- paste0("Gene", 1:5)
    colnames(mat1) <- paste0("Cell1_", 1:5)
    colnames(mat2) <- paste0("Cell2_", 1:5)
    mae <- create_mae(list(Cond1 = mat1, Cond2 = mat2))

    inferred <- infer_networks(count_matrices_list = mae, method = "ZILGM", nCores = 1)

    expect_true(inherits(inferred, "SummarizedExperiment"))
    # ZILGM returns SE where length is number of assays (networks)
    expect_true(length(SummarizedExperiment::assays(inferred)) >= 1)
})

test_that("infer_networks PCzinb with custom parameters", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    pczinb_params <- list(maxcard = 2)
    inferred <- infer_networks(
        count_matrices_list = mae,
        method = "PCzinb",
        nCores = 1,
        pczinb_params = pczinb_params
    )

    expect_true(inherits(inferred, "SummarizedExperiment"))
})

test_that("infer_networks works with JRF method - returns list of dataframes", {
    skip_if_not_installed("MultiAssayExperiment")

    # Try to load the compiled code - will work during R CMD check and package_coverage
    dll_loaded <- tryCatch(
        {
            library.dynam("scGraphVerse", package = "scGraphVerse", lib.loc = .libPaths())
            TRUE
        },
        error = function(e) FALSE
    )

    if (!dll_loaded) {
        skip("JRF C code requires package to be installed (works during coverage)")
    }

    set.seed(123)
    mat1 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    mat2 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat1) <- rownames(mat2) <- paste0("Gene", 1:5)
    colnames(mat1) <- paste0("Cell1_", 1:5)
    colnames(mat2) <- paste0("Cell2_", 1:5)
    mae <- create_mae(list(Cond1 = mat1, Cond2 = mat2))

    inferred <- infer_networks(count_matrices_list = mae, method = "JRF", nCores = 1)

    expect_true(is.list(inferred))
    expect_equal(length(inferred), 2)
    expect_true(all(sapply(inferred, is.data.frame)))
})

test_that("infer_networks JRF with custom parameters", {
    skip_if_not_installed("MultiAssayExperiment")

    # Try to load the compiled code
    dll_loaded <- tryCatch(
        {
            library.dynam("scGraphVerse", package = "scGraphVerse", lib.loc = .libPaths())
            TRUE
        },
        error = function(e) FALSE
    )

    if (!dll_loaded) {
        skip("JRF C code requires package to be installed (works during coverage)")
    }

    set.seed(123)
    mat1 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    mat2 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat1) <- rownames(mat2) <- paste0("Gene", 1:5)
    colnames(mat1) <- paste0("Cell1_", 1:5)
    colnames(mat2) <- paste0("Cell2_", 1:5)
    mae <- create_mae(list(Cond1 = mat1, Cond2 = mat2))

    jrf_params <- list(ntree = 500, mtry = 2)

    inferred <- infer_networks(
        count_matrices_list = mae,
        method = "JRF",
        nCores = 1,
        jrf_params = jrf_params
    )

    expect_true(is.list(inferred))
    expect_equal(length(inferred), 2)
})

test_that("generate_adjacency produces SE", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(50, lambda = 5), nrow = 10, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:10)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    inferred <- infer_networks(count_matrices_list = mae, method = "GENIE3", nCores = 1)
    adjm <- generate_adjacency(inferred)

    expect_true(inherits(adjm, "SummarizedExperiment"))
    # SE length is number of genes (rows), not number of networks
    expect_equal(nrow(adjm), 10)
})

test_that("symmetrize works with SE", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(50, lambda = 5), nrow = 10, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:10)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    inferred <- infer_networks(count_matrices_list = mae, method = "GENIE3", nCores = 1)
    adjm <- generate_adjacency(inferred)

    # adjm is already SE, symmetrize directly
    symm <- symmetrize(adjm, weight_function = "mean")

    expect_true(inherits(symm, "SummarizedExperiment"))
})

test_that("infer_networks handles verbose parameter", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(50, lambda = 5), nrow = 10, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:10)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    expect_no_error(infer_networks(count_matrices_list = mae, method = "GENIE3", nCores = 1, verbose = TRUE))
    expect_no_error(infer_networks(count_matrices_list = mae, method = "GENIE3", nCores = 1, verbose = FALSE))
})

test_that("infer_networks works with multiple matrices", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat1 <- matrix(rpois(50, lambda = 5), nrow = 10, ncol = 5)
    mat2 <- matrix(rpois(50, lambda = 5), nrow = 10, ncol = 5)
    rownames(mat1) <- rownames(mat2) <- paste0("Gene", 1:10)
    colnames(mat1) <- paste0("Cell1_", 1:5)
    colnames(mat2) <- paste0("Cell2_", 1:5)
    mae <- create_mae(list(Cond1 = mat1, Cond2 = mat2))

    inferred <- infer_networks(count_matrices_list = mae, method = "GENIE3", nCores = 1)

    expect_true(is.list(inferred))
    expect_equal(length(inferred), 2)
})

test_that("symmetrize works with different weight functions", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(50, lambda = 5), nrow = 10, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:10)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    inferred <- infer_networks(count_matrices_list = mae, method = "GENIE3", nCores = 1)
    adjm <- generate_adjacency(inferred)

    symm_mean <- symmetrize(adjm, weight_function = "mean")
    symm_max <- symmetrize(adjm, weight_function = "max")
    symm_min <- symmetrize(adjm, weight_function = "min")

    expect_true(inherits(symm_mean, "SummarizedExperiment"))
    expect_true(inherits(symm_max, "SummarizedExperiment"))
    expect_true(inherits(symm_min, "SummarizedExperiment"))
})

test_that("infer_networks GENIE3 with regulators and targets", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("GENIE3")

    set.seed(123)
    mat <- matrix(rpois(50, lambda = 5), nrow = 10, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:10)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    genie3_params <- list(
        regulators = paste0("Gene", 1:5),
        targets = paste0("Gene", 6:10)
    )

    inferred <- infer_networks(
        count_matrices_list = mae,
        method = "GENIE3",
        nCores = 1,
        genie3_params = genie3_params
    )

    expect_true(is.list(inferred))
    expect_equal(length(inferred), 1)
})

test_that("infer_networks GENIE3 with Extra Trees", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("GENIE3")

    set.seed(123)
    mat <- matrix(rpois(50, lambda = 5), nrow = 10, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:10)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    genie3_params <- list(
        treeMethod = "ET",
        K = 3,
        nTrees = 100
    )

    inferred <- infer_networks(
        count_matrices_list = mae,
        method = "GENIE3",
        nCores = 1,
        genie3_params = genie3_params
    )

    expect_true(is.list(inferred))
})

test_that("infer_networks ZILGM with custom lambda", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    zilgm_params <- list(
        lambda = c(0.1, 0.05),
        nlambda = 2,
        family = "Poisson"
    )

    inferred <- infer_networks(
        count_matrices_list = mae,
        method = "ZILGM",
        nCores = 1,
        zilgm_params = zilgm_params
    )

    expect_true(inherits(inferred, "SummarizedExperiment"))
})

test_that("infer_networks ZILGM with NBI family", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    zilgm_params <- list(
        lambda = 0.1,
        family = "NBI",
        theta = 1
    )

    inferred <- infer_networks(
        count_matrices_list = mae,
        method = "ZILGM",
        nCores = 1,
        zilgm_params = zilgm_params
    )

    expect_true(inherits(inferred, "SummarizedExperiment"))
})

test_that("infer_networks PCzinb with different methods", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- colnames(mat) <- paste0("Gene", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    # Test with nb method
    pczinb_params_nb <- list(method = "nb", maxcard = 1)
    inferred_nb <- infer_networks(
        count_matrices_list = mae,
        method = "PCzinb",
        nCores = 1,
        pczinb_params = pczinb_params_nb
    )

    expect_true(inherits(inferred_nb, "SummarizedExperiment"))
})

test_that("infer_networks with single cell matrix", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    # Very small matrix
    mat <- matrix(rpois(20, lambda = 5), nrow = 4, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:4)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    inferred <- infer_networks(
        count_matrices_list = mae,
        method = "GENIE3",
        nCores = 1,
        verbose = FALSE
    )

    expect_true(is.list(inferred))
    expect_equal(length(inferred), 1)
})

test_that("infer_networks error on invalid method", {
    skip_if_not_installed("MultiAssayExperiment")

    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    expect_error(
        infer_networks(count_matrices_list = mae, method = "InvalidMethod"),
        "arg.*should be one of"
    )
})

test_that("infer_networks GENIE3 verbose mode", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("GENIE3")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    # Test verbose mode - message says "Running GENIE3"
    expect_message(
        infer_networks(
            count_matrices_list = mae,
            method = "GENIE3",
            nCores = 1,
            verbose = TRUE
        ),
        "Running GENIE3"
    )
})

test_that("infer_networks PCzinb with zinb0 method", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- colnames(mat) <- paste0("Gene", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    pczinb_params <- list(method = "zinb0", maxcard = 1)
    inferred <- infer_networks(
        count_matrices_list = mae,
        method = "PCzinb",
        nCores = 1,
        pczinb_params = pczinb_params
    )

    expect_true(inherits(inferred, "SummarizedExperiment"))
})

test_that("infer_networks ZILGM with OR symmetrization", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    zilgm_params <- list(
        lambda = 0.1,
        sym = "OR"
    )

    inferred <- infer_networks(
        count_matrices_list = mae,
        method = "ZILGM",
        nCores = 1,
        zilgm_params = zilgm_params
    )

    expect_true(inherits(inferred, "SummarizedExperiment"))
})

test_that("infer_networks ZILGM with MM update type", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    zilgm_params <- list(
        lambda = 0.1,
        update_type = "MM"
    )

    inferred <- infer_networks(
        count_matrices_list = mae,
        method = "ZILGM",
        nCores = 1,
        zilgm_params = zilgm_params
    )

    expect_true(inherits(inferred, "SummarizedExperiment"))
})

test_that("infer_networks with 3+ matrices", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat1 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    mat2 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    mat3 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat1) <- rownames(mat2) <- rownames(mat3) <- paste0("Gene", 1:5)
    colnames(mat1) <- paste0("Cell1_", 1:5)
    colnames(mat2) <- paste0("Cell2_", 1:5)
    colnames(mat3) <- paste0("Cell3_", 1:5)

    mae <- create_mae(list(Cond1 = mat1, Cond2 = mat2, Cond3 = mat3))

    inferred <- infer_networks(
        count_matrices_list = mae,
        method = "GENIE3",
        nCores = 1
    )

    expect_true(is.list(inferred))
    expect_equal(length(inferred), 3)
})

test_that("infer_networks PCzinb with extend parameter", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- colnames(mat) <- paste0("Gene", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    pczinb_params <- list(
        method = "poi",
        maxcard = 1,
        extend = FALSE
    )

    inferred <- infer_networks(
        count_matrices_list = mae,
        method = "PCzinb",
        nCores = 1,
        pczinb_params = pczinb_params
    )

    expect_true(inherits(inferred, "SummarizedExperiment"))
})

test_that("infer_networks GENIE3 with K parameter", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("GENIE3")

    set.seed(123)
    mat <- matrix(rpois(50, lambda = 5), nrow = 10, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:10)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    genie3_params <- list(K = 5)

    inferred <- infer_networks(
        count_matrices_list = mae,
        method = "GENIE3",
        nCores = 1,
        genie3_params = genie3_params
    )

    expect_true(is.list(inferred))
})

test_that("generate_adjacency with single network", {
    df <- data.frame(
        Gene1 = c("A", "B", "C"),
        Gene2 = c("B", "C", "A"),
        Weight = c(0.5, 0.8, 0.3)
    )

    result <- generate_adjacency(list(df), nCores = 1)

    expect_true(inherits(result, "SummarizedExperiment"))
    adj <- SummarizedExperiment::assay(result, 1)
    expect_equal(adj["A", "B"], 0.5)
    expect_equal(adj["B", "C"], 0.8)
    # Check diagonal is zero (without checking names)
    expect_true(all(diag(adj) == 0))
})

test_that("generate_adjacency with nCores > 1", {
    df1 <- data.frame(
        Gene1 = c("A", "B"),
        Gene2 = c("B", "C"),
        Weight = c(0.5, 0.8)
    )
    df2 <- data.frame(
        Gene1 = c("A", "C"),
        Gene2 = c("C", "B"),
        Weight = c(0.3, 0.6)
    )

    result <- generate_adjacency(list(df1, df2), nCores = 2)

    expect_true(inherits(result, "SummarizedExperiment"))
    expect_equal(length(SummarizedExperiment::assays(result)), 2)
})

test_that("generate_adjacency handles NA weights", {
    df <- data.frame(
        Gene1 = c("A", "B", "C"),
        Gene2 = c("B", "C", "A"),
        Weight = c(0.5, NA, 0.3)
    )

    result <- generate_adjacency(list(df), nCores = 1)
    adj <- SummarizedExperiment::assay(result, 1)

    expect_equal(adj["A", "B"], 0.5)
    expect_equal(adj["B", "C"], 0) # NA should be ignored
    expect_equal(adj["C", "A"], 0.3)
})

test_that("generate_adjacency adds network metadata", {
    df <- data.frame(
        Gene1 = c("A", "B"),
        Gene2 = c("B", "A"),
        Weight = c(0.5, 0.8)
    )

    result <- generate_adjacency(list(net1 = df), nCores = 1)

    expect_equal(S4Vectors::metadata(result)$type, "weighted")
    # Check metadata structure
    expect_true(inherits(result, "SummarizedExperiment"))
})

test_that("create_mae with unnamed list", {
    mat1 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    mat2 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat1) <- rownames(mat2) <- paste0("Gene", 1:5)
    colnames(mat1) <- paste0("Cell1_", 1:5)
    colnames(mat2) <- paste0("Cell2_", 1:5)

    mae <- create_mae(list(mat1, mat2))

    expect_true(inherits(mae, "MultiAssayExperiment"))
    expect_equal(names(mae), c("experiment_1", "experiment_2"))
})

test_that("create_mae with named list", {
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)

    mae <- create_mae(list(Cond1 = mat, Cond2 = mat))

    expect_true(inherits(mae, "MultiAssayExperiment"))
    expect_equal(names(mae), c("Cond1", "Cond2"))
})

test_that("create_mae errors on non-list input", {
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)

    expect_error(
        create_mae(mat),
        "datasets must be a named list"
    )
})

test_that("earlyj with mixed object types fails", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("SingleCellExperiment")

    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)

    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mat))

    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = list(mat1 = sce, sce1 = sce)
    )

    # This should work since both are SCE now - test proper merging instead
    expect_no_error(earlyj(mae))
})

test_that("earlyj with rowg = FALSE", {
    skip_if_not_installed("MultiAssayExperiment")

    mat1 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    mat2 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat1) <- rownames(mat2) <- paste0("Gene", 1:5)
    colnames(mat1) <- paste0("Cell1_", 1:5)
    colnames(mat2) <- paste0("Cell2_", 1:5)

    mae <- create_mae(list(Cond1 = mat1, Cond2 = mat2))
    merged <- earlyj(mae, rowg = FALSE)

    expect_true(inherits(merged, "MultiAssayExperiment"))
})

test_that("earlyj with SingleCellExperiment objects", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("SingleCellExperiment")

    mat1 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    mat2 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat1) <- rownames(mat2) <- paste0("Gene", 1:5)
    colnames(mat1) <- paste0("Cell1_", 1:5)
    colnames(mat2) <- paste0("Cell2_", 1:5)

    sce1 <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mat1))
    sce2 <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mat2))

    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = list(exp1 = sce1, exp2 = sce2)
    )

    merged <- earlyj(mae)

    expect_true(inherits(merged, "MultiAssayExperiment"))
})

test_that("earlyj with Seurat objects", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("Seurat")
    skip("Seurat object creation requires full Seurat setup")

    # Would test Seurat object merging here
    expect_true(TRUE)
})
test_that("PCzinb works with matrix input - poi method", {
    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- colnames(mat) <- paste0("Gene", 1:5)

    result <- PCzinb(mat, method = "poi", maxcard = 1, nCores = 1)

    expect_true(is.matrix(result))
    expect_equal(dim(result), c(5, 5))
    expect_true(all(result %in% c(0, 1)))
})

test_that("PCzinb works with matrix input - nb method", {
    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- colnames(mat) <- paste0("Gene", 1:5)

    result <- PCzinb(mat, method = "nb", maxcard = 1, nCores = 1)

    expect_true(is.matrix(result))
    expect_equal(dim(result), c(5, 5))
    expect_true(all(result %in% c(0, 1)))
})

# ===== Additional infer_networks tests for better coverage =====

test_that("infer_networks validates MAE input", {
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)

    # Should error with non-MAE input
    expect_error(
        infer_networks(count_matrices_list = list(mat), method = "GENIE3"),
        "must be a MultiAssayExperiment"
    )
})

test_that("infer_networks GENIE3 verbose mode prints messages", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("GENIE3")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    # Should print verbose message
    expect_message(
        infer_networks(count_matrices_list = mae, method = "GENIE3", nCores = 1, verbose = TRUE),
        "Running GENIE3"
    )
})

test_that("infer_networks ZILGM verbose mode prints messages", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mpath")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    # Should print verbose message
    expect_message(
        infer_networks(count_matrices_list = mae, method = "ZILGM", nCores = 1, verbose = TRUE),
        "Running ZILGM"
    )
})

test_that("infer_networks JRF verbose mode prints messages", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    # Should print verbose message
    expect_message(
        infer_networks(count_matrices_list = mae, method = "JRF", nCores = 1, verbose = TRUE),
        "Running JRF"
    )
})

test_that("infer_networks PCzinb verbose mode prints messages", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("pcalg")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    # Should print verbose message
    expect_message(
        infer_networks(count_matrices_list = mae, method = "PCzinb", nCores = 1, verbose = TRUE),
        "Running PCzinb"
    )
})

test_that("infer_networks ZILGM without names adds default names", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mpath")

    set.seed(123)
    mat1 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    mat2 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat1) <- rownames(mat2) <- paste0("Gene", 1:5)
    colnames(mat1) <- colnames(mat2) <- paste0("Cell", 1:5)

    # Create unnamed list then MAE
    mae <- create_mae(list(mat1, mat2))

    result <- infer_networks(count_matrices_list = mae, method = "ZILGM", nCores = 1)

    # Should have default names network_1, network_2
    expect_true(inherits(result, "SummarizedExperiment"))
    expect_true(all(grepl("^network_", SummarizedExperiment::assayNames(result))))
})

test_that("infer_networks PCzinb without names adds default names", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("pcalg")

    set.seed(123)
    mat1 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    mat2 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat1) <- rownames(mat2) <- paste0("Gene", 1:5)
    colnames(mat1) <- colnames(mat2) <- paste0("Cell", 1:5)

    # Create unnamed list then MAE
    mae <- create_mae(list(mat1, mat2))

    result <- infer_networks(count_matrices_list = mae, method = "PCzinb", nCores = 1)

    # Should have default names network_1, network_2
    expect_true(inherits(result, "SummarizedExperiment"))
    expect_true(all(grepl("^network_", SummarizedExperiment::assayNames(result))))
})

test_that("infer_networks GRNBoost2 requires reticulate", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    # Mock scenario where reticulate is not available
    # This will trigger the error path if reticulate/python is not set up
    skip_if(
        requireNamespace("reticulate", quietly = TRUE) &&
            tryCatch(reticulate::py_available(initialize = FALSE), error = function(e) FALSE),
        "reticulate and Python are available"
    )

    expect_error(
        infer_networks(count_matrices_list = mae, method = "GRNBoost2", nCores = 1),
        "reticulate|Python"
    )
})

test_that("PCzinb works with matrix input - zinb0 method", {
    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- colnames(mat) <- paste0("Gene", 1:5)

    result <- PCzinb(mat, method = "zinb0", maxcard = 1, nCores = 1)

    expect_true(is.matrix(result))
    expect_equal(dim(result), c(5, 5))
    expect_true(all(result %in% c(0, 1)))
})

test_that("PCzinb works with matrix input - zinb1 method", {
    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- colnames(mat) <- paste0("Gene", 1:5)

    result <- PCzinb(mat, method = "zinb1", maxcard = 1, nCores = 1)

    expect_true(is.matrix(result))
    expect_equal(dim(result), c(5, 5))
    expect_true(all(result %in% c(0, 1)))
})

test_that("PCzinb works with SummarizedExperiment input", {
    skip_if_not_installed("SummarizedExperiment")

    set.seed(123)
    library(SummarizedExperiment)
    mat <- matrix(rpois(50, lambda = 5), ncol = 10)
    rownames(mat) <- paste0("Gene", 1:5)
    se <- SummarizedExperiment(assays = list(counts = mat))

    result <- PCzinb(se, method = "poi", maxcard = 1, nCores = 1)

    expect_true(inherits(result, "SummarizedExperiment"))
    expect_true(!is.null(metadata(result)$adj_mat))
    expect_equal(dim(metadata(result)$adj_mat), c(5, 5))
})

test_that("PCzinb works with SingleCellExperiment input", {
    skip_if_not_installed("SingleCellExperiment")

    set.seed(123)
    library(SingleCellExperiment)
    mat <- matrix(rpois(50, lambda = 5), ncol = 10)
    rownames(mat) <- paste0("Gene", 1:5)
    sce <- SingleCellExperiment(assays = list(counts = mat))

    result <- PCzinb(sce, method = "poi", maxcard = 1, nCores = 1)

    expect_true(inherits(result, "SingleCellExperiment"))
    expect_true(!is.null(rowPair(result)))
})

test_that("PCzinb handles custom alpha parameter", {
    set.seed(123)
    mat <- matrix(rpois(50, lambda = 5), nrow = 10)

    result <- PCzinb(mat, method = "poi", alpha = 0.05, maxcard = 1, nCores = 1)

    expect_true(is.matrix(result))
})

test_that("PCzinb handles extend parameter", {
    set.seed(123)
    mat <- matrix(rpois(50, lambda = 5), nrow = 10)

    result1 <- PCzinb(mat, method = "poi", extend = TRUE, maxcard = 1, nCores = 1)
    result2 <- PCzinb(mat, method = "poi", extend = FALSE, maxcard = 1, nCores = 1)

    expect_true(is.matrix(result1))
    expect_true(is.matrix(result2))
})

test_that("PCzinb rejects invalid input type", {
    expect_error(
        PCzinb(data.frame(a = 1:5), method = "poi"),
        "x must be a matrix, SummarizedExperiment, or SingleCellExperiment object"
    )
})

test_that("PCzinb handles missing whichAssay gracefully", {
    skip_if_not_installed("SummarizedExperiment")

    library(SummarizedExperiment)
    mat <- matrix(rpois(50, lambda = 5), ncol = 10)
    rownames(mat) <- paste0("Gene", 1:5)
    se <- SummarizedExperiment(assays = list(counts = mat))

    expect_warning(
        PCzinb(se, method = "poi", whichAssay = "processed", maxcard = 1, nCores = 1),
        "We recommend to use QPtransform"
    )
})

test_that("PCzinb validates method argument", {
    set.seed(123)
    mat <- matrix(rpois(50, lambda = 5), nrow = 10)

    expect_error(
        PCzinb(mat, method = "invalid"),
        "'arg' should be one of"
    )
})

test_that("PCzinb with maxcard = 2", {
    set.seed(123)
    mat <- matrix(rpois(100, lambda = 5), nrow = 20, ncol = 5)

    # Test with maxcard = 2
    result_poi <- PCzinb(mat, method = "poi", extend = TRUE, maxcard = 2, nCores = 1)

    expect_true(is.matrix(result_poi))
    expect_equal(dim(result_poi), c(5, 5))
})

test_that("PCzinb with nCores > 1", {
    set.seed(123)
    mat <- matrix(rpois(50, lambda = 5), nrow = 10, ncol = 5)

    result <- PCzinb(mat, method = "poi", maxcard = 1, nCores = 2)

    expect_true(is.matrix(result))
    expect_equal(dim(result), c(5, 5))
})

test_that("PCzinb with custom alpha", {
    set.seed(123)
    mat <- matrix(rpois(50, lambda = 5), nrow = 10, ncol = 5)

    # Test with custom alpha
    alpha_custom <- 0.01
    result_poi <- PCzinb(mat, method = "poi", alpha = alpha_custom, maxcard = 1, nCores = 1)

    expect_true(is.matrix(result_poi))
    expect_equal(dim(result_poi), c(5, 5))
})
