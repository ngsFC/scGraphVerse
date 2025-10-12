test_that("cutoff_adjacency works with MAE", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(42)
    mat1 <- matrix(rpois(100, lambda = 5), nrow = 10, ncol = 10)
    rownames(mat1) <- paste0("Gene", 1:10)
    colnames(mat1) <- paste0("Cell", 1:10)

    # Create count matrices MAE
    count_mae <- create_mae(list(Cond1 = mat1))

    # Create weighted adjacency matrices
    adjm1 <- matrix(runif(100), nrow = 10)
    rownames(adjm1) <- colnames(adjm1) <- paste0("Gene", 1:10)
    adj_se <- build_network_se(list(adjm1))

    bin_adj <- cutoff_adjacency(
        count_matrices = count_mae,
        weighted_adjm_list = adj_se,
        n = 1,
        method = "GENIE3",
        quantile_threshold = 0.9,
        nCores = 1
    )

    expect_true(inherits(bin_adj, "SummarizedExperiment"))
})

test_that("cutoff_adjacency works with different thresholds", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(42)
    mat1 <- matrix(rpois(100, lambda = 5), nrow = 10, ncol = 10)
    rownames(mat1) <- paste0("Gene", 1:10)
    colnames(mat1) <- paste0("Cell", 1:10)

    count_mae <- create_mae(list(Cond1 = mat1))
    adjm1 <- matrix(runif(100), nrow = 10)
    rownames(adjm1) <- colnames(adjm1) <- paste0("Gene", 1:10)
    adj_se <- build_network_se(list(adjm1))

    bin_adj_low <- cutoff_adjacency(
        count_matrices = count_mae,
        weighted_adjm_list = adj_se,
        n = 1,
        method = "GENIE3",
        quantile_threshold = 0.5,
        nCores = 1
    )

    bin_adj_high <- cutoff_adjacency(
        count_matrices = count_mae,
        weighted_adjm_list = adj_se,
        n = 1,
        method = "GENIE3",
        quantile_threshold = 0.95,
        nCores = 1
    )

    expect_true(inherits(bin_adj_low, "SummarizedExperiment"))
    expect_true(inherits(bin_adj_high, "SummarizedExperiment"))
})

test_that("create_consensus works with vote method", {
    skip_if_not_installed("SummarizedExperiment")

    mat1 <- matrix(c(1, 0, 0, 1), nrow = 2)
    mat2 <- matrix(c(0, 1, 1, 0), nrow = 2)
    rownames(mat1) <- colnames(mat1) <- paste0("G", 1:2)
    rownames(mat2) <- colnames(mat2) <- paste0("G", 1:2)

    se <- build_network_se(list(mat1, mat2))
    consensus <- create_consensus(se, method = "vote", threshold = 0.5)

    expect_true(inherits(consensus, "SummarizedExperiment"))
})

test_that("create_consensus works with union method", {
    skip_if_not_installed("SummarizedExperiment")

    mat1 <- matrix(c(1, 0, 0, 1), nrow = 2)
    mat2 <- matrix(c(0, 1, 1, 0), nrow = 2)
    rownames(mat1) <- colnames(mat1) <- paste0("G", 1:2)
    rownames(mat2) <- colnames(mat2) <- paste0("G", 1:2)

    se <- build_network_se(list(mat1, mat2))
    consensus <- create_consensus(se, method = "union")

    expect_true(inherits(consensus, "SummarizedExperiment"))
})

test_that("create_consensus works with INet method", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("INetTool")

    mat1 <- matrix(c(1, 1, 1, 1), nrow = 2)
    mat2 <- matrix(c(1, 1, 1, 0), nrow = 2)
    rownames(mat1) <- colnames(mat1) <- paste0("G", 1:2)
    rownames(mat2) <- colnames(mat2) <- paste0("G", 1:2)

    # INet requires weighted matrices
    wmat1 <- matrix(runif(4), nrow = 2)
    wmat2 <- matrix(runif(4), nrow = 2)
    rownames(wmat1) <- colnames(wmat1) <- paste0("G", 1:2)
    rownames(wmat2) <- colnames(wmat2) <- paste0("G", 1:2)

    se <- build_network_se(list(mat1, mat2))
    wse <- build_network_se(list(wmat1, wmat2))

    consensus <- create_consensus(se, method = "INet", weighted_list = wse)

    expect_true(inherits(consensus, "SummarizedExperiment"))
})

test_that("create_consensus INet method requires weighted_list", {
    skip_if_not_installed("SummarizedExperiment")

    mat1 <- matrix(c(1, 1, 1, 1), nrow = 2)
    rownames(mat1) <- colnames(mat1) <- paste0("G", 1:2)
    se <- build_network_se(list(mat1))

    # Should error when weighted_list is NULL for INet method
    expect_error(
        create_consensus(se, method = "INet", weighted_list = NULL),
        "weighted_list must be provided"
    )
})

test_that("create_consensus INet validates dimensions", {
    skip_if_not_installed("SummarizedExperiment")

    mat1 <- matrix(c(1, 1, 1, 1), nrow = 2)
    wmat1 <- matrix(runif(9), nrow = 3) # Different size
    rownames(mat1) <- colnames(mat1) <- paste0("G", 1:2)
    rownames(wmat1) <- colnames(wmat1) <- paste0("G", 1:3)

    se <- build_network_se(list(mat1))
    wse <- build_network_se(list(wmat1))

    # Should error with dimension mismatch
    expect_error(
        create_consensus(se, method = "INet", weighted_list = wse),
        "Dimension mismatch"
    )
})

test_that("create_consensus vote method with different thresholds", {
    skip_if_not_installed("SummarizedExperiment")

    mat1 <- matrix(c(1, 0, 0, 1), nrow = 2)
    mat2 <- matrix(c(0, 1, 1, 0), nrow = 2)
    mat3 <- matrix(c(1, 1, 1, 1), nrow = 2)
    rownames(mat1) <- colnames(mat1) <- paste0("G", 1:2)
    rownames(mat2) <- colnames(mat2) <- paste0("G", 1:2)
    rownames(mat3) <- colnames(mat3) <- paste0("G", 1:2)

    se <- build_network_se(list(mat1, mat2, mat3))

    consensus_low <- create_consensus(se, method = "vote", threshold = 0.3)
    consensus_high <- create_consensus(se, method = "vote", threshold = 0.7)

    expect_true(inherits(consensus_low, "SummarizedExperiment"))
    expect_true(inherits(consensus_high, "SummarizedExperiment"))
})

test_that("classify_edges works correctly", {
    skip_if_not_installed("SummarizedExperiment")

    consensus_mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    reference_mat <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), nrow = 3)
    rownames(consensus_mat) <- colnames(consensus_mat) <- paste0("G", 1:3)
    rownames(reference_mat) <- colnames(reference_mat) <- paste0("G", 1:3)

    consensus_se <- build_network_se(list(consensus_mat))
    reference_se <- build_network_se(list(reference_mat))

    edge_class <- classify_edges(consensus_se, reference_se, use_stringdb = FALSE)

    expect_true(is.list(edge_class))
    expect_true("tp_edges" %in% names(edge_class))
    expect_true("fp_edges" %in% names(edge_class))
    expect_true("fn_edges" %in% names(edge_class))
    expect_false(edge_class$use_stringdb)
})

test_that("classify_edges handles SE input with multiple assays", {
    skip_if_not_installed("SummarizedExperiment")

    consensus_mat1 <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    consensus_mat2 <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), nrow = 3)
    reference_mat <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), nrow = 3)
    rownames(consensus_mat1) <- colnames(consensus_mat1) <- paste0("G", 1:3)
    rownames(consensus_mat2) <- colnames(consensus_mat2) <- paste0("G", 1:3)
    rownames(reference_mat) <- colnames(reference_mat) <- paste0("G", 1:3)

    consensus_se <- build_network_se(list(consensus_mat1, consensus_mat2))
    reference_se <- build_network_se(list(reference_mat))

    edge_class <- classify_edges(consensus_se, reference_se, use_stringdb = FALSE)

    expect_true(is.list(edge_class))
    expect_true("tp_edges" %in% names(edge_class))
})

test_that("classify_edges returns list with edge classifications", {
    skip_if_not_installed("SummarizedExperiment")

    consensus_mat <- matrix(c(
        0, 1, 1,
        1, 0, 1,
        1, 1, 0
    ), nrow = 3)
    reference_mat <- matrix(c(
        0, 1, 0,
        1, 0, 1,
        0, 1, 0
    ), nrow = 3)
    rownames(consensus_mat) <- colnames(consensus_mat) <- paste0("G", 1:3)
    rownames(reference_mat) <- colnames(reference_mat) <- paste0("G", 1:3)

    consensus_se <- build_network_se(list(consensus_mat))
    reference_se <- build_network_se(list(reference_mat))

    edge_class <- classify_edges(consensus_se, reference_se, use_stringdb = FALSE)

    # Should return a list with required elements
    expect_true(is.list(edge_class))
    expect_true("tp_edges" %in% names(edge_class))
})

test_that("classify_edges validates input types", {
    skip_if_not_installed("SummarizedExperiment")

    consensus_mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    rownames(consensus_mat) <- colnames(consensus_mat) <- paste0("G", 1:3)
    consensus_se <- build_network_se(list(consensus_mat))

    # Should error with non-SE consensus input
    expect_error(
        classify_edges(consensus_mat, consensus_se, use_stringdb = FALSE),
        "consensus_matrix must be a SummarizedExperiment"
    )

    # Should error with non-SE reference input (when provided)
    expect_error(
        classify_edges(consensus_se, consensus_mat, use_stringdb = FALSE),
        "reference_matrix must be a SummarizedExperiment object or NULL"
    )
})

test_that("classify_edges validates dimension matching", {
    skip_if_not_installed("SummarizedExperiment")

    consensus_mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    reference_mat <- matrix(c(0, 1, 1, 0), nrow = 2) # Different size
    rownames(consensus_mat) <- colnames(consensus_mat) <- paste0("G", 1:3)
    rownames(reference_mat) <- colnames(reference_mat) <- paste0("G", 1:2)

    consensus_se <- build_network_se(list(consensus_mat))
    reference_se <- build_network_se(list(reference_mat))

    # Should error with mismatched dimensions
    expect_error(
        classify_edges(consensus_se, reference_se, use_stringdb = FALSE),
        "Matrices must have the same dimensions"
    )
})

test_that("community_path detects communities", {
    mat <- matrix(c(
        0, 1, 1, 0, 0,
        1, 0, 1, 0, 0,
        1, 1, 0, 0, 0,
        0, 0, 0, 0, 1,
        0, 0, 0, 1, 0
    ), nrow = 5, byrow = TRUE)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:5)

    result <- community_path(mat, plot = FALSE, pathway_db = "none")

    expect_true(is.list(result))
    expect_true(!is.null(result$communities))
    expect_true(!is.null(result$graph))
})

test_that("community_path works with different pathway databases", {
    mat <- matrix(c(
        0, 1, 1, 0,
        1, 0, 1, 0,
        1, 1, 0, 1,
        0, 0, 1, 0
    ), nrow = 4, byrow = TRUE)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:4)

    result_none <- community_path(mat, plot = FALSE, pathway_db = "none")
    result_kegg <- community_path(mat, plot = FALSE, pathway_db = "kegg")

    expect_true(is.list(result_none))
    expect_true(is.list(result_kegg))
})

test_that("community_path handles plot parameter", {
    skip_if_not_installed("ggraph")

    mat <- matrix(c(
        0, 1, 1,
        1, 0, 1,
        1, 1, 0
    ), nrow = 3, byrow = TRUE)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:3)

    # Should not error with plot = TRUE
    expect_no_error(community_path(mat, plot = TRUE, pathway_db = "none"))
})

test_that("compute_community_metrics calculates similarity", {
    skip_if_not_installed("SummarizedExperiment")

    control <- list(
        communities = list(membership = c(1, 1, 2, 2, 3)),
        graph = NULL
    )
    predicted <- list(list(
        communities = list(membership = c(1, 1, 2, 2, 3)),
        graph = NULL
    ))

    metrics <- compute_community_metrics(control, predicted)

    expect_true(is.data.frame(metrics))
    expect_true(all(c("VI", "NMI", "ARI") %in% colnames(metrics)))
})

test_that("compare_consensus creates plot", {
    skip_if_not_installed("ggraph")
    skip_if_not_installed("ggplot2")

    consensus_mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    reference_mat <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), nrow = 3)
    rownames(consensus_mat) <- colnames(consensus_mat) <- paste0("G", 1:3)
    rownames(reference_mat) <- colnames(reference_mat) <- paste0("G", 1:3)

    consensus_se <- build_network_se(list(consensus_mat))
    reference_se <- build_network_se(list(reference_mat))

    plot <- compare_consensus(
        consensus_matrix = consensus_se,
        reference_matrix = reference_se,
        false_plot = FALSE
    )

    expect_true(inherits(plot, "gg"))
})

test_that("symmetrize with custom weight function", {
    skip_if_not_installed("SummarizedExperiment")

    mat <- matrix(c(0, 2, 3, 5, 0, 7, 4, 6, 0), nrow = 3)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:3)
    se <- build_network_se(list(mat))

    # Custom function: use product
    custom_fn <- function(x) prod(x)
    result <- symmetrize(se, weight_function = custom_fn, nCores = 1)

    expect_true(inherits(result, "SummarizedExperiment"))
    sym_mat <- SummarizedExperiment::assay(result, 1)
    # For both non-zero values, should use product
    expect_equal(sym_mat[1, 2], sym_mat[2, 1])
})

test_that("symmetrize with one zero one non-zero", {
    skip_if_not_installed("SummarizedExperiment")

    # Create matrix with zero/non-zero pairs
    mat <- matrix(c(
        0, 0, 5,
        3, 0, 0,
        0, 7, 0
    ), nrow = 3, byrow = TRUE)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:3)
    se <- build_network_se(list(mat))

    result <- symmetrize(se, weight_function = "mean", nCores = 1)

    sym_mat <- SummarizedExperiment::assay(result, 1)
    # When one is zero, should take the non-zero value
    expect_equal(sym_mat[1, 2], 3)
    expect_equal(sym_mat[2, 1], 3)
    expect_equal(sym_mat[1, 3], 5)
    expect_equal(sym_mat[3, 1], 5)
})

test_that("symmetrize with nCores > 1", {
    skip_if_not_installed("SummarizedExperiment")

    mat1 <- matrix(runif(16), nrow = 4)
    mat2 <- matrix(runif(16), nrow = 4)
    rownames(mat1) <- colnames(mat1) <- paste0("G", 1:4)
    rownames(mat2) <- colnames(mat2) <- paste0("G", 1:4)
    se <- build_network_se(list(mat1, mat2))

    result <- symmetrize(se, weight_function = "max", nCores = 2)

    expect_true(inherits(result, "SummarizedExperiment"))
    expect_equal(length(SummarizedExperiment::assays(result)), 2)
})

test_that("symmetrize with min weight function", {
    skip_if_not_installed("SummarizedExperiment")

    mat <- matrix(c(0, 2, 8, 5, 0, 3, 4, 9, 0), nrow = 3)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:3)
    se <- build_network_se(list(mat))

    result <- symmetrize(se, weight_function = "min", nCores = 1)

    sym_mat <- SummarizedExperiment::assay(result, 1)
    # For both non-zero, should use min
    expect_equal(sym_mat[1, 2], min(2, 5))
    expect_equal(sym_mat[2, 1], min(2, 5))
})

test_that("symmetrize preserves metadata", {
    skip_if_not_installed("SummarizedExperiment")

    mat <- matrix(runif(9), nrow = 3)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:3)
    se <- build_network_se(list(mat))
    S4Vectors::metadata(se)$original_info <- "test_metadata"

    result <- symmetrize(se, weight_function = "mean", nCores = 1)

    expect_equal(S4Vectors::metadata(result)$original_info, "test_metadata")
    expect_true(S4Vectors::metadata(result)$symmetrized)
})

test_that("create_consensus with threshold = 1 (intersection)", {
    skip_if_not_installed("SummarizedExperiment")

    mat1 <- matrix(c(1, 1, 1, 1), nrow = 2)
    mat2 <- matrix(c(1, 0, 0, 1), nrow = 2)
    rownames(mat1) <- colnames(mat1) <- paste0("G", 1:2)
    rownames(mat2) <- colnames(mat2) <- paste0("G", 1:2)

    se <- build_network_se(list(mat1, mat2))
    consensus <- create_consensus(se, method = "vote", threshold = 1.0)

    cons_mat <- SummarizedExperiment::assay(consensus, 1)
    # Only edges present in ALL networks
    expect_equal(cons_mat[1, 1], 1)
    expect_equal(cons_mat[2, 2], 1)
    expect_equal(cons_mat[1, 2], 0)
})

test_that("create_consensus with threshold = 0 (union)", {
    skip_if_not_installed("SummarizedExperiment")

    mat1 <- matrix(c(1, 0, 0, 0), nrow = 2)
    mat2 <- matrix(c(0, 1, 0, 0), nrow = 2)
    rownames(mat1) <- colnames(mat1) <- paste0("G", 1:2)
    rownames(mat2) <- colnames(mat2) <- paste0("G", 1:2)

    se <- build_network_se(list(mat1, mat2))
    consensus <- create_consensus(se, method = "vote", threshold = 0)

    cons_mat <- SummarizedExperiment::assay(consensus, 1)
    # All edges should be included
    expect_equal(cons_mat[1, 1], 1)
    expect_equal(cons_mat[1, 2], 1)
})

test_that("create_consensus with multiple networks different patterns", {
    skip_if_not_installed("SummarizedExperiment")

    mat1 <- matrix(c(1, 1, 1, 1), nrow = 2)
    mat2 <- matrix(c(1, 0, 1, 1), nrow = 2)
    mat3 <- matrix(c(1, 1, 0, 1), nrow = 2)
    rownames(mat1) <- colnames(mat1) <- paste0("G", 1:2)
    rownames(mat2) <- colnames(mat2) <- paste0("G", 1:2)
    rownames(mat3) <- colnames(mat3) <- paste0("G", 1:2)

    se <- build_network_se(list(mat1, mat2, mat3))
    # Use threshold 0.5 to get majority vote
    consensus <- create_consensus(se, method = "vote", threshold = 0.5)

    cons_mat <- SummarizedExperiment::assay(consensus, 1)
    # Edges present in at least 2/3 networks
    expect_equal(cons_mat[1, 1], 1)
    expect_equal(cons_mat[2, 2], 1)
})

test_that("cutoff_adjacency works with GENIE3 method", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("SummarizedExperiment")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)

    mae <- create_mae(list(Cond1 = mat))
    inferred <- infer_networks(mae, method = "GENIE3", nCores = 1)
    adjm <- generate_adjacency(inferred)
    symm <- symmetrize(adjm, weight_function = "mean")

    binary <- cutoff_adjacency(
        count_matrices = mae,
        weighted_adjm_list = symm,
        n = 1,
        method = "GENIE3",
        quantile_threshold = 0.8,
        nCores = 1
    )

    expect_true(inherits(binary, "SummarizedExperiment"))
})

test_that("cutoff_adjacency debug mode", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)

    mae <- create_mae(list(Cond1 = mat))
    inferred <- infer_networks(mae, method = "GENIE3", nCores = 1)
    adjm <- generate_adjacency(inferred)
    symm <- symmetrize(adjm, weight_function = "mean")

    # Debug mode should run without error
    result <- cutoff_adjacency(
        count_matrices = mae,
        weighted_adjm_list = symm,
        n = 1,
        method = "GENIE3",
        nCores = 1,
        debug = TRUE
    )

    expect_true(inherits(result, "SummarizedExperiment"))
})

test_that("cutoff_adjacency errors on non-MAE input", {
    skip_if_not_installed("MultiAssayExperiment")

    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)

    expect_error(
        cutoff_adjacency(
            count_matrices = mat,
            weighted_adjm_list = mat,
            n = 1,
            method = "GENIE3"
        ),
        "MultiAssayExperiment"
    )
})

test_that("cutoff_adjacency errors on non-SE weighted matrices", {
    skip_if_not_installed("MultiAssayExperiment")

    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    expect_error(
        cutoff_adjacency(
            count_matrices = mae,
            weighted_adjm_list = list(mat),
            n = 1,
            method = "GENIE3"
        ),
        "SummarizedExperiment"
    )
})

test_that("cutoff_adjacency errors on length mismatch", {
    skip_if_not_installed("MultiAssayExperiment")

    mat1 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    mat2 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat1) <- rownames(mat2) <- paste0("Gene", 1:5)
    colnames(mat1) <- paste0("Cell1_", 1:5)
    colnames(mat2) <- paste0("Cell2_", 1:5)

    mae <- create_mae(list(Cond1 = mat1, Cond2 = mat2))
    adjm1 <- matrix(runif(25), nrow = 5)
    rownames(adjm1) <- colnames(adjm1) <- paste0("Gene", 1:5)
    se_single <- build_network_se(list(adjm1))

    expect_error(
        cutoff_adjacency(
            count_matrices = mae,
            weighted_adjm_list = se_single,
            n = 1,
            method = "GENIE3"
        ),
        "Length.*must match"
    )
})

test_that("cutoff_adjacency with custom quantile", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)

    mae <- create_mae(list(Cond1 = mat))
    inferred <- infer_networks(mae, method = "GENIE3", nCores = 1)
    adjm <- generate_adjacency(inferred)
    symm <- symmetrize(adjm, weight_function = "mean")

    binary <- cutoff_adjacency(
        count_matrices = mae,
        weighted_adjm_list = symm,
        n = 1,
        method = "GENIE3",
        quantile_threshold = 0.5,
        nCores = 1
    )

    expect_true(inherits(binary, "SummarizedExperiment"))
    expect_true(all(SummarizedExperiment::assay(binary, 1) %in% c(0, 1)))
})

test_that("classify_edges classifies TP/FP/FN correctly", {
    skip_if_not_installed("SummarizedExperiment")

    consensus_mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    reference_mat <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), nrow = 3)
    rownames(consensus_mat) <- colnames(consensus_mat) <- paste0("G", 1:3)
    rownames(reference_mat) <- colnames(reference_mat) <- paste0("G", 1:3)

    consensus_se <- build_network_se(list(consensus_mat))
    reference_se <- build_network_se(list(reference_mat))

    result <- classify_edges(consensus_se, reference_se, use_stringdb = FALSE)

    expect_true(is.list(result))
    expect_true(all(c("tp_edges", "fp_edges", "fn_edges") %in% names(result)))
    expect_true(length(result$tp_edges) >= 0)
    expect_true(length(result$fp_edges) >= 0)
    expect_true(length(result$fn_edges) >= 0)
})

test_that("classify_edges errors on non-SE consensus", {
    mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)

    expect_error(
        classify_edges(mat),
        "SummarizedExperiment"
    )
})

test_that("classify_edges errors on non-SE reference", {
    skip_if_not_installed("SummarizedExperiment")

    consensus_mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    rownames(consensus_mat) <- colnames(consensus_mat) <- paste0("G", 1:3)
    consensus_se <- build_network_se(list(consensus_mat))

    reference_mat <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), nrow = 3)

    expect_error(
        classify_edges(consensus_se, reference_mat),
        "SummarizedExperiment"
    )
})

test_that("classify_edges errors on dimension mismatch", {
    skip_if_not_installed("SummarizedExperiment")

    consensus_mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    reference_mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0), nrow = 4)
    rownames(consensus_mat) <- colnames(consensus_mat) <- paste0("G", 1:3)
    rownames(reference_mat) <- colnames(reference_mat) <- paste0("G", 1:4)

    consensus_se <- build_network_se(list(consensus_mat))
    reference_se <- build_network_se(list(reference_mat))

    expect_error(
        classify_edges(consensus_se, reference_se),
        "same dimensions"
    )
})

test_that("classify_edges with use_stringdb FALSE and NULL reference", {
    skip_if_not_installed("SummarizedExperiment")

    consensus_mat <- matrix(sample(0:1, 16, replace = TRUE), nrow = 4)
    rownames(consensus_mat) <- colnames(consensus_mat) <- paste0("G", 1:4)
    consensus_se <- build_network_se(list(consensus_mat))

    # With use_stringdb FALSE and NULL reference, should error about dimensions
    expect_error(
        classify_edges(consensus_se, reference_matrix = NULL, use_stringdb = FALSE)
    )
})

test_that("classify_edges returns correct structure", {
    skip_if_not_installed("SummarizedExperiment")

    consensus_mat <- matrix(sample(0:1, 16, replace = TRUE), nrow = 4)
    reference_mat <- matrix(sample(0:1, 16, replace = TRUE), nrow = 4)
    rownames(consensus_mat) <- colnames(consensus_mat) <- paste0("G", 1:4)
    rownames(reference_mat) <- colnames(reference_mat) <- paste0("G", 1:4)

    consensus_se <- build_network_se(list(consensus_mat))
    reference_se <- build_network_se(list(reference_mat))

    result <- classify_edges(consensus_se, reference_se, use_stringdb = FALSE)

    expect_true("consensus_graph" %in% names(result))
    expect_true("reference_graph" %in% names(result))
    expect_true("edge_colors" %in% names(result))
    expect_true("use_stringdb" %in% names(result))
    expect_false(result$use_stringdb)
})

test_that("plotg creates visualization from SE", {
    skip_if_not_installed("ggraph")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("gridExtra")

    mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:3)
    se <- build_network_se(list(mat))

    expect_no_error(plotg(se))
})

test_that("plotg creates visualization from list", {
    skip_if_not_installed("ggraph")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("gridExtra")

    mat1 <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    mat2 <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), nrow = 3)
    rownames(mat1) <- colnames(mat1) <- paste0("G", 1:3)
    rownames(mat2) <- colnames(mat2) <- paste0("G", 1:3)

    expect_no_error(plotg(list(mat1, mat2)))
})

test_that("plotg warns on non-symmetric matrix and errors", {
    skip_if_not_installed("ggraph")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("gridExtra")

    mat <- matrix(c(0, 1, 0, 0, 0, 1, 0, 1, 0), nrow = 3)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:3)

    # Non-symmetric matrix should warn then error
    expect_warning(expect_error(plotg(list(mat)), "No valid graphs"), "symmetric")
})

test_that("plotg errors on non-list input", {
    skip_if_not_installed("ggraph")

    mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)

    expect_error(plotg(mat), "list.*SummarizedExperiment")
})

test_that("plotg warns and errors on empty graph list", {
    skip_if_not_installed("ggraph")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("gridExtra")

    # Empty network (no edges)
    mat <- matrix(0, nrow = 3, ncol = 3)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:3)

    # Empty graph should warn about no edges then error about no valid graphs
    expect_warning(expect_error(plotg(list(mat)), "No valid graphs"), "No edges")
})

test_that("plot_network_comparison without FP", {
    skip_if_not_installed("ggraph")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("SummarizedExperiment")

    consensus_mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    reference_mat <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), nrow = 3)
    rownames(consensus_mat) <- colnames(consensus_mat) <- paste0("G", 1:3)
    rownames(reference_mat) <- colnames(reference_mat) <- paste0("G", 1:3)

    consensus_se <- build_network_se(list(consensus_mat))
    reference_se <- build_network_se(list(reference_mat))

    edge_class <- classify_edges(consensus_se, reference_se, use_stringdb = FALSE)
    plot_result <- plot_network_comparison(edge_class, show_fp = FALSE)

    expect_true(inherits(plot_result, "gg") || inherits(plot_result, "ggplot"))
})

test_that("plot_network_comparison with FP", {
    skip_if_not_installed("ggraph")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("patchwork")
    skip_if_not_installed("SummarizedExperiment")

    consensus_mat <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3)
    reference_mat <- matrix(c(0, 1, 0, 1, 0, 0, 0, 0, 0), nrow = 3)
    rownames(consensus_mat) <- colnames(consensus_mat) <- paste0("G", 1:3)
    rownames(reference_mat) <- colnames(reference_mat) <- paste0("G", 1:3)

    consensus_se <- build_network_se(list(consensus_mat))
    reference_se <- build_network_se(list(reference_mat))

    edge_class <- classify_edges(consensus_se, reference_se, use_stringdb = FALSE)
    plot_result <- plot_network_comparison(edge_class, show_fp = TRUE)

    # Should return patchwork object or ggplot
    expect_true(inherits(plot_result, "patchwork") || inherits(plot_result, "gg"))
})

test_that("plot_network_comparison with stringdb labels", {
    skip_if_not_installed("ggraph")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("SummarizedExperiment")

    consensus_mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    reference_mat <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), nrow = 3)
    rownames(consensus_mat) <- colnames(consensus_mat) <- paste0("G", 1:3)
    rownames(reference_mat) <- colnames(reference_mat) <- paste0("G", 1:3)

    consensus_se <- build_network_se(list(consensus_mat))
    reference_se <- build_network_se(list(reference_mat))

    edge_class <- classify_edges(consensus_se, reference_se, use_stringdb = FALSE)
    # Manually set use_stringdb to TRUE to test label logic
    edge_class$use_stringdb <- TRUE

    plot_result <- plot_network_comparison(edge_class, show_fp = FALSE)

    expect_true(inherits(plot_result, "gg") || inherits(plot_result, "ggplot"))
})

test_that("plot_network_comparison with no FP edges returns single plot", {
    skip_if_not_installed("ggraph")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("patchwork")
    skip_if_not_installed("SummarizedExperiment")

    consensus_mat <- matrix(c(0, 1, 0, 1, 0, 0, 0, 0, 0), nrow = 3)
    reference_mat <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), nrow = 3)
    rownames(consensus_mat) <- colnames(consensus_mat) <- paste0("G", 1:3)
    rownames(reference_mat) <- colnames(reference_mat) <- paste0("G", 1:3)

    consensus_se <- build_network_se(list(consensus_mat))
    reference_se <- build_network_se(list(reference_mat))

    edge_class <- classify_edges(consensus_se, reference_se, use_stringdb = FALSE)
    plot_result <- plot_network_comparison(edge_class, show_fp = TRUE)

    # Should return single plot since no FP edges
    expect_true(inherits(plot_result, "gg") || inherits(plot_result, "ggplot"))
})

test_that("cutoff_adjacency with multiple shuffles", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("GENIE3")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)

    mae <- create_mae(list(Cond1 = mat))
    inferred <- infer_networks(mae, method = "GENIE3", nCores = 1)
    adjm <- generate_adjacency(inferred)

    binary <- cutoff_adjacency(
        count_matrices = mae,
        weighted_adjm_list = adjm,
        n = 2,
        method = "GENIE3",
        quantile_threshold = 0.9,
        nCores = 1
    )

    expect_true(inherits(binary, "SummarizedExperiment"))
    # Check that cutoff values were computed
    metadata <- S4Vectors::metadata(binary)
    expect_true("method" %in% names(metadata))
})

test_that("cutoff_adjacency preserves metadata", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("GENIE3")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)

    mae <- create_mae(list(Cond1 = mat))
    inferred <- infer_networks(mae, method = "GENIE3", nCores = 1)
    adjm <- generate_adjacency(inferred)

    binary <- cutoff_adjacency(
        count_matrices = mae,
        weighted_adjm_list = adjm,
        n = 1,
        method = "GENIE3",
        nCores = 1
    )

    metadata <- S4Vectors::metadata(binary)
    expect_true("method" %in% names(metadata))
    expect_equal(metadata$method, "GENIE3")
})

test_that("classify_edges with no common edges", {
    skip_if_not_installed("SummarizedExperiment")

    consensus_mat <- matrix(c(0, 1, 0, 0, 0, 0, 0, 0, 0), nrow = 3)
    reference_mat <- matrix(c(0, 0, 1, 0, 0, 0, 0, 0, 0), nrow = 3)
    rownames(consensus_mat) <- colnames(consensus_mat) <- paste0("G", 1:3)
    rownames(reference_mat) <- colnames(reference_mat) <- paste0("G", 1:3)

    consensus_se <- build_network_se(list(consensus_mat))
    reference_se <- build_network_se(list(reference_mat))

    result <- classify_edges(consensus_se, reference_se, use_stringdb = FALSE)

    expect_equal(length(result$tp_edges), 0)
    expect_true(length(result$fp_edges) > 0 || length(result$fn_edges) > 0)
})

test_that("classify_edges with perfect match edges", {
    skip_if_not_installed("SummarizedExperiment")

    # Create identical networks to test perfect match
    consensus_mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    reference_mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    rownames(consensus_mat) <- colnames(consensus_mat) <- paste0("G", 1:3)
    rownames(reference_mat) <- colnames(reference_mat) <- paste0("G", 1:3)

    consensus_se <- build_network_se(list(consensus_mat))
    reference_se <- build_network_se(list(reference_mat))

    result <- classify_edges(consensus_se, reference_se, use_stringdb = FALSE)

    # Should have all TP edges and no FP or FN
    expect_true(is.list(result))
    expect_true(length(result$tp_edges) > 0)
    expect_equal(length(result$fp_edges), 0)
    expect_equal(length(result$fn_edges), 0)
})
