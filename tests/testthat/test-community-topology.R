test_that("compute_topology_metrics calculates network metrics", {
    skip_if_not_installed("igraph")

    mat <- matrix(c(
        0, 1, 1, 0, 0,
        1, 0, 1, 0, 0,
        1, 1, 0, 0, 0,
        0, 0, 0, 0, 1,
        0, 0, 0, 1, 0
    ), nrow = 5, byrow = TRUE)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:5)

    control <- community_path(mat, plot = FALSE, pathway_db = "none")
    predicted <- list(community_path(mat, plot = FALSE, pathway_db = "none"))

    result <- compute_topology_metrics(control, predicted)

    expect_true(is.list(result))
    expect_true("topology_measures" %in% names(result))
    expect_true("control_topology" %in% names(result))
    expect_true(is.data.frame(result$topology_measures))
    expect_true(all(c("Modularity", "Communities", "Density", "Transitivity") %in% colnames(result$topology_measures)))
})

test_that("compute_topology_metrics handles multiple predictions", {
    skip_if_not_installed("igraph")

    mat1 <- matrix(c(
        0, 1, 1, 0,
        1, 0, 1, 0,
        1, 1, 0, 1,
        0, 0, 1, 0
    ), nrow = 4, byrow = TRUE)
    rownames(mat1) <- colnames(mat1) <- paste0("G", 1:4)

    mat2 <- matrix(c(
        0, 1, 0, 0,
        1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0
    ), nrow = 4, byrow = TRUE)
    rownames(mat2) <- colnames(mat2) <- paste0("G", 1:4)

    control <- community_path(mat1, plot = FALSE, pathway_db = "none")
    pred1 <- community_path(mat1, plot = FALSE, pathway_db = "none")
    pred2 <- community_path(mat2, plot = FALSE, pathway_db = "none")

    result <- compute_topology_metrics(control, list(pred1, pred2))

    expect_equal(nrow(result$topology_measures), 2)
    expect_true(all(!is.na(result$control_topology)))
})

test_that("plot_community_comparison creates visualizations", {
    skip_if_not_installed("fmsb")

    control <- list(
        communities = list(membership = c(1, 1, 2, 2, 3)),
        graph = NULL
    )
    predicted <- list(list(
        communities = list(membership = c(1, 1, 2, 2, 3)),
        graph = NULL
    ))

    comm_metrics <- compute_community_metrics(control, predicted)

    # Create minimal topology data
    topology_measures <- data.frame(
        Modularity = 0.5,
        Communities = 3,
        Density = 0.4,
        Transitivity = 0.3
    )
    control_topology <- c(Modularity = 0.5, Communities = 3, Density = 0.4, Transitivity = 0.3)

    # Should not error
    expect_no_error(plot_community_comparison(comm_metrics, topology_measures, control_topology))
})

test_that("community_similarity combines all metrics", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("fmsb")

    mat <- matrix(c(
        0, 1, 1, 0, 0,
        1, 0, 1, 0, 0,
        1, 1, 0, 0, 0,
        0, 0, 0, 0, 1,
        0, 0, 0, 1, 0
    ), nrow = 5, byrow = TRUE)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:5)

    control <- community_path(mat, plot = FALSE, pathway_db = "none")
    predicted <- list(community_path(mat, plot = FALSE, pathway_db = "none"))

    result <- community_similarity(control, predicted, plot = FALSE)

    expect_true(is.list(result))
    expect_true("community_metrics" %in% names(result))
    expect_true("topology_measures" %in% names(result))
    expect_true("control_topology" %in% names(result))
    expect_true(is.data.frame(result$community_metrics))
    expect_true(is.data.frame(result$topology_measures))
})

test_that("community_similarity works with plot = TRUE", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("fmsb")

    mat <- matrix(c(
        0, 1, 1, 0,
        1, 0, 1, 0,
        1, 1, 0, 1,
        0, 0, 1, 0
    ), nrow = 4, byrow = TRUE)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:4)

    control <- community_path(mat, plot = FALSE, pathway_db = "none")
    predicted <- list(community_path(mat, plot = FALSE, pathway_db = "none"))

    # Should create plots without error
    expect_no_error(community_similarity(control, predicted, plot = TRUE))
})

test_that("community_path with verbose mode", {
    skip_if_not_installed("igraph")

    mat <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:3)

    expect_message(
        community_path(mat, plot = FALSE, pathway_db = "none", verbose = TRUE),
        "communities"
    )
})

test_that("community_path with custom method params", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("robin")

    mat <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:3)

    result <- community_path(
        mat,
        methods = "walktrap",
        plot = FALSE,
        pathway_db = "none",
        method_params = list(steps = 3)
    )

    expect_true(is.list(result))
    expect_true("communities" %in% names(result))
})

test_that("community_path with two methods for comparison", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("robin")
    skip("robinCompare parameter compatibility varies by version")

    mat <- matrix(c(
        0, 1, 1, 0, 0, 0,
        1, 0, 1, 0, 0, 0,
        1, 1, 0, 1, 0, 0,
        0, 0, 1, 0, 1, 1,
        0, 0, 0, 1, 0, 1,
        0, 0, 0, 1, 1, 0
    ), nrow = 6, byrow = TRUE)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:6)

    result <- community_path(
        mat,
        methods = c("louvain", "walktrap"),
        plot = FALSE,
        pathway_db = "none",
        verbose = FALSE
    )

    expect_true(is.list(result))
    expect_true("best_method" %in% names(result$communities))
})

test_that("community_path with custom genes_path threshold", {
    skip_if_not_installed("igraph")

    mat <- matrix(c(
        0, 1, 1, 0, 0, 0,
        1, 0, 1, 0, 0, 0,
        1, 1, 0, 1, 0, 0,
        0, 0, 1, 0, 1, 1,
        0, 0, 0, 1, 0, 1,
        0, 0, 0, 1, 1, 0
    ), nrow = 6, byrow = TRUE)
    # Use valid gene symbols to avoid annotation errors
    rownames(mat) <- colnames(mat) <- c("BRCA1", "TP53", "EGFR", "MYC", "KRAS", "AKT1")

    result <- community_path(
        mat,
        plot = FALSE,
        pathway_db = "none",
        genes_path = 2
    )

    expect_true(is.list(result))
})

test_that("community_path with plot = TRUE", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("ggraph")
    skip_if_not_installed("ggplot2")

    mat <- matrix(c(
        0, 1, 1, 0,
        1, 0, 1, 0,
        1, 1, 0, 1,
        0, 0, 1, 0
    ), nrow = 4, byrow = TRUE)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:4)

    expect_no_error(
        community_path(mat, plot = TRUE, pathway_db = "none", verbose = FALSE)
    )
})

test_that("community_path with kegg pathway", {
    skip_if_not_installed("igraph")
    skip("KEGG requires internet and annotation packages")

    mat <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3)
    rownames(mat) <- colnames(mat) <- c("TP53", "BRCA1", "EGFR")

    result <- community_path(
        mat,
        plot = FALSE,
        pathway_db = "kegg",
        genes_path = 1
    )

    expect_true(is.list(result))
})

test_that("community_path with reactome pathway", {
    skip_if_not_installed("igraph")
    skip("Reactome requires internet and annotation packages")

    mat <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3)
    rownames(mat) <- colnames(mat) <- c("TP53", "BRCA1", "EGFR")

    result <- community_path(
        mat,
        plot = FALSE,
        pathway_db = "reactome",
        genes_path = 1
    )

    expect_true(is.list(result))
})

test_that("community_path with custom comparison params", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("robin")
    skip("robinCompare parameter compatibility varies by version")

    mat <- matrix(c(
        0, 1, 1, 0, 0, 0,
        1, 0, 1, 0, 0, 0,
        1, 1, 0, 1, 0, 0,
        0, 0, 1, 0, 1, 1,
        0, 0, 0, 1, 0, 1,
        0, 0, 0, 1, 1, 0
    ), nrow = 6, byrow = TRUE)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:6)

    result <- community_path(
        mat,
        methods = c("louvain", "walktrap"),
        plot = FALSE,
        pathway_db = "none",
        comparison_params = list(measure = "nmi", type = "independent")
    )

    expect_true(is.list(result))
})

test_that("compute_topology_metrics with NULL graph (warning)", {
    skip_if_not_installed("igraph")

    mat <- matrix(c(
        0, 1, 1, 0,
        1, 0, 1, 0,
        1, 1, 0, 1,
        0, 0, 1, 0
    ), nrow = 4, byrow = TRUE)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:4)

    control <- community_path(mat, plot = FALSE, pathway_db = "none")

    # Create predicted with NULL graph
    predicted_null <- list(
        communities = list(membership = c(1, 1, 2, 2)),
        graph = NULL
    )

    expect_warning(
        result <- compute_topology_metrics(control, list(predicted_null)),
        "no valid graph"
    )

    expect_true(all(is.na(result$topology_measures[1, ])))
})

test_that("compute_topology_metrics with single prediction", {
    skip_if_not_installed("igraph")

    mat <- matrix(c(
        0, 1, 1, 0,
        1, 0, 1, 0,
        1, 1, 0, 1,
        0, 0, 1, 0
    ), nrow = 4, byrow = TRUE)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:4)

    control <- community_path(mat, plot = FALSE, pathway_db = "none")
    predicted <- list(community_path(mat, plot = FALSE, pathway_db = "none"))

    result <- compute_topology_metrics(control, predicted)

    expect_equal(nrow(result$topology_measures), 1)
    expect_true(all(!is.na(result$topology_measures[1, ])))
    expect_equal(length(result$control_topology), 4)
})

test_that("compute_topology_metrics with disconnected graph", {
    skip_if_not_installed("igraph")

    # Create disconnected graph
    mat <- matrix(c(
        0, 1, 0, 0,
        1, 0, 0, 0,
        0, 0, 0, 1,
        0, 0, 1, 0
    ), nrow = 4, byrow = TRUE)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:4)

    control <- community_path(mat, plot = FALSE, pathway_db = "none")
    predicted <- list(community_path(mat, plot = FALSE, pathway_db = "none"))

    result <- compute_topology_metrics(control, predicted)

    expect_true(is.data.frame(result$topology_measures))
    expect_true("Density" %in% colnames(result$topology_measures))
})

test_that("community_path with fastgreedy method", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("robin")

    mat <- matrix(c(
        0, 1, 1, 0, 0, 0,
        1, 0, 1, 0, 0, 0,
        1, 1, 0, 1, 0, 0,
        0, 0, 1, 0, 1, 1,
        0, 0, 0, 1, 0, 1,
        0, 0, 0, 1, 1, 0
    ), nrow = 6, byrow = TRUE)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:6)

    result <- community_path(
        mat,
        methods = "fastGreedy", # Note capital G
        plot = FALSE,
        pathway_db = "none"
    )

    expect_true(is.list(result))
    expect_true("communities" %in% names(result))
})

test_that("community_path with spinglass method", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("robin")
    skip("spinglass can be slow and unstable in tests")

    mat <- matrix(c(
        0, 1, 1, 0,
        1, 0, 1, 0,
        1, 1, 0, 1,
        0, 0, 1, 0
    ), nrow = 4, byrow = TRUE)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:4)

    result <- community_path(
        mat,
        methods = "spinglass",
        plot = FALSE,
        pathway_db = "none"
    )

    expect_true(is.list(result))
})

test_that("community_path with infomap method", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("robin")

    mat <- matrix(c(
        0, 1, 1, 0, 0, 0,
        1, 0, 1, 0, 0, 0,
        1, 1, 0, 1, 0, 0,
        0, 0, 1, 0, 1, 1,
        0, 0, 0, 1, 0, 1,
        0, 0, 0, 1, 1, 0
    ), nrow = 6, byrow = TRUE)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:6)

    result <- community_path(
        mat,
        methods = "infomap",
        plot = FALSE,
        pathway_db = "none"
    )

    expect_true(is.list(result))
    expect_true("communities" %in% names(result))
})

test_that("community_path with edge_betweenness method", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("robin")

    mat <- matrix(c(
        0, 1, 1, 0,
        1, 0, 1, 0,
        1, 1, 0, 1,
        0, 0, 1, 0
    ), nrow = 4, byrow = TRUE)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:4)

    result <- community_path(
        mat,
        methods = "edgeBetweenness", # Note camelCase
        plot = FALSE,
        pathway_db = "none"
    )

    expect_true(is.list(result))
    expect_true("communities" %in% names(result))
})
