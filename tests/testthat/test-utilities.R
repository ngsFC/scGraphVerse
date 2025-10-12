test_that("selgene selects genes from matrix", {
    set.seed(123)
    mat <- matrix(rpois(1000, lambda = 5), nrow = 100)
    rownames(mat) <- paste0("Gene", 1:100)
    colnames(mat) <- paste0("Cell", 1:10)

    selected <- selgene(
        object = mat,
        top_n = 10,
        remove_rib = FALSE,
        remove_mt = FALSE
    )

    expect_true(is.character(selected))
    expect_equal(length(selected), 10)
})

test_that("selgene removes ribosomal genes", {
    set.seed(123)
    mat <- matrix(rpois(1000, lambda = 5), nrow = 100)
    rownames(mat) <- c(paste0("RPS", 1:10), paste0("RPL", 1:10), paste0("Gene", 1:80))
    colnames(mat) <- paste0("Cell", 1:10)

    selected <- selgene(
        object = mat,
        top_n = 50,
        remove_rib = TRUE,
        remove_mt = FALSE
    )

    expect_true(all(!grepl("^RPS|^RPL", selected)))
})

test_that("selgene removes mitochondrial genes", {
    set.seed(123)
    mat <- matrix(rpois(1000, lambda = 5), nrow = 100)
    rownames(mat) <- c(paste0("MT-", 1:10), paste0("Gene", 1:90))
    colnames(mat) <- paste0("Cell", 1:10)

    selected <- selgene(
        object = mat,
        top_n = 50,
        remove_mt = TRUE,
        remove_rib = FALSE
    )

    expect_true(all(!grepl("^MT-", selected)))
})

test_that("selgene works with SingleCellExperiment", {
    skip_if_not_installed("SingleCellExperiment")

    set.seed(123)
    library(SingleCellExperiment)
    mat <- matrix(rpois(500, lambda = 10), nrow = 50, ncol = 10)
    rownames(mat) <- paste0("Gene", 1:50)
    # Add logcounts with sufficient variance
    logcounts <- log1p(mat) + matrix(rnorm(500, sd = 0.5), nrow = 50)
    sce <- SingleCellExperiment(assays = list(counts = mat, logcounts = logcounts))

    selected <- selgene(
        object = sce,
        top_n = 10,
        remove_rib = FALSE,
        remove_mt = FALSE
    )

    expect_true(is.character(selected))
    expect_true(length(selected) <= 10)
})

test_that("plotg creates visualization", {
    skip_if_not_installed("ggraph")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("gridExtra")

    mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:3)
    se <- build_network_se(list(mat))

    expect_no_error(plotg(se))
})

test_that("pscores calculates precision metrics", {
    skip_if_not_installed("fmsb")
    skip_if_not_installed("SummarizedExperiment")

    pred_mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    true_mat <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), nrow = 3)
    rownames(pred_mat) <- colnames(pred_mat) <- paste0("G", 1:3)
    rownames(true_mat) <- colnames(true_mat) <- paste0("G", 1:3)

    pred_se <- build_network_se(list(pred_mat))
    scores <- pscores(ground_truth = true_mat, predicted_list = pred_se, zero_diag = TRUE)

    expect_true(is.list(scores))
    expect_true("Statistics" %in% names(scores))
    expect_true(is.data.frame(scores$Statistics))
    expect_true("Precision" %in% colnames(scores$Statistics))
})

test_that("plotROC creates ROC plot", {
    skip_if_not_installed("ggplot2")

    pred_mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    true_mat <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), nrow = 3)
    rownames(pred_mat) <- colnames(pred_mat) <- paste0("G", 1:3)
    rownames(true_mat) <- colnames(true_mat) <- paste0("G", 1:3)

    pred_se <- build_network_se(list(pred_mat))
    result <- plotROC(predicted_se = pred_se, ground_truth = true_mat, plot_title = "Test ROC", is_binary = TRUE)

    expect_true(is.list(result))
    expect_true("plot" %in% names(result))
    expect_true("auc" %in% names(result))
    expect_true(inherits(result$plot, "gg"))
})

test_that("plotROC rejects non-SummarizedExperiment input", {
    pred_mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    true_mat <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), nrow = 3)

    # Should error with list input
    expect_error(
        plotROC(predicted_se = list(pred_mat), ground_truth = true_mat, plot_title = "Test", is_binary = TRUE),
        "must be a SummarizedExperiment"
    )

    # Should error with matrix input
    expect_error(
        plotROC(predicted_se = pred_mat, ground_truth = true_mat, plot_title = "Test", is_binary = TRUE),
        "must be a SummarizedExperiment"
    )
})

test_that("selgene validates top_n parameter", {
    mat <- matrix(rpois(100, lambda = 5), nrow = 10)
    rownames(mat) <- paste0("Gene", 1:10)

    expect_error(selgene(object = mat), "valid 'top_n'")
    expect_error(selgene(object = mat, top_n = -1), "valid 'top_n'")
    expect_error(selgene(object = mat, top_n = 0), "valid 'top_n'")
    expect_error(selgene(object = mat, top_n = c(5, 10)), "valid 'top_n'")
    expect_error(selgene(object = mat, top_n = "five"), "valid 'top_n'")
})

test_that("selgene removes both MT and ribosomal genes", {
    mat <- matrix(rpois(1000, lambda = 5), nrow = 100)
    rownames(mat) <- c(
        paste0("MT-", 1:10),
        paste0("RPS", 1:10),
        paste0("RPL", 1:10),
        paste0("Gene", 1:70)
    )

    selected <- selgene(
        object = mat,
        top_n = 50,
        remove_mt = TRUE,
        remove_rib = TRUE
    )

    expect_true(all(!grepl("^MT-", selected)))
    expect_true(all(!grepl("^RPS", selected)))
    expect_true(all(!grepl("^RPL", selected)))
})

test_that("selgene handles top_n larger than available genes", {
    mat <- matrix(rpois(50, lambda = 5), nrow = 5)
    rownames(mat) <- paste0("Gene", 1:5)

    selected <- selgene(object = mat, top_n = 100)

    expect_true(length(selected) <= 5)
})

test_that("selgene case-insensitive gene filtering", {
    # Test that MT and ribosomal filtering works case-insensitively
    mat <- matrix(rpois(400, lambda = 5), nrow = 40)
    rownames(mat) <- c(
        paste0("mt-", 1:5), # lowercase mt
        paste0("MT-", 6:10), # uppercase MT
        paste0("rps", 11:15), # lowercase rps
        paste0("RPS", 16:20), # uppercase RPS
        paste0("Gene", 1:20)
    )

    selected_mt <- selgene(object = mat, top_n = 30, remove_mt = TRUE, remove_rib = FALSE)
    selected_rib <- selgene(object = mat, top_n = 30, remove_mt = FALSE, remove_rib = TRUE)

    # Should remove both lowercase and uppercase versions
    expect_true(all(!grepl("^mt-|^MT-", selected_mt, ignore.case = FALSE)))
    expect_true(all(!grepl("^rps|^RPS|^rpl|^RPL", selected_rib, ignore.case = FALSE)))
})

test_that("selgene with Seurat object and cell_type filtering", {
    skip_if_not_installed("Seurat")
    skip("Seurat objects require actual Seurat setup")

    # This would need actual Seurat object with cell type annotations
    # Skipped to avoid dependency on full Seurat infrastructure
})

test_that("selgene with cell_type parameter on matrix errors", {
    mat <- matrix(rpois(50, lambda = 5), nrow = 5)
    rownames(mat) <- paste0("Gene", 1:5)

    # Cell type filtering not supported for matrices - should error
    expect_error(
        selgene(object = mat, top_n = 3, cell_type = "T_cells"),
        "not supported for raw matrices"
    )
})

test_that("selgene with invalid top_n parameter", {
    mat <- matrix(rpois(50, lambda = 5), nrow = 5)
    rownames(mat) <- paste0("Gene", 1:5)

    expect_error(
        selgene(object = mat, top_n = -5),
        "valid 'top_n'"
    )

    expect_error(
        selgene(object = mat, top_n = c(3, 5)),
        "valid 'top_n'"
    )
})

test_that("selgene with SingleCellExperiment validates assay parameter", {
    skip_if_not_installed("SingleCellExperiment")

    # Create SCE with logcounts which selgene expects
    set.seed(123)
    mat <- matrix(rpois(50, lambda = 10), nrow = 5, ncol = 10)
    rownames(mat) <- paste0("Gene", 1:5)
    # Add logcounts assay which is what selgene expects
    logcounts_mat <- log1p(mat)
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = mat, logcounts = logcounts_mat)
    )

    # Test that selgene works with SCE - it uses logcounts by default
    selected <- selgene(object = sce, top_n = 3)

    expect_true(is.character(selected))
})

test_that(".compute_binary_roc_point calculates TPR and FPR", {
    # Needs binary predictions
    pred <- c(1, 1, 0, 0)
    truth <- c(1, 1, 0, 0)

    result <- scGraphVerse:::.compute_binary_roc_point(pred, truth)

    expect_true(is.list(result))
    expect_true(all(c("TPR", "FPR") %in% names(result)))
})

test_that(".compute_weighted_roc_curve generates ROC curve", {
    skip_if_not_installed("pROC")

    set.seed(123)
    pred <- runif(100)
    truth <- sample(0:1, 100, replace = TRUE)

    result <- scGraphVerse:::.compute_weighted_roc_curve(pred, truth, "test")

    expect_true(is.list(result))
    expect_true(all(c("df", "auc") %in% names(result)))
    expect_true(is.data.frame(result$df))
    expect_true(all(c("TPR", "FPR") %in% colnames(result$df)))
})

test_that(".identify_edges classifies edge types", {
    predicted <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    ground_truth <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), nrow = 3)
    rownames(predicted) <- colnames(predicted) <- paste0("G", 1:3)
    rownames(ground_truth) <- colnames(ground_truth) <- paste0("G", 1:3)

    result <- scGraphVerse:::.identify_edges(predicted, ground_truth, c("TP", "FP", "FN"))

    # Returns a data frame with edge classifications
    expect_true(is.data.frame(result) || is.matrix(result))
})

test_that(".aggregate_cutoffs computes quantile cutoffs", {
    skip("Complex internal function tested via cutoff_adjacency")
    # This function requires specific structure from .run_network_on_shuffled
    # It's already tested indirectly through cutoff_adjacency()
})

test_that(".binarize_adjacency thresholds weighted matrices", {
    set.seed(123)
    weighted <- list(net1 = matrix(runif(9, 0, 1), 3, 3))
    rownames(weighted[[1]]) <- colnames(weighted[[1]]) <- paste0("G", 1:3)
    cutoffs <- c(0.5)
    names(cutoffs) <- "net1"

    binary <- scGraphVerse:::.binarize_adjacency(weighted, cutoffs, "GENIE3", FALSE)

    expect_true(is.list(binary))
    expect_true(all(binary[[1]] %in% c(0, 1)))
})

test_that(".create_adjacency_expansion expands matrix", {
    # Binary adjacency matrix for expansion
    B <- matrix(c(0, 1, 1, 0), nrow = 2)

    result <- scGraphVerse:::.create_adjacency_expansion(B)

    # Returns an expanded adjacency matrix
    expect_true(!is.null(result))
})

test_that(".simulate_counts_ZINB generates count data", {
    skip_if_not_installed("mpath")
    set.seed(123)
    n <- 100
    values <- c(5, 10, 15)
    theta <- 1
    pi <- 0.1

    counts <- scGraphVerse:::.simulate_counts_ZINB(n, values, theta, pi)

    expect_true(is.matrix(counts))
    expect_true(all(counts >= 0))
})

test_that(".add_technical_noise adds noise to data", {
    skip_if_not_installed("mpath")
    set.seed(123)
    n <- 50
    p <- 10
    mu <- 5 # Single value
    pi <- 0.1

    noisy <- scGraphVerse:::.add_technical_noise(n, p, mu, pi)

    expect_true(is.matrix(noisy))
    expect_equal(dim(noisy), c(p, n)) # Returns p x n matrix
})

test_that(".compare_communities calculates similarity metrics", {
    control_comm <- c(1, 1, 2, 2, 3, 3)
    pred_comm <- c(1, 1, 2, 2, 3, 3)

    result <- scGraphVerse:::.compare_communities(control_comm, pred_comm)

    expect_true(is.numeric(result))
    expect_true(length(result) == 3) # VI, NMI, ARI
    expect_true(all(names(result) %in% c("VI", "NMI", "ARI")))
})

test_that(".compute_topo_metrics calculates network topology", {
    skip_if_not_installed("igraph")

    # Create a simple graph
    mat <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3)
    g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected")
    comm <- c(1, 1, 2)

    result <- scGraphVerse:::.compute_topo_metrics(g, comm)

    expect_true(is.numeric(result))
    expect_equal(length(result), 4) # Modularity, Communities, Density, Transitivity
})
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
        depth_range = c(1000, 5000)
    )

    expect_true(is.list(sim_matrices))
    expect_equal(length(sim_matrices), 1)
    expect_equal(dim(sim_matrices[[1]]), c(100, 5))
    expect_true(all(sim_matrices[[1]] >= 0))
})

test_that("zinb_simdata generates multiple matrices", {
    B <- diag(3)
    rownames(B) <- colnames(B) <- paste0("Gene", 1:3)

    sim_matrices <- zinb_simdata(
        n = 50, p = 3, B = B,
        mu_range = list(c(1, 3), c(2, 4)),
        mu_noise = c(0.3, 0.3),
        theta = c(1.5, 1.5),
        pi = c(0.1, 0.1),
        kmat = 2,
        depth_range = c(500, 2000)
    )

    expect_true(is.list(sim_matrices))
    expect_equal(length(sim_matrices), 2)
    expect_equal(dim(sim_matrices[[1]]), c(50, 3))
    expect_equal(dim(sim_matrices[[2]]), c(50, 3))
})
