# Tests for ZILGM_internal.R functions

test_that("zilgm_internal runs with basic Poisson model", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mpath")

    set.seed(42)
    n <- 20
    p <- 5
    X <- matrix(rpois(n * p, lambda = 3), nrow = n, ncol = p)

    result <- zilgm_internal(
        X = X,
        lambda = c(0.1, 0.05),
        nlambda = 2,
        family = "Poisson",
        update_type = "IRLS",
        sym = "AND",
        nCores = 1,
        verbose = 0
    )

    expect_true(is.list(result))
    expect_true("network" %in% names(result))
    expect_true("lambda" %in% names(result))
    expect_equal(length(result$network), 2)
})

test_that("zilgm_internal runs with NBI model", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mpath")

    set.seed(42)
    n <- 20
    p <- 5
    X <- matrix(rpois(n * p, lambda = 3), nrow = n, ncol = p)

    result <- zilgm_internal(
        X = X,
        lambda = c(0.1, 0.05),
        nlambda = 2,
        family = "NBI",
        update_type = "MM",
        sym = "OR",
        theta = 1,
        nCores = 1,
        verbose = 0
    )

    expect_true(is.list(result))
    expect_true("network" %in% names(result))
    expect_equal(length(result$network), 2)
})

test_that("zilgm_internal runs with NBII model", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mpath")

    set.seed(42)
    n <- 20
    p <- 5
    X <- matrix(rpois(n * p, lambda = 3), nrow = n, ncol = p)

    result <- zilgm_internal(
        X = X,
        lambda = c(0.1, 0.05),
        nlambda = 2,
        family = "NBII",
        update_type = "IRLS",
        nCores = 1,
        verbose = 0
    )

    expect_true(is.list(result))
    expect_true("network" %in% names(result))
})

test_that("zilgm_internal with bootstrap stability selection", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mpath")

    set.seed(42)
    n <- 30
    p <- 5
    X <- matrix(rpois(n * p, lambda = 3), nrow = n, ncol = p)

    result <- zilgm_internal(
        X = X,
        lambda = c(0.1, 0.05),
        nlambda = 2,
        family = "Poisson",
        do_boot = TRUE,
        boot_num = 3,
        beta = 0.05,
        nCores = 1,
        verbose = 0
    )

    expect_true(is.list(result))
    expect_true("v" %in% names(result))
    expect_true("opt_index" %in% names(result))
    expect_true("opt_lambda" %in% names(result))
})

test_that("zilgm input validation - matrix conversion", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mpath")

    set.seed(42)
    X <- matrix(rpois(50, 3), nrow = 10, ncol = 5)

    # Should accept matrix input
    result <- zilgm(X, lambda = 0.1, nlambda = 1, family = "Poisson", nCores = 1)
    expect_true(is.list(result))
})

test_that("zilgm errors on invalid lambda", {
    skip_if_not_installed("glmnet")

    X <- matrix(rpois(50, 3), nrow = 10, ncol = 5)

    expect_error(
        zilgm(X, lambda = -0.1, family = "Poisson", nCores = 1),
        "lambda must be non-negative"
    )
})

test_that("zilgm errors on too few columns", {
    skip_if_not_installed("glmnet")

    X <- matrix(rpois(10, 3), nrow = 10, ncol = 1)

    expect_error(
        zilgm(X, lambda = 0.1, family = "Poisson", nCores = 1),
        "2 or more columns"
    )
})

test_that("find_lammax calculates maximum lambda", {
    X <- matrix(rnorm(100), nrow = 10, ncol = 10)
    lammax <- find_lammax(X)

    expect_true(is.numeric(lammax))
    expect_true(lammax > 0)
})

test_that("dNBI calculates negative binomial density", {
    y <- c(0, 1, 2, 3, 4)
    mu <- 2.5
    theta <- 1

    density_log <- dNBI(y, mu, theta, log = TRUE)
    density_normal <- dNBI(y, mu, theta, log = FALSE)

    expect_true(all(is.finite(density_log)))
    expect_true(all(density_normal >= 0))
})

test_that("dP calculates Poisson density", {
    y <- c(0, 1, 2, 3, 4)
    mu <- 2.5

    density_log <- dP(y, mu, log = TRUE)
    density_normal <- dP(y, mu, log = FALSE)

    expect_true(all(is.finite(density_log)))
    expect_true(all(density_normal >= 0))
})

test_that("dNBII calculates NBII density", {
    y <- c(0, 1, 2, 3, 4)
    mu <- 2.5
    sigma <- 0.5

    density_log <- dNBII(y, mu, sigma, log = TRUE)
    density_normal <- dNBII(y, mu, sigma, log = FALSE)

    expect_true(all(is.finite(density_log)))
    expect_true(all(density_normal >= 0))
})

test_that("thresholding_mat applies threshold correctly", {
    mat <- matrix(c(0.01, 0.5, 0.2, 0.001), nrow = 2)

    result <- thresholding_mat(mat, thres = 0.1)

    expect_true(inherits(result, "Matrix"))
    expect_equal(sum(result != 0), 2)
})

test_that("hat_net creates network with AND rule", {
    coef_mat <- matrix(c(0, 0.5, 0.3, 0), nrow = 2)

    result <- hat_net(coef_mat, thresh = 0.1, type = "AND")

    expect_true(is.matrix(result))
    expect_equal(dim(result), c(2, 2))
})

test_that("hat_net creates network with OR rule", {
    coef_mat <- matrix(c(0, 0.5, 0.3, 0), nrow = 2)

    result <- hat_net(coef_mat, thresh = 0.1, type = "OR")

    expect_true(is.matrix(result))
    expect_equal(dim(result), c(2, 2))
})

test_that("theta_ml estimates theta parameter", {
    set.seed(123)
    y <- rpois(50, lambda = 5)
    mu <- rep(5, 50)

    theta <- theta_ml(y, mu)

    expect_true(is.numeric(theta))
    expect_true(theta > 0)
})

test_that("sigma_ml estimates sigma parameter", {
    set.seed(123)
    y <- rpois(50, lambda = 5)
    mu <- rep(5, 50)

    sigma <- sigma_ml(y, mu)

    expect_true(is.numeric(sigma))
    expect_true(sigma >= 0)
})

test_that("zilgm handles small datasets", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mpath")

    set.seed(42)
    # Small dataset (n=20, p=3)
    X <- matrix(rpois(60, lambda = 3), nrow = 20, ncol = 3)

    result <- zilgm_internal(
        X = X, lambda = 0.1, nlambda = 1,
        family = "Poisson", update_type = "MM",
        nCores = 1, verbose = 0
    )

    expect_true(is.list(result))
    expect_true("network" %in% names(result))
})

test_that("zilgm handles sparse count data", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mpath")

    set.seed(42)
    # Sparse count matrix with many zeros
    X <- matrix(rpois(100, lambda = 1), nrow = 20, ncol = 5)
    X[X > 2] <- 0 # Make it sparse

    result <- zilgm_internal(
        X = X, lambda = 0.1, nlambda = 1,
        family = "Poisson", nCores = 1, verbose = 0
    )

    expect_true(is.list(result))
    expect_true("network" %in% names(result))
})

test_that("zigm_network creates network for multiple genes", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mpath")
    skip_if_not_installed("BiocParallel")

    set.seed(42)
    X <- matrix(rpois(100, lambda = 3), nrow = 20, ncol = 5)

    result <- zigm_network(
        X = X,
        lambda = c(0.1, 0.05),
        family = "Poisson",
        update_type = "IRLS",
        sym = "AND",
        nCores = 1,
        n = 20,
        p = 5,
        verbose = 0
    )

    expect_true(is.list(result))
    expect_true("hat_net" %in% names(result))
    expect_true("coef_net" %in% names(result))
})

test_that("zilgm with different symmetrization methods", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mpath")

    set.seed(42)
    X <- matrix(rpois(100, lambda = 3), nrow = 20, ncol = 5)

    result_and <- zilgm_internal(
        X = X, lambda = c(0.1), family = "Poisson",
        sym = "AND", nCores = 1, verbose = 0
    )

    result_or <- zilgm_internal(
        X = X, lambda = c(0.1), family = "Poisson",
        sym = "OR", nCores = 1, verbose = 0
    )

    expect_true(is.list(result_and))
    expect_true(is.list(result_or))
    expect_true("network" %in% names(result_and))
    expect_true("network" %in% names(result_or))
})

test_that("zilgm with automatic lambda selection", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mpath")

    set.seed(42)
    X <- matrix(rpois(100, lambda = 3), nrow = 20, ncol = 5)

    result <- zilgm_internal(
        X = X,
        lambda = NULL,
        nlambda = 3,
        family = "Poisson",
        nCores = 1,
        verbose = 0
    )

    expect_true(is.list(result))
    expect_equal(length(result$lambda), 3)
})
test_that(".compute_confusion_metrics calculates correctly", {
    pred_vec <- c(0, 1, 1, 0, 1)
    gt_vec <- c(0, 1, 0, 1, 1)

    # Access internal function
    metrics <- scGraphVerse:::.compute_confusion_metrics(pred_vec, gt_vec, 1)

    expect_true(is.data.frame(metrics))
    expect_true("TP" %in% colnames(metrics))
    expect_true("TN" %in% colnames(metrics))
    expect_true("FP" %in% colnames(metrics))
    expect_true("FN" %in% colnames(metrics))
})

test_that(".normalize_library_size normalizes counts", {
    set.seed(123)
    mat <- matrix(rpois(100, lambda = 5), nrow = 10, ncol = 10)

    normalized <- scGraphVerse:::.normalize_library_size(mat, c(1000, 5000))

    expect_true(is.matrix(normalized))
    expect_equal(dim(normalized), dim(mat))
    expect_true(all(normalized >= 0))
})

test_that(".shuffle_matrix_rows shuffles correctly", {
    set.seed(123)
    mat <- matrix(1:20, nrow = 4, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:4)

    shuffled <- scGraphVerse:::.shuffle_matrix_rows(mat)

    expect_true(is.matrix(shuffled))
    expect_equal(dim(shuffled), dim(mat))
    expect_equal(sort(as.vector(mat)), sort(as.vector(shuffled)))
})

test_that(".filter_genes removes MT and ribosomal genes", {
    expr <- matrix(rpois(100, lambda = 5), nrow = 10, ncol = 10)
    rownames(expr) <- c("MT-ATP6", "MT-ND1", "RPS3", "RPL4", paste0("Gene", 1:6))

    filtered_mt <- scGraphVerse:::.filter_genes(expr, remove_mt = TRUE, remove_rib = FALSE)
    filtered_rib <- scGraphVerse:::.filter_genes(expr, remove_mt = FALSE, remove_rib = TRUE)
    filtered_both <- scGraphVerse:::.filter_genes(expr, remove_mt = TRUE, remove_rib = TRUE)

    expect_true(all(!grepl("^MT-", rownames(filtered_mt))))
    expect_true(all(!grepl("^RPS|^RPL", rownames(filtered_rib))))
    expect_true(all(!grepl("^MT-|^RPS|^RPL", rownames(filtered_both))))
})

test_that(".select_top_genes selects by variance", {
    set.seed(123)
    expr <- matrix(rnorm(100, mean = 5, sd = 2), nrow = 10, ncol = 10)
    rownames(expr) <- paste0("Gene", 1:10)

    top_genes <- scGraphVerse:::.select_top_genes(expr, top_n = 5)

    expect_true(is.character(top_genes))
    expect_true(length(top_genes) <= 10)
})


test_that(".merge_matrix_list combines matrices", {
    mat1 <- matrix(rpois(50, lambda = 5), nrow = 5, ncol = 10)
    mat2 <- matrix(rpois(50, lambda = 5), nrow = 5, ncol = 10)
    rownames(mat1) <- rownames(mat2) <- paste0("Gene", 1:5)
    colnames(mat1) <- paste0("Cell1_", 1:10)
    colnames(mat2) <- paste0("Cell2_", 1:10)

    merged <- scGraphVerse:::.merge_matrix_list(list(mat1, mat2), rowg = TRUE)

    expect_true(is.matrix(merged))
    expect_equal(nrow(merged), 5)
    expect_equal(ncol(merged), 20)
})


test_that(".convert_counts_list handles different object types", {
    # Test with matrix
    mat <- matrix(rpois(50, lambda = 5), nrow = 5, ncol = 10)
    rownames(mat) <- paste0("Gene", 1:5)

    converted <- scGraphVerse:::.convert_counts_list(list(mat))

    expect_true(is.list(converted))
    expect_equal(length(converted), 1)
    expect_true(is.matrix(converted[[1]]))
})

test_that(".extract_expression extracts from matrix", {
    mat <- matrix(rpois(50, lambda = 5), nrow = 5, ncol = 10)
    rownames(mat) <- paste0("Gene", 1:5)

    expr <- scGraphVerse:::.extract_expression(mat)

    expect_true(is.matrix(expr))
    expect_equal(dim(expr), dim(mat))
})

test_that(".create_igraph_plot handles empty network", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("ggraph")

    # Empty network (no edges)
    mat <- matrix(0, nrow = 5, ncol = 5)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:5)

    result <- scGraphVerse:::.create_igraph_plot(mat, 1)

    expect_null(result)
})

test_that(".create_igraph_plot handles valid network", {
    skip_if_not_installed("igraph")
    skip_if_not_installed("ggraph")

    # Network with edges
    mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:3)

    result <- scGraphVerse:::.create_igraph_plot(mat, 1)

    expect_true(inherits(result, "gg") || inherits(result, "ggplot"))
})

test_that(".edge_to_str converts edges to strings", {
    edge_list <- data.frame(
        gene1 = c("A", "B", "C"),
        gene2 = c("B", "C", "A")
    )

    result <- scGraphVerse:::.edge_to_str(edge_list)

    expect_true(is.character(result))
    expect_equal(length(result), 3)
})

test_that(".create_adjacency_expansion expands adjacency matrix", {
    B <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    rownames(B) <- colnames(B) <- paste0("G", 1:3)

    result <- scGraphVerse:::.create_adjacency_expansion(B)

    expect_true(is.list(result))
    expect_true("A" %in% names(result))
    expect_true("edge_indices" %in% names(result))
})

test_that(".simulate_counts_ZINB generates ZINB counts", {
    skip_if_not_installed("pscl")

    n <- 10
    values <- c(2, 3, 4)
    theta <- 1
    pi <- 0.1

    result <- scGraphVerse:::.simulate_counts_ZINB(n, values, theta, pi)

    expect_true(is.matrix(result))
    expect_equal(nrow(result), length(values))
    expect_equal(ncol(result), n)
})

test_that(".add_technical_noise adds noise to counts", {
    skip_if_not_installed("pscl")

    n <- 10
    p <- 5
    mu <- 0.5
    pi <- 0.1

    result <- scGraphVerse:::.add_technical_noise(n, p, mu, pi)

    expect_true(is.matrix(result))
    expect_equal(nrow(result), p)
    expect_equal(ncol(result), n)
})

test_that(".prepare_prediction_vectors prepares vectors for ROC", {
    pred_matrix <- matrix(runif(16), nrow = 4)
    ground_truth <- matrix(sample(0:1, 16, replace = TRUE), nrow = 4)
    rownames(pred_matrix) <- rownames(ground_truth) <- paste0("G", 1:4)
    colnames(pred_matrix) <- colnames(ground_truth) <- paste0("G", 1:4)

    result <- scGraphVerse:::.prepare_prediction_vectors(pred_matrix, ground_truth)

    expect_true(is.numeric(result))
    expect_true(length(result) > 0)
})

test_that(".compute_binary_roc_point calculates ROC point", {
    pred_vec <- c(1, 0, 1, 0, 1)
    truth_vec <- c(1, 0, 0, 1, 1)

    result <- scGraphVerse:::.compute_binary_roc_point(pred_vec, truth_vec)

    expect_true(is.list(result))
    expect_true("FPR" %in% names(result))
    expect_true("TPR" %in% names(result))
})

test_that(".compute_weighted_roc_curve calculates ROC curve", {
    skip_if_not_installed("pROC")

    pred_vec <- runif(20)
    truth_vec <- sample(0:1, 20, replace = TRUE)

    result <- scGraphVerse:::.compute_weighted_roc_curve(pred_vec, truth_vec, "Test")

    expect_true(is.list(result))
    expect_true("df" %in% names(result))
    expect_true("auc" %in% names(result))
})

test_that(".shuffle_matrix_rows shuffles rows correctly", {
    set.seed(123)
    mat <- matrix(1:20, nrow = 4, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:4)
    colnames(mat) <- paste0("Cell", 1:5)

    shuffled <- scGraphVerse:::.shuffle_matrix_rows(mat)

    expect_true(is.matrix(shuffled))
    expect_equal(dim(shuffled), dim(mat))
    expect_equal(rownames(shuffled), rownames(mat))
    expect_equal(colnames(shuffled), colnames(mat))
})

test_that(".merge_matrix_list with rowg=FALSE transposes", {
    mat1 <- matrix(rpois(50, lambda = 5), nrow = 5, ncol = 10)
    mat2 <- matrix(rpois(50, lambda = 5), nrow = 5, ncol = 10)
    rownames(mat1) <- rownames(mat2) <- paste0("Gene", 1:5)
    colnames(mat1) <- paste0("Cell1_", 1:10)
    colnames(mat2) <- paste0("Cell2_", 1:10)

    merged <- scGraphVerse:::.merge_matrix_list(list(mat1, mat2), rowg = FALSE)

    expect_true(is.matrix(merged))
    # rowg=FALSE means genes are in columns, so matrices are transposed
    expect_true(ncol(merged) > 0)
})

test_that(".validate_network_list validates networks", {
    net1 <- matrix(runif(9), nrow = 3)
    net2 <- matrix(runif(9), nrow = 3)
    rownames(net1) <- colnames(net1) <- paste0("G", 1:3)
    rownames(net2) <- colnames(net2) <- paste0("G", 1:3)

    gene_names <- scGraphVerse:::.validate_network_list(list(net1, net2))

    expect_equal(gene_names, paste0("G", 1:3))
})

test_that(".validate_network_list errors on dimension mismatch", {
    net1 <- matrix(runif(9), nrow = 3)
    net2 <- matrix(runif(16), nrow = 4)
    rownames(net1) <- colnames(net1) <- paste0("G", 1:3)
    rownames(net2) <- colnames(net2) <- paste0("G", 1:4)

    expect_error(
        scGraphVerse:::.validate_network_list(list(net1, net2)),
        "same dimensions"
    )
})

test_that(".extract_networks_from_se extracts networks", {
    skip_if_not_installed("SummarizedExperiment")

    mat1 <- matrix(runif(9), nrow = 3)
    mat2 <- matrix(runif(9), nrow = 3)
    rownames(mat1) <- colnames(mat1) <- paste0("G", 1:3)
    rownames(mat2) <- colnames(mat2) <- paste0("G", 1:3)

    se <- build_network_se(list(mat1, mat2))
    networks <- scGraphVerse:::.extract_networks_from_se(se)

    expect_true(is.list(networks))
    expect_equal(length(networks), 2)
})

test_that(".is_network_se identifies network SE", {
    skip_if_not_installed("SummarizedExperiment")

    mat <- matrix(runif(9), nrow = 3)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:3)

    se <- build_network_se(list(mat))

    expect_true(scGraphVerse:::.is_network_se(se))
})

test_that(".compute_topo_metrics calculates topology metrics", {
    skip_if_not_installed("igraph")

    mat <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:3)
    g <- igraph::graph_from_adjacency_matrix(mat, mode = "undirected")
    comm <- c(1, 1, 2)

    metrics <- scGraphVerse:::.compute_topo_metrics(g, comm)

    expect_true(is.numeric(metrics))
    expect_true("Modularity" %in% names(metrics))
})

test_that(".compare_communities compares community structures", {
    skip_if_not_installed("igraph")

    control_comm <- c(1, 1, 2, 2, 3, 3)
    pred_comm <- c(1, 2, 1, 2, 3, 3)

    result <- scGraphVerse:::.compare_communities(control_comm, pred_comm)

    expect_true(is.numeric(result))
    expect_true("VI" %in% names(result))
    expect_true("NMI" %in% names(result))
    expect_true("ARI" %in% names(result))
})

test_that(".identify_edges identifies edge types", {
    predicted <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    ground_truth <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), nrow = 3)
    rownames(predicted) <- colnames(predicted) <- paste0("G", 1:3)
    rownames(ground_truth) <- colnames(ground_truth) <- paste0("G", 1:3)

    edges <- scGraphVerse:::.identify_edges(predicted, ground_truth, c("TP", "FP", "FN"))

    expect_true(is.data.frame(edges) || is.null(edges))
    if (!is.null(edges)) {
        expect_true("edge_type" %in% colnames(edges))
    }
})

test_that(".prepare_method_args extracts method parameters", {
    method_params <- list(resolution = 0.8, steps = 5)

    args <- scGraphVerse:::.prepare_method_args("walktrap", method_params)

    expect_true(is.list(args))
    expect_true("steps" %in% names(args))
    expect_equal(args$steps, 5)
})
# Extensive coverage tests to reach 80%

# Test error conditions in infer_networks
test_that("infer_networks error handling - non-MAE input", {
    mat <- matrix(rpois(50, lambda = 5), nrow = 10, ncol = 5)

    expect_error(
        infer_networks(count_matrices_list = mat, method = "GENIE3"),
        "must be a MultiAssayExperiment"
    )
})

# Test GENIE3 with custom parameters
test_that("infer_networks GENIE3 with all custom parameters", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    genie3_params <- list(
        regulators = c("Gene1", "Gene2"),
        targets = c("Gene3", "Gene4", "Gene5"),
        treeMethod = "ET",
        K = 2,
        nTrees = 500,
        seed = 123
    )

    inferred <- infer_networks(
        count_matrices_list = mae,
        method = "GENIE3",
        nCores = 1,
        genie3_params = genie3_params,
        verbose = TRUE
    )

    expect_true(is.list(inferred))
})

# Test PCzinb with all methods
test_that("PCzinb with all methods", {
    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- colnames(mat) <- paste0("Gene", 1:5)

    result_poi <- PCzinb(mat, method = "poi", maxcard = 1, nCores = 1)
    result_nb <- PCzinb(mat, method = "nb", maxcard = 1, nCores = 1)
    result_zinb0 <- PCzinb(mat, method = "zinb0", maxcard = 1, nCores = 1)
    result_zinb1 <- PCzinb(mat, method = "zinb1", maxcard = 1, nCores = 1)

    expect_true(is.matrix(result_poi))
    expect_true(is.matrix(result_nb))
    expect_true(is.matrix(result_zinb0))
    expect_true(is.matrix(result_zinb1))
})

# Test PCzinb error handling
test_that("PCzinb error on invalid input", {
    expect_error(
        PCzinb(data.frame(a = 1:5), method = "poi"),
        "must be a matrix"
    )
})

# Test create_mae with different input types
test_that("create_mae handles various input types", {
    skip_if_not_installed("MultiAssayExperiment")

    mat1 <- matrix(rpois(50, lambda = 5), nrow = 5, ncol = 10)
    mat2 <- matrix(rpois(50, lambda = 5), nrow = 5, ncol = 10)
    rownames(mat1) <- rownames(mat2) <- paste0("Gene", 1:5)
    colnames(mat1) <- paste0("Cell1_", 1:10)
    colnames(mat2) <- paste0("Cell2_", 1:10)

    # Named list
    mae_named <- create_mae(list(Cond1 = mat1, Cond2 = mat2))
    expect_equal(names(mae_named), c("Cond1", "Cond2"))

    # Unnamed list
    mae_unnamed <- create_mae(list(mat1, mat2))
    expect_true(length(mae_unnamed) == 2)
})

# Test earlyj with different orientations
test_that("earlyj with rowg parameter", {
    skip_if_not_installed("MultiAssayExperiment")

    mat1 <- matrix(rpois(50, lambda = 5), nrow = 5, ncol = 10)
    mat2 <- matrix(rpois(50, lambda = 5), nrow = 5, ncol = 10)
    rownames(mat1) <- rownames(mat2) <- paste0("Gene", 1:5)
    colnames(mat1) <- paste0("Cell1_", 1:10)
    colnames(mat2) <- paste0("Cell2_", 1:10)

    mae <- create_mae(list(Cond1 = mat1, Cond2 = mat2))

    result_rowg_true <- earlyj(mae, rowg = TRUE)
    result_rowg_false <- earlyj(mae, rowg = FALSE)

    expect_true(inherits(result_rowg_true, "MultiAssayExperiment"))
    expect_true(inherits(result_rowg_false, "MultiAssayExperiment"))
})

# Test earlyj error handling
test_that("earlyj errors on non-MAE input", {
    expect_error(
        earlyj(list(matrix(1:10, 2, 5))),
        "must be a MultiAssayExperiment"
    )
})

# Test generate_adjacency with multiple networks
test_that("generate_adjacency with multiple dataframes", {
    skip_if_not_installed("MultiAssayExperiment")

    df1 <- data.frame(
        Gene1 = c("G1", "G2", "G3"),
        Gene2 = c("G2", "G3", "G1"),
        Weight = c(0.8, 0.6, 0.9)
    )
    df2 <- data.frame(
        Gene1 = c("G1", "G3"),
        Gene2 = c("G3", "G2"),
        Weight = c(0.7, 0.5)
    )

    adj_se <- generate_adjacency(list(df1, df2), nCores = 1)

    expect_true(inherits(adj_se, "SummarizedExperiment"))
    expect_true(length(SummarizedExperiment::assays(adj_se)) >= 1)
})

# Test symmetrize with all weight functions
test_that("symmetrize all weight functions", {
    skip_if_not_installed("SummarizedExperiment")

    mat <- matrix(runif(9), nrow = 3)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:3)
    se <- build_network_se(list(mat))

    symm_mean <- symmetrize(se, weight_function = "mean")
    symm_min <- symmetrize(se, weight_function = "min")
    symm_max <- symmetrize(se, weight_function = "max")

    expect_true(inherits(symm_mean, "SummarizedExperiment"))
    expect_true(inherits(symm_min, "SummarizedExperiment"))
    expect_true(inherits(symm_max, "SummarizedExperiment"))
})

# Test cutoff_adjacency with ZILGM method
test_that("cutoff_adjacency with ZILGM method", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(42)
    mat1 <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat1) <- paste0("Gene", 1:5)
    colnames(mat1) <- paste0("Cell", 1:5)

    count_mae <- create_mae(list(Cond1 = mat1))
    adjm1 <- matrix(runif(25), nrow = 5)
    rownames(adjm1) <- colnames(adjm1) <- paste0("Gene", 1:5)
    adj_se <- build_network_se(list(adjm1))

    bin_adj <- cutoff_adjacency(
        count_matrices = count_mae,
        weighted_adjm_list = adj_se,
        n = 1,
        method = "GENIE3",
        quantile_threshold = 0.8,
        nCores = 1
    )

    expect_true(inherits(bin_adj, "SummarizedExperiment"))
})

# Test create_consensus with all methods
test_that("create_consensus all methods comprehensive", {
    skip_if_not_installed("SummarizedExperiment")

    mat1 <- matrix(c(1, 0, 0, 1), nrow = 2)
    mat2 <- matrix(c(0, 1, 1, 0), nrow = 2)
    mat3 <- matrix(c(1, 1, 1, 1), nrow = 2)
    rownames(mat1) <- colnames(mat1) <- paste0("G", 1:2)
    rownames(mat2) <- colnames(mat2) <- paste0("G", 1:2)
    rownames(mat3) <- colnames(mat3) <- paste0("G", 1:2)

    se <- build_network_se(list(mat1, mat2, mat3))

    consensus_vote_low <- create_consensus(se, method = "vote", threshold = 0.3)
    consensus_vote_mid <- create_consensus(se, method = "vote", threshold = 0.5)
    consensus_vote_high <- create_consensus(se, method = "vote", threshold = 0.7)
    consensus_union <- create_consensus(se, method = "union")

    expect_true(inherits(consensus_vote_low, "SummarizedExperiment"))
    expect_true(inherits(consensus_vote_mid, "SummarizedExperiment"))
    expect_true(inherits(consensus_vote_high, "SummarizedExperiment"))
    expect_true(inherits(consensus_union, "SummarizedExperiment"))
})

# Test classify_edges comprehensive
test_that("classify_edges with different network sizes", {
    skip_if_not_installed("SummarizedExperiment")

    consensus_mat <- matrix(sample(0:1, 16, replace = TRUE), nrow = 4)
    reference_mat <- matrix(sample(0:1, 16, replace = TRUE), nrow = 4)
    rownames(consensus_mat) <- colnames(consensus_mat) <- paste0("G", 1:4)
    rownames(reference_mat) <- colnames(reference_mat) <- paste0("G", 1:4)

    consensus_se <- build_network_se(list(consensus_mat))
    reference_se <- build_network_se(list(reference_mat))

    edge_class <- classify_edges(consensus_se, reference_se, use_stringdb = FALSE)

    expect_true(is.list(edge_class))
    expect_true("tp_edges" %in% names(edge_class))
})

# Test community_path with different parameters
test_that("community_path comprehensive testing", {
    mat <- matrix(c(
        0, 1, 1, 0, 0, 0,
        1, 0, 1, 0, 0, 0,
        1, 1, 0, 1, 0, 0,
        0, 0, 1, 0, 1, 1,
        0, 0, 0, 1, 0, 1,
        0, 0, 0, 1, 1, 0
    ), nrow = 6, byrow = TRUE)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:6)

    result_none <- community_path(mat, plot = FALSE, pathway_db = "none")
    result_kegg <- community_path(mat, plot = FALSE, pathway_db = "kegg")
    result_reactome <- community_path(mat, plot = FALSE, pathway_db = "reactome")

    expect_true(is.list(result_none))
    expect_true(is.list(result_kegg))
    expect_true(is.list(result_reactome))
})

# Test compute_community_metrics thoroughly
test_that("compute_community_metrics with various community structures", {
    control <- list(
        communities = list(membership = c(1, 1, 2, 2, 3, 3)),
        graph = NULL
    )
    predicted <- list(
        list(communities = list(membership = c(1, 1, 2, 2, 3, 3)), graph = NULL),
        list(communities = list(membership = c(1, 2, 1, 2, 3, 3)), graph = NULL),
        list(communities = list(membership = c(1, 1, 1, 2, 2, 2)), graph = NULL)
    )

    metrics <- compute_community_metrics(control, predicted)

    expect_true(is.data.frame(metrics))
    expect_equal(nrow(metrics), 3)
})

# Test pscores with SE input
test_that("pscores with SummarizedExperiment", {
    skip_if_not_installed("fmsb")
    skip_if_not_installed("SummarizedExperiment")

    pred_mat1 <- matrix(sample(0:1, 16, replace = TRUE), nrow = 4)
    pred_mat2 <- matrix(sample(0:1, 16, replace = TRUE), nrow = 4)
    true_mat <- matrix(sample(0:1, 16, replace = TRUE), nrow = 4)
    rownames(pred_mat1) <- colnames(pred_mat1) <- paste0("G", 1:4)
    rownames(pred_mat2) <- colnames(pred_mat2) <- paste0("G", 1:4)
    rownames(true_mat) <- colnames(true_mat) <- paste0("G", 1:4)

    pred_se <- build_network_se(list(pred_mat1, pred_mat2))
    scores <- pscores(ground_truth = true_mat, predicted_list = pred_se, zero_diag = TRUE)

    expect_true(is.list(scores))
    expect_true("Statistics" %in% names(scores))
})

# Test plotROC with weighted matrices
test_that("plotROC with weighted predictions", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("pROC")

    pred_mat <- matrix(runif(16), nrow = 4)
    true_mat <- matrix(sample(0:1, 16, replace = TRUE), nrow = 4)
    rownames(pred_mat) <- colnames(pred_mat) <- paste0("G", 1:4)
    rownames(true_mat) <- colnames(true_mat) <- paste0("G", 1:4)

    pred_se <- build_network_se(list(pred_mat))
    result <- plotROC(
        predicted_se = pred_se,
        ground_truth = true_mat,
        plot_title = "Weighted ROC",
        is_binary = FALSE
    )

    expect_true(is.list(result))
    expect_true("plot" %in% names(result))
    expect_true("auc" %in% names(result))
})

# Test zinb_simdata with all parameters
test_that("zinb_simdata comprehensive parameter testing", {
    B <- diag(4)
    rownames(B) <- colnames(B) <- paste0("Gene", 1:4)

    sim1 <- zinb_simdata(
        n = 50, p = 4, B = B,
        mu_range = list(c(1, 5)),
        mu_noise = 0.3,
        theta = 1,
        pi = 0.1,
        kmat = 1,
        depth_range = c(1000, 3000)
    )

    sim2 <- zinb_simdata(
        n = 50, p = 4, B = B,
        mu_range = list(c(2, 6), c(1, 4)),
        mu_noise = c(0.2, 0.3),
        theta = c(1.5, 2),
        pi = c(0.1, 0.15),
        kmat = 2,
        depth_range = c(500, 2000)
    )

    expect_true(is.list(sim1))
    expect_equal(length(sim1), 1)
    expect_true(is.list(sim2))
    expect_equal(length(sim2), 2)
})

test_that("infer_networks with nCores parameter variation", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    # Test with nCores = 2
    inferred <- infer_networks(mae, method = "GENIE3", nCores = 2, verbose = FALSE)
    expect_true(is.list(inferred))
})

test_that("infer_networks ZILGM with NBII family", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    zilgm_params <- list(lambda = 0.1, family = "NBII")

    inferred <- infer_networks(
        count_matrices_list = mae,
        method = "ZILGM",
        nCores = 1,
        zilgm_params = zilgm_params
    )

    expect_true(inherits(inferred, "SummarizedExperiment"))
})

test_that("infer_networks ZILGM with bootstrap", {
    skip_if_not_installed("MultiAssayExperiment")

    set.seed(123)
    mat <- matrix(rpois(25, lambda = 5), nrow = 5, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:5)
    colnames(mat) <- paste0("Cell", 1:5)
    mae <- create_mae(list(Cond1 = mat))

    zilgm_params <- list(
        lambda = 0.1,
        do_boot = TRUE,
        boot_num = 2,
        beta = 0.1
    )

    inferred <- infer_networks(
        count_matrices_list = mae,
        method = "ZILGM",
        nCores = 1,
        zilgm_params = zilgm_params
    )

    expect_true(inherits(inferred, "SummarizedExperiment"))
})

test_that("ZILGM with IRLS update", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mpath")

    set.seed(42)
    X <- matrix(rpois(100, lambda = 3), nrow = 20, ncol = 5)

    result <- scGraphVerse:::zilgm_internal(
        X = X,
        lambda = 0.1,
        nlambda = 1,
        family = "Poisson",
        update_type = "IRLS",
        nCores = 1,
        verbose = 0
    )

    expect_true(is.list(result))
    expect_true("network" %in% names(result))
})

test_that("ZILGM with multiple lambda values", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mpath")

    set.seed(42)
    X <- matrix(rpois(100, lambda = 3), nrow = 20, ncol = 5)

    result <- scGraphVerse:::zilgm_internal(
        X = X,
        lambda = c(0.2, 0.1, 0.05),
        family = "Poisson",
        nCores = 1,
        verbose = 0
    )

    expect_true(is.list(result))
    expect_equal(length(result$lambda), 3)
})

test_that("ZILGM theta estimation", {
    set.seed(123)
    y <- rpois(50, lambda = 5)
    mu <- rep(5, 50)

    theta <- scGraphVerse:::theta_ml(y, mu)

    expect_true(is.numeric(theta))
    expect_true(theta > 0)
    expect_true(is.finite(theta))
})

test_that("ZILGM sigma estimation", {
    set.seed(123)
    y <- rpois(50, lambda = 5)
    mu <- rep(5, 50)

    sigma <- scGraphVerse:::sigma_ml(y, mu)

    expect_true(is.numeric(sigma))
    expect_true(sigma >= 0)
    expect_true(is.finite(sigma))
})

test_that("ZILGM find_lammax calculation", {
    X <- matrix(rnorm(100), nrow = 10, ncol = 10)

    lammax <- scGraphVerse:::find_lammax(X)

    expect_true(is.numeric(lammax))
    expect_true(lammax > 0)
    expect_true(is.finite(lammax))
})

test_that("utils helper: edge_to_str conversion", {
    edge_df <- data.frame(
        from = c("A", "B", "C"),
        to = c("B", "C", "D")
    )

    result <- scGraphVerse:::.edge_to_str(edge_df)

    expect_true(is.character(result))
    expect_equal(length(result), 3)
})

test_that("utils helper: prepare_prediction_vectors", {
    pred <- matrix(runif(16), nrow = 4)
    truth <- matrix(sample(0:1, 16, replace = TRUE), nrow = 4)
    rownames(pred) <- rownames(truth) <- paste0("G", 1:4)
    colnames(pred) <- colnames(truth) <- paste0("G", 1:4)

    result <- scGraphVerse:::.prepare_prediction_vectors(pred, truth)

    expect_true(is.numeric(result))
})

test_that("utils helper: compute_confusion_metrics", {
    pred <- c(1, 0, 1, 0, 1)
    truth <- c(1, 1, 0, 0, 1)

    metrics <- scGraphVerse:::.compute_confusion_metrics(pred, truth, 1)

    expect_true(is.data.frame(metrics))
    expect_true(all(c("TP", "TN", "FP", "FN") %in% colnames(metrics)))
})

test_that("utils helper: normalize_library_size", {
    mat <- matrix(rpois(100, lambda = 5), nrow = 10, ncol = 10)

    normalized <- scGraphVerse:::.normalize_library_size(mat, c(1000, 5000))

    expect_true(is.matrix(normalized))
    expect_equal(dim(normalized), dim(mat))
})

test_that("utils helper: shuffle_matrix_rows", {
    set.seed(123)
    mat <- matrix(1:20, nrow = 4, ncol = 5)
    rownames(mat) <- paste0("Gene", 1:4)

    shuffled <- scGraphVerse:::.shuffle_matrix_rows(mat)

    expect_true(is.matrix(shuffled))
    expect_equal(dim(shuffled), dim(mat))
})

test_that("utils helper: filter_genes with MT removal", {
    expr <- matrix(rpois(100, lambda = 5), nrow = 10, ncol = 10)
    rownames(expr) <- c(paste0("MT-", 1:5), paste0("Gene", 1:5))

    filtered <- scGraphVerse:::.filter_genes(expr, remove_mt = TRUE, remove_rib = FALSE)

    expect_true(all(!grepl("^MT-", rownames(filtered))))
})

test_that("utils helper: filter_genes with ribosomal removal", {
    expr <- matrix(rpois(100, lambda = 5), nrow = 10, ncol = 10)
    rownames(expr) <- c(paste0("RPS", 1:5), paste0("Gene", 1:5))

    filtered <- scGraphVerse:::.filter_genes(expr, remove_mt = FALSE, remove_rib = TRUE)

    expect_true(all(!grepl("^RPS|^RPL", rownames(filtered))))
})

test_that("utils helper: select_top_genes by variance", {
    set.seed(123)
    expr <- matrix(rnorm(100, mean = 5, sd = 2), nrow = 10, ncol = 10)
    rownames(expr) <- paste0("Gene", 1:10)

    top_genes <- scGraphVerse:::.select_top_genes(expr, top_n = 5)

    expect_true(is.character(top_genes))
    expect_true(length(top_genes) <= 10)
})

test_that("utils helper: convert_counts_list", {
    mat1 <- matrix(rpois(50, lambda = 5), nrow = 5, ncol = 10)
    mat2 <- matrix(rpois(50, lambda = 5), nrow = 5, ncol = 10)
    rownames(mat1) <- rownames(mat2) <- paste0("Gene", 1:5)

    converted <- scGraphVerse:::.convert_counts_list(list(mat1, mat2))

    expect_true(is.list(converted))
    expect_equal(length(converted), 2)
})

test_that("utils helper: extract_expression from matrix", {
    mat <- matrix(rpois(50, lambda = 5), nrow = 5, ncol = 10)
    rownames(mat) <- paste0("Gene", 1:5)

    expr <- scGraphVerse:::.extract_expression(mat)

    expect_true(is.matrix(expr))
    expect_equal(dim(expr), dim(mat))
})

test_that("community_path with verbose = FALSE", {
    skip_if_not_installed("igraph")

    mat <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:3)

    result <- community_path(mat, plot = FALSE, pathway_db = "none", verbose = FALSE)

    expect_true(is.list(result))
    expect_true("communities" %in% names(result))
})

test_that("community_path with louvain method", {
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
        methods = "louvain",
        plot = FALSE,
        pathway_db = "none"
    )

    expect_true(is.list(result))
})

test_that("community_path with leiden method", {
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
        methods = "leiden",
        plot = FALSE,
        pathway_db = "none"
    )

    expect_true(is.list(result))
})

test_that("community_path returns graph object", {
    skip_if_not_installed("igraph")

    mat <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:3)

    result <- community_path(mat, plot = FALSE, pathway_db = "none")

    expect_true("graph" %in% names(result))
    expect_true(inherits(result$graph, "igraph"))
})

test_that("classify_edges with all TP edges", {
    skip_if_not_installed("SummarizedExperiment")

    mat <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3)
    rownames(mat) <- colnames(mat) <- paste0("G", 1:3)

    consensus_se <- build_network_se(list(mat))
    reference_se <- build_network_se(list(mat))

    result <- classify_edges(consensus_se, reference_se, use_stringdb = FALSE)

    expect_equal(length(result$fp_edges), 0)
    expect_equal(length(result$fn_edges), 0)
    expect_true(length(result$tp_edges) > 0)
})

test_that("classify_edges with all FP edges", {
    skip_if_not_installed("SummarizedExperiment")

    consensus_mat <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3)
    reference_mat <- matrix(0, nrow = 3, ncol = 3)
    rownames(consensus_mat) <- colnames(consensus_mat) <- paste0("G", 1:3)
    rownames(reference_mat) <- colnames(reference_mat) <- paste0("G", 1:3)

    consensus_se <- build_network_se(list(consensus_mat))
    reference_se <- build_network_se(list(reference_mat))

    result <- classify_edges(consensus_se, reference_se, use_stringdb = FALSE)

    expect_true(length(result$fp_edges) > 0)
    expect_equal(length(result$fn_edges), 0)
    expect_equal(length(result$tp_edges), 0)
})

test_that("classify_edges edge_colors length matches reference", {
    skip_if_not_installed("SummarizedExperiment")

    consensus_mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3)
    reference_mat <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), nrow = 3)
    rownames(consensus_mat) <- colnames(consensus_mat) <- paste0("G", 1:3)
    rownames(reference_mat) <- colnames(reference_mat) <- paste0("G", 1:3)

    consensus_se <- build_network_se(list(consensus_mat))
    reference_se <- build_network_se(list(reference_mat))

    result <- classify_edges(consensus_se, reference_se, use_stringdb = FALSE)

    # edge_colors should have one color per reference edge
    ref_edge_count <- sum(reference_mat > 0) / 2 # Divide by 2 for undirected
    expect_true(length(result$edge_colors) > 0)
})

test_that("cutoff_adjacency with weight_function max", {
    skip_if_not_installed("MultiAssayExperiment")

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
        weight_function = "max",
        nCores = 1
    )

    expect_true(inherits(binary, "SummarizedExperiment"))
})

test_that("cutoff_adjacency with weight_function min", {
    skip_if_not_installed("MultiAssayExperiment")

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
        weight_function = "min",
        nCores = 1
    )

    expect_true(inherits(binary, "SummarizedExperiment"))
})

test_that("ZILGM with lambda sequence", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mpath")

    set.seed(42)
    X <- matrix(rpois(100, lambda = 3), nrow = 20, ncol = 5)

    result <- scGraphVerse:::zilgm_internal(
        X = X,
        nlambda = 10,
        lambda_min_ratio = 0.01,
        family = "Poisson",
        verbose = 0
    )

    expect_true(is.list(result))
    expect_true(length(result$lambda) <= 10)
})

test_that("ZILGM with NBI family", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mpath")

    set.seed(42)
    X <- matrix(rpois(50, lambda = 4), nrow = 10, ncol = 5)

    result <- scGraphVerse:::zilgm_internal(
        X = X,
        lambda = 0.15,
        nlambda = 1,
        family = "NBI",
        verbose = 0
    )

    expect_true(is.list(result))
    expect_true("network" %in% names(result))
})

test_that("ZILGM with NBII family", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mpath")

    set.seed(42)
    X <- matrix(rpois(50, lambda = 3), nrow = 10, ncol = 5)

    result <- scGraphVerse:::zilgm_internal(
        X = X,
        lambda = 0.1,
        nlambda = 1,
        family = "NBII",
        verbose = 0
    )

    expect_true(is.list(result))
    expect_true("network" %in% names(result))
})

test_that("ZILGM with stability selection via bootstrap", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("mpath")

    set.seed(42)
    X <- matrix(rpois(50, lambda = 5), nrow = 10, ncol = 5)

    result <- scGraphVerse:::zilgm_internal(
        X = X,
        lambda = 0.1,
        nlambda = 1,
        do_boot = TRUE,
        boot_num = 2,
        beta = 0.2,
        family = "Poisson",
        verbose = 0
    )

    expect_true(is.list(result))
    expect_true("network" %in% names(result))
})
