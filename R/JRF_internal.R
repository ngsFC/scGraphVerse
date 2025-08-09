# =========================
# Internal helper for JRF
# =========================

.jrf_onetarget <- function(x, y, ntree, mtry, sampsize, nclasses,
                        importance = TRUE, nCores = 1) {
    totsize <- sum(sampsize)
    p_per_class <- nrow(x) / nclasses

    joint_importance <- array(0, dim = c(p_per_class, nclasses))

    use_parallel <- (nCores > 1) &&
        requireNamespace("BiocParallel", quietly = TRUE) && (ntree > 1)

    if (use_parallel) {
        # PARALLEL TREE BUILDING with BiocParallel
        # Setup parallel backend
        if (nCores > 1) {
            bp_param <- BiocParallel::MulticoreParam(workers = nCores)
        } else {
            bp_param <- BiocParallel::SerialParam()
        }

        # Build trees in parallel
        tree_results <- BiocParallel::bplapply(seq_len(ntree),
            function(tree_idx) {
                # Joint bootstrap sampling for this tree
                bootstrap_samples <- .joint_bootstrap_sample(sampsize, nclasses)

                # Build one joint tree with full recursive splitting
                tree_importance <- .build_joint_tree(
                    x, y, mtry, sampsize, nclasses,
                    bootstrap_samples, p_per_class,
                    nodesize = 5
                )

                return(tree_importance)
            },
            BPPARAM = bp_param
        )

        # Aggregate results from parallel workers
        for (tree_result in tree_results) {
            joint_importance <- joint_importance + tree_result
        }
    } else {
        for (tree_idx in seq_len(ntree)) {
            # Joint bootstrap sampling across all classes
            bootstrap_samples <- .joint_bootstrap_sample(sampsize, nclasses)

            # Build one joint tree with full recursive splitting
            tree_importance <- .build_joint_tree(
                x, y, mtry, sampsize, nclasses,
                bootstrap_samples, p_per_class,
                nodesize = 5
            )

            # Accumulate importance scores
            joint_importance <- joint_importance + tree_importance
        }
    }

    joint_importance <- joint_importance / ntree

    # Flatten importance matrix
    flattened_importance <- numeric(p_per_class * nclasses)
    for (s in seq_len(nclasses)) {
        start_idx <- (s - 1) * p_per_class + 1
        end_idx <- s * p_per_class
        flattened_importance[start_idx:end_idx] <- joint_importance[, s]
    }

    imp_matrix <- matrix(0, nrow = length(flattened_importance), ncol = 2)
    imp_matrix[, 1] <- flattened_importance
    imp_matrix[, 2] <- flattened_importance

    out <- list(importance = imp_matrix)
    class(out) <- "randomForest"
    return(out)
}

# bootstrap sampling
.joint_bootstrap_sample <- function(sampsize, nclasses) {
    bootstrap_idx <- vector("list", nclasses)

    # bootstrap logic
    for (s in seq_len(nclasses)) {
        bootstrap_idx[[s]] <- sample(seq_len(sampsize[s]), sampsize[s],
            replace = TRUE
        )
    }

    return(bootstrap_idx)
}

# Build single joint tree with full recursive splitting
.build_joint_tree <- function(x, y, mtry, sampsize, nclasses,
                        bootstrap_samples, p_per_class, nodesize = 5,
                        maxnodes = NULL) {
    # Initialize importance accumulator
    tree_importance <- array(0, dim = c(p_per_class, nclasses))

    total_samples <- sum(sampsize)

    # Calculate maximum nodes if not specified
    if (is.null(maxnodes)) {
        maxnodes <- 2 * floor(total_samples / max(1, nodesize - 4)) + 1
    } else {
        maxnodes <- 2 * maxnodes - 1 # Convert terminal nodes to total nodes
    }

    # Initialize node tracking structures
    node_queue <- list()
    node_counter <- 1

    # Create root node for each class
    root_nodes <- list()
    for (s in seq_len(nclasses)) {
        root_nodes[[s]] <- list(
            node_id = 1,
            class_id = s,
            start_idx = 1,
            end_idx = sampsize[s],
            sample_indices = bootstrap_samples[[s]],
            status = "to_split",
            depth = 0
        )
    }

    # Add root nodes to queue
    node_queue[[1]] <- root_nodes

    # Main tree building loop
    while (length(node_queue) > 0 && node_counter < maxnodes) {
        current_node_set <- node_queue[[1]]
        node_queue <- node_queue[-1]

        # Check if any nodes in this set need splitting
        nodes_to_split <- vapply(
            current_node_set,
            function(node) node$status == "to_split",
            logical(1)
        )
        if (!any(nodes_to_split)) next

        # Find best variable across all classes for joint splitting
        best_split <- .find_best_joint_variable(
            x, y, current_node_set, mtry,
            sampsize, nclasses, p_per_class,
            nodesize
        )

        if (is.null(best_split) || best_split$improvement <= 0) {
            # Mark all nodes as terminal
            for (s in seq_len(nclasses)) {
                if (current_node_set[[s]]$status == "to_split") {
                    current_node_set[[s]]$status <- "terminal"
                }
            }
            next
        }

        # Accumulate variable importance
        for (s in seq_len(nclasses)) {
            if (best_split$valid_classes[s] &&
                best_split$class_contributions[s] > 0) {
                tree_importance[best_split$variable, s] <-
                    tree_importance[best_split$variable, s] +
                    best_split$class_contributions[s]
            }
        }

        # Create child nodes
        left_nodes <- list()
        right_nodes <- list()

        for (s in seq_len(nclasses)) {
            node <- current_node_set[[s]]
            if (node$status == "to_split" && best_split$valid_classes[s]) {
                # Get split threshold for this class
                threshold <- best_split$thresholds[s]

                # Split samples based on threshold
                split_result <- .split_node_samples(
                    x, y, node,
                    best_split$variable,
                    threshold, s, nclasses,
                    p_per_class
                )

                # Create left child node
                left_nodes[[s]] <- list(
                    node_id = node_counter + 1,
                    class_id = s,
                    sample_indices = split_result$left_samples,
                    status = if (length(split_result$left_samples)<=nodesize) {
                        "terminal"
                    } else {
                        "to_split"
                    },
                    depth = node$depth + 1,
                    parent_id = node$node_id
                )

                # Create right child node
                right_nodes[[s]] <- list(
                    node_id = node_counter + 2,
                    class_id = s,
                    sample_indices = split_result$right_samples,
                    status = if (length(split_result$right_samples)<=nodesize) {
                        "terminal"
                    } else {
                        "to_split"
                    },
                    depth = node$depth + 1,
                    parent_id = node$node_id
                )
            } else {
                # Node doesn't split - create terminal copies
                left_nodes[[s]] <- node
                right_nodes[[s]] <- node
                left_nodes[[s]]$status <- "terminal"
                right_nodes[[s]]$status <- "terminal"
            }
        }

        # Add child nodes to queue if they need splitting
        left_need_split <- any(vapply(
            left_nodes, function(n) n$status == "to_split", logical(1)
        ))
        if (left_need_split) {
            node_queue <- c(node_queue, list(left_nodes))
        }

        right_need_split <- any(vapply(
            right_nodes, function(n) n$status == "to_split", logical(1)
        ))
        if (right_need_split) {
            node_queue <- c(node_queue, list(right_nodes))
        }

        node_counter <- node_counter + 2
    }

    return(tree_importance)
}

# Evaluate joint split criterion for a variable across all classes
.evaluate_joint_split_criterion <- function(x, y, var_idx, sampsize, nclasses,
                                            bootstrap_samples, p_per_class) {
    class_contributions <- numeric(nclasses)
    valid_classes <- logical(nclasses)

    # Class-specific split criteria (like critvar[s] in original C)
    class_criteria <- numeric(nclasses)
    class_weights <- numeric(nclasses)
    raw_class_scores <- numeric(nclasses)

    total_samples <- sum(sampsize)

    for (s in seq_len(nclasses)) {
        # Extract data for this class and variable
        var_row_idx <- (s - 1) * p_per_class + var_idx

        if (var_row_idx > nrow(x)) {
            valid_classes[s] <- FALSE
            next
        }

        # Get predictor and response values for this class
        class_start_col <- if (s == 1) 1 else sum(sampsize[seq_len(s - 1)]) + 1
        class_end_col <- sum(sampsize[seq_len(s)])

        predictor_vals <- x[var_row_idx, class_start_col:class_end_col]
        response_vals <- y[s, class_start_col:class_end_col]

        # Apply bootstrap sampling
        boot_pred <- predictor_vals[bootstrap_samples[[s]]]
        boot_resp <- response_vals[bootstrap_samples[[s]]]

        # Remove invalid values
        valid_idx <- is.finite(boot_pred) & is.finite(boot_resp)
        if (sum(valid_idx) < 3) {
            valid_classes[s] <- FALSE
            next
        }

        boot_pred <- boot_pred[valid_idx]
        boot_resp <- boot_resp[valid_idx]
        valid_classes[s] <- TRUE

        # Calculate split criterion for this class
        raw_score <- .calculate_split_criterion(boot_pred, boot_resp)
        raw_class_scores[s] <- raw_score
        class_criteria[s] <- raw_score
        class_weights[s] <- sampsize[s] / total_samples # class weighting
    }

    valid_idx <- which(valid_classes)
    joint_score <- 0
    if (length(valid_idx) > 0) {
        joint_score <- (sum(class_weights[valid_idx] *
            class_criteria[valid_idx]) / nclasses)

        class_contributions[valid_idx] <- class_criteria[valid_idx]
    }

    return(list(
        class_contributions = class_contributions,
        valid_classes = valid_classes,
        joint_score = joint_score,
        raw_class_scores = raw_class_scores,
        class_weights = class_weights
    ))
}

.calculate_split_criterion <- function(predictor_vals, response_vals) {
    if (length(unique(predictor_vals)) < 2) {
        return(0)
    }

    # Sort values to find best split point
    sorted_idx <- order(predictor_vals)
    sorted_pred <- predictor_vals[sorted_idx]
    sorted_resp <- response_vals[sorted_idx]

    # Calculate parent node variance
    n_total <- length(sorted_resp)
    if (n_total < 2) {
        return(0)
    }

    parent_var <- var(sorted_resp)
    best_criterion <- 0

    for (i in seq_len(length(sorted_pred) - 1)) {
        if (sorted_pred[i] >= sorted_pred[i + 1]) next

        # Left and right node data
        left_vals <- sorted_resp[seq_len(i)]
        right_vals <- sorted_resp[seq(i + 1, n_total)]

        # Require minimum node size (at least 1 observation per node)
        if (length(left_vals) == 0 || length(right_vals) == 0) next

        # Calculate variances for each child node
        left_var <- if (length(left_vals) == 1) 0 else var(left_vals)
        right_var <- if (length(right_vals) == 1) 0 else var(right_vals)

        # Calculate node weights
        left_weight <- length(left_vals) / n_total
        right_weight <- length(right_vals) / n_total

        # Variance reduction criterion (standard Random Forest)
        weighted_var <- left_weight * left_var + right_weight * right_var
        split_criterion <- parent_var - weighted_var

        if (split_criterion > best_criterion) {
            best_criterion <- split_criterion
        }
    }

    return(best_criterion)
}

# Find best variable for joint splitting across all classes (full tree version)
.find_best_joint_variable <- function(x, y, node_set, mtry, sampsize, nclasses,
                                p_per_class, nodesize) {
    # Variable selection with shuffling 
    available_vars <- seq_len(p_per_class)
    last <- p_per_class - 1

    best_improvement <- 0
    best_variable <- NULL
    best_thresholds <- numeric(nclasses)
    best_contributions <- numeric(nclasses)
    best_valid_classes <- logical(nclasses)

    # Try mtry variables
    for (var_iter in seq_len(min(mtry, p_per_class))) {
        if (last < 0) break

        # Random variable selection with array swapping
        j <- sample(0:last, 1)
        var_idx <- available_vars[j + 1]
        available_vars[j + 1] <- available_vars[last + 1]
        last <- last - 1

        # Evaluate this variable for joint splitting
        joint_result <- .evaluate_joint_variable_split(
            x, y, node_set, var_idx,
            sampsize, nclasses, p_per_class
        )

        if (joint_result$joint_improvement > best_improvement) {
            best_improvement <- joint_result$joint_improvement
            best_variable <- var_idx
            best_thresholds <- joint_result$thresholds
            best_contributions <- joint_result$class_contributions
            best_valid_classes <- joint_result$valid_classes
        }
    }

    if (is.null(best_variable) || best_improvement <= 0) {
        return(NULL)
    }

    return(list(
        variable = best_variable,
        improvement = best_improvement,
        thresholds = best_thresholds,
        class_contributions = best_contributions,
        valid_classes = best_valid_classes
    ))
}

# Evaluate a specific variable for joint splitting
.evaluate_joint_variable_split <- function(x, y, node_set, var_idx,
                                sampsize, nclasses, p_per_class) {
    class_contributions <- numeric(nclasses)
    valid_classes <- logical(nclasses)
    thresholds <- numeric(nclasses)
    class_weights <- sampsize / sum(sampsize)

    # Evaluate split for each class
    for (s in seq_len(nclasses)) {
        node <- node_set[[s]]

        if (node$status != "to_split" || length(node$sample_indices) <= 1) {
            valid_classes[s] <- FALSE
            next
        }

        # Extract data for this class and variable
        var_row_idx <- (s - 1) * p_per_class + var_idx

        if (var_row_idx > nrow(x)) {
            valid_classes[s] <- FALSE
            next
        }

        # Get predictor and response values using node's sample indices
        class_start_col <- if (s == 1) {
            1
        } else {
            sum(sampsize[seq_len(s - 1)]) + 1
        }
        class_end_col <- sum(sampsize[seq_len(s)])

        # Apply node's bootstrap sampling
        class_cols <- class_start_col:class_end_col
        predictor_vals <- x[var_row_idx, class_cols][node$sample_indices]
        response_vals <- y[s, class_cols][node$sample_indices]

        # Remove invalid values
        valid_idx <- is.finite(predictor_vals) & is.finite(response_vals)
        if (sum(valid_idx) < 3) {
            valid_classes[s] <- FALSE
            next
        }

        predictor_vals <- predictor_vals[valid_idx]
        response_vals <- response_vals[valid_idx]

        # Find best threshold for this class (matching original C logic)
        split_result <- .find_best_threshold_for_class(
            predictor_vals, response_vals
        )

        if (split_result$improvement > 0) {
            class_contributions[s] <- split_result$improvement
            thresholds[s] <- split_result$threshold
            valid_classes[s] <- TRUE
        } else {
            valid_classes[s] <- FALSE
        }
    }

    # Calculate joint improvement (weighted by class size)
    valid_idx <- which(valid_classes)
    joint_improvement <- 0

    if (length(valid_idx) > 0) {
        weighted_contributions <- class_weights[valid_idx] *
            class_contributions[valid_idx]
        joint_improvement <- sum(weighted_contributions) / nclasses
    }

    return(list(
        joint_improvement = joint_improvement,
        class_contributions = class_contributions,
        valid_classes = valid_classes,
        thresholds = thresholds
    ))
}

# Find best threshold for a single class (matching original C algorithm)
.find_best_threshold_for_class <- function(predictor_vals, response_vals) {
    if (length(unique(predictor_vals)) < 2) {
        return(list(improvement = 0, threshold = mean(predictor_vals)))
    }

    # Sort values (matching original C qsort)
    sorted_idx <- order(predictor_vals)
    sorted_pred <- predictor_vals[sorted_idx]
    sorted_resp <- response_vals[sorted_idx]

    n_total <- length(sorted_resp)
    sum_total <- sum(sorted_resp)

    # Parent criterion (sum of squares)
    parent_criterion <- sum_total * sum_total / n_total

    best_improvement <- 0
    best_threshold <- mean(range(sorted_pred))

    sum_left <- 0
    n_left <- 0

    # Search through gaps (matching original C algorithm)
    for (i in seq_len(n_total - 1)) {
        sum_left <- sum_left + sorted_resp[i]
        n_left <- n_left + 1

        sum_right <- sum_total - sum_left
        n_right <- n_total - n_left

        # Only consider if there's a gap between values
        if (sorted_pred[i] < sorted_pred[i + 1] && n_left > 0 && n_right > 0) {
            # Calculate improvement (matching original C formula)
            left_criterion <- sum_left * sum_left / n_left
            right_criterion <- sum_right * sum_right / n_right
            improvement <- left_criterion + right_criterion - parent_criterion

            if (improvement > best_improvement) {
                best_improvement <- improvement
                best_threshold <- (sorted_pred[i] + sorted_pred[i + 1]) / 2
            }
        }
    }

    return(list(
        improvement = best_improvement,
        threshold = best_threshold
    ))
}

# Split node samples based on threshold
.split_node_samples <- function(x, y, node, var_idx, threshold,
                                class_id, nclasses, p_per_class) {
    # Get predictor values for this variable and class
    var_row_idx <- (class_id - 1) * p_per_class + var_idx

    max_indices <- max(node$sample_indices)
    sampsize_cumsum <- c(0, cumsum(rep(max_indices, nclasses)))
    class_start_col <- sampsize_cumsum[class_id] + 1
    class_end_col <- sampsize_cumsum[class_id + 1]

    # Get predictor values using node's sample indices
    class_cols <- class_start_col:class_end_col
    predictor_vals <- x[var_row_idx, class_cols][node$sample_indices]

    # Split samples based on threshold
    left_mask <- predictor_vals <= threshold
    right_mask <- !left_mask

    left_samples <- node$sample_indices[left_mask]
    right_samples <- node$sample_indices[right_mask]

    return(list(
        left_samples = left_samples,
        right_samples = right_samples
    ))
}

# Main JRF network inference
.jrf_network <- function(data_list, ntree = 1000, mtry = NULL, nCores = 1) {
    nclasses <- length(data_list)
    sampsize <- vapply(data_list, ncol, FUN.VALUE = integer(1))
    p <- nrow(data_list[[1]])
    if (is.null(mtry)) mtry <- max(floor(sqrt(p - 1)), 1)

    genes.name <- rownames(data_list[[1]])
    if (is.null(genes.name)) genes.name <- paste0("G", seq_len(p))

    imp <- array(0, dim = c(p, length(genes.name), nclasses))

    for (j in seq_len(p)) {
        # Build covariate matrix: (p-1) * nclasses rows, totsize cols
        totsize <- sum(sampsize)
        covar <- matrix(0, (p - 1) * nclasses, totsize)
        y <- matrix(0, nclasses, totsize)

        # Fill data matrices
        for (c in seq_len(nclasses)) {
            idx_start <- if (c == 1) 1 else sum(sampsize[seq_len(c - 1)]) + 1
            idx_end <- sum(sampsize[seq_len(c)])

            # Target gene response for class c
            y[c, idx_start:idx_end] <- data_list[[c]][j, ]

            # Predictor genes for class c (excluding target gene j)
            covar[
                ((c - 1) * (p - 1) + 1):(c * (p - 1)),
                idx_start:idx_end
            ] <- data_list[[c]][-j, ]
        }

        # Run TRUE joint RF - joint modeling with full trees
        jrf.out <- .jrf_onetarget(
            x = covar, y = y, ntree = ntree, mtry = mtry,
            sampsize = sampsize, nclasses = nclasses, nCores = nCores
        )

        # Extract importance scores per class
        for (s in seq_len(nclasses)) {
            start_idx <- (p - 1) * (s - 1) + 1
            end_idx <- (p - 1) * s

            if (end_idx <= nrow(jrf.out$importance)) {
                imp[-j, j, s] <- jrf.out$importance[start_idx:end_idx, 2]
            } else {
                # Handle edge case - fill with small random values
                imp[-j, j, s] <- runif(p - 1, 0, 0.01)
            }
        }
    }

    # Create symmetric importance matrix
    imp.final <- matrix(0, nrow = p * (p - 1) / 2, ncol = nclasses)
    vec1 <- matrix(rep(genes.name, p), p, p)
    vec2 <- t(vec1)
    vec1 <- vec1[lower.tri(vec1, diag = FALSE)]
    vec2 <- vec2[lower.tri(vec2, diag = FALSE)]

    for (s in seq_len(nclasses)) {
        imp.s <- imp[, , s]
        t.imp <- t(imp.s)
        imp.final[, s] <- (imp.s[lower.tri(imp.s, diag = FALSE)] +
            t.imp[lower.tri(t.imp, diag = FALSE)]) / 2
    }

    result_list <- vector("list", nclasses)

    for (s in seq_len(nclasses)) {
        result_list[[s]] <- data.frame(
            gene1 = as.character(vec1),
            gene2 = as.character(vec2),
            importance = imp.final[, s],
            stringsAsFactors = FALSE
        )
    }

    return(result_list)
}
