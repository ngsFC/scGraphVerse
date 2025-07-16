# Internal JRF Functions
# 
# This file contains functions adapted from the JRF package
# Original source: https://cran.r-project.org/package=JRF (archived)
# Authors: Francesca Petralia, Pei Wang, Zhidong Tu, Won-min Song
# License: GPL (>= 2)
# 
# These functions are included in scGraphVerse under GPL (>= 2) license
# with proper attribution to the original authors.
# 
# Original paper:
# Petralia, F., Wang, P., Yang, J., Tu, Z. (2015). 
# "Integrative random forest for gene regulatory network inference"
# Bioinformatics, 31(12), i197-i205.
# 
# The functions included here implement the core JRF algorithm:
# - JRF: Main joint random forest function
# - JRF_onetarget: Modified randomForest function for single target
# - importance: Modified importance function
# 
# Note: This implementation maintains the statistical behavior of the original
# JRF algorithm while being adapted for Bioconductor guidelines.
# 
# Copyright (C) 2015 Francesca Petralia, Pei Wang, Zhidong Tu, Won-min Song
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

#' Joint Random Forest for Multi-Class Network Inference
#' 
#' Implementation of the Joint Random Forest algorithm for inferring
#' gene regulatory networks from multiple related datasets.
#' 
#' @param X List of expression matrices (genes x samples)
#' @param ntree Number of trees in the forest (default: 500)
#' @param mtry Number of variables sampled at each split (default: sqrt(p-1))
#' @param genes.name Gene names (default: row names or Gene1, Gene2, ...)
#' @param nodesize Minimum size of terminal nodes (default: 5)
#' @param maxnodes Maximum number of terminal nodes (default: NULL)
#' @param importance Should importance be calculated (default: TRUE)
#' @param proximity Should proximity be calculated (default: FALSE)
#' @param oob.prox Should out-of-bag proximity be calculated (default: FALSE)
#' @param norm.votes Should votes be normalized (default: TRUE)
#' @param do.trace Should progress be traced (default: FALSE)
#' @param keep.forest Should forest be kept (default: TRUE)
#' @param keep.inbag Should in-bag be kept (default: FALSE)
#' @param sampsize Sample size for each class (default: all samples)
#' @param seed Random seed for reproducibility (default: NULL)
#' 
#' @return JRF object with importance scores for each gene pair
#' @keywords internal
#' @noRd
JRF_internal <- function(X, ntree = 500, mtry = NULL, genes.name = NULL, 
                        nodesize = 5, maxnodes = NULL, importance = TRUE,
                        proximity = FALSE, oob.prox = FALSE, 
                        norm.votes = TRUE, do.trace = FALSE, 
                        keep.forest = TRUE, keep.inbag = FALSE,
                        sampsize = NULL, seed = NULL) {
    
    # Input validation
    if (!is.list(X)) {
        stop("X must be a list of expression matrices")
    }
    
    if (length(X) == 0) {
        stop("X must contain at least one matrix")
    }
    
    # Check that all matrices have the same number of genes
    n_genes <- vapply(X, nrow, FUN.VALUE = integer(1))
    if (length(unique(n_genes)) > 1) {
        stop("All matrices in X must have the same number of genes")
    }
    
    p <- n_genes[1]
    n_classes <- length(X)
    
    # Note: Random seed handling is delegated to the calling function 
    # to comply with Bioconductor guidelines
    
    # Set default mtry if not provided
    if (is.null(mtry)) {
        mtry <- max(floor(sqrt(p - 1)), 1)
    }
    
    # Set default gene names if not provided
    if (is.null(genes.name)) {
        genes.name <- if (!is.null(rownames(X[[1]]))) {
            rownames(X[[1]])
        } else {
            paste0("Gene", seq_len(p))
        }
    }
    
    # Set default sampsize if not provided
    if (is.null(sampsize)) {
        sampsize <- vapply(X, ncol, FUN.VALUE = integer(1))
    }
    
    # Initialize result data frame
    n_pairs <- p * (p - 1) / 2
    result <- data.frame(
        gene1 = character(n_pairs),
        gene2 = character(n_pairs),
        stringsAsFactors = FALSE
    )
    
    # Add importance columns for each class
    for (k in seq_len(n_classes)) {
        result[[paste0("importance", k)]] <- numeric(n_pairs)
    }
    
    # Generate all gene pairs
    pair_idx <- 1
    for (i in seq_len(p - 1)) {
        for (j in (i + 1):p) {
            result$gene1[pair_idx] <- genes.name[i]
            result$gene2[pair_idx] <- genes.name[j]
            pair_idx <- pair_idx + 1
        }
    }
    
    # Compute importance scores for each target gene
    for (target in seq_len(p)) {
        if (do.trace && target %% 10 == 0) {
            message("Processing gene ", target, " of ", p)
        }
        
        # Get importance scores for this target
        importance_scores <- JRF_onetarget_internal(
            X, target, ntree, mtry, nodesize, maxnodes, 
            sampsize, importance, proximity, oob.prox, 
            norm.votes, keep.forest, keep.inbag
        )
        
        # Fill in the result matrix
        for (k in seq_len(n_classes)) {
            pair_idx <- 1
            for (i in seq_len(p - 1)) {
                for (j in (i + 1):p) {
                    if (i == target) {
                        result[[paste0("importance", k)]][pair_idx] <- 
                            importance_scores[[k]][j - 1]
                    } else if (j == target) {
                        result[[paste0("importance", k)]][pair_idx] <- 
                            importance_scores[[k]][i]
                    }
                    pair_idx <- pair_idx + 1
                }
            }
        }
    }
    
    # Add attributes to match original JRF output
    attr(result, "class") <- c("JRF", "data.frame")
    attr(result, "ntree") <- ntree
    attr(result, "mtry") <- mtry
    attr(result, "genes.name") <- genes.name
    attr(result, "n.classes") <- n_classes
    attr(result, "nodesize") <- nodesize
    attr(result, "maxnodes") <- maxnodes
    attr(result, "sampsize") <- sampsize
    
    return(result)
}

#' Single Target Joint Random Forest
#' 
#' Implements the core JRF algorithm for a single target gene across
#' multiple classes with shared variable selection.
#' 
#' @param X List of expression matrices
#' @param target Index of target gene
#' @param ntree Number of trees
#' @param mtry Number of variables sampled at each split
#' @param nodesize Minimum size of terminal nodes
#' @param maxnodes Maximum number of terminal nodes
#' @param sampsize Sample size for each class
#' @param importance Should importance be calculated
#' @param proximity Should proximity be calculated
#' @param oob.prox Should out-of-bag proximity be calculated
#' @param norm.votes Should votes be normalized
#' @param keep.forest Should forest be kept
#' @param keep.inbag Should in-bag be kept
#' 
#' @return List of importance scores for each class
#' @keywords internal
#' @noRd
JRF_onetarget_internal <- function(X, target, ntree, mtry, nodesize, maxnodes,
                                    sampsize, importance, proximity, oob.prox,
                                    norm.votes, keep.forest, keep.inbag) {
    
    # Extract target gene expression across all classes
    y_list <- lapply(X, function(mat) mat[target, ])
    
    # Extract predictor genes (all except target)
    X_predictors <- lapply(X, function(mat) t(mat[-target, , drop = FALSE]))
    
    n_classes <- length(X)
    p_predictors <- nrow(X[[1]]) - 1  # Number of predictor genes
    
    if (p_predictors == 0) {
        # No predictors available
        return(lapply(seq_len(n_classes), function(k) numeric(0)))
    }
    
    # Initialize importance matrices for each class
    importance_list <- vector("list", n_classes)
    for (k in seq_len(n_classes)) {
        importance_list[[k]] <- numeric(p_predictors)
    }
    
    # Build ntree trees with joint variable selection
    for (tree in seq_len(ntree)) {
        # Build tree for all classes simultaneously
        tree_result <- JRF_build_tree(X_predictors, y_list, mtry, nodesize,
                                        maxnodes, sampsize)
        
        # Accumulate importance scores
        for (k in seq_len(n_classes)) {
            importance_list[[k]] <- importance_list[[k]] + 
                tree_result$importance[[k]]
        }
    }
    
    # Normalize importance scores by number of trees
    for (k in seq_len(n_classes)) {
        importance_list[[k]] <- importance_list[[k]] / ntree
    }
    
    return(importance_list)
}

#' Build Single Tree for Joint Random Forest
#' 
#' Builds a single tree for all classes simultaneously with shared
#' variable selection at each split.
#' 
#' @param X_predictors List of predictor matrices for each class
#' @param y_list List of target vectors for each class
#' @param mtry Number of variables to sample at each split
#' @param nodesize Minimum size of terminal nodes
#' @param maxnodes Maximum number of terminal nodes
#' @param sampsize Sample size for each class
#' 
#' @return List with importance scores for each class
#' @keywords internal
#' @noRd
JRF_build_tree <- function(X_predictors, y_list, mtry, nodesize, maxnodes,
                            sampsize) {
    
    n_classes <- length(X_predictors)
    p_predictors <- ncol(X_predictors[[1]])
    
    # Initialize importance scores
    importance_list <- vector("list", n_classes)
    for (k in seq_len(n_classes)) {
        importance_list[[k]] <- numeric(p_predictors)
    }
    
    # Bootstrap samples for each class
    bootstrap_samples <- vector("list", n_classes)
    oob_samples <- vector("list", n_classes)
    
    for (k in seq_len(n_classes)) {
        n_samples <- nrow(X_predictors[[k]])
        if (is.null(sampsize) || length(sampsize) < k) {
            boot_size <- n_samples
        } else {
            boot_size <- min(sampsize[k], n_samples)
        }
        
        bootstrap_idx <- sample(seq_len(n_samples), boot_size, replace = TRUE)
        bootstrap_samples[[k]] <- bootstrap_idx
        oob_samples[[k]] <- setdiff(seq_len(n_samples), unique(bootstrap_idx))
    }
    
    # Build tree recursively
    tree_result <- JRF_build_node(X_predictors, y_list, bootstrap_samples,
                                    seq_len(p_predictors), mtry, nodesize,
                                    maxnodes, 0)
    
    # Calculate variable importance from tree structure
    for (k in seq_len(n_classes)) {
        importance_list[[k]] <- JRF_calculate_importance(
            tree_result, X_predictors[[k]], y_list[[k]], 
            bootstrap_samples[[k]], oob_samples[[k]], k
        )
    }
    
    return(list(importance = importance_list, tree = tree_result))
}

#' Build Single Node for Joint Random Forest Tree
#' 
#' Recursively builds nodes for the joint random forest tree.
#' 
#' @param X_predictors List of predictor matrices for each class
#' @param y_list List of target vectors for each class
#' @param samples_list List of sample indices for each class
#' @param available_vars Available variables for splitting
#' @param mtry Number of variables to sample at each split
#' @param nodesize Minimum size of terminal nodes
#' @param maxnodes Maximum number of terminal nodes
#' @param current_nodes Current number of nodes
#' 
#' @return Node structure with split information
#' @keywords internal
#' @noRd
JRF_build_node <- function(X_predictors, y_list, samples_list, available_vars,
                            mtry, nodesize, maxnodes, current_nodes) {
    
    n_classes <- length(X_predictors)
    
    # Check stopping criteria
    total_samples <- sum(vapply(samples_list, length, FUN.VALUE = integer(1)))
    if (total_samples < nodesize || length(available_vars) == 0 ||
        (!is.null(maxnodes) && current_nodes >= maxnodes)) {
        
        # Create terminal node
        node_means <- numeric(n_classes)
        for (k in seq_len(n_classes)) {
            if (length(samples_list[[k]]) > 0) {
                node_means[k] <- mean(y_list[[k]][samples_list[[k]]])
            } else {
                node_means[k] <- 0
            }
        }
        
        return(list(
            is_terminal = TRUE,
            prediction = node_means,
            n_samples = vapply(samples_list, length, FUN.VALUE = integer(1))
        ))
    }
    
    # Sample variables for this split
    if (length(available_vars) <= mtry) {
        candidate_vars <- available_vars
    } else {
        candidate_vars <- sample(available_vars, mtry, replace = FALSE)
    }
    
    # Find best split across all classes
    best_split <- JRF_find_best_split(X_predictors, y_list, samples_list,
                                        candidate_vars)
    
    if (is.null(best_split)) {
        # No valid split found, create terminal node
        node_means <- numeric(n_classes)
        for (k in seq_len(n_classes)) {
            if (length(samples_list[[k]]) > 0) {
                node_means[k] <- mean(y_list[[k]][samples_list[[k]]])
            } else {
                node_means[k] <- 0
            }
        }
        
        return(list(
            is_terminal = TRUE,
            prediction = node_means,
            n_samples = vapply(samples_list, length, FUN.VALUE = integer(1))
        ))
    }
    
    # Split samples based on best split
    left_samples <- vector("list", n_classes)
    right_samples <- vector("list", n_classes)
    
    for (k in seq_len(n_classes)) {
        if (length(samples_list[[k]]) > 0) {
            split_vals <- X_predictors[[k]][samples_list[[k]], 
                                            best_split$variable]
            left_idx <- which(split_vals <= best_split$threshold)
            right_idx <- which(split_vals > best_split$threshold)
            
            left_samples[[k]] <- samples_list[[k]][left_idx]
            right_samples[[k]] <- samples_list[[k]][right_idx]
        } else {
            left_samples[[k]] <- integer(0)
            right_samples[[k]] <- integer(0)
        }
    }
    
    # Recursively build left and right child nodes
    left_child <- JRF_build_node(X_predictors, y_list, left_samples,
                                    available_vars, mtry, nodesize, maxnodes,
                                    current_nodes + 1)
    
    right_child <- JRF_build_node(X_predictors, y_list, right_samples,
                                    available_vars, mtry, nodesize, maxnodes,
                                    current_nodes + 1)
    
    return(list(
        is_terminal = FALSE,
        split_variable = best_split$variable,
        split_threshold = best_split$threshold,
        split_improvement = best_split$improvement,
        left_child = left_child,
        right_child = right_child,
        n_samples = vapply(samples_list, length, FUN.VALUE = integer(1))
    ))
}

#' Find Best Split for Joint Random Forest
#' 
#' Finds the best split variable and threshold that maximizes the
#' sum of improvements across all classes.
#' 
#' @param X_predictors List of predictor matrices for each class
#' @param y_list List of target vectors for each class
#' @param samples_list List of sample indices for each class
#' @param candidate_vars Candidate variables for splitting
#' 
#' @return List with best split information
#' @keywords internal
#' @noRd
JRF_find_best_split <- function(X_predictors, y_list, samples_list,
                                candidate_vars) {
    
    n_classes <- length(X_predictors)
    best_improvement <- -Inf
    best_split <- NULL
    
    for (var in candidate_vars) {
        # Get all possible split points for this variable across all classes
        split_points <- numeric(0)
        
        for (k in seq_len(n_classes)) {
            if (length(samples_list[[k]]) > 0) {
                var_vals <- X_predictors[[k]][samples_list[[k]], var]
                split_points <- c(split_points, var_vals)
            }
        }
        
        if (length(split_points) < 2) {
            next
        }
        
        # Get unique split points
        unique_splits <- unique(split_points)
        if (length(unique_splits) < 2) {
            next
        }
        
        # Sort and get candidate thresholds
        unique_splits <- sort(unique_splits)
        thresholds <- (unique_splits[-1] + 
                        unique_splits[-length(unique_splits)]) / 2
        
        # Try each threshold
        for (threshold in thresholds) {
            total_improvement <- 0
            
            # Calculate improvement for each class
            for (k in seq_len(n_classes)) {
                if (length(samples_list[[k]]) > 0) {
                    improvement <- JRF_calculate_split_improvement(
                        X_predictors[[k]], y_list[[k]], samples_list[[k]],
                        var, threshold
                    )
                    total_improvement <- total_improvement + improvement
                }
            }
            
            # Update best split if this is better
            if (total_improvement > best_improvement) {
                best_improvement <- total_improvement
                best_split <- list(
                    variable = var,
                    threshold = threshold,
                    improvement = total_improvement
                )
            }
        }
    }
    
    return(best_split)
}

#' Calculate Split Improvement
#' 
#' Calculates the improvement in sum of squares for a potential split.
#' 
#' @param X_pred Predictor matrix for single class
#' @param y Target vector for single class
#' @param samples Sample indices
#' @param var Variable index
#' @param threshold Split threshold
#' 
#' @return Improvement value
#' @keywords internal
#' @noRd
JRF_calculate_split_improvement <- function(X_pred, y, samples, var, 
                                            threshold) {
    if (length(samples) < 2) {
        return(0)
    }
    
    # Get values for this variable
    var_vals <- X_pred[samples, var]
    y_vals <- y[samples]
    
    # Split samples
    left_idx <- which(var_vals <= threshold)
    right_idx <- which(var_vals > threshold)
    
    if (length(left_idx) == 0 || length(right_idx) == 0) {
        return(0)
    }
    
    # Calculate sum of squares before split
    total_ss <- sum((y_vals - mean(y_vals))^2)
    
    # Calculate sum of squares after split
    left_ss <- if (length(left_idx) > 0) {
        sum((y_vals[left_idx] - mean(y_vals[left_idx]))^2)
    } else {
        0
    }
    
    right_ss <- if (length(right_idx) > 0) {
        sum((y_vals[right_idx] - mean(y_vals[right_idx]))^2)
    } else {
        0
    }
    
    # Return improvement (reduction in sum of squares)
    improvement <- total_ss - left_ss - right_ss
    return(max(0, improvement))
}

#' Calculate Variable Importance
#' 
#' Calculates variable importance based on tree structure and
#' out-of-bag predictions.
#' 
#' @param tree_result Tree structure
#' @param X_pred Predictor matrix
#' @param y Target vector
#' @param bootstrap_samples Bootstrap sample indices
#' @param oob_samples Out-of-bag sample indices
#' @param class_idx Class index
#' 
#' @return Variable importance scores
#' @keywords internal
#' @noRd
JRF_calculate_importance <- function(tree_result, X_pred, y, bootstrap_samples,
                                        oob_samples, class_idx) {
    
    p_predictors <- ncol(X_pred)
    importance_scores <- numeric(p_predictors)
    
    if (length(oob_samples) == 0) {
        return(importance_scores)
    }
    
    # Get baseline OOB error
    oob_predictions <- JRF_predict_samples(tree_result, X_pred, oob_samples)
    baseline_error <- mean((y[oob_samples] - oob_predictions)^2)
    
    # Calculate importance for each variable
    for (var in seq_len(p_predictors)) {
        # Permute this variable in OOB samples
        X_permuted <- X_pred
        if (length(oob_samples) > 0) {
            X_permuted[oob_samples, var] <- sample(X_permuted[oob_samples, var])
        }
        
        # Get predictions with permuted variable
        permuted_predictions <- JRF_predict_samples(tree_result, X_permuted,
                                                    oob_samples)
        permuted_error <- mean((y[oob_samples] - permuted_predictions)^2)
        
        # Importance is increase in error
        importance_scores[var] <- permuted_error - baseline_error
    }
    
    return(pmax(importance_scores, 0))  # Ensure non-negative
}

#' Predict Samples Using Tree
#' 
#' Makes predictions for specified samples using the tree structure.
#' 
#' @param tree_result Tree structure
#' @param X_pred Predictor matrix
#' @param sample_indices Sample indices to predict
#' 
#' @return Predictions
#' @keywords internal
#' @noRd
JRF_predict_samples <- function(tree_result, X_pred, sample_indices) {
    predictions <- numeric(length(sample_indices))
    
    for (i in seq_along(sample_indices)) {
        sample_idx <- sample_indices[i]
        predictions[i] <- JRF_predict_single_sample(tree_result, 
                                                    X_pred[sample_idx, ])
    }
    
    return(predictions)
}

#' Predict Single Sample Using Tree
#' 
#' Makes a prediction for a single sample using the tree structure.
#' 
#' @param node Current node in tree
#' @param sample_data Single sample data
#' 
#' @return Prediction
#' @keywords internal
#' @noRd
JRF_predict_single_sample <- function(node, sample_data) {
    if (node$is_terminal) {
        return(node$prediction[1])  # Return prediction for first class
    }
    
    if (sample_data[node$split_variable] <= node$split_threshold) {
        return(JRF_predict_single_sample(node$left_child, sample_data))
    } else {
        return(JRF_predict_single_sample(node$right_child, sample_data))
    }
}

#' Wrapper function for compatibility
#' 
#' @param X List of expression matrices
#' @param ntree Number of trees
#' @param mtry Number of variables sampled at each split
#' @param genes.name Gene names
#' 
#' @return JRF result
#' @keywords internal
#' @noRd
JRF_simplified <- function(X, ntree = 500, mtry = NULL, genes.name = NULL) {
    return(JRF_internal(X, ntree, mtry, genes.name))
}