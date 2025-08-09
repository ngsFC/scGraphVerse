#' Generic functions for scGraphVerse
#'
#' This file contains all generic function definitions for the scGraphVerse.
#' These must be defined before any setMethod calls.
#'
#' @importFrom methods setGeneric
#' @keywords internal

#' Structure learning for count data using PC algorithms
#'
#' This function performs structure learning for count data using various 
#' PC algorithms adapted for different distributional assumptions including 
#' Poisson, Negative Binomial, and Zero-Inflated Negative Binomial models.
#'
#' @param x A matrix of count data or SummarizedExperiment object
#' @param ... Additional arguments passed to specific methods
#' @return An adjacency matrix or SummarizedExperiment with adjacency matrix
#' @export
#' @rdname PCzinb
setGeneric("PCzinb", function(x, ...) standardGeneric("PCzinb"))