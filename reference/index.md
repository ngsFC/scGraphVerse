# Package index

## Network Inference

Core functions for inferring gene regulatory networks from single-cell
data

- [`infer_networks()`](https://ngsfc.github.io/scGraphVerse/reference/infer_networks.md)
  : Infer Gene Regulatory Networks from Expression Matrices
- [`create_consensus()`](https://ngsfc.github.io/scGraphVerse/reference/create_consensus.md)
  : Create a Consensus Adjacency Matrix from Multiple Networks
- [`compare_consensus()`](https://ngsfc.github.io/scGraphVerse/reference/compare_consensus.md)
  : Compare Consensus and Reference Graphs or STRINGdb Networks
- [`PCzinb()`](https://ngsfc.github.io/scGraphVerse/reference/PCzinb.md)
  : Structure learning for count data using PC algorithms

## Data Preparation

Functions for preparing and organizing input data

- [`create_mae()`](https://ngsfc.github.io/scGraphVerse/reference/create_mae.md)
  : Create MultiAssayExperiment from Multiple Single-Cell Datasets
- [`build_network_se()`](https://ngsfc.github.io/scGraphVerse/reference/build_network_se.md)
  : Create a SummarizedExperiment for Network Storage

## Network Analysis

Tools for analyzing network structure and properties

- [`compute_topology_metrics()`](https://ngsfc.github.io/scGraphVerse/reference/compute_topology_metrics.md)
  : Compute Network Topological Properties
- [`compute_community_metrics()`](https://ngsfc.github.io/scGraphVerse/reference/compute_community_metrics.md)
  : Compute Community Assignment Similarity Metrics
- [`community_similarity()`](https://ngsfc.github.io/scGraphVerse/reference/community_similarity.md)
  : Compare Community Assignments and Topological Properties
- [`classify_edges()`](https://ngsfc.github.io/scGraphVerse/reference/classify_edges.md)
  : Classify Edges as TP, FP, or FN
- [`edge_mining()`](https://ngsfc.github.io/scGraphVerse/reference/edge_mining.md)
  : Edge Mining of Gene Interactions Using PubMed
- [`community_path()`](https://ngsfc.github.io/scGraphVerse/reference/community_path.md)
  : Community Detection and Pathway Enrichment Analysis

## Visualization

Plotting functions for networks and analysis results

- [`plotg()`](https://ngsfc.github.io/scGraphVerse/reference/plotg.md) :
  Visualize Graphs from Adjacency Matrices
- [`plot_network_comparison()`](https://ngsfc.github.io/scGraphVerse/reference/plot_network_comparison.md)
  : Visualize Network Comparison
- [`plot_community_comparison()`](https://ngsfc.github.io/scGraphVerse/reference/plot_community_comparison.md)
  : Visualize Community and Topology Comparison
- [`plotROC()`](https://ngsfc.github.io/scGraphVerse/reference/plotROC.md)
  : Plot ROC Curves for Inferred Networks

## Network Utilities

Helper functions for network manipulation and processing

- [`symmetrize()`](https://ngsfc.github.io/scGraphVerse/reference/symmetrize.md)
  : Symmetrize Adjacency Matrices in a SummarizedExperiment
- [`generate_adjacency()`](https://ngsfc.github.io/scGraphVerse/reference/generate_adjacency.md)
  : Generate Adjacency Matrices from Gene Interaction Tables
- [`cutoff_adjacency()`](https://ngsfc.github.io/scGraphVerse/reference/cutoff_adjacency.md)
  : Threshold Adjacency Matrices Based on Shuffled Network Quantiles
- [`stringdb_adjacency()`](https://ngsfc.github.io/scGraphVerse/reference/stringdb_adjacency.md)
  : Build Adjacency Matrices for Physical Interactions from STRING (POST
  API)
- [`selgene()`](https://ngsfc.github.io/scGraphVerse/reference/selgene.md)
  : Select Top Expressed Genes from Single-Cell Data
- [`pscores()`](https://ngsfc.github.io/scGraphVerse/reference/pscores.md)
  : Compute Performance Scores for Predicted Adjacency Matrices

## Setup and Configuration

Functions for initializing Python dependencies and downloading reference
data

- [`init_py()`](https://ngsfc.github.io/scGraphVerse/reference/init_py.md)
  : Initialize Python Environment for GRNBoost2
- [`download_Atlas()`](https://ngsfc.github.io/scGraphVerse/reference/download_Atlas.md)
  : Download and Load an RDS File from a URL

## Example Data

Toy datasets for testing and examples

- [`toy_counts`](https://ngsfc.github.io/scGraphVerse/reference/toy_counts.md)
  : Toy MultiAssayExperiment for Network Inference
- [`toy_adj_matrix`](https://ngsfc.github.io/scGraphVerse/reference/toy_adj_matrix.md)
  : Toy adjacency matrix for examples
- [`zinb_simdata()`](https://ngsfc.github.io/scGraphVerse/reference/zinb_simdata.md)
  : Simulate Zero-Inflated Negative Binomial (ZINB) Count Matrices with
  Sequencing Depth

## Internal Functions

Internal helper functions (not typically called directly by users)

- [`earlyj()`](https://ngsfc.github.io/scGraphVerse/reference/earlyj.md)
  : Modify Cell Names and Combine Datasets
- [`nb.loglik()`](https://ngsfc.github.io/scGraphVerse/reference/nb.loglik.md)
  : Log-likelihood of the negative binomial model Given a vector of
  counts, this function computes the sum of the log-probabilities of the
  counts under a negative binomial (NB) model. The NB distribution is
  parametrized by two parameters: the mean value and the dispersion of
  the negative binomial distribution
- [`nb.loglik.dispersion()`](https://ngsfc.github.io/scGraphVerse/reference/nb.loglik.dispersion.md)
  : Log-likelihood of negative binomial model, for a fixed dispersion
  parameter
- [`nb.loglik.regression()`](https://ngsfc.github.io/scGraphVerse/reference/nb.loglik.regression.md)
  : log-likelihood of the NB regression model
- [`nb.loglik.regression.gradient()`](https://ngsfc.github.io/scGraphVerse/reference/nb.loglik.regression.gradient.md)
  : Gradient of the log-likelihood of the NB regression model
- [`nb.OptimizeDispersion()`](https://ngsfc.github.io/scGraphVerse/reference/nb.OptimizeDispersion.md)
  : (NB) model. The NB distribution is parametrized by two parameters:
  the mean value and the dispersion of the negative binomial
  distribution
- [`nb.regression.parseModel()`](https://ngsfc.github.io/scGraphVerse/reference/nb.regression.parseModel.md)
  : Parse ZINB regression model
- [`zinb0.noT()`](https://ngsfc.github.io/scGraphVerse/reference/zinb0.noT.md)
  : Structure learning with zero-inflated negative binomial model (mean
  only)
- [`zinb1.noT()`](https://ngsfc.github.io/scGraphVerse/reference/zinb1.noT.md)
  : Structure learning with zero-inflated negative binomial model
- [`zinbOptimizeDispersion()`](https://ngsfc.github.io/scGraphVerse/reference/zinbOptimizeDispersion.md)
  : (ZINB) model. The ZINB distribution is parametrized by three
  parameters: the mean value and the dispersion of the negative binomial
  distribution, and the probability of the zero component.
- [`zinb.regression.parseModel()`](https://ngsfc.github.io/scGraphVerse/reference/zinb.regression.parseModel.md)
  : Parse ZINB regression model
