
# Gene Regulatory Network (GRN) Analysis

This repository contains a collection of R scripts and markdown files used to perform various analyses related to network generation, adjacency matrices, and ROC (Receiver Operating Characteristic) curve plotting. The files in this repository are designed to work together for the analysis and visualization of network properties.

## Overview of Files

- **grnet.Rmd**: This R Markdown file likely serves as the main document that ties together all the scripts in this repository. It likely contains explanations, visualizations, and code executions for generating results from the other scripts.

- **cutoff_adjacency.R**: This script is responsible for creating an adjacency matrix based on a specific cutoff threshold. This is a key component of network generation, where the connections between nodes are defined based on their relationships and a cutoff value.

- **pscores.R**: This script calculate and plot TPR, FPR, PRECISION, F1, ACCURACY.

- **plotROC.R**: This script handles the plotting of ROC curves, which are used to evaluate the performance of the models. In a network analysis context, ROC curves might be used to assess the quality of link prediction models.

- **generate_adjacency.R**: This script generates the adjacency matrix for a network. It may either rely on input data or be part of a process that constructs a network from scratch based on node and edge criteria.

- **simmetric.R**: This script might calculate symmetric properties of an adjacency matrix or network, which is important in undirected network analysis where relationships between nodes are mutual.

- **plotg.R**: A script for plotting graphs (networks). This file likely provides functions to visualize the networks generated or analyzed by the other scripts.

- **earlyj.R**: concatenate and create a single matrix from a list of n matrices for early integration steps.

- **dropo.R**: dropout script.

- **compare_consensus.R**: This script compares different networks. Consensus matrices are used to summarize the agreement between different network structures results.

## Usage

To reproduce the analysis or explore the network generation process, you can run the provided scripts in R. Here's a general workflow to get started:

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/your-repo.git
   cd your-repo
   ```

2. Install the necessary R packages:
   Ensure you have all the required packages installed. You can install the necessary libraries by adding the following to your R session:
   ```R
   install.packages(c("ggplot2", "igraph", "pROC"))  # Add more as needed
   ```

3. Run the R scripts individually or use the `grnet.Rmd` file to generate the full analysis:
   ```R
   rmarkdown::render("grnet.Rmd")
   ```

4. Customization:
   Modify the cutoff thresholds, input data, or network parameters in the individual scripts to suit your needs.

## Contributions

If you would like to contribute to this project, feel free to submit a pull request or report issues. Contributions such as bug fixes, enhancements, and suggestions are welcome!

## License

This project is licensed under the MIT License.

## Founding
"National Centre for HPC, Big Data and Quantum Computing" - CN00000013 - CUP B93C22999629996
