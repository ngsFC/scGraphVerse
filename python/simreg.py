import pandas as pd
import numpy as np
from scipy.stats import nbinom

# Step 1: Load gene interaction data from the file
def load_interaction_data(file):
    # Load interaction file
    data = pd.read_csv(file, sep='\t')
    return data

# Step 2: Simulate mean expression based on regression models
def simulate_mean_expression(interactions, num_cells=1000):
    genes = list(set(interactions['Gene 1']).union(set(interactions['Gene 2'])))
    num_genes = len(genes)
    gene_index = {gene: idx for idx, gene in enumerate(genes)}  # Map gene names to indices

    # Initialize an empty mean expression matrix (cells x genes)
    mean_expression_matrix = np.zeros((num_cells, num_genes))

    # Iterate through each cell and simulate the mean expression of each gene
    for cell in range(num_cells):
        mean_expression_levels = np.random.uniform(0.1, 1, num_genes)  # Random initial expression for genes
        for gene_idx, gene in enumerate(genes):
            interacting_genes = interactions[interactions['Gene 1'] == gene]
            if not interacting_genes.empty:
                total_influence = 0
                for _, row in interacting_genes.iterrows():
                    regulator_gene = row['Gene 2']
                    weight = row['Weight']
                    regulator_idx = gene_index[regulator_gene]
                    total_influence += weight * mean_expression_levels[regulator_idx]
                
                # Add noise and update mean expression level
                mean_expression_levels[gene_idx] = total_influence + np.random.normal(0, 0.1)
        
        # Store the simulated mean expression levels for this cell
        mean_expression_matrix[cell, :] = mean_expression_levels
    
    return mean_expression_matrix, genes

# Step 3: Simulate scRNA-seq count data using Negative Binomial distribution
def simulate_count_data(mean_expression_matrix, num_cells, library_size_mean=5000, library_size_sd=1000, dispersion=1.0):
    num_genes = mean_expression_matrix.shape[1]
    
    # Generate a library size for each cell (total counts per cell)
    library_sizes = np.random.normal(library_size_mean, library_size_sd, num_cells).astype(int)
    library_sizes[library_sizes < 0] = library_size_mean  # Ensure no negative library sizes

    # Initialize an empty count matrix (cells x genes)
    count_matrix = np.zeros((num_cells, num_genes), dtype=int)

    # For each cell, generate gene expression counts using Negative Binomial
    for cell in range(num_cells):
        total_expression = mean_expression_matrix[cell, :]
        
        # Scale the mean expression by the library size
        scaled_mean_expression = total_expression / np.sum(total_expression) * library_sizes[cell]
        
        # Simulate counts for each gene using Negative Binomial distribution
        for gene_idx in range(num_genes):
            mean_exp = scaled_mean_expression[gene_idx]
            if mean_exp > 0:
                # Variance = mean + mean^2 / dispersion
                size = 1 / dispersion  # Inverse dispersion parameter
                prob = size / (size + mean_exp)
                count = nbinom.rvs(size, prob)
                count_matrix[cell, gene_idx] = count
    
    return count_matrix

# Step 4: Save simulated count data to a file
def save_simulated_data(count_matrix, genes, output_file='simulated_scRNAseq_counts.csv'):
    # Create DataFrame
    num_cells = count_matrix.shape[0]
    cell_names = [f"Cell_{i+1}" for i in range(num_cells)]
    expression_df = pd.DataFrame(count_matrix, columns=genes, index=cell_names)
    
    # Save to CSV
    expression_df.to_csv(output_file)
    print(f"Simulated scRNA-seq count data saved to {output_file}")

# Main: Load interaction data and run simulation
interaction_file = './../data/genemania-interactions.txt'
interactions = load_interaction_data(interaction_file)

# Step 2: Simulate mean expression levels using the regression model
mean_expression_matrix, gene_names = simulate_mean_expression(interactions, num_cells=500)

# Step 3: Simulate scRNA-seq count data using Negative Binomial distribution
count_matrix = simulate_count_data(mean_expression_matrix, num_cells=500)

# Step 4: Save the simulated count data
save_simulated_data(count_matrix, gene_names)

