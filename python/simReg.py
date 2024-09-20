import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy.stats import nbinom
import matplotlib.pyplot as plt
import seaborn as sns

# Step 1: Load gene interaction data from the file
def load_interaction_data(file):
    # Load interaction file as tab-separated
    data = pd.read_csv(file, sep='\t')
    return data

# Step 2: Simulate gene expression using a regression-based model
def simulate_gene_expression_with_regression(interactions, num_cells=1000, num_iterations=10, noise_level=0.01):
    # Extract unique genes from the interaction data
    genes = list(set(interactions['Gene 1']).union(set(interactions['Gene 2'])))
    num_genes = len(genes)
    
    # Create a mapping from gene name to index
    gene_index = {gene: idx for idx, gene in enumerate(genes)}
    
    # Initialize an expression matrix (cells x genes) with random initial expression levels using Gamma distribution
    expression_matrix = np.random.gamma(shape=2.0, scale=2.0, size=(num_cells, num_genes))
    
    # Step 3: Iteratively update the expression based on interactions
    for iteration in range(num_iterations):
        for cell in range(num_cells):
            new_expression_levels = np.copy(expression_matrix[cell, :])  # Keep the current expression levels
            
            # For each gene, update its expression based on its interacting genes using linear regression
            for gene_idx, gene in enumerate(genes):
                interacting_genes = interactions[interactions['Gene 1'] == gene]
                
                # If the gene has no interacting neighbors, skip
                if interacting_genes.empty:
                    continue
                
                # Prepare the input for regression: neighbor expression levels and their weights
                X = []
                y = []
                for _, row in interacting_genes.iterrows():
                    neighbor_gene = row['Gene 2']
                    weight = row['Weight']
                    neighbor_idx = gene_index[neighbor_gene]
                    X.append(expression_matrix[cell, neighbor_idx])
                    y.append(weight)
                
                # Convert lists to numpy arrays for regression
                X = np.array(X).reshape(-1, 1)  # Features (neighbor expression levels)
                y = np.array(y)  # Targets (weights)
                
                # Perform linear regression (predict the gene expression from its neighbors)
                if len(X) > 1:  # Ensure there is more than one neighbor
                    model = LinearRegression()
                    model.fit(X, y)
                    predicted_expression = model.predict(X).sum()
                else:
                    predicted_expression = X[0] * y[0]  # Handle the case with only one neighbor
                
                # Update the gene expression level with some added noise
                new_expression_levels[gene_idx] = predicted_expression + noise_level * np.random.normal()
                new_expression_levels[gene_idx] = max(0, new_expression_levels[gene_idx])  # Ensure no negative values
            
            # Update the expression matrix with new expression levels
            expression_matrix[cell, :] = new_expression_levels
    
    return expression_matrix, genes

# Step 4: Generate counts using Negative Binomial distribution based on regression model predictions
def generate_counts_from_nb(expression_matrix, dispersion=0.5):
    num_cells, num_genes = expression_matrix.shape
    count_matrix = np.zeros_like(expression_matrix, dtype=int)
    
    for cell in range(num_cells):
        for gene_idx in range(num_genes):
            # Mean expression predicted by the regression model
            mu = expression_matrix[cell, gene_idx]
            
            # Convert mean expression to Negative Binomial count
            if mu > 0:
                r = dispersion  # Dispersion parameter
                p = mu / (mu + r)  # Calculate probability for NB
                count_matrix[cell, gene_idx] = nbinom.rvs(r, p)
    
    return count_matrix

# Step 5: Simulate dropout
def apply_dropout(count_matrix, dropout_rate=0.2):
    dropout_mask = np.random.binomial(1, 1 - dropout_rate, size=count_matrix.shape)  # Generate dropout mask
    count_matrix *= dropout_mask  # Apply dropout: set some counts to 0
    return count_matrix

# Step 6: Save simulated count data to a file
def save_simulated_data(count_matrix, genes, output_file='simulated_scRNAseq_counts.csv'):
    num_cells = count_matrix.shape[0]
    cell_names = [f"Cell_{i+1}" for i in range(num_cells)]
    expression_df = pd.DataFrame(count_matrix, columns=genes, index=cell_names)
    expression_df.to_csv(output_file)
    print(f"Simulated scRNA-seq count data saved to {output_file}")

# Step 7: Plot results for visualization
def plot_simulation_results(count_matrix):
    counts = count_matrix.flatten()
    
    # Plot a histogram of count values
    plt.figure(figsize=(10, 6))
    plt.hist(counts[counts > 0], bins=50, color='blue', alpha=0.7)
    plt.title('Distribution of Gene Expression Counts')
    plt.xlabel('Count')
    plt.ylabel('Frequency')
    plt.show()
    
    # Plot a heatmap of a sample of the count matrix (e.g., first 50 cells and genes)
    plt.figure(figsize=(12, 8))
    sns.heatmap(count_matrix[:50, :50], cmap='viridis')
    plt.title('Heatmap of Expression Counts (First 50 Cells and Genes)')
    plt.xlabel('Genes')
    plt.ylabel('Cells')
    plt.show()

# Main: Load interaction data and run the simulation
interaction_file = './../data/genemania-interactions.txt'  # Replace with the correct path
interactions = load_interaction_data(interaction_file)

# Step 2: Simulate gene expression levels based on interactions using regression
expression_matrix, gene_names = simulate_gene_expression_with_regression(interactions, num_cells=500, num_iterations=10, noise_level=0.01)

# Step 4: Generate counts using Negative Binomial distribution
count_matrix = generate_counts_from_nb(expression_matrix, dispersion=0.5)

# Step 5: Apply dropout to simulate scRNA-seq technical noise
count_matrix = apply_dropout(count_matrix, dropout_rate=0.2)

# Step 6: Save the simulated count data
save_simulated_data(count_matrix, gene_names)

# Step 7: Plot the simulation results
plot_simulation_results(count_matrix)
