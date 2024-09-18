import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression

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
    
    # Initialize an expression matrix (cells x genes) with random initial expression levels
    expression_matrix = np.random.uniform(0.1, 1.0, (num_cells, num_genes))
    
    # Step 3: Iteratively update the expression based on interactions
    for iteration in range(num_iterations):
        print(f"--- Iteration {iteration + 1} ---")
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
                    
                    # Show the model parameters (intercept and coefficients)
                    print(f"Gene: {gene}")
                    print(f"Model intercept: {model.intercept_}")
                    print(f"Model coefficients: {model.coef_}\n")
                else:
                    predicted_expression = X[0] * y[0]  # Handle the case with only one neighbor
                
                # Update the gene expression level with some added noise
                new_expression_levels[gene_idx] = predicted_expression + noise_level * np.random.normal()
                new_expression_levels[gene_idx] = max(0, new_expression_levels[gene_idx])  # Ensure no negative values
            
            # Update the expression matrix with new expression levels
            expression_matrix[cell, :] = new_expression_levels
    
    return expression_matrix, genes

# Step 4: Save simulated expression data to a file
def save_simulated_data(expression_matrix, genes, output_file='simulated_scRNAseq_expression.csv'):
    # Create DataFrame with cells as rows and genes as columns
    num_cells = expression_matrix.shape[0]
    cell_names = [f"Cell_{i+1}" for i in range(num_cells)]
    expression_df = pd.DataFrame(expression_matrix, columns=genes, index=cell_names)
    
    # Save to CSV
    expression_df.to_csv(output_file)
    print(f"Simulated scRNA-seq expression data saved to {output_file}")

# Main: Load interaction data and run regression-based simulation
interaction_file = './../data/genemania-interactions.txt'  # Replace with the correct path
interactions = load_interaction_data(interaction_file)

# Step 2: Simulate gene expression levels based on interactions
expression_matrix, gene_names = simulate_gene_expression_with_regression(interactions, num_cells=500, num_iterations=10, noise_level=0.01)

# Step 4: Save the simulated expression data
save_simulated_data(expression_matrix, gene_names)

