import pandas as pd
import numpy as np

# Step 1: Load gene interaction data from the file
def load_interaction_data(file):
    # Load interaction file
    data = pd.read_csv(file, sep='\t')
    return data

# Step 2: Simulate gene expression using regression-based model
def simulate_gene_expression(interactions, num_cells=1000, num_iterations=10, noise_level=0.01):
    # Extract unique genes from the interaction data
    genes = list(set(interactions['Gene 1']).union(set(interactions['Gene 2'])))
    num_genes = len(genes)
    
    # Create a mapping from gene name to index
    gene_index = {gene: idx for idx, gene in enumerate(genes)}
    
    # Initialize an expression matrix (cells x genes) with random initial expression levels
    expression_matrix = np.random.uniform(0.1, 1.0, (num_cells, num_genes))
    
    # Step 3: Iteratively update the expression based on interactions
    for iteration in range(num_iterations):
        for cell in range(num_cells):
            new_expression_levels = np.copy(expression_matrix[cell, :])  # Keep the current expression levels
            
            # For each gene, update its expression based on its interacting genes
            for gene_idx, gene in enumerate(genes):
                interacting_genes = interactions[interactions['Gene 1'] == gene]
                total_influence = 0
                
                # Apply interaction weights from neighbors
                for _, row in interacting_genes.iterrows():
                    neighbor_gene = row['Gene 2']
                    weight = row['Weight']
                    neighbor_idx = gene_index[neighbor_gene]
                    total_influence += weight * expression_matrix[cell, neighbor_idx]
                
                # Update expression using the weighted sum of neighbors plus a small noise term
                new_expression_levels[gene_idx] = total_influence + noise_level * np.random.normal()
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
expression_matrix, gene_names = simulate_gene_expression(interactions, num_cells=500, num_iterations=10, noise_level=0.01)

# Step 4: Save the simulated expression data
save_simulated_data(expression_matrix, gene_names)

