import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

# Step 1: Create a Gene Regulatory Network (GRN) from interaction data
def load_grn_from_file(interaction_file):
    # Load interaction file
    data = pd.read_csv(interaction_file, sep='\t')  # Adjust separator if needed (e.g., ',' for CSV)
    
    # Create a graph and add edges based on the interaction data
    G = nx.Graph()
    
    for index, row in data.iterrows():
        gene1 = row['Gene 1']
        gene2 = row['Gene 2']
        weight = row['Weight'] if 'Weight' in row else np.random.uniform(-1, 1)  # Assign weights if available
        G.add_edge(gene1, gene2, weight=weight)
    
    return G

# Step 2: Simulate Gene Expression using an ODE-based Model
def simulate_expression(grn, num_cells, num_timepoints, noise_level=0.1):
    gene_names = list(grn.nodes)  # Get the list of gene names
    num_genes = len(gene_names)
    
    # Initialize gene expression levels (all start with a baseline expression)
    expression = np.zeros((num_cells, num_genes))
    time_series = np.linspace(0, 10, num_timepoints)
    
    # Simulate dynamics for each cell independently
    for cell in range(num_cells):
        expression_levels = np.random.uniform(0.1, 1, num_genes)  # Random initial expression levels
        for t in time_series:
            for gene_idx in range(num_genes):
                gene = gene_names[gene_idx]  # Use the gene name instead of an index
                # Get regulation from neighboring genes in the GRN
                regulation = sum(grn.edges[gene, nbr]['weight'] * expression_levels[gene_names.index(nbr)]
                                 for nbr in grn.neighbors(gene) if grn.has_edge(gene, nbr))
                
                # Update expression level using a simple ODE
                expression_levels[gene_idx] += 0.01 * regulation + noise_level * np.random.normal()
                expression_levels[gene_idx] = max(0, expression_levels[gene_idx])  # No negative expression
                
        expression[cell, :] = expression_levels
    
    return expression, gene_names

# Step 3: Introduce scRNA-seq Noise (Dropout and Technical Variability)
def introduce_scrnaseq_noise(expression, dropout_rate=0.3, technical_noise_level=0.1):
    # Apply dropout (zero out some gene expressions)
    dropout_mask = np.random.binomial(1, 1 - dropout_rate, size=expression.shape)
    expression_with_dropout = expression * dropout_mask
    
    # Add technical noise (Gaussian noise)
    noise = np.random.normal(0, technical_noise_level, size=expression.shape)
    expression_with_noise = expression_with_dropout + noise
    
    # Ensure no negative values after adding noise
    expression_with_noise[expression_with_noise < 0] = 0
    
    return expression_with_noise

# Step 4: Run the Simulation with Interaction Data
def run_simulation(interaction_file, num_cells=1000, num_timepoints=5):
    # Load the GRN from the provided interaction file
    grn = load_grn_from_file(interaction_file)
    
    # Simulate gene expression dynamics
    expression, gene_names = simulate_expression(grn, num_cells, num_timepoints)
    
    # Introduce scRNA-seq-like noise (dropout and technical noise)
    expression_noisy = introduce_scrnaseq_noise(expression)
    
    return expression_noisy, gene_names

# Step 5: Save and Visualize the Data
def save_and_plot_expression(expression, gene_names, output_file="simulated_scrna_seq.csv"):
    # Create row names as 'cell1', 'cell2', etc.
    cell_names = [f'cell{i+1}' for i in range(expression.shape[0])]
    
    # Create a DataFrame with cell names as rows and gene names as columns
    expression_df = pd.DataFrame(expression, index=cell_names, columns=gene_names)
    
    # Save the expression data to a CSV file
    expression_df.to_csv(output_file)
    
    print(f"Simulation complete. Data saved to '{output_file}'")
    
    # Plot a heatmap of a subset of the data for visualization
    plt.figure(figsize=(10, 8))
    plt.imshow(expression_df.values[:50, :50], aspect='auto', cmap='viridis')
    plt.colorbar(label='Gene Expression')
    plt.title('Heatmap of Simulated scRNA-seq Data (First 50 cells and genes)')
    plt.xlabel('Genes')
    plt.ylabel('Cells')
    plt.show()

# Run the simulation with your interaction data file
interaction_file = './../data/genemania-interactions.txt'  # Replace with your interaction file path
expression_data, gene_names = run_simulation(interaction_file, num_cells=500, num_timepoints=10)
save_and_plot_expression(expression_data, gene_names)

