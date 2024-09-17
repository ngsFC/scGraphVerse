import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# Load the T-cell interaction network file
tcell_interactions_path = './../data/Tcell_interactions.tsv'
tcell_interactions = pd.read_csv(tcell_interactions_path, sep="\t")

# Create a directed graph using the interaction data
G = nx.DiGraph()

# Add edges with combined score as edge weight
for index, row in tcell_interactions.iterrows():
    G.add_edge(row['#node1'], row['node2'], weight=row['combined_score'])

# Function to simulate gene expression based on the network
def simulate_expression(G, num_cells=100, num_iterations=10, noise_level=0.1):
    nodes = list(G.nodes())
    num_genes = len(nodes)
    
    # Initialize random gene expression matrix (cells x genes)
    expression_matrix = np.random.rand(num_cells, num_genes)

    # Simulate gene expression based on network interactions
    for iteration in range(num_iterations):
        for i in range(num_genes):
            gene = nodes[i]
            # Get regulators of the gene (predecessors in the graph)
            regulators = list(G.predecessors(gene))
            if regulators:
                regulator_indices = [nodes.index(r) for r in regulators]
                weights = [G.edges[r, gene]['weight'] for r in regulators]
                
                # Update gene expression based on regulators
                expression_matrix[:, i] = np.dot(expression_matrix[:, regulator_indices], weights) + \
                                          np.random.normal(0, noise_level, num_cells)
    
    # Ensure non-negative values (mimicking expression values)
    expression_matrix = np.clip(expression_matrix, 0, None)
    
    return expression_matrix

# Step 1: Simulate gene expression for 100 cells and 10 iterations
simulated_expression = simulate_expression(G, num_cells=100, num_iterations=10)

# Step 2: Apply dropout to simulate scRNA-seq-like data
def apply_dropout(expression_matrix, dropout_rate=0.3):
    dropout_mask = np.random.rand(*expression_matrix.shape) > dropout_rate
    return expression_matrix * dropout_mask

# Step 3: Apply noise to the simulated gene expression
def add_noise(expression_matrix, noise_level=0.1):
    noise = np.random.normal(0, noise_level, expression_matrix.shape)
    return expression_matrix + noise

# Apply dropout
expression_with_dropout = apply_dropout(simulated_expression, dropout_rate=0.3)

# Apply noise
expression_noisy = add_noise(expression_with_dropout, noise_level=0.1)

# Convert to pandas dataframe for easier manipulation
nodes = list(G.nodes())
expression_noisy_df = pd.DataFrame(expression_noisy, columns=nodes)

# Add row names as 'cell1', 'cell2', ..., 'celln'
expression_noisy_df.index = [f'cell{i+1}' for i in range(expression_noisy_df.shape[0])]

# Save the noisy expression data with row names
expression_noisy_df.to_csv("simulated_expression_tcell_with_dropout.csv", index=True)

# Display the first few rows of the simulated expression data
print(expression_noisy_df.head())

